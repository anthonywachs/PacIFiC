#include <DS_HeatTransfer.hh>
#include <FV_DomainAndFields.hh>
#include <FV_DiscreteField.hh>
#include <FV_DomainBuilder.hh>
#include <DS_HeatTransferSystem.hh>
#include <FV_SystemNumbering.hh>
#include <FV_Mesh.hh>
#include <FV_TimeIterator.hh>
#include <MAC.hh>
#include <MAC_Root.hh>
#include <MAC_Error.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_Vector.hh>
#include <MAC_BoolArray2D.hh>
#include <MAC_Communicator.hh>
#include <MAC_Exec.hh>
#include <MAC_Application.hh>
#include <intVector.hh>
#include <LA_Vector.hh>
#include <LA_SeqMatrix.hh>
#include <LA_MatrixIterator.hh>
#include <PAC_Misc.hh>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <DS_AllRigidBodies.hh>

//---------------------------------------------------------------------------
DS_HeatTransfer*
DS_HeatTransfer:: create( MAC_Object* a_owner,
		MAC_ModuleExplorer const* exp,
                struct DS2HE const& transfer )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_HeatTransfer:: create" ) ;
   MAC_CHECK_PRE( exp != 0 ) ;

   DS_HeatTransfer* result =
                        new DS_HeatTransfer( a_owner, exp, transfer ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;

   return( result ) ;

}

//---------------------------------------------------------------------------
DS_HeatTransfer:: DS_HeatTransfer( MAC_Object* a_owner,
											  MAC_ModuleExplorer const* exp,
				    			  			  struct DS2HE const& fromDS )
//---------------------------------------------------------------------------
   : MAC_Object( a_owner )
   , PAC_ComputingTime("Solver")
   , TF ( fromDS.dom_->discrete_field( "temperature" ) )
   , UF ( 0 )
   , TF_DS_ERROR( 0 )
   , GLOBAL_EQ( 0 )
   , rho( fromDS.rho_ )
   , AdvectionScheme ( fromDS.AdvectionScheme_ )
   , AdvectionTimeAccuracy ( fromDS.AdvectionTimeAccuracy_ )
	, ViscousStressOrder ( fromDS.ViscousStressOrder_ )
   , heat_capacity( exp->double_data( "Heat_capacity") )
   , thermal_conductivity( exp->double_data( "Thermal_conductivity") )
   , b_bodyterm ( false )
   , is_solids ( fromDS.is_solids_ )
	, is_NSwithHE (fromDS.is_NSwithHE_ )
	, is_stressCal (fromDS.is_stressCal_ )
	, stressCalFreq ( fromDS.stressCalFreq_ )
	, b_restart ( fromDS.b_restart_ )
	, is_par_motion ( fromDS.is_par_motion_ )
	, allrigidbodies ( fromDS.allrigidbodies_ )
{
   MAC_LABEL( "DS_HeatTransfer:: DS_HeatTransfer" ) ;

   MAC_ASSERT( TF->discretization_type() == "centered" ) ;
   MAC_ASSERT( TF->storage_depth() == 5 ) ;

   // Call of MAC_Communicator routine to set the rank of each proces and
   // the number of processes during execution
   macCOMM = MAC_Exec::communicator();
   my_rank = macCOMM->rank();
   nb_procs = macCOMM->nb_ranks();
   is_master = 0;
   is_iperiodic[0] = false;
   is_iperiodic[1] = false;
   is_iperiodic[2] = false;

   // Timing routines
   if ( my_rank == is_master ) {
     CT_set_start();
     SCT_insert_app("Objects_Creation");
     SCT_set_start("Objects_Creation");
   }

   // Get space dimension
   dim = TF->primary_grid()->nb_space_dimensions() ;
   nb_comps = TF->nb_components() ;

   if ( dim == 1 ) {
     string error_message="Space dimension should either 2 or 3";
     MAC_Error::object()->raise_bad_data_value(exp,
		  													 "nb_space_dimensions",
															 error_message );
   }

   // Create the Direction Splitting subcommunicators
   create_DS_subcommunicators();

   // Read with or without body term
   if ( exp->has_entry( "BodyTerm" ) )
      b_bodyterm = exp->bool_data( "BodyTerm" ) ;

   // Periodic boundary condition check
   periodic_comp = TF->primary_grid()->get_periodic_directions();
   is_iperiodic[0] = periodic_comp->operator()( 0 );
   is_iperiodic[1] = periodic_comp->operator()( 1 );
   if(dim >2) {
      is_iperiodic[2] = periodic_comp->operator()( 2 );
   }

   is_par_motion = fromDS.is_par_motion_;

	if (is_NSwithHE) {
		UF = fromDS.dom_->discrete_field( "velocity" );
	}

   // Create structure to input in the solver system

   // Build the matrix system
   MAC_ModuleExplorer* se =
								exp->create_subexplorer( 0,"DS_HeatTransferSystem" ) ;
   GLOBAL_EQ = DS_HeatTransferSystem::create( this, se, TF ) ;
   se->destroy() ;

   // Timing routines
   if ( my_rank == is_master ) {
     SCT_insert_app("Matrix_Assembly&Initialization");
     SCT_insert_app("Solver first step");
     SCT_insert_app("Solver x solution");
     SCT_insert_app("Solver y solution");
     SCT_insert_app("Solver z solution");
     SCT_get_elapsed_time("Objects_Creation");
   }
}




//---------------------------------------------------------------------------
DS_HeatTransfer:: ~DS_HeatTransfer( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_HeatTransfer:: ~DS_HeatTransfer" ) ;

   free_DS_subcommunicators() ;

}

//---------------------------------------------------------------------------
void
DS_HeatTransfer:: do_before_time_stepping( FV_TimeIterator const* t_it,
      	std::string const& basename )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_HeatTransfer:: do_before_time_stepping" ) ;

   if ( my_rank == is_master ) SCT_set_start("Matrix_Assembly&Initialization");

   allocate_mpi_variables();

   // Necessary especially for cases with non-zero field initiallization
   if (b_restart == false) ugradu_initialization ( );

	// Calculate row index for each field in each direction`
	calculate_row_indexes ( );

   // Generate solid particles if required
	if (is_solids) {
		// Build void frac and intersection variable
		allrigidbodies->build_solid_variables_on_fluid_grid(TF);
		// Compute void fraction for temperature field
		allrigidbodies->compute_void_fraction_on_grid(TF, false);
		// Compute intersection with RB for temperature field
		allrigidbodies->compute_grid_intersection_with_rigidbody(TF, false);
		if (my_rank == 0)
			cout << "HE: Finished void fraction and grid intersection... \n" << endl;
	}

   if (is_solids) {
      nodes_temperature_initialization(0);
      nodes_temperature_initialization(1);
      nodes_temperature_initialization(3);
      if (dim == 3) nodes_temperature_initialization(4);
   }

   // Assemble 1D tridiagonal matrices and schur complement calculation
   assemble_temperature_and_schur(t_it);

   if (my_rank == 0)
		cout << "HE: Finished assembling pre-coefficient matrix... \n" << endl;

   if ( my_rank == is_master ) SCT_get_elapsed_time( "Matrix_Assembly&Initialization" );

}

//---------------------------------------------------------------------------
void
DS_HeatTransfer:: do_before_inner_iterations_stage(
	FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_HeatTransfer:: do_before_inner_iterations_stage" ) ;
	//
   // if ((is_par_motion) && (is_solids)) {
   //    nodes_temperature_initialization(0);
   //    nodes_temperature_initialization(1);
   //    nodes_temperature_initialization(3);
   //    if (dim == 3) nodes_temperature_initialization(4);
	//
   //    // Assemble 1D tridiagonal matrices and schur complement calculation
   //    assemble_temperature_and_schur(t_it);
   // }

}

//---------------------------------------------------------------------------
void
DS_HeatTransfer:: do_one_inner_iteration( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_HeatTransfer:: do_one_inner_iteration" ) ;

   // Solve heat equation using direction splitting
   HeatEquation_DirectionSplittingSolver(t_it);

}

//---------------------------------------------------------------------------
void
DS_HeatTransfer:: do_after_inner_iterations_stage(
	FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_HeatTransfer:: do_after_inner_iterations_stage" ) ;

	// Compute fluid particle interaction
	if (is_stressCal && (t_it->iteration_number() % stressCalFreq == 0)) {
		allrigidbodies->compute_temperature_gradient_for_allRB(ViscousStressOrder);
	}

   // Compute temperature change over the time step
   double temperature_time_change = compute_DS_temperature_change()
   	/ t_it->time_step() ;
   if ( my_rank == is_master )
     cout << "Temperature change = " <<
     	MAC::doubleToString( ios::scientific, 5, temperature_time_change )
	<< endl;

}

//---------------------------------------------------------------------------
void
DS_HeatTransfer:: do_after_time_stepping( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_HeatTransfer:: do_after_time_stepping" ) ;

	GLOBAL_EQ->display_debug();
	output_l2norm();

	// Elapsed time by sub-problems
   if ( my_rank == is_master ) {
     double cputime = CT_get_elapsed_time();
	  cout << endl
	  		 << "========================================================" << endl
			 << "                Heat Transfer Problem                   " << endl
			 << "========================================================" << endl;
     write_elapsed_time_smhd(cout,cputime,"Computation time");
     SCT_get_summary(cout,cputime);
   }

	deallocate_mpi_variables();
}

//---------------------------------------------------------------------------
void
DS_HeatTransfer:: do_additional_savings( FV_TimeIterator const* t_it,
      	int const& cycleNumber )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_HeatTransfer:: do_additional_savings" ) ;

}




//----------------------------------------------------------------------
void
DS_HeatTransfer::write_output_field()
//----------------------------------------------------------------------
{

  ofstream outputFile ;

  std::ostringstream os2;
  os2 << "./DS_results/intersection_data.csv";
  std::string filename = os2.str();
  outputFile.open(filename.c_str());

  size_t i,j,k;
  outputFile << "x,y,z,void_frac,left,lv"
             << ",right,rv"
             << ",bottom,bov"
             << ",top,tv"
             << ",behind,bev"
             << ",front,fv" << endl;

  size_t_vector min_index(dim,0);
  size_t_vector max_index(dim,0);

  size_t_array2D* void_frac = (is_solids) ?
                    allrigidbodies->get_void_fraction_on_grid(TF) : 0;
  size_t_array2D* intersect_vector = (is_solids) ?
                    allrigidbodies->get_intersect_vector_on_grid(TF) : 0;
  doubleArray2D* intersect_distance = (is_solids) ?
                    allrigidbodies->get_intersect_distance_on_grid(TF) : 0;

  for (size_t comp=0;comp<nb_comps;comp++) {
     // Get local min and max indices
     for (size_t l=0;l<dim;++l) {
		  min_index(l) = 0;
        max_index(l) = TF->get_local_nb_dof( comp, l );
     }

     size_t local_min_k = 0;
     size_t local_max_k = 1;

     if (dim == 3) {
        local_min_k = min_index(2);
        local_max_k = max_index(2);
     }

     for (i=min_index(0);i<=max_index(0);++i) {
        double xC = TF->get_DOF_coordinate( i, comp, 0 ) ;
        for (j=min_index(1);j<=max_index(1);++j) {
           double yC = TF->get_DOF_coordinate( j, comp, 1 ) ;
           for (k=local_min_k;k<=local_max_k;++k) {
              double zC = 0.;
              if (dim == 3) zC = TF->get_DOF_coordinate( k, comp, 2 ) ;
				  size_t p = TF->DOF_local_number(i,j,k,comp);

              outputFile << xC << "," << yC << "," << zC
              << "," << void_frac->operator()(p,0)
              << "," << intersect_vector->operator()(p,0)
              << "," << intersect_distance->operator()(p,0)
              << "," << intersect_vector->operator()(p,1)
              << "," << intersect_distance->operator()(p,1)
              << "," << intersect_vector->operator()(p,2)
              << "," << intersect_distance->operator()(p,2)
              << "," << intersect_vector->operator()(p,3)
              << "," << intersect_distance->operator()(p,3)
              << "," << intersect_vector->operator()(p,4)
              << "," << intersect_distance->operator()(p,4)
              << "," << intersect_vector->operator()(p,5)
              << "," << intersect_distance->operator()(p,5)
              << endl;
           }
        }
     }
  }
  outputFile.close();
}

//---------------------------------------------------------------------------
double
DS_HeatTransfer:: bodyterm_value ( double const& xC
											, double const& yC
											, double const& zC)
//---------------------------------------------------------------------------
{
   double bodyterm = 0.;
   if (b_bodyterm) {
      if (dim == 2) {
         bodyterm = 2. * MAC::pi() * MAC::pi() * MAC::sin( MAC::pi() * xC )
                                               * MAC::sin( MAC::pi() * yC );
      } else if (dim == 3) {
         bodyterm = 3. * MAC::pi() * MAC::pi() * MAC::sin( MAC::pi() * xC )
                                               * MAC::sin( MAC::pi() * yC )
                                               * MAC::sin( MAC::pi() * zC );
      }
   }

   return (bodyterm);
}

//---------------------------------------------------------------------------
void
DS_HeatTransfer:: assemble_DS_un_at_rhs ( FV_TimeIterator const* t_it
													 , double const& gamma)
//---------------------------------------------------------------------------
{
  // Assemble the diffusive components of velocity once in each iteration
  assemble_temperature_diffusion_terms();

  size_t_vector min_unknown_index(3,0);
  size_t_vector max_unknown_index(3,0);

  size_t_array2D* void_frac = (is_solids) ?
						  allrigidbodies->get_void_fraction_on_grid(TF) : 0;

  vector<doubleVector*> T_diffusion = GLOBAL_EQ->get_temperature_diffusion();

  for (size_t comp=0;comp<nb_comps;comp++) {
     // Get local min and max indices
     for (size_t l=0;l<dim;++l) {
        min_unknown_index(l) =
		  						TF->get_min_index_unknown_handled_by_proc( comp, l ) ;
        max_unknown_index(l) =
		  						TF->get_max_index_unknown_handled_by_proc( comp, l ) ;
     }

     for (size_t i = min_unknown_index(0); i <= max_unknown_index(0); ++i) {
        double dxC = TF->get_cell_size( i, comp, 0 ) ;
        double xC = TF->get_DOF_coordinate( i, comp, 0 ) ;
        for (size_t j = min_unknown_index(1); j <= max_unknown_index(1); ++j) {
           double dyC = TF->get_cell_size( j, comp, 1 ) ;
	   	  double yC = TF->get_DOF_coordinate( j, comp, 1 ) ;
			  for (size_t k = min_unknown_index(2);k <= max_unknown_index(2);++k) {
				  double dzC = (dim == 2) ? 1. : TF->get_cell_size( k, comp, 2 ) ;
				  double zC = (dim == 2) ? 0. : TF->get_DOF_coordinate( k, comp, 2 ) ;

				  size_t p = TF->DOF_local_number(i,j,k,comp);
              // Dxx for un
              double xvalue = T_diffusion[0]->operator()(p);
              // Dyy for un
              double yvalue = T_diffusion[1]->operator()(p);
				  // Dzz for un
				  double zvalue = (dim == 2) ? 0. : T_diffusion[2]->operator()(p);
		        // Bodyterm for rhs
		        double bodyterm = bodyterm_value(xC,yC,zC);
              // Advection term
              double adv_value = (is_NSwithHE) ? compute_adv_component(comp,i,j,k)
				  									 		  : 0.;

              if (is_solids) {
                 if (void_frac->operator()(p,0) != 0) {
                    adv_value = 0.;
                 }
              }

				  double rhs = gamma*(xvalue*dyC*dzC
					  					  + yvalue*dxC*dzC
										  + zvalue*dxC*dyC)
			  				  	 - adv_value
				  		 		 + (TF->DOF_value( i, j, k, comp, 1 )*dxC*dyC*dzC)
				  		 		  / t_it -> time_step();

              TF->set_DOF_value( i, j, k, comp, 0,
				  					rhs*(t_it -> time_step())/(dxC*dyC*dzC)
								 + gamma*bodyterm*(t_it -> time_step()));
           }
        }
     }
  }
}

//---------------------------------------------------------------------------
double DS_HeatTransfer:: divergence_of_U ( size_t const& comp
													  , size_t const& i
													  , size_t const& j
													  , size_t const& k
													  , size_t const& level)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_HeatTransfer:: divergence_of_U" ) ;

   FV_SHIFT_TRIPLET shift = TF->shift_staggeredToCentered() ;

   double xvalue = 0.,yvalue=0.,zvalue=0.,value=0.;

   double dxC = TF->get_cell_size( i, comp, 0 ) ;
   double dyC = TF->get_cell_size( j, comp, 1 ) ;
   double dzC = 0.;

   // du/dx
   double xhr = UF->get_DOF_coordinate( shift.i+i,0, 0 )
				  - UF->get_DOF_coordinate( shift.i+i-1, 0, 0 ) ;
   double xright = UF->DOF_value( shift.i+i, j, k, 0, level )
					  - UF->DOF_value( shift.i+i-1, j, k, 0, level ) ;

   xvalue = xright/xhr;

   // dv/dy
   double yhr = UF->get_DOF_coordinate( shift.j+j,1, 1 )
				  - UF->get_DOF_coordinate( shift.j+j-1, 1, 1 ) ;
   double yright = UF->DOF_value( i, shift.j+j, k, 1, level )
					  - UF->DOF_value( i, shift.j+j-1, k, 1, level ) ;

   yvalue = yright/yhr;

   if (dim == 3) {
      // dw/dz
      dzC = TF->get_cell_size( k, comp, 2 ) ;
      double zhr = UF->get_DOF_coordinate( shift.k+k,2, 2 )
					  - UF->get_DOF_coordinate( shift.k+k-1, 2, 2 ) ;
      double zright = UF->DOF_value( i, j, shift.k+k, 2, level )
						  - UF->DOF_value( i, j, shift.k+k-1, 2, level ) ;

      zvalue = zright/zhr;
   }

   value = (dim == 2) ? (xvalue + yvalue)*dxC*dyC :
                        (xvalue + yvalue + zvalue)*dxC*dyC*dzC ;

   return value;
}




//---------------------------------------------------------------------------
void DS_HeatTransfer:: assemble_temperature_diffusion_terms ( )
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_HeatTransfer:: assemble_temperature_diffusion_terms" ) ;

   size_t_vector min_unknown_index(3,0);
   size_t_vector max_unknown_index(3,0);

   vector<doubleVector*> T_diffusion = GLOBAL_EQ->get_temperature_diffusion();

   for (size_t comp = 0; comp < nb_comps; comp++) {
      // Get local min and max indices
      for (size_t l=0;l<dim;++l) {
         min_unknown_index(l) =
                        TF->get_min_index_unknown_handled_by_proc( comp, l ) ;
         max_unknown_index(l) =
                        TF->get_max_index_unknown_handled_by_proc( comp, l ) ;
      }

      for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i) {
         for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j) {
            for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k) {
               size_t p = TF->DOF_local_number(i,j,k,comp);
               // dxx of level3
               T_diffusion[0]->operator()(p) =
                                    compute_un_component(comp,i,j,k,0,3);
               // dyy of level1(2D) or level4(3D)
               if (dim == 2) {
                  T_diffusion[1]->operator()(p) =
                                    compute_un_component(comp,i,j,k,1,1);
               } else {
                  T_diffusion[1]->operator()(p) =
                                    compute_un_component(comp,i,j,k,1,4);
               }
               // dzz of level3
               if (dim == 3)
                  T_diffusion[2]->operator()(p) =
                                    compute_un_component(comp,i,j,k,2,1);

            }
         }
      }
   }
}




//---------------------------------------------------------------------------
double
DS_HeatTransfer:: compute_un_component ( size_t const& comp
													, size_t const& i
													, size_t const& j
													, size_t const& k
													, size_t const& dir
													, size_t const& level)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_HeatTransfer:: compute_un_component" ) ;

   double xhr,xhl,xright,xleft,yhr,yhl,yright,yleft;
   double zhr,zhl,zright,zleft, value=0.;

	size_t_array2D* intersect_vector = (is_solids) ?
									allrigidbodies->get_intersect_vector_on_grid(TF)
									: 0;
	doubleArray2D* intersect_distance = (is_solids) ?
									allrigidbodies->get_intersect_distance_on_grid(TF)
									: 0;
	doubleArray2D* intersect_fieldVal = (is_solids) ?
									allrigidbodies->get_intersect_fieldValue_on_grid(TF)
									: 0;

	size_t p = TF->DOF_local_number(i,j,k,comp);

   if (dir == 0) {
      xhr = TF->get_DOF_coordinate( i+1,comp, 0 )
			 - TF->get_DOF_coordinate( i, comp, 0 ) ;
      xhl = TF->get_DOF_coordinate( i, comp, 0 )
		    - TF->get_DOF_coordinate( i-1, comp, 0 ) ;
      xright = TF->DOF_value( i+1, j, k, comp, level )
				 - TF->DOF_value( i, j, k, comp, level ) ;
      xleft = TF->DOF_value( i, j, k, comp, level )
				- TF->DOF_value( i-1, j, k, comp, level ) ;

      if (is_solids) {
			size_t_array2D* void_frac = allrigidbodies
												->get_void_fraction_on_grid(TF);
         if (void_frac->operator()(p,0) == 0) {
            if (intersect_vector->operator()(p,2*dir+0) == 1) {
               xleft = TF->DOF_value( i, j, k, comp, level )
					      - intersect_fieldVal->operator()(p,2*dir+0);
               xhl = intersect_distance->operator()(p,2*dir+0);
            }
            if (intersect_vector->operator()(p,2*dir+1) == 1) {
               xright = intersect_fieldVal->operator()(p,2*dir+1)
							 - TF->DOF_value( i, j, k, comp, level );
               xhr = intersect_distance->operator()(p,2*dir+1);
            }
         } else {
            xright = 0.; xleft = 0.;
         }
      }

      //xvalue = xright/xhr - xleft/xhl;
      if (TF->DOF_in_domain( (int)i-1, (int)j, (int)k, comp)
		 && TF->DOF_in_domain( (int)i+1, (int)j, (int)k, comp))
         value = xright/xhr - xleft/xhl;
      else if (TF->DOF_in_domain( (int)i-1, (int)j, (int)k, comp))
         value = - xleft/xhl;
      else
         value = xright/xhr;
   } else if (dir == 1) {
      yhr = TF->get_DOF_coordinate( j+1,comp, 1 )
		    - TF->get_DOF_coordinate( j, comp, 1 ) ;
      yhl = TF->get_DOF_coordinate( j, comp, 1 )
			 - TF->get_DOF_coordinate( j-1, comp, 1 ) ;
      yright = TF->DOF_value( i, j+1, k, comp, level )
				 - TF->DOF_value( i, j, k, comp, level ) ;
      yleft = TF->DOF_value( i, j, k, comp, level )
				- TF->DOF_value( i, j-1, k, comp, level ) ;

      if (is_solids) {
			size_t_array2D* void_frac = allrigidbodies->get_void_fraction_on_grid(TF);
         if (void_frac->operator()(p,0) == 0) {
            if (intersect_vector->operator()(p,2*dir+0) == 1) {
               yleft = TF->DOF_value( i, j, k, comp, level )
					      - intersect_fieldVal->operator()(p,2*dir+0);
               yhl = intersect_distance->operator()(p,2*dir+0);
            }
            if (intersect_vector->operator()(p,2*dir+1) == 1) {
               yright = intersect_fieldVal->operator()(p,2*dir+1)
					       - TF->DOF_value( i, j, k, comp, level );
               yhr = intersect_distance->operator()(p,2*dir+1);
            }
         } else {
            yleft = 0.; yright = 0.;
         }
      }

      //yvalue = yright/yhr - yleft/yhl;
      if (TF->DOF_in_domain((int)i, (int)j-1, (int)k, comp)
		 && TF->DOF_in_domain((int)i, (int)j+1, (int)k, comp))
         value = yright/yhr - yleft/yhl;
      else if(TF->DOF_in_domain((int)i, (int)j-1, (int)k, comp))
         value = - yleft/yhl;
      else
         value = yright/yhr;
   } else if (dir == 2) {
      zhr = TF->get_DOF_coordinate( k+1,comp, 2 )
		    - TF->get_DOF_coordinate( k, comp, 2 ) ;
      zhl = TF->get_DOF_coordinate( k, comp, 2 )
		    - TF->get_DOF_coordinate( k-1, comp, 2 ) ;
      zright = TF->DOF_value( i, j, k+1, comp, level )
				 - TF->DOF_value( i, j, k, comp, level ) ;
      zleft = TF->DOF_value( i, j, k, comp, level )
			   - TF->DOF_value( i, j, k-1, comp, level ) ;

      if (is_solids) {
         size_t_array2D* void_frac = allrigidbodies->get_void_fraction_on_grid(TF);
         if (void_frac->operator()(p,0) == 0) {
            if (intersect_vector->operator()(p,2*dir+0) == 1) {
               zleft = TF->DOF_value( i, j, k, comp, level )
					 		- intersect_fieldVal->operator()(p,2*dir+0);
               zhl = intersect_distance->operator()(p,2*dir+0);
            }
            if (intersect_vector->operator()(p,2*dir+1) == 1) {
               zright = intersect_fieldVal->operator()(p,2*dir+1)
							 - TF->DOF_value( i, j, k, comp, level );
               zhr = intersect_distance->operator()(p,2*dir+1);
            }
         } else {
            zleft = 0.; zright = 0.;
         }
      }

      //zvalue = zright/zhr - zleft/zhl;
      if (TF->DOF_in_domain((int)i, (int)j, (int)k-1, comp)
		 && TF->DOF_in_domain((int)i, (int)j, (int)k+1, comp))
         value = zright/zhr - zleft/zhl;
      else if (TF->DOF_in_domain((int)i, (int)j, (int)k-1, comp))
         value = - zleft/zhl;
      else
         value = zright/zhr;
   }

   return(value);

}

//---------------------------------------------------------------------------
double
DS_HeatTransfer:: compute_adv_component ( size_t const& comp
													 , size_t const& i
													 , size_t const& j
													 , size_t const& k)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_HeatTransfer:: compute_adv_component" ) ;
   double ugradu = 0., value = 0.;

   if ( AdvectionScheme == "TVD" ) {
      ugradu = assemble_advection_TVD(UF,1,1.,i,j,k,1)
				 - TF->DOF_value(i,j,k,comp,1)
				 * divergence_of_U(comp,i,j,k,1);
   } else if ( AdvectionScheme == "Upwind" ) {
      ugradu = assemble_advection_Upwind(UF,1,1.,i,j,k,1)
				 - TF->DOF_value(i,j,k,comp,1)
				 * divergence_of_U(comp,i,j,k,1);
   } else if ( AdvectionScheme == "Centered" ) {
      ugradu = assemble_advection_Centered(UF,1,1.,i,j,k,1)
				 - TF->DOF_value(i,j,k,comp,1)
				 * divergence_of_U(comp,i,j,k,1);
   }

   if ( AdvectionTimeAccuracy == 1 ) {
      value = ugradu;
   } else {
      value = 1.5*ugradu - 0.5*TF->DOF_value(i,j,k,comp,2);
      TF->set_DOF_value(i,j,k,comp,2,ugradu);
   }

   return(value);
}




//---------------------------------------------------------------------------
void
DS_HeatTransfer:: calculate_row_indexes ( )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_HeatTransfer:: calculate_row_indexes" ) ;

   for (size_t comp = 0; comp < nb_comps; comp++) {
      // Get local min and max indices
      size_t_vector min_unknown_index(3,0);
      size_t_vector max_unknown_index(3,0);

      for (size_t l = 0; l < dim; ++l) {
         min_unknown_index(l) =
                           TF->get_min_index_unknown_handled_by_proc( comp, l ) ;
         max_unknown_index(l) =
                           TF->get_max_index_unknown_handled_by_proc( comp, l ) ;
      }

      for (size_t dir = 0; dir < dim; dir++) {
         size_t_array2D* row_index = GLOBAL_EQ->get_row_indexes(0,dir,comp);
         switch (dir) {
            case 0:
             for (size_t j=min_unknown_index(1); j<=max_unknown_index(1); j++) {
              for (size_t k=min_unknown_index(2); k<=max_unknown_index(2); k++) {
                  size_t p = (j - min_unknown_index(1))
                           + (1 + max_unknown_index(1) - min_unknown_index(1))
                           * (k - min_unknown_index(2));
                  row_index->operator()(j,k) = p;
              }
             }
             break;
            case 1:
             for (size_t i=min_unknown_index(0); i<=max_unknown_index(0); i++) {
              for (size_t k=min_unknown_index(2); k<=max_unknown_index(2); k++) {
                  size_t p = (i - min_unknown_index(0))
                           + (1 + max_unknown_index(0) - min_unknown_index(0))
                           * (k - min_unknown_index(2));
						row_index->operator()(i,k) = p;
              }
             }
             break;
            case 2:
             for (size_t i=min_unknown_index(0); i<=max_unknown_index(0); i++) {
              for (size_t j=min_unknown_index(1); j<=max_unknown_index(1); j++) {
                  size_t p = (i - min_unknown_index(0))
                           + (1 + max_unknown_index(0) - min_unknown_index(0))
                           * (j - min_unknown_index(1));
						row_index->operator()(i,j) = p;
              }
             }
             break;
         }
      }
   }
}




//----------------------------------------------------------------------
double DS_HeatTransfer:: compute_DS_temperature_change( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "DS_HeatTransfer:: compute_DS_temperature_change" ) ;

	size_t_vector min_unknown_index(3,0);
	size_t_vector max_unknown_index(3,0);

	size_t comp = 0;
	double sum_sq_U=0.,sum_sq_dU=0.;

	for (size_t l = 0; l < dim; ++l) {
		min_unknown_index(l) =
						 TF->get_min_index_unknown_handled_by_proc( comp, l ) ;
		max_unknown_index(l) =
						 TF->get_max_index_unknown_handled_by_proc( comp, l ) ;
	}

	for (size_t i = min_unknown_index(0); i <= max_unknown_index(0); ++i) {
		for (size_t j = min_unknown_index(1); j <= max_unknown_index(1); ++j) {
			for (size_t k = min_unknown_index(2); k <= max_unknown_index(2); ++k) {
				sum_sq_U += pow(TF->DOF_value( i, j, k, comp, 0 ),2.);
			   sum_sq_dU += pow(TF->DOF_value( i, j, k, comp, 0 )
								   - TF->DOF_value( i, j, k, comp, 1 ),2.);

			}
		}
	}

	sum_sq_U = macCOMM->sum(sum_sq_U);
	sum_sq_dU = macCOMM->sum(sum_sq_dU);

   return ( MAC::sqrt(sum_sq_dU/sum_sq_U) ) ;
}




//---------------------------------------------------------------------------
void
DS_HeatTransfer:: assemble_temperature_matrix (FV_DiscreteField const* FF
															, FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_HeatTransfer:: assemble_temperature_matrix" ) ;

	double gamma = (1.0/2.0)*(thermal_conductivity/rho/heat_capacity);
	TDMatrix* A = GLOBAL_EQ-> get_A();

   size_t_vector min_unknown_index(3,0);
   size_t_vector max_unknown_index(3,0);

	for (size_t comp = 0; comp < nb_comps; comp++) {
		// Get local min and max indices
		for (size_t l=0;l<dim;++l) {
			min_unknown_index(l) =
							FF->get_min_index_unknown_handled_by_proc( comp, l ) ;
			max_unknown_index(l) =
							FF->get_max_index_unknown_handled_by_proc( comp, l ) ;
	   }

		for (size_t dir = 0; dir < dim; dir++) {

			bool r_bound = false;
         bool l_bound = false;
         // All the proc will have open right bound,
  		 	// except last proc for non periodic systems
         if ((is_iperiodic[dir] != 1)
  		    && (rank_in_i[dir] == nb_ranks_comm_i[dir]-1))
			   r_bound = true;
         // All the proc will have open left bound,
  		   // except first proc for non periodic systems
         if ((is_iperiodic[dir] != 1)
			 && (rank_in_i[dir] == 0))
			 	l_bound = true;

			size_t dir_j = (dir == 0) ? 1 : 0;
			size_t dir_k = (dir == 2) ? 1 : 2;

			size_t local_min_k = (dim == 2) ? 0 : min_unknown_index(dir_k);
			size_t local_max_k = (dim == 2) ? 0 : max_unknown_index(dir_k);

			size_t_array2D* row_index = GLOBAL_EQ->get_row_indexes(0,dir,comp);

			for (size_t j = min_unknown_index(dir_j);
						  j <= max_unknown_index(dir_j);++j) {
				for (size_t k = local_min_k; k <= local_max_k; ++k) {
					size_t r_index = row_index->operator()(j,k);
					// Perform assembling
   				double Aee_diagcoef=0.;
   				for (size_t m = 0, i = min_unknown_index(dir);
											 i <= max_unknown_index(dir); ++i, ++m) {
    					double xC = FF->get_DOF_coordinate( i, comp, dir ) ;
       				double xR = FF->get_DOF_coordinate( i+1, comp, dir ) ;
       				double xL = FF->get_DOF_coordinate( i-1, comp, dir ) ;

       				double dxr = xR - xC;
       				double dxl = xC - xL;

       				double right = -gamma/(dxr);
       				double left = -gamma/(dxl);

						double center = -(right+left);

				      if (is_solids) {
				         size_t p=0;
				         if (dir == 0) {
				            p = FF->DOF_local_number(i,j,k,comp);
				         } else if (dir == 1) {
				            p = FF->DOF_local_number(j,i,k,comp);
				         } else if (dir == 2) {
				            p = FF->DOF_local_number(j,k,i,comp);
				         }

							size_t_array2D* intersect_vector =
				                   allrigidbodies->get_intersect_vector_on_grid(FF);
				         doubleArray2D* intersect_distance =
				                 allrigidbodies->get_intersect_distance_on_grid(FF);
				         size_t_array2D* void_frac =
				                      allrigidbodies->get_void_fraction_on_grid(FF);

				         // if left node is inside the solid particle
				         if (intersect_vector->operator()(p,2*dir+0) == 1) {
				            left = -gamma/intersect_distance->operator()(p,2*dir+0);
				         }
				         // if right node is inside the solid particle
				         if (intersect_vector->operator()(p,2*dir+1) == 1) {
				            right = -gamma/intersect_distance->operator()(p,2*dir+1);
				         }
				         // if center node is inside the solid particle
				         if (void_frac->operator()(p,0) != 0) {
				            left = 0.;
				            right = 0.;
				         }

          				center = -(right+left);

				         if (intersect_vector->operator()(p,2*dir+0) == 1)
								left = 0.;
				         if (intersect_vector->operator()(p,2*dir+1) == 1)
								right = 0.;
				      }

			         // add unsteady term
			         double value = center;
       				double unsteady_term =
									(FF->get_cell_size(i,comp,dir))/(t_it->time_step());

				      // Condition for handling the pressure neumann conditions at wall
				      if (i==min_unknown_index(dir) && l_bound) {
							int ii = (int)((dir == 0) ? i-1 : min_unknown_index(0));
							int jj = (int)((dir == 1) ? i-1 : min_unknown_index(1));
							int kk = (int)((dir == 2) ? i-1 : min_unknown_index(2));

			 				if (FF->DOF_in_domain(ii,jj,kk,comp)
			  				 && FF->DOF_has_imposed_Dirichlet_value
			  									((size_t)ii,(size_t)jj,(size_t)kk,comp)) {
             				// For Dirichlet boundary condition
             				value = center;
          				} else {
          					// For Neumann homogeneous boundary condition
             				value = -right;
          				}
       				} else if (i==max_unknown_index(dir) && r_bound) {
							int ii = (int)((dir == 0) ? i+1 : max_unknown_index(0));
							int jj = (int)((dir == 1) ? i+1 : max_unknown_index(1));
							int kk = (int)((dir == 2) ? i+1 : max_unknown_index(2));

          				if (FF->DOF_in_domain(ii,jj,kk,comp)
			  				 && FF->DOF_has_imposed_Dirichlet_value
		  										((size_t)ii,(size_t)jj,(size_t)kk,comp)) {
          					// For Dirichlet boundary condition
             				value = center;
          				} else {
             				// For Neumann homogeneous boundary condition
             				value = -left;
          				}
       				}

       				value = value + unsteady_term;

				      // Set Aie, Aei and Aee
				      if ((!l_bound) && (i == min_unknown_index(dir))) {
			            // Periodic boundary condition at minimum unknown index
			            // First proc has non zero value in Aie,Aei for first & last index
	  			  		   if (rank_in_i[dir] == 0) {
             				A[dir].ie[comp][r_index]->set_item(m,nb_ranks_comm_i[dir]-1,left);
      	     				A[dir].ei[comp][r_index]->set_item(nb_ranks_comm_i[dir]-1,m,left);
	  						} else {
             				A[dir].ie[comp][r_index]->set_item(m,rank_in_i[dir]-1,left);
	     						A[dir].ei[comp][r_index]->set_item(rank_in_i[dir]-1,m,left);
	  						}
       				}

       				if ((!r_bound) && (i == max_unknown_index(dir))) {
          				// Periodic boundary condition at maximum unknown index
          				// For last index, Aee comes from this proc as it
							// is interface unknown wrt this proc
          				A[dir].ie[comp][r_index]->set_item(m-1,rank_in_i[dir],left);
          				Aee_diagcoef = value;
          				A[dir].ei[comp][r_index]->set_item(rank_in_i[dir],m-1,left);
       				}

				      // Set Aii_sub_diagonal
				      if ((rank_in_i[dir] == nb_ranks_comm_i[dir]-1)
						 && (is_iperiodic[dir] != 1)) {
				         if (i > min_unknown_index(dir))
								A[dir].ii_sub[comp][r_index]->set_item(m-1,left);
				       } else {
				          if (i<max_unknown_index(dir))
				             if (i>min_unknown_index(dir))
				                A[dir].ii_sub[comp][r_index]->set_item(m-1,left);
				       }

       				 // Set Aii_super_diagonal
       				 if ((rank_in_i[dir] == nb_ranks_comm_i[dir]-1)
						  && (is_iperiodic[dir] != 1)) {
          				 if (i < max_unknown_index(dir))
							 	 A[dir].ii_super[comp][r_index]->set_item(m,right);
       				 } else {
          				 if (i < max_unknown_index(dir)-1)
             				 A[dir].ii_super[comp][r_index]->set_item(m,right);
       				 }

			        	 // Set Aii_main_diagonal
			       	 if ((rank_in_i[dir] == nb_ranks_comm_i[dir]-1)
						  && (is_iperiodic[dir] != 1)) {
          				 A[dir].ii_main[comp][r_index]->set_item(m,value);
       				 } else {
       					 if (i<max_unknown_index(dir))
          					 A[dir].ii_main[comp][r_index]->set_item(m,value);
       				 }
   		 		 }  // End of for loop

   		 		 GLOBAL_EQ->pre_thomas_treatment(comp,dir,A,r_index);

					 // Storing Aee for MPI communication
 					 ProdMatrix* Ap = GLOBAL_EQ->get_Ap();
 					 double* local_coeff = data_for_S[dir].send[comp][0];
 					 size_t nbrow = Ap[dir].ei_ii_ie[comp]->nb_rows();
 					 size_t ii = (nbrow*nbrow + 1)*r_index;
 					 local_coeff[ii] = Aee_diagcoef;
				 }
			 }
		 }
	 }
}




//---------------------------------------------------------------------------
void
DS_HeatTransfer:: assemble_schur_matrix ( FV_DiscreteField const* FF )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_HeatTransfer:: assemble_schur_matrix" ) ;

   TDMatrix* A = GLOBAL_EQ-> get_A();
	ProdMatrix* Ap = GLOBAL_EQ->get_Ap();

	size_t_vector min_unknown_index(dim,0);
	size_t_vector max_unknown_index(dim,0);

	for (size_t comp = 0; comp < nb_comps; comp++) {
		// Get local min and max indices
		for (size_t l=0;l<dim;++l) {
			min_unknown_index(l) =
								FF->get_min_index_unknown_handled_by_proc( comp, l ) ;
			max_unknown_index(l) =
								FF->get_max_index_unknown_handled_by_proc( comp, l ) ;
		}

		for (size_t dir = 0; dir < dim; dir++) {
			size_t dir_j = (dir == 0) ? 1 : 0;
			size_t dir_k = (dir == 2) ? 1 : 2;

			size_t_array2D* row_index = GLOBAL_EQ->get_row_indexes(0,dir,comp);

			size_t local_min_k = (dim == 2) ? 0 : min_unknown_index(dir_k);
			size_t local_max_k = (dim == 2) ? 0 : max_unknown_index(dir_k);

   		if (nb_ranks_comm_i[dir] > 1) {
				size_t nbrow = Ap[dir].ei_ii_ie[comp]->nb_rows();
				double* local_packet = data_for_S[dir].send[comp][0];

				// Calculating the product matrix and creating the container for MPI
	         for (size_t j = min_unknown_index(dir_j);
	                    j <= max_unknown_index(dir_j);++j) {
	            for (size_t k = local_min_k; k <= local_max_k; ++k) {
						size_t r_index = row_index->operator()(j,k);

      				GLOBAL_EQ->compute_product_matrix(A,Ap,comp,dir,r_index);

						// Create the data container to send and in the master as well
						for (size_t kk = 0; kk < nbrow; kk++) {
							for (size_t jj = 0; jj < nbrow; jj++) {
								size_t ii = (nbrow*nbrow + 1)*r_index
											 + nbrow*kk + jj + 1;
							   local_packet[ii] = Ap[dir].ei_ii_ie[comp]->item(kk,jj);
							}
						}
					}
				}

				if (rank_in_i[dir] != 0 ) {
					// Send the packed data to master
					MPI_Send( data_for_S[dir].send[comp][0],
						 (int) data_for_S[dir].size[comp],
							MPI_DOUBLE, 0, 0, DS_Comm_i[dir] ) ;
				}

				// Assemble the global product matrix by adding contribution from
				// all procs
      		if ( rank_in_i[dir] == 0 ) {
					for (size_t i = 1; i < (size_t)nb_ranks_comm_i[dir]; ++i) {
						// Recieve the data packet sent by master
						static MPI_Status status;
						MPI_Recv( data_for_S[dir].receive[comp][0],
							(int) data_for_S[dir].size[comp],
							MPI_DOUBLE, (int) i, 0, DS_Comm_i[dir], &status ) ;

		         	for (size_t j = min_unknown_index(dir_j);
		                    		j <= max_unknown_index(dir_j);++j) {
		            	for (size_t k = local_min_k; k <= local_max_k; ++k) {
								size_t r_index = row_index->operator()(j,k);
								size_t p = (nbrow*nbrow + 1)*r_index;
								double ee_proc0 = local_packet[p];
         					A[dir].ee[comp][r_index]->set_item(0,0,ee_proc0);

								for (size_t kk = 0;kk < nbrow; kk++) {
									for (size_t jj = 0;jj < nbrow; jj++) {
										size_t ii = (nbrow*nbrow + 1)*r_index
													 + nbrow*kk + jj + 1;
										double value =
												data_for_S[dir].receive[comp][0][ii];
										local_packet[ii] += value;
									}
								}

								size_t ii = (nbrow*nbrow + 1)*r_index;
								double value_Aee =
												data_for_S[dir].receive[comp][0][ii];

								// Assemble the global Aee matrix
								// Only for (nb_proc-1) in case no PBC
								// Otherwise for (nb_proc)
								if (!is_iperiodic[dir]) {
									if (i<(size_t)nb_ranks_comm_i[dir]-1) {
										A[dir].ee[comp][r_index]->set_item(i,i,value_Aee);
									}
								} else {
									A[dir].ee[comp][r_index]->set_item(i,i,value_Aee);
								}
							}
						}
					}
				}

      		// Assemble the schlur complement in the master proc
				if (rank_in_i[dir] == 0) {
					TDMatrix* Schur = GLOBAL_EQ-> get_Schur();
					for (size_t j = min_unknown_index(dir_j);
		                    j <= max_unknown_index(dir_j);++j) {
		            for (size_t k = local_min_k; k <= local_max_k; ++k) {
							size_t r_index = row_index->operator()(j,k);
							size_t schur_size = Schur[dir].ii_main[comp][r_index]
													  ->nb_rows();
							for (int p = 0; p < (int)schur_size; p++) {
								size_t ii = (nbrow*nbrow + 1)*r_index
											 + nbrow*p + p + 1;
								Schur[dir].ii_main[comp][r_index]
										->set_item(p,A[dir].ee[comp][r_index]->item(p,p)
											     								-local_packet[ii]);
								if (p < (int)schur_size-1)
									Schur[dir].ii_super[comp][r_index]
										->set_item(p,-local_packet[ii+1]);
								if (p > 0)
									Schur[dir].ii_sub[comp][r_index]
										->set_item(p-1,-local_packet[ii-1]);
								// In case of periodic and multi-processor,
								// there will be a variant of Tridiagonal matrix
								// instead of normal format
								if (is_iperiodic[dir] == 1) {
									ii = (nbrow*nbrow + 1)*r_index
										+ nbrow*p + schur_size + 1;
									Schur[dir].ie[comp][r_index]->set_item(p,0,
																				-local_packet[ii]);
									ii = (nbrow*nbrow + 1)*r_index
										+ nbrow*schur_size + p + 1;
									Schur[dir].ei[comp][r_index]->set_item(0,p,
																				-local_packet[ii]);
								}
							}
							// Pre-thomas treatment on Schur complement
							GLOBAL_EQ->pre_thomas_treatment(comp,dir,Schur,r_index);

							// In case of periodic and multi-processor, there will be
							// a variant of Tridiagonal matrix instead of normal format
							// So, Schur complement of Schur complement is calculated
							if (is_iperiodic[dir] == 1) {
								size_t ii = ((size_t)pow(nbrow,2) + 1)*r_index
											 + nbrow*schur_size
											 + schur_size + 1;
								Schur[dir].ee[comp][r_index]->set_item(0,0,
										A[dir].ee[comp][r_index]->item(schur_size,schur_size)
									 										- local_packet[ii]);

								ProdMatrix* SchurP = GLOBAL_EQ->get_SchurP();
								GLOBAL_EQ->compute_product_matrix_interior(Schur
															 ,SchurP,comp,0,dir,r_index);

								TDMatrix* DoubleSchur = GLOBAL_EQ->get_DoubleSchur();
								DoubleSchur[dir].ii_main[comp][r_index]
								->set_item(0,Schur[dir].ee[comp][r_index]->item(0,0)
								-SchurP[dir].ei_ii_ie[comp]->item(0,0));
							}
						}
					}
				}
			// Condition for single processor in any
			// direction with periodic boundary conditions
			} else if (is_iperiodic[dir] == 1) {
				for (size_t j = min_unknown_index(dir_j);
								j <= max_unknown_index(dir_j);++j) {
					for (size_t k = local_min_k; k <= local_max_k; ++k) {
						size_t r_index = row_index->operator()(j,k);

						GLOBAL_EQ->compute_product_matrix(A,Ap,comp,dir,r_index);

						LA_SeqMatrix* product_matrix = Ap[dir].ei_ii_ie[comp];

						A[dir].ee[comp][r_index]->set_item(0,0
								,data_for_S[dir].send[comp][0][0]);

						TDMatrix* Schur = GLOBAL_EQ-> get_Schur();
						size_t nb_row = Schur[dir].ii_main[comp][r_index]->nb_rows();
						for (int p = 0; p < (int)nb_row; p++) {
							Schur[dir].ii_main[comp][r_index]
										->set_item(p,A[dir].ee[comp][r_index]->item(p,p)
														-product_matrix->item(p,p));
							if (p < (int)nb_row-1)
								Schur[dir].ii_super[comp][r_index]
										->set_item(p,-product_matrix->item(p,p+1));
							if (p > 0)
							Schur[dir].ii_sub[comp][r_index]
							->set_item(p-1,-product_matrix->item(p,p-1));
						}
						GLOBAL_EQ->pre_thomas_treatment(comp,dir,Schur,r_index);
					}
				}
			}
      }
   }
}

//---------------------------------------------------------------------------
void
DS_HeatTransfer:: assemble_temperature_and_schur ( FV_TimeIterator const* t_it)
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_HeatTransfer:: assemble_temperature_and_schur" ) ;
	// Assemble temperature matrix
	assemble_temperature_matrix (TF,t_it);
	// Calculate and assemble Schur complement
	assemble_schur_matrix (TF);
}

//---------------------------------------------------------------------------
double DS_HeatTransfer:: assemble_local_rhs ( size_t const& j
														  , size_t const& k
														  , double const& gamma
														  , FV_TimeIterator const* t_it
														  , size_t const& comp
														  , size_t const& dir)
//---------------------------------------------------------------------------
{
   // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l) {
     min_unknown_index(l) = TF->get_min_index_unknown_handled_by_proc(comp,l) ;
     max_unknown_index(l) = TF->get_max_index_unknown_handled_by_proc(comp,l) ;
   }

   // Compute VEC_rhs_x = rhs in x
   double fe = 0.;

	// Since, this function is used in all directions;
   // ii, jj, and kk are used to convert the passed
   // arguments corresponding to correct direction
   size_t ii=0,jj=0,kk=0;
	size_t level = 0;

   // Vector for fi
   LocalVector* VEC = GLOBAL_EQ->get_VEC();

	vector<doubleVector*> T_diffusion = GLOBAL_EQ->get_temperature_diffusion();
	size_t_array2D* intersect_vector = (is_solids) ?
                           allrigidbodies->get_intersect_vector_on_grid(TF)
									: 0;
   doubleArray2D* intersect_distance = (is_solids) ?
                           allrigidbodies->get_intersect_distance_on_grid(TF)
									: 0;
   doubleArray2D* intersect_fieldVal = (is_solids) ?
                           allrigidbodies->get_intersect_fieldValue_on_grid(TF)
									: 0;

   for (size_t i = min_unknown_index(dir);i <= max_unknown_index(dir); ++i) {
		if (dir == 0) {
         ii = i; jj = j; kk = k; level = 0;
      } else if (dir == 1) {
         ii = j; jj = i; kk = k; level = 3;
      } else if (dir == 2) {
         ii = j; jj = k; kk = i; level = 4;
      }
      size_t pos = i - min_unknown_index(dir);
		size_t p = TF->DOF_local_number(ii,jj,kk,comp);

		double value = T_diffusion[dir]->operator()(p);

		if (is_solids) {
			if (intersect_vector->operator()(p,2*dir+0) == 1) {
				value = value - intersect_fieldVal->operator()(p,2*dir+0)
								  / intersect_distance->operator()(p,2*dir+0);
			}
			if (intersect_vector->operator()(p,2*dir+1) == 1) {
				value = value - intersect_fieldVal->operator()(p,2*dir+1)
								  / intersect_distance->operator()(p,2*dir+1);
			}
		}

      double dC = TF->get_cell_size(i,comp,dir) ;

      double temp_val = (TF->DOF_value(ii,jj,kk,comp,level)*dC)
							 / (t_it->time_step()) - gamma*value;

      if (is_iperiodic[dir] == 0) {
         if (rank_in_i[dir] == nb_ranks_comm_i[dir]-1) {
            VEC[dir].local_T[comp]->set_item( pos,temp_val);
         } else {
            if (i == max_unknown_index(dir))
               fe = temp_val;
            else
               VEC[dir].local_T[comp]->set_item( pos,temp_val);
         }
      } else {
         if (i == max_unknown_index(dir))
            fe = temp_val;
         else
            VEC[dir].local_T[comp]->set_item( pos,temp_val);
      }

   }

   // Effect of boundary conditions in case of non-periodic direction
   int m = int(min_unknown_index(dir)) - 1;

   if (dir == 0) {
      ii = m; jj = j; kk = k;
   } else if (dir == 1) {
      ii = j; jj = m; kk = k;
   } else if (dir == 2) {
      ii = j; jj = k; kk = m;
   }

   if ( TF->DOF_in_domain((int)ii,(int)jj,(int)kk,comp))
      if ( TF->DOF_has_imposed_Dirichlet_value(ii,jj,kk,comp)) {
         double ai = 1/(TF->get_DOF_coordinate(m+1,comp,dir)
							 - TF->get_DOF_coordinate(m,comp,dir));
         double dirichlet_value = TF->DOF_value(ii,jj,kk,comp,1) ;
         VEC[dir].local_T[comp]->add_to_item( 0
														, + gamma * ai * dirichlet_value );
      }

   m = int(max_unknown_index(dir)) + 1;

   if (dir == 0) {
      ii = m; jj = j; kk = k;
   } else if (dir == 1) {
      ii = j; jj = m; kk = k;
   } else if (dir == 2) {
      ii = j; jj = k; kk = m;
   }

   if ( TF->DOF_in_domain((int)ii,(int)jj,(int)kk,comp))
      if ( TF->DOF_has_imposed_Dirichlet_value(ii,jj,kk,comp)) {
      double ai = 1/(TF->get_DOF_coordinate(m,comp,dir)
						 - TF->get_DOF_coordinate(m-1,comp,dir));
      double dirichlet_value = TF->DOF_value(ii,jj,kk,comp,1) ;
      VEC[dir].local_T[comp]->add_to_item( VEC[dir].local_T[comp]->nb_rows()-1
													, + gamma * ai * dirichlet_value );
   }

   return fe;
}

//---------------------------------------------------------------------------
void
DS_HeatTransfer:: compute_Aei_ui (struct TDMatrix* arr
										  , struct LocalVector* VEC
										  , size_t const& comp
										  , size_t const& dir
										  , size_t const& r_index)
//---------------------------------------------------------------------------
{
   // create a replica of local rhs vector in local solution vector
   for (size_t i=0;i<VEC[dir].local_T[comp]->nb_rows();i++){
      VEC[dir].local_solution_T[comp]->set_item(i
															,VEC[dir].local_T[comp]->item(i));
   }

   // Solve for ui locally and put it in local solution vector
   GLOBAL_EQ->mod_thomas_algorithm(arr, VEC[dir].local_solution_T[comp]
												  , comp, dir,r_index);

   for (size_t i=0;i<VEC[dir].T[comp]->nb_rows();i++){
          VEC[dir].T[comp]->set_item(i,0);
   }

   // Calculate Aei*ui in each proc locally and put it in T vector
   arr[dir].ei[comp][r_index]->multiply_vec_then_add(
																VEC[dir].local_solution_T[comp]
															  ,VEC[dir].T[comp]);

}

//---------------------------------------------------------------------------
void
DS_HeatTransfer:: data_packing ( double const& fe
										 , size_t const& comp
										 , size_t const& dir
										 , size_t const& p)
//---------------------------------------------------------------------------
{
   LocalVector* VEC = GLOBAL_EQ->get_VEC() ;

   double *packed_data = first_pass[dir].send[comp][0];

   if (rank_in_i[dir] == 0) {
      // Check if bc is periodic in x
      // If it is, we need to pack two elements apart from fe
      if(is_iperiodic[dir])
          packed_data[3*p+0] = VEC[dir].T[comp]->item(nb_ranks_comm_i[dir]-1);
      else
          packed_data[3*p+0] = 0;

      packed_data[3*p+1] = VEC[dir].T[comp]->item(rank_in_i[dir]);

   } else if (rank_in_i[dir] == nb_ranks_comm_i[dir]-1) {
      // Check if bc is periodic in x
      // If it is, we need to pack two elements apart from fe
      if(is_iperiodic[dir])
          packed_data[3*p+1] = VEC[dir].T[comp]->item(rank_in_i[dir]);
      else
          packed_data[3*p+1]=0;

      packed_data[3*p+0] = VEC[dir].T[comp]->item(rank_in_i[dir]-1);

   } else {
      packed_data[3*p+0] = VEC[dir].T[comp]->item(rank_in_i[dir]-1);
      packed_data[3*p+1] = VEC[dir].T[comp]->item(rank_in_i[dir]);
   }

   packed_data[3*p+2] = fe; // Send the fe values and 0 for last proc

}

//---------------------------------------------------------------------------
void
DS_HeatTransfer:: unpack_compute_ue_pack(size_t const& comp
													, size_t const& dir
													, size_t const& p)
//---------------------------------------------------------------------------
{
   LocalVector* VEC = GLOBAL_EQ->get_VEC() ;

   size_t nb_interface_unknowns = VEC[dir].T[comp]->nb_rows();

	VEC[dir].T[comp]->set(0.);
   VEC[dir].interface_T[comp]->set(0.);

   // If periodic in x, first proc contributes to last interface unknown
   if (is_iperiodic[dir])
      VEC[dir].T[comp]->set_item(nb_ranks_comm_i[dir]-1,
							first_pass[dir].send[comp][0][3*p]);
   VEC[dir].T[comp]->set_item(0,
							first_pass[dir].send[comp][0][3*p+1]);
   VEC[dir].interface_T[comp]->set_item(0,
							first_pass[dir].send[comp][0][3*p+2]);

   // Vec_temp might contain previous values
   for (size_t i = 1; i < (size_t)nb_ranks_comm_i[dir];i++) {
      if (i!=(size_t)nb_ranks_comm_i[dir]-1) {
         VEC[dir].T[comp]->add_to_item(i-1,first_pass[dir].receive[comp][i][3*p]);
         VEC[dir].T[comp]->add_to_item(i,first_pass[dir].receive[comp][i][3*p+1]);
         VEC[dir].interface_T[comp]->set_item(i,first_pass[dir].receive[comp][i][3*p+2]);  // Assemble the interface rhs fe
      } else {
         if (is_iperiodic[dir] ==0) {
            VEC[dir].T[comp]->add_to_item(i-1,first_pass[dir].receive[comp][i][3*p]);
         } else {
            VEC[dir].T[comp]->add_to_item(i-1,first_pass[dir].receive[comp][i][3*p]);
            // If periodic in x, last proc has an interface unknown
            VEC[dir].T[comp]->add_to_item(i,first_pass[dir].receive[comp][i][3*p+1]);
            VEC[dir].interface_T[comp]->set_item(i,first_pass[dir].receive[comp][i][3*p+2]);
         }
      }
   }

   for (size_t i=0;i<nb_interface_unknowns;i++) {
      VEC[dir].interface_T[comp]->set_item(i,VEC[dir].interface_T[comp]->item(i)-VEC[dir].T[comp]->item(i)); // Get fe - Aei*xi to solve for ue
   }

   // Solve for ue (interface unknowns) in the master proc
   DS_interface_unknown_solver(VEC[dir].interface_T[comp], comp, dir,p);

   // Pack the interface_rhs_x into the appropriate send_data
   for (size_t i = 1; i < (size_t)nb_ranks_comm_i[dir]; ++i) {
      if (i!=(size_t)nb_ranks_comm_i[dir]-1) {
         second_pass[dir].send[comp][i][2*p+0] = VEC[dir].interface_T[comp]->item(i-1);
         second_pass[dir].send[comp][i][2*p+1] = VEC[dir].interface_T[comp]->item(i);
      } else {
         second_pass[dir].send[comp][i][2*p+0] = VEC[dir].interface_T[comp]->item(i-1);
         if (is_iperiodic[dir])
            second_pass[dir].send[comp][i][2*p+1] = VEC[dir].interface_T[comp]->item(i);
         else
            second_pass[dir].send[comp][i][2*p+1] = 0;
      }
   }

}

//---------------------------------------------------------------------------
void
DS_HeatTransfer:: unpack_ue(size_t const& comp
								  , double * received_data
								  , size_t const& dir
								  , size_t const& p)
//---------------------------------------------------------------------------
{
   LocalVector* VEC = GLOBAL_EQ->get_VEC() ;

   if (rank_in_i[dir] != nb_ranks_comm_i[dir]-1) {
      VEC[dir].interface_T[comp]->set_item(rank_in_i[dir]-1,received_data[2*p]);
      VEC[dir].interface_T[comp]->set_item(rank_in_i[dir],received_data[2*p+1]);
   } else {
      if (is_iperiodic[dir] ==0) {
         VEC[dir].interface_T[comp]->set_item(rank_in_i[dir]-1,received_data[2*p]);
      } else {
         VEC[dir].interface_T[comp]->set_item(rank_in_i[dir]-1,received_data[2*p]);
         VEC[dir].interface_T[comp]->set_item(rank_in_i[dir],received_data[2*p+1]);
      }
   }
}

//---------------------------------------------------------------------------
void
DS_HeatTransfer:: solve_interface_unknowns ( double const& gamma
														 , FV_TimeIterator const* t_it
														 , size_t const& comp
														 , size_t const& dir
													 	 , size_t const& level )
//---------------------------------------------------------------------------
{
   // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l) {
      min_unknown_index(l) = TF->get_min_index_unknown_handled_by_proc( comp, l ) ;
      max_unknown_index(l) = TF->get_max_index_unknown_handled_by_proc( comp, l ) ;
   }

   TDMatrix* A = GLOBAL_EQ->get_A();
   LocalVector* VEC = GLOBAL_EQ->get_VEC() ;

	// Array declaration for sending data from master to all slaves
	size_t local_min_j = (dir == 0) ? min_unknown_index(1) : min_unknown_index(0);
	size_t local_max_j = (dir == 0) ? max_unknown_index(1) : max_unknown_index(0);
	size_t local_min_k = (dim == 2) ? 0 :
							  ((dir == 2) ? min_unknown_index(1) : min_unknown_index(2));
	size_t local_max_k = (dim == 2) ? 0 :
							  ((dir == 2) ? max_unknown_index(1) : max_unknown_index(2));

   size_t_array2D* row_index = GLOBAL_EQ->get_row_indexes(0,dir,comp);

   // Send and receive the data first pass
   if ( rank_in_i[dir] == 0 ) {
      // Receiving data from all the slave procs iff multi processors are used
      if (nb_ranks_comm_i[dir] != 1) {
         for (size_t i = 1; i < (size_t)nb_ranks_comm_i[dir]; ++i) {
            static MPI_Status status;
            MPI_Recv( first_pass[dir].receive[comp][i],
				    (int) first_pass[dir].size[comp],
				  MPI_DOUBLE, (int) i, 0, DS_Comm_i[dir], &status ) ;
         }
      }

		for (size_t j = local_min_j; j <= local_max_j; j++) {
         for (size_t k = local_min_k; k <= local_max_k; k++) {

            size_t p = row_index->operator()(j,k);

            unpack_compute_ue_pack(comp,dir,p);

  	         // Need to have the original rhs function
            // assembled for corrosponding j,k pair
            assemble_local_rhs(j,k,gamma,t_it,comp,dir);

            // Setup RHS = fi - Aie*xe for solving ui
            A[dir].ie[comp][p]->multiply_vec_then_add(VEC[dir].interface_T[comp]
                                             ,VEC[dir].local_T[comp],-1.0,1.0);

            // Solve ui and transfer solution into distributed vector
            GLOBAL_EQ->DS_HeatEquation_solver(j,k,min_unknown_index(dir)
                                                            ,comp,dir,p,level);
         }
      }
   } else {
      // Send the packed data to master
      MPI_Send( first_pass[dir].send[comp][0],
			 (int) first_pass[dir].size[comp],
			 MPI_DOUBLE, 0, 0, DS_Comm_i[dir] ) ;
   }

   // Send the data from master iff multi processor are used
   if (nb_ranks_comm_i[dir] != 1) {
      if ( rank_in_i[dir] == 0 ) {
         for (size_t i = 1;i < (size_t)nb_ranks_comm_i[dir]; ++i) {
            MPI_Send( second_pass[dir].send[comp][i],
					 (int) second_pass[dir].size[comp],
					 MPI_DOUBLE,(int) i, 0, DS_Comm_i[dir] ) ;
         }
      } else {
         // Create the container to receive the ue
         static MPI_Status status ;
         MPI_Recv( second_pass[dir].receive[comp][0],
				 (int) second_pass[dir].size[comp],
				 MPI_DOUBLE, 0, 0, DS_Comm_i[dir], &status ) ;

			// Solve the system of equations in each proc
			for (size_t j = local_min_j; j <= local_max_j; j++) {
				for (size_t k = local_min_k; k <= local_max_k; k++) {
					size_t p = row_index->operator()(j,k);

					unpack_ue(comp,second_pass[dir].receive[comp][0],dir,p);

					// Need to have the original rhs function
					// assembled for corrosponding j,k pair
					assemble_local_rhs(j,k,gamma,t_it,comp,dir);

					// Setup RHS = fi - Aie*xe for solving ui
					A[dir].ie[comp][p]
							  ->multiply_vec_then_add(VEC[dir].interface_T[comp]
															 ,VEC[dir].local_T[comp],-1.0,1.0);

					// Solve ui and transfer solution into distributed vector
					GLOBAL_EQ->DS_HeatEquation_solver(j,k,min_unknown_index(dir)
																				,comp,dir,p,level);
				}
			}
      }
   }
}

//---------------------------------------------------------------------------
void
DS_HeatTransfer:: Solve_i_in_jk ( FV_TimeIterator const* t_it
										  , double const& gamma
										  , size_t const& dir_i
										  , size_t const& dir_j
										  , size_t const& dir_k
									     , size_t const& level)
//---------------------------------------------------------------------------
{
  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);

  for (size_t comp=0;comp<nb_comps;comp++) {
     // Get local min and max indices
     for (size_t l=0;l<dim;++l) {
        min_unknown_index(l) = TF->get_min_index_unknown_handled_by_proc( comp, l ) ;
        max_unknown_index(l) = TF->get_max_index_unknown_handled_by_proc( comp, l ) ;
     }

     size_t local_min_k = 0;
     size_t local_max_k = 0;

     if (dim == 3) {
        local_min_k = min_unknown_index(dir_k);
        local_max_k = max_unknown_index(dir_k);
     }

     LocalVector* VEC = GLOBAL_EQ->get_VEC() ;
     TDMatrix* A = GLOBAL_EQ->get_A();
	  size_t_array2D* row_index = GLOBAL_EQ->get_row_indexes(0,dir_i,comp);

     // Solve in i
     if ((nb_ranks_comm_i[dir_i]>1)||(is_iperiodic[dir_i] == 1)) {
        for (size_t j=min_unknown_index(dir_j);j<=max_unknown_index(dir_j);++j) {
	   	  for (size_t k=local_min_k; k <= local_max_k; ++k) {
				  size_t r_index = row_index->operator()(j,k);
              // Assemble fi and return fe for each proc locally
              double fe = assemble_local_rhs(j,k,gamma,t_it,comp,dir_i);
              // Calculate Aei*ui in each proc locally
              compute_Aei_ui(A,VEC,comp,dir_i,r_index);
              // Pack Aei_ui and fe for sending it to master
              data_packing (fe,comp,dir_i,r_index);
	   	  }
        }
        solve_interface_unknowns (gamma,t_it,comp,dir_i,level);

     } else if (is_iperiodic[dir_i] == 0) {  // Serial mode with non-periodic condition
        for (size_t j=min_unknown_index(dir_j);j<=max_unknown_index(dir_j);++j) {
           for (size_t k=local_min_k; k <= local_max_k; ++k) {
				  size_t r_index = row_index->operator()(j,k);
              assemble_local_rhs(j,k,gamma,t_it,comp,dir_i);
              GLOBAL_EQ->DS_HeatEquation_solver(j,k,min_unknown_index(dir_i)
				  													,comp,dir_i,r_index,level);
           }
        }
     }
  }
}

//----------------------------------------------------------------------
void
DS_HeatTransfer::DS_interface_unknown_solver(LA_SeqVector* interface_rhs
														 , size_t const& comp
														 , size_t const& dir
														 , size_t const& r_index)
//----------------------------------------------------------------------
{
   MAC_LABEL( "DS_HeatTransfer:: DS_interface_unknown_solver" ) ;

   TDMatrix* Schur = GLOBAL_EQ-> get_Schur();

   // Condition for variant of Tridiagonal Schur complement
	// in Perioidic direction with multi-processor
   if ((is_iperiodic[dir] == 1) && (nb_ranks_comm_i[dir] != 1)) {
      LocalVector* Schur_VEC = GLOBAL_EQ->get_Schur_VEC();
      TDMatrix* DoubleSchur = GLOBAL_EQ-> get_DoubleSchur();

      // Transfer interface_rhs to Schur VEC (i.e. S_fi and S_fe)
      size_t nrows = Schur_VEC[dir].local_T[comp]->nb_rows();
      for (size_t i = 0; i < nrows; i++) {
          Schur_VEC[dir].local_T[comp]->set_item(i,interface_rhs->item(i));
      }
      Schur_VEC[dir].interface_T[comp]->set_item(0,interface_rhs->item(nrows));

      // Calculate Sei*(Sii)-1*S_fi
      compute_Aei_ui(Schur,Schur_VEC,comp,dir,r_index);

      // Calculate S_fe - Sei*(Sii)-1*S_fi
      Schur_VEC[dir].interface_T[comp]->set_item(0,
												Schur_VEC[dir].interface_T[comp]->item(0)
											 - Schur_VEC[dir].T[comp]->item(0));

      // Calculate S_ue, using Schur complement of Schur complement
      GLOBAL_EQ->mod_thomas_algorithm(DoubleSchur
												, Schur_VEC[dir].interface_T[comp]
												, comp
												, dir
												, r_index);

      // Calculate S_fi-Sie*S_ue
      Schur[dir].ie[comp][r_index]->multiply_vec_then_add(
													Schur_VEC[dir].interface_T[comp]
												 , Schur_VEC[dir].local_T[comp],-1.0,1.0);

      // Calculate S_ui
      GLOBAL_EQ->mod_thomas_algorithm(Schur
												, Schur_VEC[dir].local_T[comp]
												, comp
												, dir
												, r_index);

      // Transfer back the solution to interface_rhs
      for (size_t i = 0; i < nrows; i++) {
          interface_rhs->set_item(i,Schur_VEC[dir].local_T[comp]->item(i));
      }
      interface_rhs->set_item(nrows,Schur_VEC[dir].interface_T[comp]->item(0));
   } else {
      GLOBAL_EQ->mod_thomas_algorithm(Schur, interface_rhs, comp, dir,r_index);
   }

}




//---------------------------------------------------------------------------
void
DS_HeatTransfer:: ugradu_initialization (  )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_HeatTransfer:: ugradu_initialization" ) ;

  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);

  for (size_t comp=0;comp<nb_comps;comp++) {
     // Get local min and max indices
     for (size_t l=0;l<dim;++l) {
        min_unknown_index(l) =
		  					TF->get_min_index_unknown_handled_by_proc( comp, l ) ;
        max_unknown_index(l) =
		  					TF->get_max_index_unknown_handled_by_proc( comp, l ) ;
     }

     size_t local_min_k = 0;
     size_t local_max_k = 0;

     if (dim == 3) {
        local_min_k = min_unknown_index(2);
        local_max_k = max_unknown_index(2);
     }

     for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i) {
        for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j) {
           for (size_t k=local_min_k;k<=local_max_k;++k) {
              TF->set_DOF_value( i, j, k, comp, 2, 0.);
           }
        }
     }
  }
}

//---------------------------------------------------------------------------
void
DS_HeatTransfer:: nodes_temperature_initialization ( size_t const& level )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_HeatTransfer:: Solids_flux_correction" ) ;

  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);

  // Vector for solid presence
  size_t_array2D* void_frac = allrigidbodies->get_void_fraction_on_grid(TF);

  for (size_t comp=0;comp<nb_comps;comp++) {
     // Get local min and max indices
     for (size_t l=0;l<dim;++l) {
        if (is_iperiodic[l]) {
           min_unknown_index(l) =
			  				TF->get_min_index_unknown_handled_by_proc( comp, l ) - 1;
           max_unknown_index(l) =
			  				TF->get_max_index_unknown_handled_by_proc( comp, l ) + 1;
        } else {
           min_unknown_index(l) =
			  				TF->get_min_index_unknown_handled_by_proc( comp, l );
           max_unknown_index(l) =
			  				TF->get_max_index_unknown_handled_by_proc( comp, l );
        }
     }

     size_t local_min_k = 0;
     size_t local_max_k = 0;

     if (dim == 3) {
        local_min_k = min_unknown_index(2);
        local_max_k = max_unknown_index(2);
     }

     for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i) {
		  double xC = TF->get_DOF_coordinate( i, comp, 0 ) ;
        for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j) {
			  double yC = TF->get_DOF_coordinate( j, comp, 1 ) ;
           for (size_t k=local_min_k;k<=local_max_k;++k) {
				  double zC = (dim == 2) ? 0 : TF->get_DOF_coordinate( k, comp, 2 );
				  geomVector pt(xC,yC,zC);
              size_t p = TF->DOF_local_number(i,j,k,comp);
              if (void_frac->operator()(p,0) != 0) {
                 size_t par_id = void_frac->operator()(p,0) - 1;
					  geomVector Tpart = allrigidbodies->rigid_body_temperature(par_id,pt);
                 TF->set_DOF_value( i, j, k, comp, level,Tpart(comp));
              }
           }
        }
     }
  }
}




//---------------------------------------------------------------------------
void
DS_HeatTransfer:: HeatEquation_DirectionSplittingSolver ( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_HeatTransfer:: HeatEquation_DirectionSplittingSolver" ) ;

  double gamma= thermal_conductivity/rho/heat_capacity;

  TF->copy_DOFs_value( 0, 1 );

  // First Equation
  if ( my_rank == is_master ) SCT_set_start("Solver first step");
  assemble_DS_un_at_rhs (t_it,gamma);
  if ( my_rank == is_master ) SCT_get_elapsed_time("Solver first step");
  // Update gamma based for invidual direction
  gamma = (1.0/2.0)*(thermal_conductivity/rho/heat_capacity);


  if ( my_rank == is_master ) SCT_set_start("Solver x solution");
  // Solve x-direction(i.e. 0) in y(i.e. 1) and z(i.e. 2)
  Solve_i_in_jk (t_it,gamma,0,1,2,3);
  // Synchronize the temperature field
  TF->synchronize( 3 );
  if (is_solids) nodes_temperature_initialization(3);
  if ( my_rank == is_master ) SCT_get_elapsed_time("Solver x solution");

  size_t level = (dim == 2) ? 0 : 4 ;

  if ( my_rank == is_master ) SCT_set_start("Solver y solution");
  // Solve y-direction(i.e. 1) in x(i.e. 0) and z(i.e. 2)
  Solve_i_in_jk (t_it,gamma,1,0,2,level);
  // Synchronize the temperature field
  TF->synchronize( level );
  if (is_solids) nodes_temperature_initialization(level);
  if ( my_rank == is_master ) SCT_get_elapsed_time("Solver y solution");


  if (dim == 3) {
     if ( my_rank == is_master ) SCT_set_start("Solver z solution");
     // Solve z-direction(i.e. 2) in x(i.e. 0) and y(i.e. 1)
     Solve_i_in_jk (t_it,gamma,2,0,1,0);
     // Synchronize the temperature field
	  TF->synchronize( 0 );
     if (is_solids) nodes_temperature_initialization(0);
	  if ( my_rank == is_master ) SCT_get_elapsed_time("Solver z solution");
  }
}

//---------------------------------------------------------------------------
void
DS_HeatTransfer:: output_l2norm ( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL("DS_HeatTransfer:: output_l2norm" ) ;

   // Parameters
   size_t k=0;

   size_t_vector min_unknown_index(dim,0);
   size_t_vector max_unknown_index(dim,0);

   for (size_t comp=0;comp<nb_comps;comp++) {
      double computed_DS_field=0., computed_DS_L2 = 0.;

      // Get local min and max indices
      for (size_t l=0;l<dim;++l) {
         min_unknown_index(l) = TF->get_min_index_unknown_handled_by_proc( comp, l ) ;
         max_unknown_index(l) = TF->get_max_index_unknown_handled_by_proc( comp, l ) ;
      }

      for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i) {
         for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j) {

            if (dim == 2) {
               k=0;
               computed_DS_field = TF->DOF_value( i, j, k, comp, 0 ) ;

               if ( TF->DOF_is_unknown_handled_by_proc( i, j, k, comp ) ) {
                  computed_DS_L2 += computed_DS_field*computed_DS_field
                               * TF->get_cell_measure( i, j, k, comp ) ;
               }
            } else {
               for (k=min_unknown_index(2);k<=max_unknown_index(2);++k) {
                  computed_DS_field = TF->DOF_value( i, j, k, comp, 0 ) ;

                  if ( TF->DOF_is_unknown_handled_by_proc( i, j, k, comp ) ) {
                     computed_DS_L2 += computed_DS_field*computed_DS_field
                                  * TF->get_cell_measure( i, j, k, comp ) ;
                  }
               }
            }
         }
      }

      computed_DS_L2 = macCOMM->sum( computed_DS_L2 ) ;
      computed_DS_L2 = MAC::sqrt(computed_DS_L2);

      if (my_rank == is_master) {
			MAC::out() << "Norm L2 T = "
	                 << MAC::doubleToString( ios::scientific, 12, computed_DS_L2 )
	                 << endl;
      }
   }
}




//---------------------------------------------------------------------------
void
DS_HeatTransfer:: create_DS_subcommunicators ( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_HeatTransfer:: create_DS_subcommunicators" ) ;

   int color = 0, key = 0;
   //int const* number_of_subdomains_per_direction = TF->primary_grid()->get_domain_decomposition() ;
   int const* MPI_coordinates_world =
									TF->primary_grid()->get_MPI_coordinates() ;
   int const* MPI_number_of_coordinates =
									TF->primary_grid()->get_domain_decomposition() ;

   if (dim == 2) {
      // Assign color and key for splitting in x
      color = MPI_coordinates_world[1];
      key = MPI_coordinates_world[0];
      // Split by direction in x
      processor_splitting (color,key,0);

      // Assign color and key for splitting in y
      color = MPI_coordinates_world[0];
      key = MPI_coordinates_world[1];
      // Split by direction in y
      processor_splitting (color,key,1);
   } else {
      // Assign color and key for splitting in x
      color = MPI_coordinates_world[1]
				+ MPI_coordinates_world[2]*MPI_number_of_coordinates[1] ;
      key = MPI_coordinates_world[0];
      // Split by direction in x
      processor_splitting (color,key,0);

      // Assign color and key for splitting in y
      color = MPI_coordinates_world[2]
				+ MPI_coordinates_world[0]*MPI_number_of_coordinates[2];
      key = MPI_coordinates_world[1];
      // Split by direction in y
      processor_splitting (color,key,1);

      // Assign color and key for splitting in y
      color = MPI_coordinates_world[0]
				+ MPI_coordinates_world[1]*MPI_number_of_coordinates[0];;
      key = MPI_coordinates_world[2];

      // Split by direction in y
      processor_splitting (color,key,2);
   }
}

//---------------------------------------------------------------------------
void
DS_HeatTransfer:: processor_splitting ( int const& color
												  , int const& key
												  , size_t const& dir )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_HeatTransfer:: processor_splitting" ) ;

   MPI_Comm_split(MPI_COMM_WORLD, color, key, &DS_Comm_i[dir]);
   MPI_Comm_size( DS_Comm_i[dir], &nb_ranks_comm_i[dir] ) ;
   MPI_Comm_rank( DS_Comm_i[dir], &rank_in_i[dir] ) ;

}

//---------------------------------------------------------------------------
void
DS_HeatTransfer:: allocate_mpi_variables ( void )
//---------------------------------------------------------------------------
{

   for (size_t dir = 0; dir < dim; dir++) {
      first_pass[dir].size = new size_t [nb_comps];
      second_pass[dir].size = new size_t [nb_comps];
		data_for_S[dir].size = new size_t [nb_comps];
      for (size_t comp = 0; comp < nb_comps; comp++) {
         // Get local min and max indices
         size_t_vector min_unknown_index(dim,0);
         size_t_vector max_unknown_index(dim,0);
         for (size_t l=0;l<dim;++l) {
            min_unknown_index(l) =
								TF->get_min_index_unknown_handled_by_proc( comp, l ) ;
            max_unknown_index(l) =
								TF->get_max_index_unknown_handled_by_proc( comp, l ) ;
         }

			size_t local_min_j = (dir == 0) ? min_unknown_index(1)
													  : min_unknown_index(0);
			size_t local_max_j = (dir == 0) ? max_unknown_index(1)
													  : max_unknown_index(0);
			size_t local_min_k = (dim == 2) ? 0 :
									  ((dir == 2) ? min_unknown_index(1)
													  : min_unknown_index(2));
			size_t local_max_k = (dim == 2) ? 0 :
									  ((dir == 2) ? max_unknown_index(1)
													  : max_unknown_index(2));

			size_t local_length_j = (local_max_j-local_min_j+1);
			size_t local_length_k = (local_max_k-local_min_k+1);

         if (dim != 3) {
            first_pass[dir].size[comp] = 3*local_length_j;
            second_pass[dir].size[comp] = 2*local_length_j;
				data_for_S[dir].size[comp] = is_iperiodic[dir] ?
						 (size_t)(pow(nb_ranks_comm_i[dir],2) + 1)*local_length_j
					  : (size_t)(pow(nb_ranks_comm_i[dir]-1,2) + 1)*local_length_j ;
         } else if (dim == 3) {
            first_pass[dir].size[comp] = 3*local_length_j*local_length_k;
            second_pass[dir].size[comp] = 2*local_length_j*local_length_k;
				data_for_S[dir].size[comp] = is_iperiodic[dir] ?
			    (size_t)(pow(nb_ranks_comm_i[dir],2) + 1)*local_length_j*local_length_k
			  : (size_t)(pow(nb_ranks_comm_i[dir]-1,2) + 1)*local_length_j*local_length_k ;
         }
      }
   }

   // Array declarations
   for (size_t dir = 0; dir < dim; dir++) {
      first_pass[dir].send = new double** [nb_comps];
      first_pass[dir].receive = new double** [nb_comps];
      second_pass[dir].send = new double** [nb_comps];
      second_pass[dir].receive = new double** [nb_comps];
		data_for_S[dir].send = new double** [nb_comps];
		data_for_S[dir].receive = new double** [nb_comps];
      for (size_t comp = 0; comp < nb_comps; comp++) {
         first_pass[dir].send[comp] = new double* [1];
         first_pass[dir].receive[comp] = new double* [nb_ranks_comm_i[dir]];
         second_pass[dir].send[comp] = new double* [nb_ranks_comm_i[dir]];
         second_pass[dir].receive[comp] = new double* [1];
			data_for_S[dir].send[comp] = new double* [1];
	      data_for_S[dir].receive[comp] = new double* [1];
			data_for_S[dir].send[comp][0] =
								new double[data_for_S[dir].size[comp]];
			data_for_S[dir].receive[comp][0] =
								new double[data_for_S[dir].size[comp]];
			first_pass[dir].send[comp][0] =
								new double[first_pass[dir].size[comp]];
			second_pass[dir].receive[comp][0] =
								new double[second_pass[dir].size[comp]];
         for (size_t i = 0; i < (size_t)nb_ranks_comm_i[dir]; i++) {
            first_pass[dir].receive[comp][i] =
								new double[first_pass[dir].size[comp]];
            second_pass[dir].send[comp][i] =
								new double[second_pass[dir].size[comp]];
         }
      }
   }
}

//---------------------------------------------------------------------------
void
DS_HeatTransfer:: deallocate_mpi_variables ( void )
//---------------------------------------------------------------------------
{
   // Array declarations
   for (size_t dir = 0; dir < dim; dir++) {
      for (size_t comp = 0; comp < nb_comps; comp++) {
			delete [] data_for_S[dir].send[comp][0];
			delete [] data_for_S[dir].receive[comp][0];
			delete [] first_pass[dir].send[comp][0];
			delete [] second_pass[dir].receive[comp][0];
			for (size_t i = 0; i < (size_t) nb_ranks_comm_i[dir]; i++) {
            delete [] first_pass[dir].receive[comp][i];
            delete [] second_pass[dir].send[comp][i];
         }
         delete [] first_pass[dir].send[comp];
         delete [] first_pass[dir].receive[comp];
         delete [] second_pass[dir].send[comp];
         delete [] second_pass[dir].receive[comp];
			delete [] data_for_S[dir].send[comp];
			delete [] data_for_S[dir].receive[comp];
      }
      delete [] first_pass[dir].send;
      delete [] first_pass[dir].receive;
      delete [] second_pass[dir].send;
      delete [] second_pass[dir].receive;
		delete [] data_for_S[dir].send;
		delete [] data_for_S[dir].receive;
      delete [] first_pass[dir].size;
      delete [] second_pass[dir].size;
		delete [] data_for_S[dir].size;
   }
}

//---------------------------------------------------------------------------
void
DS_HeatTransfer:: free_DS_subcommunicators ( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DS_HeatTransfer:: free_DS_subcommunicators" ) ;


}

//----------------------------------------------------------------------
double DS_HeatTransfer:: assemble_advection_TVD(
										FV_DiscreteField const* AdvectingField
									 ,	size_t advecting_level
									 , double const& coef
									 , size_t const& i
									 , size_t const& j
									 , size_t const& k
									 , size_t advected_level) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "DS_HeatTransfer:: assemble_advection_TVD" );
   MAC_CHECK_PRE( advecting_level < AdvectingField->storage_depth() ) ;
   MAC_ASSERT( AdvectingField->discretization_type() == "staggered" ) ;

   // Parameters
   size_t component = 0 ;
   double xC = 0., yC = 0., zC = 0.,
   	xr = 0., xR = 0., xl = 0., xL = 0., yt = 0., yT = 0., yb = 0., yB = 0.,
	zf = 0., zF = 0., zb = 0., zB = 0.;
   double dxC = 0., dyC = 0., dzC = 0.,
   	dxr = 0., dxl = 0., dxCr = 0., dxCl = 0., dxRr = 0., dxR = 0.,
	dxLl = 0., dyt = 0., dyb = 0., dyCt = 0., dyCb = 0., dyTt = 0.,
	dyT = 0., dyBb = 0., dzf = 0., dzb = 0., dzCf = 0., dzCb = 0.,
	dzFf = 0., dzF = 0., dzBb = 0.;

   double AdvectedvalueC = 0., AdvectedvalueRi = 0., AdvectedvalueRiRi = 0.,
   	AdvectedvalueLe = 0.,  AdvectedvalueLeLe = 0., AdvectedvalueTo = 0.,
	AdvectedvalueToTo = 0., AdvectedvalueBo = 0., AdvectedvalueBoBo = 0.,
   	AdvectedvalueFr = 0., AdvectedvalueFrFr = 0., AdvectedvalueBe = 0.,
	AdvectedvalueBeBe = 0.;
   double ur = 0., ul = 0., vt = 0., vb = 0., wf = 0., wb = 0.,
	fri = 0., fle = 0., fto = 0., fbo = 0., ffr = 0., fbe = 0., flux = 0.;
   double cRip12 = 0., cLip12 = 0., cRim12 = 0., cLim12 = 0.,
   	thetaC = 0., thetaRi = 0., thetaLe = 0., thetaTo = 0., thetaBo = 0.,
	thetaFr = 0., thetaBe = 0.;

   FV_SHIFT_TRIPLET shift = TF->shift_staggeredToCentered() ;

   // Perform assembling

   xC = TF->get_DOF_coordinate( i, component, 0 );
   dxC = TF->get_cell_size( i, component, 0 ) ;

   yC = TF->get_DOF_coordinate( j, component, 1 );
   dyC = TF->get_cell_size( j, component, 1 ) ;

   if ( dim == 2 ) {
      AdvectedvalueC = TF->DOF_value( i, j, k, component, advected_level );

      // Right and Left
      // --------------
      AdvectedvalueRi = TF->DOF_value( i+1, j, k, component, advected_level );
      AdvectedvalueLe = TF->DOF_value( i-1, j, k, component, advected_level );

      thetaC = fabs( AdvectedvalueRi - AdvectedvalueC ) > 1.e-20  ?
                   ( AdvectedvalueC - AdvectedvalueLe )
		 / ( AdvectedvalueRi - AdvectedvalueC ) : 1e20 ;

      // Right (X)
      ur = AdvectingField->DOF_value( i+shift.i, j, k, 0, advecting_level );

      if ( TF->DOF_color( i+1, j, k, component ) == FV_BC_RIGHT ) {
         if ( ur > 0. ) fri = ur * AdvectedvalueC;
         else fri = ur * AdvectedvalueRi;
      } else {
         xr = AdvectingField->get_DOF_coordinate( i+shift.i, 0, 0 );
	 xR = TF->get_DOF_coordinate( i+1, component, 0 );
         dxCr = xr - xC;
         dxr  = xR - xC;
         dxRr = xR - xr;
         dxR = TF->get_cell_size( i+1, component, 0 );
         AdvectedvalueRiRi = TF->DOF_value( i+2, j, k, component, advected_level );

         thetaRi = fabs( AdvectedvalueRiRi - AdvectedvalueRi) > 1.e-20 ?
 	               ( AdvectedvalueRi - AdvectedvalueC)
		     / ( AdvectedvalueRiRi - AdvectedvalueRi) : 1e20 ;
         cRip12 = AdvectedvalueRi
		- ( dxRr / dxR ) * FV_DiscreteField::SuperBee_phi(thetaRi)
		* ( AdvectedvalueRiRi - AdvectedvalueRi );
         cLip12 = AdvectedvalueC + ( dxCr / dxr )
	   	* FV_DiscreteField::SuperBee_phi( thetaC )
		* ( AdvectedvalueRi - AdvectedvalueC );

         fri = 0.5 * ( ur * ( cRip12 + cLip12 )
		- fabs(ur) * ( cRip12 - cLip12 ) ) ;
      }

      // Left (X)
      ul = AdvectingField->DOF_value( i+shift.i-1, j, k, 0, advecting_level );

      if ( TF->DOF_color( i-1, j, k, component ) == FV_BC_LEFT ) {
         if ( ul > 0. ) fle = ul * AdvectedvalueLe;
         else fle = ul * AdvectedvalueC;
      } else {
         xl = AdvectingField->get_DOF_coordinate( i+shift.i-1, 0, 0 );
         xL = TF->get_DOF_coordinate( i-1, component, 0 );
         dxCl = xC - xl;
         dxl  = xC - xL;
         dxLl = xl - xL;
         AdvectedvalueLeLe = TF->DOF_value( i-2, j, k, component, advected_level );

         thetaLe = fabs( AdvectedvalueC - AdvectedvalueLe ) > 1.e-20 ?
 		       ( AdvectedvalueLe - AdvectedvalueLeLe )
		     / ( AdvectedvalueC - AdvectedvalueLe ) : 1e20 ;
         cLim12 = AdvectedvalueLe
		+ ( dxLl / dxl ) * FV_DiscreteField::SuperBee_phi( thetaLe )
		* ( AdvectedvalueC - AdvectedvalueLe ) ;
         if ( TF->DOF_color( i, j, k, component ) == FV_BC_RIGHT )
           cRim12 = AdvectedvalueC;
         else {
           xR = TF->get_DOF_coordinate( i+1, component, 0 );
           dxr  = xR - xC;
           cRim12 = AdvectedvalueC - ( dxCl / dxr )
	     	* FV_DiscreteField::SuperBee_phi( thetaC )
		* ( AdvectedvalueRi - AdvectedvalueC ) ;
         }

	   fle = 0.5 * ( ul * ( cRim12 + cLim12 )
			- fabs(ul) * ( cRim12 - cLim12 ) ) ;
      }

      // Top and Bottom
      // --------------
      AdvectedvalueTo = TF->DOF_value( i, j+1, k, component, advected_level );
      AdvectedvalueBo = TF->DOF_value( i, j-1, k, component, advected_level );

      thetaC = fabs( AdvectedvalueTo - AdvectedvalueC ) > 1.e-20 ?
    	           ( AdvectedvalueC - AdvectedvalueBo )
		 / ( AdvectedvalueTo - AdvectedvalueC) : 1e20 ;

      // Top (Y)
      vt = AdvectingField->DOF_value( i, j+shift.j, k, 1, advecting_level );

      if ( TF->DOF_color( i, j+1, k, component ) == FV_BC_TOP ) {
         if ( vt > 0. ) fto = vt * AdvectedvalueC;
         else fto = vt * AdvectedvalueTo;
      } else {
         yt = AdvectingField->get_DOF_coordinate( j+shift.j, 1, 1 );
         yT = TF->get_DOF_coordinate( j+1, component, 1 );
         dyCt = yt - yC;
         dyt  = yT - yC;
         dyTt = yT - yt;
	 dyT = TF->get_cell_size( j+1, component, 1 );
         AdvectedvalueToTo = TF->DOF_value( i, j+2, k, component, advected_level );

         thetaTo = fabs( AdvectedvalueToTo - AdvectedvalueTo ) > 1.e-20 ?
                       ( AdvectedvalueTo - AdvectedvalueC )
                     / ( AdvectedvalueToTo - AdvectedvalueTo) : 1e20 ;
         cRip12 = AdvectedvalueTo
		- ( dyTt / dyT ) * FV_DiscreteField::SuperBee_phi( thetaTo )
		* ( AdvectedvalueToTo - AdvectedvalueTo );
         cLip12 = AdvectedvalueC + ( dyCt / dyt )
	   	* FV_DiscreteField::SuperBee_phi( thetaC )
		* ( AdvectedvalueTo - AdvectedvalueC );

         fto = 0.5 * ( vt * ( cRip12 + cLip12 ) - fabs(vt) * ( cRip12 - cLip12 ) ) ;
      }

      // Bottom (Y)
      vb = AdvectingField->DOF_value( i, j+shift.j-1, k, 1, advecting_level );

      if ( TF->DOF_color( i, j-1, k, component ) == FV_BC_BOTTOM ) {
         if ( vb > 0. ) fbo = vb * AdvectedvalueBo;
         else fbo = vb * AdvectedvalueC;
      } else {
         yb = AdvectingField->get_DOF_coordinate( j+shift.j-1, 1, 1 );
         yB = TF->get_DOF_coordinate( j-1, component, 1 );
         dyCb = yC - yb;
         dyb  = yC - yB;
         dyBb = yb - yB;
         AdvectedvalueBoBo = TF->DOF_value( i, j-2, k, component, advected_level );

         thetaBo = fabs( AdvectedvalueC - AdvectedvalueBo ) > 1.e-20 ?
                       ( AdvectedvalueBo - AdvectedvalueBoBo )
		     / ( AdvectedvalueC - AdvectedvalueBo ) : 1e20 ;
         cLim12 = AdvectedvalueBo
		+ ( dyBb / dyb ) * FV_DiscreteField::SuperBee_phi( thetaBo )
		* ( AdvectedvalueC - AdvectedvalueBo );
         if ( TF->DOF_color( i, j, k, component ) == FV_BC_TOP )
            cRim12 = AdvectedvalueC;
         else {
            yT = TF->get_DOF_coordinate( j+1, component, 1 );
            dyt  = yT - yC;
            cRim12 = AdvectedvalueC - ( dyCb / dyt )
	     	* FV_DiscreteField::SuperBee_phi( thetaC )
        	* ( AdvectedvalueTo - AdvectedvalueC );
	 }

	 fbo = 0.5 * ( vb * ( cRim12 + cLim12 ) - fabs(vb) * ( cRim12 - cLim12 ) ) ;
      }

      flux = ( fto - fbo ) * dxC + ( fri - fle ) * dyC;

   } else {
      zC = TF->get_DOF_coordinate( k, component, 2 );
      dzC = TF->get_cell_size( k, component, 2 ) ;

      AdvectedvalueC = TF->DOF_value( i, j, k, component, advected_level );

      // Right and Left
      // --------------
      AdvectedvalueRi = TF->DOF_value( i+1, j, k, component, advected_level );
      AdvectedvalueLe = TF->DOF_value( i-1, j, k, component, advected_level );

      thetaC = fabs( AdvectedvalueRi - AdvectedvalueC) > 1.e-20 ?
                   ( AdvectedvalueC - AdvectedvalueLe )
                 / ( AdvectedvalueRi - AdvectedvalueC) : 1e20 ;


      // Right (X)
      ur = AdvectingField->DOF_value( i+shift.i, j, k, 0, advecting_level );


      if ( TF->DOF_color( i+1, j, k, component ) == FV_BC_RIGHT ) {
          if ( ur > 0. ) fri = ur * AdvectedvalueC;
          else fri = ur * AdvectedvalueRi;
      } else {
          xr = AdvectingField->get_DOF_coordinate( i+shift.i, 0, 0 );
          xR = TF->get_DOF_coordinate( i+1, component, 0 );
          dxCr = xr - xC;
          dxr  = xR - xC;
          dxRr = xR - xr;
          dxR = TF->get_cell_size( i+1, component, 0 );

          AdvectedvalueRiRi = TF->DOF_value( i+2, j, k, component, advected_level );

          thetaRi = fabs( AdvectedvalueRiRi - AdvectedvalueRi ) > 1.e-20 ?
                        ( AdvectedvalueRi - AdvectedvalueC ) \
                      / ( AdvectedvalueRiRi - AdvectedvalueRi ) : 1e20 ;
          cRip12 = AdvectedvalueRi
		- ( dxRr / dxR ) * FV_DiscreteField::SuperBee_phi( thetaRi )
		* ( AdvectedvalueRiRi - AdvectedvalueRi );
          cLip12 = AdvectedvalueC + ( dxCr / dxr )
	     	* FV_DiscreteField::SuperBee_phi( thetaC )
		* ( AdvectedvalueRi - AdvectedvalueC );

          fri = 0.5 * ( ur * ( cRip12 + cLip12 ) - fabs(ur) * ( cRip12 - cLip12 ) ) ;
      }

      // Left (X)
      ul = AdvectingField->DOF_value( i+shift.i-1, j, k, 0, advecting_level );

      if ( TF->DOF_color( i-1, j, k, component ) == FV_BC_LEFT ) {
         if ( ul > 0. ) fle = ul * AdvectedvalueLe;
         else fle = ul * AdvectedvalueC;
      } else {
         xl = AdvectingField->get_DOF_coordinate( i+shift.i-1, 0, 0 );
         xL = TF->get_DOF_coordinate( i-1, component, 0 );
         dxCl = xC - xl;
         dxl  = xC - xL;
         dxLl = xl - xL;
         AdvectedvalueLeLe = TF->DOF_value( i-2, j, k, component, advected_level );

         thetaLe = fabs( AdvectedvalueC - AdvectedvalueLe) > 1.e-20 ?
	               ( AdvectedvalueLe - AdvectedvalueLeLe )
		     / ( AdvectedvalueC - AdvectedvalueLe) : 1e20 ;
	 cLim12 = AdvectedvalueLe
		+ ( dxLl / dxl ) * FV_DiscreteField::SuperBee_phi( thetaLe )
		* ( AdvectedvalueC - AdvectedvalueLe ) ;
         if ( TF->DOF_color( i, j, k, component ) == FV_BC_RIGHT )
            cRim12 = AdvectedvalueC;
         else {
//            xR = TF->get_DOF_coordinate( i+1, advected_level, 0 );
            xR = TF->get_DOF_coordinate( i+1, component, 0 );
            dxr  = xR - xC;
            cRim12 = AdvectedvalueC - ( dxCl / dxr )
			* FV_DiscreteField::SuperBee_phi( thetaC )
			* ( AdvectedvalueRi - AdvectedvalueC );
         }

         fle = 0.5 * ( ul * ( cRim12 + cLim12 )	- fabs(ul) * ( cRim12 - cLim12 ) ) ;
      }

      // Top and Bottom
      // --------------
      AdvectedvalueTo = TF->DOF_value( i, j+1, k, component, advected_level );
      AdvectedvalueBo = TF->DOF_value( i, j-1, k, component, advected_level );

      thetaC = fabs( AdvectedvalueTo - AdvectedvalueC ) > 1.e-20 ?
	   	   ( AdvectedvalueC - AdvectedvalueBo )
		 / ( AdvectedvalueTo - AdvectedvalueC ) : 1e20 ;

      // Top (Y)
      vt = AdvectingField->DOF_value( i, j+shift.j, k, 1, advecting_level );

      if ( TF->DOF_color( i, j+1, k, component ) == FV_BC_TOP ) {
         if ( vt > 0. ) fto = vt * AdvectedvalueC;
         else fto = vt * AdvectedvalueTo;
      } else {
         yt = AdvectingField->get_DOF_coordinate( j+shift.j, 1, 1 );
         yT = TF->get_DOF_coordinate( j+1, component, 1 );
         dyCt = yt - yC;
         dyt  = yT - yC;
         dyTt = yT - yt;
         dyT = TF->get_cell_size( j+1, component, 1 );
         AdvectedvalueToTo = TF->DOF_value( i, j+2, k, component, advected_level );

         thetaTo = fabs( AdvectedvalueToTo - AdvectedvalueTo) > 1.e-20 ?
	               ( AdvectedvalueTo - AdvectedvalueC )
		     / ( AdvectedvalueToTo - AdvectedvalueTo ) : 1e20 ;
         cRip12 = AdvectedvalueTo
		- ( dyTt / dyT ) * FV_DiscreteField::SuperBee_phi( thetaTo )
		* ( AdvectedvalueToTo - AdvectedvalueTo ) ;
         cLip12 = AdvectedvalueC + ( dyCt / dyt )
		* FV_DiscreteField::SuperBee_phi(thetaC)
		* ( AdvectedvalueTo - AdvectedvalueC ) ;

         fto = 0.5 * ( vt * ( cRip12 + cLip12 ) - fabs(vt) * ( cRip12 - cLip12 ) ) ;
      }

      // Bottom (Y)
      vb = AdvectingField->DOF_value( i, j+shift.j-1, k, 1, advecting_level );

      if ( TF->DOF_color( i, j-1, k, component ) == FV_BC_BOTTOM ) {
         if ( vb > 0. ) fbo = vb * AdvectedvalueBo;
         else fbo = vb * AdvectedvalueC;
      } else {
          yb = AdvectingField->get_DOF_coordinate( j+shift.j-1, 1, 1 );
          yB = TF->get_DOF_coordinate( j-1, component, 1 );
          dyCb = yC - yb;
          dyb  = yC - yB;
          dyBb = yb - yB;
          AdvectedvalueBoBo = TF->DOF_value( i, j-2, k, component, advected_level );

          thetaBo = fabs( AdvectedvalueC - AdvectedvalueBo) > 1.e-20 ?
                        ( AdvectedvalueBo - AdvectedvalueBoBo )
		      / ( AdvectedvalueC - AdvectedvalueBo ) : 1e20 ;
          cLim12 = AdvectedvalueBo
		+ ( dyBb / dyb ) * FV_DiscreteField::SuperBee_phi( thetaBo )
		* ( AdvectedvalueC - AdvectedvalueBo ) ;
          if ( TF->DOF_color( i, j, k, component ) == FV_BC_TOP )
             cRim12 = AdvectedvalueC;
          else {
             yT = TF->get_DOF_coordinate( j+1, component, 1 );
             dyt  = yT - yC;
             cRim12 = AdvectedvalueC - ( dyCb / dyt )
			* FV_DiscreteField::SuperBee_phi( thetaC )
			* ( AdvectedvalueTo - AdvectedvalueC ) ;
          }

          fbo = 0.5 * ( vb * ( cRim12 + cLim12 ) - fabs(vb) * ( cRim12 - cLim12 ) ) ;
      }

      // Front and Behind
      // ----------------
      AdvectedvalueFr = TF->DOF_value( i, j, k+1, component, advected_level );
      AdvectedvalueBe = TF->DOF_value( i, j, k-1, component, advected_level );

      thetaC = fabs( AdvectedvalueFr - AdvectedvalueC) > 1.e-20 ?
    	           ( AdvectedvalueC - AdvectedvalueBe )
		 / ( AdvectedvalueFr - AdvectedvalueC ) : 1e20 ;

      // Front (Z)
      wf = AdvectingField->DOF_value( i, j, k+shift.k, 2, advecting_level );

      if ( TF->DOF_color( i, j, k+1, component ) == FV_BC_FRONT ) {
         if ( wf > 0. ) ffr = wf * AdvectedvalueC;
         else ffr = wf * AdvectedvalueFr;
      } else {
         zf = AdvectingField->get_DOF_coordinate( k+shift.k, 2, 2 );
         zF = TF->get_DOF_coordinate( k+1, component, 2 );
         dzCf = zf - zC;
         dzf  = zF - zC;
         dzFf = zF - zf;
         dzF = TF->get_cell_size( k+1, component, 2 );
         AdvectedvalueFrFr = TF->DOF_value( i, j, k+2, component, advected_level );

         thetaFr = fabs( AdvectedvalueFrFr - AdvectedvalueFr) > 1.e-20 ?
		       ( AdvectedvalueFr - AdvectedvalueC )
		     / ( AdvectedvalueFrFr - AdvectedvalueFr ) : 1e20 ;
	 cRip12 = AdvectedvalueFr
		- ( dzFf / dzF ) * FV_DiscreteField::SuperBee_phi( thetaFr )
		* ( AdvectedvalueFrFr - AdvectedvalueFr ) ;
         cLip12 = AdvectedvalueC + ( dzCf / dzf )
		* FV_DiscreteField::SuperBee_phi( thetaC )
		* ( AdvectedvalueFr - AdvectedvalueC ) ;

	 ffr = 0.5 * ( wf * ( cRip12 + cLip12 ) - fabs(wf) * ( cRip12 - cLip12 ) ) ;
      }

      // Behind (Z)
      wb = AdvectingField->DOF_value( i, j, k+shift.k-1, 2, advecting_level );

      if ( TF->DOF_color( i, j, k-1, component ) == FV_BC_BEHIND ) {
         if ( wb > 0. ) fbe = wb * AdvectedvalueBe;
         else fbe = wb * AdvectedvalueC;
      } else {
          zb = AdvectingField->get_DOF_coordinate( k+shift.k-1, 2, 2 );
          zB = TF->get_DOF_coordinate( k-1, component, 2 );
          dzCb = zC - zb;
          dzb  = zC - zB;
          dzBb = zb - zB;
          AdvectedvalueBeBe = TF->DOF_value( i, j, k-2, component, advected_level );

          thetaBe = fabs( AdvectedvalueC - AdvectedvalueBe) > 1.e-20 ?
		        ( AdvectedvalueBe - AdvectedvalueBeBe )
		      / ( AdvectedvalueC - AdvectedvalueBe ) : 1e20 ;
          cLim12 = AdvectedvalueBe
		+ ( dzBb / dzb ) * FV_DiscreteField::SuperBee_phi( thetaBe )
		* ( AdvectedvalueC - AdvectedvalueBe ) ;
          if ( TF->DOF_color( i, j, k, component ) == FV_BC_FRONT )
	     cRim12 = AdvectedvalueC;
	  else {
             zF = TF->get_DOF_coordinate( k+1, component, 2 );
             dzf  = zF - zC;
	     cRim12 = AdvectedvalueC - ( dzCb / dzf )
	       		* FV_DiscreteField::SuperBee_phi( thetaC )
			* ( AdvectedvalueFr - AdvectedvalueC ) ;
	  }

	  fbe = 0.5 * ( wb * ( cRim12 + cLim12 ) - fabs(wb) * ( cRim12 - cLim12 ) ) ;
      }

      flux = ( fto - fbo ) * dxC * dzC + ( fri - fle ) * dyC * dzC + ( ffr - fbe ) * dxC * dyC;
   }

   return (coef * flux);

}




//----------------------------------------------------------------------
double DS_HeatTransfer:: assemble_advection_Centered(
											FV_DiscreteField const* AdvectingField
										 , size_t advecting_level
										 , double const& coef
										 , size_t const& i
										 , size_t const& j
										 , size_t const& k
										 , size_t advected_level) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "DS_HeatTransfer:: assemble_advection_Centered" );
   MAC_CHECK_PRE( advecting_level < AdvectingField->storage_depth() ) ;
   MAC_ASSERT( AdvectingField->discretization_type() == "staggered" ) ;

   // Parameters
   size_t component = 0 ;
   double dxC = 0., dyC = 0., dzC = 0.;

   double AdvectedvalueC = 0., AdvectedvalueRi = 0.,
   	AdvectedvalueLe = 0., AdvectedvalueTo = 0.,
	AdvectedvalueBo = 0., AdvectedvalueFr = 0., AdvectedvalueBe = 0.;
   double ur = 0., ul = 0., vt = 0., vb = 0., wf = 0., wb = 0.,
	fri = 0., fle = 0., fto = 0., fbo = 0., ffr = 0., fbe = 0., flux = 0.;

   FV_SHIFT_TRIPLET shift = TF->shift_staggeredToCentered() ;

   // Perform assembling

   dxC = TF->get_cell_size( i, component, 0 ) ;

   dyC = TF->get_cell_size( j, component, 1 ) ;

   if ( dim == 2 ) {
      AdvectedvalueC = TF->DOF_value( i, j, k, component, advected_level );

      // Right and Left
      // --------------
      AdvectedvalueRi = TF->DOF_value( i+1, j, k, component, advected_level );
      AdvectedvalueLe = TF->DOF_value( i-1, j, k, component, advected_level );

      // Right (X)
      ur = AdvectingField->DOF_value( i+shift.i, j, k, 0, advecting_level );

      if ( TF->DOF_color( i+1, j, k, component ) == FV_BC_RIGHT ) {
         fri = ur * AdvectedvalueRi;
      } else {
         fri = ur * 0.5*(AdvectedvalueRi + AdvectedvalueC);
      }

      // Left (X)
      ul = AdvectingField->DOF_value( i+shift.i-1, j, k, 0, advecting_level );

      if ( TF->DOF_color( i-1, j, k, component ) == FV_BC_LEFT ) {
         fle = ul * AdvectedvalueLe;
      } else {
         fle = ul * 0.5*(AdvectedvalueLe + AdvectedvalueC);
      }

      // Top and Bottom
      // --------------
      AdvectedvalueTo = TF->DOF_value( i, j+1, k, component, advected_level );
      AdvectedvalueBo = TF->DOF_value( i, j-1, k, component, advected_level );

      // Top (Y)
      vt = AdvectingField->DOF_value( i, j+shift.j, k, 1, advecting_level );

      if ( TF->DOF_color( i, j+1, k, component ) == FV_BC_TOP ) {
         fto = vt * AdvectedvalueTo;
      } else {
         fto = vt * 0.5*(AdvectedvalueTo + AdvectedvalueC);
      }

      // Bottom (Y)
      vb = AdvectingField->DOF_value( i, j+shift.j-1, k, 1, advecting_level );

      if ( TF->DOF_color( i, j-1, k, component ) == FV_BC_BOTTOM ) {
         fbo = vb * AdvectedvalueBo;
      } else {
         fbo = vb * 0.5*(AdvectedvalueBo + AdvectedvalueC);
      }

      flux = ( fto - fbo ) * dxC + ( fri - fle ) * dyC;

   } else {
      dzC = TF->get_cell_size( k, component, 2 ) ;

      AdvectedvalueC = TF->DOF_value( i, j, k, component, advected_level );

      // Right and Left
      // --------------
      AdvectedvalueRi = TF->DOF_value( i+1, j, k, component, advected_level );
      AdvectedvalueLe = TF->DOF_value( i-1, j, k, component, advected_level );

      // Right (X)
      ur = AdvectingField->DOF_value( i+shift.i, j, k, 0, advecting_level );

      if ( TF->DOF_color( i+1, j, k, component ) == FV_BC_RIGHT ) {
         fri = ur * AdvectedvalueRi;
      } else {
         fri = ur * 0.5*(AdvectedvalueRi + AdvectedvalueC);
      }

      // Left (X)
      ul = AdvectingField->DOF_value( i+shift.i-1, j, k, 0, advecting_level );

      if ( TF->DOF_color( i-1, j, k, component ) == FV_BC_LEFT ) {
         fle = ul * AdvectedvalueLe;
      } else {
         fle = ul * 0.5*(AdvectedvalueLe + AdvectedvalueC);
      }

      // Top and Bottom
      // --------------
      AdvectedvalueTo = TF->DOF_value( i, j+1, k, component, advected_level );
      AdvectedvalueBo = TF->DOF_value( i, j-1, k, component, advected_level );

      // Top (Y)
      vt = AdvectingField->DOF_value( i, j+shift.j, k, 1, advecting_level );

      if ( TF->DOF_color( i, j+1, k, component ) == FV_BC_TOP ) {
         fto = vt * AdvectedvalueTo;
      } else {
         fto = vt * 0.5*(AdvectedvalueTo + AdvectedvalueC);
      }

      // Bottom (Y)
      vb = AdvectingField->DOF_value( i, j+shift.j-1, k, 1, advecting_level );

      if ( TF->DOF_color( i, j-1, k, component ) == FV_BC_BOTTOM ) {
         fbo = vb * AdvectedvalueBo;
      } else {
         fbo = vb * 0.5*(AdvectedvalueBo + AdvectedvalueC);
      }

      // Front and Behind
      // ----------------
      AdvectedvalueFr = TF->DOF_value( i, j, k+1, component, advected_level );
      AdvectedvalueBe = TF->DOF_value( i, j, k-1, component, advected_level );

      // Front (Z)
      wf = AdvectingField->DOF_value( i, j, k+shift.k, 2, advecting_level );

      if ( TF->DOF_color( i, j, k+1, component ) == FV_BC_FRONT ) {
         ffr = wf * AdvectedvalueFr;
      } else {
         ffr = wf * 0.5*(AdvectedvalueFr + AdvectedvalueC);
      }

      // Behind (Z)
      wb = AdvectingField->DOF_value( i, j, k+shift.k-1, 2, advecting_level );

      if ( TF->DOF_color( i, j, k-1, component ) == FV_BC_BEHIND ) {
         fbe = wb * AdvectedvalueBe;
      } else {
         fbe = wb * 0.5*(AdvectedvalueBe + AdvectedvalueC);
      }

      flux = ( fto - fbo ) * dxC * dzC + ( fri - fle ) * dyC * dzC + ( ffr - fbe ) * dxC * dyC;
   }

   return (coef * flux);

}

//----------------------------------------------------------------------
double DS_HeatTransfer:: assemble_advection_Upwind(
											FV_DiscreteField const* AdvectingField
										 , size_t advecting_level
										 , double const& coef
										 , size_t const& i
										 , size_t const& j
										 , size_t const& k
										 , size_t advected_level) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "DS_HeatTransfer:: assemble_advection_Upwind" );
   MAC_CHECK_PRE( advecting_level < AdvectingField->storage_depth() ) ;
   MAC_ASSERT( AdvectingField->discretization_type() == "staggered" ) ;

   // Parameters
   size_t component = 0 ;
   double dxC = 0., dyC = 0., dzC = 0.;
   double AdvectedvalueC = 0., AdvectedvalueRi = 0., AdvectedvalueLe = 0.,
   	AdvectedvalueTo = 0., AdvectedvalueBo = 0.,
   	AdvectedvalueFr = 0., AdvectedvalueBe = 0,
	ur = 0., ul = 0., vt = 0., vb = 0., wf = 0., wb = 0.,
	fri = 0., fle = 0., fto = 0., fbo = 0., ffr = 0., fbe = 0., flux = 0.;

   // Comment: cell centered unknowns always have a defined value at +1/-1
   // indices in all 3 directions. Whether one of the +1/-1 DOF values is on a
   // boundary or not, and whether that boundary has a Dirichlet or Neumann
   // condition is irrelevant, this +1/-1 DOF always has the right value.
   // For Neumann, this is guaranted by
   // FV_BoundaryCondition:: set_free_DOF_values in
   // FV_DiscreteField:: update_free_DOFs_value or
   // FV_DiscreteField:: add_to_free_DOFs_value

   FV_SHIFT_TRIPLET shift = TF->shift_staggeredToCentered() ;

   dxC = TF->get_cell_size( i, component, 0 ) ;
   dyC = TF->get_cell_size( j, component, 1 ) ;

   if ( dim == 2 ) {
      AdvectedvalueC = TF->DOF_value( i, j, k, component, advected_level );

      // Right (X)
      AdvectedvalueRi = TF->DOF_value( i+1, j, k, component, advected_level );
      ur = AdvectingField->DOF_value( i+shift.i, j, k, 0, advecting_level );
      if ( ur > 0. ) fri = ur * AdvectedvalueC;
      else fri = ur * AdvectedvalueRi;

      // Left (X)
      AdvectedvalueLe = TF->DOF_value( i-1, j, k, component, advected_level );
      ul = AdvectingField->DOF_value( i+shift.i-1, j, k, 0, advecting_level );
      if ( ul > 0. ) fle = ul * AdvectedvalueLe;
      else fle = ul * AdvectedvalueC;

      // Top (Y)
      AdvectedvalueTo = TF->DOF_value( i, j+1, k, component, advected_level );
      vt = AdvectingField->DOF_value( i, j+shift.j, k, 1, advecting_level );
      if ( vt > 0. ) fto = vt * AdvectedvalueC;
      else fto = vt * AdvectedvalueTo;

      // Bottom (Y)
      AdvectedvalueBo = TF->DOF_value( i, j-1, k, component, advected_level );
      vb = AdvectingField->DOF_value( i, j+shift.j-1, k, 1, advecting_level );
      if ( vb > 0. ) fbo = vb * AdvectedvalueBo;
      else fbo = vb * AdvectedvalueC;

      flux = (fto - fbo) * dxC + (fri - fle) * dyC;

   } else {
      dzC = TF->get_cell_size( k, component, 2);
      AdvectedvalueC = TF->DOF_value( i, j, k, component, advected_level );

      // Right (X)
      AdvectedvalueRi = TF->DOF_value( i+1, j, k, component, advected_level );
      ur = AdvectingField->DOF_value( i+shift.i, j, k, 0, advecting_level );
      if ( ur > 0. ) fri = ur * AdvectedvalueC;
      else fri = ur * AdvectedvalueRi;

      // Left (X)
      AdvectedvalueLe = TF->DOF_value( i-1, j, k, component, advected_level );
      ul = AdvectingField->DOF_value( i+shift.i-1, j, k, 0, advecting_level );
      if ( ul > 0. ) fle = ul * AdvectedvalueLe;
      else fle = ul * AdvectedvalueC;

      // Top (Y)
      AdvectedvalueTo = TF->DOF_value( i, j+1, k, component, advected_level );
      vt = AdvectingField->DOF_value( i, j+shift.j, k, 1, advecting_level );
      if ( vt > 0. ) fto = vt * AdvectedvalueC;
      else fto = vt * AdvectedvalueTo;

      // Bottom (Y)
      AdvectedvalueBo = TF->DOF_value( i, j-1, k, component, advected_level );
      vb = AdvectingField->DOF_value( i, j+shift.j-1, k, 1, advecting_level );
      if ( vb > 0. ) fbo = vb * AdvectedvalueBo;
      else fbo = vb * AdvectedvalueC;

      // Front (Z)
      AdvectedvalueFr = TF->DOF_value( i, j, k+1, component, advected_level );
      wf = AdvectingField->DOF_value( i, j, k+shift.k, 2, advecting_level );
      if ( wf > 0. ) ffr = wf * AdvectedvalueC;
      else ffr = wf * AdvectedvalueFr;

      // Behind (Z)
      AdvectedvalueBe = TF->DOF_value( i, j, k-1, component, advected_level );
      wb = AdvectingField->DOF_value( i, j, k+shift.k-1, 2, advecting_level );
      if ( wb > 0. ) fbe = wb * AdvectedvalueBe;
      else fbe = wb * AdvectedvalueC;

      flux = (fto - fbo) * dxC * dzC + (fri - fle) * dyC * dzC + (ffr - fbe) * dxC * dyC;
   }

   return (coef * flux);
}
