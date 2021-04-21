#include <DDS_HeatTransfer.hh>
#include <FV_DomainAndFields.hh>
#include <FV_DiscreteField.hh>
#include <FV_DomainBuilder.hh>
#include <DDS_HeatTransferSystem.hh>
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

//---------------------------------------------------------------------------
DDS_HeatTransfer*
DDS_HeatTransfer:: create( MAC_Object* a_owner,
		MAC_ModuleExplorer const* exp,
                struct NavierStokes2Temperature const& transfer )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransfer:: create" ) ;
   MAC_CHECK_PRE( exp != 0 ) ;

   DDS_HeatTransfer* result =
                        new DDS_HeatTransfer( a_owner, exp, transfer ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;

   return( result ) ;

}

//---------------------------------------------------------------------------
DDS_HeatTransfer:: DDS_HeatTransfer( MAC_Object* a_owner,
		MAC_ModuleExplorer const* exp,
                struct NavierStokes2Temperature const& fromNS )
//---------------------------------------------------------------------------
   : MAC_Object( a_owner )
   , ComputingTime("Solver")
   , TF ( fromNS.dom_->discrete_field( "temperature" ) )
   , UF (fromNS.UF_)
   , TF_DS_ERROR( 0 )
   , GLOBAL_EQ( 0 )
   , rho( fromNS.rho_ )
   , AdvectionScheme ( fromNS.AdvectionScheme_ )
   , AdvectionTimeAccuracy ( fromNS.AdvectionTimeAccuracy_ ) 
   , ViscousStressOrder ( fromNS.ViscousStressOrder_ )
   , heat_capacity( exp->double_data( "Heat_capacity") )
   , thermal_conductivity( exp->double_data( "Thermal_conductivity") )
   , b_bodyterm( false )
   , is_solids ( fromNS.is_solids_ )
{
   MAC_LABEL( "DDS_HeatTransfer:: DDS_HeatTransfer" ) ;

   MAC_ASSERT( TF->discretization_type() == "centered" ) ;
   MAC_ASSERT( TF->storage_depth() == 5 ) ;

   // Call of MAC_Communicator routine to set the rank of each proces and
   // the number of processes during execution
   pelCOMM = MAC_Exec::communicator();
   my_rank = pelCOMM->rank();
   nb_procs = pelCOMM->nb_ranks();
   is_master = 0;
   is_iperiodic[0] = false;
   is_iperiodic[1] = false;
   is_iperiodic[2] = false;

   // Timing routines
   if ( my_rank == is_master )
   {
     CT_set_start();
     SCT_insert_app("Objects_Creation");
     SCT_set_start("Objects_Creation");
   }


   // Is the run a follow up of a previous job
//   b_restart = MAC_Application::is_follow();
   b_restart = fromNS.b_restart_ ;
   // Clear results directory in case of a new run
   if( !b_restart ) PAC_Misc::clearAllFiles( "Res", "Savings", my_rank ) ;

   // Get space dimension
   dim = TF->primary_grid()->nb_space_dimensions() ;
   nb_comps = TF->nb_components() ;

   if ( dim == 1 )
   {
     string error_message="Space dimension should either 2 or 3";
     MAC_Error::object()->raise_bad_data_value(exp,
	"nb_space_dimensions",
	error_message );
   }


   // Create the Direction Splitting subcommunicators
   create_DDS_subcommunicators();

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

   is_stressCal = fromNS.is_stressCal_;
   is_par_motion = fromNS.is_par_motion_;
   particle_information = fromNS.particle_information_;
   insertion_type = fromNS.insertion_type_;

   if (is_solids) {
      Npart = fromNS.Npart_ ;
      loc_thres = fromNS.loc_thres_ ;
      level_set_type = fromNS.level_set_type_ ;

      if (is_stressCal) {
         Npoints = fromNS.Npoints_ ;
         if (dim == 3) {
            if ((level_set_type == "Sphere") || (level_set_type == "Cylinder")) {
               Pmin = fromNS.Pmin_ ;
               ar = fromNS.ar_ ;
               pole_loc = fromNS.pole_loc_ ;
            }
         }
      }         
   }

   // Create structure to input in the solver system
   struct HeatTransfer2System inputDataHE;
   inputDataHE.is_solids_ = is_solids ; 
   inputDataHE.is_stressCal_ = is_stressCal ;
   inputDataHE.Npart_ = Npart ;
   inputDataHE.level_set_type_ = level_set_type ;
   inputDataHE.Npoints_ = Npoints ;
   inputDataHE.ar_ = ar ;

   // Build the matrix system
   MAC_ModuleExplorer* se =
	exp->create_subexplorer( 0,"DDS_HeatTransferSystem" ) ;
   GLOBAL_EQ = DDS_HeatTransferSystem::create( this, se, TF, inputDataHE ) ;
   se->destroy() ;


   // Duplicates field for error calculation in case of sinusoidal
   // solution with body term 3.pi^2.sin(pi.x).sin(pi.y).sin(pi.z)
   // on a [0:1]x[0:1]x[0:1] domain
   if ( b_bodyterm )
   {
     const_cast<FV_DomainAndFields*>(fromNS.dom_)->duplicate_field(
   	"temperature", "tf_ds_error" ) ;

     TF_DS_ERROR = fromNS.dom_->discrete_field( "tf_ds_error" ) ;
     TF_DS_ERROR->set_BC_values_modif_status( true ) ;
   }

   // Timing routines
   if ( my_rank == is_master )
   {
     SCT_insert_app("Matrix_Assembly&Initialization");
     SCT_insert_app("Matrix_Solution");
     SCT_insert_app("DS_Solution");
     SCT_insert_app("Solver first step");
     SCT_insert_app("Solver x solution");
     SCT_insert_app("Transfer x solution");
     SCT_insert_app("Solver y solution");
     SCT_insert_app("Transfer y solution");
     SCT_insert_app("Solver z solution");
     SCT_insert_app("Transfer z solution");
     SCT_get_elapsed_time("Objects_Creation");
   }
}




//---------------------------------------------------------------------------
DDS_HeatTransfer:: ~DDS_HeatTransfer( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransfer:: ~DDS_HeatTransfer" ) ;

   free_DDS_subcommunicators() ;

}

//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: do_before_time_stepping( FV_TimeIterator const* t_it,
      	std::string const& basename )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransfer:: do_before_time_stepping" ) ;

   if ( my_rank == is_master ) SCT_set_start("Matrix_Assembly&Initialization");

//   FV_OneStepIteration::do_before_time_stepping( t_it, basename ) ;

   allocate_mpi_variables();

   // Initialize temperature vector at the matrix level
   GLOBAL_EQ->initialize_temperature();

   // Necessary especially for cases with non-zero field initiallization
   if (b_restart == false) ugradu_initialization ( );

   // Generate solid particles if required
   if (is_solids) {
      Solids_generation();
      node_property_calculation();
      nodes_temperature_initialization(0);
      nodes_temperature_initialization(1);
      nodes_temperature_initialization(3);
      if (dim == 3) nodes_temperature_initialization(4);
      if (is_stressCal) {
         // Generate discretization of surface in approximate equal area
         generate_surface_discretization ();
      }     
   }

   // Assemble 1D tridiagonal matrices and schur complement calculation
   assemble_temperature_and_schur(t_it);

   if ( my_rank == is_master ) SCT_get_elapsed_time( "Matrix_Assembly&Initialization" );

}

//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: do_before_inner_iterations_stage(
	FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransfer:: do_before_inner_iterations_stage" ) ;

//   FV_OneStepIteration::do_before_inner_iterations_stage( t_it ) ;

   if ((is_par_motion) && (is_solids)) {
      Solids_generation();
      node_property_calculation();
      nodes_temperature_initialization(0);
      nodes_temperature_initialization(1);
      nodes_temperature_initialization(3);
      if (dim == 3) nodes_temperature_initialization(4);

      // Assemble 1D tridiagonal matrices and schur complement calculation
      assemble_temperature_and_schur(t_it);
   }

   // Perform matrix level operations before each time step
   GLOBAL_EQ->at_each_time_step( );
}

//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: do_one_inner_iteration( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransfer:: do_one_inner_iteration" ) ;

   if ( my_rank == is_master ) SCT_set_start("DS_Solution");

   // Solve heat equation using direction splitting 
   HeatEquation_DirectionSplittingSolver(t_it);

   if ( my_rank == is_master ) SCT_get_elapsed_time( "DS_Solution" );

}

//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: do_after_inner_iterations_stage(
	FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransfer:: do_after_inner_iterations_stage" ) ;

//   FV_OneStepIteration::do_after_inner_iterations_stage( t_it ) ;

   if (is_stressCal) {
      compute_fluid_particle_interaction(t_it);
   }

   // Compute temperature change over the time step
   double temperature_time_change = GLOBAL_EQ->compute_temperature_change()
   	/ t_it->time_step() ;
   if ( my_rank == is_master )
     cout << "Temperature change = " <<
     	MAC::doubleToString( ios::scientific, 5, temperature_time_change )
	<< endl;

}

//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: do_after_time_stepping( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransfer:: do_after_time_stepping" ) ;

//   write_output_field();

   // Elapsed time by sub-problems
   if ( my_rank == is_master )
   {
     double cputime = CT_get_elapsed_time();
     cout << endl << "Full problem" << endl;
     write_elapsed_time_smhd(cout,cputime,"Computation time");
     SCT_get_summary(cout,cputime);
   }
//   DS_error_with_analytical_solution(TF,TF_DS_ERROR);
   GLOBAL_EQ->display_debug();
   output_l2norm();

   deallocate_mpi_variables();
}

//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: do_additional_savings( FV_TimeIterator const* t_it,
      	int const& cycleNumber )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransfer:: do_additional_savings" ) ;

}




//----------------------------------------------------------------------
void
DDS_HeatTransfer::write_output_field()
//----------------------------------------------------------------------
{

  ofstream outputFile ;

  std::ostringstream os2;
  os2 << "./DS_results/outputT_" << my_rank << ".csv";
  std::string filename = os2.str();
  outputFile.open(filename.c_str());


//  outputFile.open("./DS_results/output_" << my_rank << ".txt", std::ios_base::out | std::ios_base::trunc) ;
  size_t i,j,k;
  outputFile << "x,y,z,par_ID,void_frac,left,lv,right,rv,bottom,bov,top,tv,behind,bev,front,fv" << endl;
//  outputFile << "x,y,z,par_ID,void_frac,ugradu,error" << endl;

  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);

  NodeProp node = GLOBAL_EQ->get_node_property();
  BoundaryBisec* b_intersect = GLOBAL_EQ->get_b_intersect(0);

  for (size_t comp=0;comp<nb_comps;comp++) {
     // Get local min and max indices
     for (size_t l=0;l<dim;++l) {
        min_unknown_index(l) = TF->get_min_index_unknown_handled_by_proc( comp, l );// - 1;
        max_unknown_index(l) = TF->get_max_index_unknown_handled_by_proc( comp, l );// + 1;
     }

     size_t local_min_k = 0;
     size_t local_max_k = 0;

     if (dim == 3) {
        local_min_k = min_unknown_index(2);
        local_max_k = max_unknown_index(2);
     }

     for (i=min_unknown_index(0);i<=max_unknown_index(0);++i) {
        double xC = TF->get_DOF_coordinate( i, comp, 0 ) ;
        for (j=min_unknown_index(1);j<=max_unknown_index(1);++j) {
           double yC = TF->get_DOF_coordinate( j, comp, 1 ) ;
           for (k=local_min_k;k<=local_max_k;++k) {
              double zC = 0.;
              if (dim == 3) zC = TF->get_DOF_coordinate( k, comp, 2 ) ;
              size_t p = return_node_index(TF,comp,i,j,k);
              size_t id = (size_t) node.parID[comp]->item(p);
              double voidf = node.void_frac[comp]->item(p);

              outputFile << xC << "," << yC << "," << zC << "," << id << "," << voidf;
              for (size_t dir = 0; dir < dim; dir++) {
                  for (size_t off = 0; off < 2; off++) {
                      outputFile << "," << b_intersect[dir].offset[comp]->item(p,off) << "," << b_intersect[dir].value[comp]->item(p,off); 
//                      outputFile << "," << compute_adv_component(comp,i,j,k) << "," << TF->DOF_value(i,j,k,comp,0)*divergence_of_U(comp,i,j,k,0); 
                  }
              }
              outputFile << endl;
           } 
        }
     }
  }
  outputFile.close();
}

//---------------------------------------------------------------------------
double
DDS_HeatTransfer:: bodyterm_value ( double const& xC, double const& yC, double const& zC) 
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
DDS_HeatTransfer:: assemble_DS_un_at_rhs (
        FV_TimeIterator const* t_it, double const& gamma)
//---------------------------------------------------------------------------
{
  double dxC, dyC, dzC, xC, yC, zC=0.;
  double xvalue=0.,yvalue=0.,zvalue=0.,rhs=0., bodyterm=0., adv_value=0.;

  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);

  NodeProp node = GLOBAL_EQ->get_node_property();

  for (size_t comp=0;comp<nb_comps;comp++) {
     // Get local min and max indices
     for (size_t l=0;l<dim;++l) {
        min_unknown_index(l) = TF->get_min_index_unknown_handled_by_proc( comp, l ) ;
        max_unknown_index(l) = TF->get_max_index_unknown_handled_by_proc( comp, l ) ;
     }

     size_t i, j, k;
     for (i=min_unknown_index(0);i<=max_unknown_index(0);++i) {
        // Compute VEC_rhs_x = rhs in x
        dxC = TF->get_cell_size( i, comp, 0 ) ;
        xC = TF->get_DOF_coordinate( i, comp, 0 ) ;
        for (j=min_unknown_index(1);j<=max_unknown_index(1);++j) {
           dyC = TF->get_cell_size( j, comp, 1 ) ;
	   yC = TF->get_DOF_coordinate( j, comp, 1 ) ;
           if (dim ==2 ) {
              k = 0;
              // Dxx for un
              xvalue = compute_un_component(comp,i,j,k,0,3);
              // Dyy for un
              yvalue = compute_un_component(comp,i,j,k,1,1);
	      // Bodyterm for rhs
	      bodyterm = bodyterm_value(xC,yC,zC);
              // Advection term
              adv_value = compute_adv_component(comp,i,j,k);

              if (is_solids) {
                 size_t p = return_node_index(TF,comp,i,j,k);
                 if (node.void_frac[comp]->item(p) == 1) {
                    adv_value = 0.;
                 }
              } 

              rhs = gamma*(xvalue*dyC + yvalue*dxC) - adv_value + (TF->DOF_value( i, j, k, comp, 1 )*dxC*dyC)/(t_it -> time_step());
              TF->set_DOF_value( i, j, k, comp, 0, rhs*(t_it -> time_step())/(dxC*dyC) + gamma*bodyterm*(t_it -> time_step()));

           } else {
              for (k=min_unknown_index(2);k<=max_unknown_index(2);++k) {
                 dzC = TF->get_cell_size( k, comp, 2 ) ;
	         zC = TF->get_DOF_coordinate( k, comp, 2 ) ;
                 // Dxx for un
                 xvalue = compute_un_component(comp,i,j,k,0,3);
                 // Dyy for un
                 yvalue = compute_un_component(comp,i,j,k,1,4);
                 // Dzz for un
                 zvalue = compute_un_component(comp,i,j,k,2,1);
	         // Bodyterm for rhs
	         bodyterm = bodyterm_value(xC,yC,zC);
                 // Advection term
                 adv_value = compute_adv_component(comp,i,j,k);

                 if (is_solids) {
                    size_t p = return_node_index(TF,comp,i,j,k);
                    if (node.void_frac[comp]->item(p) == 1) {
                       adv_value = 0.;
                    }
                 } 

                 rhs = gamma*(xvalue*dyC*dzC + yvalue*dxC*dzC + zvalue*dxC*dyC) - adv_value
                               + (TF->DOF_value( i, j, k, comp, 1 )*dxC*dyC*dzC)/(t_it -> time_step());
                 TF->set_DOF_value( i, j, k, comp, 0, rhs*(t_it -> time_step())/(dxC*dyC*dzC) + gamma*bodyterm*(t_it -> time_step()));
              }
           }
        }
     }
  }
}

//---------------------------------------------------------------------------
double
DDS_HeatTransfer:: divergence_of_U ( size_t const& comp, size_t const& i, size_t const& j, size_t const& k, size_t const& level)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_HeatTransfer:: divergence_of_U" ) ;

   FV_SHIFT_TRIPLET shift = TF->shift_staggeredToCentered() ;

   double xvalue = 0.,yvalue=0.,zvalue=0.,value=0.;

   double dxC = TF->get_cell_size( i, comp, 0 ) ;    
   double dyC = TF->get_cell_size( j, comp, 1 ) ;    
   double dzC = 0.;

   // du/dx
   double xhr= UF->get_DOF_coordinate( shift.i+i,0, 0 ) - UF->get_DOF_coordinate( shift.i+i-1, 0, 0 ) ;
   double xright = UF->DOF_value( shift.i+i, j, k, 0, level ) - UF->DOF_value( shift.i+i-1, j, k, 0, level ) ;

   xvalue = xright/xhr;

   // dv/dy
   double yhr= UF->get_DOF_coordinate( shift.j+j,1, 1 ) - UF->get_DOF_coordinate( shift.j+j-1, 1, 1 ) ;
   double yright = UF->DOF_value( i, shift.j+j, k, 1, level ) - UF->DOF_value( i, shift.j+j-1, k, 1, level ) ;

   yvalue = yright/yhr;

   if (dim == 3) {
      // dw/dz
      dzC = TF->get_cell_size( k, comp, 2 ) ;    
      double zhr= UF->get_DOF_coordinate( shift.k+k,2, 2 ) - UF->get_DOF_coordinate( shift.k+k-1, 2, 2 ) ;
      double zright = UF->DOF_value( i, j, shift.k+k, 2, level ) - UF->DOF_value( i, j, shift.k+k-1, 2, level ) ;

      zvalue = zright/zhr;
   }

   value = (dim == 2) ? (xvalue + yvalue)*dxC*dyC :
                        (xvalue + yvalue + zvalue)*dxC*dyC*dzC ;

   return value;
}

//---------------------------------------------------------------------------
double
DDS_HeatTransfer:: compute_un_component ( size_t const& comp, size_t const& i, size_t const& j, size_t const& k, size_t const& dir, size_t const& level)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_HeatTransfer:: compute_un_component" ) ;

   double xhr,xhl,xright,xleft,yhr,yhl,yright,yleft;
   double zhr,zhl,zright,zleft, value=0.;

   BoundaryBisec* b_intersect = GLOBAL_EQ->get_b_intersect(0);
   NodeProp node = GLOBAL_EQ->get_node_property();

   if (dir == 0) {
      xhr= TF->get_DOF_coordinate( i+1,comp, 0 ) - TF->get_DOF_coordinate( i, comp, 0 ) ;
      xhl= TF->get_DOF_coordinate( i, comp, 0 ) - TF->get_DOF_coordinate( i-1, comp, 0 ) ;
      xright = TF->DOF_value( i+1, j, k, comp, level ) - TF->DOF_value( i, j, k, comp, level ) ;
      xleft = TF->DOF_value( i, j, k, comp, level ) - TF->DOF_value( i-1, j, k, comp, level ) ;

      if (is_solids) {
         size_t p = return_node_index(TF,comp,i,j,k); 
         if (node.void_frac[comp]->item(p) == 0) {
            if ((b_intersect[dir].offset[comp]->item(p,0) == 1)) {
               xleft = TF->DOF_value( i, j, k, comp, level ) - b_intersect[dir].field[comp]->item(p,0);
               xhl = b_intersect[dir].value[comp]->item(p,0);
            } 
            if ((b_intersect[dir].offset[comp]->item(p,1) == 1)) {
               xright = b_intersect[dir].field[comp]->item(p,1) - TF->DOF_value( i, j, k, comp, level );
               xhr = b_intersect[dir].value[comp]->item(p,1);
            }
         } else {
            xright = 0.; xleft = 0.;
         }
      }

      //xvalue = xright/xhr - xleft/xhl;
      if (TF->DOF_in_domain( (int)i-1, (int)j, (int)k, comp) && TF->DOF_in_domain( (int)i+1, (int)j, (int)k, comp))
         value = xright/xhr - xleft/xhl;
      else if (TF->DOF_in_domain( (int)i-1, (int)j, (int)k, comp))
         value = - xleft/xhl;
      else
         value = xright/xhr;
   } else if (dir == 1) {
      yhr= TF->get_DOF_coordinate( j+1,comp, 1 ) - TF->get_DOF_coordinate( j, comp, 1 ) ;
      yhl= TF->get_DOF_coordinate( j, comp, 1 ) - TF->get_DOF_coordinate( j-1, comp, 1 ) ;
      yright = TF->DOF_value( i, j+1, k, comp, level ) - TF->DOF_value( i, j, k, comp, level ) ;
      yleft = TF->DOF_value( i, j, k, comp, level ) - TF->DOF_value( i, j-1, k, comp, level ) ;

      if (is_solids) {
         size_t p = return_node_index(TF,comp,i,j,k);           
         if (node.void_frac[comp]->item(p) == 0) {
            if ((b_intersect[dir].offset[comp]->item(p,0) == 1)) {
               yleft = TF->DOF_value( i, j, k, comp, level ) - b_intersect[dir].field[comp]->item(p,0);
               yhl = b_intersect[dir].value[comp]->item(p,0);
            } 
            if ((b_intersect[dir].offset[comp]->item(p,1) == 1)) {
               yright = b_intersect[dir].field[comp]->item(p,1) - TF->DOF_value( i, j, k, comp, level );
               yhr = b_intersect[dir].value[comp]->item(p,1);
            }
         } else {
            yleft = 0.; yright = 0.;
         }
      }

      //yvalue = yright/yhr - yleft/yhl;
      if (TF->DOF_in_domain((int)i, (int)j-1, (int)k, comp) && TF->DOF_in_domain((int)i, (int)j+1, (int)k, comp))
         value = yright/yhr - yleft/yhl;
      else if(TF->DOF_in_domain((int)i, (int)j-1, (int)k, comp))
         value = - yleft/yhl;
      else
         value = yright/yhr;
   } else if (dir == 2) {
      zhr= TF->get_DOF_coordinate( k+1,comp, 2 ) - TF->get_DOF_coordinate( k, comp, 2 ) ;
      zhl= TF->get_DOF_coordinate( k, comp, 2 ) - TF->get_DOF_coordinate( k-1, comp, 2 ) ;
      zright = TF->DOF_value( i, j, k+1, comp, level ) - TF->DOF_value( i, j, k, comp, level ) ;
      zleft = TF->DOF_value( i, j, k, comp, level ) - TF->DOF_value( i, j, k-1, comp, level ) ;

      if (is_solids) {
         size_t p = return_node_index(TF,comp,i,j,k);
         if (node.void_frac[comp]->item(p) == 0) {
            if ((b_intersect[dir].offset[comp]->item(p,0) == 1)) {
               zleft = TF->DOF_value( i, j, k, comp, level ) - b_intersect[dir].field[comp]->item(p,0);
               zhl = b_intersect[dir].value[comp]->item(p,0);
            }
            if ((b_intersect[dir].offset[comp]->item(p,1) == 1)) {
               zright = b_intersect[dir].field[comp]->item(p,1) - TF->DOF_value( i, j, k, comp, level );
               zhr = b_intersect[dir].value[comp]->item(p,1);
            }
         } else {
            zleft = 0.; zright = 0.;
         }
      }
      
      //zvalue = zright/zhr - zleft/zhl;
      if (TF->DOF_in_domain((int)i, (int)j, (int)k-1, comp) && TF->DOF_in_domain((int)i, (int)j, (int)k+1, comp))
         value = zright/zhr - zleft/zhl;
      else if(TF->DOF_in_domain((int)i, (int)j, (int)k-1, comp))
         value = - zleft/zhl;
      else
         value = zright/zhr;
   }

   return(value);
	   
}

//---------------------------------------------------------------------------
double
DDS_HeatTransfer:: compute_adv_component ( size_t const& comp, size_t const& i, size_t const& j, size_t const& k)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_HeatTransfer:: compute_adv_component" ) ;
   double ugradu = 0., value = 0.;

   if ( AdvectionScheme == "TVD" ) {
      ugradu = assemble_advection_TVD(UF,1,1.,i,j,k,1) - TF->DOF_value(i,j,k,comp,1)*divergence_of_U(comp,i,j,k,1);
   } else if ( AdvectionScheme == "Upwind" ) {
      ugradu = assemble_advection_Upwind(UF,1,1.,i,j,k,1) - TF->DOF_value(i,j,k,comp,1)*divergence_of_U(comp,i,j,k,1);
   } else if ( AdvectionScheme == "Centered" ) {
      ugradu = assemble_advection_Centered(UF,1,1.,i,j,k,1) - TF->DOF_value(i,j,k,comp,1)*divergence_of_U(comp,i,j,k,1);
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
size_t
DDS_HeatTransfer:: return_node_index (
  FV_DiscreteField const* FF,
  size_t const& comp,
  size_t const& i,
  size_t const& j,
  size_t const& k )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransfer:: return_node_index" ) ;

   // Get local min and max indices
   size_t_vector min_index(dim,0);
   size_t_vector max_index(dim,0);
   size_t_vector i_length(dim,0);
   for (size_t l=0;l<dim;++l) {
      // To include knowns at dirichlet boundary in the indexing as well, wherever required pow(2,64)
//     min_unknown_index(l) = FF->get_min_index_unknown_on_proc( comp, l ) - 1;
//     max_unknown_index(l) = FF->get_max_index_unknown_on_proc( comp, l ) + 1;
      min_index(l) = 0;
      max_index(l) = FF->get_local_nb_dof( comp, l );
      i_length(l) = 1 + max_index(l) - min_index(l);
   }

   size_t local_min_k = 0;
   if (dim == 3) local_min_k = min_index(2);

   size_t p = (j-min_index(1)) + i_length(1)*(i-min_index(0)) + i_length(0)*i_length(1)*(k-local_min_k);

   return(p);

}

//---------------------------------------------------------------------------
size_t
DDS_HeatTransfer:: return_row_index (
  FV_DiscreteField const* FF,
  size_t const& comp,
  size_t const& dir,
  size_t const& j,
  size_t const& k )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransfer:: return_row_index" ) ;

   // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l) {
	min_unknown_index(l) = FF->get_min_index_unknown_handled_by_proc( comp, l ) ;
	max_unknown_index(l) = FF->get_max_index_unknown_handled_by_proc( comp, l ) ;
   }

   size_t p=0;

   if (dim == 2) {
      if (dir == 0) {
         p = j - min_unknown_index(1);
      } else if (dir == 1) {
         p = j - min_unknown_index(0);
      }
   } else if (dim == 3) {
      if (dir == 0) {
         p = (j-min_unknown_index(1))+(1+max_unknown_index(1)-min_unknown_index(1))*(k-min_unknown_index(2));
      } else if (dir == 1) {
         p = (j-min_unknown_index(0))+(1+max_unknown_index(0)-min_unknown_index(0))*(k-min_unknown_index(2));
      } else if (dir == 2) {
         p = (j-min_unknown_index(0))+(1+max_unknown_index(0)-min_unknown_index(0))*(k-min_unknown_index(1));
      }
   }

   return(p); 
}

//---------------------------------------------------------------------------
double
DDS_HeatTransfer:: assemble_temperature_matrix (
  FV_DiscreteField const* FF,
  FV_TimeIterator const* t_it,
  double const& gamma,
  size_t const& comp,
  size_t const& dir,
  size_t const& j,
  size_t const& k,
  size_t const& r_index )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransfer:: assemble_temperature_matrix" ) ;

   // Parameters
   double dxr,dxl,xR,xL,xC,right,left,center = 0. ;

   // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l) {
	min_unknown_index(l) = FF->get_min_index_unknown_handled_by_proc( comp, l ) ;
	max_unknown_index(l) = FF->get_max_index_unknown_handled_by_proc( comp, l ) ;
   }

   // Perform assembling
   size_t m, i;
   TDMatrix* A = GLOBAL_EQ-> get_A();

   double Aee_diagcoef=0.;

   for (m=0,i=min_unknown_index(dir);i<=max_unknown_index(dir);++i,++m) {
       xC = FF->get_DOF_coordinate( i, comp, dir ) ;
       xR = FF->get_DOF_coordinate( i+1, comp, dir ) ;
       xL = FF->get_DOF_coordinate( i-1, comp, dir ) ;

       dxr = xR - xC;
       dxl = xC - xL;

       right = -gamma/(dxr);
       left = -gamma/(dxl);

       if (is_solids) {
          size_t p=0;
          BoundaryBisec* b_intersect = GLOBAL_EQ->get_b_intersect(0);
          NodeProp node = GLOBAL_EQ->get_node_property();
          if (dir == 0) {
             p = return_node_index(FF,comp,i,j,k);
          } else if (dir == 1) {
             p = return_node_index(FF,comp,j,i,k);
          } else if (dir == 2) {
             p = return_node_index(FF,comp,j,k,i);
          }
          // if left node is inside the solid particle
          if ((b_intersect[dir].offset[comp]->item(p,0) == 1)) {
             left = -gamma/(b_intersect[dir].value[comp]->item(p,0));
          }
          // if right node is inside the solid particle
          if ((b_intersect[dir].offset[comp]->item(p,1) == 1)) {
             right = -gamma/(b_intersect[dir].value[comp]->item(p,1));
          }
          // if center node is inside the solid particle
          if (node.void_frac[comp]->item(p) == 1.) {
             left = 0.;
             right = 0.;
          }

          center = -(right+left);

          if ((b_intersect[dir].offset[comp]->item(p,0) == 1)) left = 0.;
          if ((b_intersect[dir].offset[comp]->item(p,1) == 1)) right = 0.;
       } else {
          center = - (right+left);
       }

       bool r_bound = false;
       bool l_bound = false;
       // All the proc will have open right bound, except last proc for non periodic systems
       if ((is_iperiodic[dir] != 1) && (rank_in_i[dir] == nb_ranks_comm_i[dir]-1)) r_bound = true;
       // All the proc will have open left bound, except first proc for non periodic systems
       if ((is_iperiodic[dir] != 1) && (rank_in_i[dir] == 0)) l_bound = true;

       // add unsteady term
       double value;
       size_t k_min, k_max;
       double unsteady_term = (FF->get_cell_size(i,comp,dir))/(t_it->time_step());

       if (dim == 2) {
          k_min = 0; k_max = 0;
       } else {
          k_min = min_unknown_index(2); k_max = max_unknown_index(2);
       }

       // Since, this function is used in all directions; 
       // ii, jj, and kk are used to convert the passed arguments corresponding to correct direction
       int ii=0,jj=0,kk=0;

       // Condition for handling the pressure neumann conditions at wall
       if (i==min_unknown_index(dir) && l_bound) {
          if (dir == 0) {
             ii = (int)i-1; jj = (int)min_unknown_index(1); kk = (int)k_min;
          } else if (dir == 1) {
             ii = (int)min_unknown_index(0); jj = (int)i-1; kk = (int)k_min;
          } else if (dir == 2) {
             ii = (int)min_unknown_index(0); jj = (int)min_unknown_index(1); kk = (int)i-1;
          }
          if (FF->DOF_in_domain(ii,jj,kk,comp) && FF->DOF_has_imposed_Dirichlet_value((size_t)ii,(size_t)jj,(size_t)kk,comp)) {
             // For Dirichlet boundary condition
             value = center;
          } else {
             // For Neumann homogeneous boundary condition
             value = -right;
          }
       } else if (i==max_unknown_index(dir) && r_bound) {
          if (dir == 0) {
             ii = (int)i+1; jj = (int)max_unknown_index(1); kk = (int)k_max;
          } else if (dir == 1) {
             ii = (int)max_unknown_index(0); jj = (int)i+1; kk = (int)k_max;
          } else if (dir == 2) {
             ii = (int)max_unknown_index(0); jj = (int)max_unknown_index(1); kk = (int)i+1;
          }
          if (FF->DOF_in_domain(ii,jj,kk,comp) && FF->DOF_has_imposed_Dirichlet_value((size_t)ii,(size_t)jj,(size_t)kk,comp)) {
             // For Dirichlet boundary condition
             value = center;
          } else {
             // For Neumann homogeneous boundary condition
             value = -left;
          }
       } else {
          value = center;
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
          // For last index, Aee comes from this proc as it is interface unknown wrt this proc
          A[dir].ie[comp][r_index]->set_item(m-1,rank_in_i[dir],left);			
          Aee_diagcoef = value;
          A[dir].ei[comp][r_index]->set_item(rank_in_i[dir],m-1,left);
       }

       // Set Aii_sub_diagonal
       if ((rank_in_i[dir] == nb_ranks_comm_i[dir]-1) && (is_iperiodic[dir] != 1)) {
          if (i > min_unknown_index(dir)) A[dir].ii_sub[comp][r_index]->set_item(m-1,left);
       } else {
          if (i<max_unknown_index(dir)) {
             if (i>min_unknown_index(dir)) {
                A[dir].ii_sub[comp][r_index]->set_item(m-1,left);
	     }
	  }
       }

       // Set Aii_super_diagonal
       if ((rank_in_i[dir] == nb_ranks_comm_i[dir]-1) && (is_iperiodic[dir] != 1)) {
          if (i < max_unknown_index(dir)) A[dir].ii_super[comp][r_index]->set_item(m,right);
       } else {
          if (i < max_unknown_index(dir)-1) {
             A[dir].ii_super[comp][r_index]->set_item(m,right);
          }
       }

       // Set Aii_main_diagonal
       if ((rank_in_i[dir] == nb_ranks_comm_i[dir]-1) && (is_iperiodic[dir] != 1)) {
          A[dir].ii_main[comp][r_index]->set_item(m,value);
       } else {
          if (i<max_unknown_index(dir)) {
             A[dir].ii_main[comp][r_index]->set_item(m,value);
          }
       }
   }  // End of for loop

   GLOBAL_EQ->pre_thomas_treatment(comp,dir,A,r_index);

   return(Aee_diagcoef);
}

//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: assemble_schur_matrix (size_t const& comp, size_t const& dir, double const& Aee_diagcoef, size_t const& r_index)
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransfer:: assemble_schur_matrix" ) ;

   TDMatrix* A = GLOBAL_EQ-> get_A();

   if (nb_ranks_comm_i[dir]>1) {

      ProdMatrix* Ap = GLOBAL_EQ->get_Ap();
      ProdMatrix* Ap_proc0 = GLOBAL_EQ->get_Ap_proc0();

      GLOBAL_EQ->compute_product_matrix(A,Ap,comp,dir,r_index);

      LA_SeqMatrix* product_matrix = Ap[dir].ei_ii_ie[comp];

      size_t nbrow = product_matrix->nb_rows();
      // Create a copy of product matrix to receive matrix, this will eliminate the memory leak issue which caused by "create_copy" command
      for (size_t k=0;k<nbrow;k++) {
         for (size_t j=0;j<nbrow;j++) {
            Ap_proc0[dir].ei_ii_ie[comp]->set_item(k,j,product_matrix->item(k,j));
         }
      }

      LA_SeqMatrix* receive_matrix = Ap_proc0[dir].ei_ii_ie[comp];

      if ( rank_in_i[dir] == 0 ) {
         A[dir].ee[comp][r_index]->set_item(0,0,Aee_diagcoef);
         for (size_t i=1;i<(size_t)nb_ranks_comm_i[dir];++i) {

             // Create the container to receive
             size_t nbrows = product_matrix->nb_rows();
             size_t nb_received_data = nbrows*nbrows+1;
             double * received_data = new double [nb_received_data];

             // Receive the data
             static MPI_Status status ;
             MPI_Recv( received_data, (int)nb_received_data, MPI_DOUBLE, (int) i, 0, DDS_Comm_i[dir], &status ) ;

             // Transfer the received data to the receive matrix
             for (size_t k=0;k<(size_t)nbrows;k++) {
                for (size_t j=0;j<(size_t)nbrows;j++) {
                   // Assemble the global product matrix by adding contributions from all the procs
                   receive_matrix->add_to_item(k,j,received_data[k*(nbrows)+j]);
                }
             }

   	     if (is_iperiodic[dir] == 0) {
                if (i<(size_t)nb_ranks_comm_i[dir]-1) {
                   // Assemble the global Aee matrix
                   // No periodic condition in x. So no fe contribution from last proc
                   A[dir].ee[comp][r_index]->set_item(i,i,received_data[nb_received_data-1]);
                }
             } else{
                // Assemble the global Aee matrix
                // Periodic condition in x. So there is fe contribution from last proc
                A[dir].ee[comp][r_index]->set_item(i,i,received_data[nb_received_data-1]);
             }
	     delete [] received_data;

         }
      } else {
         // Create the packed data container

         size_t nbrows = product_matrix->nb_rows();
         size_t nb_send_data = nbrows*nbrows+1;
         double * packed_data = new double [nb_send_data];

         // Fill the packed data container with Aie
         // Iterator only fetches the values present. Zeros are not fetched.

         for (size_t i=0 ; i<nbrows ; i++ ) {
            for ( size_t j=0 ; j<nbrows ; j++ ) {
               // Packing rule
               // Pack the product matrix into a vector
               packed_data[i*nbrows+j]=product_matrix->item(i,j);
            }
         }

         // Fill the last element of packed data with the diagonal coefficient Aee
         packed_data[nb_send_data-1] = Aee_diagcoef;

         // Send the data
         MPI_Send( packed_data, (int)nb_send_data, MPI_DOUBLE, 0, 0, DDS_Comm_i[dir] ) ;

         delete [] packed_data;

      }

      // Assemble the schlur complement in the master proc

      if (rank_in_i[dir] == 0) {
	 TDMatrix* Schur = GLOBAL_EQ-> get_Schur();
         size_t nb_row = Schur[dir].ii_main[comp][r_index]->nb_rows();
         for (int p = 0; p < (int)nb_row; p++) {
            Schur[dir].ii_main[comp][r_index]->set_item(p,A[dir].ee[comp][r_index]->item(p,p)-receive_matrix->item(p,p));
            if (p < (int)nb_row-1) Schur[dir].ii_super[comp][r_index]->set_item(p,-receive_matrix->item(p,p+1));
            if (p > 0) Schur[dir].ii_sub[comp][r_index]->set_item(p-1,-receive_matrix->item(p,p-1));
	    // In case of periodic and multi-processor, there will be a variant of Tridiagonal matrix instead of normal format
            if (is_iperiodic[dir] == 1) {
               Schur[dir].ie[comp][r_index]->set_item(p,0,-receive_matrix->item(p,nb_row)); 
               Schur[dir].ei[comp][r_index]->set_item(0,p,-receive_matrix->item(nb_row,p)); 
	    }
	 }

	 // Pre-thomas treatment on Schur complement
	 GLOBAL_EQ->pre_thomas_treatment(comp,dir,Schur,r_index);

	 // In case of periodic and multi-processor, there will be a variant of Tridiagonal matrix instead of normal format
	 // So, Schur complement of Schur complement is calculated
	 if (is_iperiodic[dir] == 1) {
            Schur[dir].ee[comp][r_index]->set_item(0,0,A[dir].ee[comp][r_index]->item(nb_row,nb_row)-receive_matrix->item(nb_row,nb_row));

	    ProdMatrix* SchurP = GLOBAL_EQ->get_SchurP();
            GLOBAL_EQ->compute_product_matrix_interior(Schur,SchurP,comp,0,dir,r_index);

            TDMatrix* DoubleSchur = GLOBAL_EQ-> get_DoubleSchur();
            DoubleSchur[dir].ii_main[comp][r_index]->set_item(0,Schur[dir].ee[comp][r_index]->item(0,0)-SchurP[dir].ei_ii_ie[comp]->item(0,0));
	 }

      }
   } else if (is_iperiodic[dir] == 1) {
      // Condition for single processor in any direction with periodic boundary conditions
      ProdMatrix* Ap = GLOBAL_EQ->get_Ap();
      ProdMatrix* Ap_proc0 = GLOBAL_EQ->get_Ap_proc0();
      GLOBAL_EQ->compute_product_matrix(A,Ap,comp,dir,r_index);

      LA_SeqMatrix* product_matrix = Ap[dir].ei_ii_ie[comp];

      size_t nbrow = product_matrix->nb_rows();
      // Create a copy of product matrix to receive matrix, this will eliminate the memory leak issue which caused by "create_copy" command
      for (size_t k=0;k<nbrow;k++) {
         for (size_t j=0;j<nbrow;j++) {
            Ap_proc0[dir].ei_ii_ie[comp]->set_item(k,j,product_matrix->item(k,j));
         }
      }

      LA_SeqMatrix* receive_matrix = Ap_proc0[dir].ei_ii_ie[comp];

      A[dir].ee[comp][r_index]->set_item(0,0,Aee_diagcoef);

      TDMatrix* Schur = GLOBAL_EQ-> get_Schur();
      size_t nb_row = Schur[dir].ii_main[comp][r_index]->nb_rows();
      for (int p = 0; p < (int)nb_row; p++) {
         Schur[dir].ii_main[comp][r_index]->set_item(p,A[dir].ee[comp][r_index]->item(p,p)-receive_matrix->item(p,p));
         if (p < (int)nb_row-1) Schur[dir].ii_super[comp][r_index]->set_item(p,-receive_matrix->item(p,p+1));
         if (p > 0) Schur[dir].ii_sub[comp][r_index]->set_item(p-1,-receive_matrix->item(p,p-1));
      }
      GLOBAL_EQ->pre_thomas_treatment(comp,dir,Schur,r_index);
   }
}

//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: assemble_temperature_and_schur ( FV_TimeIterator const* t_it)
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransfer:: assemble_temperature_and_schur" ) ;

   double gamma = (1.0/2.0)*(thermal_conductivity/rho/heat_capacity);

   size_t_vector min_unknown_index(dim,0);
   size_t_vector max_unknown_index(dim,0);

   // Assemble temperature matrix and schur complement for each component
   for (size_t comp = 0; comp < nb_comps; comp++) {
      // Get local min and max indices
      for (size_t l=0;l<dim;++l) {
         min_unknown_index(l) = TF->get_min_index_unknown_handled_by_proc( comp, l ) ;
         max_unknown_index(l) = TF->get_max_index_unknown_handled_by_proc( comp, l ) ;
      }

      for (size_t dir = 0; dir < dim; dir++) {
         size_t dir_j, dir_k;
         size_t local_min_k = 0;
         size_t local_max_k = 0;

         if (dir == 0) {
            dir_j = 1; dir_k = 2;
         } else if (dir == 1) {
            dir_j = 0; dir_k = 2;
         } else if (dir == 2) {
            dir_j = 0; dir_k = 1;
         }

         if (dim == 3) {
            local_min_k = min_unknown_index(dir_k);
            local_max_k = max_unknown_index(dir_k);
         }

         for (size_t j=min_unknown_index(dir_j);j<=max_unknown_index(dir_j);++j) {
            for (size_t k=local_min_k; k <= local_max_k; ++k) {
               size_t r_index = return_row_index (TF,comp,dir,j,k);
               double Aee_diagcoef = assemble_temperature_matrix (TF,t_it,gamma,comp,dir,j,k,r_index);
               assemble_schur_matrix(comp,dir,Aee_diagcoef,r_index);
            }
         }
      }
   }
}

//---------------------------------------------------------------------------
double
DDS_HeatTransfer:: assemble_local_rhs ( size_t const& j, size_t const& k, double const& gamma, FV_TimeIterator const* t_it, size_t const& comp, size_t const& dir)
//---------------------------------------------------------------------------
{
   // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l) {
     min_unknown_index(l) = TF->get_min_index_unknown_handled_by_proc(comp,l) ;
     max_unknown_index(l) = TF->get_max_index_unknown_handled_by_proc(comp,l) ;
   }

   size_t i,pos;
   int m;

   // Compute VEC_rhs_x = rhs in x
   double dC=0, fe=0.;

   // Vector for fi
   LocalVector* VEC = GLOBAL_EQ->get_VEC();

   for (i=min_unknown_index(dir);i<=max_unknown_index(dir);++i) {
     double value=0.;
     pos = i - min_unknown_index(dir);

     // Get contribution of un
     dC = TF->get_cell_size(i,comp,dir) ;

     // x direction
     if (dir == 0) {
        value = compute_un_component(comp,i,j,k,dir,3);
        if (is_solids) {
           BoundaryBisec* b_intersect = GLOBAL_EQ->get_b_intersect(0);
           size_t p = return_node_index(TF,comp,i,j,k);
           if ((b_intersect[dir].offset[comp]->item(p,0) == 1)) {
              value = value - b_intersect[dir].field[comp]->item(p,0)/b_intersect[dir].value[comp]->item(p,0);
           }
           if ((b_intersect[dir].offset[comp]->item(p,1) == 1)) {
              value = value - b_intersect[dir].field[comp]->item(p,1)/b_intersect[dir].value[comp]->item(p,1);
           }
        }
     // y direction
     } else if (dir == 1) {
        if (dim == 2) {
           value = compute_un_component(comp,j,i,k,dir,1);
           if (is_solids) {
              BoundaryBisec* b_intersect = GLOBAL_EQ->get_b_intersect(0);
              size_t p = return_node_index(TF,comp,j,i,k);
              if ((b_intersect[dir].offset[comp]->item(p,0) == 1)) {
                 value = value - b_intersect[dir].field[comp]->item(p,0)/b_intersect[dir].value[comp]->item(p,0);
              }
              if ((b_intersect[dir].offset[comp]->item(p,1) == 1)) {
                 value = value - b_intersect[dir].field[comp]->item(p,1)/b_intersect[dir].value[comp]->item(p,1);
              }
           }
        } else if (dim == 3) {
           value = compute_un_component(comp,j,i,k,dir,4);
           if (is_solids) {
              BoundaryBisec* b_intersect = GLOBAL_EQ->get_b_intersect(0);
              size_t p = return_node_index(TF,comp,j,i,k);
              if ((b_intersect[dir].offset[comp]->item(p,0) == 1)) {
                 value = value - b_intersect[dir].field[comp]->item(p,0)/b_intersect[dir].value[comp]->item(p,0);
              }
              if ((b_intersect[dir].offset[comp]->item(p,1) == 1)) {
                 value = value - b_intersect[dir].field[comp]->item(p,1)/b_intersect[dir].value[comp]->item(p,1);
              }
           }
        }
     // z direction
     } else if (dir == 2) {
        value = compute_un_component(comp,j,k,i,dir,1);
        if (is_solids) {
           BoundaryBisec* b_intersect = GLOBAL_EQ->get_b_intersect(0);
           size_t p = return_node_index(TF,comp,j,k,i);
           if ((b_intersect[dir].offset[comp]->item(p,0) == 1)) {
              value = value - b_intersect[dir].field[comp]->item(p,0)/b_intersect[dir].value[comp]->item(p,0);
           }
           if ((b_intersect[dir].offset[comp]->item(p,1) == 1)) {
              value = value - b_intersect[dir].field[comp]->item(p,1)/b_intersect[dir].value[comp]->item(p,1);
           }
        }       
     }

     double temp_val=0.;
     if (dir == 0) {
        temp_val = (TF->DOF_value(i,j,k,comp,0)*dC)/(t_it->time_step()) - gamma*value;
     } else if (dir == 1) {
        temp_val = (TF->DOF_value(j,i,k,comp,3)*dC)/(t_it->time_step()) - gamma*value;
     } else if (dir == 2) {
        temp_val = (TF->DOF_value(j,k,i,comp,4)*dC)/(t_it->time_step()) - gamma*value;
     }

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

   // Since, this function is used in all directions; 
   // ii, jj, and kk are used to convert the passed arguments corresponding to correct direction
   size_t ii=0,jj=0,kk=0;


   // Effect of boundary conditions in case of non-periodic direction
   m = int(min_unknown_index(dir)) - 1;

   if (dir == 0) {
      ii = m; jj = j; kk = k;
   } else if (dir == 1) {
      ii = j; jj = m; kk = k;
   } else if (dir == 2) {
      ii = j; jj = k; kk = m;
   }

   if ( TF->DOF_in_domain((int)ii,(int)jj,(int)kk,comp))
      if ( TF->DOF_has_imposed_Dirichlet_value(ii,jj,kk,comp)) {
         double ai = 1/(TF->get_DOF_coordinate(m+1,comp,dir) - TF->get_DOF_coordinate(m,comp,dir));
         double dirichlet_value = TF->DOF_value(ii,jj,kk,comp,1) ;
         VEC[dir].local_T[comp]->add_to_item( 0, + gamma * ai * dirichlet_value );
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
      double ai = 1/(TF->get_DOF_coordinate(m,comp,dir) - TF->get_DOF_coordinate(m-1,comp,dir));
      double dirichlet_value = TF->DOF_value(ii,jj,kk,comp,1) ;
      VEC[dir].local_T[comp]->add_to_item( VEC[dir].local_T[comp]->nb_rows()-1 , + gamma * ai * dirichlet_value );
   }

   return fe;
}

//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: compute_Aei_ui (struct TDMatrix* arr, struct LocalVector* VEC, size_t const& comp, size_t const& dir, size_t const& r_index)
//---------------------------------------------------------------------------
{
   // create a replica of local rhs vector in local solution vector
   for (size_t i=0;i<VEC[dir].local_T[comp]->nb_rows();i++){
      VEC[dir].local_solution_T[comp]->set_item(i,VEC[dir].local_T[comp]->item(i));
   }

   // Solve for ui locally and put it in local solution vector
   GLOBAL_EQ->mod_thomas_algorithm(arr, VEC[dir].local_solution_T[comp], comp, dir,r_index);

   for (size_t i=0;i<VEC[dir].T[comp]->nb_rows();i++){
          VEC[dir].T[comp]->set_item(i,0);
   }

   // Calculate Aei*ui in each proc locally and put it in T vector
   arr[dir].ei[comp][r_index]->multiply_vec_then_add(VEC[dir].local_solution_T[comp],VEC[dir].T[comp]);

}

//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: data_packing ( double const& fe, size_t const& comp, size_t const& dir, size_t const& vec_pos)
//---------------------------------------------------------------------------
{
   LocalVector* VEC = GLOBAL_EQ->get_VEC() ;

   double *packed_data = first_pass[dir].send[comp][rank_in_i[dir]]; 

   if (rank_in_i[dir] == 0) {
      // Check if bc is periodic in x
      // If it is, we need to pack two elements apart from fe
      if(is_iperiodic[dir])
          packed_data[3*vec_pos+0] = VEC[dir].T[comp]->item(nb_ranks_comm_i[dir]-1);
      else
          packed_data[3*vec_pos+0] = 0;

      packed_data[3*vec_pos+1] = VEC[dir].T[comp]->item(rank_in_i[dir]);

   } else if (rank_in_i[dir] == nb_ranks_comm_i[dir]-1) {
      // Check if bc is periodic in x
      // If it is, we need to pack two elements apart from fe
      if(is_iperiodic[dir])
          packed_data[3*vec_pos+1] = VEC[dir].T[comp]->item(rank_in_i[dir]);
      else
          packed_data[3*vec_pos+1]=0;

      packed_data[3*vec_pos+0] = VEC[dir].T[comp]->item(rank_in_i[dir]-1);

   } else {
      packed_data[3*vec_pos+0] = VEC[dir].T[comp]->item(rank_in_i[dir]-1);
      packed_data[3*vec_pos+1] = VEC[dir].T[comp]->item(rank_in_i[dir]);
   }

   packed_data[3*vec_pos+2] = fe; // Send the fe values and 0 for last proc
   
}

//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: unpack_compute_ue_pack(size_t const& comp, size_t const& dir, size_t const& p)
//---------------------------------------------------------------------------
{
   LocalVector* VEC = GLOBAL_EQ->get_VEC() ;

   size_t nb_interface_unknowns = VEC[dir].T[comp]->nb_rows();
   for (size_t i=0;i<nb_interface_unknowns;i++) {
       VEC[dir].T[comp]->set_item(i,0);
       VEC[dir].interface_T[comp]->set_item(i,0);
   }

   // If periodic in x, first proc contributes to last interface unknown
   if (is_iperiodic[dir])
      VEC[dir].T[comp]->set_item(nb_ranks_comm_i[dir]-1,first_pass[dir].send[comp][rank_in_i[dir]][3*p]);
   VEC[dir].T[comp]->set_item(0,first_pass[dir].send[comp][rank_in_i[dir]][3*p+1]);
   VEC[dir].interface_T[comp]->set_item(0,first_pass[dir].send[comp][rank_in_i[dir]][3*p+2]);

   // Vec_temp might contain previous values
   for (size_t i=1;i<(size_t)nb_ranks_comm_i[dir];i++) {
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
   for (size_t i=1;i<(size_t)nb_ranks_comm_i[dir];++i) {
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
DDS_HeatTransfer:: unpack_ue(size_t const& comp, double * received_data, size_t const& dir, size_t const& p)
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
DDS_HeatTransfer:: solve_interface_unknowns ( double const& gamma,  FV_TimeIterator const* t_it, size_t const& comp, size_t const& dir)
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
   size_t local_length_j=0;
   size_t local_min_j=0, local_max_j=0;
   size_t local_min_k=0, local_max_k=0;

   if (dir == 0) {
      local_min_j = min_unknown_index(1);
      local_max_j = max_unknown_index(1);
      if (dim == 3) {
         local_min_k = min_unknown_index(2);
         local_max_k = max_unknown_index(2);
      }
   } else if (dir == 1) {
      local_min_j = min_unknown_index(0);
      local_max_j = max_unknown_index(0);
      if (dim == 3) {
         local_min_k = min_unknown_index(2);
         local_max_k = max_unknown_index(2);
      }
   } else if (dir == 2) {
      local_min_j = min_unknown_index(0);
      local_max_j = max_unknown_index(0);
      local_min_k = min_unknown_index(1);
      local_max_k = max_unknown_index(1);
   }

   local_length_j = (local_max_j-local_min_j+1);

   // Send and receive the data first pass
   if ( rank_in_i[dir] == 0 ) {
      // Receiving data from all the slave procs iff multi processors are used
      if (nb_ranks_comm_i[dir] != 1) {
         for (size_t i=1;i<(size_t)nb_ranks_comm_i[dir];++i) {
            static MPI_Status status;
            MPI_Recv( first_pass[dir].receive[comp][i], (int) first_pass[dir].size[comp], MPI_DOUBLE, (int) i, 0, DDS_Comm_i[dir], &status ) ;
         }
      }

      // Solve system of interface unknowns for each y
      if (dim == 2) {
         size_t k = 0;
         for (size_t j=local_min_j;j<=local_max_j;j++) {

     	    size_t p = j-local_min_j;

	    unpack_compute_ue_pack(comp,dir,p); 

            // Need to have the original rhs function assembled for corrosponding j,k pair
            assemble_local_rhs(j,k,gamma,t_it,comp,dir);

            // Setup RHS = fi - Aie*ue for solving ui
            A[dir].ie[comp][p]->multiply_vec_then_add(VEC[dir].interface_T[comp],VEC[dir].local_T[comp],-1.0,1.0);

            // Solve ui and transfer solution into distributed vector
            GLOBAL_EQ->DS_HeatEquation_solver(j,k,min_unknown_index(dir),comp,dir,p);
         }
      } else {
         for (size_t k=local_min_k;k<=local_max_k;k++) {
            for (size_t j=local_min_j;j<=local_max_j;j++) {

   	       size_t p = (j-local_min_j)+local_length_j*(k-local_min_k);

	       unpack_compute_ue_pack(comp,dir,p); 

               // Need to have the original rhs function assembled for corrosponding j,k pair
               assemble_local_rhs(j,k,gamma,t_it,comp,dir);

               // Setup RHS = fi - Aie*ue for solving ui
               A[dir].ie[comp][p]->multiply_vec_then_add(VEC[dir].interface_T[comp],VEC[dir].local_T[comp],-1.0,1.0);

               // Solve ui and transfer solution into distributed vector
               GLOBAL_EQ->DS_HeatEquation_solver(j,k,min_unknown_index(dir),comp,dir,p);
            }
         }
      }
   } else {
      // Send the packed data to master
      MPI_Send( first_pass[dir].send[comp][rank_in_i[dir]], (int) first_pass[dir].size[comp], MPI_DOUBLE, 0, 0, DDS_Comm_i[dir] ) ;
   }

   // Send the data from master iff multi processor are used
   if (nb_ranks_comm_i[dir] != 1) {
      if ( rank_in_i[dir] == 0 ) {
         for (size_t i=1;i<(size_t)nb_ranks_comm_i[dir];++i) {
            MPI_Send( second_pass[dir].send[comp][i], (int) second_pass[dir].size[comp], MPI_DOUBLE,(int) i, 0, DDS_Comm_i[dir] ) ;
         }
      } else {
         // Create the container to receive the ue
         static MPI_Status status ;
         MPI_Recv( second_pass[dir].receive[comp][rank_in_i[dir]], (int) second_pass[dir].size[comp], MPI_DOUBLE, 0, 0, DDS_Comm_i[dir], &status ) ;

         // Solve the system of equations in each proc

         if (dim == 2) {
            size_t k = 0;
            for (size_t j = local_min_j;j<=local_max_j;j++) {
               size_t p = j-local_min_j;

               unpack_ue(comp,second_pass[dir].receive[comp][rank_in_i[dir]],dir,p);

               // Need to have the original rhs function assembled for corrosponding j,k pair
               assemble_local_rhs(j,k,gamma,t_it,comp,dir);

               // Setup RHS = fi - Aie*xe for solving ui
               A[dir].ie[comp][p]->multiply_vec_then_add(VEC[dir].interface_T[comp],VEC[dir].local_T[comp],-1.0,1.0);

               // Solve ui and transfer solution into distributed vector
               GLOBAL_EQ->DS_HeatEquation_solver(j,k,min_unknown_index(dir),comp,dir,p);
            }
         } else {
            for (size_t k = local_min_k;k<=local_max_k;k++) {
               for (size_t j = local_min_j;j<=local_max_j;j++) {
                  size_t p = (j-local_min_j)+local_length_j*(k-local_min_k);

	          unpack_ue(comp,second_pass[dir].receive[comp][rank_in_i[dir]],dir,p);

                  // Need to have the original rhs function assembled for corrosponding j,k pair
                  assemble_local_rhs(j,k,gamma,t_it,comp,dir);

                  // Setup RHS = fi - Aie*xe for solving ui
                  A[dir].ie[comp][p]->multiply_vec_then_add(VEC[dir].interface_T[comp],VEC[dir].local_T[comp],-1.0,1.0);

                  // Solve ui and transfer solution into distributed vector
                  GLOBAL_EQ->DS_HeatEquation_solver(j,k,min_unknown_index(dir),comp,dir,p);
               }
            }
         }
      }
   }
}

//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: Solve_i_in_jk ( FV_TimeIterator const* t_it, double const& gamma, size_t const& dir_i, size_t const& dir_j, size_t const& dir_k )
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

     // Solve in i
     if ((nb_ranks_comm_i[dir_i]>1)||(is_iperiodic[dir_i] == 1)) {
        for (size_t j=min_unknown_index(dir_j);j<=max_unknown_index(dir_j);++j) {
	   for (size_t k=local_min_k; k <= local_max_k; ++k) {
              size_t r_index = return_row_index (TF,comp,dir_i,j,k);
              // Assemble fi and return fe for each proc locally        
              double fe = assemble_local_rhs(j,k,gamma,t_it,comp,dir_i);
              // Calculate Aei*ui in each proc locally
              compute_Aei_ui(A,VEC,comp,dir_i,r_index);
              // Pack Aei_ui and fe for sending it to master
              data_packing (fe,comp,dir_i,r_index);
	   }
        }
        solve_interface_unknowns (gamma,t_it,comp,dir_i);

     } else if (is_iperiodic[dir_i] == 0) {  // Serial mode with non-periodic condition
        for (size_t j=min_unknown_index(dir_j);j<=max_unknown_index(dir_j);++j) {
           for (size_t k=local_min_k; k <= local_max_k; ++k) {
              size_t r_index = return_row_index (TF,comp,dir_i,j,k);
              assemble_local_rhs(j,k,gamma,t_it,comp,dir_i);
              GLOBAL_EQ->DS_HeatEquation_solver(j,k,min_unknown_index(dir_i),comp,dir_i,r_index);
           }
        }
     }
  }
}

//----------------------------------------------------------------------
void
DDS_HeatTransfer::DS_interface_unknown_solver(LA_SeqVector* interface_rhs, size_t const& comp, size_t const& dir, size_t const& r_index)
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransfer:: DS_interface_unknown_solver" ) ;

   TDMatrix* Schur = GLOBAL_EQ-> get_Schur();

   // Condition for variant of Tridiagonal Schur complement in Perioidic direction with multi-processor 
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
      Schur_VEC[dir].interface_T[comp]->set_item(0,Schur_VEC[dir].interface_T[comp]->item(0)-Schur_VEC[dir].T[comp]->item(0));

      // Calculate S_ue, using Schur complement of Schur complement
      GLOBAL_EQ->mod_thomas_algorithm(DoubleSchur, Schur_VEC[dir].interface_T[comp], comp, dir,r_index);

      // Calculate S_fi-Sie*S_ue
      Schur[dir].ie[comp][r_index]->multiply_vec_then_add(Schur_VEC[dir].interface_T[comp],Schur_VEC[dir].local_T[comp],-1.0,1.0);

      // Calculate S_ui
      GLOBAL_EQ->mod_thomas_algorithm(Schur, Schur_VEC[dir].local_T[comp], comp, dir, r_index);

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
DDS_HeatTransfer:: Solids_generation ()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_HeatTransfer:: Solids_generation" ) ;

  // Convert string to istringstream
  istringstream global_par_info(*particle_information);
  // Import particle information in Heat Transfer solver
  ReadStringofSolids(global_par_info);

} 

//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: ReadStringofSolids( istringstream &is )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NSWithHeatTransfer:: ReadStringofSolids" ) ;

   // Structure of particle input data
   PartInput solid = GLOBAL_EQ->get_solid();

   string line;
   istringstream lineStream;
   string cell;
   int cntr = -1;

   // Ensure that the getline ALWAYS reads from start of the string
   is.clear();
   is.seekg(0);

   // read lines
   while (getline(is, line)) {
      lineStream.clear();
      lineStream.str(line);

      // read first word in a line
      getline(lineStream, cell, '\t');

      if (cell == "P") {
         cntr++;

         // Extracting position
         getline(lineStream, cell, '\t');
         double xp = stod(cell);
         getline(lineStream, cell, '\t');
         double yp = stod(cell);
         getline(lineStream, cell, '\t');
         double zp = stod(cell);
         // Extracting Radius
         getline(lineStream, cell, '\t');
         double Rp = stod(cell);
         // Extracting velocity
         getline(lineStream, cell, '\t');
         double vx = stod(cell);
         getline(lineStream, cell, '\t');
         double vy = stod(cell);
         getline(lineStream, cell, '\t');
         double vz = stod(cell);
         // Extracting rotational velocity
         getline(lineStream, cell, '\t');
         double wx = stod(cell);
         getline(lineStream, cell, '\t');
         double wy = stod(cell);
         getline(lineStream, cell, '\t');
         double wz = stod(cell);
	 // Extracting temperature
         getline(lineStream, cell, '\t');
         double Tp = stod(cell);
         getline(lineStream, cell, '\t');
         double off = stod(cell);
/*
	 cout << "Particle info: " << xp << "," << yp << "," << zp << "," 
		 		   << vx << "," << vy << "," << vz << "," 
		 		   << wx << "," << wy << "," << wz << "," 
		 		   << Rp << "," << Tp << endl; 
*/
	 // Storing the information in particle structure
         for (size_t comp=0;comp<nb_comps;comp++) {
            solid.coord[comp]->set_item(cntr,0,xp);
            solid.coord[comp]->set_item(cntr,1,yp);
            solid.coord[comp]->set_item(cntr,2,zp);
            solid.size[comp]->set_item(cntr,Rp);
            solid.vel[comp]->set_item(cntr,0,vx);
            solid.vel[comp]->set_item(cntr,1,vy);
            solid.vel[comp]->set_item(cntr,2,vz);
            solid.ang_vel[comp]->set_item(cntr,0,wx);
            solid.ang_vel[comp]->set_item(cntr,1,wy);
            solid.ang_vel[comp]->set_item(cntr,2,wz);
            solid.temp[comp]->set_item(cntr,Tp);
            solid.inside[comp]->set_item(cntr,off);
         }
      } else if (cell == "O") {
         // Extracting Orientation Matrix
         getline(lineStream, cell, '\t');
         double txx = stod(cell);
         getline(lineStream, cell, '\t');
         double txy = stod(cell);
         getline(lineStream, cell, '\t');
         double txz = stod(cell);
         getline(lineStream, cell, '\t');
         double tyx = stod(cell);
         getline(lineStream, cell, '\t');
         double tyy = stod(cell);
         getline(lineStream, cell, '\t');
         double tyz = stod(cell);
         getline(lineStream, cell, '\t');
         double tzx = stod(cell);
         getline(lineStream, cell, '\t');
         double tzy = stod(cell);
         getline(lineStream, cell, '\t');
         double tzz = stod(cell);
/*
	 cout << "Orientation info: " << txx << "," << txy << "," << txz << "," 
		 	              << tyx << "," << tyy << "," << tyz << "," 
		 	              << tzx << "," << tzy << "," << tzz << endl; 
*/
         // Storing the information in particle structure
         for (size_t comp=0;comp<nb_comps;comp++) {
            solid.thetap[comp]->set_item(cntr,0,txx);
            solid.thetap[comp]->set_item(cntr,1,txy);
            solid.thetap[comp]->set_item(cntr,2,txz);
            solid.thetap[comp]->set_item(cntr,3,tyx);
            solid.thetap[comp]->set_item(cntr,4,tyy);
            solid.thetap[comp]->set_item(cntr,5,tyz);
            solid.thetap[comp]->set_item(cntr,6,tzx);
            solid.thetap[comp]->set_item(cntr,7,tzy);
            solid.thetap[comp]->set_item(cntr,8,tzz);
         }
      }
   }

   // Make sure the number of particle in GRAINS and FLUID are same
   MAC_ASSERT( (size_t)(cntr+1) == Npart ) ;
}


//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: ugradu_initialization (  )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_HeatTransfer:: ugradu_initialization" ) ;

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
DDS_HeatTransfer:: nodes_temperature_initialization ( size_t const& level )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_HeatTransfer:: Solids_flux_correction" ) ;

  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);

  // Vector for solid presence
  NodeProp node = GLOBAL_EQ->get_node_property();

  for (size_t comp=0;comp<nb_comps;comp++) {
     // Get local min and max indices
     for (size_t l=0;l<dim;++l) {
        if (is_iperiodic[l]) {
           min_unknown_index(l) = TF->get_min_index_unknown_handled_by_proc( comp, l ) - 1;
           max_unknown_index(l) = TF->get_max_index_unknown_handled_by_proc( comp, l ) + 1;
        } else {
           min_unknown_index(l) = TF->get_min_index_unknown_handled_by_proc( comp, l );
           max_unknown_index(l) = TF->get_max_index_unknown_handled_by_proc( comp, l );
        }
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
              size_t p = return_node_index(TF,comp,i,j,k);
              if (node.void_frac[comp]->item(p) == 1.) {
                 size_t par_id = (size_t) node.parID[comp]->item(p);
                 double Tpart = impose_solid_temperature (comp,0,10,i,j,k,0.,par_id);
                 TF->set_DOF_value( i, j, k, comp, level,Tpart);
              }
           }
        }
     }
  }
}

//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: compute_fluid_particle_interaction( FV_TimeIterator const* t_it)
//---------------------------------------------------------------------------
{
  MAC_LABEL("DDS_HeatTransfer:: compute_fluid_particle_interaction" ) ;

  string fileName = "./DS_results/particle_Tstress.csv" ;

  doubleVector temp_force(Npart,0);
  double avg_force = 0.;

  size_t Nmax = 0.;
  if (level_set_type == "Sphere") {
     Nmax = (dim == 2) ? ((size_t) Npoints) : ((size_t) (2*Npoints)) ;
  } else if (level_set_type == "Cube") {
     Nmax = (dim == 2) ? ((size_t) (4*Npoints)) : ((size_t) (6*pow(Npoints,2))) ;
  } else if (level_set_type == "Cylinder") {
     double Npm1 = round(pow(MAC::sqrt(Npoints) - MAC::sqrt(MAC::pi()/ar),2.));
     double Ncyl = (Npoints - Npm1);
     double dh = 1. - MAC::sqrt(Npm1/Npoints);
     double Nr = round(2./dh);
     Nmax = (dim == 3) ? ((size_t)(2*Npoints + Nr*Ncyl)) : 0 ;
  }

  for (size_t parID = 0; parID < Npart; parID++) {
     // Contribution of stress tensor
     compute_temperature_gradient_on_particle(temp_force, parID, Nmax);
     // Gathering information from all procs
     temp_force(parID) = pelCOMM->sum(temp_force(parID)) ;

     if (my_rank == 0) {
	avg_force += temp_force(parID);

        ofstream MyFile( fileName.c_str(), ios::app ) ;
        MyFile << t_it -> time() << "," << parID << "," << Nmax << "," << temp_force(parID) << endl;
        MyFile.close( ) ;
     }
  }

  if (my_rank == 0) cout << "Average dimensional Nusselt number for Np " << Nmax << " : " << avg_force/double(Npart) <<endl;
}

//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: compute_temperature_gradient_on_particle(class doubleVector& force, size_t const& parID, size_t const& Np )
//---------------------------------------------------------------------------
{
  MAC_LABEL("DDS_HeatTransfer:: compute_temperature_gradient_on_particle" ) ;

  if (ViscousStressOrder == "first") {
     first_order_temperature_gradient(force, parID, Np );
  } else if (ViscousStressOrder == "second") {
     second_order_temperature_gradient(force, parID, Np );
  }
}

//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: second_order_temperature_gradient(class doubleVector& force, size_t const& parID, size_t const& Np )
//---------------------------------------------------------------------------
{
  MAC_LABEL("DDS_HeatTransfer:: second_order_temperature_gradient" ) ;

  size_t i0_temp;
  size_t comp = 0;

  // Structure of particle input data
  PartInput solid = GLOBAL_EQ->get_solid();
  // Structure of particle surface input data
  SurfaceDiscretize surface = GLOBAL_EQ->get_surface();

  // comp won't matter as the particle position is independent of comp
  double xp = solid.coord[0]->item(parID,0);
  double yp = solid.coord[0]->item(parID,1);
  double zp = solid.coord[0]->item(parID,2);
  double ri = solid.size[0]->item(parID);
/*  
  ofstream outputFile ;
  std::ostringstream os2;
  os2 << "./DS_results/temp_grad_" << my_rank << "_" << parID << ".csv";
  std::string filename = os2.str();
  outputFile.open(filename.c_str());
//  outputFile << "x,y,z,s_xx,s_yy,s_xy" << endl;
  outputFile << "x,y,z,id" << endl;
*/
  doubleArray2D point(3,3,0);
  doubleVector fini(3,0);
  doubleVector level_set(2,1.);          
  boolVector point_in_domain(2,true);
  size_t_vector in_parID(2,0);         //Store particle ID if level_set becomes negative
  boolArray2D found(dim,3,false);
  size_t_array2D i0(3,3,0);

  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);

  doubleVector Dmin(dim,0);
  doubleVector Dmax(dim,0);
  doubleVector rotated_coord(dim,0);
  doubleVector rotated_normal(dim,0);
  intVector sign(dim,0);

  for (size_t i=0;i<Np;i++) {
     // Get local min and max indices
     // Get min and max coordinates in the current processor
     for (size_t l=0;l<dim;++l) {
         min_unknown_index(l) = TF->get_min_index_unknown_handled_by_proc( comp, l );
         max_unknown_index(l) = TF->get_max_index_unknown_handled_by_proc( comp, l );
         Dmin(l) = TF->primary_grid()->get_min_coordinate_on_current_processor(l);
         Dmax(l) = TF->primary_grid()->get_max_coordinate_on_current_processor(l);
     }

     // Rotating surface points
     rotated_coord(0) = ri*surface.coordinate->item(i,0);
     rotated_coord(1) = ri*surface.coordinate->item(i,1);
     rotated_coord(2) = ri*surface.coordinate->item(i,2);

     rotation_matrix(parID,rotated_coord);

     point(0,0) = xp + rotated_coord(0);
     point(0,1) = yp + rotated_coord(1);
     point(0,2) = zp + rotated_coord(2);

     // Rotating surface normal
     rotated_normal(0) = surface.normal->item(i,0);
     rotated_normal(1) = surface.normal->item(i,1);
     rotated_normal(2) = surface.normal->item(i,2);

     rotation_matrix(parID,rotated_normal);

     for (size_t dir=0;dir<dim;dir++) {
        // PBC on rotated surface points
        if (is_iperiodic[dir]) {
           double isize = TF->primary_grid()->get_main_domain_max_coordinate(dir) - TF->primary_grid()->get_main_domain_min_coordinate(dir);
           double imin = TF->primary_grid()->get_main_domain_min_coordinate(dir);
           point(0,dir) = point(0,dir) - MAC::floor((point(0,dir)-imin)/isize)*isize;
        }
        // Finding the grid indexes next to ghost points
        found(dir,0) = FV_Mesh::between(TF->get_DOF_coordinates_vector(comp,dir), point(0,dir), i0_temp);
        if (found(dir,0) == 1) i0(0,dir) = i0_temp;
     }

     // Accessing the smallest grid size in domain
     double dh = TF->primary_grid()->get_smallest_grid_size();

     bool status = (dim==2) ? ((point(0,0) > Dmin(0)) && (point(0,0) <= Dmax(0)) && (point(0,1) > Dmin(1)) && (point(0,1) <= Dmax(1))) :
                              ((point(0,0) > Dmin(0)) && (point(0,0) <= Dmax(0)) && (point(0,1) > Dmin(1)) && (point(0,1) <= Dmax(1))
                                                                                 && (point(0,2) > Dmin(2)) && (point(0,2) <= Dmax(2)));
     double threshold = pow(loc_thres,0.5)*dh;
     double dfdi=0.;
     size_t major_dir=4;

     if (status) {
        for (size_t dir=0;dir<dim;dir++) {
           sign(dir) = (rotated_normal(dir) > 0.) ? 1 : -1;
	   if (dim == 2) {
              if (MAC::abs(rotated_normal(dir)) == MAC::max(MAC::abs(rotated_normal(0)),MAC::abs(rotated_normal(1)))) 
                 major_dir = dir;
	   } else if (dim == 3) {
	      if (MAC::abs(rotated_normal(dir)) == MAC::max(MAC::abs(rotated_normal(0)),
				                   MAC::max(MAC::abs(rotated_normal(1)),MAC::abs(rotated_normal(2)))))
		 major_dir = dir;
	   }
	}

	// Ghost points in i for the calculation of i-derivative of field
        ghost_points_generation( point, i0, sign(major_dir), major_dir, point_in_domain, rotated_normal);

        // Assuming all ghost points are in fluid
        level_set(0) = 1.; level_set(1) = 1.;

        // Checking all the ghost points in the solid/fluid, and storing the parID if present in solid
        for (size_t m=0;m<Npart;m++) {
           // In Normal direction
           if (level_set(0) > threshold) {
              level_set(0) = level_set_function(m,comp,point(1,0),point(1,1),point(1,2),level_set_type);
              level_set(0) *= solid.inside[comp]->item(m);
              if (level_set(0) < threshold) in_parID(0) = m;
           }
           if (level_set(1) > threshold) {
              level_set(1) = level_set_function(m,comp,point(2,0),point(2,1),point(2,2),level_set_type);
              level_set(1) *= solid.inside[comp]->item(m);
              if (level_set(1) < threshold) in_parID(1) = m;
           }
        }

        // Calculation of field variable on ghost point(0)
        fini(0) = impose_solid_temperature_for_ghost(comp,point(0,0),point(0,1),point(0,2),parID);

        // Calculation of field variable on ghost point(1)
        if (level_set(0) > threshold) {
           fini(1) = third_order_ghost_field_estimate(TF, comp, point(1,0), point(1,1), point(1,2), i0(1,0), i0(1,1), i0(1,2), major_dir, sign,0);
        } else if (level_set(0) <= threshold) {
           fini(1) = impose_solid_temperature_for_ghost(comp,point(1,0),point(1,1),point(1,2),in_parID(0));
        }
        // Calculation of field variable on ghost point(2)
        if (level_set(1) > threshold) {
           fini(2) = third_order_ghost_field_estimate(TF, comp, point(2,0), point(2,1), point(2,2), i0(2,0), i0(2,1), i0(2,2), major_dir, sign,0);
        } else if (level_set(1) <= threshold) {
           fini(2) = impose_solid_temperature_for_ghost(comp,point(2,0),point(2,1),point(2,2),in_parID(1));
        }

        // Derivative
        // Point 1 and 2 in computational domain
        if (point_in_domain(0) && point_in_domain(1)) {
           if ((level_set(0) > threshold) && (level_set(1) > threshold)) {
              double dx1 = pow(pow(point(1,0)-point(0,0),2) + pow(point(1,1)-point(0,1),2) + pow(point(1,2)-point(0,2),2),0.5);
              double dx2 = pow(pow(point(2,0)-point(0,0),2) + pow(point(2,1)-point(0,1),2) + pow(point(2,2)-point(0,2),2),0.5);
              dfdi = ((fini(1) - fini(0))*dx2/dx1 - (fini(2) - fini(0))*dx1/dx2)/(dx2-dx1);
           // Point 1 in fluid and 2 in the solid
           } else if ((level_set(0) > threshold) && (level_set(1) <= threshold)) {
              double dx1 = pow(pow(point(1,0)-point(0,0),2) + pow(point(1,1)-point(0,1),2) + pow(point(1,2)-point(0,2),2),0.5);
              dfdi = (fini(1) - fini(0))/dx1;
           // Point 1 is present in solid 
           } else if (level_set(0) <= threshold) {
              double dx1 = pow(pow(point(1,0)-point(0,0),2) + pow(point(1,1)-point(0,1),2) + pow(point(1,2)-point(0,2),2),0.5);
              dfdi = (fini(1) - fini(0))/dx1;
           }
        // Point 1 in computational domain
        } else if (point_in_domain(0) && !point_in_domain(1)) {
           double dx1 = pow(pow(point(1,0)-point(0,0),2) + pow(point(1,1)-point(0,1),2) + pow(point(1,2)-point(0,2),2),0.5);
           dfdi = (fini(1) - fini(0))/dx1;
        // Particle close to wall
        } else if (!point_in_domain(0)) {
           // Creating a new ghost point in domain boundary or wall
           i0(1,major_dir) = (sign(major_dir) == 1) ? (i0(0,major_dir) + 1*sign(major_dir)) : 
		                                      (i0(0,major_dir) + 0*sign(major_dir)) ;
           point(1,major_dir) = TF->get_DOF_coordinate(i0(1,major_dir), comp, major_dir);
           double t1 = (point(1,major_dir) - point(0,major_dir))/rotated_normal(major_dir);

           for (size_t dir=0; dir<dim; dir++) {
              if (dir != major_dir) {
                 point(1,dir) = point(0,dir) + rotated_normal(dir)*t1;
		 size_t i0_t = 0;
                 bool temp = FV_Mesh::between(TF->get_DOF_coordinates_vector(0,dir), point(1,dir), i0_t);
                 if (temp) i0(1,dir) = i0_t; 
              }
	   }

	   // Estimating the field value on new ghost point
           fini(1) = third_order_ghost_field_estimate(TF,comp,point(1,0),point(1,1),point(1,2),i0(1,0),i0(1,1),i0(1,2),major_dir,sign,0);
           double dx1 = pow(pow(point(1,0)-point(0,0),2) + pow(point(1,1)-point(0,1),2) + pow(point(1,2)-point(0,2),2),0.5);
           dfdi = (fini(1) - fini(0))/dx1;
//           outputFile << point(1,0) << "," << point(1,1) << "," << point(1,2) << "," << fini(1) << endl;
        }

//        outputFile << point(0,0) << "," << point(0,1) << "," << point(0,2) << "," << fini(0) << endl;
//        if (point_in_domain(0)) outputFile << point(1,0) << "," << point(1,1) << "," << point(1,2) << "," << fini(1) << endl;
//        if (point_in_domain(1)) outputFile << point(2,0) << "," << point(2,1) << "," << point(2,2) << "," << fini(2) << endl;
     }

     double scale = (dim == 2) ? (ri/(2.*MAC::pi()*ri)) : (ri*ri/(4.*MAC::pi()*ri*ri));

     // Ref: Keating thesis Pg-85
     // point_coord*(area) --> Component of area in particular direction
     force(parID) = force(parID) + dfdi*(surface.area->item(i)*scale);
  }
//  outputFile.close();
}

//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: ghost_points_generation(class doubleArray2D& point, class size_t_array2D& i0, int const& sign, size_t const& major_dir, class boolVector& point_in_domain, class doubleVector& rotated_vector )
//---------------------------------------------------------------------------
{
  MAC_LABEL("DDS_HeatTransfer:: ghost_points_generation" ) ;

  intVector i0_temp(2,0);
  size_t i0_t;

  // Ghost points in i for the calculation of i-derivative of field
  i0_temp(0) = (sign == 1) ? (int(i0(0,major_dir)) + 1*sign) : (int(i0(0,major_dir)) + 0*sign);
  i0_temp(1) = (sign == 1) ? (int(i0(0,major_dir)) + 2*sign) : (int(i0(0,major_dir)) + 1*sign);

  point(1,major_dir) = TF->get_DOF_coordinate(i0_temp(0), 0, major_dir);
  point(2,major_dir) = TF->get_DOF_coordinate(i0_temp(1), 0, major_dir);

  if (MAC::abs(point(0,major_dir)-point(1,major_dir)) < MAC::abs(point(1,major_dir)-point(2,major_dir))) {
     i0_temp(0) = (sign == 1) ? (int(i0(0,major_dir)) + 2*sign) : (int(i0(0,major_dir)) + 1*sign);
     i0_temp(1) = (sign == 1) ? (int(i0(0,major_dir)) + 3*sign) : (int(i0(0,major_dir)) + 2*sign);

     point(1,major_dir) = TF->get_DOF_coordinate(i0_temp(0), 0, major_dir);
     point(2,major_dir) = TF->get_DOF_coordinate(i0_temp(1), 0, major_dir);
  }

  double t1 = (point(1,major_dir) - point(0,major_dir))/rotated_vector(major_dir);
  double t2 = (point(2,major_dir) - point(0,major_dir))/rotated_vector(major_dir);

  for (size_t dir=0; dir<dim; dir++) {
      if (dir != major_dir) {
         point(1,dir) = point(0,dir) + rotated_vector(dir)*t1;
         point(2,dir) = point(0,dir) + rotated_vector(dir)*t2;
      }
  }

  boolArray2D found(2,dim,true);

  for (size_t dir=0; dir<dim; dir++) {
     if (dir == major_dir) {

       // Checking the points in domain or not	     
       if ((i0_temp(0) < 0) || (i0_temp(0) >= (int)TF->get_local_nb_dof(0,dir))) {
          found(0,dir) = 0;
       }

       if ((i0_temp(1) < 0) || (i0_temp(1) >= (int)TF->get_local_nb_dof(0,dir))) {
          found(1,dir) = 0;
       }

       i0(1,dir) = i0_temp(0);
       i0(2,dir) = i0_temp(1);
     } else {
        found(0,dir) = FV_Mesh::between(TF->get_DOF_coordinates_vector(0,dir), point(1,dir), i0_t);
        if (found(0,dir)) i0(1,dir) = i0_t; 
        found(1,dir) = FV_Mesh::between(TF->get_DOF_coordinates_vector(0,dir), point(2,dir), i0_t);
        if (found(1,dir)) i0(2,dir) = i0_t; 
     }
  }

  // Checking the ghost points in domain or not
  point_in_domain(0) = (dim == 2) ? found(0,0) && found(0,1) :
	                            found(0,0) && found(0,1) && found(0,2) ;

  point_in_domain(1) = (dim == 2) ? found(1,0) && found(1,1) :
	                            found(1,0) && found(1,1) && found(1,2) ;

}

//---------------------------------------------------------------------------
double
DDS_HeatTransfer:: third_order_ghost_field_estimate ( FV_DiscreteField* FF, size_t const& comp, double const& xp, double const& yp, double const& zp, size_t const& ii, size_t const& ji, size_t const& ki, size_t const& ghost_points_dir, class intVector& sign, size_t const& level)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_HeatTransfer:: third_order_ghost_field_estimate" ) ;

// Call respective functions based on 2D or 3D domain

   double result = 0.;

   if (dim == 2) {
      result = (ghost_points_dir == 0) ? quadratic_interpolation2D(FF, comp, xp, yp, zp, ii, ji, ki, 1, sign, level) :
                                         quadratic_interpolation2D(FF, comp, xp, yp, zp, ii, ji, ki, 0, sign, level) ;
   } else if (dim == 3) {
      result = quadratic_interpolation3D(FF, comp, xp, yp, zp, ii, ji, ki, ghost_points_dir, sign, level) ;
   }

   return(result);
}

//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: gen_dir_index_of_secondary_ghost_points ( class size_t_vector& index, class intVector& sign, size_t const& interpol_dir, class size_t_array2D& index_g, class boolVector& point_in_domain, size_t const& comp)
//---------------------------------------------------------------------------
{
  MAC_LABEL("DDS_HeatTransfer:: gen_dir_index_of_secondary_ghost_points" ) ;

  intVector i0_temp(3,0);

  for (size_t dir=0;dir<dim;dir++) {
     index_g(0,dir) = index(dir);
     index_g(1,dir) = index(dir);
     index_g(2,dir) = index(dir);
  }

  if (sign(interpol_dir) > 0.) {
     i0_temp(0) = int(index(interpol_dir));	
     i0_temp(1) = int(index(interpol_dir)) + 1;
     i0_temp(2) = int(index(interpol_dir)) + 2;
  } else if (sign(interpol_dir) <= 0.) {
     i0_temp(0) = int(index(interpol_dir)) - 1;	
     i0_temp(1) = int(index(interpol_dir));
     i0_temp(2) = int(index(interpol_dir)) + 1;
  }

  // Checking the ghost points in domain or not
  if ((i0_temp(0) < 0) || (i0_temp(0) >= (int)TF->get_local_nb_dof(comp,interpol_dir))) {
     point_in_domain(0) = 0;
  } else {
     point_in_domain(0) = 1;
  }

  if ((i0_temp(1) < 0) || (i0_temp(1) >= (int)TF->get_local_nb_dof(comp,interpol_dir))) {
     point_in_domain(1) = 0;
  } else {
     point_in_domain(1) = 1;
  }

  if ((i0_temp(2) < 0) || (i0_temp(2) >= (int)TF->get_local_nb_dof(comp,interpol_dir))) {
     point_in_domain(2) = 0;
  } else {
     point_in_domain(2) = 1;
  }

  index_g(0,interpol_dir) = i0_temp(0);
  index_g(1,interpol_dir) = i0_temp(1);
  index_g(2,interpol_dir) = i0_temp(2);

}
//---------------------------------------------------------------------------
double
DDS_HeatTransfer:: quadratic_interpolation2D ( FV_DiscreteField* FF, size_t const& comp, double const& xp, double const& yp, double const& zp, size_t const& i0, size_t const& j0, size_t const& k0, size_t const& interpol_dir, class intVector& sign, size_t const& level)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_HeatTransfer:: quadratic_interpolation2D" ) ;

// Calculates the field value at the ghost points 
// near the particle boundary using the quadratic interpolation 
// inspired from Johansen 1998;
// xp,yp,zp are the ghost point coordinated; interpol_dir is the direction 
// in which the additional points will be used for quadratic interpolation
/*
   ofstream outputFile ;
   std::ostringstream os2;
   os2 << "./DS_results/temp_drag_extras.csv";
   std::string filename = os2.str();
   outputFile.open(filename.c_str(), ios::app);
*/
   // Node information for field(1)
   NodeProp node = GLOBAL_EQ->get_node_property();
   // Intersection information for field(1) in fluid(0)
   BoundaryBisec* bf_intersect = GLOBAL_EQ->get_b_intersect(0);

   // Directional index of point
   size_t_vector index(3,0);
   boolVector point_in_solid(3,0);
   boolVector point_in_domain(3,1);
   doubleVector xi(3,0.);
   // Directional indexes of ghost points
   size_t_array2D index_g(3,3,0);
   // Local node index of ghost points
   size_t_vector node_index(3,0);
   // Decide which scheme to use
   string scheme = "quadratic";

   xi(0) = xp; xi(1) = yp; xi(2) = zp;
   index(0) = i0; index(1) = j0; index(2) = k0;

   double f0=0.,f1=0.,f2=0.;
   double x0=0.,x1=0.,x2=0.;
   double l0=0.,l1=0.,l2=0.;

   // Creating ghost points for quadratic interpolation
   gen_dir_index_of_secondary_ghost_points(index, sign, interpol_dir, index_g, point_in_domain, comp);

   // Check weather the ghost points are in solid or not; TRUE if they are   
   node_index(0) = return_node_index(FF,comp,index_g(0,0),index_g(0,1),index_g(0,2));
   point_in_solid(0) = node.void_frac[comp]->item(node_index(0));
   node_index(1) = return_node_index(FF,comp,index_g(1,0),index_g(1,1),index_g(1,2));
   point_in_solid(1) = node.void_frac[comp]->item(node_index(1));
   node_index(2) = return_node_index(FF,comp,index_g(2,0),index_g(2,1),index_g(2,2));
   point_in_solid(2) = node.void_frac[comp]->item(node_index(2));

   // Assume all the ghost points in fluid
   x0 = FF->get_DOF_coordinate(index_g(0,interpol_dir), comp, interpol_dir);
   x1 = FF->get_DOF_coordinate(index_g(1,interpol_dir), comp, interpol_dir);
   x2 = FF->get_DOF_coordinate(index_g(2,interpol_dir), comp, interpol_dir);

   // Storing the field values assuming all ghost points in fluid and domain
   f0 = FF->DOF_value( index_g(0,0), index_g(0,1), index_g(0,2), comp, level );
   f1 = FF->DOF_value( index_g(1,0), index_g(1,1), index_g(1,2), comp, level );
   f2 = FF->DOF_value( index_g(2,0), index_g(2,1), index_g(2,2), comp, level );

   // Ghost points corrections
   // All points in domain
   if (point_in_domain(0) && point_in_domain(1) && point_in_domain(2)) {
      // 0 in solid, rest in fluid
      if (point_in_solid(0) && !point_in_solid(1) && !point_in_solid(2)) {
         x0 = FF->get_DOF_coordinate(index_g(1,interpol_dir), comp, interpol_dir) - bf_intersect[interpol_dir].value[comp]->item(node_index(1),0);
         f0 = bf_intersect[interpol_dir].field[comp]->item(node_index(1),0);
      // 2 in solid, rest in fluid
      } else if (!point_in_solid(0) && !point_in_solid(1) && point_in_solid(2)) {
         x2 = FF->get_DOF_coordinate(index_g(1,interpol_dir), comp, interpol_dir) + bf_intersect[interpol_dir].value[comp]->item(node_index(1),1);
         f2 = bf_intersect[interpol_dir].field[comp]->item(node_index(1),1);
      // 0, 2 in solid; 1 in fluid
      } else if (point_in_solid(0) && !point_in_solid(1) && point_in_solid(2)) {
         x0 = FF->get_DOF_coordinate(index_g(1,interpol_dir), comp, interpol_dir) - bf_intersect[interpol_dir].value[comp]->item(node_index(1),0);
         f0 = bf_intersect[interpol_dir].field[comp]->item(node_index(1),0);
         x2 = FF->get_DOF_coordinate(index_g(1,interpol_dir), comp, interpol_dir) + bf_intersect[interpol_dir].value[comp]->item(node_index(1),1);
         f2 = bf_intersect[interpol_dir].field[comp]->item(node_index(1),1);
      // 0, 1 in solid; 2 in fluid
      } else if (point_in_solid(0) && point_in_solid(1) && !point_in_solid(2)) {
         x1 = FF->get_DOF_coordinate(index_g(2,interpol_dir), comp, interpol_dir) - bf_intersect[interpol_dir].value[comp]->item(node_index(2),0);
         f1 = bf_intersect[interpol_dir].field[comp]->item(node_index(2),0);
         scheme = "linear12";
      // 1, 2 in solid; 0 in fluid
      } else if (!point_in_solid(0) && point_in_solid(1) && point_in_solid(2)) {
         x1 = FF->get_DOF_coordinate(index_g(0,interpol_dir), comp, interpol_dir) + bf_intersect[interpol_dir].value[comp]->item(node_index(0),1);
         f1 = bf_intersect[interpol_dir].field[comp]->item(node_index(0),1);
         scheme = "linear01";
      }
   // Point 0 and 1 are in domain, 2 not in domain
   } else if (point_in_domain(0) && point_in_domain(1) && !point_in_domain(2)) {
      scheme = "linear01";
      // 0 in fluid; 1 in solid
      if (!point_in_solid(0) && point_in_solid(1)) {
         x1 = FF->get_DOF_coordinate(index_g(0,interpol_dir), comp, interpol_dir) + bf_intersect[interpol_dir].value[comp]->item(node_index(0),1);
         f1 = bf_intersect[interpol_dir].field[comp]->item(node_index(0),1);
      // 0 in solid, 1 in fluid
      } else if (point_in_solid(0) && !point_in_solid(1)) {
         x0 = FF->get_DOF_coordinate(index_g(1,interpol_dir), comp, interpol_dir) - bf_intersect[interpol_dir].value[comp]->item(node_index(1),0);
         f0 = bf_intersect[interpol_dir].field[comp]->item(node_index(1),0);
      }
   // Point 1 and 2 are in domain, 0 not in domain
   } else if (!point_in_domain(0) && point_in_domain(1) && point_in_domain(2)) {
      scheme = "linear12";
      // 1 in fluid; 2 in solid
      if (!point_in_solid(1) && point_in_solid(2)) {
         x2 = FF->get_DOF_coordinate(index_g(1,interpol_dir), comp, interpol_dir) + bf_intersect[interpol_dir].value[comp]->item(node_index(1),1);
         f2 = bf_intersect[interpol_dir].field[comp]->item(node_index(1),1);
      // 1 in solid, 2 in fluid
      } else if (point_in_solid(1) && !point_in_solid(2)) {
         x1 = FF->get_DOF_coordinate(index_g(2,interpol_dir), comp, interpol_dir) - bf_intersect[interpol_dir].value[comp]->item(node_index(2),0);
         f1 = bf_intersect[interpol_dir].field[comp]->item(node_index(2),0);
      }
   }

/*   
   if (comp == 0) {
      if (interpol_dir == 0) {
         double yt = FF->get_DOF_coordinate(index(1), comp, 1); 
         outputFile << x0 << "," << yt << "," << 0. << "," << f0 << endl;
         outputFile << x1 << "," << yt << "," << 0. << "," << f1 << endl;
         outputFile << x2 << "," << yt << "," << 0. << "," << f2 << endl;
      } else if (interpol_dir == 1) {
         double xt = FF->get_DOF_coordinate(index(0), comp, 0);
         outputFile << xt << "," << x0 << "," << 0. << "," << f0 << endl;
         outputFile << xt << "," << x1 << "," << 0. << "," << f1 << endl;
         outputFile << xt << "," << x2 << "," << 0. << "," << f2 << endl;
      }
   }
*/
   double result = 0.;

   if (scheme == "quadratic") {
      l0 = (xi(interpol_dir) - x1)*(xi(interpol_dir) - x2)/(x0 - x1)/(x0 - x2);
      l1 = (xi(interpol_dir) - x0)*(xi(interpol_dir) - x2)/(x1 - x0)/(x1 - x2);
      l2 = (xi(interpol_dir) - x0)*(xi(interpol_dir) - x1)/(x2 - x0)/(x2 - x1);
      result = f0*l0 + f1*l1 + f2*l2;
   } else if (scheme == "linear01") {
      l0 = (xi(interpol_dir) - x1)/(x0 - x1);
      l1 = (xi(interpol_dir) - x0)/(x1 - x0);
      result = f0*l0 + f1*l1;
   } else if (scheme == "linear12") {
      l1 = (xi(interpol_dir) - x2)/(x1 - x2);
      l2 = (xi(interpol_dir) - x1)/(x2 - x1);
      result = f1*l1 + f2*l2;
   }

   return(result);
}

//---------------------------------------------------------------------------
double
DDS_HeatTransfer:: quadratic_interpolation3D ( FV_DiscreteField* FF, size_t const& comp, double const& xp, double const& yp, double const& zp, size_t const& ii, size_t const& ji, size_t const& ki, size_t const& ghost_points_dir, class intVector& sign, size_t const& level)
//---------------------------------------------------------------------------
{
  MAC_LABEL("DDS_HeatTransfer:: quadratic_interpolation3D" ) ;

  // Structure of particle input data
  PartInput solid = GLOBAL_EQ->get_solid();

  doubleArray2D point(3,3,0);
  // Directional indexes of ghost points
  size_t_vector index(3,0);
  // Directional indexes of secondary ghost points
  size_t_array2D index_g(3,3,0);
  // Coordinates of secondary ghost points
  doubleArray2D coord_g(3,3,0.);
  // level set value of the secondary ghost points
  doubleVector level_set(3,1.);          
  // Store particle ID if level_set becomes negative
  size_t_vector in_parID(3,0);         
  // Presence in solid or not
  boolVector point_in_solid(3,0);
  // Presence in domain or not
  boolVector point_in_domain(3,1);
  // Decide which scheme to use
  string scheme = "quadratic";

  size_t sec_ghost_dir = 0;
  size_t sec_interpol_dir = 0;

  point(0,0) = xp; point(0,1) = yp; point(0,2) = zp;
  index(0) = ii; index(1) = ji; index(2) = ki;

  // Ghost points generated in y and then quadratic interpolation in z will generate the same
  // stencil if ghost points are generated in z and the quadratic interpolation done in y
  if (ghost_points_dir == 0) {
     sec_ghost_dir = 1;
     sec_interpol_dir = 2;
  } else if (ghost_points_dir == 1) {
     sec_ghost_dir = 0;
     sec_interpol_dir = 2;
  } else if (ghost_points_dir == 2) {
     sec_ghost_dir = 0;
     sec_interpol_dir = 1;
  }

  // Generate indexes of secondary ghost points in sec_ghost_dir direction 
  gen_dir_index_of_secondary_ghost_points(index, sign, sec_ghost_dir, index_g, point_in_domain, comp);

  // Assume all secondary ghost points in fluid
  double x0 = FF->get_DOF_coordinate(index_g(0,sec_ghost_dir), comp, sec_ghost_dir);
  double x1 = FF->get_DOF_coordinate(index_g(1,sec_ghost_dir), comp, sec_ghost_dir);
  double x2 = FF->get_DOF_coordinate(index_g(2,sec_ghost_dir), comp, sec_ghost_dir);

  if (sec_ghost_dir == 0) {
     coord_g(0,0) = x0; coord_g(0,1) = point(0,1); coord_g(0,2) = point(0,2);
     coord_g(1,0) = x1; coord_g(1,1) = point(0,1); coord_g(1,2) = point(0,2);
     coord_g(2,0) = x2; coord_g(2,1) = point(0,1); coord_g(2,2) = point(0,2);
  } else if (sec_ghost_dir == 1) {
     coord_g(0,0) = point(0,0); coord_g(0,1) = x0; coord_g(0,2) = point(0,2);
     coord_g(1,0) = point(0,0); coord_g(1,1) = x1; coord_g(1,2) = point(0,2);
     coord_g(2,0) = point(0,0); coord_g(2,1) = x2; coord_g(2,2) = point(0,2);
  } else if (sec_ghost_dir == 2) {
     coord_g(0,0) = point(0,0); coord_g(0,1) = point(0,1); coord_g(0,2) = x0;
     coord_g(1,0) = point(0,0); coord_g(1,1) = point(0,1); coord_g(1,2) = x1;
     coord_g(2,0) = point(0,0); coord_g(2,1) = point(0,1); coord_g(2,2) = x2;
  }

  double dh = FF->primary_grid()->get_smallest_grid_size();

  double threshold = pow(loc_thres,0.5)*dh;
  // Checking the secondary ghost points in the solid/fluid, and storing the parID if present in solid
  for (size_t m=0;m<Npart;m++) {
     // x0
     if (level_set(0) > threshold) {
        level_set(0) = level_set_function(m,comp,coord_g(0,0),coord_g(0,1),coord_g(0,2),level_set_type);
        level_set(0) *= solid.inside[comp]->item(m);
        if (level_set(0) < threshold) { in_parID(0) = m; point_in_solid(0) = 1; }
     }
     // x1
     if (level_set(1) > threshold) {
        level_set(1) = level_set_function(m,comp,coord_g(1,0),coord_g(1,1),coord_g(1,2),level_set_type);
        level_set(1) *= solid.inside[comp]->item(m);
        if (level_set(1) < threshold) { in_parID(1) = m; point_in_solid(1) = 1; }
     }
     // x2
     if (level_set(2) > threshold) {
        level_set(2) = level_set_function(m,comp,coord_g(2,0),coord_g(2,1),coord_g(2,2),level_set_type);
        level_set(2) *= solid.inside[comp]->item(m);
        if (level_set(2) < threshold) { in_parID(2) = m; point_in_solid(2) = 1; }
     }
  }

  double f0 = 0., f1 = 0., f2 = 0., del = 0.;

  // Estimate the field values at the secondary ghost points 
  f0=quadratic_interpolation2D(FF,comp,coord_g(0,0),coord_g(0,1),coord_g(0,2),index_g(0,0),index_g(0,1),index_g(0,2),sec_interpol_dir,sign,level);
  f1=quadratic_interpolation2D(FF,comp,coord_g(1,0),coord_g(1,1),coord_g(1,2),index_g(1,0),index_g(1,1),index_g(1,2),sec_interpol_dir,sign,level);
  f2=quadratic_interpolation2D(FF,comp,coord_g(2,0),coord_g(2,1),coord_g(2,2),index_g(2,0),index_g(2,1),index_g(2,2),sec_interpol_dir,sign,level);

  // Ghost points corrections
  if (point_in_domain(0) && point_in_domain(1) && point_in_domain(2)) {
     // 0 in solid, rest in fluid
     if (point_in_solid(0) && !point_in_solid(1) && !point_in_solid(2)) {
        if (sec_ghost_dir == 0) {
           del = find_intersection_for_ghost(x0, x1, yp, zp, in_parID(0), comp, sec_ghost_dir, dh, 0, 0);    
           f0 = impose_solid_temperature_for_ghost(comp,x1-del,yp,zp,in_parID(0));
        } else if (sec_ghost_dir == 1) {
           del = find_intersection_for_ghost(x0, x1, xp, zp, in_parID(0), comp, sec_ghost_dir, dh, 0, 0);    
           f0 = impose_solid_temperature_for_ghost(comp,xp,x1-del,zp,in_parID(0));
        } else if (sec_ghost_dir == 2) {
           del = find_intersection_for_ghost(x0, x1, xp, yp, in_parID(0), comp, sec_ghost_dir, dh, 0, 0);    
           f0 = impose_solid_temperature_for_ghost(comp,xp,yp,x1-del,in_parID(0));
        }
        x0 = x1 - del;
     // 2 in solid, rest in fluid
     } else if (!point_in_solid(0) && !point_in_solid(1) && point_in_solid(2)) {
        if (sec_ghost_dir == 0) {
           del = find_intersection_for_ghost(x1, x2, yp, zp, in_parID(2), comp, sec_ghost_dir, dh, 0, 1);    
           f2 = impose_solid_temperature_for_ghost(comp,x1+del,yp,zp,in_parID(2));
        } else if (sec_ghost_dir == 1) {
           del = find_intersection_for_ghost(x1, x2, xp, zp, in_parID(2), comp, sec_ghost_dir, dh, 0, 1);    
           f2 = impose_solid_temperature_for_ghost(comp,xp,x1+del,zp,in_parID(2));
        } else if (sec_ghost_dir == 2) {
           del = find_intersection_for_ghost(x1, x2, xp, yp, in_parID(2), comp, sec_ghost_dir, dh, 0, 1);    
           f2 = impose_solid_temperature_for_ghost(comp,xp,yp,x1+del,in_parID(2));
        }
        x2 = x1 + del;
     // 0, 2 in solid; 1 in fluid
     } else if (point_in_solid(0) && !point_in_solid(1) && point_in_solid(2)) {
        if (sec_ghost_dir == 0) {
           del = find_intersection_for_ghost(x0, x1, yp, zp, in_parID(0), comp, sec_ghost_dir, dh, 0, 0);    
           f0 = impose_solid_temperature_for_ghost(comp,x1-del,yp,zp,in_parID(0));
        } else if (sec_ghost_dir == 1) {
           del = find_intersection_for_ghost(x0, x1, xp, zp, in_parID(0), comp, sec_ghost_dir, dh, 0, 0);    
           f0 = impose_solid_temperature_for_ghost(comp,xp,x1-del,zp,in_parID(0));
        } else if (sec_ghost_dir == 2) {
           del = find_intersection_for_ghost(x0, x1, xp, yp, in_parID(0), comp, sec_ghost_dir, dh, 0, 0);    
           f0 = impose_solid_temperature_for_ghost(comp,xp,yp,x1-del,in_parID(0));
        }
        x0 = x1 - del;

        if (sec_ghost_dir == 0) {
           del = find_intersection_for_ghost(x1, x2, yp, zp, in_parID(2), comp, sec_ghost_dir, dh, 0, 1);    
           f2 = impose_solid_temperature_for_ghost(comp,x1+del,yp,zp,in_parID(2));
        } else if (sec_ghost_dir == 1) {
           del = find_intersection_for_ghost(x1, x2, xp, zp, in_parID(2), comp, sec_ghost_dir, dh, 0, 1);    
           f2 = impose_solid_temperature_for_ghost(comp,xp,x1+del,zp,in_parID(2));
        } else if (sec_ghost_dir == 2) {
           del = find_intersection_for_ghost(x1, x2, xp, yp, in_parID(2), comp, sec_ghost_dir, dh, 0, 1);    
           f2 = impose_solid_temperature_for_ghost(comp,xp,yp,x1+del,in_parID(2));
        }
        x2 = x1 + del;

     // 0, 1 in solid; 2 in fluid
     } else if (point_in_solid(0) && point_in_solid(1) && !point_in_solid(2)) {
        if (sec_ghost_dir == 0) {
           del = find_intersection_for_ghost(x1, x2, yp, zp, in_parID(1), comp, sec_ghost_dir, dh, 0, 0);    
           f1 = impose_solid_temperature_for_ghost(comp,x2-del,yp,zp,in_parID(1));
        } else if (sec_ghost_dir == 1) {
           del = find_intersection_for_ghost(x1, x2, xp, zp, in_parID(1), comp, sec_ghost_dir, dh, 0, 0);    
           f1 = impose_solid_temperature_for_ghost(comp,xp,x2-del,zp,in_parID(1));
        } else if (sec_ghost_dir == 2) {
           del = find_intersection_for_ghost(x1, x2, xp, yp, in_parID(1), comp, sec_ghost_dir, dh, 0, 0);    
           f1 = impose_solid_temperature_for_ghost(comp,xp,yp,x2-del,in_parID(1));
        }
        x1 = x2 - del;
        scheme = "linear12";
     // 1, 2 in solid; 0 in fluid
     } else if (!point_in_solid(0) && point_in_solid(1) && point_in_solid(2)) {
        if (sec_ghost_dir == 0) {
           del = find_intersection_for_ghost(x0, x1, yp, zp, in_parID(1), comp, sec_ghost_dir, dh, 0, 1);    
           f1 = impose_solid_temperature_for_ghost(comp,x0+del,yp,zp,in_parID(1));
        } else if (sec_ghost_dir == 1) {
           del = find_intersection_for_ghost(x0, x1, xp, zp, in_parID(1), comp, sec_ghost_dir, dh, 0, 1);    
           f1 = impose_solid_temperature_for_ghost(comp,xp,x0+del,zp,in_parID(1));
        } else if (sec_ghost_dir == 2) {
           del = find_intersection_for_ghost(x0, x1, xp, yp, in_parID(1), comp, sec_ghost_dir, dh, 0, 1);    
           f1 = impose_solid_temperature_for_ghost(comp,xp,yp,x0+del,in_parID(1));
        }
        x1 = x0 + del;
        scheme = "linear01";
     }
  // Point 0 and 1 are in domain, 2 not in domain
  } else if (point_in_domain(0) && point_in_domain(1) && !point_in_domain(2)) {
     scheme = "linear01";
     // 0 in fluid; 1 in solid
     if (!point_in_solid(0) && point_in_solid(1)) {
        if (sec_ghost_dir == 0) {
           del = find_intersection_for_ghost(x0, x1, yp, zp, in_parID(1), comp, sec_ghost_dir, dh, 0, 1);    
           f1 = impose_solid_temperature_for_ghost(comp,x0+del,yp,zp,in_parID(1));
        } else if (sec_ghost_dir == 1) {
           del = find_intersection_for_ghost(x0, x1, xp, zp, in_parID(1), comp, sec_ghost_dir, dh, 0, 1);    
           f1 = impose_solid_temperature_for_ghost(comp,xp,x0+del,zp,in_parID(1));
        } else if (sec_ghost_dir == 2) {
           del = find_intersection_for_ghost(x0, x1, xp, yp, in_parID(1), comp, sec_ghost_dir, dh, 0, 1);    
           f1 = impose_solid_temperature_for_ghost(comp,xp,yp,x0+del,in_parID(1));
        }
        x1 = x0 + del;
     // 0 in solid, 1 in fluid
     } else if (point_in_solid(0) && !point_in_solid(1)) {
        if (sec_ghost_dir == 0) {
           del = find_intersection_for_ghost(x0, x1, yp, zp, in_parID(0), comp, sec_ghost_dir, dh, 0, 0);    
           f0 = impose_solid_temperature_for_ghost(comp,x1-del,yp,zp,in_parID(0));
        } else if (sec_ghost_dir == 1) {
           del = find_intersection_for_ghost(x0, x1, xp, zp, in_parID(0), comp, sec_ghost_dir, dh, 0, 0);    
           f0 = impose_solid_temperature_for_ghost(comp,xp,x1-del,zp,in_parID(0));
        } else if (sec_ghost_dir == 2) {
           del = find_intersection_for_ghost(x0, x1, xp, yp, in_parID(0), comp, sec_ghost_dir, dh, 0, 0);    
           f0 = impose_solid_temperature_for_ghost(comp,xp,yp,x1-del,in_parID(0));
        }
        x0 = x1 - del;
     }
  // Point 1 and 2 are in domain, 0 not in domain
  } else if (!point_in_domain(0) && point_in_domain(1) && point_in_domain(2)) {
     scheme = "linear12";
     // 1 in fluid; 2 in solid
     if (!point_in_solid(1) && point_in_solid(2)) {
        if (sec_ghost_dir == 0) {
           del = find_intersection_for_ghost(x1, x2, yp, zp, in_parID(2), comp, sec_ghost_dir, dh, 0, 1);    
           f2 = impose_solid_temperature_for_ghost(comp,x1+del,yp,zp,in_parID(2));
        } else if (sec_ghost_dir == 1) {
           del = find_intersection_for_ghost(x1, x2, xp, zp, in_parID(2), comp, sec_ghost_dir, dh, 0, 1);    
           f2 = impose_solid_temperature_for_ghost(comp,xp,x1+del,zp,in_parID(2));
        } else if (sec_ghost_dir == 2) {
           del = find_intersection_for_ghost(x1, x2, xp, yp, in_parID(2), comp, sec_ghost_dir, dh, 0, 1);    
           f2 = impose_solid_temperature_for_ghost(comp,xp,yp,x1+del,in_parID(2));
        }
        x2 = x1 + del;
     // 1 in solid, 2 in fluid
     } else if (point_in_solid(1) && !point_in_solid(2)) {
        if (sec_ghost_dir == 0) {
           del = find_intersection_for_ghost(x1, x2, yp, zp, in_parID(1), comp, sec_ghost_dir, dh, 0, 0);    
           f1 = impose_solid_temperature_for_ghost(comp,x2-del,yp,zp,in_parID(1));
        } else if (sec_ghost_dir == 1) {
           del = find_intersection_for_ghost(x1, x2, xp, zp, in_parID(1), comp, sec_ghost_dir, dh, 0, 0);    
           f1 = impose_solid_temperature_for_ghost(comp,xp,x2-del,zp,in_parID(1));
        } else if (sec_ghost_dir == 2) {
           del = find_intersection_for_ghost(x1, x2, xp, yp, in_parID(1), comp, sec_ghost_dir, dh, 0, 0);    
           f1 = impose_solid_temperature_for_ghost(comp,xp,yp,x2-del,in_parID(1));
        }
        x1 = x2 - del;
     }
  }

  double l0 = 0., l1 = 0., l2 = 0.;
  double result = 0.;

  if (scheme == "quadratic") {
     l0 = (point(0,sec_ghost_dir) - x1)*(point(0,sec_ghost_dir) - x2)/(x0 - x1)/(x0 - x2);
     l1 = (point(0,sec_ghost_dir) - x0)*(point(0,sec_ghost_dir) - x2)/(x1 - x0)/(x1 - x2);
     l2 = (point(0,sec_ghost_dir) - x0)*(point(0,sec_ghost_dir) - x1)/(x2 - x0)/(x2 - x1);
     result = f0*l0 + f1*l1 + f2*l2;
  } else if (scheme == "linear01") {
     l0 = (point(0,sec_ghost_dir) - x1)/(x0 - x1);
     l1 = (point(0,sec_ghost_dir) - x0)/(x1 - x0);
     result = f0*l0 + f1*l1;
  } else if (scheme == "linear12") {
     l1 = (point(0,sec_ghost_dir) - x2)/(x1 - x2);
     l2 = (point(0,sec_ghost_dir) - x1)/(x2 - x1);
     result = f1*l1 + f2*l2;
  }

  return(result);
}



//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: first_order_temperature_gradient(class doubleVector& force, size_t const& parID, size_t const& Np )
//---------------------------------------------------------------------------
{
  MAC_LABEL("DDS_HeatTransfer:: first_order_temperature_gradient" ) ;

  size_t i0_temp;
  size_t comp = 0;

  // Structure of particle input data
  PartInput solid = GLOBAL_EQ->get_solid();
  // Structure of particle surface input data
  SurfaceDiscretize surface = GLOBAL_EQ->get_surface();

  // comp won't matter as the particle position is independent of comp
  double xp = solid.coord[comp]->item(parID,0);
  double yp = solid.coord[comp]->item(parID,1);
  double zp = solid.coord[comp]->item(parID,2);
  double ri = solid.size[comp]->item(parID);
/*  
  ofstream outputFile ;
  std::ostringstream os2;
  os2 << "./DS_results/temp_grad_" << my_rank << "_" << parID << ".csv";
  std::string filename = os2.str();
  outputFile.open(filename.c_str());
  outputFile << "x,y,z,id" << endl;*/
//  outputFile << "i,Nu" << endl;

  doubleArray2D ipoint(3,3,0.);         
  doubleVector fini(3,0);
  doubleVector level_set(2,1.);          
  boolVector in_domain(2,true);        //true if ghost point in the computational domain
  size_t_vector in_parID(2,0);         //Store particle ID if level_set becomes negative
  boolArray2D found(3,dim,false);
  size_t_array2D i0(3,dim,0);

  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);

  doubleVector Dmin(dim,0);
  doubleVector Dmax(dim,0);
  doubleVector rotated_coord(dim,0);
  doubleVector rotated_normal(dim,0);

  for (size_t i=0;i<Np;i++) {
     // Get local min and max indices
     for (size_t l=0;l<dim;++l) {
        min_unknown_index(l) = TF->get_min_index_unknown_handled_by_proc( comp, l );
        max_unknown_index(l) = TF->get_max_index_unknown_handled_by_proc( comp, l );
        Dmin(l) = TF->primary_grid()->get_min_coordinate_on_current_processor(l);
        Dmax(l) = TF->primary_grid()->get_max_coordinate_on_current_processor(l);
     }

     // Rotating surface points
     rotated_coord(0) = ri*surface.coordinate->item(i,0);
     rotated_coord(1) = ri*surface.coordinate->item(i,1);
     rotated_coord(2) = ri*surface.coordinate->item(i,2);

     rotation_matrix(parID,rotated_coord);

     ipoint(0,0) = xp + rotated_coord(0);
     ipoint(0,1) = yp + rotated_coord(1);
     ipoint(0,2) = zp + rotated_coord(2);

     // Rotating surface normal
     rotated_normal(0) = surface.normal->item(i,0);
     rotated_normal(1) = surface.normal->item(i,1);
     rotated_normal(2) = surface.normal->item(i,2);

     rotation_matrix(parID,rotated_normal);

     // Correction in case of periodic boundary condition in any direction
     for (size_t dir=0;dir<dim;dir++) {
        // PBC on rotated surface points
        if (is_iperiodic[dir]) {
           double isize = TF->primary_grid()->get_main_domain_max_coordinate(dir) - TF->primary_grid()->get_main_domain_min_coordinate(dir);
           double imin = TF->primary_grid()->get_main_domain_min_coordinate(dir);
           ipoint(0,dir) = ipoint(0,dir) - MAC::floor((ipoint(0,dir)-imin)/isize)*isize;
        }
        // Finding the grid indexes next to ghost points
        found(0,dir) = FV_Mesh::between(TF->get_DOF_coordinates_vector(comp,dir), ipoint(0,dir), i0_temp) ;
        if (found(0,dir) == 1) i0(0,dir) = i0_temp ;
     }

     // Accessing the smallest grid size in domain
     double dh = TF->primary_grid()->get_smallest_grid_size();

     bool status = (dim==2) ? ((ipoint(0,0) > Dmin(0)) && (ipoint(0,0) <= Dmax(0)) && (ipoint(0,1) > Dmin(1)) && (ipoint(0,1) <= Dmax(1))) :
                              ((ipoint(0,0) > Dmin(0)) && (ipoint(0,0) <= Dmax(0)) && (ipoint(0,1) > Dmin(1)) && (ipoint(0,1) <= Dmax(1))
                                                                                   && (ipoint(0,2) > Dmin(2)) && (ipoint(0,2) <= Dmax(2)));

     double threshold = pow(loc_thres,0.5)*dh;
     double dfdi=0.;

     if (status) {
        for (size_t l=0;l<dim;++l) {
           // Ghost points in normal direction at the particle surface
           ipoint(1,l) = ipoint(0,l) + dh*rotated_normal(l);
           ipoint(2,l) = ipoint(0,l) + 2.*dh*rotated_normal(l);

           // Periocid boundary conditions
           if (is_iperiodic[l]) {
              double isize = TF->primary_grid()->get_main_domain_max_coordinate(l) - TF->primary_grid()->get_main_domain_min_coordinate(l);
              double imin = TF->primary_grid()->get_main_domain_min_coordinate(l);
              ipoint(1,l) = ipoint(1,l) - MAC::floor((ipoint(1,l)-imin)/isize)*isize;
              ipoint(2,l) = ipoint(2,l) - MAC::floor((ipoint(2,l)-imin)/isize)*isize;
           }

           // Finding the grid indexes next to ghost points
           found(1,l) = FV_Mesh::between(TF->get_DOF_coordinates_vector(comp,l), ipoint(1,l), i0_temp);
           if (found(1,l) == 1) i0(1,l) = i0_temp;
           found(2,l) = FV_Mesh::between(TF->get_DOF_coordinates_vector(comp,l), ipoint(2,l), i0_temp);
           if (found(2,l) == 1) i0(2,l) = i0_temp;

        }

        // Assuming all ghost points are in fluid
        level_set(0) = 1.; level_set(1) = 1.;

        // Checking all the ghost points in the solid/fluid, and storing the parID if present in solid
        for (size_t m=0;m<Npart;m++) {
           if (level_set(0) > threshold) {
              level_set(0) = level_set_function(m,comp,ipoint(1,0),ipoint(1,1),ipoint(1,2),level_set_type);
              level_set(0) *= solid.inside[comp]->item(m);
              if (level_set(0) < threshold) in_parID(0) = m;
           }
           if (level_set(1) > threshold) {
              level_set(1) = level_set_function(m,comp,ipoint(2,0),ipoint(2,1),ipoint(2,2),level_set_type);
              level_set(1) *= solid.inside[comp]->item(m);
              if (level_set(1) < threshold) in_parID(1) = m;
           }
        }

        // Calculation of field variable on ghost point(0)
        fini(0) = impose_solid_temperature_for_ghost(comp,ipoint(0,0),ipoint(0,1),ipoint(0,2),parID);

        if (dim == 2) {
           in_domain(0) = found(1,0) && found(1,1);
           in_domain(1) = found(2,0) && found(2,1);
           // Calculation of field variable on ghost point(1)
           if ((level_set(0) > threshold) && in_domain(0)) 
              fini(1) = ghost_field_estimate_on_face (TF,comp,i0(1,0),i0(1,1),0, ipoint(1,0), ipoint(1,1),0, dh,2,0);
           // Calculation of field variable on ghost point(2)
           if ((level_set(1) > threshold) && in_domain(1)) 
              fini(2) = ghost_field_estimate_on_face (TF,comp,i0(2,0),i0(2,1),0, ipoint(2,0), ipoint(2,1),0, dh,2,0);

        } else if (dim == 3) {
           in_domain(0) = found(1,0) && found(1,1) && found(1,2);
           in_domain(1) = found(2,0) && found(2,1) && found(2,2);

           // Calculation of field variable on ghost point(1)
           if ((level_set(0) > threshold) && in_domain(0)) 
              fini(1) = ghost_field_estimate_in_box (TF,comp,i0(1,0),i0(1,1),i0(1,2),ipoint(1,0),ipoint(1,1),ipoint(1,2),dh,0,parID);
           // Calculation of field variable on ghost point(2)
           if ((level_set(1) > threshold) && in_domain(1))
              fini(2) = ghost_field_estimate_in_box (TF,comp,i0(2,0),i0(2,1),i0(2,2),ipoint(2,0),ipoint(2,1),ipoint(2,2),dh,0,parID);
        }

        // Derivative
        // Both points 1 and 2 are in fluid, and both in the computational domain
        if ((level_set(0) > threshold) && (level_set(1) > threshold) && (in_domain(0) && in_domain(1))) {
           dfdi = (-fini(2) + 4.*fini(1) - 3.*fini(0))/2./dh;
        // Point 1 in fluid and 2 is either in the solid or out of the computational domain
        } else if ((level_set(0) > threshold) && ((level_set(1) <= threshold) || ((in_domain(1) == 0) && (in_domain(0) == 1)))) {
           dfdi = (fini(1) - fini(0))/dh;
        // Point 1 is present in solid 
        } else if (level_set(0) <= threshold) {
           dfdi = (impose_solid_temperature_for_ghost(comp,ipoint(1,0),ipoint(1,1),ipoint(1,2),in_parID(0)) - fini(0))/dh;
        // Point 1 is out of the computational domain 
        } else if (!in_domain(0)) { 
           intVector sign(dim,0);
	   size_t major_dir = 4;
           for (size_t l=0;l<dim;l++) {
              sign(l) = (rotated_normal(l) > 0.) ? 1 : -1;
              if (dim == 2) {
                 if (MAC::abs(rotated_normal(l)) == MAC::max(MAC::abs(rotated_normal(0)),MAC::abs(rotated_normal(1)))) 
                    major_dir = l;
              } else if (dim == 3) {
	         if (MAC::abs(rotated_normal(l)) == MAC::max(MAC::abs(rotated_normal(0)),
	                                              MAC::max(MAC::abs(rotated_normal(1)),MAC::abs(rotated_normal(2)))))
                    major_dir = l;
	      }
	   }
           // Creating a new ghost point in domain boundary or wall
           i0(1,major_dir) = (sign(major_dir) == 1) ? (i0(0,major_dir) + 1*sign(major_dir)) : 
		                                      (i0(0,major_dir) + 0*sign(major_dir)) ;
           ipoint(1,major_dir) = TF->get_DOF_coordinate(i0(1,major_dir), comp, major_dir);
           double t1 = (ipoint(1,major_dir) - ipoint(0,major_dir))/rotated_normal(major_dir);

           for (size_t dir=0; dir<dim; dir++) {
              if (dir != major_dir) {
                 ipoint(1,dir) = ipoint(0,dir) + rotated_normal(dir)*t1;
		 size_t i0_t = 0;
                 bool temp = FV_Mesh::between(TF->get_DOF_coordinates_vector(0,dir), ipoint(1,dir), i0_t);
                 if (temp) i0(1,dir) = i0_t; 
              }
	   }

	   // Estimating the field value on new ghost point
           fini(1) = third_order_ghost_field_estimate(TF,comp,ipoint(1,0),ipoint(1,1),ipoint(1,2),i0(1,0),i0(1,1),i0(1,2),major_dir,sign,0);
           double dx1 = pow(pow(ipoint(1,0)-ipoint(0,0),2) + pow(ipoint(1,1)-ipoint(0,1),2) + pow(ipoint(1,2)-ipoint(0,2),2),0.5);
           dfdi = (fini(1) - fini(0))/dx1;
        }
/*
        outputFile << ipoint(0,0) << "," << ipoint(0,1) << "," << ipoint(0,2) << "," << fini(0) << endl;
        outputFile << ipoint(1,0) << "," << ipoint(1,1) << "," << ipoint(1,2) << "," << fini(1) << endl;
        outputFile << ipoint(2,0) << "," << ipoint(2,1) << "," << ipoint(2,2) << "," << fini(2) << endl;*/
     }

//     outputFile << i << "," << dfdi << endl;
     double scale = (dim == 2) ? (ri/(2.*MAC::pi()*ri)) : (ri*ri/(4.*MAC::pi()*ri*ri));

     // Ref: Keating thesis Pg-85
     // point_coord*(area) --> Component of area in particular direction
     force(parID) = force(parID) + dfdi*(surface.area->item(i)*scale);
  }
//  outputFile.close();  
}

//---------------------------------------------------------------------------
double
DDS_HeatTransfer:: ghost_field_estimate_on_face ( FV_DiscreteField* FF, size_t const& comp, size_t const& i0, size_t const& j0, size_t const& k0, double const& x0, double const& y0, double const& z0, double const& dh, size_t const& face_vec, size_t const& level)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_HeatTransfer:: ghost_field_estimate_on_face" ) ;

// Calculates the field value on a face at the ghost points 
// near the particle boundary considering boundary affects;
// x0,y0,z0 are the ghost point coordinated; i0,j0,k0 is the
// bottom left index of the face cell; face_vec is the 
// normal vector of the face (i.e. 0 is x,1 is y, 2 is z) 

   BoundaryBisec* bf_intersect = GLOBAL_EQ->get_b_intersect(0);    
   NodeProp node = GLOBAL_EQ->get_node_property();                
   PartInput solid = GLOBAL_EQ->get_solid();

   size_t_array2D p(dim,dim,0);
   // arrays of vertex indexes of the cube/square
   size_t_array2D ix(dim,dim,0);
   size_t_array2D iy(dim,dim,0);
   size_t_array2D iz(dim,dim,0);
   // Min and max of the cell containing ghost point
   doubleArray2D extents(dim,2,0);
   // Field value at vertex of face
   doubleArray2D f(2,2,0);
   // Interpolated values at the walls
   doubleArray2D fwall(2,2,0);
   // Distance of grid/particle walls from the ghost point
   doubleArray2D del_wall(2,2,0);
   // Ghost point coordinate
   doubleVector xghost(3,0);

   // Direction in the plane of face
   size_t dir1=0, dir2=0;

   if (face_vec == 0) {
      ix(0,0) = i0,     iy(0,0) = j0,    iz(0,0) = k0;
      ix(1,0) = i0,     iy(1,0) = j0,    iz(1,0) = k0+1;
      ix(0,1) = i0,     iy(0,1) = j0+1,  iz(0,1) = k0;
      ix(1,1) = i0,     iy(1,1) = j0+1,  iz(1,1) = k0+1;
      dir1 = 2, dir2 = 1;
   } else if (face_vec == 1) {
      ix(0,0) = i0,     iy(0,0) = j0,    iz(0,0) = k0;
      ix(1,0) = i0+1,   iy(1,0) = j0,    iz(1,0) = k0;
      ix(0,1) = i0,     iy(0,1) = j0,    iz(0,1) = k0+1;
      ix(1,1) = i0+1,   iy(1,1) = j0,    iz(1,1) = k0+1;
      dir1 = 0, dir2 = 2;
   } else if (face_vec == 2) {
      ix(0,0) = i0,     iy(0,0) = j0,    iz(0,0) = k0;
      ix(1,0) = i0+1,   iy(1,0) = j0,    iz(1,0) = k0;
      ix(0,1) = i0,     iy(0,1) = j0+1,  iz(0,1) = k0;
      ix(1,1) = i0+1,   iy(1,1) = j0+1,  iz(1,1) = k0;
      dir1 = 0, dir2 = 1;
   }

   xghost(0) = x0;
   xghost(1) = y0;
   xghost(2) = z0;

   for (size_t i=0; i < 2; i++) {
      for (size_t j=0; j < 2; j++) {
         // Face vertex index
         p(i,j) = return_node_index(FF,comp,ix(i,j),iy(i,j),iz(i,j));
         // Vertex field values
         f(i,j) = FF->DOF_value( ix(i,j), iy(i,j), iz(i,j), comp, level ); 
      }
   }

   // Min x-coordinate in the grid cell
   extents(0,0) = FF->get_DOF_coordinate( ix(0,0), comp, 0 ) ;
   // Max x-coordinate in the grid cell
   extents(0,1) = FF->get_DOF_coordinate( ix(1,1), comp, 0 ) ;
   // Min y-coordinate in the grid cell
   extents(1,0) = FF->get_DOF_coordinate( iy(0,0), comp, 1 ) ;
   // Max y-coordinate in the grid cell
   extents(1,1) = FF->get_DOF_coordinate( iy(1,1), comp, 1 ) ;

   if (dim == 3) {
      // Min z-coordinate in the grid cell
      extents(2,0) = FF->get_DOF_coordinate( iz(0,0), comp, 2 ) ;
      // Max z-coordinate in the grid cell
      extents(2,1) = FF->get_DOF_coordinate( iz(1,1), comp, 2 ) ;
   }

   // Contribution from left and right wall
   for (size_t i = 0; i < 2; i++) {     // 0 --> left; 1 --> right
      if ((node.void_frac[comp]->item(p(i,0)) == 0) && (node.void_frac[comp]->item(p(i,1)) == 0)) {
         fwall(0,i) = ((extents(dir2,1) - xghost(dir2))*f(i,0) + (xghost(dir2) - extents(dir2,0))*f(i,1))/(extents(dir2,1)-extents(dir2,0));
         del_wall(0,i) = MAC::abs(extents(dir1,i) - xghost(dir1));
      // if bottom vertex is in fluid domain
      } else if ((node.void_frac[comp]->item(p(i,0)) == 0) && (bf_intersect[dir2].offset[comp]->item(p(i,0),1) == 1)) {
         double yint = bf_intersect[dir2].value[comp]->item(p(i,0),1);
         // Condition where intersection distance is more than ghost point distance, it means that the ghost 
         // point can be projected on the wall
         if (yint >= (xghost(dir2)-extents(dir2,0))) {
            fwall(0,i) = ((extents(dir2,0)+yint-xghost(dir2))*f(i,0)+(xghost(dir2)-extents(dir2,0))*bf_intersect[dir2].field[comp]->item(p(i,0),1))/yint;
            del_wall(0,i) = MAC::abs(extents(dir1,i) - xghost(dir1));
         // Ghost point cannot be projected on the wall, as the solid surface come first
         } else {
            size_t id = (size_t) node.parID[comp]->item(p(i,1));
            if (face_vec > dir2) {
               del_wall(0,i) = (i==1) ? 
                   find_intersection_for_ghost(xghost(dir1), extents(dir1,i), xghost(dir2), xghost(face_vec), id, comp, dir1, dh, level, i) :
                   find_intersection_for_ghost(extents(dir1,i), xghost(dir1), xghost(dir2), xghost(face_vec), id, comp, dir1, dh, level, i) ;
            } else {
               del_wall(0,i) = (i==1) ? 
                   find_intersection_for_ghost(xghost(dir1), extents(dir1,i), xghost(face_vec), xghost(dir2), id, comp, dir1, dh, level, i) :
                   find_intersection_for_ghost(extents(dir1,i), xghost(dir1), xghost(face_vec), xghost(dir2), id, comp, dir1, dh, level, i) ;
            }

            double ghost_temp = solid.temp[comp]->item(id);

            if (dir1 == 0) {
               (i == 1) ? (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0)+del_wall(0,i),xghost(1),xghost(2),id)) :
                          (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0)-del_wall(0,i),xghost(1),xghost(2),id)) ;
            } else if (dir1 == 1) {
               (i == 1) ? (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0),xghost(1)+del_wall(0,i),xghost(2),id)) :
                          (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0),xghost(1)-del_wall(0,i),xghost(2),id)) ;
            } else if (dir1 == 2) {
               (i == 1) ? (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0),xghost(1),xghost(2)+del_wall(0,i),id)) :
                          (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0),xghost(1),xghost(2)-del_wall(0,i),id)) ;
            }
            fwall(0,i) = ghost_temp;
         }
      // if top vertex is in fluid domain
      } else if ((node.void_frac[comp]->item(p(i,1)) == 0) && (bf_intersect[dir2].offset[comp]->item(p(i,1),0) == 1)) {
         double yint = bf_intersect[dir2].value[comp]->item(p(i,1),0);
         // Condition where intersection distance is more than ghost point distance, it means that the ghost 
         // point can be projected on the wall
         if (yint >= (extents(dir2,1)-xghost(dir2))) {
            fwall(0,i) = ((xghost(dir2)+yint-extents(dir2,1))*f(i,1)+(extents(dir2,1)-xghost(dir2))*bf_intersect[dir2].field[comp]->item(p(i,1),0))/yint;
            del_wall(0,i) = MAC::abs(extents(dir1,i) - xghost(dir1));
         // Ghost point cannot be projected on the wall, as the solid surface come first
         } else {
            size_t id = (size_t) node.parID[comp]->item(p(i,0));
            if (face_vec > dir2) {
               del_wall(0,i) = (i==1) ? 
                   find_intersection_for_ghost(xghost(dir1), extents(dir1,i), xghost(dir2), xghost(face_vec), id, comp, dir1, dh, level, i) :
                   find_intersection_for_ghost(extents(dir1,i), xghost(dir1), xghost(dir2), xghost(face_vec), id, comp, dir1, dh, level, i) ;
            } else {
               del_wall(0,i) = (i==1) ? 
                   find_intersection_for_ghost(xghost(dir1), extents(dir1,i), xghost(face_vec), xghost(dir2), id, comp, dir1, dh, level, i) :
                   find_intersection_for_ghost(extents(dir1,i), xghost(dir1), xghost(face_vec), xghost(dir2), id, comp, dir1, dh, level, i) ;
            }

            double ghost_temp = solid.temp[comp]->item(id);

            if (dir1 == 0) {
               (i == 1) ? (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0)+del_wall(0,i),xghost(1),xghost(2),id)) :
                          (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0)-del_wall(0,i),xghost(1),xghost(2),id)) ;
            } else if (dir1 == 1) {
               (i == 1) ? (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0),xghost(1)+del_wall(0,i),xghost(2),id)) :
                          (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0),xghost(1)-del_wall(0,i),xghost(2),id)) ;
            } else if (dir1 == 2) {
               (i == 1) ? (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0),xghost(1),xghost(2)+del_wall(0,i),id)) :
                          (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0),xghost(1),xghost(2)-del_wall(0,i),id)) ;
            }
            fwall(0,i) = ghost_temp;
         }
      // if both vertex's are in solid domain
      } else if ((node.void_frac[comp]->item(p(i,0)) == 1) && (node.void_frac[comp]->item(p(i,1)) == 1)) {
         size_t id = (size_t) node.parID[comp]->item(p(i,0));

         if (face_vec > dir2) {
            del_wall(0,i) = (i==1) ? 
                find_intersection_for_ghost(xghost(dir1), extents(dir1,i), xghost(dir2), xghost(face_vec), id, comp, dir1, dh, level, i) :
                find_intersection_for_ghost(extents(dir1,i), xghost(dir1), xghost(dir2), xghost(face_vec), id, comp, dir1, dh, level, i) ;
         } else {
            del_wall(0,i) = (i==1) ? 
                find_intersection_for_ghost(xghost(dir1), extents(dir1,i), xghost(face_vec), xghost(dir2), id, comp, dir1, dh, level, i) :
                find_intersection_for_ghost(extents(dir1,i), xghost(dir1), xghost(face_vec), xghost(dir2), id, comp, dir1, dh, level, i) ;
         }

         double ghost_temp = solid.temp[comp]->item(id);

         if (dir1 == 0) {
            (i == 1) ? (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0)+del_wall(0,i),xghost(1),xghost(2),id)) :
                       (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0)-del_wall(0,i),xghost(1),xghost(2),id)) ;
         } else if (dir1 == 1) {
            (i == 1) ? (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0),xghost(1)+del_wall(0,i),xghost(2),id)) :
                       (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0),xghost(1)-del_wall(0,i),xghost(2),id)) ;
         } else if (dir1 == 2) {
            (i == 1) ? (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0),xghost(1),xghost(2)+del_wall(0,i),id)) :
                       (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0),xghost(1),xghost(2)-del_wall(0,i),id)) ;
         }
         fwall(0,i) = ghost_temp;
      }
   }

   // Contribution from top and bottom wall
   for (size_t j = 0; j < 2; j++) {         // 0 --> bottom; 1 --> top
      if ((node.void_frac[comp]->item(p(0,j)) == 0) && (node.void_frac[comp]->item(p(1,j)) == 0)) {
         fwall(1,j) = ((extents(dir1,1) - xghost(dir1))*f(0,j) + (xghost(dir1) - extents(dir1,0))*f(1,j))/(extents(dir1,1)-extents(dir1,0));
         del_wall(1,j) = MAC::abs(extents(dir2,j) - xghost(dir2));
      // if left vertex is in fluid domain
      } else if ((node.void_frac[comp]->item(p(0,j)) == 0) && (bf_intersect[dir1].offset[comp]->item(p(0,j),1) == 1)) {
         double xint = bf_intersect[dir1].value[comp]->item(p(0,j),1);
         if (xint >= (xghost(dir1)-extents(dir1,0))) {
            fwall(1,j) = ((extents(dir1,0)+xint-xghost(dir1))*f(0,j)+(xghost(dir1)-extents(dir1,0))*bf_intersect[dir1].field[comp]->item(p(0,j),1))/xint;
            del_wall(1,j) = MAC::abs(extents(dir2,j) - xghost(dir2));
         } else {
            size_t id = (size_t) node.parID[comp]->item(p(1,j));
            if (face_vec > dir1) {
               del_wall(1,j) = (j==1) ? 
                   find_intersection_for_ghost(xghost(dir2), extents(dir2,j), xghost(dir1), xghost(face_vec), id, comp, dir2, dh, level, j) :
                   find_intersection_for_ghost(extents(dir2,j), xghost(dir2), xghost(dir1), xghost(face_vec), id, comp, dir2, dh, level, j) ;
            } else {
               del_wall(1,j) = (j==1) ? 
                   find_intersection_for_ghost(xghost(dir2), extents(dir2,j), xghost(face_vec), xghost(dir1), id, comp, dir2, dh, level, j) :
                   find_intersection_for_ghost(extents(dir2,j), xghost(dir2), xghost(face_vec), xghost(dir1), id, comp, dir2, dh, level, j) ;
            }

            double ghost_temp = solid.temp[comp]->item(id);

            if (dir2 == 0) {
               (j == 1) ? (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0)+del_wall(1,j),xghost(1),xghost(2),id)) :
                          (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0)-del_wall(1,j),xghost(1),xghost(2),id)) ;
            } else if (dir2 == 1) {
               (j == 1) ? (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0),xghost(1)+del_wall(1,j),xghost(2),id)) :
                          (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0),xghost(1)-del_wall(1,j),xghost(2),id)) ;
            } else if (dir2 == 2) {
               (j == 1) ? (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0),xghost(1),xghost(2)+del_wall(1,j),id)) :
                          (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0),xghost(1),xghost(2)-del_wall(1,j),id)) ;
            }
            fwall(1,j) = ghost_temp;
  
         }
      // if right vertex is in fluid domain
      } else if ((node.void_frac[comp]->item(p(1,j)) == 0) && (bf_intersect[dir1].offset[comp]->item(p(1,j),0) == 1)) {
         double xint = bf_intersect[dir1].value[comp]->item(p(1,j),0);
         if (xint >= (extents(dir1,1)-xghost(dir1))) {
            fwall(1,j) = ((xghost(dir1)+xint-extents(dir1,1))*f(1,j)+(extents(dir1,1)-xghost(dir1))*bf_intersect[dir1].field[comp]->item(p(1,j),0))/xint;
            del_wall(1,j) = MAC::abs(extents(dir2,j) - xghost(dir2));
         } else {
            size_t id = (size_t) node.parID[comp]->item(p(0,j));
            if (face_vec > dir1) {
               del_wall(1,j) = (j==1) ? 
                   find_intersection_for_ghost(xghost(dir2), extents(dir2,j), xghost(dir1), xghost(face_vec), id, comp, dir2, dh, level, j) :
                   find_intersection_for_ghost(extents(dir2,j), xghost(dir2), xghost(dir1), xghost(face_vec), id, comp, dir2, dh, level, j) ;
            } else {
               del_wall(1,j) = (j==1) ? 
                   find_intersection_for_ghost(xghost(dir2), extents(dir2,j), xghost(face_vec), xghost(dir1), id, comp, dir2, dh, level, j) :
                   find_intersection_for_ghost(extents(dir2,j), xghost(dir2), xghost(face_vec), xghost(dir1), id, comp, dir2, dh, level, j) ;
            }

            double ghost_temp = solid.temp[comp]->item(id);

            if (dir2 == 0) {
               (j == 1) ? (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0)+del_wall(1,j),xghost(1),xghost(2),id)) :
                          (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0)-del_wall(1,j),xghost(1),xghost(2),id)) ;
            } else if (dir2 == 1) {
               (j == 1) ? (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0),xghost(1)+del_wall(1,j),xghost(2),id)) :
                          (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0),xghost(1)-del_wall(1,j),xghost(2),id)) ;
            } else if (dir2 == 2) {
               (j == 1) ? (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0),xghost(1),xghost(2)+del_wall(1,j),id)) :
                          (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0),xghost(1),xghost(2)-del_wall(1,j),id)) ;
            }
            fwall(1,j) = ghost_temp;
         }
      // if both vertex's are in solid domain
      } else if ((node.void_frac[comp]->item(p(0,j)) == 1) && (node.void_frac[comp]->item(p(1,j)) == 1)) {
         size_t id = (size_t) node.parID[comp]->item(p(0,j));
         if (face_vec > dir1) {
            del_wall(1,j) = (j==1) ? 
                find_intersection_for_ghost(xghost(dir2), extents(dir2,j), xghost(dir1), xghost(face_vec), id, comp, dir2, dh, level, j) :
                find_intersection_for_ghost(extents(dir2,j), xghost(dir2), xghost(dir1), xghost(face_vec), id, comp, dir2, dh, level, j) ;
         } else {
            del_wall(1,j) = (j==1) ? 
                find_intersection_for_ghost(xghost(dir2), extents(dir2,j), xghost(face_vec), xghost(dir1), id, comp, dir2, dh, level, j) :
                find_intersection_for_ghost(extents(dir2,j), xghost(dir2), xghost(face_vec), xghost(dir1), id, comp, dir2, dh, level, j) ;
         }

         double ghost_temp = solid.temp[comp]->item(id);

         if (dir2 == 0) {
            (j == 1) ? (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0)+del_wall(1,j),xghost(1),xghost(2),id)) :
                       (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0)-del_wall(1,j),xghost(1),xghost(2),id)) ;
         } else if (dir2 == 1) {
            (j == 1) ? (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0),xghost(1)+del_wall(1,j),xghost(2),id)) :
                       (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0),xghost(1)-del_wall(1,j),xghost(2),id)) ;
         } else if (dir2 == 2) {
            (j == 1) ? (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0),xghost(1),xghost(2)+del_wall(1,j),id)) :
                       (ghost_temp = impose_solid_temperature_for_ghost(comp,xghost(0),xghost(1),xghost(2)-del_wall(1,j),id)) ;
         }
         fwall(1,j) = ghost_temp;
      }
   }

   double field_value = (1./2.)*((del_wall(0,1)*fwall(0,0) + del_wall(0,0)*fwall(0,1))/(del_wall(0,1)+del_wall(0,0)) + 
                                 (del_wall(1,0)*fwall(1,1) + del_wall(1,1)*fwall(1,0))/(del_wall(1,0)+del_wall(1,1)));

   return (field_value);

}

//---------------------------------------------------------------------------
double
DDS_HeatTransfer:: ghost_field_estimate_in_box ( FV_DiscreteField* FF, size_t const& comp, size_t const& i0, size_t const& j0, size_t const& k0, double const& x0, double const& y0, double const& z0, double const& dh, size_t const& level, size_t const& parID)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_NSWithHeatTransfer:: ghost_field_estimate_in_box" ) ;

// Calculates the field value at the ghost points in the box
// near the particle boundary considering boundary affects;
// x0,y0,z0 are the ghost point coordinated; i0,j0,k0 is the
// bottom left index of the grid coordinate  

   doubleArray2D vel(dim,2,0);
   doubleArray2D del(dim,2,0);

   // Behind face
   double temp = TF->get_DOF_coordinate(k0,comp , 2); 
   double face_solid = level_set_function (parID,comp,x0,y0,z0,level_set_type)*
                       level_set_function (parID,comp,x0,y0,temp,level_set_type);

   if (face_solid > 0) {
      vel(2,0) = ghost_field_estimate_on_face (TF,comp,i0,j0,k0,x0,y0,temp,dh,2,0);
      del(2,0) = MAC::abs(temp - z0);
   } else {
      del(2,0) = find_intersection_for_ghost(temp, z0, x0, y0, parID, comp, 2, dh, 0, 0);    
      vel(2,0) = impose_solid_temperature_for_ghost(comp,x0,y0,z0-del(2,0),parID);
   }

   // Front face
   temp = TF->get_DOF_coordinate(k0+1,comp , 2); 
   face_solid = level_set_function (parID,comp,x0,y0,z0,level_set_type)*
                level_set_function (parID,comp,x0,y0,temp,level_set_type);

   if (face_solid > 0) {
      vel(2,1) = ghost_field_estimate_on_face (TF,comp,i0,j0,k0+1,x0,y0,temp,dh,2,0);
      del(2,1) = MAC::abs(temp - z0);
   } else {
      del(2,1) = find_intersection_for_ghost(z0, temp, x0, y0, parID, comp, 2, dh, 0, 1);    
      vel(2,1) = impose_solid_temperature_for_ghost(comp,x0,y0,z0+del(2,1),parID);
   }

   // Left face
   temp = TF->get_DOF_coordinate(i0,comp, 0); 
   face_solid = level_set_function (parID,comp,x0,y0,z0,level_set_type)*
                level_set_function (parID,comp,temp,y0,z0,level_set_type);

   if (face_solid > 0) {
      vel(0,0) = ghost_field_estimate_on_face (TF,comp,i0,j0,k0,temp,y0,z0,dh,0,0);
      del(0,0) = MAC::abs(temp - x0);
   } else {
      del(0,0) = find_intersection_for_ghost(temp, x0, y0, z0, parID, comp, 0, dh, 0, 0);    
      vel(0,0) = impose_solid_temperature_for_ghost(comp,x0-del(0,0),y0,z0,parID);
   }

   // Right face
   temp = TF->get_DOF_coordinate(i0+1,comp, 0); 
   face_solid = level_set_function (parID,comp,x0,y0,z0,level_set_type)*
                level_set_function (parID,comp,temp,y0,z0,level_set_type);

   if (face_solid > 0) {
      vel(0,1) = ghost_field_estimate_on_face (TF,comp,i0+1,j0,k0,temp,y0,z0,dh,0,0);
      del(0,1) = MAC::abs(temp - x0);
   } else {
      del(0,1) = find_intersection_for_ghost(x0, temp, y0, z0, parID, comp, 0, dh, 0, 1);    
      vel(0,1) = impose_solid_temperature_for_ghost(comp,x0+del(0,1),y0,z0,parID);
   }

   // Bottom face
   temp = TF->get_DOF_coordinate(j0,comp, 1); 
   face_solid = level_set_function (parID,comp,x0,y0,z0,level_set_type)*
                level_set_function (parID,comp,x0,temp,z0,level_set_type);

   if (face_solid > 0) {
      vel(1,0) = ghost_field_estimate_on_face (TF,comp,i0,j0,k0,x0,temp,z0,dh,1,0);
      del(1,0) = MAC::abs(temp - y0);
   } else {
      del(1,0) = find_intersection_for_ghost(temp, y0, x0, z0, parID, comp, 1, dh, 0, 0);    
      vel(1,0) = impose_solid_temperature_for_ghost(comp,x0,y0-del(1,0),z0,parID);
   }

   // Top face
   temp = TF->get_DOF_coordinate(j0+1,comp, 1); 
   face_solid = level_set_function (parID,comp,x0,y0,z0,level_set_type)*
                level_set_function (parID,comp,x0,temp,z0,level_set_type);

   if (face_solid > 0) {
      vel(1,1) = ghost_field_estimate_on_face (TF,comp,i0,j0+1,k0,x0,temp,z0,dh,1,0);
      del(1,1) = MAC::abs(temp - y0);
   } else {
      del(1,1) = find_intersection_for_ghost(y0, temp, x0, z0, parID, comp, 1, dh, 0, 1);    
      vel(1,1) = impose_solid_temperature_for_ghost(comp,x0,y0+del(1,1),z0,parID);
   }

   double value = (1./3.)*((vel(0,1)*del(0,0)+vel(0,0)*del(0,1))/(del(0,0)+del(0,1)) + 
                           (vel(1,1)*del(1,0)+vel(1,0)*del(1,1))/(del(1,0)+del(1,1)) +
                           (vel(2,1)*del(2,0)+vel(2,0)*del(2,1))/(del(2,0)+del(2,1)));

   return(value);
}

//---------------------------------------------------------------------------
double
DDS_HeatTransfer:: impose_solid_temperature (size_t const& comp, size_t const& dir, size_t const& off, size_t const& i, size_t const& j, size_t const& k, double const& xb, size_t const& parID )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_NSWithHeatTransfer:: impose_solid_velocity" ) ;

  doubleVector delta(3,0.);
  doubleVector grid_coord(3,0.);
  doubleVector par_coord(3,0.);

  PartInput solid = GLOBAL_EQ->get_solid();

  grid_coord(0) = TF->get_DOF_coordinate( i, comp, 0 ) ;
  grid_coord(1) = TF->get_DOF_coordinate( j, comp, 1 ) ;

  par_coord(0) = solid.coord[comp]->item(parID,0);
  par_coord(1) = solid.coord[comp]->item(parID,1);

  if (dim == 3) {
     grid_coord(2) = TF->get_DOF_coordinate( k, comp, 2 ) ;
     par_coord(2) = solid.coord[comp]->item(parID,2);
  }

  double sign = 0.;

  if (off == 0) {
     sign = -1.;
  } else if (off == 1) {
     sign = +1.;
  }

  grid_coord(dir) = grid_coord(dir) + sign*xb;

  for (size_t m = 0; m < 3; m++) {
     delta(m) = grid_coord(m) - par_coord(m);
  }

//  double ghost_temp = solid.temp[comp]->item(parID);
//  double ghost_temp = pow((grid_coord(0)-0.5)*(grid_coord(1)-0.5)*(grid_coord(2)-0.5),2.) + pow((grid_coord(0)-0.5),4.)*pow((grid_coord(1)-0.5),4.)*pow((grid_coord(2)-0.5),4.);
  double ghost_temp = 3.*MAC::sin(MAC::pi()*grid_coord(0))*MAC::sin(MAC::pi()*grid_coord(1))*MAC::sin(MAC::pi()*grid_coord(2));

  return(ghost_temp);
}

//---------------------------------------------------------------------------
double
DDS_HeatTransfer:: impose_solid_temperature_for_ghost (size_t const& comp, double const& xg, double const& yg, double const& zg, size_t const& parID )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_NSWithHeatTransfer:: impose_solid_temperature_for_ghost" ) ;

  doubleVector delta(3,0.);
  doubleVector grid_coord(3,0.);
  doubleVector par_coord(3,0.);

  PartInput solid = GLOBAL_EQ->get_solid();

  grid_coord(0) = xg;
  grid_coord(1) = yg;

  par_coord(0) = solid.coord[comp]->item(parID,0);
  par_coord(1) = solid.coord[comp]->item(parID,1);

  if (dim == 3) {
     grid_coord(2) = zg;
     par_coord(2) = solid.coord[comp]->item(parID,2);
  }

  for (size_t m = 0; m < 3; m++) {
     delta(m) = grid_coord(m) - par_coord(m);
  }

//  double ghost_temp = solid.temp[comp]->item(parID);
//  double ghost_temp = pow((grid_coord(0)-0.5)*(grid_coord(1)-0.5)*(grid_coord(2)-0.5),2.) + pow((grid_coord(0)-0.5),4.)*pow((grid_coord(1)-0.5),4.)*pow((grid_coord(2)-0.5),4.);
  double ghost_temp = 3.*MAC::sin(MAC::pi()*grid_coord(0))*MAC::sin(MAC::pi()*grid_coord(1))*MAC::sin(MAC::pi()*grid_coord(2));

  return(ghost_temp);
}

//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: compute_surface_points_on_cube(size_t const& Np)
//---------------------------------------------------------------------------
{
  MAC_LABEL("DDS_HeatTransfer:: compute_surface_points_on_cube" ) ;
/*
  ofstream outputFile ;
  std::ostringstream os2;
  os2 << "./DS_results/point_data_" << my_rank << ".csv";
  std::string filename = os2.str();
  outputFile.open(filename.c_str());
  outputFile << "x,y,z,area,nx,ny,nz" << endl;
*/
  // Structure of particle input data
  SurfaceDiscretize surface = GLOBAL_EQ->get_surface();

  double dp = 2./((double)Np);

  if (dim == 3) {
     // Generating discretization on surface
     doubleVector lsp(Np,0.);
     for (size_t i=0; i<Np; i++) lsp(i) = -1. + dp*(double(i)+0.5);

     size_t cntr = 0;

     for (size_t i=0; i<Np; i++) {
	for (size_t j=0; j<Np; j++) {
           //Front
           surface.coordinate->set_item(cntr,0,lsp(i));
           surface.coordinate->set_item(cntr,1,lsp(j));
           surface.coordinate->set_item(cntr,2,1.);
           surface.area->set_item(cntr,dp*dp);
           surface.normal->set_item(cntr,2,1.);
           //Behind
           surface.coordinate->set_item(Np*Np+cntr,0,lsp(i));
           surface.coordinate->set_item(Np*Np+cntr,1,lsp(j));
           surface.coordinate->set_item(Np*Np+cntr,2,-1.);
           surface.area->set_item(Np*Np+cntr,dp*dp);
           surface.normal->set_item(Np*Np+cntr,2,-1.);
           //Top
           surface.coordinate->set_item(2*Np*Np+cntr,0,lsp(i));
           surface.coordinate->set_item(2*Np*Np+cntr,2,lsp(j));
           surface.coordinate->set_item(2*Np*Np+cntr,1,1.);
           surface.area->set_item(2*Np*Np+cntr,dp*dp);
           surface.normal->set_item(2*Np*Np+cntr,1,1.);
           //Bottom
           surface.coordinate->set_item(3*Np*Np+cntr,0,lsp(i));
           surface.coordinate->set_item(3*Np*Np+cntr,2,lsp(j));
           surface.coordinate->set_item(3*Np*Np+cntr,1,-1.);
           surface.area->set_item(3*Np*Np+cntr,dp*dp);
           surface.normal->set_item(3*Np*Np+cntr,1,-1.);
           //Right
           surface.coordinate->set_item(4*Np*Np+cntr,1,lsp(i));
           surface.coordinate->set_item(4*Np*Np+cntr,2,lsp(j));
           surface.coordinate->set_item(4*Np*Np+cntr,0,1.);
           surface.area->set_item(4*Np*Np+cntr,dp*dp);
           surface.normal->set_item(4*Np*Np+cntr,0,1.);
           //Left
           surface.coordinate->set_item(5*Np*Np+cntr,1,lsp(i));
           surface.coordinate->set_item(5*Np*Np+cntr,2,lsp(j));
           surface.coordinate->set_item(5*Np*Np+cntr,0,-1.);
           surface.area->set_item(5*Np*Np+cntr,dp*dp);
           surface.normal->set_item(5*Np*Np+cntr,0,-1.);

/*	   outputFile << surface.coordinate->item(cntr,0) << "," << surface.coordinate->item(cntr,1) << "," << surface.coordinate->item(cntr,2) << "," << surface.area->item(cntr) << "," << surface.normal->item(cntr,0) << "," << surface.normal->item(cntr,1) << "," << surface.normal->item(cntr,2) << endl; 
           outputFile << surface.coordinate->item(Np*Np+cntr,0) << "," << surface.coordinate->item(Np*Np+cntr,1) << "," << surface.coordinate->item(Np*Np+cntr,2) << "," << surface.area->item(Np*Np+cntr) << "," << surface.normal->item(Np*Np+cntr,0) << "," << surface.normal->item(Np*Np+cntr,1) << "," << surface.normal->item(Np*Np+cntr,2) << endl; 
           outputFile << surface.coordinate->item(2*Np*Np+cntr,0) << "," << surface.coordinate->item(2*Np*Np+cntr,1) << "," << surface.coordinate->item(2*Np*Np+cntr,2) << "," << surface.area->item(2*Np*Np+cntr) << "," << surface.normal->item(2*Np*Np+cntr,0) << "," << surface.normal->item(2*Np*Np+cntr,1) << "," << surface.normal->item(2*Np*Np+cntr,2) << endl; 
           outputFile << surface.coordinate->item(3*Np*Np+cntr,0) << "," << surface.coordinate->item(3*Np*Np+cntr,1) << "," << surface.coordinate->item(3*Np*Np+cntr,2) << "," << surface.area->item(3*Np*Np+cntr) << "," << surface.normal->item(3*Np*Np+cntr,0) << "," << surface.normal->item(3*Np*Np+cntr,1) << "," << surface.normal->item(3*Np*Np+cntr,2) << endl; 
           outputFile << surface.coordinate->item(4*Np*Np+cntr,0) << "," << surface.coordinate->item(4*Np*Np+cntr,1) << "," << surface.coordinate->item(4*Np*Np+cntr,2) << "," << surface.area->item(4*Np*Np+cntr) << "," << surface.normal->item(4*Np*Np+cntr,0) << "," << surface.normal->item(4*Np*Np+cntr,1) << "," << surface.normal->item(4*Np*Np+cntr,2) << endl; 
           outputFile << surface.coordinate->item(5*Np*Np+cntr,0) << "," << surface.coordinate->item(5*Np*Np+cntr,1) << "," << surface.coordinate->item(5*Np*Np+cntr,2) << "," << surface.area->item(5*Np*Np+cntr) << "," << surface.normal->item(5*Np*Np+cntr,0) << "," << surface.normal->item(5*Np*Np+cntr,1) << "," << surface.normal->item(5*Np*Np+cntr,2) << endl; */
           cntr++;
	}
     }
  } else if (dim == 2) {
     // Generating discretization on surface
     double lsp=0.;
     for (size_t i=0; i<Np; i++) {
        lsp = -1. + dp*((double)i+0.5);
	//Bottom
	surface.coordinate->set_item(i,0,lsp);
	surface.coordinate->set_item(i,1,-1.);
	surface.area->set_item(i,dp);
	surface.normal->set_item(i,1,-1.);
	//Top
	surface.coordinate->set_item(Np+i,0,lsp);
	surface.coordinate->set_item(Np+i,1,1.);
	surface.area->set_item(Np+i,dp);
	surface.normal->set_item(Np+i,1,1.);
	//Left
	surface.coordinate->set_item(2*Np+i,0,-1.);
	surface.coordinate->set_item(2*Np+i,1,lsp);
	surface.area->set_item(2*Np+i,dp);
	surface.normal->set_item(2*Np+i,0,-1.);
	//Right
	surface.coordinate->set_item(3*Np+i,0,1.);
	surface.coordinate->set_item(3*Np+i,1,lsp);
	surface.area->set_item(3*Np+i,dp);
	surface.normal->set_item(3*Np+i,0,1.);
/*  
  	outputFile << surface.coordinate->item(i,0) << "," << surface.coordinate->item(i,1) << "," << surface.area->item(i) << "," << surface.normal->item(i,0) << "," << surface.normal->item(i,1) << endl; 
  	outputFile << surface.coordinate->item(1*Np+i,0) << "," << surface.coordinate->item(1*Np+i,1) << "," << surface.area->item(1*Np+i) << "," << surface.normal->item(1*Np+i,0) << "," << surface.normal->item(1*Np+i,1) << endl; 
  	outputFile << surface.coordinate->item(2*Np+i,0) << "," << surface.coordinate->item(2*Np+i,1) << "," << surface.area->item(2*Np+i) << "," << surface.normal->item(2*Np+i,0) << "," << surface.normal->item(2*Np+i,1) << endl; 
  	outputFile << surface.coordinate->item(3*Np+i,0) << "," << surface.coordinate->item(3*Np+i,1) << "," << surface.area->item(3*Np+i) << "," << surface.normal->item(3*Np+i,0) << "," << surface.normal->item(3*Np+i,1) << endl; */
     }
  }
//  outputFile.close();
}

//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: compute_surface_points_on_cylinder(class doubleVector& k, size_t const& Nring)
//---------------------------------------------------------------------------
{
  MAC_LABEL("DDS_HeatTransfer:: compute_surface_points_on_cylinder" ) ;
/*
  ofstream outputFile ;
  std::ostringstream os2;
  os2 << "./DS_results/point_data_" << my_rank << ".csv";
  std::string filename = os2.str();
  outputFile.open(filename.c_str());
  outputFile << "x,y,z,area,nx,ny,nz" << endl;
*/

  // Structure of particle input data
  SurfaceDiscretize surface = GLOBAL_EQ->get_surface();

  // Radius of the rings in lamber projection plane
  doubleVector Rring(Nring,0.);

  Rring(Nring-1) = 1.;

  size_t maxby2 = (size_t) k(Nring-1);

  if (dim == 3) {
     // Calculation for all rings except at the pole
     for (int i=(int)Nring-1; i>0; --i) {
        double Ri = Rring(i);
        Rring(i-1) = MAC::sqrt(k(i-1)/k(i))*Rring(i);
        Rring(i) = (Rring(i) + Rring(i-1))/2.;
        double d_theta = 2.*MAC::pi()/(k(i)-k(i-1));
        // Theta initialize as 1% of the d_theta, so there would be no chance of point overlap with mesh gridlines
        double theta = 0.01*d_theta;
        for (int j=(int)k(i-1); j<k(i); j++) {
	   // For top disk
           theta = theta + d_theta;
           surface.coordinate->set_item(j,0,Rring(i)*MAC::cos(theta));
           surface.coordinate->set_item(j,1,Rring(i)*MAC::sin(theta));
           surface.coordinate->set_item(j,2,1.);
           surface.area->set_item(j,0.5*d_theta*(pow(Ri,2)-pow(Rring(i-1),2)));
	   // For bottom disk
	   surface.coordinate->set_item(maxby2+j,0,Rring(i)*MAC::cos(theta));
           surface.coordinate->set_item(maxby2+j,1,Rring(i)*MAC::sin(theta));
           surface.coordinate->set_item(maxby2+j,2,-1.);
           surface.area->set_item(maxby2+j,0.5*d_theta*(pow(Ri,2)-pow(Rring(i-1),2)));
	   // Create surface normal vectors
	   surface.normal->set_item(j,0,0.);
	   surface.normal->set_item(j,1,0.);
	   surface.normal->set_item(j,2,1.);
	   surface.normal->set_item(maxby2+j,0,0.);
	   surface.normal->set_item(maxby2+j,1,0.);
	   surface.normal->set_item(maxby2+j,2,-1.);

/*           outputFile << surface.coordinate->item(j,0) << "," << surface.coordinate->item(j,1) << "," << surface.coordinate->item(j,2) << "," 
		      << surface.area->item(j) << "," 
		      << surface.normal->item(j,0) << "," << surface.normal->item(j,1) << "," << surface.normal->item(j,2) << endl;
           outputFile << surface.coordinate->item(maxby2+j,0) << "," << surface.coordinate->item(maxby2+j,1) << "," << surface.coordinate->item(maxby2+j,2) << "," << surface.area->item(maxby2+j) << "," << surface.normal->item(maxby2+j,0) << "," << surface.normal->item(maxby2+j,1) << "," << surface.normal->item(maxby2+j,2) << endl;*/
        }
     } 

     // Calculation at the ring on pole (i=0)
     double Ri = Rring(0);
     Rring(0) = Rring(0)/2.;
     double d_theta = 2.*MAC::pi()/(k(0));
     // Theta initialize as 1% of the d_theta, so there would be no chance of point overlap with mesh gridlines
     double theta = 0.01*d_theta;
     if (k(0)>1) {
        for (int j=0; j < k(0); j++) {
	   // For top disk
           theta = theta + d_theta;
           surface.coordinate->set_item(j,0,Rring(0)*MAC::cos(theta));
           surface.coordinate->set_item(j,1,Rring(0)*MAC::sin(theta));
           surface.coordinate->set_item(j,2,1.);
           surface.area->set_item(j,0.5*d_theta*pow(Ri,2));
           // For bottom disk
           surface.coordinate->set_item(maxby2+j,0,Rring(0)*MAC::cos(theta));
           surface.coordinate->set_item(maxby2+j,1,Rring(0)*MAC::sin(theta));
           surface.coordinate->set_item(maxby2+j,2,-1.);
           surface.area->set_item(maxby2+j,0.5*d_theta*pow(Ri,2.));
	   // Create surface normal vectors
	   surface.normal->set_item(j,0,0.);
	   surface.normal->set_item(j,1,0.);
	   surface.normal->set_item(j,2,1.);
	   surface.normal->set_item(maxby2+j,0,0.);
	   surface.normal->set_item(maxby2+j,1,0.);
	   surface.normal->set_item(maxby2+j,2,-1.);
/*           outputFile << surface.coordinate->item(j,0) << "," << surface.coordinate->item(j,1) << "," << surface.coordinate->item(j,2) << "," << surface.area->item(j) << "," << surface.normal->item(j,0) << "," << surface.normal->item(j,1) << "," << surface.normal->item(j,2) << endl;
           outputFile << surface.coordinate->item(maxby2+j,0) << "," << surface.coordinate->item(maxby2+j,1) << "," << surface.coordinate->item(maxby2+j,2) << "," << surface.area->item(maxby2+j) << "," << surface.normal->item(maxby2+j,0) << "," << surface.normal->item(maxby2+j,1) << "," << surface.normal->item(maxby2+j,2) << endl;*/
        }
     } else {
	// For top disk
        surface.coordinate->set_item(0,0,0.);
        surface.coordinate->set_item(0,1,0.);
        surface.coordinate->set_item(0,2,1.);
        surface.area->set_item(0,0.5*d_theta*pow(Ri,2.));
	// For bottom disk
        surface.coordinate->set_item(maxby2,0,0.);
        surface.coordinate->set_item(maxby2,1,0.);
        surface.coordinate->set_item(maxby2,2,-1.);
        surface.area->set_item(maxby2,0.5*d_theta*pow(Ri,2.));
        // Create surface normal vectors
        surface.normal->set_item(0,0,0.);
        surface.normal->set_item(0,1,0.);
        surface.normal->set_item(0,2,1.);
        surface.normal->set_item(maxby2,0,0.);
        surface.normal->set_item(maxby2,1,0.);
        surface.normal->set_item(maxby2,2,-1.);
/*        outputFile << surface.coordinate->item(0,0) << "," << surface.coordinate->item(0,1) << "," << surface.coordinate->item(0,2) << "," << surface.area->item(0) << "," << surface.normal->item(0,0) << "," << surface.normal->item(0,1) << "," << surface.normal->item(0,2) << endl;
        outputFile << surface.coordinate->item(maxby2,0) << "," << surface.coordinate->item(maxby2,1) << "," << surface.coordinate->item(maxby2,2) << "," << surface.area->item(maxby2) << "," << surface.normal->item(maxby2,0) << "," << surface.normal->item(maxby2,1) << "," << surface.normal->item(maxby2,2) << endl;*/
     }

     // Generating one ring of points on cylindrical surface
     // Can be used to calculate stress on whole surface by a constant shift of points

     // Estimating number of points on cylindrical surface
     int pts_1_ring = (int)(k(Nring-1) - k(Nring-2));
     int cyl_rings = (int)round(2./(1-MAC::sqrt(k(Nring-2)/k(Nring-1))));
     double cell_area = 2.*MAC::pi()/(double(pts_1_ring))*(2./(double(cyl_rings)));

     d_theta = 2.*MAC::pi()/pts_1_ring;
     for (int j=0; j<cyl_rings; j++) {
        theta = 0.01*d_theta;
	for (int ij=0; ij<pts_1_ring; ij++) {
           theta = theta + d_theta;
	   int n = 2*(int)k(Nring-1) + j*pts_1_ring + ij;
           surface.coordinate->set_item(n,0,MAC::cos(theta));
           surface.coordinate->set_item(n,1,MAC::sin(theta));
           surface.coordinate->set_item(n,2,-1.+ 2.*(j+0.5)/(double(cyl_rings)));
           surface.area->set_item(n,cell_area);
           surface.normal->set_item(n,0,MAC::cos(theta));
           surface.normal->set_item(n,1,MAC::sin(theta));
           surface.normal->set_item(n,2,0.);

//           outputFile << surface.coordinate->item(2*k(Nring-1)+j*ij,0) << "," << surface.coordinate->item(2*k(Nring-1)+j*ij,1) << "," << surface.coordinate->item(2*k(Nring-1)+j*ij,2) << "," << surface.area->item(2*k(Nring-1)+j*ij) << "," << surface.normal->item(2*k(Nring-1)+j*ij,0) << "," << surface.normal->item(2*k(Nring-1)+j*ij,1) << "," << surface.normal->item(2*k(Nring-1)+j*ij,2) << endl;
	}
     }
  }
//  outputFile.close();
}

//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: compute_surface_points_on_sphere(class doubleVector& eta, class doubleVector& k, class doubleVector& Rring, size_t const& Nring)
//---------------------------------------------------------------------------
{
  MAC_LABEL("DDS_HeatTransfer:: compute_surface_points_on_sphere" ) ;
/*
  ofstream outputFile ;
  std::ostringstream os2;
  os2 << "./DS_results/point_data_" << my_rank << ".csv";
  std::string filename = os2.str();
  outputFile.open(filename.c_str());
  outputFile << "x,y,z,area,nx,ny,nz" << endl;
*/

  // Structure of particle input data
  SurfaceDiscretize surface = GLOBAL_EQ->get_surface();

  size_t maxby2 = (size_t) k(Nring-1);

  if (dim == 3) {
     // Calculation for all rings except at the pole
     for (int i=(int)Nring-1; i>0; --i) {
        double Ri = Rring(i);
        Rring(i) = (Rring(i) + Rring(i-1))/2.;
        eta(i) = (eta(i) + eta(i-1))/2.;
        double d_theta = 2.*MAC::pi()/(k(i)-k(i-1));
        // Theta initialize as 1% of the d_theta, so there would be no chance of point overlap with mesh gridlines
        double theta = 0.01*d_theta;
        for (int j=(int)k(i-1); j<k(i); j++) {
           theta = theta + d_theta;
           if (pole_loc == 2) {
              surface.coordinate->set_item(j,0,MAC::cos(theta)*MAC::sin(eta(i)));
              surface.coordinate->set_item(j,1,MAC::sin(theta)*MAC::sin(eta(i)));
              surface.coordinate->set_item(j,2,MAC::cos(eta(i)));
              surface.area->set_item(j,0.5*d_theta*(pow(Ri,2.)-pow(Rring(i-1),2.)));
              // For second half of sphere
              surface.coordinate->set_item(maxby2+j,0,MAC::cos(theta)*MAC::sin(eta(i)));
              surface.coordinate->set_item(maxby2+j,1,MAC::sin(theta)*MAC::sin(eta(i)));
              surface.coordinate->set_item(maxby2+j,2,-MAC::cos(eta(i)));
              surface.area->set_item(maxby2+j,0.5*d_theta*(pow(Ri,2.)-pow(Rring(i-1),2.)));
           } else if (pole_loc == 1) {
              surface.coordinate->set_item(j,2,MAC::cos(theta)*MAC::sin(eta(i)));
              surface.coordinate->set_item(j,0,MAC::sin(theta)*MAC::sin(eta(i)));
              surface.coordinate->set_item(j,1,MAC::cos(eta(i)));
              surface.area->set_item(j,0.5*d_theta*(pow(Ri,2.)-pow(Rring(i-1),2.)));
              // For second half of sphere
              surface.coordinate->set_item(maxby2+j,2,MAC::cos(theta)*MAC::sin(eta(i)));
              surface.coordinate->set_item(maxby2+j,0,MAC::sin(theta)*MAC::sin(eta(i)));
              surface.coordinate->set_item(maxby2+j,1,-MAC::cos(eta(i)));
              surface.area->set_item(maxby2+j,0.5*d_theta*(pow(Ri,2.)-pow(Rring(i-1),2.)));
           } else if (pole_loc == 0) {
              surface.coordinate->set_item(j,1,MAC::cos(theta)*MAC::sin(eta(i)));
              surface.coordinate->set_item(j,2,MAC::sin(theta)*MAC::sin(eta(i)));
              surface.coordinate->set_item(j,0,MAC::cos(eta(i)));
              surface.area->set_item(j,0.5*d_theta*(pow(Ri,2.)-pow(Rring(i-1),2.)));
              // For second half of sphere
              surface.coordinate->set_item(maxby2+j,1,MAC::cos(theta)*MAC::sin(eta(i)));
              surface.coordinate->set_item(maxby2+j,2,MAC::sin(theta)*MAC::sin(eta(i)));
              surface.coordinate->set_item(maxby2+j,0,-MAC::cos(eta(i)));
              surface.area->set_item(maxby2+j,0.5*d_theta*(pow(Ri,2.)-pow(Rring(i-1),2.)));
           }
	   // Create surface normal vectors
	   surface.normal->set_item(j,0,surface.coordinate->item(j,0));
	   surface.normal->set_item(j,1,surface.coordinate->item(j,1));
	   surface.normal->set_item(j,2,surface.coordinate->item(j,2));
	   surface.normal->set_item(maxby2+j,0,surface.coordinate->item(maxby2+j,0));
	   surface.normal->set_item(maxby2+j,1,surface.coordinate->item(maxby2+j,1));
	   surface.normal->set_item(maxby2+j,2,surface.coordinate->item(maxby2+j,2));

/*           outputFile << surface.coordinate->item(j,0) << "," << surface.coordinate->item(j,1) << "," << surface.coordinate->item(j,2) << "," 
		      << surface.area->item(j) << "," 
		      << surface.normal->item(j,0) << "," << surface.normal->item(j,1) << "," << surface.normal->item(j,2) << endl;
           outputFile << surface.coordinate->item(maxby2+j,0) << "," << surface.coordinate->item(maxby2+j,1) << "," << surface.coordinate->item(maxby2+j,2) << "," << surface.area->item(maxby2+j) << "," << surface.normal->item(maxby2+j,0) << "," << surface.normal->item(maxby2+j,1) << "," << surface.normal->item(maxby2+j,2) << endl;*/
        }
     } 

     // Calculation at the ring on pole (i=0)
     double Ri = Rring(0);
     Rring(0) = Rring(0)/2.;
     eta(0) = eta(0)/2.;
     double d_theta = 2.*MAC::pi()/(k(0));
     // Theta initialize as 1% of the d_theta, so there would be no chance of point overlap with mesh gridlines
     double theta = 0.01*d_theta;
     if (k(0)>1) {
        for (int j=0; j < k(0); j++) {
           theta = theta + d_theta;
           if (pole_loc == 2) {
              surface.coordinate->set_item(j,0,MAC::cos(theta)*MAC::sin(eta(0)));
              surface.coordinate->set_item(j,1,MAC::sin(theta)*MAC::sin(eta(0)));
              surface.coordinate->set_item(j,2,MAC::cos(eta(0)));
              surface.area->set_item(j,0.5*d_theta*pow(Ri,2.));
              // For second half of sphere
              surface.coordinate->set_item(maxby2+j,0,MAC::cos(theta)*MAC::sin(eta(0)));
              surface.coordinate->set_item(maxby2+j,1,MAC::sin(theta)*MAC::sin(eta(0)));
              surface.coordinate->set_item(maxby2+j,2,-MAC::cos(eta(0)));
              surface.area->set_item(maxby2+j,0.5*d_theta*pow(Ri,2.));
           } else if (pole_loc == 1) {
              surface.coordinate->set_item(j,2,MAC::cos(theta)*MAC::sin(eta(0)));
              surface.coordinate->set_item(j,0,MAC::sin(theta)*MAC::sin(eta(0)));
              surface.coordinate->set_item(j,1,MAC::cos(eta(0)));
              surface.area->set_item(j,0.5*d_theta*pow(Ri,2.));
              // For second half of sphere
              surface.coordinate->set_item(maxby2+j,2,MAC::cos(theta)*MAC::sin(eta(0)));
              surface.coordinate->set_item(maxby2+j,0,MAC::sin(theta)*MAC::sin(eta(0)));
              surface.coordinate->set_item(maxby2+j,1,-MAC::cos(eta(0)));
              surface.area->set_item(maxby2+j,0.5*d_theta*pow(Ri,2.));
           } else if (pole_loc == 0) {
              surface.coordinate->set_item(j,1,MAC::cos(theta)*MAC::sin(eta(0)));
              surface.coordinate->set_item(j,2,MAC::sin(theta)*MAC::sin(eta(0)));
              surface.coordinate->set_item(j,0,MAC::cos(eta(0)));
              surface.area->set_item(j,0.5*d_theta*pow(Ri,2.));
              // For second half of sphere
              surface.coordinate->set_item(maxby2+j,1,MAC::cos(theta)*MAC::sin(eta(0)));
              surface.coordinate->set_item(maxby2+j,2,MAC::sin(theta)*MAC::sin(eta(0)));
              surface.coordinate->set_item(maxby2+j,0,-MAC::cos(eta(0)));
              surface.area->set_item(maxby2+j,0.5*d_theta*pow(Ri,2.));
           } 
	   // Create surface normal vectors
	   surface.normal->set_item(j,0,surface.coordinate->item(j,0));
	   surface.normal->set_item(j,1,surface.coordinate->item(j,1));
	   surface.normal->set_item(j,2,surface.coordinate->item(j,2));
	   surface.normal->set_item(maxby2+j,0,surface.coordinate->item(maxby2+j,0));
	   surface.normal->set_item(maxby2+j,1,surface.coordinate->item(maxby2+j,1));
	   surface.normal->set_item(maxby2+j,2,surface.coordinate->item(maxby2+j,2));
/*           outputFile << surface.coordinate->item(j,0) << "," << surface.coordinate->item(j,1) << "," << surface.coordinate->item(j,2) << "," << surface.area->item(j) << "," << surface.normal->item(j,0) << "," << surface.normal->item(j,1) << "," << surface.normal->item(j,2) << endl;
           outputFile << surface.coordinate->item(maxby2+j,0) << "," << surface.coordinate->item(maxby2+j,1) << "," << surface.coordinate->item(maxby2+j,2) << "," << surface.area->item(maxby2+j) << "," << surface.normal->item(maxby2+j,0) << "," << surface.normal->item(maxby2+j,1) << "," << surface.normal->item(maxby2+j,2) << endl;*/
        }
     } else {
        if (pole_loc == 2) { 
           surface.coordinate->set_item(0,0,0.);
           surface.coordinate->set_item(0,1,0.);
           surface.coordinate->set_item(0,2,1.);
           surface.area->set_item(0,0.5*d_theta*pow(Ri,2.));
           // For second half of sphere
           surface.coordinate->set_item(maxby2,0,0.);
           surface.coordinate->set_item(maxby2,1,0.);
           surface.coordinate->set_item(maxby2,2,-1.);
           surface.area->set_item(maxby2,0.5*d_theta*pow(Ri,2.));
        } else if (pole_loc == 1) {
           surface.coordinate->set_item(0,2,0.);
           surface.coordinate->set_item(0,0,0.);
           surface.coordinate->set_item(0,1,1.);
           surface.area->set_item(0,0.5*d_theta*pow(Ri,2.));
           // For second half of sphere
           surface.coordinate->set_item(maxby2,2,0.);
           surface.coordinate->set_item(maxby2,0,0.);
           surface.coordinate->set_item(maxby2,1,-1.);
           surface.area->set_item(maxby2,0.5*d_theta*pow(Ri,2.));
        } else if (pole_loc == 0) {
           surface.coordinate->set_item(0,1,0.);
           surface.coordinate->set_item(0,2,0.);
           surface.coordinate->set_item(0,0,1.);
           surface.area->set_item(0,0.5*d_theta*pow(Ri,2.));
           // For second half of sphere
           surface.coordinate->set_item(maxby2,1,0.);
           surface.coordinate->set_item(maxby2,2,0.);
           surface.coordinate->set_item(maxby2,0,-1.);
           surface.area->set_item(maxby2,0.5*d_theta*pow(Ri,2.));
        } 
        // Create surface normal vectors
        surface.normal->set_item(0,0,surface.coordinate->item(0,0));
        surface.normal->set_item(0,1,surface.coordinate->item(0,1));
        surface.normal->set_item(0,2,surface.coordinate->item(0,2));
        surface.normal->set_item(maxby2,0,surface.coordinate->item(maxby2,0));
        surface.normal->set_item(maxby2,1,surface.coordinate->item(maxby2,1));
        surface.normal->set_item(maxby2,2,surface.coordinate->item(maxby2,2));
/*        outputFile << surface.coordinate->item(0,0) << "," << surface.coordinate->item(0,1) << "," << surface.coordinate->item(0,2) << "," << surface.area->item(0) << "," << surface.normal->item(0,0) << "," << surface.normal->item(0,1) << "," << surface.normal->item(0,2) << endl;
        outputFile << surface.coordinate->item(maxby2,0) << "," << surface.coordinate->item(maxby2,1) << "," << surface.coordinate->item(maxby2,2) << "," << surface.area->item(maxby2) << "," << surface.normal->item(maxby2,0) << "," << surface.normal->item(maxby2,1) << "," << surface.normal->item(maxby2,2) << endl;*/
     }
  } else if (dim == 2) {
     double d_theta = 2.*MAC::pi()/(double(Nring));
     double theta = 0.01*d_theta;
     for (int j=0; j < (int) Nring; j++) {
        theta = theta + d_theta;
        surface.coordinate->set_item(j,0,MAC::cos(theta));
        surface.coordinate->set_item(j,1,MAC::sin(theta));
        surface.area->set_item(j,d_theta);
        // Create surface normal vectors
        surface.normal->set_item(j,0,surface.coordinate->item(j,0));
        surface.normal->set_item(j,1,surface.coordinate->item(j,1));
//        outputFile << surface.coordinate->item(j,0) << "," << surface.coordinate->item(j,1) << "," << surface.coordinate->item(j,2) << "," << surface.area->item(j) << "," << surface.normal->item(j,0) << "," << surface.normal->item(j,1) << "," << surface.normal->item(j,2) << endl;
     }
  }
//  outputFile.close();
}
//---------------------------------------------------------------------------
void DDS_HeatTransfer:: generate_surface_discretization()
//---------------------------------------------------------------------------
{
  MAC_LABEL("DDS_HeatTransfer:: generate_surface_discretization" ) ;

  size_t kmax = (int) Npoints;

  if ( level_set_type == "Sphere" ) {
     // Reference paper: Becker and Becker, A general rule for disk and hemisphere partition into 
     // equal-area cells, Computational Geometry 45 (2012) 275-283

     double eta_temp = MAC::pi()/2.;
     double k_temp = (double) kmax;
     double Ro_temp = MAC::sqrt(2);
     double Rn_temp = MAC::sqrt(2);
     size_t cntr = 0; 

     // Estimating the number of rings on the hemisphere
     while (k_temp > double(Pmin+2)) {
        eta_temp = eta_temp - 2./ar*MAC::sqrt(MAC::pi()/k_temp)*MAC::sin(eta_temp/2.);
        Rn_temp = 2.*MAC::sin(eta_temp/2.);
        k_temp = round(k_temp*pow(Rn_temp/Ro_temp,2.));
        Ro_temp = Rn_temp;
        cntr++;
     }

     size_t Nrings = cntr+1;

     // Summation of total discretized points with increase in number of rings radially
     doubleVector k(Nrings,0.);
     // Zenithal angle for the sphere
     doubleVector eta(Nrings,0.);
     // Radius of the rings in lamber projection plane
     doubleVector Rring(Nrings,0.);

     // Assigning the maximum number of discretized points to the last element of the array
     k(Nrings-1) = (double) kmax;
     // Zenithal angle for the last must be pi/2.
     eta(Nrings-1) = MAC::pi()/2.;
     // Radius of last ring in lamber projection plane
     Rring(Nrings-1) = MAC::sqrt(2.);

     for (int i=int(Nrings)-2; i>=0; --i) {
        eta(i) = eta(i+1) - 2./ar*MAC::sqrt(MAC::pi()/k(i+1))*MAC::sin(eta(i+1)/2.);
        Rring(i) = 2.*MAC::sin(eta(i)/2.);
        k(i) = round(k(i+1)*pow(Rring(i)/Rring(i+1),2.));
        if (i==0) k(0) = (double) Pmin;
     } 

     // Discretize the particle surface into approximate equal area cells
     if (dim == 3) {
        compute_surface_points_on_sphere(eta, k, Rring, Nrings);
     } else {
        compute_surface_points_on_sphere(eta, k, Rring, kmax);
     }
  } else if (level_set_type == "Cube") {
     compute_surface_points_on_cube(kmax);
  } else if (level_set_type == "Cylinder") {
     // Reference paper: Becker and Becker, A general rule for disk and hemisphere partition into 
     // equal-area cells, Computational Geometry 45 (2012) 275-283

     double p = MAC::pi()/ar;
     double k_temp = (double) kmax;
     size_t cntr = 0; 

     // Estimating the number of rings on either of the disc
     while (k_temp > double(Pmin+2)) {
        k_temp = round(pow(MAC::sqrt(k_temp) - MAC::sqrt(p),2.));
        cntr++;
     }

     size_t Nrings = cntr+1;

     // Summation of total discretized points with increase in number of rings radially
     doubleVector k(Nrings,0.);
     // Assigning the maximum number of discretized points to the last element of the array
     k(Nrings-1) = (double) kmax;

     for (int i=(int)Nrings-2; i>=0; --i) {
	k(i) = round(pow(MAC::sqrt(k(i+1)) - MAC::sqrt(p),2.));
        if (i==0) k(0) = (double) Pmin;
     } 

     if (dim == 3) {
        compute_surface_points_on_cylinder(k, Nrings);
     } 
  }
}

//---------------------------------------------------------------------------
double
DDS_HeatTransfer:: level_set_function (size_t const& m, size_t const& comp, double const& xC, double const& yC, double const& zC, string const& type)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_HeatTransfer:: level_set_solids" ) ;

  PartInput solid = GLOBAL_EQ->get_solid();

  double xp = solid.coord[comp]->item(m,0);
  double yp = solid.coord[comp]->item(m,1);
  double zp = solid.coord[comp]->item(m,2);
  double Rp = solid.size[comp]->item(m);

  doubleVector delta(3,0);

  delta(0) = xC-xp;
  delta(1) = yC-yp;
  delta(2) = 0;
  if (dim == 3) delta(2) = zC-zp;

  // Displacement correction in case of periodic boundary condition in any or all directions
  for (size_t dir=0;dir<dim;dir++) {
     if (is_iperiodic[dir]) {
        double isize = TF->primary_grid()->get_main_domain_max_coordinate(dir) - TF->primary_grid()->get_main_domain_min_coordinate(dir);
        delta(dir) = delta(dir) - round(delta(dir)/isize)*isize;
     }
  }

  // Try to add continuous level set function; solver performs better in this case for nodes at interface
  double level_set = 0.;
  if (type == "Sphere") {
     level_set = pow(pow(delta(0),2.)+pow(delta(1),2.)+pow(delta(2),2.),0.5)-Rp;
  } else if (type == "Ellipsoid") {
     // Solid object rotation, if any     
     trans_rotation_matrix(m,delta);
     level_set = pow(delta(0)/1.,2.)+pow(delta(1)/0.5,2.)+pow(delta(2)/0.5,2.)-Rp;
  } else if (type == "Superquadric") {
     // Solid object rotation, if any     
     trans_rotation_matrix(m,delta);
     level_set = pow(pow(delta(0),4.)+pow(delta(1),4.)+pow(delta(2),4.),0.25)-Rp;
  } else if (type == "PipeX") {
     level_set = pow(pow(delta(1),2.)+pow(delta(2),2.),0.5)-Rp;
  } else if (type == "Cube") {
     // Solid object rotation, if any     
     trans_rotation_matrix(m,delta);
     delta(0) = MAC::abs(delta(0)) - Rp;
     delta(1) = MAC::abs(delta(1)) - Rp;
     delta(2) = MAC::abs(delta(2)) - Rp;

     if ((delta(0) < 0.) && (delta(1) < 0.) && (delta(2) < 0.)) {
        level_set = MAC::min(delta(0),MAC::min(delta(1),delta(2)));
     } else {
        level_set = MAC::max(delta(0),MAC::max(delta(1),delta(2)));
     }
  } else if (type == "Cylinder") {
     // Solid object rotation, if any     
     trans_rotation_matrix(m,delta);

     level_set = pow(pow(delta(0),2.)+pow(delta(1),2.),0.5)-Rp;
     if ((MAC::abs(delta(2)) < Rp) && (level_set < 0.)) {
        level_set = MAC::abs(MAC::abs(delta(2))-Rp)*level_set;
     } else {
        level_set = MAC::max(MAC::abs(MAC::abs(delta(2))-Rp),MAC::abs(level_set));
     }
  }

  return(level_set);
}

//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: trans_rotation_matrix (size_t const& m, class doubleVector& delta)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_HeatTransfer:: trans_rotation_matrix" ) ;

  PartInput solid = GLOBAL_EQ->get_solid();

  // yaw along z-axis; pitch along y-axis; roll along x-axis
  doubleArray2D rot_matrix(3,3,0);

  // Rotation matrix assemble
  if (insertion_type == "file") {
     double roll = (MAC::pi()/180.)*solid.thetap[0]->item(m,0);
     double pitch = (MAC::pi()/180.)*solid.thetap[0]->item(m,1);
     double yaw = (MAC::pi()/180.)*solid.thetap[0]->item(m,2);
     rot_matrix(0,0) = MAC::cos(yaw)*MAC::cos(pitch);
     rot_matrix(1,0) = MAC::cos(yaw)*MAC::sin(pitch)*MAC::sin(roll) - MAC::sin(yaw)*MAC::cos(roll);
     rot_matrix(2,0) = MAC::cos(yaw)*MAC::sin(pitch)*MAC::cos(roll) + MAC::sin(yaw)*MAC::sin(roll);
     rot_matrix(0,1) = MAC::sin(yaw)*MAC::cos(pitch);
     rot_matrix(1,1) = MAC::sin(yaw)*MAC::sin(pitch)*MAC::sin(roll) + MAC::cos(yaw)*MAC::cos(roll);
     rot_matrix(2,1) = MAC::sin(yaw)*MAC::sin(pitch)*MAC::cos(roll) - MAC::cos(yaw)*MAC::sin(roll);
     rot_matrix(0,2) = -MAC::sin(pitch);
     rot_matrix(1,2) = MAC::cos(pitch)*MAC::sin(roll);
     rot_matrix(2,2) = MAC::cos(pitch)*MAC::cos(roll);
  } else if (insertion_type == "GRAINS") {
     rot_matrix(0,0) = solid.thetap[0]->item(m,0);
     rot_matrix(1,0) = solid.thetap[0]->item(m,1);
     rot_matrix(2,0) = solid.thetap[0]->item(m,2);
     rot_matrix(0,1) = solid.thetap[0]->item(m,3);
     rot_matrix(1,1) = solid.thetap[0]->item(m,4);
     rot_matrix(2,1) = solid.thetap[0]->item(m,5);
     rot_matrix(0,2) = solid.thetap[0]->item(m,6);
     rot_matrix(1,2) = solid.thetap[0]->item(m,7);
     rot_matrix(2,2) = solid.thetap[0]->item(m,8);
  }

  double delta_x = delta(0)*rot_matrix(0,0) + delta(1)*rot_matrix(0,1) + delta(2)*rot_matrix(0,2);
  double delta_y = delta(0)*rot_matrix(1,0) + delta(1)*rot_matrix(1,1) + delta(2)*rot_matrix(1,2);
  double delta_z = delta(0)*rot_matrix(2,0) + delta(1)*rot_matrix(2,1) + delta(2)*rot_matrix(2,2);

  delta(0) = delta_x;
  delta(1) = delta_y;
  delta(2) = delta_z;
}

//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: rotation_matrix (size_t const& m, class doubleVector& delta)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_HeatTransfer:: rotation_matrix" ) ;

  PartInput solid = GLOBAL_EQ->get_solid();

  // yaw along z-axis; pitch along y-axis; roll along x-axis
  doubleArray2D rot_matrix(3,3,0);

  // Rotation matrix assemble
  if (insertion_type == "file") {
     double roll = (MAC::pi()/180.)*solid.thetap[0]->item(m,0);
     double pitch = (MAC::pi()/180.)*solid.thetap[0]->item(m,1);
     double yaw = (MAC::pi()/180.)*solid.thetap[0]->item(m,2);
     rot_matrix(0,0) = MAC::cos(yaw)*MAC::cos(pitch);
     rot_matrix(0,1) = MAC::cos(yaw)*MAC::sin(pitch)*MAC::sin(roll) - MAC::sin(yaw)*MAC::cos(roll);
     rot_matrix(0,2) = MAC::cos(yaw)*MAC::sin(pitch)*MAC::cos(roll) + MAC::sin(yaw)*MAC::sin(roll);
     rot_matrix(1,0) = MAC::sin(yaw)*MAC::cos(pitch);
     rot_matrix(1,1) = MAC::sin(yaw)*MAC::sin(pitch)*MAC::sin(roll) + MAC::cos(yaw)*MAC::cos(roll);
     rot_matrix(1,2) = MAC::sin(yaw)*MAC::sin(pitch)*MAC::cos(roll) - MAC::cos(yaw)*MAC::sin(roll);
     rot_matrix(2,0) = -MAC::sin(pitch);
     rot_matrix(2,1) = MAC::cos(pitch)*MAC::sin(roll);
     rot_matrix(2,2) = MAC::cos(pitch)*MAC::cos(roll);
  } else if (insertion_type == "GRAINS") {
     rot_matrix(0,0) = solid.thetap[0]->item(m,0);
     rot_matrix(0,1) = solid.thetap[0]->item(m,1);
     rot_matrix(0,2) = solid.thetap[0]->item(m,2);
     rot_matrix(1,0) = solid.thetap[0]->item(m,3);
     rot_matrix(1,1) = solid.thetap[0]->item(m,4);
     rot_matrix(1,2) = solid.thetap[0]->item(m,5);
     rot_matrix(2,0) = solid.thetap[0]->item(m,6);
     rot_matrix(2,1) = solid.thetap[0]->item(m,7);
     rot_matrix(2,2) = solid.thetap[0]->item(m,8);
  }
     
  double delta_x = delta(0)*rot_matrix(0,0) + delta(1)*rot_matrix(0,1) + delta(2)*rot_matrix(0,2);
  double delta_y = delta(0)*rot_matrix(1,0) + delta(1)*rot_matrix(1,1) + delta(2)*rot_matrix(1,2);
  double delta_z = delta(0)*rot_matrix(2,0) + delta(1)*rot_matrix(2,1) + delta(2)*rot_matrix(2,2);

  delta(0) = delta_x;
  delta(1) = delta_y;
  delta(2) = delta_z;
}

//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: node_property_calculation ( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_HeatTransfer:: node_property_calculation" ) ;

  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);

  PartInput solid = GLOBAL_EQ->get_solid();
  NodeProp node = GLOBAL_EQ->get_node_property();

  for (size_t comp=0;comp<nb_comps;comp++) {
     // Get local min and max indices; Calculation on the rows next to the proc as well
     for (size_t l=0;l<dim;++l) {
/*        min_unknown_index(l) = TF->get_min_index_unknown_handled_by_proc( comp, l ) - 1 ;
        max_unknown_index(l) = TF->get_max_index_unknown_handled_by_proc( comp, l ) + 1 ;*/
        min_unknown_index(l) = TF->get_min_index_unknown_on_proc( comp, l );
        max_unknown_index(l) = TF->get_max_index_unknown_on_proc( comp, l );
     }

     size_t local_min_k = 0;
     size_t local_max_k = 0;

     if (dim == 3) {
        local_min_k = min_unknown_index(2);
        local_max_k = max_unknown_index(2);
     }

     // Clearing the matrices to store new time step value
     node.void_frac[comp]->nullify();
     node.parID[comp]->nullify();

     for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i) {
        double xC = TF->get_DOF_coordinate( i, comp, 0 ) ;
        double dx = TF->get_cell_size(i,comp,0) ;
        for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j) {
           double yC = TF->get_DOF_coordinate( j, comp, 1 ) ;
           double dy = TF->get_cell_size(j,comp,1) ;
           for (size_t k=local_min_k;k<=local_max_k;++k) {
              double zC = 0.;
              double dC = min(dx,dy);
              if (dim == 3) {
                 zC = TF->get_DOF_coordinate( k, comp, 2 ) ;
                 double dz = TF->get_cell_size(k,comp,2) ;
                 dC = min(dC,dz);
              }
              size_t p = return_node_index(TF,comp,i,j,k);
              for (size_t m=0;m<Npart;m++) {
                 double level_set = level_set_function(m,comp,xC,yC,zC,level_set_type);
                 level_set *= solid.inside[comp]->item(m);  

                 // level_set is xb, if local critical time scale is 0.01 of the global time scale 
                 // then the node is considered inside the solid object
                 // (xb/dC)^2 = 0.01 --> (xb/xC) = 0.1
                 if (level_set <= pow(loc_thres,0.5)*dC) {
                 //if (level_set <= 1.E-1*dC) {
                    node.void_frac[comp]->set_item(p,1.);
                    node.parID[comp]->set_item(p,(double)m);
                    break;
                 }
              }
           }
        }
     }

     // Level 0 is for the intersection matrix corresponding to fluid side
     assemble_intersection_matrix(comp,0);
     // Level 1 is for the intersection matrix corresponding to solids side
     assemble_intersection_matrix(comp,1);

//     BoundaryBisec* b_intersect = GLOBAL_EQ->get_b_intersect(0);
//     b_intersect[0].value[comp]->print_items(MAC::out(),0);
  }
}

//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: assemble_intersection_matrix ( size_t const& comp, size_t const& level)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_HeatTransfer:: assemble_intersection_matrix" ) ;

  size_t_vector min_index(dim,0);
  size_t_vector max_index(dim,0);
//  size_t_vector min_unknown_index(dim,0);
//  size_t_vector max_unknown_index(dim,0);
  size_t_vector ipos(3,0);
  size_t_array2D local_extents(dim,2,0);
  size_t_array2D node_neigh(dim,2,0);

  NodeProp node = GLOBAL_EQ->get_node_property();
  BoundaryBisec* b_intersect = GLOBAL_EQ->get_b_intersect(level);

  for (size_t l=0;l<dim;++l) {
      // To include knowns at dirichlet boundary in the intersection calculation as well, important in cases where the particle is close to domain boundary pow(2,64)
//     min_unknown_index(l) = TF->get_min_index_unknown_on_proc( comp, l );
//     max_unknown_index(l) = TF->get_max_index_unknown_on_proc( comp, l );
     min_index(l) = 0 ; 
     max_index(l) = TF->get_local_nb_dof( comp, l ) ;
     local_extents(l,0) = 0;
     local_extents(l,1) = (max_index(l)-min_index(l));
  }

  size_t local_min_k = 0;
  size_t local_max_k = 1;

  if (dim == 3) {
     local_min_k = min_index(2);
     local_max_k = max_index(2);
  }

  // Clearing the matrices to store new time step value
  for (size_t dir=0;dir<dim;dir++) {
     b_intersect[dir].offset[comp]->nullify();
     b_intersect[dir].value[comp]->nullify();
     b_intersect[dir].field[comp]->nullify();            
  }

  for (size_t i=min_index(0);i<max_index(0);++i) {
     ipos(0) = i - min_index(0);
     for (size_t j=min_index(1);j<max_index(1);++j) {
        ipos(1) = j - min_index(1);
        for (size_t k=local_min_k;k<local_max_k;++k) {
           ipos(2) = k - local_min_k;
           size_t p = return_node_index(TF,comp,i,j,k);

           double center_void_frac = 0.;
           if (level == 0) {
              center_void_frac = 1.;
           } else if (level == 1) {
              center_void_frac = 0.;
           }

           if (node.void_frac[comp]->item(p) != center_void_frac) {
              node_neigh(0,0) = return_node_index(TF,comp,i-1,j,k);
              node_neigh(0,1) = return_node_index(TF,comp,i+1,j,k);
              node_neigh(1,0) = return_node_index(TF,comp,i,j-1,k);
              node_neigh(1,1) = return_node_index(TF,comp,i,j+1,k);
              if (dim == 3) {
                 node_neigh(2,0) = return_node_index(TF,comp,i,j,k-1);
                 node_neigh(2,1) = return_node_index(TF,comp,i,j,k+1);
              }

              for (size_t dir=0;dir<dim;dir++) {
                 size_t ii,jj,kk;
                 if (dir == 0) {
                    ii = i;jj = j; kk = k;
                 } else if (dir == 1) {
                    ii = j;jj = i; kk = k;
                 } else if (dir == 2) {
                    ii = k;jj = i; kk = j;
                 }
                 for (size_t off=0;off<2;off++) {
		    // Checking if the nodes are on domain boundary or not, 
		    // if so, the check the intersection only on one side
		    if (ipos(dir) != local_extents(dir,off)) {
                       size_t left, right;
                       if (off == 0) {
                          left = ii-1; right = ii;
                       } else if (off == 1) {
                          left = ii; right = ii+1;
                       }

                       if (node.void_frac[comp]->item(node_neigh(dir,off)) != node.void_frac[comp]->item(p)) {
                          double xb = find_intersection(left,right,jj,kk,comp,dir,off,0);
                          b_intersect[dir].offset[comp]->set_item(p,off,1);
                          b_intersect[dir].value[comp]->set_item(p,off,xb);
                          // ID of particle having the intersection with node i
                          // If level==0, then the neighbour node is present in the particle
                          // If level==1, then the reference node is present in the particle
                          size_t par_id;
                          if (level == 0) {
                             par_id = (size_t)node.parID[comp]->item(node_neigh(dir,off));
                          } else if ( level == 1) {
                             par_id = (size_t)node.parID[comp]->item(p);
                          }
                          double field_value = impose_solid_temperature(comp,dir,off,i,j,k,xb,par_id);
                          b_intersect[dir].field[comp]->set_item(p,off,field_value);
		       }
                    } 
                 }
              }
           }
        }
     }
  }
}

//---------------------------------------------------------------------------
double
DDS_HeatTransfer:: find_intersection ( size_t const& left, size_t const& right, size_t const& yconst, size_t const& zconst, size_t const& comp, size_t const& dir, size_t const& off, size_t const& level)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_HeatTransfer:: find_intersection" ) ;

  NodeProp node = GLOBAL_EQ->get_node_property();

  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);
  size_t_vector side(2,0);

  side(0) = left;
  side(1) = right;

  for (size_t l=0;l<dim;++l) {
     min_unknown_index(l) = TF->get_min_index_unknown_handled_by_proc( comp, l ) - 1;
     max_unknown_index(l) = TF->get_max_index_unknown_handled_by_proc( comp, l ) + 1;
  }

  double funl=0., func=0., funr=0.;

  double xleft = TF->get_DOF_coordinate( side(0), comp, dir ) ;
  double xright = TF->get_DOF_coordinate( side(1), comp, dir ) ;

  double yvalue=0.,zvalue=0.;
  size_t p=0;

  if (dir == 0) {
     yvalue = TF->get_DOF_coordinate( yconst, comp, 1 ) ;
     if (dim == 3) zvalue = TF->get_DOF_coordinate( zconst, comp, 2 ) ;
     // If level==0, then the neighbour node is present in the particle
     // If level==1, then the reference node is present in the particle
//     p = return_node_index(TF,comp,side(off),yconst,zconst);
     if (off != level) {
        p = return_node_index(TF,comp,side(1),yconst,zconst);
     } else if (off == level) {
        p = return_node_index(TF,comp,side(0),yconst,zconst);
     }
  } else if (dir == 1) {
     yvalue = TF->get_DOF_coordinate( yconst, comp, 0 ) ;
     if (dim == 3) zvalue = TF->get_DOF_coordinate( zconst, comp, 2 ) ;
     // If level==0, then the neighbour node is present in the particle
     // If level==1, then the reference node is present in the particle
//     p = return_node_index(TF,comp,yconst,side(off),zconst);
     if (off != level) {
        p = return_node_index(TF,comp,yconst,side(1),zconst);
     } else if ( off == level) {
        p = return_node_index(TF,comp,yconst,side(0),zconst);
     }
  } else if (dir == 2) {
     yvalue = TF->get_DOF_coordinate( yconst, comp, 0 ) ;
     if (dim == 3) zvalue = TF->get_DOF_coordinate( zconst, comp, 1 ) ;
     // If level==0, then the neighbour node is present in the particle
     // If level==1, then the reference node is present in the particle
//     p = return_node_index(TF,comp,yconst,zconst,side(off));
     if (off != level) {
        p = return_node_index(TF,comp,yconst,zconst,side(1));
     } else if ( off == level) {
        p = return_node_index(TF,comp,yconst,zconst,side(0));
     }
  }

  size_t id = (size_t)node.parID[comp]->item(p);

  double xcenter;

  if (dir == 0) {
     funl = level_set_function(id,comp,xleft,yvalue,zvalue,level_set_type);
     funr = level_set_function(id,comp,xright,yvalue,zvalue,level_set_type);
  } else if (dir == 1) {
     funl = level_set_function(id,comp,yvalue,xleft,zvalue,level_set_type);
     funr = level_set_function(id,comp,yvalue,xright,zvalue,level_set_type);
  } else if (dir == 2) {
     funl = level_set_function(id,comp,yvalue,zvalue,xleft,level_set_type);
     funr = level_set_function(id,comp,yvalue,zvalue,xright,level_set_type);
  }

  // In case both the points are on the same side of solid interface
  // This will occur when the point just outside the solid interface will be considered inside the solid
  // This condition enables the intersection with the interface using the point in fluid and the ACTUAL node in the solid 
  // by shifting the point by 5% of grid size 
  if (funl*funr > 0.) {
     double dx = TF->get_cell_size(side(off),comp,dir) ;
     if (off == level) {
        xleft = xleft - 0.05*dx;
     } else {
        xright = xright + 0.05*dx; 
     }
  }

  if (dir == 0) {
     funl = level_set_function(id,comp,xleft,yvalue,zvalue,level_set_type);
     funr = level_set_function(id,comp,xright,yvalue,zvalue,level_set_type);
  } else if (dir == 1) {
     funl = level_set_function(id,comp,yvalue,xleft,zvalue,level_set_type);
     funr = level_set_function(id,comp,yvalue,xright,zvalue,level_set_type);
  } else if (dir == 2) {
     funl = level_set_function(id,comp,yvalue,zvalue,xleft,level_set_type);
     funr = level_set_function(id,comp,yvalue,zvalue,xright,level_set_type);
  }

  // If the shifted point is also physically outside the solid then xb = dx
  if (funl*funr > 0.) {
     xcenter = TF->get_DOF_coordinate( side(off), comp, dir ) ;
  } else {
     // Bisection method algorithm
     while (MAC::abs(xright-xleft) > 1.E-14) {
        xcenter = (xleft+xright)/2.;
        if (dir == 0) {
           funl = level_set_function(id,comp,xleft,yvalue,zvalue,level_set_type);
           func = level_set_function(id,comp,xcenter,yvalue,zvalue,level_set_type);
        } else if (dir == 1) {
           funl = level_set_function(id,comp,yvalue,xleft,zvalue,level_set_type);
           func = level_set_function(id,comp,yvalue,xcenter,zvalue,level_set_type);
        } else if (dir == 2) {
           funl = level_set_function(id,comp,yvalue,zvalue,xleft,level_set_type);
           func = level_set_function(id,comp,yvalue,zvalue,xcenter,level_set_type);
        }

        if ((func == 1.E-16) || ((xcenter-xleft)/2. <= 1.E-16)) break;

        if (func*funl >= 1.E-16) {
           xleft = xcenter;
        } else {
           xright = xcenter;
        }
     }
  }

  if (off == 0) {
     xcenter = MAC::abs(xcenter - TF->get_DOF_coordinate( side(1), comp, dir ));
  } else if (off == 1) {
     xcenter = MAC::abs(xcenter - TF->get_DOF_coordinate( side(0), comp, dir ));
  }

  return (xcenter);
}

//---------------------------------------------------------------------------
double
DDS_HeatTransfer:: find_intersection_for_ghost ( double const& xl, double const& xr, double const& yvalue, double const& zvalue, size_t const& id, size_t const& comp, size_t const& dir, double const& dx, size_t const& level, size_t const& off)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_HeatTransfer:: find_intersection_for_ghost" ) ;

  doubleVector side(2,0);

  double xleft = xl;
  double xright = xr;

  side(0) = xleft;
  side(1) = xright;

  double funl=0., func=0., funr=0.;

  double xcenter;

  if (dir == 0) {
     funl = level_set_function(id,comp,xleft,yvalue,zvalue,level_set_type);
     funr = level_set_function(id,comp,xright,yvalue,zvalue,level_set_type);
  } else if (dir == 1) {
     funl = level_set_function(id,comp,yvalue,xleft,zvalue,level_set_type);
     funr = level_set_function(id,comp,yvalue,xright,zvalue,level_set_type);
  } else if (dir == 2) {
     funl = level_set_function(id,comp,yvalue,zvalue,xleft,level_set_type);
     funr = level_set_function(id,comp,yvalue,zvalue,xright,level_set_type);
  }

  // In case both the points are on the same side of solid interface
  // This will occur when the point just outside the solid interface will be considered inside the solid
  // This condition enables the intersection with the interface using the point in fluid and the ACTUAL node in the solid 
  // by shifting the point by 5% of grid size 
  if (funl*funr > 0.) {
     if (off == level) {
        xleft = xleft - 0.05*dx;
     } else {
        xright = xright + 0.05*dx;
     }
  }

  if (dir == 0) {
     funl = level_set_function(id,comp,xleft,yvalue,zvalue,level_set_type);
     funr = level_set_function(id,comp,xright,yvalue,zvalue,level_set_type);
  } else if (dir == 1) {
     funl = level_set_function(id,comp,yvalue,xleft,zvalue,level_set_type);
     funr = level_set_function(id,comp,yvalue,xright,zvalue,level_set_type);
  } else if (dir == 2) {
     funl = level_set_function(id,comp,yvalue,zvalue,xleft,level_set_type);
     funr = level_set_function(id,comp,yvalue,zvalue,xright,level_set_type);
  }

  // If the shifted point is also physically outside the solid then xb = dx
  if (funl*funr > 0.) {
     xcenter = side(off) ;
  } else {
     // Bisection method algorithm
     while (MAC::abs(xright-xleft) > 1.E-14) {
        xcenter = (xleft+xright)/2.;
        if (dir == 0) {
           funl = level_set_function(id,comp,xleft,yvalue,zvalue,level_set_type);
           func = level_set_function(id,comp,xcenter,yvalue,zvalue,level_set_type);
        } else if (dir == 1) {
           funl = level_set_function(id,comp,yvalue,xleft,zvalue,level_set_type);
           func = level_set_function(id,comp,yvalue,xcenter,zvalue,level_set_type);
        } else if (dir == 2) {
           funl = level_set_function(id,comp,yvalue,zvalue,xleft,level_set_type);
           func = level_set_function(id,comp,yvalue,zvalue,xcenter,level_set_type);
        }

        if ((func == 1.E-16) || ((xcenter-xleft)/2. <= 1.E-16)) break;

        if (func*funl >= 1.E-16) {
           xleft = xcenter;
        } else {
           xright = xcenter;
        }
     }
  }

  if (off == 0) {
     xcenter = MAC::abs(xcenter - side(1));
  } else if (off == 1) {
     xcenter = MAC::abs(xcenter - side(0));
  }

  return (xcenter);
}

//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: HeatEquation_DirectionSplittingSolver ( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_HeatTransfer:: HeatEquation_DirectionSplittingSolver" ) ;

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
  Solve_i_in_jk (t_it,gamma,0,1,2);
  if ( my_rank == is_master ) SCT_get_elapsed_time("Solver x solution");
  if ( my_rank == is_master ) SCT_set_start("Transfer x solution");
  // Synchronize the distributed DS solution vector
  GLOBAL_EQ->synchronize_DS_solution_vec();
  // Tranfer back to field
  TF->update_free_DOFs_value( 3, GLOBAL_EQ->get_solution_DS_temperature() ) ;
  if (is_solids) nodes_temperature_initialization(3);
  if ( my_rank == is_master ) SCT_get_elapsed_time("Transfer x solution");

  if ( my_rank == is_master ) SCT_set_start("Solver y solution");
  // Solve y-direction(i.e. 1) in x(i.e. 0) and z(i.e. 2)
  Solve_i_in_jk (t_it,gamma,1,0,2);
  if ( my_rank == is_master ) SCT_get_elapsed_time("Solver y solution");
  if ( my_rank == is_master ) SCT_set_start("Transfer y solution");
  // Synchronize the distributed DS solution vector
  GLOBAL_EQ->synchronize_DS_solution_vec();
  // Tranfer back to field
  if (dim == 2) {
     TF->update_free_DOFs_value( 0 , GLOBAL_EQ->get_solution_DS_temperature() ) ;
     if (is_solids) nodes_temperature_initialization(0);
  } else if (dim == 3) {
     TF->update_free_DOFs_value( 4 , GLOBAL_EQ->get_solution_DS_temperature() ) ;
     if (is_solids) nodes_temperature_initialization(4);
  }
  if ( my_rank == is_master ) SCT_get_elapsed_time("Transfer y solution");

  if (dim == 3) {
     if ( my_rank == is_master ) SCT_set_start("Solver z solution");
     // Solve z-direction(i.e. 2) in x(i.e. 0) and y(i.e. 1)
     Solve_i_in_jk (t_it,gamma,2,0,1);
     if ( my_rank == is_master ) SCT_get_elapsed_time("Solver z solution");
     if ( my_rank == is_master ) SCT_set_start("Transfer z solution");
     // Synchronize the distributed DS solution vector
     GLOBAL_EQ->synchronize_DS_solution_vec();
     // Tranfer back to field
     TF->update_free_DOFs_value( 0 , GLOBAL_EQ->get_solution_DS_temperature() ) ;
     if (is_solids) nodes_temperature_initialization(0);
     if ( my_rank == is_master ) SCT_get_elapsed_time("Transfer z solution");
  }
}

//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: output_l2norm ( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_HeatTransfer:: output_l2norm" ) ;

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

      computed_DS_L2 = pelCOMM->sum( computed_DS_L2 ) ;
      computed_DS_L2 = MAC::sqrt(computed_DS_L2);

      if (my_rank == is_master) {
         cout << "L2 Norm for component "<<comp<<" with DS = " << std::fixed << std::setprecision(10) << computed_DS_L2<<endl;
      }
   }
}

//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: DS_error_with_analytical_solution ( FV_DiscreteField const* FF,FV_DiscreteField* FF_ERROR )
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_HeatTransfer:: DS_error_with_analytical_solution" ) ;

   // Parameters
   size_t i,j,k;

   size_t_vector min_unknown_index(dim,0);
   size_t_vector max_unknown_index(dim,0);

   for (size_t comp=0;comp<nb_comps;comp++) {

      double x,y,z,analytical_solution=0., computed_DS_field=0.,error_L2 = 0. ,denom = 0., norm_error = 0.;
      // Get local min and max indices
      for (size_t l=0;l<dim;++l) {
	 min_unknown_index(l) = FF->get_min_index_unknown_handled_by_proc( comp, l ) ;
         max_unknown_index(l) = FF->get_max_index_unknown_handled_by_proc( comp, l ) ;
      }

      for (i=min_unknown_index(0);i<=max_unknown_index(0);++i) {
         x = FF->get_DOF_coordinate( i, comp, 0 ) ;
         for (j=min_unknown_index(1);j<=max_unknown_index(1);++j) {
            y = FF->get_DOF_coordinate( j, comp, 1 ) ;
            if (dim == 2) {
               k=0;

               computed_DS_field = FF->DOF_value( i, j, k, comp, 0 ) ;

               // Presence of solids; analytical solution of hollow cylindrical shell for 2D
               if (is_solids) {
                  PartInput solid = GLOBAL_EQ->get_solid();
                  NodeProp node = GLOBAL_EQ->get_node_property();
                  size_t p = return_node_index(FF,comp,i,j,0);
                  if (node.void_frac[comp]->item(p) != 1.) {
                     double xp = solid.coord[comp]->item(0,0);
                     double yp = solid.coord[comp]->item(0,1);
                     double r1 = solid.size[comp]->item(0);
                     double t1 = solid.temp[comp]->item(0);
                     double r2 = solid.size[comp]->item(1);
                     double t2 = solid.temp[comp]->item(1);
                     double r = pow(pow(x-xp,2.)+pow(y-yp,2.),0.5);
                     analytical_solution = t2 - (t2-t1)*(log(r/r2)/log(r1/r2));
                     FF_ERROR->set_DOF_value( i, j, k, comp, 0, MAC::abs( computed_DS_field - analytical_solution) ) ;
                     if ( FF->DOF_is_unknown_handled_by_proc( i, j, k, comp ) ) {
                        error_L2 += MAC::sqr( computed_DS_field - analytical_solution)
                                               * FF->get_cell_measure( i, j, k, comp ) ;
                     }
                     denom += MAC::sqr( analytical_solution) * FF->get_cell_measure( i, j, k, comp ) ;
                  }
               } else {
                  analytical_solution = MAC::sin( MAC::pi() * x )* MAC::sin( MAC::pi() * y ) ;
                  //Get error between computed solutions with Direction splitting
                  FF_ERROR->set_DOF_value( i, j, k, comp, 0, MAC::abs( computed_DS_field - analytical_solution) ) ;
                  if ( FF->DOF_is_unknown_handled_by_proc( i, j, k, comp ) ) {
                     error_L2 += MAC::sqr( computed_DS_field - analytical_solution)
                                            * FF->get_cell_measure( i, j, k, comp ) ;
                  }
                  denom += MAC::sqr( analytical_solution) * FF->get_cell_measure( i, j, k, comp ) ;
               }
            } else {
               for (k=min_unknown_index(2);k<=max_unknown_index(2);++k) {
                  z = FF->get_DOF_coordinate( k, comp, 2 ) ;
                  computed_DS_field = FF->DOF_value( i, j, k, comp, 0 ) ;
                  // Presence of solids; analytical solution for the hollow spherical shell in 3D
                  if (is_solids) {
                     PartInput solid = GLOBAL_EQ->get_solid();
                     NodeProp node = GLOBAL_EQ->get_node_property();
                     size_t p = return_node_index(FF,comp,i,j,k);
                     if (node.void_frac[comp]->item(p) != 1.) {
                        double xp = solid.coord[comp]->item(0,0);
                        double yp = solid.coord[comp]->item(0,1);
                        double zp = solid.coord[comp]->item(0,2);
                        double r1 = solid.size[comp]->item(0);
                        double t1 = solid.temp[comp]->item(0);
                        double r2 = solid.size[comp]->item(1);
                        double t2 = solid.temp[comp]->item(1);
                        double r = pow(pow(x-xp,2.)+pow(y-yp,2.)+pow(z-zp,2.),0.5);
                        analytical_solution = t1 - (t1-t2)*(1-r1/r)/(1-r1/r2);
                        FF_ERROR->set_DOF_value( i, j, k, comp, 0, MAC::abs( computed_DS_field - analytical_solution) ) ;
                        if ( FF->DOF_is_unknown_handled_by_proc( i, j, k, comp ) ) {
                           error_L2 += MAC::sqr( computed_DS_field - analytical_solution)
                                                  * FF->get_cell_measure( i, j, k, comp ) ;
                        }
                        denom += MAC::sqr( analytical_solution) * FF->get_cell_measure( i, j, k, comp ) ;
                     }
                  } else {
                     analytical_solution = MAC::sin( MAC::pi() * x )* MAC::sin( MAC::pi() * y ) * MAC::sin( MAC::pi() * z ) ;
                     //Get error between computed solutions with Direction splitting
                     FF_ERROR->set_DOF_value( i, j, k, comp, 0, MAC::abs( computed_DS_field- analytical_solution ) ) ;
                     if ( FF->DOF_is_unknown_handled_by_proc( i, j, k, comp ) ) {
                        error_L2 += MAC::sqr( computed_DS_field - analytical_solution)
                                               * FF->get_cell_measure( i, j, k, comp ) ;
	             }
                     denom += MAC::sqr( analytical_solution) * FF->get_cell_measure( i, j, k, comp ) ;
                  }
               }
            }
         }
      }

      error_L2 = pelCOMM->sum( error_L2 ) ;
      error_L2 = MAC::sqrt(error_L2);
      denom = pelCOMM->sum( denom ) ;
      denom = MAC::sqrt(denom);
      norm_error = error_L2/denom;

      if (my_rank == is_master) {
         cout << "L2 Error with analytical solution for "<<comp<<" DS (Normalized) = " << norm_error<<endl;
         cout << "L2 Error with analytical solution for "<<comp<<" DS = " << error_L2<<endl;
      }
   }
}




//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: create_DDS_subcommunicators ( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransfer:: create_DDS_subcommunicators" ) ;

   int color = 0, key = 0;
   //int const* number_of_subdomains_per_direction = TF->primary_grid()->get_domain_decomposition() ;
   int const* MPI_coordinates_world = TF->primary_grid()->get_MPI_coordinates() ;
   int const* MPI_number_of_coordinates = TF->primary_grid()->get_domain_decomposition() ;

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
      color = MPI_coordinates_world[1] + MPI_coordinates_world[2]*MPI_number_of_coordinates[1] ;
      key = MPI_coordinates_world[0];
      // Split by direction in x
      processor_splitting (color,key,0);

      // Assign color and key for splitting in y
      color = MPI_coordinates_world[2] + MPI_coordinates_world[0]*MPI_number_of_coordinates[2];
      key = MPI_coordinates_world[1];
      // Split by direction in y
      processor_splitting (color,key,1);

      // Assign color and key for splitting in y
      color = MPI_coordinates_world[0] + MPI_coordinates_world[1]*MPI_number_of_coordinates[0];;
      key = MPI_coordinates_world[2];

      // Split by direction in y
      processor_splitting (color,key,2);
   }
}

//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: processor_splitting ( int const& color, int const& key, size_t const& dir )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransfer:: processor_splitting" ) ;

   MPI_Comm_split(MPI_COMM_WORLD, color, key, &DDS_Comm_i[dir]);
   MPI_Comm_size( DDS_Comm_i[dir], &nb_ranks_comm_i[dir] ) ;
   MPI_Comm_rank( DDS_Comm_i[dir], &rank_in_i[dir] ) ;

}

//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: allocate_mpi_variables ( void )
//---------------------------------------------------------------------------
{

   for (size_t dir = 0; dir < dim; dir++) {  
      first_pass[dir].size = new size_t [nb_comps];
      second_pass[dir].size = new size_t [nb_comps];
      for (size_t comp = 0; comp < nb_comps; comp++) {
         size_t local_min_j=0, local_max_j=0;
         size_t local_min_k=0, local_max_k=0;

         // Get local min and max indices
         size_t_vector min_unknown_index(dim,0);
         size_t_vector max_unknown_index(dim,0);
         for (size_t l=0;l<dim;++l) {
            min_unknown_index(l) = TF->get_min_index_unknown_handled_by_proc( comp, l ) ;
            max_unknown_index(l) = TF->get_max_index_unknown_handled_by_proc( comp, l ) ;
         }

         if (dir == 0) {
            local_min_j = min_unknown_index(1);
            local_max_j = max_unknown_index(1);
            if (dim == 3) {
               local_min_k = min_unknown_index(2);
               local_max_k = max_unknown_index(2);
            }
         } else if (dir == 1) {
            local_min_j = min_unknown_index(0);
            local_max_j = max_unknown_index(0);
            if (dim == 3) {
               local_min_k = min_unknown_index(2);
               local_max_k = max_unknown_index(2);
            }
         } else if (dir == 2) {
            local_min_j = min_unknown_index(0);
            local_max_j = max_unknown_index(0);
            local_min_k = min_unknown_index(1);
            local_max_k = max_unknown_index(1);
         }

         size_t local_length_j = (local_max_j-local_min_j+1);
         size_t local_length_k = (local_max_k-local_min_k+1);

         if (dim != 3) {
            first_pass[dir].size[comp] = 3*local_length_j;
            second_pass[dir].size[comp] = 2*local_length_j;
         } else if (dim == 3) {
            first_pass[dir].size[comp] = 3*local_length_j*local_length_k;
            second_pass[dir].size[comp] = 2*local_length_j*local_length_k;
         }
      }
   }

   // Array declarations
   for (size_t dir = 0; dir < dim; dir++) {  
      first_pass[dir].send = new double** [nb_comps];
      first_pass[dir].receive = new double** [nb_comps];
      second_pass[dir].send = new double** [nb_comps];
      second_pass[dir].receive = new double** [nb_comps];
      for (size_t comp = 0; comp < nb_comps; comp++) {
         first_pass[dir].send[comp] = new double* [nb_ranks_comm_i[dir]];
         first_pass[dir].receive[comp] = new double* [nb_ranks_comm_i[dir]];
         second_pass[dir].send[comp] = new double* [nb_ranks_comm_i[dir]];
         second_pass[dir].receive[comp] = new double* [nb_ranks_comm_i[dir]];
         for (size_t i = 0; i < (size_t)nb_ranks_comm_i[dir]; i++) {
            first_pass[dir].send[comp][i] = new double[first_pass[dir].size[comp]];
            first_pass[dir].receive[comp][i] = new double[first_pass[dir].size[comp]];
            second_pass[dir].send[comp][i] = new double[second_pass[dir].size[comp]];
            second_pass[dir].receive[comp][i] = new double[second_pass[dir].size[comp]];
         }
      }
   }
}

//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: deallocate_mpi_variables ( void )
//---------------------------------------------------------------------------
{
   // Array declarations
   for (size_t dir = 0; dir < dim; dir++) {
      for (size_t comp = 0; comp < nb_comps; comp++) {
         for (size_t i = 0; i < (size_t) nb_ranks_comm_i[dir]; i++) {
            delete [] first_pass[dir].send[comp][i];
            delete [] first_pass[dir].receive[comp][i];
            delete [] second_pass[dir].send[comp][i];
            delete [] second_pass[dir].receive[comp][i];
         }
         delete [] first_pass[dir].send[comp];
         delete [] first_pass[dir].receive[comp];
         delete [] second_pass[dir].send[comp];
         delete [] second_pass[dir].receive[comp];
      }
      delete [] first_pass[dir].send;
      delete [] first_pass[dir].receive;
      delete [] second_pass[dir].send;
      delete [] second_pass[dir].receive;
      delete [] first_pass[dir].size;
      delete [] second_pass[dir].size;
   }
}

//---------------------------------------------------------------------------
void
DDS_HeatTransfer:: free_DDS_subcommunicators ( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransfer:: free_DDS_subcommunicators" ) ;


}

//----------------------------------------------------------------------
double DDS_HeatTransfer:: assemble_advection_TVD( FV_DiscreteField const* AdvectingField, 
	size_t advecting_level, double const& coef, size_t const& i, size_t const& j, size_t const& k, size_t advected_level) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransfer:: assemble_advection_TVD" );   
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
            xR = TF->get_DOF_coordinate( i+1, advected_level, 0 );
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
double DDS_HeatTransfer:: assemble_advection_Centered_new( FV_DiscreteField const* AdvectingField, 
	size_t advecting_level, double const& coef, size_t const& i, size_t const& j, size_t const& k, size_t advected_level) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransfer:: assemble_advection_Centered" );   
   MAC_CHECK_PRE( advecting_level < AdvectingField->storage_depth() ) ;
   MAC_ASSERT( AdvectingField->discretization_type() == "staggered" ) ; 
   
   // Parameters
   size_t component = 0 ;  
   double dxC = 0., dyC = 0., dzC = 0.; 

   double AdvectedvalueC = 0., AdvectedvalueRi = 0., 
   	AdvectedvalueLe = 0., AdvectedvalueTo = 0., 
	AdvectedvalueBo = 0., AdvectedvalueFr = 0., 
	AdvectedvalueBe = 0.;
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
	 	 
      ur = AdvectingField->DOF_value( i+shift.i, j, k, 0, advecting_level );
      ul = AdvectingField->DOF_value( i+shift.i-1, j, k, 0, advecting_level );
      // Right (X)

      if ( TF->DOF_color( i+1, j, k, component ) == FV_BC_RIGHT ) {
         fri = 0.5*(ur+ul) * AdvectedvalueRi;
      } else {
         fri = 0.5*(ur+ul) * 0.5*(AdvectedvalueRi + AdvectedvalueC);
      }
	 
      // Left (X)
	   
      if ( TF->DOF_color( i-1, j, k, component ) == FV_BC_LEFT ) {
         fle = 0.5*(ur+ul) * AdvectedvalueLe;
      } else {
         fle = 0.5*(ur+ul) * 0.5*(AdvectedvalueLe + AdvectedvalueC);
      }
	 
      // Top and Bottom
      // --------------
      AdvectedvalueTo = TF->DOF_value( i, j+1, k, component, advected_level );
      AdvectedvalueBo = TF->DOF_value( i, j-1, k, component, advected_level );

      vt = AdvectingField->DOF_value( i, j+shift.j, k, 1, advecting_level );
      vb = AdvectingField->DOF_value( i, j+shift.j-1, k, 1, advecting_level );
      // Top (Y)
	   
      if ( TF->DOF_color( i, j+1, k, component ) == FV_BC_TOP ) {
         fto = 0.5*(vt+vb) * AdvectedvalueTo;
      } else {
         fto = 0.5*(vt+vb) * 0.5*(AdvectedvalueTo + AdvectedvalueC);
      }

      // Bottom (Y)
	   
      if ( TF->DOF_color( i, j-1, k, component ) == FV_BC_BOTTOM ) {
         fbo = 0.5*(vt+vb) * AdvectedvalueBo;
      } else {
         fbo = 0.5*(vt+vb) * 0.5*(AdvectedvalueBo + AdvectedvalueC);
      }

      flux = ( fto - fbo ) * dxC + ( fri - fle ) * dyC;	 

   } else {
      dzC = TF->get_cell_size( k, component, 2 ) ;    

      AdvectedvalueC = TF->DOF_value( i, j, k, component, advected_level );
	 
      // Right and Left
      // --------------
      AdvectedvalueRi = TF->DOF_value( i+1, j, k, component, advected_level );
      AdvectedvalueLe = TF->DOF_value( i-1, j, k, component, advected_level );
	 	 
      ul = AdvectingField->DOF_value( i+shift.i-1, j, k, 0, advecting_level );
      ur = AdvectingField->DOF_value( i+shift.i, j, k, 0, advecting_level );
      // Right (X)

      if ( TF->DOF_color( i+1, j, k, component ) == FV_BC_RIGHT ) {
         fri = 0.5*(ur+ul) * AdvectedvalueRi;
      } else {
         fri = 0.5*(ur+ul) * 0.5*(AdvectedvalueRi + AdvectedvalueC);
      }
	 
      // Left (X)
	   
      if ( TF->DOF_color( i-1, j, k, component ) == FV_BC_LEFT ) {
         fle = 0.5*(ur+ul) * AdvectedvalueLe;
      } else {
         fle = 0.5*(ur+ul) * 0.5*(AdvectedvalueLe + AdvectedvalueC);
      }
	 
      // Top and Bottom
      // --------------
      AdvectedvalueTo = TF->DOF_value( i, j+1, k, component, advected_level );
      AdvectedvalueBo = TF->DOF_value( i, j-1, k, component, advected_level );

      vb = AdvectingField->DOF_value( i, j+shift.j-1, k, 1, advecting_level );
      vt = AdvectingField->DOF_value( i, j+shift.j, k, 1, advecting_level );
      // Top (Y)
	   
      if ( TF->DOF_color( i, j+1, k, component ) == FV_BC_TOP ) {
         fto = 0.5*(vt+vb) * AdvectedvalueTo;
      } else {   
         fto = 0.5*(vt+vb) * 0.5*(AdvectedvalueTo + AdvectedvalueC);
      }

      // Bottom (Y)
	   
      if ( TF->DOF_color( i, j-1, k, component ) == FV_BC_BOTTOM ) {
         fbo = 0.5*(vt+vb) * AdvectedvalueBo;
      } else {
         fbo = 0.5*(vt+vb) * 0.5*(AdvectedvalueBo + AdvectedvalueC);
      }
	 
      // Front and Behind
      // ----------------
      AdvectedvalueFr = TF->DOF_value( i, j, k+1, component, advected_level );
      AdvectedvalueBe = TF->DOF_value( i, j, k-1, component, advected_level );

      wf = AdvectingField->DOF_value( i, j, k+shift.k, 2, advecting_level );
      wb = AdvectingField->DOF_value( i, j, k+shift.k-1, 2, advecting_level );
      // Front (Z)
	   
      if ( TF->DOF_color( i, j, k+1, component ) == FV_BC_FRONT ) {
         ffr = 0.5*(wf+wb) * AdvectedvalueFr;
      } else {
         ffr = 0.5*(wf+wb) * 0.5*(AdvectedvalueFr + AdvectedvalueC);
      }

      // Behind (Z)
	   
      if ( TF->DOF_color( i, j, k-1, component ) == FV_BC_BEHIND ) {
         fbe = 0.5*(wf+wb) * AdvectedvalueBe;
      } else {
         fbe = 0.5*(wf+wb) * 0.5*(AdvectedvalueBe + AdvectedvalueC);
      }

      flux = ( fto - fbo ) * dxC * dzC + ( fri - fle ) * dyC * dzC + ( ffr - fbe ) * dxC * dyC;
   } 
   
   return (coef * flux);
         
}


//----------------------------------------------------------------------
double DDS_HeatTransfer:: assemble_advection_Centered( FV_DiscreteField const* AdvectingField, 
	size_t advecting_level, double const& coef, size_t const& i, size_t const& j, size_t const& k, size_t advected_level) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransfer:: assemble_advection_Centered" );   
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
double DDS_HeatTransfer:: assemble_advection_Upwind( FV_DiscreteField const* AdvectingField,
	size_t advecting_level, double const& coef, size_t const& i, size_t const& j, size_t const& k, size_t advected_level) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransfer:: assemble_advection_Upwind" );
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
