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
   b_restart = MAC_Application::is_follow();
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

   if (is_solids) {
      Npart = fromNS.Npart_ ;
      solid_filename = fromNS.solid_filename_ ;
      loc_thres = fromNS.loc_thres_ ;
      level_set_type = fromNS.level_set_type_ ; 
   }

   // Build the matrix system
   MAC_ModuleExplorer* se =
	exp->create_subexplorer( 0,"DDS_HeatTransferSystem" ) ;
   GLOBAL_EQ = DDS_HeatTransferSystem::create( this, se, TF ) ;
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

   // Generate solid particles if required
   if (is_solids) {
      Solids_generation();
      node_property_calculation();
      nodes_temperature_initialization(0);
      nodes_temperature_initialization(1);
      nodes_temperature_initialization(3);
      if (dim == 3) nodes_temperature_initialization(4);
   }

   GLOBAL_EQ->initialize_temperature();

   // Necessary especially for cases with non-zero field initiallization
   if (b_restart == false) ugradu_initialization ( );

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

   // Perform matrix level operations before each time step
   GLOBAL_EQ->at_each_time_step( );

   // Assemble 1D tridiagonal matrices and schur complement calculation
   assemble_temperature_and_schur(t_it);

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
  os2 << "/home/goyal001/Documents/Computing/MAC-Test/DS_results/output_" << my_rank << ".csv";
  std::string filename = os2.str();
  outputFile.open(filename.c_str());


//  outputFile.open("/home/goyal001/Documents/Computing/MAC-Test/DS_results/output_" << my_rank << ".txt", std::ios_base::out | std::ios_base::trunc) ;
  size_t i,j,k;
  outputFile << "x,y,z,par_ID,void_frac,left,lv,right,rv,bottom,bov,top,tv,behind,bev,front,fv" << endl;

  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);

  PartInput solid = GLOBAL_EQ->get_solid();
  NodeProp node = GLOBAL_EQ->get_node_property();
  BoundaryBisec* b_intersect = GLOBAL_EQ->get_b_intersect(0);

  for (size_t comp=0;comp<nb_comps;comp++) {
     // Get local min and max indices
     for (size_t l=0;l<dim;++l) {
        min_unknown_index(l) = TF->get_min_index_unknown_handled_by_proc( comp, l ) - 1;
        max_unknown_index(l) = TF->get_max_index_unknown_handled_by_proc( comp, l ) + 1;
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
              size_t id = node.parID[comp]->item(p);
              size_t voidf = node.void_frac[comp]->item(p);

              outputFile << xC << "," << yC << "," << zC << "," << id << "," << voidf;
              for (size_t dir = 0; dir < dim; dir++) {
                  for (size_t off = 0; off < 2; off++) {
                      outputFile << "," << b_intersect[dir].offset[comp]->item(p,off) << "," << b_intersect[dir].value[comp]->item(p,off); 
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
      if (TF->DOF_in_domain( i-1, j, k, comp) && TF->DOF_in_domain( i+1, j, k, comp))
         value = xright/xhr - xleft/xhl;
      else if (TF->DOF_in_domain( i-1, j, k, comp))
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
      if (TF->DOF_in_domain(i, j-1, k, comp) && TF->DOF_in_domain(i, j+1, k, comp))
         value = yright/yhr - yleft/yhl;
      else if(TF->DOF_in_domain(i, j-1, k, comp))
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
      if (TF->DOF_in_domain(i, j, k-1, comp) && TF->DOF_in_domain(i, j, k+1, comp))
         value = zright/zhr - zleft/zhl;
      else if(TF->DOF_in_domain(i, j, k-1, comp))
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
      ugradu = assemble_advection_TVD(UF,1,1.,i,j,k,1);
   } else if ( AdvectionScheme == "Upwind" ) {
      ugradu = assemble_advection_Upwind_new(UF,1,1.,i,j,k,1);
/*   } else if ( AdvectionScheme == "Centered" ) {
      ugradu = assemble_advection_Centered(1,rho,1,i,j,k,comp);*/
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
   size_t_vector min_unknown_index(dim,0);
   size_t_vector max_unknown_index(dim,0);
   size_t_vector i_length(dim,0);
   for (size_t l=0;l<dim;++l) {
	min_unknown_index(l) = FF->get_min_index_unknown_handled_by_proc( comp, l ) - 1;
	max_unknown_index(l) = FF->get_max_index_unknown_handled_by_proc( comp, l ) + 1;
        i_length(l) = 1 + max_unknown_index(l) - min_unknown_index(l);
   }

   size_t local_min_k = 0;
   if (dim == 3) local_min_k = min_unknown_index(2);

   size_t p = (j-min_unknown_index(1)) + i_length(1)*(i-min_unknown_index(0)) + i_length(0)*i_length(1)*(k-local_min_k);

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
       size_t ii=0,jj=0,kk=0;

       // Condition for handling the pressure neumann conditions at wall
       if (i==min_unknown_index(dir) && l_bound) {
          if (dir == 0) {
             ii = i-1; jj = min_unknown_index(1); kk = k_min;
          } else if (dir == 1) {
             ii = min_unknown_index(0); jj = i-1; kk = k_min;
          } else if (dir == 2) {
             ii = min_unknown_index(0); jj = min_unknown_index(1); kk = i-1;
          }
          if (FF->DOF_in_domain(ii,jj,kk,comp) && FF->DOF_has_imposed_Dirichlet_value(ii,jj,kk,comp)) {
             // For Dirichlet boundary condition
             value = center;
          } else {
             // For Neumann homogeneous boundary condition
             value = -right;
          }
       } else if (i==max_unknown_index(dir) && r_bound) {
          if (dir == 0) {
             ii = i+1; jj = max_unknown_index(1); kk = k_max;
          } else if (dir == 1) {
             ii = max_unknown_index(0); jj = i+1; kk = k_max;
          } else if (dir == 2) {
             ii = max_unknown_index(0); jj = max_unknown_index(1); kk = i+1;
          }
          if (FF->DOF_in_domain(ii,jj,kk,comp) && FF->DOF_has_imposed_Dirichlet_value(ii,jj,kk,comp)) {
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
         for (size_t i=1;i<nb_ranks_comm_i[dir];++i) {

             // Create the container to receive
             size_t nbrows = product_matrix->nb_rows();
             size_t nb_received_data = nbrows*nbrows+1;
             double * received_data = new double [nb_received_data];

             // Receive the data
             static MPI_Status status ;
             MPI_Recv( received_data, nb_received_data, MPI_DOUBLE, i, 0,
                                       DDS_Comm_i[dir], &status ) ;

             // Transfer the received data to the receive matrix
             for (size_t k=0;k<nbrows;k++) {
                for (size_t j=0;j<nbrows;j++) {
                   // Assemble the global product matrix by adding contributions from all the procs
                   receive_matrix->add_to_item(k,j,received_data[k*(nbrows)+j]);
                }
             }

   	     if (is_iperiodic[dir] == 0) {
                if (i<nb_ranks_comm_i[dir]-1) {
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
         MPI_Send( packed_data, nb_send_data, MPI_DOUBLE, 0, 0, DDS_Comm_i[dir] ) ;

         delete [] packed_data;

      }

      // Assemble the schlur complement in the master proc

      if (rank_in_i[dir] == 0) {
	 TDMatrix* Schur = GLOBAL_EQ-> get_Schur();
         size_t nb_row = Schur[dir].ii_main[comp][r_index]->nb_rows();
         for (int p = 0; p < nb_row; p++) {
            Schur[dir].ii_main[comp][r_index]->set_item(p,A[dir].ee[comp][r_index]->item(p,p)-receive_matrix->item(p,p));
            if (p < nb_row-1) Schur[dir].ii_super[comp][r_index]->set_item(p,-receive_matrix->item(p,p+1));
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
            size_t nb_row = DoubleSchur[dir].ii_main[comp][r_index]->nb_rows();
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
      for (int p = 0; p < nb_row; p++) {
         Schur[dir].ii_main[comp][r_index]->set_item(p,A[dir].ee[comp][r_index]->item(p,p)-receive_matrix->item(p,p));
         if (p < nb_row-1) Schur[dir].ii_super[comp][r_index]->set_item(p,-receive_matrix->item(p,p+1));
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
   double dC,xC,yC,zC=0;
   double fe=0.;

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
           NodeProp node = GLOBAL_EQ->get_node_property();
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
              NodeProp node = GLOBAL_EQ->get_node_property();
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
              NodeProp node = GLOBAL_EQ->get_node_property();
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
           NodeProp node = GLOBAL_EQ->get_node_property();
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

   if ( TF->DOF_in_domain(ii,jj,kk,comp))
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

   if ( TF->DOF_in_domain(ii,jj,kk,comp))
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
   for (size_t i=1;i<nb_ranks_comm_i[dir];i++) {
      if (i!=nb_ranks_comm_i[dir]-1) {
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
   for (size_t i=1;i<nb_ranks_comm_i[dir];++i) {
      if (i!=nb_ranks_comm_i[dir]-1) {
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
DDS_HeatTransfer:: unpack_ue(size_t const& comp, double * received_data, size_t const& dir, int const& p)
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
         for (size_t i=1;i<nb_ranks_comm_i[dir];++i) {
            static MPI_Status status;
            MPI_Recv( first_pass[dir].receive[comp][i], first_pass[dir].size[comp], MPI_DOUBLE, i, 0,
                                             DDS_Comm_i[dir], &status ) ;
         }
      }

      // Solve system of interface unknowns for each y
      if (dim == 2) {
         size_t k = 0;
         for (size_t j=local_min_j;j<=local_max_j;j++) {

     	    size_t p = j-local_min_j;

	    unpack_compute_ue_pack(comp,dir,p); 

            // Need to have the original rhs function assembled for corrosponding j,k pair
            double fe = assemble_local_rhs(j,k,gamma,t_it,comp,dir);

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
               double fe = assemble_local_rhs(j,k,gamma,t_it,comp,dir);

               // Setup RHS = fi - Aie*ue for solving ui
               A[dir].ie[comp][p]->multiply_vec_then_add(VEC[dir].interface_T[comp],VEC[dir].local_T[comp],-1.0,1.0);

               // Solve ui and transfer solution into distributed vector
               GLOBAL_EQ->DS_HeatEquation_solver(j,k,min_unknown_index(dir),comp,dir,p);
            }
         }
      }
   } else {
      // Send the packed data to master
      MPI_Send( first_pass[dir].send[comp][rank_in_i[dir]], first_pass[dir].size[comp], MPI_DOUBLE, 0, 0, DDS_Comm_i[dir] ) ;
   }

   // Send the data from master iff multi processor are used
   if (nb_ranks_comm_i[dir] != 1) {
      if ( rank_in_i[dir] == 0 ) {
         for (size_t i=1;i<nb_ranks_comm_i[dir];++i) {
            MPI_Send( second_pass[dir].send[comp][i], second_pass[dir].size[comp], MPI_DOUBLE, i, 0, DDS_Comm_i[dir] ) ;
         }
      } else {
         // Create the container to receive the ue
         static MPI_Status status ;
         MPI_Recv( second_pass[dir].receive[comp][rank_in_i[dir]], second_pass[dir].size[comp], MPI_DOUBLE, 0, 0,
                             DDS_Comm_i[dir], &status ) ;

         // Solve the system of equations in each proc

         if (dim == 2) {
            size_t k = 0;
            for (size_t j = local_min_j;j<=local_max_j;j++) {
               size_t p = j-local_min_j;

               unpack_ue(comp,second_pass[dir].receive[comp][rank_in_i[dir]],dir,p);

               // Need to have the original rhs function assembled for corrosponding j,k pair
               double fe = assemble_local_rhs(j,k,gamma,t_it,comp,dir);

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
                  double fe = assemble_local_rhs(j,k,gamma,t_it,comp,dir);

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
              double fe = assemble_local_rhs(j,k,gamma,t_it,comp,dir_i);
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

  // Structure of particle input data
  PartInput solid = GLOBAL_EQ->get_solid();

  double xp,yp,zp,Rp,vx,vy,vz,wx,wy,wz,Tp,off;

  for (size_t comp=0;comp<nb_comps;comp++) {
     ifstream inFile;
     std::ostringstream os2;
     os2 << "./InputFiles/" << solid_filename;
     std::string filename = os2.str();

     inFile.open(filename.c_str());
     string line;
     getline(inFile,line);
     for (size_t i=0;i<Npart;i++) {
        inFile >> xp >> yp >> zp >> Rp >> vx >> vy >> vz >> wx >> wy >> wz >> Tp >> off;
        solid.coord[comp]->set_item(i,0,xp);
        solid.coord[comp]->set_item(i,1,yp);
        solid.coord[comp]->set_item(i,2,zp);
        solid.size[comp]->set_item(i,Rp);
        solid.temp[comp]->set_item(i,Tp);
        solid.inside[comp]->set_item(i,off);
     }
     inFile.close();
  }
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
  PartInput solid = GLOBAL_EQ->get_solid();

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
              size_t p = return_node_index(TF,comp,i,j,k);
              if (node.void_frac[comp]->item(p) == 1.) {
                 size_t par_id = node.parID[comp]->item(p);
                 double Tpart = solid.temp[comp]->item(par_id);
                 TF->set_DOF_value( i, j, k, comp, level,Tpart);
              }
           }
        }
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

  double level_set = 0.;
  // Type 0 is for circular/spherical solids in 2D/3D system
  if (type == "Sphere") {
     level_set = pow(pow(delta(0),2.)+pow(delta(1),2.)+pow(delta(2),2.),0.5)-Rp;
  } else if (type == "PipeX") {
     level_set = pow(pow(delta(1),2.)+pow(delta(2),2.),0.5)-Rp;
  }

  return(level_set);

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
        min_unknown_index(l) = TF->get_min_index_unknown_handled_by_proc( comp, l ) - 1 ;
        max_unknown_index(l) = TF->get_max_index_unknown_handled_by_proc( comp, l ) + 1 ;
     }

     size_t local_min_k = 0;
     size_t local_max_k = 0;

     if (dim == 3) {
        local_min_k = min_unknown_index(2);
        local_max_k = max_unknown_index(2);
     }

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
                    node.parID[comp]->set_item(p,m);
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

  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);
  size_t_vector ipos(3,0);
  size_t_array2D local_unknown_extents(dim,2,0);
  size_t_array2D node_neigh(dim,2,0);

  PartInput solid = GLOBAL_EQ->get_solid();
  NodeProp node = GLOBAL_EQ->get_node_property();
  BoundaryBisec* b_intersect = GLOBAL_EQ->get_b_intersect(level);

  for (size_t l=0;l<dim;++l) {
     min_unknown_index(l) = TF->get_min_index_unknown_handled_by_proc( comp, l ) - 1 ;
     max_unknown_index(l) = TF->get_max_index_unknown_handled_by_proc( comp, l ) + 1 ;
     local_unknown_extents(l,0) = 0;
     local_unknown_extents(l,1) = (max_unknown_index(l)-min_unknown_index(l));
  }

  size_t local_min_k = 0;
  size_t local_max_k = 0;

  if (dim == 3) {
     local_min_k = min_unknown_index(2);
     local_max_k = max_unknown_index(2);
  }

  for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i) {
     ipos(0) = i - min_unknown_index(0);
     for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j) {
        ipos(1) = j - min_unknown_index(1);
        for (size_t k=local_min_k;k<=local_max_k;++k) {
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
                    size_t left, right;
                    if (off == 0) {
                       left = ii-1; right = ii;
                    } else if (off == 1) {
                       left = ii; right = ii+1;
                    }

                    if ((node.void_frac[comp]->item(node_neigh(dir,off)) != node.void_frac[comp]->item(p)) && (ipos(dir) != local_unknown_extents(dir,off))) {
                       double xb = find_intersection(left,right,jj,kk,comp,dir,off);
                       b_intersect[dir].offset[comp]->set_item(p,off,1);
                       b_intersect[dir].value[comp]->set_item(p,off,xb);
                       size_t par_id = node.parID[comp]->item(node_neigh(dir,off));
                       double field_value = solid.temp[comp]->item(par_id);
                       b_intersect[dir].field[comp]->set_item(p,off,field_value);
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
DDS_HeatTransfer:: find_intersection ( size_t const& left, size_t const& right, size_t const& yconst, size_t const& zconst, size_t const& comp, size_t const& dir, size_t const& off)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_HeatTransfer:: find_intersection" ) ;

  PartInput solid = GLOBAL_EQ->get_solid();
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
  size_t p;

  if (dir == 0) {
     yvalue = TF->get_DOF_coordinate( yconst, comp, 1 ) ;
     if (dim == 3) zvalue = TF->get_DOF_coordinate( zconst, comp, 2 ) ;
     p = return_node_index(TF,comp,side(off),yconst,zconst);
  } else if (dir == 1) {
     yvalue = TF->get_DOF_coordinate( yconst, comp, 0 ) ;
     if (dim == 3) zvalue = TF->get_DOF_coordinate( zconst, comp, 2 ) ;
     p = return_node_index(TF,comp,yconst,side(off),zconst);
  } else if (dir == 2) {
     yvalue = TF->get_DOF_coordinate( yconst, comp, 0 ) ;
     if (dim == 3) zvalue = TF->get_DOF_coordinate( zconst, comp, 1 ) ;
     p = return_node_index(TF,comp,yconst,zconst,side(off));
  }

  size_t id = node.parID[comp]->item(p);

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
  if (funl*funr > 0.) {
     double dx = TF->get_cell_size(side(off),comp,dir) ;
     if (off == 0) {
        xleft = xleft - dx;
     } else {
        xright = xright + dx; 
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
   size_t local_number=0,k=0;

   LA_Vector const* DS_Solution = GLOBAL_EQ->get_solution_DS_temperature();

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
               local_number = TF->DOF_local_number( i, j, k, comp );
               computed_DS_field = TF->DOF_value( i, j, k, comp, 0 ) ;

               if ( TF->DOF_is_unknown_handled_by_proc( i, j, k, comp ) ) {
                  computed_DS_L2 += computed_DS_field*computed_DS_field
                               * TF->get_cell_measure( i, j, k, comp ) ;
               }
            } else {
               for (k=min_unknown_index(2);k<=max_unknown_index(2);++k) {
                  local_number = TF->DOF_local_number( i, j, k, comp );
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
      first_pass[dir].size = new int [nb_comps];
      second_pass[dir].size = new int [nb_comps];
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
         for (size_t i = 0; i < nb_ranks_comm_i[dir]; i++) {
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
         for (size_t i = 0; i < nb_ranks_comm_i[dir]; i++) {
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
   size_t center_pos_in_matrix = 0, component = 0 ;  
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
      size_t k = 0 ;
      center_pos_in_matrix = TF->DOF_global_number( i, j, k, component );
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

      center_pos_in_matrix = TF->DOF_global_number( i, j, k, component );
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
double DDS_HeatTransfer:: assemble_advection_Upwind( FV_DiscreteField const* AdvectingField,
	size_t advecting_level, double const& coef, size_t const& i, size_t const& j, size_t const& k, size_t advected_level) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransfer:: assemble_advection_Upwind" );
   MAC_CHECK_PRE( advecting_level < AdvectingField->storage_depth() ) ;
   MAC_ASSERT( AdvectingField->discretization_type() == "staggered" ) ;

   // Parameters
   size_t center_pos_in_matrix = 0, component = 0 ;
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
      size_t k = 0;
      center_pos_in_matrix = TF->DOF_global_number( i, j, k, component );
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
      center_pos_in_matrix = TF->DOF_global_number( i, j, k, component );
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

//----------------------------------------------------------------------
double DDS_HeatTransfer:: assemble_advection_Upwind_new( FV_DiscreteField const* AdvectingField,
	size_t advecting_level, double const& coef, size_t const& i, size_t const& j, size_t const& k, size_t advected_level) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatTransfer:: assemble_advection_Upwind_new" );
   MAC_CHECK_PRE( advecting_level < AdvectingField->storage_depth() ) ;
   MAC_ASSERT( AdvectingField->discretization_type() == "staggered" ) ;

   // Parameters
   size_t component = 0 ;
   double dxC = 0., dyC = 0., dzC = 0.;
   double AdvectedvalueC = 0., AdvectedvalueRi = 0., AdvectedvalueLe = 0.,
   	AdvectedvalueTo = 0., AdvectedvalueBo = 0.,
   	AdvectedvalueFr = 0., AdvectedvalueBe = 0,
	ur = 0., ul = 0., vt = 0., vb = 0., wf = 0., wb = 0.,
        uavg = 0., vavg = 0., wavg = 0.,
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
      size_t k = 0;
      AdvectedvalueC = TF->DOF_value( i, j, k, component, advected_level );

      // Right (X) and Left (X)
      AdvectedvalueRi = TF->DOF_value( i+1, j, k, component, advected_level );
      AdvectedvalueLe = TF->DOF_value( i-1, j, k, component, advected_level );
      ur = AdvectingField->DOF_value( i+shift.i, j, k, 0, advecting_level );
      ul = AdvectingField->DOF_value( i+shift.i-1, j, k, 0, advecting_level );
      uavg = 0.5 * (ur + ul);

      if ( uavg > 0. ) fri = uavg * AdvectedvalueC;
      else fri = uavg * AdvectedvalueRi;

      if ( uavg > 0. ) fle = uavg * AdvectedvalueLe;
      else fle = uavg * AdvectedvalueC;

      // Top (Y) and Bottom (Y)
      AdvectedvalueTo = TF->DOF_value( i, j+1, k, component, advected_level );
      AdvectedvalueBo = TF->DOF_value( i, j-1, k, component, advected_level );
      vt = AdvectingField->DOF_value( i, j+shift.j, k, 1, advecting_level );
      vb = AdvectingField->DOF_value( i, j+shift.j-1, k, 1, advecting_level );
      vavg = 0.5 * (vt + vb);

      if ( vavg > 0. ) fto = vavg * AdvectedvalueC;
      else fto = vavg * AdvectedvalueTo;

      if ( vavg > 0. ) fbo = vavg * AdvectedvalueBo;
      else fbo = vavg * AdvectedvalueC;

      flux = (fto - fbo) * dxC + (fri - fle) * dyC;

   } else {
      dzC = TF->get_cell_size( k, component, 2);
      AdvectedvalueC = TF->DOF_value( i, j, k, component, advected_level );

      // Right (X) and Left (X)
      AdvectedvalueRi = TF->DOF_value( i+1, j, k, component, advected_level );
      AdvectedvalueLe = TF->DOF_value( i-1, j, k, component, advected_level );
      ur = AdvectingField->DOF_value( i+shift.i, j, k, 0, advecting_level );
      ul = AdvectingField->DOF_value( i+shift.i-1, j, k, 0, advecting_level );
      uavg = 0.5 * (ur + ul);

      if ( uavg > 0. ) fri = uavg * AdvectedvalueC;
      else fri = uavg * AdvectedvalueRi;

      if ( uavg > 0. ) fle = uavg * AdvectedvalueLe;
      else fle = uavg * AdvectedvalueC;

      // Top (Y) and Bottom (Y)
      AdvectedvalueTo = TF->DOF_value( i, j+1, k, component, advected_level );
      AdvectedvalueBo = TF->DOF_value( i, j-1, k, component, advected_level );
      vt = AdvectingField->DOF_value( i, j+shift.j, k, 1, advecting_level );
      vb = AdvectingField->DOF_value( i, j+shift.j-1, k, 1, advecting_level );
      vavg = 0.5 * (vt + vb);

      if ( vavg > 0. ) fto = vavg * AdvectedvalueC;
      else fto = vavg * AdvectedvalueTo;

      if ( vavg > 0. ) fbo = vavg * AdvectedvalueBo;
      else fbo = vavg * AdvectedvalueC;

      // Front (Z) and Behind (Z)
      AdvectedvalueFr = TF->DOF_value( i, j, k+1, component, advected_level );
      AdvectedvalueBe = TF->DOF_value( i, j, k-1, component, advected_level );
      wf = AdvectingField->DOF_value( i, j, k+shift.k, 2, advecting_level );
      wb = AdvectingField->DOF_value( i, j, k+shift.k-1, 2, advecting_level );
      wavg = 0.5 * (wf + wb);

      if ( wavg > 0. ) ffr = wavg * AdvectedvalueC;
      else ffr = wavg * AdvectedvalueFr;

      if ( wavg > 0. ) fbe = wavg * AdvectedvalueBe;
      else fbe = wavg * AdvectedvalueC;

      flux = (fto - fbo) * dxC * dzC + (fri - fle) * dyC * dzC + (ffr - fbe) * dxC * dyC;
   }

   return (coef * flux);
}
