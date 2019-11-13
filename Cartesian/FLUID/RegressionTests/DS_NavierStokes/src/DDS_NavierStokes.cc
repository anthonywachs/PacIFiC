#include <DDS_NavierStokes.hh>
#include <FV_DomainAndFields.hh>
#include <FV_DomainBuilder.hh>
#include <FV_DiscreteField.hh>
#include <DDS_NavierStokesSystem.hh>
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


DDS_NavierStokes const* DDS_NavierStokes::PROTOTYPE
                                                 = new DDS_NavierStokes() ;


//---------------------------------------------------------------------------
DDS_NavierStokes:: DDS_NavierStokes( void )
//--------------------------------------------------------------------------
   : FV_OneStepIteration( "DDS_NavierStokes" )
   , ComputingTime("Solver")
{
   MAC_LABEL( "DDS_NavierStokes:: DDS_NavierStokes" ) ;

}




//---------------------------------------------------------------------------
DDS_NavierStokes*
DDS_NavierStokes:: create_replica( MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, dom, exp ) ) ;

   DDS_NavierStokes* result =
                        new DDS_NavierStokes( a_owner, dom, exp ) ;

   MAC_CHECK( create_replica_POST( result, a_owner, dom, exp ) ) ;
   return( result ) ;

}




//---------------------------------------------------------------------------
DDS_NavierStokes:: DDS_NavierStokes( MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : FV_OneStepIteration( a_owner, dom, exp )
   , ComputingTime("Solver")
   , UF ( dom->discrete_field( "velocity" ) )
   , PF ( dom->discrete_field( "pressure" ) )
   , GLOBAL_EQ( 0 )
   , peclet( 1. )
   , mu( 1. )
   , kai( 1. )
   , AdvectionScheme( "TVD" )
   , AdvectionTimeAccuracy( 1 )   
   , rho( 1. )
{
   MAC_LABEL( "DDS_NavierStokes:: DDS_NavierStokes" ) ;
   MAC_ASSERT( UF->discretization_type() == "staggered" ) ;
   MAC_ASSERT( PF->discretization_type() == "centered" ) ;
   //MAC_ASSERT( UF->storage_depth() == 1 ) ;

   // Call of MAC_Communicator routine to set the rank of each proces and
   // the number of processes during execution
   pelCOMM = MAC_Exec::communicator();
   my_rank = pelCOMM->rank();
   nb_procs = pelCOMM->nb_ranks();
   is_master = 0;

   is_Uperiodic[0] = false;
   is_Uperiodic[1] = false;
   is_Uperiodic[2] = false;
   is_Pperiodic[0] = false;
   is_Pperiodic[1] = false;
   is_Pperiodic[2] = false;

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
   dim = UF->primary_grid()->nb_space_dimensions() ;
   nb_comps = UF->nb_components() ;

   if ( dim == 1 )
   {
     string error_message="Space dimension should either 2 or 3";
     MAC_Error::object()->raise_bad_data_value(exp,
	"nb_space_dimensions",
	error_message );
   }

   // Create the Direction Splitting subcommunicators
   create_DDS_subcommunicators();

   // Read Density
   if ( exp->has_entry( "Density" ) )
   {
     rho = exp->double_data( "Density" ) ;
     exp->test_data( "Density", "Density>0." ) ;
   }

   // Read Viscosity
   if ( exp->has_entry( "Viscosity" ) )
   {
     mu = exp->double_data( "Viscosity" ) ;
     exp->test_data( "Viscosity", "Viscosity>0." ) ;
   }

   // Read Kai
   if ( exp->has_entry( "Kai" ) )
   {
     kai = exp->double_data( "Kai" ) ;
     exp->test_data( "Kai", "Kai>=0." ) ;
   }

   // Advection scheme
   if ( exp->has_entry( "AdvectionScheme" ) )
     AdvectionScheme = exp->string_data( "AdvectionScheme" );
   if ( AdvectionScheme != "Upwind" && AdvectionScheme != "TVD" )
   {
     string error_message="   - Upwind\n   - TVD";
     MAC_Error::object()->raise_bad_data_value( exp,
        "AdvectionScheme", error_message );
   }
   if ( AdvectionScheme == "TVD"
   	&& UF->primary_grid()->get_security_bandwidth() < 2 )
   {
     string error_message="   >= 2 with TVD scheme";
     MAC_Error::object()->raise_bad_data_value( exp,
        "security_bandwidth", error_message );
   }


   // Advection term time accuracy
   if ( exp->has_entry( "AdvectionTimeAccuracy" ) )
     AdvectionTimeAccuracy = exp->int_data( "AdvectionTimeAccuracy" );
   if ( AdvectionTimeAccuracy != 1 && AdvectionTimeAccuracy != 2 )
   {
     string error_message ="   - 1\n   - 2\n   ";
     MAC_Error::object()->raise_bad_data_value( exp,
	"AdvectionTimeAccuracy", error_message );
   }

   // Periodic boundary condition check for velocity
   U_periodic_comp = UF->primary_grid()->get_periodic_directions();
   is_Uperiodic[0] = U_periodic_comp->operator()( 0 );
   is_Uperiodic[1] = U_periodic_comp->operator()( 1 );
   if(dim >2)
      is_Uperiodic[2] = U_periodic_comp->operator()( 2 ); 

   // Periodic boundary condition check for pressure
   P_periodic_comp = PF->primary_grid()->get_periodic_directions();
   is_Pperiodic[0] = P_periodic_comp->operator()( 0 );
   is_Pperiodic[1] = P_periodic_comp->operator()( 1 );
   if(dim >2)
      is_Pperiodic[2] = P_periodic_comp->operator()( 2 ); 

   // Build the matrix system
   MAC_ModuleExplorer* se = exp->create_subexplorer( 0,"DDS_NavierStokesSystem" ) ;
   GLOBAL_EQ = DDS_NavierStokesSystem::create( this, se, UF, PF ) ;
   se->destroy() ;

   // Timing routines
   if ( my_rank == is_master )
   {
     SCT_insert_app("Matrix_Assembly&Initialization");
     SCT_insert_app("Pressure predictor");
//     SCT_insert_app("Explicit Velocity step");
//     SCT_insert_app("Velocity x update");
//     SCT_insert_app("Velocity y update");
//     SCT_insert_app("Velocity z update");
     SCT_insert_app("Velocity update");
     SCT_insert_app("Penalty Step");
     SCT_insert_app("Pressure Update");
     SCT_insert_app("Writing CSV");
     SCT_get_elapsed_time("Objects_Creation");
   }
}




//---------------------------------------------------------------------------
DDS_NavierStokes:: ~DDS_NavierStokes( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: ~DDS_NavierStokes" ) ;

   free_DDS_subcommunicators() ;

   deallocate_mpi_vectors_U();

   deallocate_mpi_vectors_P();

}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: do_one_inner_iteration( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: do_one_inner_iteration" ) ;
   MAC_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;

   start_total_timer( "DDS_NavierStokes:: do_one_inner_iteration" ) ;
   start_solving_timer() ;

   if ( my_rank == is_master ) SCT_set_start("Pressure predictor");
   NS_first_step(t_it);
   if ( my_rank == is_master ) SCT_get_elapsed_time( "Pressure predictor" );

   if ( my_rank == is_master ) SCT_set_start( "Velocity update" );
   NS_velocity_update(t_it);
   if ( my_rank == is_master ) SCT_get_elapsed_time( "Velocity update" );

   if ( my_rank == is_master ) SCT_set_start( "Penalty Step" );
   NS_pressure_update(t_it);
   if ( my_rank == is_master ) SCT_get_elapsed_time( "Penalty Step" );
   
   if ( my_rank == is_master ) SCT_set_start( "Pressure Update" );
   NS_final_step(t_it);
   if ( my_rank == is_master ) SCT_get_elapsed_time( "Pressure Update" );

   stop_solving_timer() ;
   stop_total_timer() ;

}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: do_before_time_stepping( FV_TimeIterator const* t_it,
      	std::string const& basename )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: do_before_time_stepping" ) ;

   start_total_timer( "DDS_NavierStokes:: do_before_time_stepping" ) ;

   if ( my_rank == is_master ) SCT_set_start("Matrix_Assembly&Initialization");

   FV_OneStepIteration::do_before_time_stepping( t_it, basename ) ;

   // Initialize velocity vector at the matrix level
   GLOBAL_EQ->initialize_DS_velocity();
   GLOBAL_EQ->initialize_DS_pressure();

   // Direction splitting
   // Assemble 1D tridiagonal matrices
   assemble_velocity_1D_matrices(t_it);
   //GLOBAL_EQ->call_compute_LU_decomposition();

   if ( rank_in_i[0] == 0 && nb_ranks_comm_i[0]>1){
      GLOBAL_EQ->compute_schlur_x_ref_P();
   }
   if ( rank_in_i[1] == 0 && nb_ranks_comm_i[1]>1){
      GLOBAL_EQ->compute_schlur_y_ref_P();
   }
   if(dim > 2){
      if ( rank_in_i[2] == 0 && nb_ranks_comm_i[2]>1){
        GLOBAL_EQ->compute_schlur_z_ref_P();
      }
   }

   for(size_t comp = 0;comp<nb_comps;comp++){

      if ( rank_in_i[0] == 0 && nb_ranks_comm_i[0]>1){
        GLOBAL_EQ->compute_schlur_x_ref(comp);
       }
       if ( rank_in_i[1] == 0 && nb_ranks_comm_i[1]>1){
        GLOBAL_EQ->compute_schlur_y_ref(comp);
       }
       if(dim > 2){
        if ( rank_in_i[2] == 0 && nb_ranks_comm_i[2]>1){
          GLOBAL_EQ->compute_schlur_z_ref(comp);
        }
       }
   }

   allocate_mpi_vectors_U();
   allocate_mpi_vectors_P();

   if ( my_rank == is_master ) SCT_get_elapsed_time( "Matrix_Assembly&Initialization" );

   stop_total_timer() ;

}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: do_after_time_stepping( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: do_after_time_stepping" ) ;

   // Elapsed time by sub-problems
   
   // SCT_set_start( "Writing CSV" );
   // write_pressure_field(t_it);
   // write_velocity_field(t_it);
   // SCT_get_elapsed_time( "Writing CSV" );

   if ( my_rank == is_master )
   {
     double cputime = CT_get_elapsed_time();
     cout << endl << "Full problem" << endl;
     write_elapsed_time_smhd(cout,cputime,"Computation time");
     SCT_get_summary(cout,cputime);
   }
}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: do_before_inner_iterations_stage(
	FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: do_before_inner_iterations_stage" ) ;

   start_total_timer( "DDS_NavierStokes:: do_before_inner_iterations_stage" ) ;

   FV_OneStepIteration::do_before_inner_iterations_stage( t_it ) ;

   // Perform matrix level operations before each time step
   GLOBAL_EQ->at_each_time_step( );

   stop_total_timer() ;

}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: do_after_inner_iterations_stage(
	FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: do_after_inner_iterations_stage" ) ;

   start_total_timer( "DDS_NavierStokes:: do_after_inner_iterations_stage" ) ;

   FV_OneStepIteration::do_after_inner_iterations_stage( t_it ) ;

   // Compute velocity change over the time step
   double velocity_time_change = GLOBAL_EQ->compute_directionsplitting_velocity_change()
   	/ t_it->time_step() ;
   if ( my_rank == is_master )
     cout << "velocity change = " <<
     	MAC::doubleToString( ios::scientific, 5, velocity_time_change )
	<< endl;

   get_velocity_divergence();

   //compute_CFL(t_it,0);
   double cfl = UF->compute_CFL( t_it, 0 );
   if ( my_rank == is_master )
      MAC::out() << "CFL: "<< cfl <<endl;

   stop_total_timer() ;

}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: do_additional_savings( FV_TimeIterator const* t_it,
      	int const& cycleNumber )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: do_additional_savings" ) ;

   start_total_timer( "DDS_NavierStokes:: do_additional_savings" ) ;

   stop_total_timer() ;

   GLOBAL_EQ->display_debug();
}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: allocate_mpi_vectors_U( void )
//---------------------------------------------------------------------------
{
  // Allocate MPI vectors for velocity components

  mpi_packed_data_U_x = new double* [nb_comps];
  mpi_packed_data_U_y = new double* [nb_comps];
  if(dim == 3)
    mpi_packed_data_U_z = new double* [nb_comps];

  if(rank_in_i[0] == 0){
    all_receive_data_U_x = new double** [nb_comps];
    all_send_data_U_x = new double** [nb_comps];
  }
  if(rank_in_i[1] == 0){
    all_receive_data_U_y = new double** [nb_comps];
    all_send_data_U_y = new double** [nb_comps];
  }
  if(dim == 3){
    if(rank_in_i[2] == 0){
      all_receive_data_U_z = new double** [nb_comps];
      all_send_data_U_z = new double** [nb_comps];
    }
  }
    

  size_t nb_first_send_x,nb_first_send_y,nb_first_send_z;

  size_t nb_second_send_x,nb_second_send_y,nb_second_send_z;

  size_t comp,p;

  for(comp=0;comp<nb_comps;comp++)
  {
    size_t_vector min_unknown_index(dim,comp);
    size_t_vector max_unknown_index(dim,comp);

    for (size_t l=0;l<dim;++l)
      min_unknown_index(l) =
       UF->get_min_index_unknown_handled_by_proc( comp, l ) ;
    for (size_t l=0;l<dim;++l)
      max_unknown_index(l) =
       UF->get_max_index_unknown_handled_by_proc( comp, l ) ;
    
    if(dim == 2){
      nb_first_send_x = 3*(max_unknown_index(1)-min_unknown_index(1)+1);
      nb_first_send_y = 3*(max_unknown_index(0)-min_unknown_index(0)+1);
    }
    else{
      nb_first_send_x = 3*(max_unknown_index(1)-min_unknown_index(1)+1)*(max_unknown_index(2)-min_unknown_index(2)+1);
      nb_first_send_y = 3*(max_unknown_index(0)-min_unknown_index(0)+1)*(max_unknown_index(2)-min_unknown_index(2)+1);
      nb_first_send_z = 3*(max_unknown_index(0)-min_unknown_index(0)+1)*(max_unknown_index(1)-min_unknown_index(1)+1);
      mpi_packed_data_U_z[comp] = new double[nb_first_send_z];
    }

    nb_second_send_x = (2*nb_first_send_x)/3;
    nb_second_send_y = (2*nb_first_send_y)/3;
    nb_second_send_z = (2*nb_first_send_z)/3;

    mpi_packed_data_U_x[comp] = new double[nb_first_send_x];
    mpi_packed_data_U_y[comp] = new double[nb_first_send_y];

    if(rank_in_i[0] == 0){
      
      all_receive_data_U_x[comp] = new double* [nb_ranks_comm_i[0]];
      all_send_data_U_x[comp] = new double* [nb_ranks_comm_i[0]];

      for(p = 0; p < nb_ranks_comm_i[0]; ++p) {
          all_receive_data_U_x[comp][p] = new double[nb_first_send_x];
          all_send_data_U_x[comp][p] = new double[nb_second_send_x];
      }

    }
    
    if(rank_in_i[1] == 0){

      all_receive_data_U_y[comp] = new double* [nb_ranks_comm_i[1]];
      all_send_data_U_y[comp] = new double* [nb_ranks_comm_i[1]];

      for(p = 0; p < nb_ranks_comm_i[1]; ++p) {
          all_receive_data_U_y[comp][p] = new double[nb_first_send_y];
          all_send_data_U_y[comp][p] = new double[nb_second_send_y];
      }

    }

    if(dim == 3){
      if(rank_in_i[2] == 0){
        all_receive_data_U_z[comp] = new double* [nb_ranks_comm_i[2]];
        all_send_data_U_z[comp] = new double* [nb_ranks_comm_i[2]];

        for(p = 0; p < nb_ranks_comm_i[2]; ++p) {
            all_receive_data_U_z[comp][p] = new double[nb_first_send_z];
            all_send_data_U_z[comp][p] = new double[nb_second_send_z];
        }
      }  
    }
    
    
  }

}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: allocate_mpi_vectors_P( void )
//---------------------------------------------------------------------------
{
  // Allocate MPI vectors for pressure
  if(rank_in_i[0] == 0){
    all_receive_data_P_x = new double* [nb_ranks_comm_i[0]];
    all_send_data_P_x = new double* [nb_ranks_comm_i[0]];
  }
  if(rank_in_i[1] == 0){
    all_receive_data_P_y = new double* [nb_ranks_comm_i[1]];
    all_send_data_P_y = new double* [nb_ranks_comm_i[1]];
  }
  if(dim == 3){
    if(rank_in_i[2] == 0){
      all_receive_data_P_z = new double* [nb_ranks_comm_i[2]];
      all_send_data_P_z = new double* [nb_ranks_comm_i[2]];
    }
  }
    

  size_t nb_first_send_x,nb_first_send_y,nb_first_send_z;

  size_t nb_second_send_x,nb_second_send_y,nb_second_send_z;

  size_t p,comp;

  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);

  for (size_t l=0;l<dim;++l)
    min_unknown_index(l) =
     PF->get_min_index_unknown_handled_by_proc( 0, l ) ;
  for (size_t l=0;l<dim;++l)
    max_unknown_index(l) =
     PF->get_max_index_unknown_handled_by_proc( 0, l ) ;
  
  if(dim == 2){
    nb_first_send_x = 3*(max_unknown_index(1)-min_unknown_index(1)+1);
    nb_first_send_y = 3*(max_unknown_index(0)-min_unknown_index(0)+1);
  }
  else{
    nb_first_send_x = 3*(max_unknown_index(1)-min_unknown_index(1)+1)*(max_unknown_index(2)-min_unknown_index(2)+1);
    nb_first_send_y = 3*(max_unknown_index(0)-min_unknown_index(0)+1)*(max_unknown_index(2)-min_unknown_index(2)+1);
    nb_first_send_z = 3*(max_unknown_index(0)-min_unknown_index(0)+1)*(max_unknown_index(1)-min_unknown_index(1)+1);
    mpi_packed_data_P_z = new double[nb_first_send_z];
  }

  nb_second_send_x = (2*nb_first_send_x)/3;
  nb_second_send_y = (2*nb_first_send_y)/3;
  nb_second_send_z = (2*nb_first_send_z)/3;

  mpi_packed_data_P_x = new double[nb_first_send_x];
  mpi_packed_data_P_y = new double[nb_first_send_y];

  if(rank_in_i[0] == 0){

      for(p = 0; p < nb_ranks_comm_i[0]; ++p) {
          all_receive_data_P_x[p] = new double[nb_first_send_x];
          all_send_data_P_x[p] = new double[nb_second_send_x];
      }

    }
    
    if(rank_in_i[1] == 0){
      
      for(p = 0; p < nb_ranks_comm_i[1]; ++p) {
          all_receive_data_P_y[p] = new double[nb_first_send_y];
          all_send_data_P_y[p] = new double[nb_second_send_y];
      }

    }

    if(dim == 3){
      if(rank_in_i[2] == 0){

        for(p = 0; p < nb_ranks_comm_i[2]; ++p) {
            all_receive_data_P_z[p] = new double[nb_first_send_z];
            all_send_data_P_z[p] = new double[nb_second_send_z];
        }

      }
    }
      

}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: deallocate_mpi_vectors_U( void )
//---------------------------------------------------------------------------
{
  size_t p;
  for(size_t comp=0;comp<nb_comps;comp++)
  {
    delete [] mpi_packed_data_U_x[comp];
    delete [] mpi_packed_data_U_y[comp];
    if(dim == 3)
      delete [] mpi_packed_data_U_z[comp];

    if(rank_in_i[0] == 0)
    {
      for(p = 0; p < nb_ranks_comm_i[0]; ++p) {
          delete [] all_receive_data_U_x[comp][p];
          delete [] all_send_data_U_x[comp][p];
      }
      delete [] all_receive_data_U_x[comp];
      delete [] all_send_data_U_x[comp];
    }
    
    if(rank_in_i[1] == 0){
      for(p = 0; p < nb_ranks_comm_i[1]; ++p) {
        delete [] all_receive_data_U_y[comp][p];
        delete [] all_send_data_U_y[comp][p];
      }  
      delete [] all_receive_data_U_y[comp];
      delete [] all_send_data_U_y[comp];
    }
    
    if(dim == 3){
      if(rank_in_i[2] == 0)
      {
        for(p = 0; p < nb_ranks_comm_i[2]; ++p) {
          delete [] all_receive_data_U_z[comp][p];
          delete [] all_send_data_U_z[comp][p];
        }
        delete [] all_receive_data_U_z[comp];  
        delete [] all_send_data_U_z[comp];
      }
    }
      
    

  }
  delete [] mpi_packed_data_U_x;
  delete [] mpi_packed_data_U_y;
  if(dim == 3)
    delete [] mpi_packed_data_U_z;
}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: deallocate_mpi_vectors_P( void )
//---------------------------------------------------------------------------
{
  size_t p;

  if(rank_in_i[0] == 0)
  {
    for(p = 0; p < nb_ranks_comm_i[0]; ++p) {
        delete [] all_receive_data_P_x[p];
        delete [] all_send_data_P_x[p];
    }
    delete [] all_receive_data_P_x;
    delete [] all_send_data_P_x;
  }
  
  if(rank_in_i[1] == 0)
  {
    for(p = 0; p < nb_ranks_comm_i[1]; ++p) {
        delete [] all_receive_data_P_y[p];
        delete [] all_send_data_P_y[p];
    }
    delete [] all_receive_data_P_y;
    delete [] all_send_data_P_y;
  }
  
  if(dim == 3){
    if(rank_in_i[2] == 0)
    {
      for(p = 0; p < nb_ranks_comm_i[2]; ++p) {
          delete [] all_receive_data_P_z[p];
          delete [] all_send_data_P_z[p];
      }
      delete [] all_receive_data_P_z;
      delete [] all_send_data_P_z;
    }
  }
    

  delete [] mpi_packed_data_P_x;
  delete [] mpi_packed_data_P_y;
  if(dim == 3)
    delete [] mpi_packed_data_P_z;
}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: assemble_velocity_bodyterm_rhs (
	FV_DiscreteField const* FF,
	LA_Vector* VEC_rhs )
//---------------------------------------------------------------------------
{
   MAC_LABEL(
   "DDS_NavierStokes:: assemble_velocity_bodyterm_rhs" ) ;

   if ( my_rank == is_master ) cout << "velocity body term rhs "
   	<< endl;

   // Parameters
   double dxC, dyC, dzC, xC, yC, zC, bodyterm = 0. ;
   size_t center_pos_in_matrix = 0 ;

   for (size_t comp=0;comp<nb_comps;++comp)
   {
     // Get local min and max indices
     size_t_vector min_unknown_index(dim,0);
     for (size_t l=0;l<dim;++l)
       min_unknown_index(l) =
       	FF->get_min_index_unknown_handled_by_proc( comp, l ) ;
     size_t_vector max_unknown_index(dim,0);
     for (size_t l=0;l<dim;++l)
       max_unknown_index(l) =
       	FF->get_max_index_unknown_handled_by_proc( comp, l ) ;

     // Perform assembling
     for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i)
     {
       dxC = FF->get_cell_size( i, comp, 0 ) ;
       xC = FF->get_DOF_coordinate( i, comp, 0 ) ;
       for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j)
       {

         dyC = FF->get_cell_size( j, comp, 1 ) ;
         yC = FF->get_DOF_coordinate( j, comp, 1 ) ;

         if ( dim == 2 )
	 {
	   size_t k = 0 ;
	   bodyterm = 2. * MAC::pi() * MAC::pi() * MAC::sin( MAC::pi() * xC )
	   	* MAC::sin( MAC::pi() * yC ) / peclet ;
	   center_pos_in_matrix = FF->DOF_global_number( i, j, k, comp );
	   VEC_rhs->add_to_item(
	 	center_pos_in_matrix, bodyterm * dxC * dyC );
	 }
	 else
	 {
	   for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k)
	   {
	     dzC = FF->get_cell_size( k, comp, 2 ) ;
	     zC = FF->get_DOF_coordinate( k, comp, 2 ) ;
	     bodyterm = 3. * MAC::pi() * MAC::pi() * MAC::sin( MAC::pi() * xC )
	   	* MAC::sin( MAC::pi() * yC ) * MAC::sin( MAC::pi() * zC )
		/ peclet ;
	     center_pos_in_matrix = FF->DOF_global_number( i, j, k, comp );
	     VEC_rhs->add_to_item(
	 	center_pos_in_matrix, bodyterm * dxC * dyC * dzC );
	   }
	 }
       }
     }
   }

  // Synchronize vector for parallel usage
  VEC_rhs->synchronize();

}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: error_with_analytical_solution (
	FV_DiscreteField const* FF,
	FV_DiscreteField* FF_ERROR )
//---------------------------------------------------------------------------
{
   MAC_LABEL(
   "DDS_NavierStokes:: error_with_analytical_solution" ) ;

   // Parameters
   double x, y, z, computed_field, analytical_solution, error_L2 = 0. ;

   for (size_t comp=0;comp<nb_comps;++comp)
   {
     // Get nb of local dof
     size_t_vector local_dof_number( dim, 0 );
     for (size_t l=0;l<dim;++l)
       local_dof_number(l) = FF->get_local_nb_dof( comp, l ) ;

     // Compute error
     error_L2 = 0. ;
     for (size_t i=0;i<local_dof_number(0);++i)
     {
       x = FF->get_DOF_coordinate( i, comp, 0 ) ;
       for (size_t j=0;j<local_dof_number(1);++j)
       {
         y = FF->get_DOF_coordinate( j, comp, 1 ) ;

         if ( dim == 2 )

	 {
	   size_t k = 0 ;
	   computed_field = FF->DOF_value( i, j, k, comp, 0 ) ;

	   analytical_solution = MAC::sin( MAC::pi() * x )
	   	* MAC::sin( MAC::pi() * y ) ;

	   FF_ERROR->set_DOF_value( i, j, k, comp, 0, MAC::abs( computed_field
	   	- analytical_solution ) ) ;

	   if ( FF->DOF_is_unknown_handled_by_proc( i, j, k, comp ) )
	     error_L2 += MAC::sqr( computed_field - analytical_solution )
	     	* FF->get_cell_measure( i, j, k, comp ) ;
	 }
	 else
	 {
	   for (size_t k=0;k<local_dof_number(2);++k)
	   {
             z = FF->get_DOF_coordinate( k, comp, 2 ) ;
	     computed_field = FF->DOF_value( i, j, k, comp, 0 ) ;

	     analytical_solution = MAC::sin( MAC::pi() * x )
	   	* MAC::sin( MAC::pi() * y ) * MAC::sin( MAC::pi() * z ) ;

	     FF_ERROR->set_DOF_value( i, j, k, comp, 0, MAC::abs( computed_field
	   	- analytical_solution ) ) ;

             if ( FF->DOF_is_unknown_handled_by_proc( i, j, k, comp ) )
	       error_L2 += MAC::sqr( computed_field - analytical_solution )
	     	* FF->get_cell_measure( i, j, k, comp ) ;
	   }
	 }
       }
     }

     error_L2 = pelCOMM->sum( error_L2 ) ;
     error_L2 = MAC::sqrt(error_L2);
     if ( my_rank == 0 )
       cout << "L2 Error without direction splitting, field " << FF->name() << ", component " << comp
     	<< " = " << error_L2 << endl;
   }

}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: assemble_velocity_matrix_1D (
  FV_DiscreteField const* FF,
  FV_TimeIterator const* t_it,
  double gamma,
  size_t const& comp,
  size_t const dir )
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_NavierStokes:: assemble_velocity_matrix_1D_x" ) ;

   if ( my_rank == is_master ) cout << "velocity matrix 1D in x " << endl;

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
   int m;
   size_t i;
   LA_SeqVector* Aii_main_diagonal = GLOBAL_EQ-> get_aii_main_diag(comp,dir);
   LA_SeqVector* Aii_super_diagonal = GLOBAL_EQ-> get_aii_super_diag(comp,dir);
   LA_SeqVector* Aii_sub_diagonal = GLOBAL_EQ-> get_aii_sub_diag(comp,dir);

   LA_SeqMatrix* Aie = GLOBAL_EQ-> get_aie(comp,dir);
   LA_SeqMatrix* Aei = GLOBAL_EQ-> get_aei(comp,dir);

   double Aee_diagcoef=0.;

   for (m=0,i=min_unknown_index(dir);i<=max_unknown_index(dir);++i,++m) {
      xC= FF->get_DOF_coordinate( i, comp, dir) ;
      xR= FF->get_DOF_coordinate( i+1,comp, dir) ;
      xL= FF->get_DOF_coordinate( i-1, comp, dir) ;

      dxr= xR - xC;
      dxl= xC - xL;

      right = -gamma/(dxr);
      left = -gamma/(dxl);
      center = - (right+left);

      // add unsteady term
      size_t k;
      double value;
      double unsteady_term = rho*(FF->get_cell_size(i,comp,dir))/(t_it->time_step());

      if (dim == 2) {
         k = 0;
      } else {
         k = min_unknown_index(2);
      }

      value = center;

      value = value + unsteady_term;

      bool r_bound = false;
      bool l_bound = false;
      // All the proc will have open right bound, except last proc for non periodic systems
      if ((is_Uperiodic[dir] != 1) && (rank_in_i[dir] == nb_ranks_comm_i[dir]-1)) r_bound = true;
      // All the proc will have open left bound, except first proc for non periodic systems
      if ((is_Uperiodic[dir] != 1) && (rank_in_i[dir] == 0)) l_bound = true;

      // Set Aie, Aei and Ae 
      if ((!l_bound) && (i == min_unknown_index(dir))) {
         // Periodic boundary condition at minimum unknown index
         // First proc has non zero value in Aie,Aei for first & last index
         if (rank_in_i[dir] == 0) {
            Aie->set_item(m,nb_ranks_comm_i[dir]-1,left);
            Aei->set_item(nb_ranks_comm_i[dir]-1,m,right);
         } else {
            Aie->set_item(m,rank_in_i[dir]-1,left);
            Aei->set_item(rank_in_i[dir]-1,m,right);
         }
      }

      if ((!r_bound) && (i == max_unknown_index(dir))) {
         // Periodic boundary condition at maximum unknown index
         // For last index, Aee comes from this proc as it is interface unknown wrt this proc
         Aie->set_item(m-1,rank_in_i[dir],right);
         Aee_diagcoef = value;
         Aei->set_item(rank_in_i[dir],m-1,left);
      }

      // Set Aii_sub_diagonal
      if ((rank_in_i[dir] == nb_ranks_comm_i[dir]-1) && (is_Uperiodic[dir] != 1)) {
         if (i > min_unknown_index(dir)) Aii_sub_diagonal->set_item(m-1,left);
      } else {
         if (i<max_unknown_index(dir)) {
            if (i>min_unknown_index(dir)) {
               Aii_sub_diagonal->set_item(m-1,left);
            }
         }
      }

      // Set Aii_super_diagonal
      if ((rank_in_i[dir] == nb_ranks_comm_i[dir]-1) && (is_Uperiodic[dir] != 1)) {
         if (i < max_unknown_index(dir)) Aii_super_diagonal->set_item(m,right);
      } else {
         if (i < max_unknown_index(dir)-1) {
            Aii_super_diagonal->set_item(m,right);
         }
      }

      // Set Aii_main_diagonal
      if ((rank_in_i[dir] == nb_ranks_comm_i[dir]-1) && (is_Uperiodic[dir] != 1)) {
         Aii_main_diagonal->set_item(m,value);
      } else {
         if (i<max_unknown_index(dir)) {
            Aii_main_diagonal->set_item(m,value);
         }
      }
   }

   GLOBAL_EQ->compute_Aii_ref(comp,dir);
   // Compute the product matrix for each proc

   if (nb_ranks_comm_i[dir]>1) {
      GLOBAL_EQ->compute_product_matrix(comp,dir);

      LA_SeqMatrix* product_matrix = GLOBAL_EQ-> get_Aei_Aii_Aie_product(comp,dir);
      LA_SeqMatrix* receive_matrix = product_matrix->create_copy(this,product_matrix);
      LA_SeqMatrix* Aee = GLOBAL_EQ-> get_Aee_matrix(comp,dir);
      LA_SeqMatrix* schlur_complement = GLOBAL_EQ-> get_schlur_complement(comp,dir);

      if ( rank_in_i[dir] == 0 ) {
         Aee->set_item(0,0,Aee_diagcoef);
   	 for (i=1;i<nb_ranks_comm_i[dir];++i) {
            // Create the container to receive
            size_t nbrows = product_matrix->nb_rows();
            size_t nb_received_data = pow(nbrows,2)+1;
            double * received_data = new double [nb_received_data];

            // Receive the data
            static MPI_Status status ;
            MPI_Recv( received_data, nb_received_data, MPI_DOUBLE, i, 0,
                            DDS_Comm_i[dir], &status ) ;

            // Transfer the received data to the receive matrix
            for (int k=0;k<nbrows;k++) {
               for (int j=0;j<nbrows;j++) {
                  // Assemble the global product matrix by adding contributions from all the procs
                  receive_matrix->add_to_item(k,j,received_data[k*(nbrows)+j]);
               }
            }

  	    if (is_Uperiodic[dir] == 0) {
               if (i<nb_ranks_comm_i[dir]-1) {
                  // Assemble the global Aee matrix 
                  // No periodic condition in x. So no fe contribution from last proc
                  Aee->set_item(i,i,received_data[nb_received_data-1]);
               }   
            } else {
               // Assemble the global Aee matrix
               // Periodic condition in x. So there is fe contribution from last proc
               Aee->set_item(i,i,received_data[nb_received_data-1]); 
            }
         }
      } else {
         // Create the packed data container
         size_t nbrows = product_matrix->nb_rows();
         size_t nb_send_data = pow(nbrows,2)+1;
         double * packed_data = new double [nb_send_data];

         // Fill the packed data container with Aie
         // Iterator only fetches the values present. Zeros are not fetched.

         for ( i=0 ; i<nbrows ; i++ ) {
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

      if (rank_in_i[dir] == 0){
         schlur_complement->add_Mat(Aee);
         schlur_complement->add_Mat(receive_matrix,-1.0);
      }
   }
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: assemble_pressure_matrix_1D (
  FV_DiscreteField const* FF,
  FV_TimeIterator const* t_it,
  size_t const dir )
//---------------------------------------------------------------------------
{
   MAC_LABEL("DDS_NavierStokes:: assemble_pressure_matrix_1D" ) ;

   if ( my_rank == is_master ) cout << "Pressure matrix 1D in x " << endl;

   // Parameters
   double dxr,dxl,dx,xR,xL,xC,right,left,center = 0. ;

   // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l) {
      min_unknown_index(l) = FF->get_min_index_unknown_handled_by_proc( 0, l ) ;
      max_unknown_index(l) = FF->get_max_index_unknown_handled_by_proc( 0, l ) ;
   }

   // Perform assembling
   int m;
   size_t i;
   size_t comp = 0;

   LA_SeqVector* Aii_main_diagonal = GLOBAL_EQ-> get_aii_main_diag_P(dir);
   LA_SeqVector* Aii_super_diagonal = GLOBAL_EQ-> get_aii_super_diag_P(dir);
   LA_SeqVector* Aii_sub_diagonal = GLOBAL_EQ-> get_aii_sub_diag_P(dir);

   LA_SeqMatrix* Aie = GLOBAL_EQ-> get_aie_P(dir);
   LA_SeqMatrix* Aei = GLOBAL_EQ-> get_aei_P(dir);

   double Aee_diagcoef=0.;

   for (m=0,i=min_unknown_index(dir);i<=max_unknown_index(dir);++i,++m) {
      xC= FF->get_DOF_coordinate( i, comp, dir) ;
      xR= FF->get_DOF_coordinate( i+1,comp, dir) ;
      xL= FF->get_DOF_coordinate( i-1, comp, dir) ;

      dx = FF->get_cell_size( i,comp, 0 );

      dxr= xR - xC;
      dxl= xC - xL;

      right = -1.0/(dxr);
      left = -1.0/(dxl);
      center = - (right+left);

      // add unsteady term
      size_t k;
      double value;
      double unsteady_term = 1.0*dx;

      if (dim == 2) {
         k = 0;
      } else {
         k = min_unknown_index(2);
      }

      value = center;

      value = value + unsteady_term;

      bool r_bound = false;
      bool l_bound = false;
      // All the proc will have open right bound, except last proc for non periodic systems
      if ((is_Pperiodic[dir] != 1) && (rank_in_i[dir] == nb_ranks_comm_i[dir]-1)) r_bound = true;
      // All the proc will have open left bound, except first proc for non periodic systems
      if ((is_Pperiodic[dir] != 1) && (rank_in_i[dir] == 0)) l_bound = true;

      // Set Aie, Aei and Ae 
      if ((!l_bound) && (i == min_unknown_index(dir))) {
         // Periodic boundary condition at minimum unknown index
         // First proc has non zero value in Aie,Aei for first & last index
         if (rank_in_i[dir] == 0) {
            Aie->set_item(m,nb_ranks_comm_i[dir]-1,left);
            Aei->set_item(nb_ranks_comm_i[dir]-1,m,right);
         } else {
            Aie->set_item(m,rank_in_i[dir]-1,left);
            Aei->set_item(rank_in_i[dir]-1,m,right);
         }
      }

      if ((!r_bound) && (i == max_unknown_index(dir))) {
         // Periodic boundary condition at maximum unknown index
         // For last index, Aee comes from this proc as it is interface unknown wrt this proc
         Aie->set_item(m-1,rank_in_i[dir],right);
         Aee_diagcoef = value;
         Aei->set_item(rank_in_i[dir],m-1,left);
      }

      // Set Aii_sub_diagonal
      if ((rank_in_i[dir] == nb_ranks_comm_i[dir]-1) && (is_Pperiodic[dir] != 1)) {
         if (i > min_unknown_index(dir)) Aii_sub_diagonal->set_item(m-1,left);
      } else {
         if (i<max_unknown_index(dir)) {
            if (i>min_unknown_index(dir)) {
               Aii_sub_diagonal->set_item(m-1,left);
            }
         }
      }

      // Set Aii_super_diagonal
      if ((rank_in_i[dir] == nb_ranks_comm_i[dir]-1) && (is_Pperiodic[dir] != 1)) {
         if (i < max_unknown_index(dir)) Aii_super_diagonal->set_item(m,right);
      } else {
         if (i < max_unknown_index(dir)-1) {
            Aii_super_diagonal->set_item(m,right);
         }
      }

      // Set Aii_main_diagonal
      if ((rank_in_i[dir] == nb_ranks_comm_i[dir]-1) && (is_Pperiodic[dir] != 1)) {
         Aii_main_diagonal->set_item(m,value);
      } else {
         if (i<max_unknown_index(dir)) {
            Aii_main_diagonal->set_item(m,value);
         }
      }
   }

   GLOBAL_EQ->compute_Aii_ref_P(dir);
   // Compute the product matrix for each proc

   if (nb_ranks_comm_i[dir]>1) {
      GLOBAL_EQ->compute_product_matrix_P(dir);

      LA_SeqMatrix* product_matrix = GLOBAL_EQ-> get_Aei_Aii_Aie_product_P(dir);
      LA_SeqMatrix* receive_matrix = product_matrix->create_copy(this,product_matrix);
      LA_SeqMatrix* Aee = GLOBAL_EQ-> get_Aee_matrix_P(dir);
      LA_SeqMatrix* schlur_complement = GLOBAL_EQ-> get_schlur_complement_P(dir);

      if ( rank_in_i[dir] == 0 ) {
         Aee->set_item(0,0,Aee_diagcoef);
   	 for (i=1;i<nb_ranks_comm_i[dir];++i) {
            // Create the container to receive
            size_t nbrows = product_matrix->nb_rows();
            size_t nb_received_data = pow(nbrows,2)+1;
            double * received_data = new double [nb_received_data];

            // Receive the data
            static MPI_Status status ;
            MPI_Recv( received_data, nb_received_data, MPI_DOUBLE, i, 0,
                            DDS_Comm_i[dir], &status ) ;

            // Transfer the received data to the receive matrix
            for (int k=0;k<nbrows;k++) {
               for (int j=0;j<nbrows;j++) {
                  // Assemble the global product matrix by adding contributions from all the procs
                  receive_matrix->add_to_item(k,j,received_data[k*(nbrows)+j]);
               }
            }

  	    if (is_Uperiodic[dir] == 0) {
               if (i<nb_ranks_comm_i[dir]-1) {
                  // Assemble the global Aee matrix 
                  // No periodic condition in x. So no fe contribution from last proc
                  Aee->set_item(i,i,received_data[nb_received_data-1]);
               }   
            } else {
               // Assemble the global Aee matrix
               // Periodic condition in x. So there is fe contribution from last proc
               Aee->set_item(i,i,received_data[nb_received_data-1]); 
            }
         }
      } else {
         // Create the packed data container
         size_t nbrows = product_matrix->nb_rows();
         size_t nb_send_data = pow(nbrows,2)+1;
         double * packed_data = new double [nb_send_data];

         // Fill the packed data container with Aie
         // Iterator only fetches the values present. Zeros are not fetched.

         for ( i=0 ; i<nbrows ; i++ ) {
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

      if (rank_in_i[dir] == 0){
         schlur_complement->add_Mat(Aee);
         schlur_complement->add_Mat(receive_matrix,-1.0);
      }
   }
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: assemble_velocity_1D_matrices ( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: assemble_velocity_1D_matrices" ) ;

   // Assemble the matrices for each component

   for (size_t dir=0;dir<dim;dir++) {
      assemble_pressure_matrix_1D (PF,t_it,dir);
   }

   for (size_t comp=0;comp<nb_comps;comp++) {
      double gamma = mu/2.0;
      for (size_t dir=0;dir<dim;dir++) {
         assemble_velocity_matrix_1D (UF,t_it,gamma,comp,dir);
      }
   }
}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: NS_first_step ( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_NavierStokes:: NS_first_step" ) ;

   size_t i, j, k;

  double value=0.;

  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);
  // First Equation

    // Get local min and max indices
    for (size_t l=0;l<dim;++l)
      min_unknown_index(l) =
       PF->get_min_index_unknown_on_proc( 0, l ) ;
    for (size_t l=0;l<dim;++l)
      max_unknown_index(l) =
       PF->get_max_index_unknown_on_proc( 0, l ) ;
    for (i=min_unknown_index(0);i<=max_unknown_index(0);++i)
    {
      for (j=min_unknown_index(1);j<=max_unknown_index(1);++j)
      {
        if(dim ==2 )
        {
          k=0;
          // Set P*_n as sum of P_(n-1/2)+phi_(n-1/2)
          value = PF->DOF_value( i, j, k, 0, 0 ) + PF->DOF_value( i, j, k, 0, 1 );
          PF->set_DOF_value( i, j, k, 0, 1, value);
        }
        else
        {
          for (k=min_unknown_index(2);k<=max_unknown_index(2);++k)
          {
            // Set P*_n as sum of P_(n-1/2)+phi_(n-1/2)
            value = PF->DOF_value( i, j, k, 0, 0 ) + PF->DOF_value( i, j, k, 0, 1 );
            PF->set_DOF_value( i, j, k, 0, 1, value);
          }
        }
      }
    }
    PF->set_neumann_DOF_values();
    // What to do for periodic conditions in this function??????????????????
}




//---------------------------------------------------------------------------
double
DDS_NavierStokes:: assemble_local_rhs_x ( size_t const& j, size_t const& k, double gamma, FV_TimeIterator const* t_it, size_t const& comp)
//---------------------------------------------------------------------------
{
  // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     min_unknown_index(l) =
      UF->get_min_index_unknown_handled_by_proc( comp, l ) ;
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     max_unknown_index(l) =
      UF->get_max_index_unknown_handled_by_proc( comp, l ) ;

 size_t i,pos;
 int m;

  // Compute VEC_rhs_x = rhs in x
  double dxC,xhr,xhl,xright,xleft;
  double fe=0.;

  // Vector for fi
  LA_SeqVector* local_rhs_x = GLOBAL_EQ->get_local_temp_x(comp) ;

  for (i=min_unknown_index(0);i<=max_unknown_index(0);++i)
  {
    double xvalue=0.;
    // Get contribution of un
    xhl= UF->get_DOF_coordinate( i, comp, 0 ) - UF->get_DOF_coordinate( i-1, comp, 0 ) ;
    xhr= UF->get_DOF_coordinate( i+1,comp, 0 ) - UF->get_DOF_coordinate( i, comp, 0 ) ;
    xright = UF->DOF_value( i+1, j, k, comp, 1 ) - UF->DOF_value( i, j, k, comp, 1 ) ;
    xleft = UF->DOF_value( i, j, k, comp, 1 ) - UF->DOF_value( i-1, j, k, comp, 1 ) ;

    pos = i - min_unknown_index(0);

    /** Assemble rhs in x for Neumann/Dirichlet BC Conditions */
    if(UF->DOF_in_domain(i-1, j, k, comp) && UF->DOF_in_domain(i+1, j, k, comp))
      xvalue = xright/xhr - xleft/xhl;
    else if(UF->DOF_in_domain(i-1, j, k, comp))
      xvalue = - xleft/xhl;
    else
      xvalue = xright/xhr ;

    dxC = UF->get_cell_size( i, comp, 0 ) ;

    if(is_Uperiodic[0] == 0){
      if(rank_in_i[0] == nb_ranks_comm_i[0]-1){
         local_rhs_x->set_item( pos, (UF->DOF_value( i, j, k, comp, 0 )*rho*dxC)/(t_it -> time_step()) - gamma*xvalue );
      }
      else{
         if(i == max_unknown_index(0))
           fe = (UF->DOF_value( i, j, k, comp, 0 )*rho*dxC)/(t_it -> time_step()) - gamma*xvalue;
         else
           local_rhs_x->set_item( pos, (UF->DOF_value( i, j, k, comp, 0 )*rho*dxC)/(t_it -> time_step()) - gamma*xvalue );
      }
    }
    else{
      if(nb_ranks_comm_i[0] > 1){
         if(i == max_unknown_index(0))
           fe = (UF->DOF_value( i, j, k, comp, 0 )*rho*dxC)/(t_it -> time_step()) - gamma*xvalue;
         else
           local_rhs_x->set_item( pos, (UF->DOF_value( i, j, k, comp, 0 )*rho*dxC)/(t_it -> time_step()) - gamma*xvalue );
      }
      else
        local_rhs_x->set_item( pos, (UF->DOF_value( i, j, k, comp, 0 )*rho*dxC)/(t_it -> time_step()) - gamma*xvalue ); 
    }
    
  }

  m = int(min_unknown_index(0)) - 1;
  if ( UF->DOF_in_domain( m, j, k, comp ) )
    if ( UF->DOF_has_imposed_Dirichlet_value( m, j, k, comp ) )
    {
      double ai = 1/(UF->get_DOF_coordinate( m+1,comp,0) - UF->get_DOF_coordinate( m,comp,0));
      double dirichlet_value = UF->DOF_value( m, j, k, comp, 1 ) ;
      local_rhs_x->add_to_item( 0, + gamma * ai * dirichlet_value );
    }

  m = int(max_unknown_index(0)) + 1;
  if ( UF->DOF_in_domain( m, j, k, comp ) )
    if ( UF->DOF_has_imposed_Dirichlet_value( m, j, k, comp ) )
    {
      double ai = 1/(UF->get_DOF_coordinate( m,comp,0) - UF->get_DOF_coordinate( m-1,comp,0));
      double dirichlet_value = UF->DOF_value( m, j, k, comp, 1 ) ;
      local_rhs_x->add_to_item( local_rhs_x->nb_rows()-1 , + gamma * ai * dirichlet_value );
    } 
    return fe;
}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: solve_x_for_secondorder ( size_t const& j, size_t const& k, double gamma, FV_TimeIterator const* t_it,double * packed_data, size_t const& comp)
//---------------------------------------------------------------------------
{
  // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     min_unknown_index(l) =
      UF->get_min_index_unknown_handled_by_proc( comp, l ) ;
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     max_unknown_index(l) =
      UF->get_max_index_unknown_handled_by_proc( comp, l ) ;

   double fe = assemble_local_rhs_x(j,k,gamma,t_it,comp);

   // create a replica of local rhs vector in local solution vector
   LA_SeqVector* local_rhs_x = GLOBAL_EQ->get_local_temp_x(comp) ;

   LA_SeqVector* local_solution_x = GLOBAL_EQ->get_local_solution_temp_x(comp) ;
   for(size_t i=0;i<local_rhs_x->nb_rows();i++){
    local_solution_x->set_item(i,local_rhs_x->item(i));
   }

   // Solve for xi locally and put it in local solution vector
   GLOBAL_EQ->DS_NavierStokes_x_local_unknown_solver(local_solution_x,comp);

   // if(comp == 0){
   //    MAC::out()<<"J:"<<j<<endl;
   //    local_solution_x->print_items(MAC::out(),0); 
   // }

   LA_SeqMatrix* Aei_x = GLOBAL_EQ-> get_aei(comp,0);
   LA_SeqVector* Vec_temp_x = GLOBAL_EQ-> get_temp_x(comp);

   for(int i=0;i<Vec_temp_x->nb_rows();i++){
          Vec_temp_x->set_item(i,0);
    }
   // Calculate Aei*xi in each proc locally
   Aei_x->multiply_vec_then_add(local_solution_x,Vec_temp_x);

   // Pack the data
   // Pack the data
   size_t vec_pos;
   if(dim == 2){
      vec_pos=j-min_unknown_index(1);
   }
   else{
      vec_pos=(j-min_unknown_index(1))+(max_unknown_index(1)-min_unknown_index(1)+1)*(k-min_unknown_index(2));
   }

   if(rank_in_i[0] == 0){
      if(is_Uperiodic[0])
          packed_data[3*vec_pos+0]=Vec_temp_x->item(nb_ranks_comm_i[0]-1);
      else
          packed_data[3*vec_pos+0]=0;
      packed_data[3*vec_pos+1]=Vec_temp_x->item(rank_in_i[0]);
   }
   else if(rank_in_i[0] == nb_ranks_comm_i[0]-1){
      packed_data[3*vec_pos+0]=Vec_temp_x->item(rank_in_i[0]-1);
      if(is_Uperiodic[0])
          packed_data[3*vec_pos+1]=Vec_temp_x->item(rank_in_i[0]);
      else
          packed_data[3*vec_pos+1]=0;
   }
   else{
      packed_data[3*vec_pos+0]=Vec_temp_x->item(rank_in_i[0]-1);
      packed_data[3*vec_pos+1]=Vec_temp_x->item(rank_in_i[0]);
   }
   packed_data[3*vec_pos+2] = fe; // Send the fe values and 0 for last proc
   

}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: solve_interface_unknowns_x ( double * packed_data, size_t nb_received_data, double gamma,  FV_TimeIterator const* t_it, size_t const& comp)
//---------------------------------------------------------------------------
{
   size_t i,j,p;
   size_t k =0;

   // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     min_unknown_index(l) =
      UF->get_min_index_unknown_handled_by_proc( comp, l ) ;
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     max_unknown_index(l) =
      UF->get_max_index_unknown_handled_by_proc( comp, l ) ;

   // Vector for fe
   LA_SeqVector* interface_rhs_x = GLOBAL_EQ->get_interface_temp_x(comp) ;
   // Vector for fi
   LA_SeqVector* local_rhs_x = GLOBAL_EQ->get_local_temp_x(comp) ;

   LA_SeqVector* Vec_temp_x = GLOBAL_EQ-> get_temp_x(comp);

   LA_SeqMatrix* Aie_x = GLOBAL_EQ-> get_aie(comp,0);

   size_t nb_send_data = (nb_received_data*2)/3;

  // Send and receive the data first pass
   //if ( my_rank == is_master )
   if ( rank_in_i[0] == 0 )
   {
      for (i=1;i<nb_ranks_comm_i[0];++i)
      //for (i=1;i<nb_procs;++i)
      {
        // Receive the data
        //pelCOMM->receive( i, all_received_data[i], nb_received_data );
        static MPI_Status status;
        MPI_Recv( all_receive_data_U_x[comp][i], nb_received_data, MPI_DOUBLE, i, 0,
                          DDS_Comm_i[0], &status ) ;
      }

      // Solve system of interface unknowns for each y
      if(dim == 2)
      {
        for(j=min_unknown_index(1);j<=max_unknown_index(1);j++){

          size_t nb_interface_unknowns = Vec_temp_x->nb_rows();
          for(i=0;i<nb_interface_unknowns;i++){
            Vec_temp_x->set_item(i,0);
            interface_rhs_x->set_item(i,0);
          }
          p = j-min_unknown_index(1);
          if(is_Uperiodic[0])
              Vec_temp_x->set_item(nb_ranks_comm_i[0]-1,packed_data[3*p]);
          Vec_temp_x->set_item(0,packed_data[3*p+1]);
          interface_rhs_x->set_item(0,packed_data[3*p+2]);

          // Vec_temp might contain previous values

          for(i=1;i<nb_ranks_comm_i[0];i++){

            if(i!=nb_ranks_comm_i[0]-1){
              Vec_temp_x->add_to_item(i-1,all_receive_data_U_x[comp][i][3*p]);
              Vec_temp_x->add_to_item(i,all_receive_data_U_x[comp][i][3*p+1]);
              interface_rhs_x->set_item(i,all_receive_data_U_x[comp][i][3*p+2]);  // Assemble the interface rhs fe
            }
            else{
              if(is_Uperiodic[0] ==0){
                  Vec_temp_x->add_to_item(i-1,all_receive_data_U_x[comp][i][3*p]);
              }
              else{
                  Vec_temp_x->add_to_item(i-1,all_receive_data_U_x[comp][i][3*p]);
                  // If periodic in x, last proc has an interface unknown
                  Vec_temp_x->add_to_item(i,all_receive_data_U_x[comp][i][3*p+1]);
                  interface_rhs_x->set_item(i,all_receive_data_U_x[comp][i][3*p+2]);
              }
            }
          }

          for(i=0;i<nb_interface_unknowns;i++){
          interface_rhs_x->set_item(i,interface_rhs_x->item(i)-Vec_temp_x->item(i)); // Get fe - Aei*xi to solve for ue
          }

          // Solve for ue (interface unknowns) in the master proc
          GLOBAL_EQ->DS_NavierStokes_x_interface_unknown_solver(interface_rhs_x,comp);

          // Pack the interface_rhs_x into the appropriate send_data
          for (i=1;i<nb_ranks_comm_i[0];++i)
          {
            if(i!=nb_ranks_comm_i[0]-1){
              all_send_data_U_x[comp][i][2*p+0]=interface_rhs_x->item(i-1);
              all_send_data_U_x[comp][i][2*p+1]=interface_rhs_x->item(i);
            }
            else{
              all_send_data_U_x[comp][i][2*p+0]=interface_rhs_x->item(i-1);
              if(is_Uperiodic[0])
                all_send_data_U_x[comp][i][2*p+1]=interface_rhs_x->item(i);  
              else
                all_send_data_U_x[comp][i][2*p+1]=0;
            }

          }

          // Need to have the original rhs function assembled for corrosponding j,k pair
          double fe = assemble_local_rhs_x(j,k,gamma,t_it,comp);

          // Setup RHS = fi - Aie*xe for solving ui
          Aie_x->multiply_vec_then_add(interface_rhs_x,local_rhs_x,-1.0,1.0);


          // Solve ui and transfer solution into distributed vector
          GLOBAL_EQ->DS_NavierStokes_x_solver(j,k,min_unknown_index(0),local_rhs_x,interface_rhs_x,comp);

          // if(comp == 0){
          //   MAC::out()<<"J is"<<j<<endl;
          //   local_rhs_x->print_items(MAC::out(),0);  
          // }

          // if(comp == 0)
          //       MAC::out()<<"x: "<<rank_in_i[0]<<" & y: "<<rank_in_i[1]<<" ,min: "<<local_rhs_x->item(0)<<endl;

        }
      }
      else
      {
        for(k=min_unknown_index(2);k<=max_unknown_index(2);k++)
        {
          for(j=min_unknown_index(1);j<=max_unknown_index(1);j++)
          {
            size_t nb_interface_unknowns = Vec_temp_x->nb_rows();
            for(i=0;i<nb_interface_unknowns;i++){
              Vec_temp_x->set_item(i,0);
              interface_rhs_x->set_item(i,0);
            }
            p = (j-min_unknown_index(1))+(max_unknown_index(1)-min_unknown_index(1)+1)*(k-min_unknown_index(2));
            if(is_Uperiodic[0])
                Vec_temp_x->set_item(nb_ranks_comm_i[0]-1,packed_data[3*p]);
            Vec_temp_x->set_item(0,packed_data[3*p+1]);
            interface_rhs_x->set_item(0,packed_data[3*p+2]);

            // Vec_temp might contain previous values


            for(i=1;i<nb_ranks_comm_i[0];i++){

              if(i!=nb_ranks_comm_i[0]-1){
                Vec_temp_x->add_to_item(i-1,all_receive_data_U_x[comp][i][3*p]);
                Vec_temp_x->add_to_item(i,all_receive_data_U_x[comp][i][3*p+1]);
                interface_rhs_x->set_item(i,all_receive_data_U_x[comp][i][3*p+2]);  // Assemble the interface rhs fe
              }
              else{
                if(is_Uperiodic[0] ==0){
                    Vec_temp_x->add_to_item(i-1,all_receive_data_U_x[comp][i][3*p]);
                }
                else{
                    Vec_temp_x->add_to_item(i-1,all_receive_data_U_x[comp][i][3*p]);
                    // If periodic in x, last proc has an interface unknown
                    Vec_temp_x->add_to_item(i,all_receive_data_U_x[comp][i][3*p+1]);
                    interface_rhs_x->set_item(i,all_receive_data_U_x[comp][i][3*p+2]);
                }
              }
            }

            for(i=0;i<nb_interface_unknowns;i++){
            interface_rhs_x->set_item(i,interface_rhs_x->item(i)-Vec_temp_x->item(i)); // Get fe - Aei*xi to solve for ue
            }


            // Solve for ue (interface unknowns) in the master proc
            GLOBAL_EQ->DS_NavierStokes_x_interface_unknown_solver(interface_rhs_x,comp);

            // Pack the interface_rhs_x into the appropriate send_data
            for (i=1;i<nb_ranks_comm_i[0];++i)
            {
              if(i!=nb_ranks_comm_i[0]-1){
                all_send_data_U_x[comp][i][2*p+0]=interface_rhs_x->item(i-1);
                all_send_data_U_x[comp][i][2*p+1]=interface_rhs_x->item(i);
              }
              else{
                all_send_data_U_x[comp][i][2*p+0]=interface_rhs_x->item(i-1);
                if(is_Uperiodic[0])
                  all_send_data_U_x[comp][i][2*p+1]=interface_rhs_x->item(i);  
                else
                  all_send_data_U_x[comp][i][2*p+1]=0;
              }

            }

          // Need to have the original rhs function assembled for corrosponding j,k pair
          double fe = assemble_local_rhs_x(j,k,gamma,t_it,comp);

          // Setup RHS = fi - Aie*xe for solving ui
          Aie_x->multiply_vec_then_add(interface_rhs_x,local_rhs_x,-1.0,1.0);

          // Solve ui and transfer solution into distributed vector
          GLOBAL_EQ->DS_NavierStokes_x_solver(j,k,min_unknown_index(0),local_rhs_x,interface_rhs_x,comp);
          }
        }
      }


    }
   else
   {
      // Send the packed data to master
      //pelCOMM->send( is_master, packed_data, nb_received_data );
      MPI_Send( packed_data, nb_received_data, MPI_DOUBLE, 0, 0, DDS_Comm_i[0] ) ;
   }

   // Send the data from master
   if ( rank_in_i[0] == 0 )
   {
     for (i=1;i<nb_ranks_comm_i[0];++i)
      {
        //pelCOMM->send( i, all_send_data[i], nb_send_data );
        MPI_Send( all_send_data_U_x[comp][i], nb_send_data, MPI_DOUBLE, i, 0, DDS_Comm_i[0] ) ;
      }
   }
   else
   {
      // Receive the data

      //pelCOMM->receive( is_master, received_data, nb_received_data );
      static MPI_Status status ;
      MPI_Recv( packed_data, nb_received_data, MPI_DOUBLE, 0, 0,
                          DDS_Comm_i[0], &status ) ;

     // Solve the system of equations in each proc

     if(dim == 2)
     {
        for(j = min_unknown_index(1);j<=max_unknown_index(1);j++)
        {
          k=0;
          p = j-min_unknown_index(1);
          if(rank_in_i[0] != nb_ranks_comm_i[0]-1){
            interface_rhs_x->set_item(rank_in_i[0]-1,packed_data[2*p]);
            interface_rhs_x->set_item(rank_in_i[0],packed_data[2*p+1]);
          }
          else{
            if(is_Uperiodic[0] ==0){
              interface_rhs_x->set_item(rank_in_i[0]-1,packed_data[2*p]);
            }
            else{
              interface_rhs_x->set_item(rank_in_i[0]-1,packed_data[2*p]);
              interface_rhs_x->set_item(rank_in_i[0],packed_data[2*p+1]);
            }
          }
          // Need to have the original rhs function assembled for corrosponding j,k pair
          double fe = assemble_local_rhs_x(j,k,gamma,t_it,comp);

          // Setup RHS = fi - Aie*xe for solving ui
          Aie_x->multiply_vec_then_add(interface_rhs_x,local_rhs_x,-1.0,1.0);

          // Solve ui and transfer solution into distributed vector
          GLOBAL_EQ->DS_NavierStokes_x_solver(j,k,min_unknown_index(0),local_rhs_x,interface_rhs_x,comp);

          // if(comp == 0)
          //       MAC::out()<<"x: "<<rank_in_i[0]<<" & y: "<<rank_in_i[1]<<" ,max: "<<interface_rhs_x->item(1)<<endl;

        }
     }
     else
     {
        for(k = min_unknown_index(2);k<=max_unknown_index(2);k++)
        {
          for(j = min_unknown_index(1);j<=max_unknown_index(1);j++)
          {
            p = (j-min_unknown_index(1))+(max_unknown_index(1)-min_unknown_index(1)+1)*(k-min_unknown_index(2));
            if(rank_in_i[0] != nb_ranks_comm_i[0]-1){
              interface_rhs_x->set_item(rank_in_i[0]-1,packed_data[2*p]);
              interface_rhs_x->set_item(rank_in_i[0],packed_data[2*p+1]);
            }
            else{
              if(is_Uperiodic[0] ==0){
                interface_rhs_x->set_item(rank_in_i[0]-1,packed_data[2*p]);
              }
              else{
                interface_rhs_x->set_item(rank_in_i[0]-1,packed_data[2*p]);
                interface_rhs_x->set_item(rank_in_i[0],packed_data[2*p+1]);
              }
            }
            // Need to have the original rhs function assembled for corrosponding j,k pair
            double fe = assemble_local_rhs_x(j,k,gamma,t_it,comp);

            // Setup RHS = fi - Aie*xe for solving ui
            Aie_x->multiply_vec_then_add(interface_rhs_x,local_rhs_x,-1.0,1.0);

            // Solve ui and transfer solution into distributed vector
            GLOBAL_EQ->DS_NavierStokes_x_solver(j,k,min_unknown_index(0),local_rhs_x,interface_rhs_x,comp);
          }

        }
     }

   }
}




//---------------------------------------------------------------------------
double
DDS_NavierStokes:: assemble_local_rhs_y ( size_t const& i, size_t const& k, double gamma, FV_TimeIterator const* t_it, size_t const& comp )
//---------------------------------------------------------------------------
{
  // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     min_unknown_index(l) =
      UF->get_min_index_unknown_handled_by_proc( comp, l ) ;
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     max_unknown_index(l) =
      UF->get_max_index_unknown_handled_by_proc( comp, l ) ;

size_t j,pos;
int m;

// Compute VEC_rhs_x = rhs in x
double dyC,yhr,yhl,yright,yleft;
double fe=0.;

 // Vector for fi
LA_SeqVector* local_rhs_y = GLOBAL_EQ->get_local_temp_y(comp) ;

 // Compute VEC_rhs_y = rhs in y
 for (j=min_unknown_index(1);j<=max_unknown_index(1);++j)
 {
   double yvalue=0.;

   yhr= UF->get_DOF_coordinate( j+1,comp, 1 ) - UF->get_DOF_coordinate( j, comp, 1 ) ;
   yhl= UF->get_DOF_coordinate( j, comp, 1 ) - UF->get_DOF_coordinate( j-1, comp, 1 ) ;
   yright = UF->DOF_value( i, j+1, k, comp, 1 ) - UF->DOF_value( i, j, k, comp, 1 ) ;
   yleft = UF->DOF_value( i, j, k, comp, 1 ) - UF->DOF_value( i, j-1, k, comp, 1 ) ;

   pos = j - min_unknown_index(1);

   /** Assemble rhs in y for Neumann/Dirichlet BC Conditions */
   if(UF->DOF_in_domain(i, j-1, k, comp) && UF->DOF_in_domain(i, j+1, k, comp))
     yvalue = yright/yhr - yleft/yhl;
   else if(UF->DOF_in_domain(i, j-1, k, comp))
     yvalue = - yleft/yhl;
   else
     yvalue = yright/yhr;

   dyC = UF->get_cell_size( j, comp, 1 ) ;

   if(is_Uperiodic[1] == 0){
      if(rank_in_i[1] == nb_ranks_comm_i[1]-1){
         local_rhs_y->set_item( pos, (UF->DOF_value( i, j, k, comp, 0 )*rho*dyC)/(t_it -> time_step()) - gamma*yvalue );
       }
       else{
         if(j == max_unknown_index(1))
           fe = (UF->DOF_value( i, j, k, comp, 0 )*rho*dyC)/(t_it -> time_step()) - gamma*yvalue ;
         else
           local_rhs_y->set_item( pos, (UF->DOF_value( i, j, k, comp, 0 )*rho*dyC)/(t_it -> time_step()) - gamma*yvalue );
       }
   }
   else{
      if(nb_ranks_comm_i[1] > 1){
         if(j == max_unknown_index(1))
           fe = (UF->DOF_value( i, j, k, comp, 0 )*rho*dyC)/(t_it -> time_step()) - gamma*yvalue ;
         else
           local_rhs_y->set_item( pos, (UF->DOF_value( i, j, k, comp, 0 )*rho*dyC)/(t_it -> time_step()) - gamma*yvalue );
      }
      else
        local_rhs_y->set_item( pos, (UF->DOF_value( i, j, k, comp, 0 )*rho*dyC)/(t_it -> time_step()) - gamma*yvalue );
   }
   
 }

 m = int(min_unknown_index(1)) - 1;
 if ( UF->DOF_in_domain( i, m, k, comp ) )
   if ( UF->DOF_has_imposed_Dirichlet_value( i, m, k, comp ) )
   {
     double ai = 1/(UF->get_DOF_coordinate( m+1,comp,1) - UF->get_DOF_coordinate( m,comp,1));
     double dirichlet_value = UF->DOF_value( i, m, k, comp, 1 ) ;
     local_rhs_y->add_to_item( 0, + gamma * ai * dirichlet_value );
   }

 m = int(max_unknown_index(1)) + 1;
 if ( UF->DOF_in_domain( i, m, k, comp ) )
   if ( UF->DOF_has_imposed_Dirichlet_value( i, m, k, comp ) )
   {
     double ai = 1/(UF->get_DOF_coordinate( m,comp,1) - UF->get_DOF_coordinate( m-1,comp,1));
     double dirichlet_value = UF->DOF_value( i, m, k, comp, 1 ) ;
     local_rhs_y->add_to_item( local_rhs_y->nb_rows()-1, + gamma * ai * dirichlet_value );
   }

  return fe;
}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: solve_y_for_secondorder ( size_t const& i, size_t const& k, double gamma, FV_TimeIterator const* t_it,double * packed_data, size_t const& comp )
//---------------------------------------------------------------------------
{
  // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     min_unknown_index(l) =
      UF->get_min_index_unknown_handled_by_proc( comp, l ) ;
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     max_unknown_index(l) =
      UF->get_max_index_unknown_handled_by_proc( comp, l ) ;

   double fe = assemble_local_rhs_y(i,k,gamma,t_it,comp);

   size_t j;
   // create a replica of local rhs vector in local solution vector
   LA_SeqVector* local_rhs_y = GLOBAL_EQ->get_local_temp_y(comp) ;

   LA_SeqVector* local_solution_y = GLOBAL_EQ->get_local_solution_temp_y(comp) ;
   for(j=0;j<local_rhs_y->nb_rows();j++){
    local_solution_y->set_item(j,local_rhs_y->item(j));
   }

   // Solve for xi locally and put it in local solution vector
   GLOBAL_EQ->DS_NavierStokes_y_local_unknown_solver(local_solution_y,comp);

   LA_SeqMatrix* Aei = GLOBAL_EQ-> get_aei(comp,1);
   LA_SeqVector* Vec_temp = GLOBAL_EQ-> get_temp_y(comp);

   for(j=0;j<Vec_temp->nb_rows();j++){
          Vec_temp->set_item(j,0);
    }
   // Calculate Aei*xi in each proc locally
   Aei->multiply_vec_then_add(local_solution_y,Vec_temp);

   size_t vec_pos;
   if(dim == 2){
      vec_pos=i-min_unknown_index(0);
   }
   else{
      vec_pos=(k-min_unknown_index(2))+(max_unknown_index(2)-min_unknown_index(2)+1)*(i-min_unknown_index(0));
   }

   if(rank_in_i[1] == 0){
    // Check if bc is periodic in x
     // If it is, we need to pack two elements apart from fe

      if(is_Uperiodic[1])
          packed_data[3*vec_pos+0]=Vec_temp->item(nb_ranks_comm_i[1]-1);
      else
          packed_data[3*vec_pos+0]=0;
      packed_data[3*vec_pos+1]=Vec_temp->item(rank_in_i[1]);
    }
    else if(rank_in_i[1] == nb_ranks_comm_i[1]-1){
      packed_data[3*vec_pos+0]=Vec_temp->item(rank_in_i[1]-1);
      // Check if bc is periodic in x
      // If it is, we need to pack two elements apart from fe
      if(is_Uperiodic[1])
          packed_data[3*vec_pos+1]=Vec_temp->item(rank_in_i[1]);
      else
          packed_data[3*vec_pos+1]=0;
    }
    else{
      packed_data[3*vec_pos+0]=Vec_temp->item(rank_in_i[1]-1);
      packed_data[3*vec_pos+1]=Vec_temp->item(rank_in_i[1]);
    }
    packed_data[3*vec_pos+2] = fe; // Send the fe values and 0 for last proc
}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: solve_interface_unknowns_y ( double * packed_data, size_t nb_received_data, double gamma,  FV_TimeIterator const* t_it, size_t const& comp  )
//---------------------------------------------------------------------------
{
   size_t i,j,p;
   size_t k =0;

   // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     min_unknown_index(l) =
      UF->get_min_index_unknown_handled_by_proc( comp, l ) ;
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     max_unknown_index(l) =
      UF->get_max_index_unknown_handled_by_proc( comp, l ) ;

   // Vector for fe
   LA_SeqVector* interface_rhs_y = GLOBAL_EQ->get_interface_temp_y(comp) ;
   // Vector for fi
   LA_SeqVector* local_rhs_y = GLOBAL_EQ->get_local_temp_y(comp);

   LA_SeqVector* Vec_temp = GLOBAL_EQ-> get_temp_y(comp);

   LA_SeqMatrix* Aie = GLOBAL_EQ-> get_aie(comp,1);

   size_t nb_send_data = (nb_received_data*2)/3;

  // Send and receive the data first pass
   if ( rank_in_i[1] == 0 )
   {
      for (i=1;i<nb_ranks_comm_i[1];++i)
      {
        // Receive the data
        //pelCOMM->receive( i, all_received_data[i], nb_received_data );
        static MPI_Status status ;
        MPI_Recv( all_receive_data_U_y[comp][i], nb_received_data, MPI_DOUBLE, i, 0,
                          DDS_Comm_i[1], &status ) ;
      }

      // Solve system of interface unknowns for each x
      if(dim == 2)
      {
        for(i=min_unknown_index(0);i<=max_unknown_index(0);i++)
        {
          size_t nb_interface_unknowns = Vec_temp->nb_rows();
          for(j=0;j<nb_interface_unknowns;j++){
            Vec_temp->set_item(j,0);
            interface_rhs_y->set_item(j,0);
          }
          p = i-min_unknown_index(0);
          if(is_Uperiodic[1])
              Vec_temp->set_item(nb_ranks_comm_i[1]-1,packed_data[3*p]);
          Vec_temp->set_item(0,packed_data[3*p+1]);
          interface_rhs_y->set_item(0,packed_data[3*p+2]);

          // Vec_temp might contain previous values

          for(j=1;j<nb_ranks_comm_i[1];j++){

            if(j!=nb_ranks_comm_i[1]-1){
              Vec_temp->add_to_item(j-1,all_receive_data_U_y[comp][j][3*p]);
              Vec_temp->add_to_item(j,all_receive_data_U_y[comp][j][3*p+1]);
              interface_rhs_y->set_item(j,all_receive_data_U_y[comp][j][3*p+2]);  // Assemble the interface rhs fe
            }
            else{
              if(is_Uperiodic[1] ==0){
                  Vec_temp->add_to_item(j-1,all_receive_data_U_y[comp][j][3*p]);
              }
              else{
                  // If periodic in y, last proc has an interface unknown
                  Vec_temp->add_to_item(j-1,all_receive_data_U_y[comp][j][3*p]);
                  Vec_temp->add_to_item(j,all_receive_data_U_y[comp][j][3*p+1]);
                  interface_rhs_y->set_item(j,all_receive_data_U_y[comp][j][3*p+2]);
              }
            }
          }

          for(j=0;j<nb_interface_unknowns;j++){
          interface_rhs_y->set_item(j,interface_rhs_y->item(j)-Vec_temp->item(j)); // Get fe - Aei*xi to solve for ue
          }

          // Solve for ue (interface unknowns) in the master proc
          GLOBAL_EQ->DS_NavierStokes_y_interface_unknown_solver(interface_rhs_y,comp);

          //if(p==0 && rank_in_i[0] == 1)
          //  interface_rhs_y->print_items(MAC::out(),0);

          // Pack the interface_rhs_x into the appropriate send_data
          for (j=1;j<nb_ranks_comm_i[1];++j)
          {
            if(j!=nb_ranks_comm_i[1]-1){
              all_send_data_U_y[comp][j][2*p+0]=interface_rhs_y->item(j-1);
              all_send_data_U_y[comp][j][2*p+1]=interface_rhs_y->item(j);
            }
            else{
              all_send_data_U_y[comp][j][2*p+0]=interface_rhs_y->item(j-1);
              if(is_Uperiodic[1])
                all_send_data_U_y[comp][j][2*p+1]=interface_rhs_y->item(j); 
              else
                all_send_data_U_y[comp][j][2*p+1]=0;
            }

          }

          // Need to have the original rhs function assembled for corrosponding j,k pair
          double fe = assemble_local_rhs_y(i,k,gamma,t_it,comp);

          // Setup RHS = fi - Aie*xe for solving ui
          Aie->multiply_vec_then_add(interface_rhs_y,local_rhs_y,-1.0,1.0);

          // Solve ui and transfer solution into distributed vector
          GLOBAL_EQ->DS_NavierStokes_y_solver(i,k,min_unknown_index(1),local_rhs_y,interface_rhs_y,comp);

          

        }
      }
      else{
        for(i=min_unknown_index(0);i<=max_unknown_index(0);i++)
        {
          for(k=min_unknown_index(2);k<=max_unknown_index(2);k++)
          { 
            size_t nb_interface_unknowns = Vec_temp->nb_rows();
            for(j=0;j<nb_interface_unknowns;j++){
              Vec_temp->set_item(j,0);
              interface_rhs_y->set_item(j,0);
            }
            p = (k-min_unknown_index(2))+(max_unknown_index(2)-min_unknown_index(2)+1)*(i-min_unknown_index(0));
            if(is_Uperiodic[1])
                Vec_temp->set_item(nb_ranks_comm_i[1]-1,packed_data[3*p]);
            Vec_temp->set_item(0,packed_data[3*p+1]);
            interface_rhs_y->set_item(0,packed_data[3*p+2]);

            // Vec_temp might contain previous values

            for(j=1;j<nb_ranks_comm_i[1];j++){

              if(j!=nb_ranks_comm_i[1]-1){
                Vec_temp->add_to_item(j-1,all_receive_data_U_y[comp][j][3*p]);
                Vec_temp->add_to_item(j,all_receive_data_U_y[comp][j][3*p+1]);
                interface_rhs_y->set_item(j,all_receive_data_U_y[comp][j][3*p+2]);  // Assemble the interface rhs fe
              }
              else{
                if(is_Uperiodic[1] == 0){
                  Vec_temp->add_to_item(j-1,all_receive_data_U_y[comp][j][3*p]);
                }
                else{
                  Vec_temp->add_to_item(j-1,all_receive_data_U_y[comp][j][3*p]);
                  Vec_temp->add_to_item(j,all_receive_data_U_y[comp][j][3*p+1]);
                  interface_rhs_y->set_item(j,all_receive_data_U_y[comp][j][3*p+2]);
                }
              }
            }

            for(j=0;j<nb_interface_unknowns;j++){
            interface_rhs_y->set_item(j,interface_rhs_y->item(j)-Vec_temp->item(j)); // Get fe - Aei*xi to solve for ue
            }

            // Solve for ue (interface unknowns) in the master proc
            GLOBAL_EQ->DS_NavierStokes_y_interface_unknown_solver(interface_rhs_y,comp);

            //if(p==0 && rank_in_i[0] == 1)
            //  interface_rhs_y->print_items(MAC::out(),0);

            // Pack the interface_rhs_x into the appropriate send_data
            for (j=1;j<nb_ranks_comm_i[1];++j)
            {
              if(j!=nb_ranks_comm_i[1]-1){
                all_send_data_U_y[comp][j][2*p+0]=interface_rhs_y->item(j-1);
                all_send_data_U_y[comp][j][2*p+1]=interface_rhs_y->item(j);
              }
              else{
                all_send_data_U_y[comp][j][2*p+0]=interface_rhs_y->item(j-1);
                if(is_Uperiodic[1])
                  all_send_data_U_y[comp][j][2*p+1]=interface_rhs_y->item(j);
                else
                  all_send_data_U_y[comp][j][2*p+1]=0;
              }

            }

            // Need to have the original rhs function assembled for corrosponding j,k pair
            double fe = assemble_local_rhs_y(i,k,gamma,t_it,comp);

            // Setup RHS = fi - Aie*xe for solving ui
            Aie->multiply_vec_then_add(interface_rhs_y,local_rhs_y,-1.0,1.0);

            // Solve ui and transfer solution into distributed vector
            GLOBAL_EQ->DS_NavierStokes_y_solver(i,k,min_unknown_index(1),local_rhs_y,interface_rhs_y,comp);

          }
        }
      }

    }
   else
   {
      // Send the packed data to master
      //pelCOMM->send( is_master, packed_data, nb_received_data );
      MPI_Send( packed_data, nb_received_data, MPI_DOUBLE, 0, 0, DDS_Comm_i[1] ) ;
   }
   // Send the data from master
   if ( rank_in_i[1] == 0 )
   {
     for (i=1;i<nb_ranks_comm_i[1];++i)
      {
        //pelCOMM->send( i, all_send_data[i], nb_send_data );
        MPI_Send( all_send_data_U_y[comp][i], nb_send_data, MPI_DOUBLE, i, 0, DDS_Comm_i[1] ) ;
      }
   }
   else
   {
      // Receive the data
      //pelCOMM->receive( is_master, received_data, nb_received_data );
      static MPI_Status status ;
      MPI_Recv( packed_data, nb_received_data, MPI_DOUBLE, 0, 0,
                          DDS_Comm_i[1], &status ) ;

     // Solve the system of equations in each proc

     if(dim == 2)
     {
        for(i = min_unknown_index(0);i<=max_unknown_index(0);i++)
        {
          p = i-min_unknown_index(0);

          if(rank_in_i[1] != nb_ranks_comm_i[1]-1){
            interface_rhs_y->set_item(rank_in_i[1]-1,packed_data[2*p]);
            interface_rhs_y->set_item(rank_in_i[1],packed_data[2*p+1]);
          }
          else{
            if(is_Uperiodic[1] ==0){
              interface_rhs_y->set_item(rank_in_i[1]-1,packed_data[2*p]);
            }
            else{
              interface_rhs_y->set_item(rank_in_i[1]-1,packed_data[2*p]);
              interface_rhs_y->set_item(rank_in_i[1],packed_data[2*p+1]);
            }
          }
          // Need to have the original rhs function assembled for corrosponding j,k pair
          double fe = assemble_local_rhs_y(i,k,gamma,t_it,comp);

          // Setup RHS = fi - Aie*xe for solving ui
          Aie->multiply_vec_then_add(interface_rhs_y,local_rhs_y,-1.0,1.0);

          // Solve ui and transfer solution into distributed vector
          GLOBAL_EQ->DS_NavierStokes_y_solver(i,k,min_unknown_index(1),local_rhs_y,interface_rhs_y,comp);


        }
     }
     else
     {
        for(i = min_unknown_index(0);i<=max_unknown_index(0);i++)
        {
          for(k = min_unknown_index(2);k<=max_unknown_index(2);k++)
          {
            p = (k-min_unknown_index(2))+(max_unknown_index(2)-min_unknown_index(2)+1)*(i-min_unknown_index(0));
            if(rank_in_i[1] != nb_ranks_comm_i[1]-1){
              interface_rhs_y->set_item(rank_in_i[1]-1,packed_data[2*p]);
              interface_rhs_y->set_item(rank_in_i[1],packed_data[2*p+1]);
            }
            else{
              if(is_Uperiodic[1] ==0){
                interface_rhs_y->set_item(rank_in_i[1]-1,packed_data[2*p]);
              }
              else{
                interface_rhs_y->set_item(rank_in_i[1]-1,packed_data[2*p]);
                interface_rhs_y->set_item(rank_in_i[1],packed_data[2*p+1]);
              }
            }
            // Need to have the original rhs function assembled for corrosponding j,k pair
            double fe = assemble_local_rhs_y(i,k,gamma,t_it,comp);

            // Setup RHS = fi - Aie*xe for solving ui
            Aie->multiply_vec_then_add(interface_rhs_y,local_rhs_y,-1.0,1.0);

            // Solve ui and transfer solution into distributed vector
            GLOBAL_EQ->DS_NavierStokes_y_solver(i,k,min_unknown_index(1),local_rhs_y,interface_rhs_y,comp);
          }
        }
     }

   }
}




//---------------------------------------------------------------------------
double
DDS_NavierStokes:: assemble_local_rhs_z ( size_t const& i, size_t const& j, double gamma, FV_TimeIterator const* t_it, size_t const& comp )
//---------------------------------------------------------------------------
{
  // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     min_unknown_index(l) =
      UF->get_min_index_unknown_handled_by_proc( comp, l ) ;
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     max_unknown_index(l) =
      UF->get_max_index_unknown_handled_by_proc( comp, l ) ;

size_t k,pos;
int m;


double dzC,zhr,zhl,zright,zleft;
double fe=0.;

// Vector for fi
LA_SeqVector* local_rhs_z = GLOBAL_EQ->get_local_temp_z(comp);

 // Compute VEC_rhs_z = rhs in z
 for (k=min_unknown_index(2);k<=max_unknown_index(2);++k)
 {
   double zvalue=0.;

   zhr= UF->get_DOF_coordinate( k+1,comp, 2 ) - UF->get_DOF_coordinate( k, comp, 2 ) ;
   zhl= UF->get_DOF_coordinate( k, comp, 2 ) - UF->get_DOF_coordinate( k-1, comp, 2 ) ;
   zright = UF->DOF_value( i, j, k+1, comp, 1 ) - UF->DOF_value( i, j, k, comp, 1 ) ;
   zleft = UF->DOF_value( i, j, k, comp, 1 ) - UF->DOF_value( i, j, k-1, comp, 1 ) ;

   /** Assemble rhs in z for Neumann/Dirichlet BC Conditions */
   if(UF->DOF_in_domain(i, j, k-1, comp) && UF->DOF_in_domain(i, j, k+1, comp))
     zvalue = zright/zhr - zleft/zhl;
   else if(UF->DOF_in_domain(i, j, k-1, comp))
     zvalue = - zleft/zhl;
   else
     zvalue = zright/zhr;

   dzC = UF->get_cell_size( k, comp, 2 ) ;

   pos = k - min_unknown_index(2);
   if(is_Uperiodic[2] == 0){
      if(rank_in_i[2] == nb_ranks_comm_i[2]-1){
        local_rhs_z->set_item( pos, (UF->DOF_value( i, j, k, comp, 0 )*rho*dzC)/(t_it -> time_step()) - gamma*zvalue );
      }
      else{
         if(k == max_unknown_index(2))
           fe = (UF->DOF_value( i, j, k, comp, 0 )*rho*dzC)/(t_it -> time_step()) - gamma*zvalue ;
         else
           local_rhs_z->set_item( pos, (UF->DOF_value( i, j, k, comp, 0 )*rho*dzC)/(t_it -> time_step()) - gamma*zvalue );
      }
   }
   else{
      if(nb_ranks_comm_i[2] > 1){
        if(k == max_unknown_index(2))
           fe = (UF->DOF_value( i, j, k, comp, 0 )*rho*dzC)/(t_it -> time_step()) - gamma*zvalue ;
        else
           local_rhs_z->set_item( pos, (UF->DOF_value( i, j, k, comp, 0 )*rho*dzC)/(t_it -> time_step()) - gamma*zvalue );
      }
      else
        local_rhs_z->set_item( pos, (UF->DOF_value( i, j, k, comp, 0 )*rho*dzC)/(t_it -> time_step()) - gamma*zvalue );
   }

    
 }

 m = int(min_unknown_index(2)) - 1;
 if ( UF->DOF_in_domain( i, j, m, comp ) )
   if ( UF->DOF_has_imposed_Dirichlet_value( i, j, m, comp ) )
   {
     double ai = 1/(UF->get_DOF_coordinate( m+1,comp,2) - UF->get_DOF_coordinate( m,comp,2));
     double dirichlet_value = UF->DOF_value( i, j, m, comp, 1 ) ;
     local_rhs_z->add_to_item( 0, + gamma * ai * dirichlet_value );
   }

 m = int(max_unknown_index(2)) + 1;
 if ( UF->DOF_in_domain( i, j, m, comp ) )
   if ( UF->DOF_has_imposed_Dirichlet_value( i, j, m, comp ) )
   {
     double ai = 1/(UF->get_DOF_coordinate( m,comp,2) - UF->get_DOF_coordinate( m-1,comp,2));
     double dirichlet_value = UF->DOF_value( i, j, m, comp, 1 ) ;
     local_rhs_z->add_to_item( local_rhs_z->nb_rows()-1, + gamma * ai * dirichlet_value );
   }

return fe;
}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: solve_z_for_secondorder ( size_t const& i, size_t const& j, double gamma, FV_TimeIterator const* t_it,double * packed_data, size_t const& comp )
//---------------------------------------------------------------------------
{
  // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     min_unknown_index(l) =
      UF->get_min_index_unknown_handled_by_proc( comp, l ) ;
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     max_unknown_index(l) =
      UF->get_max_index_unknown_handled_by_proc( comp, l ) ;

   double fe = assemble_local_rhs_z(i,j,gamma,t_it,comp);

   size_t k;
   // create a replica of local rhs vector in local solution vector
   LA_SeqVector* local_rhs_z = GLOBAL_EQ->get_local_temp_z(comp) ;

   LA_SeqVector* local_solution_z = GLOBAL_EQ->get_local_solution_temp_z(comp) ;
   for(k=0;k<local_rhs_z->nb_rows();k++){
    local_solution_z->set_item(k,local_rhs_z->item(k));
   }

   // Solve for xi locally and put it in local solution vector
   GLOBAL_EQ->DS_NavierStokes_z_local_unknown_solver(local_solution_z,comp);

   LA_SeqMatrix* Aei = GLOBAL_EQ-> get_aei(comp,2);
   LA_SeqVector* Vec_temp = GLOBAL_EQ-> get_temp_z(comp);

   for(k=0;k<Vec_temp->nb_rows();k++){
          Vec_temp->set_item(k,0);
    }
   // Calculate Aei*xi in each proc locally
   Aei->multiply_vec_then_add(local_solution_z,Vec_temp);


    size_t vec_pos=(i-min_unknown_index(0))+(max_unknown_index(0)-min_unknown_index(0)+1)*(j-min_unknown_index(1));

    if(rank_in_i[2] == 0){
      // Check if bc is periodic in x
     // If it is, we need to pack two elements apart from fe

      if(is_Uperiodic[2])
          packed_data[3*vec_pos+0]=Vec_temp->item(nb_ranks_comm_i[2]-1);
      else
          packed_data[3*vec_pos+0]=0;
      packed_data[3*vec_pos+1]=Vec_temp->item(rank_in_i[2]);
    }
    else if(rank_in_i[2] == nb_ranks_comm_i[2]-1){
      packed_data[3*vec_pos+0]=Vec_temp->item(rank_in_i[2]-1);
      // Check if bc is periodic in x
      // If it is, we need to pack two elements apart from fe
      if(is_Uperiodic[2])
          packed_data[3*vec_pos+1]=Vec_temp->item(rank_in_i[2]);
      else
          packed_data[3*vec_pos+1]=0;
    }
    else{
      packed_data[3*vec_pos+0]=Vec_temp->item(rank_in_i[2]-1);
      packed_data[3*vec_pos+1]=Vec_temp->item(rank_in_i[2]);
    }
    packed_data[3*vec_pos+2] = fe; // Send the fe values and 0 for last proc
}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: solve_interface_unknowns_z ( double * packed_data, size_t nb_received_data, double gamma,  FV_TimeIterator const* t_it, size_t const& comp  )
//---------------------------------------------------------------------------
{
   size_t i,j,p,m;

  // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     min_unknown_index(l) =
      UF->get_min_index_unknown_handled_by_proc( comp, l ) ;
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     max_unknown_index(l) =
      UF->get_max_index_unknown_handled_by_proc( comp, l ) ;

   // Vector for fe
   LA_SeqVector* interface_rhs_z = GLOBAL_EQ->get_interface_temp_z(comp) ;
   // Vector for fi
   LA_SeqVector* local_rhs_z = GLOBAL_EQ->get_local_temp_z(comp);

   LA_SeqVector* Vec_temp = GLOBAL_EQ-> get_temp_z(comp);

   LA_SeqMatrix* Aie = GLOBAL_EQ-> get_aie(comp,2);

   size_t nb_send_data = (2*nb_received_data)/3;
   
  // Send and receive the data first pass
   if ( rank_in_i[2] == 0 )
   {
      for (p=1;p<nb_ranks_comm_i[2];++p)
      {
        // Receive the data
        //pelCOMM->receive( i, all_received_data[i], nb_received_data );
        static MPI_Status status ;
        MPI_Recv( all_receive_data_U_z[comp][p], nb_received_data, MPI_DOUBLE, p, 0,
                          DDS_Comm_i[2], &status ) ;
      }

      // Solve system of interface unknowns for each x

      for(j=min_unknown_index(1);j<=max_unknown_index(1);j++)
      {
        for(i=min_unknown_index(0);i<=max_unknown_index(0);i++)
        { 
          size_t nb_interface_unknowns = Vec_temp->nb_rows();
          for(m=0;m<nb_interface_unknowns;m++){
            Vec_temp->set_item(m,0);
            interface_rhs_z->set_item(m,0);
          }
          p = (i-min_unknown_index(0))+(max_unknown_index(0)-min_unknown_index(0)+1)*(j-min_unknown_index(1));
          if(is_Uperiodic[2])
              Vec_temp->set_item(nb_ranks_comm_i[2]-1,packed_data[3*p]);
          Vec_temp->set_item(0,packed_data[3*p+1]);
          interface_rhs_z->set_item(0,packed_data[3*p+2]);

          // Vec_temp might contain previous values

          for(m=1;m<nb_ranks_comm_i[2];m++){

            if(m!=nb_ranks_comm_i[2]-1){
              Vec_temp->add_to_item(m-1,all_receive_data_U_z[comp][m][3*p]);
              Vec_temp->add_to_item(m,all_receive_data_U_z[comp][m][3*p+1]);
              interface_rhs_z->set_item(m,all_receive_data_U_z[comp][m][3*p+2]);  // Assemble the interface rhs fe
            }
            else{
              if(is_Uperiodic[2] == 0)
                  Vec_temp->add_to_item(m-1,all_receive_data_U_z[comp][m][3*p]);
              else{
                  Vec_temp->add_to_item(m-1,all_receive_data_U_z[comp][m][3*p]);
                  Vec_temp->add_to_item(m,all_receive_data_U_z[comp][m][3*p+1]);
                  interface_rhs_z->set_item(m,all_receive_data_U_z[comp][m][3*p+2]);
              }
            }
          }

          for(m=0;m<nb_interface_unknowns;m++){
          interface_rhs_z->set_item(m,interface_rhs_z->item(m)-Vec_temp->item(m)); // Get fe - Aei*xi to solve for ue
          }

          // Solve for ue (interface unknowns) in the master proc
          GLOBAL_EQ->DS_NavierStokes_z_interface_unknown_solver(interface_rhs_z,comp);

          // Pack the interface_rhs_x into the appropriate send_data
          for (m=1;m<nb_ranks_comm_i[2];++m)
          {
            if(m!=nb_ranks_comm_i[2]-1){
              all_send_data_U_z[comp][m][2*p+0]=interface_rhs_z->item(m-1);
              all_send_data_U_z[comp][m][2*p+1]=interface_rhs_z->item(m);
            }
            else{
              all_send_data_U_z[comp][m][2*p+0]=interface_rhs_z->item(m-1);
              if(is_Uperiodic[2])
                all_send_data_U_z[comp][m][2*p+1]=interface_rhs_z->item(m);
              else
                all_send_data_U_z[comp][m][2*p+1]=0;
            }

          }

          // Need to have the original rhs function assembled for corrosponding j,k pair
          double fe = assemble_local_rhs_z(i,j,gamma,t_it,comp);

          // Setup RHS = fi - Aie*xe for solving ui
          Aie->multiply_vec_then_add(interface_rhs_z,local_rhs_z,-1.0,1.0);

          // Solve ui and transfer solution into distributed vector
          GLOBAL_EQ->DS_NavierStokes_z_solver(i,j,min_unknown_index(2),local_rhs_z,interface_rhs_z,comp);
        }
      }

    }
   else
   {
      // Send the packed data to master
      //pelCOMM->send( is_master, packed_data, nb_received_data );
      MPI_Send( packed_data, nb_received_data, MPI_DOUBLE, 0, 0, DDS_Comm_i[2] ) ;
   }
   // Send the data from master
   if ( rank_in_i[2] == 0 )
   {
     for (m=1;m<nb_ranks_comm_i[2];++m)
      {
        //pelCOMM->send( i, all_send_data[i], nb_send_data );
        MPI_Send( all_send_data_U_z[comp][m], nb_send_data, MPI_DOUBLE, m, 0, DDS_Comm_i[2] ) ;
      }
   }
   else
   {
      
      // Receive the data
      //pelCOMM->receive( is_master, received_data, nb_received_data );
      static MPI_Status status ;
      MPI_Recv( packed_data, nb_received_data, MPI_DOUBLE, 0, 0,
                          DDS_Comm_i[2], &status ) ;

     // Solve the system of equations in each proc
     for(j = min_unknown_index(1);j<=max_unknown_index(1);j++){
        for(i = min_unknown_index(0);i<=max_unknown_index(0);i++){
          p = (i-min_unknown_index(0))+(max_unknown_index(0)-min_unknown_index(0)+1)*(j-min_unknown_index(1));
          if(rank_in_i[2] != nb_ranks_comm_i[2]-1){
            interface_rhs_z->set_item(rank_in_i[2]-1,packed_data[2*p]);
            interface_rhs_z->set_item(rank_in_i[2],packed_data[2*p+1]);
          }
          else{
            if(is_Uperiodic[2] == 0)
                interface_rhs_z->set_item(rank_in_i[2]-1,packed_data[2*p]);    
            else{
                interface_rhs_z->set_item(rank_in_i[2]-1,packed_data[2*p]);
                interface_rhs_z->set_item(rank_in_i[2],packed_data[2*p+1]);
            }
          }
          // Need to have the original rhs function assembled for corrosponding j,k pair
          double fe = assemble_local_rhs_z(i,j,gamma,t_it,comp);

          // Setup RHS = fi - Aie*xe for solving ui
          Aie->multiply_vec_then_add(interface_rhs_z,local_rhs_z,-1.0,1.0);

          // Solve ui and transfer solution into distributed vector
          GLOBAL_EQ->DS_NavierStokes_z_solver(i,j,min_unknown_index(2),local_rhs_z,interface_rhs_z,comp);
        }
     }
   }
}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: NS_velocity_update ( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: NS_velocity_update" ) ;

   double gamma=mu, ugradu = 0.;
   size_t i, j, k;

  double dxC,xC,dyC,yC,xhr,xhl,xright,xleft,yhr,yhl,yright,yleft;
  double pvalue =0., xvalue=0.,yvalue=0.,rhs=0.,bodyterm=0.0,advection_value = 0.;
  FV_SHIFT_TRIPLET shift ;
  int cpp=-1;

  // Periodic pressure gradient
  if ( UF->primary_grid()->is_periodic_flow() )
  {
    cpp = UF->primary_grid()->get_periodic_flow_direction() ;
    bodyterm = UF->primary_grid()->get_periodic_pressure_drop() /
            ( UF->primary_grid()->get_main_domain_max_coordinate( cpp )
            - UF->primary_grid()->get_main_domain_min_coordinate( cpp ) ) ;
  }

  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);
  // First Equation

  //if ( my_rank == is_master ) SCT_set_start("Explicit Velocity step");

  for(size_t comp=0;comp<nb_comps;comp++)
  {
    // Get local min and max indices
    for (size_t l=0;l<dim;++l)
      min_unknown_index(l) =
       UF->get_min_index_unknown_handled_by_proc( comp, l ) ;
    for (size_t l=0;l<dim;++l)
      max_unknown_index(l) =
       UF->get_max_index_unknown_handled_by_proc( comp, l ) ;

    // Here we use shift_staggeredToStaggered to shitf centered to staggered
    // for the following reason: the ith component of UU is staggered with the 
    // centered field only in the ith direction, which shares its ith location 
    // with the non-ith components of the velocity field
    // Example: 
    // For ux, it is staggered with the centered field tf in the x
    // direction only, in the y & z direction, ux and tf have the same location
    // Now, in the x direction, tf is located at the same x position as uy and
    // uz and hence the shift_staggeredToStaggered can be used for tf
    // When interpolating a centered field to a staggered field, we use 
    // shift_staggeredToStaggered for each ith component and consider the ith 
    // shift in the ith direction only, i.e.:
    // * for ux, use shift_staggeredToStaggered(0) and shift in the x direction
    // with shift.i (the xth component of the shift) only
    // * for uy, use shift_staggeredToStaggered(1) and shift in the y direction
    // with shift.j (the yth component of the shift) only
    // * for uz, use shift_staggeredToStaggered(2) and shift in the z direction
    // with shift.k (the zth component of the shift) only    
    shift = UF->shift_staggeredToStaggered( comp ) ;   

    for (i=min_unknown_index(0);i<=max_unknown_index(0);++i)
    {
      // Compute VEC_rhs_x = rhs in x
      for (j=min_unknown_index(1);j<=max_unknown_index(1);++j)
      {

        if(dim ==2 )
        {

          k=0;

          // Dxx for un
          xhr= UF->get_DOF_coordinate( i+1,comp, 0 ) - UF->get_DOF_coordinate( i, comp, 0 ) ;
          xhl= UF->get_DOF_coordinate( i, comp, 0 ) - UF->get_DOF_coordinate( i-1, comp, 0 ) ;
          xright = UF->DOF_value( i+1, j, 0, comp, 1 ) - UF->DOF_value( i, j, 0, comp, 1 ) ;
          xleft = UF->DOF_value( i, j, 0, comp, 1 ) - UF->DOF_value( i-1, j, 0, comp, 1 ) ;

          /** Assemble rhs in x for Neumann/Dirichlet BC Conditions */
          if(UF->DOF_in_domain( i-1, j, k, comp) && UF->DOF_in_domain( i+1, j, k, comp))
              xvalue = xright/xhr - xleft/xhl;
          else if (UF->DOF_in_domain( i-1, j, k, comp))
            xvalue = - xleft/xhl;
          else
            xvalue = xright/xhr;

          // Dyy for un
          yhr= UF->get_DOF_coordinate( j+1,comp, 1 ) - UF->get_DOF_coordinate( j, comp, 1 ) ;
          yhl= UF->get_DOF_coordinate( j, comp, 1 ) - UF->get_DOF_coordinate( j-1, comp, 1 ) ;
          yright = UF->DOF_value( i, j+1, 0, comp, 1 ) - UF->DOF_value( i, j, 0, comp, 1 ) ;
          yleft = UF->DOF_value( i, j, 0, comp, 1 ) - UF->DOF_value( i, j-1, 0, comp, 1 ) ;

          /** Assemble rhs in y for Neumann/Dirichlet BC Conditions */
          //yvalue = yright/yhr - yleft/yhl;
          if(UF->DOF_in_domain(i, j-1, k, comp) && UF->DOF_in_domain(i, j+1, k, comp))
            yvalue = yright/yhr - yleft/yhl;
          else if(UF->DOF_in_domain(i, j-1, k, comp))
            yvalue = - yleft/yhl;
          else
            yvalue = yright/yhr;

          // Assemble the bodyterm
          dxC = UF->get_cell_size( i, comp, 0 ) ;
          xC = UF->get_DOF_coordinate( i, comp, 0 ) ;
          dyC = UF->get_cell_size( j, comp, 1 ) ;
          yC = UF->get_DOF_coordinate( j, comp, 1 ) ;

          // if(b_bodyterm)
          // {
          //   bodyterm = 2. * MAC::pi() * MAC::pi() * MAC::sin( MAC::pi() * xC )* MAC::sin( MAC::pi() * yC );
          // }

          // Pressure contribution
          if(comp == 0){
            pvalue = (PF->DOF_value( shift.i+i, j, 0, 0, 1 ) - PF->DOF_value( shift.i+i-1, j, 0, 0, 1 ))*dyC;
          }
          else{
            pvalue = (PF->DOF_value( i, shift.j+j, 0, 0, 1 ) - PF->DOF_value( i, shift.j+j-1, 0, 0, 1 ))*dxC;
          }
          
      	  // Advection term
      	  if ( AdvectionScheme == "TVD" )
      	    ugradu = assemble_advection_TVD(1,rho,1,i,j,k,comp);
      	  else
      	    ugradu = assemble_advection_Upwind(1,rho,1,i,j,k,comp);
      	  
      	  if ( AdvectionTimeAccuracy == 1 )
      	    advection_value = ugradu;
      	  else
      	  {
      	    advection_value = 1.5*ugradu - 0.5*UF->DOF_value(i,j,k,comp,2);
                  UF->set_DOF_value(i,j,k,comp,2,ugradu);
      	  }

          rhs = gamma*(xvalue*dyC + yvalue*dxC) - pvalue - advection_value
            + (UF->DOF_value( i, j, k, comp, 1 )*dxC*dyC*rho)/(t_it -> time_step());

          if ( cpp >= 0 && cpp==comp ) rhs += - bodyterm*dxC*dyC;  

          UF->set_DOF_value( i, j, k, comp, 0, rhs*(t_it -> time_step())/(dxC*dyC*rho));

        }
        else
        {

          for (k=min_unknown_index(2);k<=max_unknown_index(2);++k)
          {

            // Dxx for un
            xhr= UF->get_DOF_coordinate( i+1,comp, 0 ) - UF->get_DOF_coordinate( i, comp, 0 ) ;
            xhl= UF->get_DOF_coordinate( i, comp, 0 ) - UF->get_DOF_coordinate( i-1, comp, 0 ) ;
            xright = UF->DOF_value( i+1, j, k, comp, 1 ) - UF->DOF_value( i, j, k, comp, 1 ) ;
            xleft = UF->DOF_value( i, j, k, comp, 1 ) - UF->DOF_value( i-1, j, k, comp, 1 ) ;

            /** Assemble rhs in x for Neumann/Dirichlet BC Conditions */
            if(UF->DOF_in_domain( i-1, j, k, comp) && UF->DOF_in_domain( i+1, j, k, comp))
                xvalue = xright/xhr - xleft/xhl;
            else if (UF->DOF_in_domain( i-1, j, k, comp))
              xvalue = - xleft/xhl;
            else
              xvalue = xright/xhr;

            // Dyy for un
            yhr= UF->get_DOF_coordinate( j+1,comp, 1 ) - UF->get_DOF_coordinate( j, comp, 1 ) ;
            yhl= UF->get_DOF_coordinate( j, comp, 1 ) - UF->get_DOF_coordinate( j-1, comp, 1 ) ;
            yright = UF->DOF_value( i, j+1, k, comp, 1 ) - UF->DOF_value( i, j, k, comp, 1 ) ;
            yleft = UF->DOF_value( i, j, k, comp, 1 ) - UF->DOF_value( i, j-1, k, comp, 1 ) ;

            /** Assemble rhs in y for Neumann/Dirichlet BC Conditions */
            if(UF->DOF_in_domain(i, j-1, k, comp) && UF->DOF_in_domain(i, j+1, k, comp))
              yvalue = yright/yhr - yleft/yhl;
            else if(UF->DOF_in_domain(i, j-1, k, comp))
              yvalue = - yleft/yhl;
            else
              yvalue = yright/yhr;

            // Assemble the bodyterm
            dxC = UF->get_cell_size( i, comp, 0 ) ;
            xC = UF->get_DOF_coordinate( i, comp, 0 ) ;
            dyC = UF->get_cell_size( j, comp, 1 ) ;
            yC = UF->get_DOF_coordinate( j, comp, 1 ) ;

            double dzC,zC,zhr,zhl,zright,zleft,zvalue=0.;
            zhr= UF->get_DOF_coordinate( k+1,comp, 2 ) - UF->get_DOF_coordinate( k, comp, 2 ) ;
            zhl= UF->get_DOF_coordinate( k, comp, 2 ) - UF->get_DOF_coordinate( k-1, comp, 2 ) ;
            zright = UF->DOF_value( i, j, k+1, comp, 1 ) - UF->DOF_value( i, j, k, comp, 1 ) ;
            zleft = UF->DOF_value( i, j, k, comp, 1 ) - UF->DOF_value( i, j, k-1, comp, 1 ) ;

            /** Assemble rhs in z for Neumann/Dirichlet BC Conditions */
            if(UF->DOF_in_domain(i, j, k-1, comp) && UF->DOF_in_domain(i, j, k+1, comp))
              zvalue = zright/zhr - zleft/zhl;
            else if(UF->DOF_in_domain(i, j, k-1, comp))
              zvalue = - zleft/zhl;
            else
              zvalue = zright/zhr;


            dzC = UF->get_cell_size( k, comp, 2 ) ;
            zC = UF->get_DOF_coordinate( k, comp, 2 ) ;

            // Pressure contribution
            if(comp == 0){
              pvalue = (PF->DOF_value( shift.i+i, j, k, 0, 1 ) - PF->DOF_value( shift.i+i-1, j, k, 0, 1 ))*dyC*dzC;
            }
            else if(comp==1){
              pvalue = (PF->DOF_value( i, shift.j+j, k, 0, 1 ) - PF->DOF_value( i, shift.j+j-1, k, 0, 1 ))*dxC*dzC;
            }
            else{
              pvalue = (PF->DOF_value( i, j, shift.k+k, 0, 1 ) - PF->DOF_value( i, j, shift.k+k-1, 0, 1 ))*dxC*dyC;
            }

      	    // Advection term
      	    if ( AdvectionScheme == "TVD" )
      	      ugradu = assemble_advection_TVD(1,rho,1,i,j,k,comp);
      	    else
      	      ugradu = assemble_advection_Upwind(1,rho,1,i,j,k,comp);
      	  
      	    if ( AdvectionTimeAccuracy == 1 )
      	      advection_value = ugradu;
      	    else
      	    {
      	      advection_value = 1.5*ugradu - 0.5*UF->DOF_value(i,j,k,comp,2);
                    UF->set_DOF_value(i,j,k,comp,2,ugradu);
      	    }

            rhs = gamma*(xvalue*dyC*dzC + yvalue*dxC*dzC + zvalue*dxC*dyC) - pvalue - advection_value + (UF->DOF_value( i, j, k, comp, 1 )*dxC*dyC*dzC*rho)/(t_it -> time_step());

            if ( cpp >= 0 && cpp==comp ) rhs += - bodyterm*dxC*dyC*dzC;

            UF->set_DOF_value( i, j, k, comp, 0, rhs*(t_it -> time_step())/(dxC*dyC*dzC*rho));
          }
        }
      }

    }
  }

  //if ( my_rank == is_master ) SCT_get_elapsed_time( "Explicit Velocity step" );

  // Update gamma based for invidual direction
  gamma = mu/2.0;

  if(dim == 2){

    for(size_t comp=0;comp<nb_comps;comp++){

      // Get local min and max indices
      for (size_t l=0;l<dim;++l)
        min_unknown_index(l) =
         UF->get_min_index_unknown_handled_by_proc( comp, l ) ;
      for (size_t l=0;l<dim;++l)
        max_unknown_index(l) =
         UF->get_max_index_unknown_handled_by_proc( comp, l ) ;

      size_t nb_send_data_x = 3*(max_unknown_index(1)-min_unknown_index(1)+1);
      // Solve in x
      if(nb_ranks_comm_i[0]>1){

        for (j=min_unknown_index(1);j<=max_unknown_index(1);++j)
         {
            k=0;
            solve_x_for_secondorder(j,k,gamma,t_it,mpi_packed_data_U_x[comp],comp);
         }

         solve_interface_unknowns_x ( mpi_packed_data_U_x[comp], nb_send_data_x, gamma, t_it,comp );

      }
      else{
        for (j=min_unknown_index(1);j<=max_unknown_index(1);++j)
         {
            k=0;

            double fe = assemble_local_rhs_x(j,k,gamma,t_it,comp);
            LA_SeqVector* local_rhs_x = GLOBAL_EQ->get_local_temp_x(comp);
            
            if(is_Uperiodic[0] == 1)
                GLOBAL_EQ->DS_NavierStokes_x_solver_periodic(j,k,min_unknown_index(0),local_rhs_x,NULL,comp); 
            else
                GLOBAL_EQ->DS_NavierStokes_x_solver(j,k,min_unknown_index(0),local_rhs_x,NULL,comp);                            
          }
      }


    }

     // Synchronize the distributed DS solution vector
     GLOBAL_EQ->synchronize_DS_solution_vec();

     // Tranfer back to field
     UF->update_free_DOFs_value( 0, GLOBAL_EQ->get_solution_DS_velocity() ) ;

     for(size_t comp=0;comp<nb_comps;comp++){

      // Get local min and max indices
      for (size_t l=0;l<dim;++l)
        min_unknown_index(l) =
         UF->get_min_index_unknown_handled_by_proc( comp, l ) ;
      for (size_t l=0;l<dim;++l)
        max_unknown_index(l) =
         UF->get_max_index_unknown_handled_by_proc( comp, l ) ;

      size_t nb_send_data_y = 3*(max_unknown_index(0)-min_unknown_index(0)+1);

      // Solve in y
       if(nb_ranks_comm_i[1]>1){

         for (i=min_unknown_index(0);i<=max_unknown_index(0);++i)
         {
            k=0;
            solve_y_for_secondorder(i,k,gamma,t_it,mpi_packed_data_U_y[comp],comp);
         }

         solve_interface_unknowns_y ( mpi_packed_data_U_y[comp], nb_send_data_y, gamma, t_it,comp );
       }
       else{
        for (i=min_unknown_index(0);i<=max_unknown_index(0);++i)
         {
            k=0;
           double fe = assemble_local_rhs_y(i,k,gamma,t_it,comp);
           LA_SeqVector* local_rhs_y = GLOBAL_EQ->get_local_temp_y(comp);

           if(is_Uperiodic[1] == 1)
              GLOBAL_EQ->DS_NavierStokes_y_solver_periodic(i,k,min_unknown_index(1),local_rhs_y,NULL,comp);
           else
              GLOBAL_EQ->DS_NavierStokes_y_solver(i,k,min_unknown_index(1),local_rhs_y,NULL,comp);
         }
       }

     }

     // // Synchronize the distributed DS solution vector
     GLOBAL_EQ->synchronize_DS_solution_vec();

     // Tranfer back to field
     UF->update_free_DOFs_value( 0 , GLOBAL_EQ->get_solution_DS_velocity() ) ; 
     
   }
  else{
    
    for(size_t comp=0;comp<nb_comps;comp++)
    {

      // Get local min and max indices
      for (size_t l=0;l<dim;++l)
        min_unknown_index(l) =
         UF->get_min_index_unknown_handled_by_proc( comp, l ) ;
      for (size_t l=0;l<dim;++l)
        max_unknown_index(l) =
         UF->get_max_index_unknown_handled_by_proc( comp, l ) ;
      size_t nb_send_data_x = 3*(max_unknown_index(1)-min_unknown_index(1)+1)*(max_unknown_index(2)-min_unknown_index(2)+1);

      // Solve in x
      if(nb_ranks_comm_i[0]>1){
        
        for (k=min_unknown_index(2);k<=max_unknown_index(2);++k){
          for (j=min_unknown_index(1);j<=max_unknown_index(1);++j){
              solve_x_for_secondorder(j,k,gamma,t_it,mpi_packed_data_U_x[comp],comp);
          }
        }
        solve_interface_unknowns_x ( mpi_packed_data_U_x[comp], nb_send_data_x, gamma, t_it,comp );
      }
      else{
        
        for (k=min_unknown_index(2);k<=max_unknown_index(2);++k){
          for (j=min_unknown_index(1);j<=max_unknown_index(1);++j){
//              if ( my_rank == is_master ) SCT_set_start("Velocity x update");
              double fe = assemble_local_rhs_x(j,k,gamma,t_it,comp);
              LA_SeqVector* local_rhs_x = GLOBAL_EQ->get_local_temp_x(comp);
              if(is_Uperiodic[0] == 1)
                  GLOBAL_EQ->DS_NavierStokes_x_solver_periodic(j,k,min_unknown_index(0),local_rhs_x,NULL,comp);
              else
                  GLOBAL_EQ->DS_NavierStokes_x_solver(j,k,min_unknown_index(0),local_rhs_x,NULL,comp);
//              if ( my_rank == is_master ) SCT_get_elapsed_time( "Velocity x update" );
          }
        }
      }
    }

    // Synchronize the distributed DS solution vector
    GLOBAL_EQ->synchronize_DS_solution_vec();

    // Tranfer back to field
    UF->update_free_DOFs_value( 0, GLOBAL_EQ->get_solution_DS_velocity() ) ;

    for(size_t comp=0;comp<nb_comps;comp++)
    {

      // Get local min and max indices
      for (size_t l=0;l<dim;++l)
        min_unknown_index(l) =
         UF->get_min_index_unknown_handled_by_proc( comp, l ) ;
      for (size_t l=0;l<dim;++l)
        max_unknown_index(l) =
         UF->get_max_index_unknown_handled_by_proc( comp, l ) ;
      size_t nb_send_data_y = 3*(max_unknown_index(0)-min_unknown_index(0)+1)*(max_unknown_index(2)-min_unknown_index(2)+1);
        
      // Solve in y
      if(nb_ranks_comm_i[1]>1){
        
        // Third equation - Solve in y
        for (i=min_unknown_index(0);i<=max_unknown_index(0);++i)
         {
           for (k=min_unknown_index(2);k<=max_unknown_index(2);++k)
           {
                solve_y_for_secondorder(i,k,gamma,t_it,mpi_packed_data_U_y[comp],comp);
           }
         }

         solve_interface_unknowns_y ( mpi_packed_data_U_y[comp], nb_send_data_y, gamma, t_it,comp );
      }
      else{
        for (i=min_unknown_index(0);i<=max_unknown_index(0);++i)
         {
           for (k=min_unknown_index(2);k<=max_unknown_index(2);++k)
           {
              //if ( my_rank == is_master ) SCT_set_start("Velocity y update");
              double fe = assemble_local_rhs_y(i,k,gamma,t_it,comp);
              LA_SeqVector* local_rhs_y = GLOBAL_EQ->get_local_temp_y(comp);
              if(is_Uperiodic[1] == 1)
                  GLOBAL_EQ->DS_NavierStokes_y_solver_periodic(i,k,min_unknown_index(1),local_rhs_y,NULL,comp);
              else
                  GLOBAL_EQ->DS_NavierStokes_y_solver(i,k,min_unknown_index(1),local_rhs_y,NULL,comp);
              //if ( my_rank == is_master ) SCT_get_elapsed_time( "Velocity y update" );
           }
         }
      }

    }

     // Synchronize the distributed DS solution vector
     GLOBAL_EQ->synchronize_DS_solution_vec();

     // Tranfer back to field
     UF->update_free_DOFs_value( 0 , GLOBAL_EQ->get_solution_DS_velocity() ) ;

    
    for(size_t comp=0;comp<nb_comps;comp++)
    {

      // Get local min and max indices
      for (size_t l=0;l<dim;++l)
        min_unknown_index(l) =
         UF->get_min_index_unknown_handled_by_proc( comp, l ) ;
      for (size_t l=0;l<dim;++l)
        max_unknown_index(l) =
         UF->get_max_index_unknown_handled_by_proc( comp, l ) ;

       // Solve in z
       size_t nb_send_data_z = 3*(max_unknown_index(0)-min_unknown_index(0)+1)*(max_unknown_index(1)-min_unknown_index(1)+1);
         
       if(nb_ranks_comm_i[2]>1){
         
         for (j=min_unknown_index(1);j<=max_unknown_index(1);++j)
         {
           for (i=min_unknown_index(0);i<=max_unknown_index(0);++i)
           {
             solve_z_for_secondorder(i,j,gamma,t_it,mpi_packed_data_U_z[comp],comp);
           }
         }

         solve_interface_unknowns_z ( mpi_packed_data_U_z[comp], nb_send_data_z, gamma, t_it,comp );
       }
       else{
         for (j=min_unknown_index(1);j<=max_unknown_index(1);++j)
         {
           for (i=min_unknown_index(0);i<=max_unknown_index(0);++i)
           {
              //if ( my_rank == is_master ) SCT_set_start("Velocity z update");
              double fe = assemble_local_rhs_z(i,j,gamma,t_it,comp);
              LA_SeqVector* local_rhs_z = GLOBAL_EQ->get_local_temp_z(comp);
              if(is_Uperiodic[2] == 1)
                  GLOBAL_EQ->DS_NavierStokes_z_solver_periodic(i,j,min_unknown_index(2),local_rhs_z,NULL,comp);
              else
                  GLOBAL_EQ->DS_NavierStokes_z_solver(i,j,min_unknown_index(2),local_rhs_z,NULL,comp);
              //if ( my_rank == is_master ) SCT_get_elapsed_time( "Velocity z update" );
           }
         }
       }

    }

     // Synchronize the distributed DS solution vector
     GLOBAL_EQ->synchronize_DS_solution_vec();

     // Tranfer back to field
     UF->update_free_DOFs_value( 0, GLOBAL_EQ->get_solution_DS_velocity() ) ;

  }

  output_L2norm_velocity(0);
  output_L2norm_velocity(1);
  output_L2norm_velocity(2);

}




//---------------------------------------------------------------------------
double
DDS_NavierStokes:: pressure_assemble_local_rhs_x ( size_t const& j, size_t const& k, FV_TimeIterator const* t_it)
//---------------------------------------------------------------------------
{
  // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     min_unknown_index(l) =
      PF->get_min_index_unknown_handled_by_proc( 0, l ) ;
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     max_unknown_index(l) =
      PF->get_max_index_unknown_handled_by_proc( 0, l ) ;

  size_t i,pos,m;
  FV_SHIFT_TRIPLET shift = PF->shift_staggeredToCentered() ;

  // Compute VEC_rhs_x = rhs in x
  double xhr,xright,yhr,yright,dx;
  double fe=0.;
  double xvalue = 0.,yvalue=0.,value=0.;

  // Vector for fi
  LA_SeqVector* local_rhs_x = GLOBAL_EQ->get_local_temp_x_P() ;

  for (i=min_unknown_index(0);i<=max_unknown_index(0);++i)
  {
    if(dim ==2 )
    {

      // Dxx for un

      xhr= UF->get_DOF_coordinate( shift.i+i,0, 0 ) - UF->get_DOF_coordinate( shift.i+i-1, 0, 0 ) ;
      xright = UF->DOF_value( shift.i+i, j, 0, 0, 0 ) - UF->DOF_value( shift.i+i-1, j, 0, 0, 0 ) ;

      xvalue = xright/xhr;
      
      // Dyy for un
      yhr= UF->get_DOF_coordinate( shift.j+j,1, 1 ) - UF->get_DOF_coordinate( shift.j+j-1, 1, 1 ) ;
      yright = UF->DOF_value( i, shift.j+j, 0, 1, 0 ) - UF->DOF_value( i, shift.j+j-1, 0, 1, 0 ) ;

      yvalue = yright/yhr;

      dx = PF->get_cell_size( i,0, 0 );

      // Assemble the bodyterm
      value = -(rho*(xvalue + yvalue)*dx)/(t_it -> time_step());


    }
    else
    {
        // Dxx for un
        xhr= UF->get_DOF_coordinate( shift.i+i,0, 0 ) - UF->get_DOF_coordinate( shift.i+i-1, 0, 0 ) ;
        xright = UF->DOF_value( shift.i+i, j, k, 0, 0 ) - UF->DOF_value( shift.i+i-1, j, k, 0, 0 ) ;

        xvalue = xright/xhr;

        // Dyy for un
        yhr= UF->get_DOF_coordinate( shift.j+j,1, 1 ) - UF->get_DOF_coordinate( shift.j+j-1, 1, 1 ) ;
        yright = UF->DOF_value( i, shift.j+j, k, 1, 0 ) - UF->DOF_value( i, shift.j+j-1, k, 1, 0 ) ;

        yvalue = yright/yhr;

        double zhr,zright,zvalue=0.;
        zhr= UF->get_DOF_coordinate( shift.k+k,2, 2 ) - UF->get_DOF_coordinate( shift.k+k-1, 2, 2 ) ;
        zright = UF->DOF_value( i, j, shift.k+k, 2, 0 ) - UF->DOF_value( i, j, shift.k+k-1, 2, 0 ) ;

        zvalue = zright/zhr;

        dx = PF->get_cell_size( i,0, 0 );
        value = -(rho*(xvalue + yvalue + zvalue)*dx)/(t_it -> time_step());
    }

    pos = i - min_unknown_index(0);

    // if(j==min_unknown_index(1))
    //   MAC::out()<<xright<<"  "<<yright<<endl;

    if(is_Pperiodic[0] == 0){
      if(rank_in_i[0] == nb_ranks_comm_i[0]-1){
         local_rhs_x->set_item( pos, value);
      }
      else{
         if(i == max_unknown_index(0))
           fe = value;
         else
           local_rhs_x->set_item( pos, value);
      }  
    }
    else{
      if(nb_ranks_comm_i[0] > 1){
         if(i == max_unknown_index(0))
           fe = value;
         else
           local_rhs_x->set_item( pos, value);
      }
      else
        local_rhs_x->set_item( pos, value);
    }
    
  }

 // No term added to rhs for neumann as flux was set to zero during assemble of lhs
 // No term added to rhs for dirichlet as boundary condition for pseudo pressure(psi) becomes zero

  /*m = int(min_unknown_index(0)) - 1;
  if ( PF->DOF_in_domain( m, j, k, 0 ) )
    if ( PF->DOF_has_imposed_Dirichlet_value( m, j, k, 0 ) )
    {
      double ai = 1/(PF->get_DOF_coordinate( m+1,0,0) - PF->get_DOF_coordinate( m,0,0));
      double dirichlet_value = PF->DOF_value( m, j, k, 0, 0 ) ;
      local_rhs_x->add_to_item( 0, + ai * dirichlet_value );
    }

  m = int(max_unknown_index(0)) + 1;
  if ( PF->DOF_in_domain( m, j, k, 0 ) )
    if ( PF->DOF_has_imposed_Dirichlet_value( m, j, k, 0 ) )
    {
      double ai = 1/(PF->get_DOF_coordinate( m,0,0) - PF->get_DOF_coordinate( m-1,0,0));
      double dirichlet_value = PF->DOF_value( m, j, k, 0, 0 ) ;
      local_rhs_x->add_to_item( local_rhs_x->nb_rows()-1 , + ai * dirichlet_value );
    }*/

  return fe;
}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: pressure_solve_x_for_secondorder ( size_t const& j, size_t const& k, FV_TimeIterator const* t_it,double * packed_data)
//---------------------------------------------------------------------------
{
  // Get local min and max indices
  size_t comp=0;
   size_t_vector min_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     min_unknown_index(l) =
      PF->get_min_index_unknown_handled_by_proc( comp, l ) ;
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     max_unknown_index(l) =
      PF->get_max_index_unknown_handled_by_proc( comp, l ) ;

   double fe = pressure_assemble_local_rhs_x(j,k,t_it);

   // create a replica of local rhs vector in local solution vector
   LA_SeqVector* local_rhs_x = GLOBAL_EQ->get_local_temp_x_P() ;
   LA_SeqVector* local_solution_x = GLOBAL_EQ->get_local_solution_temp_x_P() ;
   for(size_t i=0;i<local_rhs_x->nb_rows();i++){
    local_solution_x->set_item(i,local_rhs_x->item(i));
   }

   // Solve for xi locally and put it in local solution vector
   GLOBAL_EQ->DS_NavierStokes_x_local_unknown_solver_P(local_solution_x);

   LA_SeqMatrix* Aei_x = GLOBAL_EQ-> get_aei_P(0);
   LA_SeqVector* Vec_temp_x = GLOBAL_EQ-> get_temp_x_P();

   for(int i=0;i<Vec_temp_x->nb_rows();i++){
          Vec_temp_x->set_item(i,0);
    }
   // Calculate Aei*xi in each proc locally
   Aei_x->multiply_vec_then_add(local_solution_x,Vec_temp_x);

   // Pack the data
   size_t vec_pos;
   if(dim == 2){
      vec_pos=j-min_unknown_index(1);
   }
   else{
      vec_pos=(j-min_unknown_index(1))+(max_unknown_index(1)-min_unknown_index(1)+1)*(k-min_unknown_index(2));
   }/**/
   if(rank_in_i[0] == 0){
      // Check if bc is periodic in x
     // If it is, we need to pack two elements apart from fe

      if(is_Pperiodic[0])
          packed_data[3*vec_pos+0]=Vec_temp_x->item(nb_ranks_comm_i[0]-1);
      else
          packed_data[3*vec_pos+0]=0;

      packed_data[3*vec_pos+1]=Vec_temp_x->item(rank_in_i[0]);
   }
   else if(rank_in_i[0] == nb_ranks_comm_i[0]-1){
      packed_data[3*vec_pos+0]=Vec_temp_x->item(rank_in_i[0]-1);
      
      // Check if bc is periodic in x
      // If it is, we need to pack two elements apart from fe
      if(is_Pperiodic[0])
          packed_data[3*vec_pos+1]=Vec_temp_x->item(rank_in_i[0]);
      else
          packed_data[3*vec_pos+1]=0;
   }
   else{
      packed_data[3*vec_pos+0]=Vec_temp_x->item(rank_in_i[0]-1);
      packed_data[3*vec_pos+1]=Vec_temp_x->item(rank_in_i[0]);
   }
   packed_data[3*vec_pos+2] = fe; // Send the fe values and 0 for last proc
}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: pressure_solve_interface_unknowns_x ( double * packed_data, size_t nb_received_data, FV_TimeIterator const* t_it)
//---------------------------------------------------------------------------
{
   size_t i,j,p;
   size_t k =0;
   size_t comp=0;

   // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     min_unknown_index(l) =
      PF->get_min_index_unknown_handled_by_proc( comp, l ) ;
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     max_unknown_index(l) =
      PF->get_max_index_unknown_handled_by_proc( comp, l ) ;

   // Vector for fe
   LA_SeqVector* interface_rhs_x = GLOBAL_EQ->get_interface_temp_x_P() ;
   // Vector for fi
   LA_SeqVector* local_rhs_x = GLOBAL_EQ->get_local_temp_x_P() ;

   LA_SeqVector* Vec_temp_x = GLOBAL_EQ-> get_temp_x_P();

   LA_SeqMatrix* Aie_x = GLOBAL_EQ-> get_aie_P(0);

   size_t nb_send_data= (2*nb_received_data)/3;


  // Send and receive the data first pass
   //if ( my_rank == is_master )
   if ( rank_in_i[0] == 0 )
   {
      for (i=1;i<nb_ranks_comm_i[0];++i)
      //for (i=1;i<nb_procs;++i)
      {
        // Receive the data
        //pelCOMM->receive( i, all_received_data[i], nb_received_data );
        static MPI_Status status;
        MPI_Recv( all_receive_data_P_x[i], nb_received_data, MPI_DOUBLE, i, 0,
                          DDS_Comm_i[0], &status ) ;
      }

      // Solve system of interface unknowns for each y
      if(dim == 2)
      {
        for(j=min_unknown_index(1);j<=max_unknown_index(1);j++){
          size_t nb_interface_unknowns = Vec_temp_x->nb_rows();
          for(i=0;i<nb_interface_unknowns;i++){
            Vec_temp_x->set_item(i,0);
            interface_rhs_x->set_item(i,0);
          }
          p = j-min_unknown_index(1);
          if(is_Pperiodic[0])
              Vec_temp_x->set_item(nb_ranks_comm_i[0]-1,packed_data[3*p]);
          Vec_temp_x->set_item(0,packed_data[3*p+1]);
          interface_rhs_x->set_item(0,packed_data[3*p+2]);

          // Vec_temp might contain previous values

          for(i=1;i<nb_ranks_comm_i[0];i++){

            if(i!=nb_ranks_comm_i[0]-1){
              Vec_temp_x->add_to_item(i-1,all_receive_data_P_x[i][3*p]);
              Vec_temp_x->add_to_item(i,all_receive_data_P_x[i][3*p+1]);
              interface_rhs_x->set_item(i,all_receive_data_P_x[i][3*p+2]);  // Assemble the interface rhs fe
            }
            else{
              if(is_Pperiodic[0] ==0){
                  Vec_temp_x->add_to_item(i-1,all_receive_data_P_x[i][3*p]);
              }
              else{
                  Vec_temp_x->add_to_item(i-1,all_receive_data_P_x[i][3*p]);
                  // If periodic in x, last proc has an interface unknown
                  Vec_temp_x->add_to_item(i,all_receive_data_P_x[i][3*p+1]);
                  interface_rhs_x->set_item(i,all_receive_data_P_x[i][3*p+2]);
              }
            }
          }

          for(i=0;i<nb_interface_unknowns;i++){
          interface_rhs_x->set_item(i,interface_rhs_x->item(i)-Vec_temp_x->item(i)); // Get fe - Aei*xi to solve for ue
          }

          // Solve for ue (interface unknowns) in the master proc
          GLOBAL_EQ->DS_NavierStokes_x_interface_unknown_solver_P(interface_rhs_x);

          // Pack the interface_rhs_x into the appropriate send_data
          for (i=1;i<nb_ranks_comm_i[0];++i)
          {
            if(i!=nb_ranks_comm_i[0]-1){
              all_send_data_P_x[i][2*p+0]=interface_rhs_x->item(i-1);
              all_send_data_P_x[i][2*p+1]=interface_rhs_x->item(i);
            }
            else{
              all_send_data_P_x[i][2*p+0]=interface_rhs_x->item(i-1);
              if(is_Pperiodic[0])
                all_send_data_P_x[i][2*p+1]=interface_rhs_x->item(i);  
              else
                all_send_data_P_x[i][2*p+1]=0;
            }

          }

          // Need to have the original rhs function assembled for corrosponding j,k pair
          double fe = pressure_assemble_local_rhs_x(j,k,t_it);

          // Setup RHS = fi - Aie*xe for solving ui
          Aie_x->multiply_vec_then_add(interface_rhs_x,local_rhs_x,-1.0,1.0);

          // Solve ui and transfer solution into distributed vector
          GLOBAL_EQ->DS_NavierStokes_x_solver_P(j,k,min_unknown_index(0),local_rhs_x,interface_rhs_x);

        }
      }
      else
      {
        for(k=min_unknown_index(2);k<=max_unknown_index(2);k++)
        {
          for(j=min_unknown_index(1);j<=max_unknown_index(1);j++)
          {
            size_t nb_interface_unknowns = Vec_temp_x->nb_rows();
            for(i=0;i<nb_interface_unknowns;i++){
              Vec_temp_x->set_item(i,0);
              interface_rhs_x->set_item(i,0);
            }
            p = (j-min_unknown_index(1))+(max_unknown_index(1)-min_unknown_index(1)+1)*(k-min_unknown_index(2));
            
            if(is_Pperiodic[0])
                Vec_temp_x->set_item(nb_ranks_comm_i[0]-1,packed_data[3*p]);
            Vec_temp_x->set_item(0,packed_data[3*p+1]);
            interface_rhs_x->set_item(0,packed_data[3*p+2]);

            // Vec_temp might contain previous values


            for(i=1;i<nb_ranks_comm_i[0];i++){

              if(i!=nb_ranks_comm_i[0]-1){
                Vec_temp_x->add_to_item(i-1,all_receive_data_P_x[i][3*p]);
                Vec_temp_x->add_to_item(i,all_receive_data_P_x[i][3*p+1]);
                interface_rhs_x->set_item(i,all_receive_data_P_x[i][3*p+2]);  // Assemble the interface rhs fe
              }
              else{
                if(is_Pperiodic[0] ==0){
                    Vec_temp_x->add_to_item(i-1,all_receive_data_P_x[i][3*p]);
                }
                else{
                    Vec_temp_x->add_to_item(i-1,all_receive_data_P_x[i][3*p]);
                    // If periodic in x, last proc has an interface unknown
                    Vec_temp_x->add_to_item(i,all_receive_data_P_x[i][3*p+1]);
                    interface_rhs_x->set_item(i,all_receive_data_P_x[i][3*p+2]);
                }
              }
            }

            for(i=0;i<nb_interface_unknowns;i++){
            interface_rhs_x->set_item(i,interface_rhs_x->item(i)-Vec_temp_x->item(i)); // Get fe - Aei*xi to solve for ue
            }


            // Solve for ue (interface unknowns) in the master proc
            GLOBAL_EQ->DS_NavierStokes_x_interface_unknown_solver_P(interface_rhs_x);

            // Pack the interface_rhs_x into the appropriate send_data
            for (i=1;i<nb_ranks_comm_i[0];++i)
            {
              if(i!=nb_ranks_comm_i[0]-1){
                all_send_data_P_x[i][2*p+0]=interface_rhs_x->item(i-1);
                all_send_data_P_x[i][2*p+1]=interface_rhs_x->item(i);
              }
              else{
                all_send_data_P_x[i][2*p+0]=interface_rhs_x->item(i-1);
                if(is_Pperiodic[0])
                  all_send_data_P_x[i][2*p+1]=interface_rhs_x->item(i);  
                else
                  all_send_data_P_x[i][2*p+1]=0;
              }

            }

          // Need to have the original rhs function assembled for corrosponding j,k pair
          double fe = pressure_assemble_local_rhs_x(j,k,t_it);

          // Setup RHS = fi - Aie*xe for solving ui
          Aie_x->multiply_vec_then_add(interface_rhs_x,local_rhs_x,-1.0,1.0);

          // Solve ui and transfer solution into distributed vector
          GLOBAL_EQ->DS_NavierStokes_x_solver_P(j,k,min_unknown_index(0),local_rhs_x,interface_rhs_x);
          }
        }
      }


    }
   else
   {
      // Send the packed data to master
      //pelCOMM->send( is_master, packed_data, nb_received_data );
      MPI_Send( packed_data, nb_received_data, MPI_DOUBLE, 0, 0, DDS_Comm_i[0] ) ;
   }

   // Send the data from master
   if ( rank_in_i[0] == 0 )
   {
     for (i=1;i<nb_ranks_comm_i[0];++i)
      {
        //pelCOMM->send( i, all_send_data[i], nb_send_data );
        MPI_Send( all_send_data_P_x[i], nb_send_data, MPI_DOUBLE, i, 0, DDS_Comm_i[0] ) ;
      }
   }
   else
   {
      // Receive the data

      //pelCOMM->receive( is_master, received_data, nb_received_data );
      static MPI_Status status ;
      MPI_Recv( packed_data, nb_received_data, MPI_DOUBLE, 0, 0,
                          DDS_Comm_i[0], &status ) ;

     // Solve the system of equations in each proc

     if(dim == 2)
     {
        for(j = min_unknown_index(1);j<=max_unknown_index(1);j++)
        {

          p = j-min_unknown_index(1);
          if(rank_in_i[0] != nb_ranks_comm_i[0]-1){
            interface_rhs_x->set_item(rank_in_i[0]-1,packed_data[2*p]);
            interface_rhs_x->set_item(rank_in_i[0],packed_data[2*p+1]);
          }
          else{
            if(is_Pperiodic[0] ==0){
              interface_rhs_x->set_item(rank_in_i[0]-1,packed_data[2*p]);
            }
            else{
              interface_rhs_x->set_item(rank_in_i[0]-1,packed_data[2*p]);
              interface_rhs_x->set_item(rank_in_i[0],packed_data[2*p+1]);
            }

          }
          // Need to have the original rhs function assembled for corrosponding j,k pair
          double fe = pressure_assemble_local_rhs_x(j,k,t_it);

          // Setup RHS = fi - Aie*xe for solving ui
          Aie_x->multiply_vec_then_add(interface_rhs_x,local_rhs_x,-1.0,1.0);

          // Solve ui and transfer solution into distributed vector
          GLOBAL_EQ->DS_NavierStokes_x_solver_P(j,k,min_unknown_index(0),local_rhs_x,interface_rhs_x);
        }
     }
     else
     {
        for(k = min_unknown_index(2);k<=max_unknown_index(2);k++)
        {
          for(j = min_unknown_index(1);j<=max_unknown_index(1);j++)
          {
            p = (j-min_unknown_index(1))+(max_unknown_index(1)-min_unknown_index(1)+1)*(k-min_unknown_index(2));
            if(rank_in_i[0] != nb_ranks_comm_i[0]-1){
              interface_rhs_x->set_item(rank_in_i[0]-1,packed_data[2*p]);
              interface_rhs_x->set_item(rank_in_i[0],packed_data[2*p+1]);
            }
            else{
              if(is_Pperiodic[0] ==0){
                interface_rhs_x->set_item(rank_in_i[0]-1,packed_data[2*p]);
              }
              else{
                interface_rhs_x->set_item(rank_in_i[0]-1,packed_data[2*p]);
                interface_rhs_x->set_item(rank_in_i[0],packed_data[2*p+1]);
              }
            }
            // Need to have the original rhs function assembled for corrosponding j,k pair
            double fe = pressure_assemble_local_rhs_x(j,k,t_it);

            // Setup RHS = fi - Aie*xe for solving ui
            Aie_x->multiply_vec_then_add(interface_rhs_x,local_rhs_x,-1.0,1.0);

            // Solve ui and transfer solution into distributed vector
            GLOBAL_EQ->DS_NavierStokes_x_solver_P(j,k,min_unknown_index(0),local_rhs_x,interface_rhs_x);
          }

        }
     }

   }
}




//---------------------------------------------------------------------------
double
DDS_NavierStokes:: pressure_assemble_local_rhs_y ( size_t const& i, size_t const& k, FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
  // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     min_unknown_index(l) =
      PF->get_min_index_unknown_handled_by_proc( 0, l ) ;
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     max_unknown_index(l) =
      PF->get_max_index_unknown_handled_by_proc( 0, l ) ;

size_t j,pos,m;

double fe=0.,dy;

 // Vector for fi
LA_SeqVector* local_rhs_y = GLOBAL_EQ->get_local_temp_y_P() ;

 // Compute local_rhs_y = psi_{n} for pressure in y
 for (j=min_unknown_index(1);j<=max_unknown_index(1);++j)
 {
   
   pos = j - min_unknown_index(1);
   dy = PF->get_cell_size( j,0,1);

   if(is_Pperiodic[1] == 0){
      if(rank_in_i[1] == nb_ranks_comm_i[1]-1){
       local_rhs_y->set_item( pos, PF->DOF_value( i, j, k, 0, 1 )*dy);
      }
      else{
         if(j == max_unknown_index(1))
           fe = PF->DOF_value( i, j, k, 0, 1 )*dy;
         else
           local_rhs_y->set_item( pos, PF->DOF_value( i, j, k, 0, 1 )*dy);
      }
   }
   else{
      if(nb_ranks_comm_i[1] > 1){
         if(j == max_unknown_index(1))
           fe = PF->DOF_value( i, j, k, 0, 1 )*dy;
         else
           local_rhs_y->set_item( pos, PF->DOF_value( i, j, k, 0, 1 )*dy);
      }
      else
        local_rhs_y->set_item( pos, PF->DOF_value( i, j, k, 0, 1 )*dy);
   }
    
 }

 // No term added to rhs for neumann as flux was set to zero during assemble of lhs
 // No term added to rhs for dirichlet as boundary condition for pseudo pressure(psi) becomes zero

 /*m = int(min_unknown_index(1)) - 1;
 if ( PF->DOF_in_domain( i, m, k, 0 ) )
   if ( PF->DOF_has_imposed_Dirichlet_value( i, m, k, 0 ) )
   {
     double ai = 1/(PF->get_DOF_coordinate( m+1,0,1) - PF->get_DOF_coordinate( m,0,1));
     double dirichlet_value = PF->DOF_value( i, m, k, 0, 0 ) ;
     local_rhs_y->add_to_item( 0, + ai * dirichlet_value );
   }

 m = int(max_unknown_index(1)) + 1;
 if ( PF->DOF_in_domain( i, m, k, 0 ) )
   if ( PF->DOF_has_imposed_Dirichlet_value( i, m, k, 0 ) )
   {
     double ai = 1/(PF->get_DOF_coordinate( m,0,1) - PF->get_DOF_coordinate( m-1,0,1));
     double dirichlet_value = PF->DOF_value( i, m, k, 0, 0 ) ;
     local_rhs_y->add_to_item( local_rhs_y->nb_rows()-1, + ai * dirichlet_value );
   }*/

  return fe;
}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: pressure_solve_y_for_secondorder ( size_t const& i, size_t const& k, FV_TimeIterator const* t_it,double * packed_data )
//---------------------------------------------------------------------------
{
  // Get local min and max indices
  size_t comp=0;
   size_t_vector min_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     min_unknown_index(l) =
      PF->get_min_index_unknown_handled_by_proc( comp, l ) ;
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     max_unknown_index(l) =
      PF->get_max_index_unknown_handled_by_proc( comp, l ) ;

   double fe = pressure_assemble_local_rhs_y(i,k,t_it);

   size_t j;
   // create a replica of local rhs vector in local solution vector
   LA_SeqVector* local_rhs_y = GLOBAL_EQ->get_local_temp_y_P() ;

   LA_SeqVector* local_solution_y = GLOBAL_EQ->get_local_solution_temp_y_P() ;
   for(j=0;j<local_rhs_y->nb_rows();j++){
    local_solution_y->set_item(j,local_rhs_y->item(j));
   }

   // Solve for xi locally and put it in local solution vector
   GLOBAL_EQ->DS_NavierStokes_y_local_unknown_solver_P(local_solution_y);

   LA_SeqMatrix* Aei = GLOBAL_EQ-> get_aei_P(1);
   LA_SeqVector* Vec_temp = GLOBAL_EQ-> get_temp_y_P();

   for(j=0;j<Vec_temp->nb_rows();j++){
          Vec_temp->set_item(j,0);
    }
   // Calculate Aei*xi in each proc locally
   Aei->multiply_vec_then_add(local_solution_y,Vec_temp);
   size_t vec_pos;
   if(dim==2){
      vec_pos=i-min_unknown_index(0);
   }
   else{
      vec_pos=(k-min_unknown_index(2))+(max_unknown_index(2)-min_unknown_index(2)+1)*(i-min_unknown_index(0));
   }
   if(rank_in_i[1] == 0){
     // Check if bc is periodic in x
     // If it is, we need to pack two elements apart from fe

      if(is_Pperiodic[1])
          packed_data[3*vec_pos+0]=Vec_temp->item(nb_ranks_comm_i[1]-1);
      else
          packed_data[3*vec_pos+0]=0;
        
     packed_data[3*vec_pos+1]=Vec_temp->item(rank_in_i[1]);
   }
   else if(rank_in_i[1] == nb_ranks_comm_i[1]-1){
     packed_data[3*vec_pos+0]=Vec_temp->item(rank_in_i[1]-1);
     
     // Check if bc is periodic in x
      // If it is, we need to pack two elements apart from fe
      if(is_Pperiodic[1])
          packed_data[3*vec_pos+1]=Vec_temp->item(rank_in_i[1]);
      else
          packed_data[3*vec_pos+1]=0;
   }
   else{
     packed_data[3*vec_pos+0]=Vec_temp->item(rank_in_i[1]-1);
     packed_data[3*vec_pos+1]=Vec_temp->item(rank_in_i[1]);
   }
   packed_data[3*vec_pos+2] = fe; // Send the fe values and 0 for last proc

}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: pressure_solve_interface_unknowns_y ( double * packed_data, size_t nb_received_data, FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   size_t i,j,p;
   size_t k =0;
   size_t comp=0;

   // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     min_unknown_index(l) =
      PF->get_min_index_unknown_handled_by_proc( comp, l ) ;
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     max_unknown_index(l) =
      PF->get_max_index_unknown_handled_by_proc( comp, l ) ;

   // Vector for fe
   LA_SeqVector* interface_rhs_y = GLOBAL_EQ->get_interface_temp_y_P() ;
   // Vector for fi
   LA_SeqVector* local_rhs_y = GLOBAL_EQ->get_local_temp_y_P();

   LA_SeqVector* Vec_temp = GLOBAL_EQ-> get_temp_y_P();

   LA_SeqMatrix* Aie = GLOBAL_EQ-> get_aie_P(1);

   size_t nb_send_data = (2*nb_received_data)/3;

  // Send and receive the data first pass
   if ( rank_in_i[1] == 0 )
   {
      for (i=1;i<nb_ranks_comm_i[1];++i)
      {
        // Receive the data
        //pelCOMM->receive( i, all_received_data[i], nb_received_data );
        static MPI_Status status ;
        MPI_Recv( all_receive_data_P_y[i], nb_received_data, MPI_DOUBLE, i, 0,
                          DDS_Comm_i[1], &status ) ;
      }

      // Solve system of interface unknowns for each x
      if(dim == 2)
      {
        for(i=min_unknown_index(0);i<=max_unknown_index(0);i++)
        {
          size_t nb_interface_unknowns = Vec_temp->nb_rows();
          for(j=0;j<nb_interface_unknowns;j++){
            Vec_temp->set_item(j,0);
            interface_rhs_y->set_item(j,0);
          }
          p = i-min_unknown_index(0);
          if(is_Pperiodic[1])
              Vec_temp->set_item(nb_ranks_comm_i[1]-1,packed_data[3*p]);
          Vec_temp->set_item(0,packed_data[3*p+1]);
          interface_rhs_y->set_item(0,packed_data[3*p+2]);

          // Vec_temp might contain previous values

          for(j=1;j<nb_ranks_comm_i[1];j++){

            if(j!=nb_ranks_comm_i[1]-1){
              Vec_temp->add_to_item(j-1,all_receive_data_P_y[j][3*p]);
              Vec_temp->add_to_item(j,all_receive_data_P_y[j][3*p+1]);
              interface_rhs_y->set_item(j,all_receive_data_P_y[j][3*p+2]);  // Assemble the interface rhs fe
            }
            else{
              if(is_Pperiodic[1] ==0){
                  Vec_temp->add_to_item(j-1,all_receive_data_P_y[j][3*p]);
              }
              else{
                  // If periodic in y, last proc has an interface unknown
                  Vec_temp->add_to_item(j-1,all_receive_data_P_y[j][3*p]);
                  Vec_temp->add_to_item(j,all_receive_data_P_y[j][3*p+1]);
                  interface_rhs_y->set_item(j,all_receive_data_P_y[j][3*p+2]);
              }
            }
          }

          for(j=0;j<nb_interface_unknowns;j++){
          interface_rhs_y->set_item(j,interface_rhs_y->item(j)-Vec_temp->item(j)); // Get fe - Aei*xi to solve for ue
          }

          // Solve for ue (interface unknowns) in the master proc
          GLOBAL_EQ->DS_NavierStokes_y_interface_unknown_solver_P(interface_rhs_y);

          //if(p==0 && rank_in_i[0] == 1)
          //  interface_rhs_y->print_items(MAC::out(),0);

          // Pack the interface_rhs_x into the appropriate send_data
          for (j=1;j<nb_ranks_comm_i[1];++j)
          {
            if(j!=nb_ranks_comm_i[1]-1){
              all_send_data_P_y[j][2*p+0]=interface_rhs_y->item(j-1);
              all_send_data_P_y[j][2*p+1]=interface_rhs_y->item(j);
            }
            else{
              all_send_data_P_y[j][2*p+0]=interface_rhs_y->item(j-1);
              if(is_Pperiodic[1])
                all_send_data_P_y[j][2*p+1]=interface_rhs_y->item(j); 
              else
                all_send_data_P_y[j][2*p+1]=0;
            }

          }

          // Need to have the original rhs function assembled for corrosponding j,k pair
          double fe = pressure_assemble_local_rhs_y(i,k,t_it);

          // Setup RHS = fi - Aie*xe for solving ui
          Aie->multiply_vec_then_add(interface_rhs_y,local_rhs_y,-1.0,1.0);

          // Solve ui and transfer solution into distributed vector
          GLOBAL_EQ->DS_NavierStokes_y_solver_P(i,k,min_unknown_index(1),local_rhs_y,interface_rhs_y);

        }
      }
      else{
        for(i=min_unknown_index(0);i<=max_unknown_index(0);i++)
        {
          for(k=min_unknown_index(2);k<=max_unknown_index(2);k++)
          {
            size_t nb_interface_unknowns = Vec_temp->nb_rows();
            for(j=0;j<nb_interface_unknowns;j++){
              Vec_temp->set_item(j,0);
              interface_rhs_y->set_item(j,0);
            }
            p = (k-min_unknown_index(2))+(max_unknown_index(2)-min_unknown_index(2)+1)*(i-min_unknown_index(0));
            if(is_Pperiodic[1])
                Vec_temp->set_item(nb_ranks_comm_i[1]-1,packed_data[3*p]);
            Vec_temp->set_item(0,packed_data[3*p+1]);
            interface_rhs_y->set_item(0,packed_data[3*p+2]);

            // Vec_temp might contain previous values

            for(j=1;j<nb_ranks_comm_i[1];j++){

              if(j!=nb_ranks_comm_i[1]-1){
                Vec_temp->add_to_item(j-1,all_receive_data_P_y[j][3*p]);
                Vec_temp->add_to_item(j,all_receive_data_P_y[j][3*p+1]);
                interface_rhs_y->set_item(j,all_receive_data_P_y[j][3*p+2]);  // Assemble the interface rhs fe
              }
              else{
                if(is_Pperiodic[1] == 0){
                  Vec_temp->add_to_item(j-1,all_receive_data_P_y[j][3*p]);
                }
                else{
                  Vec_temp->add_to_item(j-1,all_receive_data_P_y[j][3*p]);
                  Vec_temp->add_to_item(j,all_receive_data_P_y[j][3*p+1]);
                  interface_rhs_y->set_item(j,all_receive_data_P_y[j][3*p+2]);
                }
              }
            }

            for(j=0;j<nb_interface_unknowns;j++){
            interface_rhs_y->set_item(j,interface_rhs_y->item(j)-Vec_temp->item(j)); // Get fe - Aei*xi to solve for ue
            }

            // Solve for ue (interface unknowns) in the master proc
            GLOBAL_EQ->DS_NavierStokes_y_interface_unknown_solver_P(interface_rhs_y);

            //if(p==0 && rank_in_i[0] == 1)
            //  interface_rhs_y->print_items(MAC::out(),0);

            // Pack the interface_rhs_x into the appropriate send_data
            for (j=1;j<nb_ranks_comm_i[1];++j)
            {
              if(j!=nb_ranks_comm_i[1]-1){
                all_send_data_P_y[j][2*p+0]=interface_rhs_y->item(j-1);
                all_send_data_P_y[j][2*p+1]=interface_rhs_y->item(j);
              }
              else{
                all_send_data_P_y[j][2*p+0]=interface_rhs_y->item(j-1);
                if(is_Pperiodic[1])
                  all_send_data_P_y[j][2*p+1]=interface_rhs_y->item(j);
                else
                  all_send_data_P_y[j][2*p+1]=0;
              }

            }

            // Need to have the original rhs function assembled for corrosponding j,k pair
            double fe = pressure_assemble_local_rhs_y(i,k,t_it);

            // Setup RHS = fi - Aie*xe for solving ui
            Aie->multiply_vec_then_add(interface_rhs_y,local_rhs_y,-1.0,1.0);

            // Solve ui and transfer solution into distributed vector
            GLOBAL_EQ->DS_NavierStokes_y_solver_P(i,k,min_unknown_index(1),local_rhs_y,interface_rhs_y);
          }
        }
      }

    }
   else
   {
      // Send the packed data to master
      //pelCOMM->send( is_master, packed_data, nb_received_data );

      MPI_Send( packed_data, nb_received_data, MPI_DOUBLE, 0, 0, DDS_Comm_i[1] ) ;
   }
   // Send the data from master
   if ( rank_in_i[1] == 0 )
   {
     for (i=1;i<nb_ranks_comm_i[1];++i)
      {
        //pelCOMM->send( i, all_send_data[i], nb_send_data );
        MPI_Send( all_send_data_P_y[i], nb_send_data, MPI_DOUBLE, i, 0, DDS_Comm_i[1] ) ;
      }
   }
   else
   {
      // Receive the data
      //pelCOMM->receive( is_master, received_data, nb_received_data );
      static MPI_Status status ;
      MPI_Recv( packed_data, nb_received_data, MPI_DOUBLE, 0, 0,
                          DDS_Comm_i[1], &status ) ;

     // Solve the system of equations in each proc

     if(dim == 2)
     {
        for(i = min_unknown_index(0);i<=max_unknown_index(0);i++)
        {
          p = i-min_unknown_index(0);
          if(rank_in_i[1] != nb_ranks_comm_i[1]-1){
            interface_rhs_y->set_item(rank_in_i[1]-1,packed_data[2*p]);
            interface_rhs_y->set_item(rank_in_i[1],packed_data[2*p+1]);
          }
          else{
            if(is_Pperiodic[1] ==0){
              interface_rhs_y->set_item(rank_in_i[1]-1,packed_data[2*p]);
            }
            else{
              interface_rhs_y->set_item(rank_in_i[1]-1,packed_data[2*p]);
              interface_rhs_y->set_item(rank_in_i[1],packed_data[2*p+1]);
            }
          }
          // Need to have the original rhs function assembled for corrosponding j,k pair
          double fe = pressure_assemble_local_rhs_y(i,k,t_it);

          // Setup RHS = fi - Aie*xe for solving ui
          Aie->multiply_vec_then_add(interface_rhs_y,local_rhs_y,-1.0,1.0);

          // Solve ui and transfer solution into distributed vector
          GLOBAL_EQ->DS_NavierStokes_y_solver_P(i,k,min_unknown_index(1),local_rhs_y,interface_rhs_y);
        }
     }
     else
     {
        for(i = min_unknown_index(0);i<=max_unknown_index(0);i++)
        {
          for(k = min_unknown_index(2);k<=max_unknown_index(2);k++)
          {
            p = (k-min_unknown_index(2))+(max_unknown_index(2)-min_unknown_index(2)+1)*(i-min_unknown_index(0));
            if(rank_in_i[1] != nb_ranks_comm_i[1]-1){
              interface_rhs_y->set_item(rank_in_i[1]-1,packed_data[2*p]);
              interface_rhs_y->set_item(rank_in_i[1],packed_data[2*p+1]);
            }
            else{
              if(is_Pperiodic[1] ==0){
                interface_rhs_y->set_item(rank_in_i[1]-1,packed_data[2*p]);
              }
              else{
                interface_rhs_y->set_item(rank_in_i[1]-1,packed_data[2*p]);
                interface_rhs_y->set_item(rank_in_i[1],packed_data[2*p+1]);
              }
            }
            // Need to have the original rhs function assembled for corrosponding j,k pair
            double fe = pressure_assemble_local_rhs_y(i,k,t_it);

            // Setup RHS = fi - Aie*xe for solving ui
            Aie->multiply_vec_then_add(interface_rhs_y,local_rhs_y,-1.0,1.0);

            // Solve ui and transfer solution into distributed vector
            GLOBAL_EQ->DS_NavierStokes_y_solver_P(i,k,min_unknown_index(1),local_rhs_y,interface_rhs_y);
          }
        }
     }

   }
}




//---------------------------------------------------------------------------
double
DDS_NavierStokes:: pressure_assemble_local_rhs_z ( size_t const& i, size_t const& j, FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
  // Get local min and max indices

 size_t_vector min_unknown_index(dim,0);
 for (size_t l=0;l<dim;++l)
   min_unknown_index(l) =
    PF->get_min_index_unknown_handled_by_proc( 0, l ) ;
 size_t_vector max_unknown_index(dim,0);
 for (size_t l=0;l<dim;++l)
   max_unknown_index(l) =
    PF->get_max_index_unknown_handled_by_proc( 0, l ) ;

size_t k,pos,m;

double fe=0.,dz;

// Vector for fi
LA_SeqVector* local_rhs_z = GLOBAL_EQ->get_local_temp_z_P();

// Compute local_rhs_z = phi_{n} for pressure in z
 for (k=min_unknown_index(2);k<=max_unknown_index(2);++k)
 {
   
   pos = k - min_unknown_index(2);
   dz = PF->get_cell_size( k,0,2);

   if(is_Pperiodic[2] == 0){
     if(rank_in_i[2] == nb_ranks_comm_i[2]-1){
        local_rhs_z->set_item( pos, PF->DOF_value( i, j, k, 0, 1 )*dz);
     }
     else{
       if(k == max_unknown_index(2))
         fe = PF->DOF_value( i, j, k, 0, 1 )*dz ;
       else
         local_rhs_z->set_item( pos, PF->DOF_value( i, j, k, 0, 1 )*dz);
     }
   }
   else{
      if(nb_ranks_comm_i[2] > 1){
       if(k == max_unknown_index(2))
         fe = PF->DOF_value( i, j, k, 0, 1 )*dz ;
       else
         local_rhs_z->set_item( pos, PF->DOF_value( i, j, k, 0, 1 )*dz);
      }
      else
        local_rhs_z->set_item( pos, PF->DOF_value( i, j, k, 0, 1 )*dz);
   }
   
 }

 // No term added to rhs for neumann as flux was set to zero during assemble of lhs
 // No term added to rhs for dirichlet as boundary condition for pseudo pressure(psi) becomes zero

 /*m = int(min_unknown_index(2)) - 1;
 if ( PF->DOF_in_domain( i, j, m, 0 ) )
   if ( PF->DOF_has_imposed_Dirichlet_value( i, j, m, 0 ) )
   {
     double ai = 1/(PF->get_DOF_coordinate( m+1,0,2) - PF->get_DOF_coordinate( m,0,2));
     double dirichlet_value = PF->DOF_value( i, j, m, 0, 0 ) ;
     local_rhs_z->add_to_item( 0, + ai * dirichlet_value );
   }

 m = int(max_unknown_index(2)) + 1;
 if ( PF->DOF_in_domain( i, j, m, 0 ) )
   if ( PF->DOF_has_imposed_Dirichlet_value( i, j, m, 0 ) )
   {
     double ai = 1/(PF->get_DOF_coordinate( m,0,2) - PF->get_DOF_coordinate( m-1,0,2));
     double dirichlet_value = PF->DOF_value( i, j, m, 0, 0
      ) ;
     local_rhs_z->add_to_item( local_rhs_z->nb_rows()-1, + ai * dirichlet_value );
   }*/

return fe;
}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: pressure_solve_z_for_secondorder ( size_t const& i, size_t const& j,FV_TimeIterator const* t_it,double * packed_data )
//---------------------------------------------------------------------------
{
  size_t comp=0;
  // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     min_unknown_index(l) =
      PF->get_min_index_unknown_handled_by_proc( comp, l ) ;
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     max_unknown_index(l) =
      PF->get_max_index_unknown_handled_by_proc( comp, l ) ;

   double fe = pressure_assemble_local_rhs_z(i,j,t_it);

   size_t k;
   // create a replica of local rhs vector in local solution vector
   LA_SeqVector* local_rhs_z = GLOBAL_EQ->get_local_temp_z_P() ;

   LA_SeqVector* local_solution_z = GLOBAL_EQ->get_local_solution_temp_z_P() ;
   for(k=0;k<local_rhs_z->nb_rows();k++){
    local_solution_z->set_item(k,local_rhs_z->item(k));
   }

   // Solve for xi locally and put it in local solution vector
   GLOBAL_EQ->DS_NavierStokes_z_local_unknown_solver_P(local_solution_z);

   LA_SeqMatrix* Aei = GLOBAL_EQ-> get_aei_P(2);
   LA_SeqVector* Vec_temp = GLOBAL_EQ-> get_temp_z_P();

   for(k=0;k<Vec_temp->nb_rows();k++){
          Vec_temp->set_item(k,0);
    }
   // Calculate Aei*xi in each proc locally
   Aei->multiply_vec_then_add(local_solution_z,Vec_temp);


    size_t vec_pos=(i-min_unknown_index(0))+(max_unknown_index(0)-min_unknown_index(0)+1)*(j-min_unknown_index(1));

    if(rank_in_i[2] == 0){
      // Check if bc is periodic in x
     // If it is, we need to pack two elements apart from fe

      if(is_Pperiodic[2])
          packed_data[3*vec_pos+0]=Vec_temp->item(nb_ranks_comm_i[2]-1);
      else
          packed_data[3*vec_pos+0]=0;
        
      packed_data[3*vec_pos+1]=Vec_temp->item(rank_in_i[2]);
    }
    else if(rank_in_i[2] == nb_ranks_comm_i[2]-1){
      packed_data[3*vec_pos+0]=Vec_temp->item(rank_in_i[2]-1);
      
      // Check if bc is periodic in x
      // If it is, we need to pack two elements apart from fe
      if(is_Pperiodic[2])
          packed_data[3*vec_pos+1]=Vec_temp->item(rank_in_i[2]);
      else
          packed_data[3*vec_pos+1]=0;

    }
    else{
      packed_data[3*vec_pos+0]=Vec_temp->item(rank_in_i[2]-1);
      packed_data[3*vec_pos+1]=Vec_temp->item(rank_in_i[2]);
    }
    packed_data[3*vec_pos+2] = fe; // Send the fe values and 0 for last proc
}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: pressure_solve_interface_unknowns_z ( double * packed_data, size_t nb_received_data, FV_TimeIterator const* t_it  )
//---------------------------------------------------------------------------
{
   size_t i,j,p,m;
   size_t comp=0;

  // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     min_unknown_index(l) =
      PF->get_min_index_unknown_handled_by_proc( comp, l ) ;
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l)
     max_unknown_index(l) =
      PF->get_max_index_unknown_handled_by_proc( comp, l ) ;

   // Vector for fe
   LA_SeqVector* interface_rhs_z = GLOBAL_EQ->get_interface_temp_z_P() ;
   // Vector for fi
   LA_SeqVector* local_rhs_z = GLOBAL_EQ->get_local_temp_z_P();

   LA_SeqVector* Vec_temp = GLOBAL_EQ-> get_temp_z_P();

   LA_SeqMatrix* Aie = GLOBAL_EQ-> get_aie_P(2);

   size_t nb_send_data = (2*nb_received_data)/3;

  // Send and receive the data first pass
   if ( rank_in_i[2] == 0 )
   {
      for (p=1;p<nb_ranks_comm_i[2];++p)
      {
        // Receive the data
        //pelCOMM->receive( i, all_received_data[i], nb_received_data );
        static MPI_Status status ;
        MPI_Recv( all_receive_data_P_z[p], nb_received_data, MPI_DOUBLE, p, 0,
                          DDS_Comm_i[2], &status ) ;
      }

      // Solve system of interface unknowns for each x

      for(j=min_unknown_index(1);j<=max_unknown_index(1);j++)
      {
        for(i=min_unknown_index(0);i<=max_unknown_index(0);i++)
        {
          size_t nb_interface_unknowns = Vec_temp->nb_rows();
          for(m=0;m<nb_interface_unknowns;m++){
            Vec_temp->set_item(m,0);
            interface_rhs_z->set_item(m,0);
          }
          p = (i-min_unknown_index(0))+(max_unknown_index(0)-min_unknown_index(0)+1)*(j-min_unknown_index(1));

          if(is_Pperiodic[2])
              Vec_temp->set_item(nb_ranks_comm_i[2]-1,packed_data[3*p]);
          Vec_temp->set_item(0,packed_data[3*p+1]);
          interface_rhs_z->set_item(0,packed_data[3*p+2]);

          // Vec_temp might contain previous values

          for(m=1;m<nb_ranks_comm_i[2];m++){

            if(m!=nb_ranks_comm_i[2]-1){
              Vec_temp->add_to_item(m-1,all_receive_data_P_z[m][3*p]);
              Vec_temp->add_to_item(m,all_receive_data_P_z[m][3*p+1]);
              interface_rhs_z->set_item(m,all_receive_data_P_z[m][3*p+2]);  // Assemble the interface rhs fe
            }
            else{
              if(is_Pperiodic[2] == 0)
                  Vec_temp->add_to_item(m-1,all_receive_data_P_z[m][3*p]);
              else{
                  Vec_temp->add_to_item(m-1,all_receive_data_P_z[m][3*p]);
                  Vec_temp->add_to_item(m,all_receive_data_P_z[m][3*p+1]);
                  interface_rhs_z->set_item(m,all_receive_data_P_z[m][3*p+2]);
              }  
            }
          }

          for(m=0;m<nb_interface_unknowns;m++){
          interface_rhs_z->set_item(m,interface_rhs_z->item(m)-Vec_temp->item(m)); // Get fe - Aei*xi to solve for ue
          }

          // Solve for ue (interface unknowns) in the master proc
          GLOBAL_EQ->DS_NavierStokes_z_interface_unknown_solver_P(interface_rhs_z);

          // Pack the interface_rhs_x into the appropriate send_data
          for (m=1;m<nb_ranks_comm_i[2];++m)
          {
            if(m!=nb_ranks_comm_i[2]-1){
              all_send_data_P_z[m][2*p+0]=interface_rhs_z->item(m-1);
              all_send_data_P_z[m][2*p+1]=interface_rhs_z->item(m);
            }
            else{
              all_send_data_P_z[m][2*p+0]=interface_rhs_z->item(m-1);
              if(is_Pperiodic[2])
                all_send_data_P_z[m][2*p+1]=interface_rhs_z->item(m);
              else
                all_send_data_P_z[m][2*p+1]=0;
            }

          }

          // Need to have the original rhs function assembled for corrosponding j,k pair
          double fe = pressure_assemble_local_rhs_z(i,j,t_it);

          // Setup RHS = fi - Aie*xe for solving ui
          Aie->multiply_vec_then_add(interface_rhs_z,local_rhs_z,-1.0,1.0);

          // Solve ui and transfer solution into distributed vector
          GLOBAL_EQ->DS_NavierStokes_z_solver_P(i,j,min_unknown_index(2),local_rhs_z,interface_rhs_z);
        }
      }

    }
   else
   {
      // Send the packed data to master
      //pelCOMM->send( is_master, packed_data, nb_received_data );
      MPI_Send( packed_data, nb_received_data, MPI_DOUBLE, 0, 0, DDS_Comm_i[2] ) ;
   }
   // Send the data from master
   if ( rank_in_i[2] == 0 )
   {
     for (m=1;m<nb_ranks_comm_i[2];++m)
      {
        //pelCOMM->send( i, all_send_data[i], nb_send_data );
        MPI_Send( all_send_data_P_z[m], nb_send_data, MPI_DOUBLE, m, 0, DDS_Comm_i[2] ) ;
      }
   }
   else
   {
      // Receive the data
      //pelCOMM->receive( is_master, received_data, nb_received_data );
      static MPI_Status status ;
      MPI_Recv( packed_data, nb_received_data, MPI_DOUBLE, 0, 0,
                          DDS_Comm_i[2], &status ) ;

     // Solve the system of equations in each proc
     for(j = min_unknown_index(1);j<=max_unknown_index(1);j++){
        for(i = min_unknown_index(0);i<=max_unknown_index(0);i++){
          p = (i-min_unknown_index(0))+(max_unknown_index(0)-min_unknown_index(0)+1)*(j-min_unknown_index(1));
          if(rank_in_i[2] != nb_ranks_comm_i[2]-1){
            interface_rhs_z->set_item(rank_in_i[2]-1,packed_data[2*p]);
            interface_rhs_z->set_item(rank_in_i[2],packed_data[2*p+1]);
          }
          else{
            if(is_Pperiodic[2] == 0)
                interface_rhs_z->set_item(rank_in_i[2]-1,packed_data[2*p]);    
            else{
                interface_rhs_z->set_item(rank_in_i[2]-1,packed_data[2*p]);
                interface_rhs_z->set_item(rank_in_i[2],packed_data[2*p+1]);
            }
          }
          // Need to have the original rhs function assembled for corrosponding j,k pair
          double fe = pressure_assemble_local_rhs_z(i,j,t_it);

          // Setup RHS = fi - Aie*xe for solving ui
          Aie->multiply_vec_then_add(interface_rhs_z,local_rhs_z,-1.0,1.0);

          // Solve ui and transfer solution into distributed vector
          GLOBAL_EQ->DS_NavierStokes_z_solver_P(i,j,min_unknown_index(2),local_rhs_z,interface_rhs_z);
        }
     }
   }
}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: NS_pressure_update ( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DDS_NavierStokes:: NS_pressure_update" ) ;

  size_t i, j, k;

  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);
  
  for (size_t l=0;l<dim;++l)
      min_unknown_index(l) =
       PF->get_min_index_unknown_handled_by_proc( 0, l ) ;
  for (size_t l=0;l<dim;++l)
      max_unknown_index(l) =
       PF->get_max_index_unknown_handled_by_proc( 0, l ) ;

  // PF->copy_DOFs_value( 1, 0 );
  // First Equation
  double m,dirichlet_value;
  
  if(dim == 2){

    // Solve in x for pressure
    if(nb_ranks_comm_i[0]>1){
      size_t nb_send_data_x = 3*(max_unknown_index(1)-min_unknown_index(1)+1);

      for (j=min_unknown_index(1);j<=max_unknown_index(1);++j)
       {
          k=0;
          pressure_solve_x_for_secondorder(j,k,t_it,mpi_packed_data_P_x);
       }

       pressure_solve_interface_unknowns_x ( mpi_packed_data_P_x, nb_send_data_x,t_it);
    }
    else{
      for (j=min_unknown_index(1);j<=max_unknown_index(1);++j)
       {
          k=0;
          double fe = pressure_assemble_local_rhs_x(j,k,t_it);
          LA_SeqVector* local_rhs_x = GLOBAL_EQ->get_local_temp_x_P();

          if(is_Pperiodic[0] == 1)
              GLOBAL_EQ->DS_NavierStokes_x_solver_P_periodic(j,k,min_unknown_index(0),local_rhs_x,NULL);
          else
              GLOBAL_EQ->DS_NavierStokes_x_solver_P(j,k,min_unknown_index(0),local_rhs_x,NULL);
        }

    }

     // Synchronize the distributed DS solution vector
     GLOBAL_EQ->synchronize_DS_solution_vec_P();

     // Tranfer back to field
     PF->update_free_DOFs_value( 1, GLOBAL_EQ->get_solution_DS_velocity_P() ) ;

     // output_L2norm_pressure( 1 ); 

    // Solve in y for pressure
     if(nb_ranks_comm_i[1]>1){
       size_t nb_send_data_y = 3*(max_unknown_index(0)-min_unknown_index(0)+1);

       for (i=min_unknown_index(0);i<=max_unknown_index(0);++i)
       {
          k=0;
          pressure_solve_y_for_secondorder(i,k,t_it,mpi_packed_data_P_y);
       }

       pressure_solve_interface_unknowns_y ( mpi_packed_data_P_y, nb_send_data_y, t_it);
     }
     else{
      for (i=min_unknown_index(0);i<=max_unknown_index(0);++i)
       {
          k=0;
         double fe = pressure_assemble_local_rhs_y(i,k,t_it);
         LA_SeqVector* local_rhs_y = GLOBAL_EQ->get_local_temp_y_P();
         //local_rhs_y->print_items(MAC::out(),0);
         if(is_Pperiodic[1] == 1 )
            GLOBAL_EQ->DS_NavierStokes_y_solver_P_periodic(i,k,min_unknown_index(1),local_rhs_y,NULL);
         else
            GLOBAL_EQ->DS_NavierStokes_y_solver_P(i,k,min_unknown_index(1),local_rhs_y,NULL);
       }
     }

     // Synchronize the distributed DS solution vector
     GLOBAL_EQ->synchronize_DS_solution_vec_P();

     // Tranfer back to field
     PF->update_free_DOFs_value( 1 , GLOBAL_EQ->get_solution_DS_velocity_P() ) ;
     
     // output_L2norm_pressure( 1 );
   }
  else{

    
    // Solve in x for pressure
    if(nb_ranks_comm_i[0]>1){
      size_t nb_send_data_x = 3*(max_unknown_index(1)-min_unknown_index(1)+1)*(max_unknown_index(2)-min_unknown_index(2)+1);
      
      for (k=min_unknown_index(2);k<=max_unknown_index(2);++k){
        for (j=min_unknown_index(1);j<=max_unknown_index(1);++j){
            pressure_solve_x_for_secondorder(j,k,t_it,mpi_packed_data_P_x);
        }
      }
      pressure_solve_interface_unknowns_x ( mpi_packed_data_P_x, nb_send_data_x, t_it);
    }
    else{
      for (k=min_unknown_index(2);k<=max_unknown_index(2);++k){
        for (j=min_unknown_index(1);j<=max_unknown_index(1);++j){
            double fe = pressure_assemble_local_rhs_x(j,k,t_it);
            LA_SeqVector* local_rhs_x = GLOBAL_EQ->get_local_temp_x_P();
            if(is_Pperiodic[0] == 1)
                GLOBAL_EQ->DS_NavierStokes_x_solver_P_periodic(j,k,min_unknown_index(0),local_rhs_x,NULL);
            else
                GLOBAL_EQ->DS_NavierStokes_x_solver_P(j,k,min_unknown_index(0),local_rhs_x,NULL);
        }
      }
    }

    // Synchronize the distributed DS solution vector
    GLOBAL_EQ->synchronize_DS_solution_vec_P();

    // Tranfer back to field
    PF->update_free_DOFs_value( 1, GLOBAL_EQ->get_solution_DS_velocity_P() ) ;

    
    // Solve in y for pressure
    if(nb_ranks_comm_i[1]>1){
      size_t nb_send_data_y = 3*(max_unknown_index(0)-min_unknown_index(0)+1)*(max_unknown_index(2)-min_unknown_index(2)+1);
      
      // Third equation - Solve in y
      for (i=min_unknown_index(0);i<=max_unknown_index(0);++i)
       {
         for (k=min_unknown_index(2);k<=max_unknown_index(2);++k)
         {
              pressure_solve_y_for_secondorder(i,k,t_it,mpi_packed_data_P_y);
         }
       }
       pressure_solve_interface_unknowns_y ( mpi_packed_data_P_y, nb_send_data_y, t_it );
    }
    else{
      for (i=min_unknown_index(0);i<=max_unknown_index(0);++i)
       {
         for (k=min_unknown_index(2);k<=max_unknown_index(2);++k)
         {
            double fe = pressure_assemble_local_rhs_y(i,k,t_it);
            LA_SeqVector* local_rhs_y = GLOBAL_EQ->get_local_temp_y_P();
            if(is_Pperiodic[1] == 1)
                GLOBAL_EQ->DS_NavierStokes_y_solver_P_periodic(i,k,min_unknown_index(1),local_rhs_y,NULL);
            else
                GLOBAL_EQ->DS_NavierStokes_y_solver_P(i,k,min_unknown_index(1),local_rhs_y,NULL);
         }
       }
    }
     // Synchronize the distributed DS solution vector
     GLOBAL_EQ->synchronize_DS_solution_vec_P();

     // Tranfer back to field
     PF->update_free_DOFs_value( 1 , GLOBAL_EQ->get_solution_DS_velocity_P() ) ;

    
    // Solve in z for pressure

     if(nb_ranks_comm_i[2]>1){
       size_t nb_send_data_z = 3*(max_unknown_index(0)-min_unknown_index(0)+1)*(max_unknown_index(1)-min_unknown_index(1)+1);
       
       for (j=min_unknown_index(1);j<=max_unknown_index(1);++j)
       {
         for (i=min_unknown_index(0);i<=max_unknown_index(0);++i)
         {
           pressure_solve_z_for_secondorder(i,j,t_it,mpi_packed_data_P_z);
         }
       }

       pressure_solve_interface_unknowns_z ( mpi_packed_data_P_z, nb_send_data_z,t_it);
     }
     else{
       for (j=min_unknown_index(1);j<=max_unknown_index(1);++j)
       {
         for (i=min_unknown_index(0);i<=max_unknown_index(0);++i)
         {
            double fe = pressure_assemble_local_rhs_z(i,j,t_it);
            LA_SeqVector* local_rhs_z = GLOBAL_EQ->get_local_temp_z_P();
            if(is_Pperiodic[2] == 1)
                GLOBAL_EQ->DS_NavierStokes_z_solver_P_periodic(i,j,min_unknown_index(2),local_rhs_z,NULL);
            else
                GLOBAL_EQ->DS_NavierStokes_z_solver_P(i,j,min_unknown_index(2),local_rhs_z,NULL);
         }
       }
     }

     // Synchronize the distributed DS solution vector
     GLOBAL_EQ->synchronize_DS_solution_vec_P();

     // Tranfer back to field
     PF->update_free_DOFs_value( 1, GLOBAL_EQ->get_solution_DS_velocity_P() ) ;
  }


  // Debug
  output_L2norm_pressure( 0 );
  output_L2norm_pressure( 1 );  
}




//---------------------------------------------------------------------------
void
DDS_NavierStokes:: NS_final_step ( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL(
   "DDS_NavierStokes:: NS_final_step" ) ;

   size_t i, j, k,m;

  double value=0.;

  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);

  double xhr,xright,xvalue1=0.,xvalue2=0.,xvalue=0.;
  double yhr,yright,yvalue1=0.,yvalue2=0.,yvalue=0.;
  double dirichlet_value;
  FV_SHIFT_TRIPLET shift = PF->shift_staggeredToCentered() ;

    // Get local min and max indices 
    // When we are running in parallel, the unknowns in the overlapping region are not solved. So we need to include 
    // them here by calling get_min_index_unknown_on_proc() instead of get_min_index_unknown_handled_by_proc().


    for (size_t l=0;l<dim;++l)
      min_unknown_index(l) =
       PF->get_min_index_unknown_on_proc( 0, l ) ;
    for (size_t l=0;l<dim;++l)
      max_unknown_index(l) =
       PF->get_max_index_unknown_on_proc( 0, l ) ;

    for (i=min_unknown_index(0);i<=max_unknown_index(0);++i)
    {
      for (j=min_unknown_index(1);j<=max_unknown_index(1);++j)
      {
        if(dim ==2 )
        {
          k=0;
          
          xhr= UF->get_DOF_coordinate( shift.i+i,0, 0 ) - UF->get_DOF_coordinate( shift.i+i-1, 0, 0 ) ;

          // Divergence of un (x component)
          xright = UF->DOF_value( shift.i+i, j, k, 0, 0 ) - UF->DOF_value( shift.i+i-1, j, k, 0, 0 ) ;
          xvalue1 = xright/xhr;

          // Divergence of un+1 (x component)
          xright = UF->DOF_value( shift.i+i, j, k, 0, 1 ) - UF->DOF_value( shift.i+i-1, j, k, 0, 1 ) ;
          xvalue2 = xright/xhr;
        
          xvalue = 0.5*(xvalue1+xvalue2);

          yhr= UF->get_DOF_coordinate( shift.j+j,1, 1 ) - UF->get_DOF_coordinate( shift.j+j-1, 1, 1 ) ;
          
          // Divergence of un (y component)
          yright = UF->DOF_value( i, shift.j+j, k, 1, 0 ) - UF->DOF_value( i, shift.j+j-1, k, 1, 0 ) ;
          yvalue1 = yright/yhr;

          // Divergence of un+1 (y component)
          yright = UF->DOF_value( i, shift.j+j, k, 1, 1 ) - UF->DOF_value( i, shift.j+j-1, k, 1, 1 ) ;
          yvalue2 = yright/yhr;

          yvalue = 0.5*(yvalue1+yvalue2);

          // Assemble the bodyterm
          value = PF->DOF_value( i, j, k, 0, 0 ) + PF->DOF_value( i, j, k, 0, 1 ) - kai*mu*(xvalue + yvalue);

          PF->set_DOF_value( i, j, k, 0, 0, value);
          //MAC::out() << PF->DOF_value( i, j, k, 0, 1 ) << endl;
        }
        else
        {
          for (k=min_unknown_index(2);k<=max_unknown_index(2);++k)
          {

            xhr= UF->get_DOF_coordinate( shift.i+i,0, 0 ) - UF->get_DOF_coordinate( shift.i+i-1, 0, 0 ) ;

            // Divergence of un (x component)
            xright = UF->DOF_value( shift.i+i, j, k, 0, 0 ) - UF->DOF_value( shift.i+i-1, j, k, 0, 0 ) ;
            xvalue1 = xright/xhr;

            // Divergence of un+1 (x component)
            xright = UF->DOF_value( shift.i+i, j, k, 0, 1 ) - UF->DOF_value( shift.i+i-1, j, k, 0, 1 ) ;
            xvalue2 = xright/xhr;
          
            xvalue = 0.5*(xvalue1+xvalue2);

            yhr= UF->get_DOF_coordinate( shift.j+j,1, 1 ) - UF->get_DOF_coordinate( shift.j+j-1, 1, 1 ) ;
  
            // Divergence of un (y component)            
            yright = UF->DOF_value( i, shift.j+j, k, 1, 0 ) - UF->DOF_value( i, shift.j+j-1, k, 1, 0 ) ;
            yvalue1 = yright/yhr;

            // Divergence of un+1 (y component)
            yright = UF->DOF_value( i, shift.j+j, k, 1, 1 ) - UF->DOF_value( i, shift.j+j-1, k, 1, 1 ) ;
            yvalue2 = yright/yhr;

            yvalue = 0.5*(yvalue1+yvalue2);

            double zhr,zright,zvalue1=0.,zvalue2=0.,zvalue=0.;

            zhr= UF->get_DOF_coordinate( shift.k+k,2, 2 ) - UF->get_DOF_coordinate( shift.k+k-1, 2, 2 ) ;
            
            // Divergence of un (z component)
            zright = UF->DOF_value( i, j, shift.k+k, 2, 0 ) - UF->DOF_value( i, j, shift.k+k-1, 2, 0 ) ;
            zvalue1 = zright/zhr;

            // Divergence of un+1 (z component)
            zright = UF->DOF_value( i, j, shift.k+k, 2, 1 ) - UF->DOF_value( i, j, shift.k+k-1, 2, 1 ) ;
            zvalue2 = zright/zhr;            
            
            zvalue = 0.5*(zvalue1+zvalue2);

            // Assemble the bodyterm
            value = PF->DOF_value( i, j, k, 0, 0 ) + PF->DOF_value( i, j, k, 0, 1 ) - kai*mu*(xvalue + yvalue+ zvalue);

            PF->set_DOF_value( i, j, k, 0, 0, value);

          }
        }
      }

    }

    // Propagate values to the boundaries depending on BC conditions
    PF->set_neumann_DOF_values();

    // Un=Un+1
    //UF->copy_DOFs_value( 1, 2 );

    UF->copy_DOFs_value( 0, 1 );

}




//----------------------------------------------------------------------
void
DDS_NavierStokes::write_pressure_field(FV_TimeIterator const* t_it)
//----------------------------------------------------------------------
{
  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);
  size_t i,j,k;
  for (size_t l=0;l<dim;++l)
    min_unknown_index(l) =
     PF->get_min_index_unknown_on_proc( 0, l ) ;
  for (size_t l=0;l<dim;++l)
    max_unknown_index(l) =
     PF->get_max_index_unknown_on_proc( 0, l ) ;

  double value;
  ofstream pressureFile ;

  ostringstream fileNameStream;                    // let this be empty

  fileNameStream << "/home/arun95/NS_2D_lidCavity_results/" << t_it->time_step() << "output_pressure.csv"; // and pass "dice_" here
  string fileName = fileNameStream.str(); 

  MAC::out()<<"The file name is "<<fileName<<endl;

  pressureFile.open(fileName.c_str(),std::ios_base::out | std::ios_base::trunc) ;

  for (i=min_unknown_index(0);i<=max_unknown_index(0);++i)
  {
    for (j=min_unknown_index(1);j<=max_unknown_index(1);++j)
    {
      if(dim ==2 )
      {
        k=0;
        value = PF->DOF_value( i, j, k, 0, 0 );
           pressureFile << value<<endl ;
      }
      else{
        for (k=min_unknown_index(2);k<=max_unknown_index(2);++k){
          value = PF->DOF_value( i, j, k, 0, 0 );
           pressureFile << value<<endl ;
        }
      }
    }
  }
  pressureFile.close();
}




//----------------------------------------------------------------------
void
DDS_NavierStokes::write_velocity_field(FV_TimeIterator const* t_it)
//----------------------------------------------------------------------
{

  double value;
  ofstream velocityFile ;

  ostringstream fileNameStream;                    // let this be empty
  fileNameStream << "/home/arun95/NS_2D_lidCavity_results/" << t_it->time_step() << "output_velocity.csv"; // and pass "dice_" here
  string fileName = fileNameStream.str(); 

  velocityFile.open(fileName.c_str(), std::ios_base::out | std::ios_base::trunc) ;
  size_t i,j,k;

  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);
  size_t nb_comps = UF->nb_components() ;

  for(size_t comp=0;comp<nb_comps;comp++)
  {
    // Get local min and max indices
    for (size_t l=0;l<dim;++l)
      min_unknown_index(l) =
       UF->get_min_index_unknown_handled_by_proc( comp, l ) ;
    for (size_t l=0;l<dim;++l)
      max_unknown_index(l) =
       UF->get_max_index_unknown_handled_by_proc( comp, l ) ;

    for (i=min_unknown_index(0);i<=max_unknown_index(0);++i)
    {
      for (j=min_unknown_index(1);j<=max_unknown_index(1);++j)
      {
        if(dim ==2 )
        {
          k=0;
          value = UF->DOF_value( i, j, k, comp, 0 );
             velocityFile << value<<endl ;
        }
        else{
          for (k=min_unknown_index(2);k<=max_unknown_index(2);++k){
            value = UF->DOF_value( i, j, k, comp, 0 );
             velocityFile << value<<endl ;
          }
        }
      }
    }
  }
  velocityFile.close();
}




//----------------------------------------------------------------------
void
DDS_NavierStokes::get_velocity_divergence(void)
//----------------------------------------------------------------------
{

  size_t i,j,k;
  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);
  double dux, duy, dx, dy, dz, duz=0.;
  double div_velocity = 0.;
  double cell_div=0.,max_divu=0.,divu = 0.;

  FV_SHIFT_TRIPLET shift = PF->shift_staggeredToCentered() ;
  for (size_t l=0;l<dim;++l)
    min_unknown_index(l) =
	PF->get_min_index_unknown_handled_by_proc( 0, l ) ;
  for (size_t l=0;l<dim;++l)
    max_unknown_index(l) =
	PF->get_max_index_unknown_handled_by_proc( 0, l ) ;

  for (i=min_unknown_index(0);i<=max_unknown_index(0);++i)
  {
    dx = PF->get_cell_size( i, 0, 0 );
    for (j=min_unknown_index(1);j<=max_unknown_index(1);++j)
    {
      dy = PF->get_cell_size( j, 0, 1 );
      if( dim == 2 )
      {
        k=0;

        // Divergence of u (x component)
        dux = UF->DOF_value( shift.i+i, j, k, 0, 0 ) 
		- UF->DOF_value( shift.i+i-1, j, k, 0, 0 ) ;

        // Divergence of u (y component)
        duy = UF->DOF_value( i, shift.j+j, k, 1, 0 ) 
		- UF->DOF_value( i, shift.j+j-1, k, 1, 0 ) ;

        cell_div = dux * dy + duy * dx;		
	max_divu = MAC::max( MAC::abs(cell_div) / ( dx * dy ), max_divu );
        div_velocity += cell_div * cell_div / ( dx * dy );
      }
      else
      {
        for (k=min_unknown_index(2);k<=max_unknown_index(2);++k)
	{
          dz = PF->get_cell_size( k, 0, 2 );
	  
	  // Divergence of u (x component)
          dux = UF->DOF_value( shift.i+i, j, k, 0, 0 ) 
	  	- UF->DOF_value( shift.i+i-1, j, k, 0, 0 ) ;

          // Divergence of u (y component)
          duy = UF->DOF_value( i, shift.j+j, k, 1, 0 ) 
	  	- UF->DOF_value( i, shift.j+j-1, k, 1, 0 ) ;
           
          // Divergence of u(z component)
          duz = UF->DOF_value( i, j, shift.k+k, 2, 0 ) 
	  	- UF->DOF_value( i, j, shift.k+k-1, 2, 0 ) ;
          

          cell_div = dux * dy * dz + duy * dx * dz + duz * dx * dy ;		
	  max_divu = MAC::max( MAC::abs(cell_div) / ( dx * dy * dz ), 
	  	max_divu );
          div_velocity += cell_div * cell_div / ( dx * dy * dz );
        }
      }
    }
  }

  FV_Mesh const* primary_mesh = UF->primary_grid() ;
  double domain_measure = dim == 2 ? 
  	primary_mesh->get_main_domain_boundary_perp_to_direction_measure( 0 )
  	* primary_mesh->get_main_domain_boundary_perp_to_direction_measure( 1 ):
	primary_mesh->get_main_domain_boundary_perp_to_direction_measure( 0 )
	* ( primary_mesh->get_main_domain_max_coordinate(2)
		- primary_mesh->get_main_domain_min_coordinate(2) );
  
  div_velocity = pelCOMM->sum( div_velocity ) ;
//  div_velocity = MAC::sqrt( div_velocity / domain_measure );
  div_velocity = MAC::sqrt( div_velocity );
  max_divu = pelCOMM->max( max_divu ) ;
  if ( my_rank == is_master )
    MAC::out() << "Norm L2 div(u) = "<< MAC::doubleToString( ios::scientific, 12, div_velocity ) << " Max div(u) = " << MAC::doubleToString( ios::scientific, 12, max_divu ) << endl;
}




//----------------------------------------------------------------------
void
DDS_NavierStokes::output_L2norm_pressure( size_t level )
//----------------------------------------------------------------------
{
  double value;
  size_t i,j,k;

  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);

  double dx,dy;
  double L2normP = 0.;
  double cell_P=0.,max_P=0.;

  for (size_t l=0;l<dim;++l)
      min_unknown_index(l) =
       PF->get_min_index_unknown_handled_by_proc( 0, l ) ;
    for (size_t l=0;l<dim;++l)
      max_unknown_index(l) =
       PF->get_max_index_unknown_handled_by_proc( 0, l ) ;

  for (i=min_unknown_index(0);i<=max_unknown_index(0);++i)
  {
    for (j=min_unknown_index(1);j<=max_unknown_index(1);++j)
    {
      if(dim ==2 )
      {
        k=0;
        dx = PF->get_cell_size( i,0, 0 );
        dy = PF->get_cell_size( j,0, 1 );
        cell_P = PF->DOF_value( i, j, k, 0, level );
	max_P = MAC::max( MAC::abs(cell_P), max_P );
        L2normP += cell_P * cell_P * dx * dy;
      }
      else
      {
        double dz=0.;
        for (k=min_unknown_index(2);k<=max_unknown_index(2);++k)
	{          
          dx = PF->get_cell_size( i,0, 0 );
          dy = PF->get_cell_size( j,0, 1 );
          dz = PF->get_cell_size( k,0, 2 );
          cell_P = PF->DOF_value( i, j, k, 0, level );
	  max_P = MAC::max( MAC::abs(cell_P), max_P );
          L2normP += cell_P * cell_P * dx * dy * dz;
        }
      }
    }
  }

  FV_Mesh const* primary_mesh = UF->primary_grid() ;
  double domain_measure = dim == 2 ? 
  	primary_mesh->get_main_domain_boundary_perp_to_direction_measure( 0 )
  	* primary_mesh->get_main_domain_boundary_perp_to_direction_measure( 1 ):
	primary_mesh->get_main_domain_boundary_perp_to_direction_measure( 0 )
	* ( primary_mesh->get_main_domain_max_coordinate(2)
		- primary_mesh->get_main_domain_min_coordinate(2) );
  
  L2normP = pelCOMM->sum( L2normP ) ;
//  L2normP = MAC::sqrt( L2normP / domain_measure );
  L2normP = MAC::sqrt( L2normP );  
  max_P = pelCOMM->max( max_P ) ;
  if ( my_rank == is_master )
      MAC::out()<< "Norm L2 P = "<< MAC::doubleToString( ios::scientific, 12, L2normP ) << " Max P = " << MAC::doubleToString( ios::scientific, 12, max_P ) << endl;
      
}




//----------------------------------------------------------------------
void
DDS_NavierStokes::output_L2norm_velocity( size_t level )
//----------------------------------------------------------------------
{
  double value;
  size_t i,j,k;

  size_t_vector min_unknown_index(dim,0);
  size_t_vector max_unknown_index(dim,0);

  double dx,dy;

  for(size_t comp=0;comp<nb_comps;comp++){
    for (size_t l=0;l<dim;++l)
      min_unknown_index(l) =
       UF->get_min_index_unknown_handled_by_proc( comp, l ) ;
    for (size_t l=0;l<dim;++l)
      max_unknown_index(l) =
       UF->get_max_index_unknown_handled_by_proc( comp, l ) ;


    double L2normU = 0.;
    double cell_U=0.,max_U=0.;

    for (i=min_unknown_index(0);i<=max_unknown_index(0);++i)
    {
      for (j=min_unknown_index(1);j<=max_unknown_index(1);++j)
      {
        if(dim ==2 )
        {
          k=0;
          dx = UF->get_cell_size( i,comp, 0 );
          dy = UF->get_cell_size( j,comp, 1 );
          cell_U = UF->DOF_value( i, j, k, comp, level );
    max_U = MAC::max( MAC::abs(cell_U), max_U );
          L2normU += cell_U * cell_U * dx * dy;
        }
        else
        {
          double dz=0.;
          for (k=min_unknown_index(2);k<=max_unknown_index(2);++k)
    {          
            dx = UF->get_cell_size( i,comp, 0 );
            dy = UF->get_cell_size( j,comp, 1 );
            dz = UF->get_cell_size( k,comp, 2 );
            cell_U = UF->DOF_value( i, j, k, comp, level );
      max_U = MAC::max( MAC::abs(cell_U), max_U );
            L2normU += cell_U * cell_U * dx * dy * dz;
          }
        }
      }
    }

    FV_Mesh const* primary_mesh = UF->primary_grid() ;
    double domain_measure = dim == 2 ? 
      primary_mesh->get_main_domain_boundary_perp_to_direction_measure( 0 )
      * primary_mesh->get_main_domain_boundary_perp_to_direction_measure( 1 ):
    primary_mesh->get_main_domain_boundary_perp_to_direction_measure( 0 )
    * ( primary_mesh->get_main_domain_max_coordinate(2)
      - primary_mesh->get_main_domain_min_coordinate(2) );
    
    L2normU = pelCOMM->sum( L2normU ) ;
  //  L2normP = MAC::sqrt( L2normP / domain_measure );
    L2normU = MAC::sqrt( L2normU );  
    max_U = pelCOMM->max( max_U ) ;
    if ( my_rank == is_master )
        MAC::out()<< "Component: "<<comp<< " Norm L2 U = "<< MAC::doubleToString( ios::scientific, 12, L2normU ) << " Max U = " << MAC::doubleToString( ios::scientific, 12, max_U ) << endl;
    
  }

      
}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: create_DDS_subcommunicators ( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatEquation:: create_DDS_subcommunicators" ) ;

   int color = 0, key = 0;
   int const* MPI_coordinates_world = UF->primary_grid()->get_MPI_coordinates() ;
   int const* MPI_number_of_coordinates = UF->primary_grid()->get_domain_decomposition() ;

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
DDS_NavierStokes:: processor_splitting ( int color, int key, size_t const dir )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_HeatEquation:: processor_splitting" ) ;

   MPI_Comm_split(MPI_COMM_WORLD, color, key, &DDS_Comm_i[dir]);
   MPI_Comm_size( DDS_Comm_i[dir], &nb_ranks_comm_i[dir] ) ;
   MPI_Comm_rank( DDS_Comm_i[dir], &rank_in_i[dir] ) ;

}

//---------------------------------------------------------------------------
void
DDS_NavierStokes:: free_DDS_subcommunicators ( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: free_DDS_subcommunicators" ) ;
}




//----------------------------------------------------------------------
double 
DDS_NavierStokes:: assemble_advection_Upwind( 
  size_t advecting_level, double const& coef, size_t advected_level,
  size_t const& i, size_t const& j, size_t const& k, size_t const& component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "DDS_NavierStokes:: assemble_advection_Upwind" );   

   // Parameters
   double dxC = 0., dyC = 0., dzC = 0.;
   double AdvectedValueC = 0., AdvectedValueRi = 0., AdvectedValueLe = 0.,
    AdvectedValueTo = 0., AdvectedValueBo = 0., AdvectedValueFr = 0., 
  AdvectedValueBe = 0,
    AdvectorValueC = 0., AdvectorValueRi = 0., AdvectorValueLe = 0.,
    AdvectorValueTo = 0., AdvectorValueBo = 0., AdvectorValueFr = 0., 
  AdvectorValueBe = 0, AdvectorValueToLe = 0., AdvectorValueToRi = 0., 
  AdvectorValueBoLe = 0., AdvectorValueBoRi = 0., AdvectorValueFrLe = 0., 
  AdvectorValueFrRi = 0., AdvectorValueBeLe = 0., AdvectorValueBeRi = 0.,
  AdvectorValueFrTo = 0., AdvectorValueFrBo = 0., AdvectorValueBeTo = 0., 
  AdvectorValueBeBo = 0.,
  ur = 0., ul = 0., vt = 0., vb = 0., wf = 0., wb = 0.,
  fri = 0., fle = 0., fto = 0., fbo = 0., ffr = 0., fbe = 0., flux = 0.;
   
  // Comment: staggered unknowns always have a defined value at +1/-1
   // indices in directions different from their component number, 
   // i.e. u in x, v in y and w in z. 
   // For instance, if u on the right or left boundary is an unknown with
   // homogeneous Neumann BC, then the flux on the right or left needs special 
   // treatment using the center value.
   // Otherwise, whether one of the +1/-1 DOF values is on a
   // boundary or not, and whether that boundary has a Dirichlet or Neumann 
   // condition is irrelevant, this +1/-1 DOF always has the right value. 
   // For Neumann, this is guaranted by 
   // FV_BoundaryCondition:: set_free_DOF_values in 
   // FV_DiscreteField:: update_free_DOFs_value or 
   // FV_DiscreteField:: add_to_free_DOFs_value
   FV_SHIFT_TRIPLET shift = UF->shift_staggeredToStaggered( component );

   dxC = UF->get_cell_size(i,component,0) ;          
   dyC = UF->get_cell_size(j,component,1) ;
 
   if( dim == 2 )
   {

     AdvectedValueC = UF->DOF_value( i, j, k, component, 
      advected_level );
     AdvectorValueC = UF->DOF_value( i, j, k, component, 
      advecting_level );
     
     // The First Component (u)
     if ( component == 0 )
     {
       // Right (U_X)
       if ( UF->DOF_color( i, j, k, component ) == FV_BC_RIGHT )
         fri = AdvectorValueC * AdvectedValueC;
       else
       {
               AdvectedValueRi = UF->DOF_value(
                   i+1, j, k, component, advected_level );
         AdvectorValueRi = UF->DOF_value( 
             i+1, j, k, component, advecting_level );
         ur = 0.5 * ( AdvectorValueC + AdvectorValueRi );
         if ( ur > 0. ) fri = ur * AdvectedValueC;
         else fri = ur * AdvectedValueRi;
       }
           
       // Left (U_X)
       if ( UF->DOF_color( i, j, k, component ) == FV_BC_LEFT )
         fle = AdvectorValueC * AdvectedValueC;
       else
       {
               AdvectedValueLe = UF->DOF_value( 
               i-1, j, k, component, advected_level );
         AdvectorValueLe = UF->DOF_value( 
             i-1, j, k, component, advecting_level );
         ul = 0.5 * ( AdvectorValueC + AdvectorValueLe );
         if ( ul > 0. ) fle = ul * AdvectedValueLe;
         else fle = ul * AdvectedValueC;
       }
       
       // Top (U_Y)
       AdvectedValueTo = UF->DOF_value(
               i, j+1, k, component, advected_level );
       AdvectorValueToLe = UF->DOF_value(
             i+shift.i-1, j+shift.j, k, 1, advecting_level );
       AdvectorValueToRi = UF->DOF_value(
             i+shift.i, j+shift.j, k, 1, advecting_level );
       vt = 0.5 * ( AdvectorValueToLe + AdvectorValueToRi );
       if ( vt > 0. ) fto = vt * AdvectedValueC;
       else fto = vt * AdvectedValueTo;
   
       // Bottom (U_Y)
       AdvectedValueBo = UF->DOF_value(
               i, j-1, k, component, advected_level );
       AdvectorValueBoLe = UF->DOF_value(
             i+shift.i-1, j+shift.j-1, k, 1, advecting_level );
       AdvectorValueBoRi = UF->DOF_value(
             i+shift.i, j+shift.j-1, k, 1, advecting_level );
       vb = 0.5 * ( AdvectorValueBoLe + AdvectorValueBoRi );
       if ( vb > 0. ) fbo = vb * AdvectedValueBo;
       else fbo = vb * AdvectedValueC;
     }
     
     // The second Component (v)
     else
     {
       // Right (V_X)
       AdvectedValueRi = UF->DOF_value(
               i+1, j, k, component, advected_level );
       AdvectorValueToRi = UF->DOF_value(
             i+shift.i, j+shift.j, k, 0, advecting_level );
       AdvectorValueBoRi = UF->DOF_value(
             i+shift.i, j+shift.j-1, k, 0, advecting_level );
       ur = 0.5 * ( AdvectorValueToRi + AdvectorValueBoRi );
       if ( ur > 0. ) fri = ur * AdvectedValueC;
       else fri = ur * AdvectedValueRi;
           
       // Left (V_X)
       AdvectedValueLe = UF->DOF_value(
               i-1, j, k, component, advected_level );
       AdvectorValueToLe = UF->DOF_value(
             i+shift.i-1, j+shift.j, k, 0, advecting_level );
       AdvectorValueBoLe = UF->DOF_value(
             i+shift.i-1, j+shift.j-1, k, 0, advecting_level );
       ul = 0.5 * ( AdvectorValueToLe + AdvectorValueBoLe );
       if ( ul > 0. ) fle = ul * AdvectedValueLe;
       else fle = ul * AdvectedValueC;
   
       // Top (V_Y)
       if ( UF->DOF_color(i, j, k, component ) == FV_BC_TOP )
         fto = AdvectorValueC * AdvectedValueC;
       else
       {
         AdvectedValueTo = UF->DOF_value(
               i, j+1, k, component, advected_level );
         AdvectorValueTo = UF->DOF_value(
             i, j+1, k, component, advecting_level );
         vt = 0.5 * ( AdvectorValueTo + AdvectorValueC );
         if ( vt > 0. ) fto = vt * AdvectedValueC;
         else fto = vt * AdvectedValueTo;
       }
   
       // Bottom (V_Y)
       if ( UF->DOF_color(i, j, k, component ) == FV_BC_BOTTOM )
         fbo = AdvectorValueC * AdvectedValueC;
       else
       {
         AdvectedValueBo = UF->DOF_value(
               i, j-1, k, component, advected_level );
         AdvectorValueBo = UF->DOF_value(
             i, j-1, k, component, advecting_level );
         vb = 0.5 * ( AdvectorValueBo + AdvectorValueC );
         if ( vb > 0. ) fbo = vb * AdvectedValueBo;
         else fbo = vb * AdvectedValueC;
       }
     }

           flux = (fto - fbo) * dxC + (fri - fle) * dyC;
   }
   else // DIM = 3
   {
     dzC = UF->get_cell_size(k,component,2) ;       
    AdvectedValueC = UF->DOF_value( i, j, k, component, 
        advected_level );
    AdvectorValueC = UF->DOF_value( i, j, k, component, 
        advecting_level );

       // The First Component (u)
       if ( component == 0 )
       {
         // Right (U_X)
         if ( UF->DOF_color( i, j, k, component ) == FV_BC_RIGHT )
           fri = AdvectorValueC * AdvectedValueC;
         else
         {
                 AdvectedValueRi = UF->DOF_value(
           i+1, j, k, component, advected_level );
     AdvectorValueRi = UF->DOF_value(
         i+1, j, k, component, advecting_level );
                 ur = 0.5 * ( AdvectorValueRi + AdvectorValueC );
                 if ( ur > 0. ) fri = ur * AdvectedValueC;
                 else fri = ur * AdvectedValueRi;
         }
           
         // Left (U_X)
         if ( UF->DOF_color( i, j, k, component ) == FV_BC_LEFT )
           fle = AdvectorValueC * AdvectedValueC;
         else
         {
                 AdvectedValueLe = UF->DOF_value(
           i-1, j, k, component, advected_level );
     AdvectorValueLe = UF->DOF_value(
         i-1, j, k, component, advecting_level );
                 ul = 0.5 * ( AdvectorValueLe + AdvectorValueC );
                 if ( ul > 0. ) fle = ul * AdvectedValueLe;
                 else fle = ul * AdvectedValueC;
         }
   
         // Top (U_Y)
         AdvectedValueTo = UF->DOF_value(
           i, j+1, k, component, advected_level );
         AdvectorValueToLe = UF->DOF_value( 
         i+shift.i-1, j+shift.j, k, 1, advecting_level );
         AdvectorValueToRi = UF->DOF_value( 
         i+shift.i, j+shift.j, k, 1, advecting_level );
         vt = 0.5 * ( AdvectorValueToLe + AdvectorValueToRi );
         if ( vt > 0. ) fto = vt * AdvectedValueC;
         else fto = vt * AdvectedValueTo;
   
         // Bottom (U_Y)
         AdvectedValueBo = UF->DOF_value(
           i, j-1, k, component, advected_level );
         AdvectorValueBoLe = UF->DOF_value(
         i+shift.i-1, j+shift.j-1, k, 1, advecting_level );
         AdvectorValueBoRi = UF->DOF_value(
         i+shift.i, j+shift.j-1, k, 1, advecting_level );
         vb = 0.5 * ( AdvectorValueBoLe + AdvectorValueBoRi );
         if ( vb > 0. ) fbo = vb * AdvectedValueBo;
         else fbo = vb * AdvectedValueC;
      
         // Front (U_Z)
         AdvectedValueFr = UF->DOF_value(
           i, j, k+1, component, advected_level );
         AdvectorValueFrLe = UF->DOF_value(
         i+shift.i-1, j, k+shift.k, 2, advecting_level );
         AdvectorValueFrRi = UF->DOF_value(
         i+shift.i, j, k+shift.k, 2, advecting_level );
         wf = 0.5 * ( AdvectorValueFrLe + AdvectorValueFrRi );
         if ( wf > 0. ) ffr = wf * AdvectedValueC;
         else ffr = wf * AdvectedValueFr;
     
         // Behind (U_Z)
         AdvectedValueBe = UF->DOF_value(
           i, j, k-1, component, advected_level );
         AdvectorValueBeLe = UF->DOF_value(
         i+shift.i-1, j, k+shift.k-1, 2, advecting_level );
         AdvectorValueBeRi = UF->DOF_value(
         i+shift.i, j, k+shift.k-1, 2, advecting_level );
         wb = 0.5 * ( AdvectorValueBeLe + AdvectorValueBeRi );
         if ( wb > 0. ) fbe = wb * AdvectedValueBe;
         else fbe = wb * AdvectedValueC;
       }
       
       // The Second Component (v)
       else if ( component == 1 )
       {
         // Right (V_X)
         AdvectedValueRi = UF->DOF_value(
           i+1, j, k, component, advected_level );
         AdvectorValueToRi = UF->DOF_value(
         i+shift.i, j+shift.j, k, 0, advecting_level );
         AdvectorValueBoRi = UF->DOF_value(
         i+shift.i, j+shift.j-1, k, 0, advecting_level );
         ur = 0.5 * ( AdvectorValueToRi + AdvectorValueBoRi );
         if ( ur > 0. ) fri = ur * AdvectedValueC;
         else fri = ur * AdvectedValueRi;
           
         // Left (V_X)
         AdvectedValueLe = UF->DOF_value(
           i-1, j, k, component, advected_level );
         AdvectorValueToLe = UF->DOF_value(
         i+shift.i-1, j+shift.j, k, 0, advecting_level );
         AdvectorValueBoLe = UF->DOF_value(
         i+shift.i-1, j+shift.j-1, k, 0, advecting_level );
         ul = 0.5 * ( AdvectorValueToLe + AdvectorValueBoLe );
         if ( ul > 0. ) fle = ul * AdvectedValueLe;
         else fle = ul * AdvectedValueC;
   
         // Top (V_Y)
         if ( UF->DOF_color(i, j, k, component ) == FV_BC_TOP )
           fto = AdvectorValueC * AdvectedValueC;
         else
         {
           AdvectedValueTo = UF->DOF_value(
           i, j+1, k, component, advected_level );
     AdvectorValueTo = UF->DOF_value(
         i, j+1, k, component, advecting_level );
           vt = 0.5 * ( AdvectorValueTo + AdvectorValueC );
           if ( vt > 0. ) fto = vt * AdvectedValueC;
           else fto = vt * AdvectedValueTo;
         }
   
         // Bottom (V_Y)
         if ( UF->DOF_color(i, j, k, component ) == FV_BC_BOTTOM )
           fbo = AdvectorValueC * AdvectedValueC;
         else
         {
           AdvectedValueBo = UF->DOF_value(
           i, j-1, k, component, advected_level );
     AdvectorValueBo = UF->DOF_value(
         i, j-1, k, component, advecting_level );
           vb = 0.5 * ( AdvectorValueBo + AdvectorValueC );
           if ( vb > 0. ) fbo = vb * AdvectedValueBo;
           else fbo = vb * AdvectedValueC;
         }
      
         // Front (V_Z)
         AdvectedValueFr = UF->DOF_value(
           i, j, k+1, component, advected_level );
         AdvectorValueFrTo = UF->DOF_value(
         i, j+shift.j, k+shift.k, 2, advecting_level );
         AdvectorValueFrBo = UF->DOF_value(
         i, j+shift.j-1, k+shift.k, 2, advecting_level );
         wf = 0.5 * ( AdvectorValueFrTo + AdvectorValueFrBo );
         if ( wf > 0. ) ffr = wf * AdvectedValueC;
         else ffr = wf * AdvectedValueFr;
     
         // Behind (V_Z)
         AdvectedValueBe = UF->DOF_value(
           i, j, k-1, component, advected_level );
         AdvectorValueBeTo = UF->DOF_value(
         i, j+shift.j, k+shift.k-1, 2, advecting_level );
         AdvectorValueBeBo = UF->DOF_value(
         i, j+shift.j-1, k+shift.k-1, 2, advecting_level );
         wb = 0.5 * ( AdvectorValueBeTo + AdvectorValueBeBo );
         if ( wb > 0. ) fbe = wb * AdvectedValueBe;
         else fbe = wb * AdvectedValueC;
       }
       
       // The Third Component (w)
       else
       {
         // Right (W_X)
         AdvectedValueRi = UF->DOF_value(
           i+1, j, k, component, advected_level );
         AdvectorValueFrRi = UF->DOF_value(
           i+shift.i, j, k+shift.k, 0, advecting_level );
         AdvectorValueBeRi = UF->DOF_value(
         i+shift.i, j, k+shift.k-1, 0, advecting_level );
         ur = 0.5 * ( AdvectorValueFrRi + AdvectorValueBeRi );
         if ( ur > 0. ) fri = ur * AdvectedValueC;
         else fri = ur * AdvectedValueRi;
           
         // Left (W_X)
         AdvectedValueLe = UF->DOF_value(
           i-1, j, k, component, advected_level );
         AdvectorValueFrLe = UF->DOF_value(
         i+shift.i-1, j, k+shift.k, 0, advecting_level );
         AdvectorValueBeLe = UF->DOF_value(
         i+shift.i-1, j, k+shift.k-1, 0, advecting_level );
         ul = 0.5 * ( AdvectorValueFrLe + AdvectorValueBeLe );
         if ( ul > 0. ) fle = ul * AdvectedValueLe;
         else fle = ul * AdvectedValueC;

         // Top (W_Y)
         AdvectedValueTo = UF->DOF_value(
           i, j+1, k, component, advected_level );
         AdvectorValueFrTo = UF->DOF_value(
        i, j+shift.j, k+shift.k, 1, advecting_level );
         AdvectorValueBeTo = UF->DOF_value(
        i, j+shift.j, k+shift.k-1, 1, advecting_level );
         vt = 0.5 * ( AdvectorValueFrTo + AdvectorValueBeTo );
         if ( vt > 0. ) fto = vt * AdvectedValueC;
         else fto = vt * AdvectedValueTo;
   
         // Bottom (W_Y)
         AdvectedValueBo = UF->DOF_value(
           i, j-1, k, component, advected_level );
         AdvectorValueFrBo = UF->DOF_value(
         i, j+shift.j-1, k+shift.k, 1, advecting_level );
         AdvectorValueBeBo = UF->DOF_value(
         i, j+shift.j-1, k+shift.k-1, 1, advecting_level );
         vb = 0.5 * ( AdvectorValueFrBo + AdvectorValueBeBo );
         if ( vb > 0. ) fbo = vb * AdvectedValueBo;
         else fbo = vb * AdvectedValueC;
      
         // Front (W_Z)
         if ( UF->DOF_color( i, j, k, component ) == FV_BC_FRONT )
           ffr = AdvectorValueC * AdvectedValueC;
         else
         {
           AdvectedValueFr = UF->DOF_value(
           i, j, k+1, component, advected_level );
           AdvectorValueFr = UF->DOF_value(
         i, j, k+1, component, advecting_level );
           wf = 0.5 * ( AdvectorValueFr + AdvectorValueC );
           if ( wf > 0. ) ffr = wf * AdvectedValueC;
           else ffr = wf * AdvectedValueFr;
         }
     
         // Behind (W_Z)
         if ( UF->DOF_color( i, j, k, component ) == FV_BC_BEHIND )
           fbe = AdvectorValueC * AdvectedValueC;
         else
         {
           AdvectedValueBe = UF->DOF_value(
           i, j, k-1, component, advected_level );
           AdvectorValueBe = UF->DOF_value(
         i, j, k-1, component, advecting_level );
           wb = 0.5 * ( AdvectorValueBe + AdvectorValueC );
           if ( wb > 0. ) fbe = wb * AdvectedValueBe;
           else fbe = wb * AdvectedValueC;
         }
       }
       
             flux = (fto - fbo) * dxC * dzC
          + (fri - fle) * dyC * dzC
      + (ffr - fbe) * dxC * dyC;

     }   
     
   return ( coef * flux ); 

}




//----------------------------------------------------------------------
double
DDS_NavierStokes:: assemble_advection_TVD( 
  size_t advecting_level, double const& coef, size_t advected_level,
         size_t const& i, size_t const& j, size_t const& k, size_t const& component ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Staggered:: assemble_advection_TVD" );   
   
   // Parameters
   size_t_vector min_unknown_index(dim,0);
   size_t_vector max_unknown_index(dim,0);
   double xC = 0., yC = 0., zC = 0., xr = 0., xR = 0., xl = 0., xL = 0., 
    yt = 0., yT = 0., yb = 0., yB = 0.,
  zf = 0., zF = 0., zb = 0., zB = 0.;
   double dxC = 0., dyC = 0., dzC = 0., dxr = 0., dxl = 0., dxCr = 0., 
    dxCl = 0., dxRr = 0., dxR = 0., dxLl = 0., dyt = 0., dyb = 0., 
  dyCt = 0., dyCb = 0., dyTt = 0., dyT = 0., dyBb = 0., dzf = 0., 
  dzb = 0., dzCf = 0., dzCb = 0., dzFf = 0., dzF = 0., dzBb = 0.;

   double AdvectedValueC = 0., AdvectedValueRi = 0., AdvectedValueLe = 0.,
    AdvectedValueTo = 0., AdvectedValueBo = 0., AdvectedValueFr = 0., 
  AdvectedValueBe = 0, AdvectedValueLeLe=0., AdvectedValueRiRi=0., 
  AdvectedValueBoBo=0., AdvectedValueToTo=0., AdvectedValueBeBe=0.,
  AdvectedValueFrFr=0.,  
    AdvectorValueC = 0., AdvectorValueRi = 0., AdvectorValueLe = 0.,
    AdvectorValueTo = 0., AdvectorValueBo = 0., AdvectorValueFr = 0., 
  AdvectorValueBe = 0, AdvectorValueToLe = 0., AdvectorValueToRi = 0., 
  AdvectorValueBoLe = 0., AdvectorValueBoRi = 0., AdvectorValueFrLe = 0., 
  AdvectorValueFrRi = 0., AdvectorValueBeLe = 0., AdvectorValueBeRi = 0.,
  AdvectorValueFrTo = 0., AdvectorValueFrBo = 0., AdvectorValueBeTo = 0., 
  AdvectorValueBeBo = 0.;
   double ur = 0., ul = 0., vt = 0., vb = 0., wf = 0., wb = 0.,
  fri = 0., fle = 0., fto = 0., fbo = 0., ffr = 0., fbe = 0., flux = 0.;
   double cRip12 = 0., cLip12 = 0., cRim12 = 0., cLim12 = 0., thetaC = 0., 
    thetaRi = 0., thetaLe = 0., thetaTo = 0., thetaBo = 0., thetaFr = 0., 
  thetaBe = 0.;
   
   FV_SHIFT_TRIPLET shift = UF->shift_staggeredToStaggered( component ) ;
        
     // Perform assembling
   dxC = UF->get_cell_size( i, component, 0 ) ;    
   dyC = UF->get_cell_size( j, component, 1 ) ; 
   xC = UF->get_DOF_coordinate( i, component, 0 );
   yC = UF->get_DOF_coordinate( j, component, 1 );
 
   if ( dim == 2 )
   {
      AdvectorValueC = UF->DOF_value( i, j, k, component, 
      advecting_level );
      AdvectedValueC =UF->DOF_value( i, j, k, component, advected_level );
     
     // The First component (u)
     if ( component == 0 )
     {               
       // Right and Left
       // --------------
       if ( UF->DOF_color( i, j, k, component ) == FV_BC_RIGHT )
       {
         AdvectorValueRi = AdvectorValueC;
         AdvectedValueRi = AdvectedValueC;
       }
       else
       {
         AdvectorValueRi = UF->DOF_value( i+1, j, k, 
            component, advecting_level );      
         AdvectedValueRi =UF->DOF_value(
            i+1, j, k, component, advected_level );      
       }
       
       if ( UF->DOF_color( i, j, k, component ) == FV_BC_LEFT )
       {
         AdvectorValueLe = AdvectorValueC;
         AdvectedValueLe = AdvectedValueC;
       }
       else
       {
         AdvectorValueLe = UF->DOF_value( i-1, j, k, 
            component, advecting_level );
         AdvectedValueLe =UF->DOF_value(
            i-1, j, k, component, advected_level );
       }

       thetaC = fabs( AdvectedValueRi - AdvectedValueC ) > 1.e-20 ? 
      ( AdvectedValueC - AdvectedValueLe ) /
    ( AdvectedValueRi - AdvectedValueC ) : 1.e20;
       
       // Right (X)
       if ( UF->DOF_color( i, j, k, component ) == FV_BC_RIGHT )
         fri = AdvectorValueC * AdvectedValueC;
       else
       {       
               ur = 0.5 * ( AdvectorValueRi + AdvectorValueC );
         if ( UF->DOF_color( i+1, j, k, component ) == FV_BC_RIGHT )
         {
                 if ( ur > 0. ) fri = ur * AdvectedValueC;
                 else fri = ur * AdvectedValueRi;
         }
         else
         {
           xr =UF->get_DOF_coordinate( i+shift.i, 1, 0 );
           xR =UF->get_DOF_coordinate( i+1, component, 0 );
           dxCr = xr - xC;
           dxr  = xR - xC;
     cLip12 = AdvectedValueC + ( dxCr / dxr )
      * FV_DiscreteField::SuperBee_phi(thetaC)
        * ( AdvectedValueRi - AdvectedValueC );
      
           dxRr = xR - xr;
           dxR =UF->get_cell_size( i+1, component, 0 );
                 AdvectedValueRiRi =UF->DOF_value( 
               i+2, j, k, component, advected_level );
       
           thetaRi = fabs( AdvectedValueRiRi - AdvectedValueRi ) > 
      1.e-20 ? ( AdvectedValueRi - AdvectedValueC ) /
      ( AdvectedValueRiRi - AdvectedValueRi ) : 1.e20;
           cRip12 = AdvectedValueRi - ( dxRr / dxR ) 
      * FV_DiscreteField::SuperBee_phi(thetaRi)
      * ( AdvectedValueRiRi - AdvectedValueRi );
           fri = 0.5 * ( ur * ( cRip12 + cLip12 )
        - fabs(ur) * ( cRip12 - cLip12 ) );  
         }
             }

       // Left (X)
       if ( UF->DOF_color( i, j, k, component ) == FV_BC_LEFT )
         fle = AdvectorValueC * AdvectedValueC;
       else
       {       
         ul = 0.5 * ( AdvectorValueLe + AdvectorValueC );
         if ( UF->DOF_color(i-1, j, k, component ) == FV_BC_LEFT )
         {
           if ( ul > 0. ) fle = ul * AdvectedValueLe;
           else fle = ul * AdvectedValueC;
         }
         else
         {
           xl =UF->get_DOF_coordinate( i+shift.i-1, 1, 0 );
           xL =UF->get_DOF_coordinate( i-1, component, 0 );
           dxl  = xC - xL;
           dxLl = xl - xL;
     
                 AdvectedValueLeLe =UF->DOF_value( 
               i-2, j, k, component, advected_level );
           
     thetaLe = fabs( AdvectedValueC - AdvectedValueLe ) > 1.e-20 ?
      ( AdvectedValueLe - AdvectedValueLeLe )/
      ( AdvectedValueC - AdvectedValueLe ) : 1.e20;
           cLim12 = AdvectedValueLe + ( dxLl / dxl ) 
      * FV_DiscreteField::SuperBee_phi(thetaLe)
      * ( AdvectedValueC - AdvectedValueLe );
           if ( UF->DOF_color( i, j, k, component ) == FV_BC_RIGHT ) 
       cRim12 = AdvectedValueC;
           else
     {
             xR =UF->get_DOF_coordinate( i+1, component, 0 );
       dxr  = xR - xC;
       dxCl = xC - xl;
       
             cRim12 = AdvectedValueC - ( dxCl / dxr ) 
        * FV_DiscreteField::SuperBee_phi(thetaC)
      * ( AdvectedValueRi - AdvectedValueC );
     }
           fle = 0.5 * ( ul * ( cRim12 + cLim12 )
             - fabs(ul) * ( cRim12 - cLim12 ) );
         }
             }
   
       // Top and Bottom
       // --------------
       if ( UF->DOF_color( i, j, k, component ) == FV_BC_TOP )
       {
         AdvectorValueTo = AdvectorValueC;
         AdvectedValueTo = AdvectedValueC;
       }
       else
       {
         AdvectorValueTo = UF->DOF_value( i, j+1, k, 
            component, advecting_level );      
               AdvectedValueTo =UF->DOF_value( 
             i, j+1, k, component, advected_level );
       }

       if ( UF->DOF_color( i, j, k, component ) == FV_BC_BOTTOM )
       {
         AdvectorValueBo = AdvectorValueC;
         AdvectedValueBo = AdvectedValueC;
       }
       else
             {
         AdvectorValueBo = UF->DOF_value( i, j-1, k, 
            component, advecting_level );
               AdvectedValueBo =UF->DOF_value( 
      i, j-1, k, component, advected_level );
       }

       thetaC = fabs( AdvectedValueTo - AdvectedValueC ) > 1.e-20 ? 
    ( AdvectedValueC - AdvectedValueBo ) /
    ( AdvectedValueTo - AdvectedValueC ) : 1.e20;

       // Top (Y)
       AdvectorValueToLe = UF->DOF_value( i+shift.i-1, 
            j+shift.j, k, 1, advecting_level );
       AdvectorValueToRi = UF->DOF_value( i+shift.i, 
            j+shift.j, k, 1, advecting_level );
       vt = 0.5 * ( AdvectorValueToLe + AdvectorValueToRi );
       if ( UF->DOF_color( i, j+1, k, component ) == FV_BC_TOP
          || UF->DOF_color( i, j+1, k, component ) == FV_BC_TOP_LEFT
          || UF->DOF_color( i, j+1, k, component ) == FV_BC_TOP_RIGHT )
       {
         if ( vt > 0. ) fto = vt * AdvectedValueC;
         else fto = vt * AdvectedValueTo;
       }
       else
       {
         yt =UF->get_DOF_coordinate( j+shift.j, 1, 1 );
         yT =UF->get_DOF_coordinate( j+1, component, 1 );
         dyCt = yt - yC;
         dyt  = yT - yC;
     
         cLip12 = AdvectedValueC + ( dyCt / dyt ) 
      * FV_DiscreteField::SuperBee_phi(thetaC)
      * ( AdvectedValueTo - AdvectedValueC );
         dyTt = yT - yt;
         dyT =UF->get_cell_size( j+1, component, 1 );
           
               AdvectedValueToTo =UF->DOF_value( 
               i, j+2, k, component, advected_level );
     
         thetaTo = fabs( AdvectedValueToTo - AdvectedValueTo ) > 
      1.e-20 ? ( AdvectedValueTo - AdvectedValueC ) /
      ( AdvectedValueToTo - AdvectedValueTo ) : 1.e20;
         cRip12 = AdvectedValueTo - ( dyTt / dyT ) 
      * FV_DiscreteField::SuperBee_phi(thetaTo)
            * ( AdvectedValueToTo - AdvectedValueTo );
     
         fto = 0.5 * ( vt * ( cRip12 + cLip12 )
          - fabs(vt) * ( cRip12 - cLip12 ) );
       }

       // Bottom (Y)
       AdvectorValueBoLe = UF->DOF_value(
             i+shift.i-1, j+shift.j-1, k, 1, advecting_level );
       AdvectorValueBoRi = UF->DOF_value( 
             i+shift.i, j+shift.j-1, k, 1, advecting_level );
       vb = 0.5 * ( AdvectorValueBoLe + AdvectorValueBoRi );
       if ( UF->DOF_color( i, j-1, k, component ) == FV_BC_BOTTOM
          || UF->DOF_color( i, j-1, k, component ) == FV_BC_BOTTOM_LEFT
          || UF->DOF_color( i, j-1, k, component ) == FV_BC_BOTTOM_RIGHT )
       {
         if ( vb > 0. ) fbo = vb * AdvectedValueBo;
         else fbo = vb * AdvectedValueC;
       }
       else
       {
         yb =UF->get_DOF_coordinate( j+shift.j-1, 1, 1 );
         yB =UF->get_DOF_coordinate( j-1, component, 1 );
         dyb  = yC - yB;
         if ( UF->DOF_color( i, j, k, component ) == FV_BC_TOP)
     cRim12 = AdvectedValueC;
         else
         {
     yT =UF->get_DOF_coordinate( j+1, component, 1 );
     dyt  = yT - yC;
     dyCb = yC - yb;
           cRim12 = AdvectedValueC - ( dyCb / dyt ) 
        * FV_DiscreteField::SuperBee_phi(thetaC)
          * ( AdvectedValueTo - AdvectedValueC );
         }
         dyBb = yb - yB;
               AdvectedValueBoBo =UF->DOF_value( 
               i, j-2, k, component, advected_level );
     
         thetaBo = fabs( AdvectedValueC - AdvectedValueBo ) > 1.e-20 ?
        ( AdvectedValueBo - AdvectedValueBoBo ) /
      ( AdvectedValueC - AdvectedValueBo ) : 1.e20;
         cLim12 = AdvectedValueBo + ( dyBb / dyb ) 
      * FV_DiscreteField::SuperBee_phi(thetaBo)
            * ( AdvectedValueC - AdvectedValueBo );
         fbo = 0.5 * ( vb * ( cRim12 + cLim12 )
          - fabs(vb) * ( cRim12 - cLim12 ) );
       }
     }
     
     // The second component (v)
     else
     {
       // Right and Left
       // --------------
       if ( UF->DOF_color( i, j, k, component ) == FV_BC_RIGHT )
       {
         AdvectorValueRi = AdvectorValueC;
         AdvectedValueRi = AdvectedValueC;
       }
       else
       {
         AdvectorValueRi = UF->DOF_value( i+1, j, k, 
            component, advecting_level );      
               AdvectedValueRi =UF->DOF_value( 
             i+1, j, k, component, advected_level );
       }

       if ( UF->DOF_color( i, j, k, component ) == FV_BC_LEFT )
       {
         AdvectorValueLe = AdvectorValueC;
         AdvectedValueLe = AdvectedValueC;
       }
       else
             {
         AdvectorValueLe = UF->DOF_value( i-1, j, k, 
            component, advecting_level );
               AdvectedValueLe =UF->DOF_value( 
             i-1, j, k, component, advected_level );
       }

       thetaC = fabs( AdvectedValueRi - AdvectedValueC ) > 1.e-20 ? 
        ( AdvectedValueC - AdvectedValueLe ) /
    ( AdvectedValueRi - AdvectedValueC ) : 1.e20;

       // Right (X)
       AdvectorValueToRi = UF->DOF_value( 
             i+shift.i, j+shift.j, k, 0, advecting_level );
       AdvectorValueBoRi = UF->DOF_value( 
             i+shift.i, j+shift.j-1, k, 0, advecting_level );
       ur = 0.5 * ( AdvectorValueToRi + AdvectorValueBoRi );
       if ( UF->DOF_color( i+1, j, k, component ) == FV_BC_RIGHT
          || UF->DOF_color( i+1, j, k, component ) == FV_BC_BOTTOM_RIGHT
          || UF->DOF_color( i+1, j, k, component ) == FV_BC_TOP_RIGHT )
       {
         if ( ur > 0. ) fri = ur * AdvectedValueC;
         else fri = ur * AdvectedValueRi;
       }
       else
       {
         xr =UF->get_DOF_coordinate( i+shift.i, 0, 0 );
         xR =UF->get_DOF_coordinate( i+1, component, 0 );
         dxCr = xr - xC;
         dxr  = xR - xC;
     
         cLip12 = AdvectedValueC + ( dxCr / dxr ) 
      * FV_DiscreteField::SuperBee_phi(thetaC)
          * ( AdvectedValueRi - AdvectedValueC );
     
         dxRr = xR - xr;
         dxR =UF->get_cell_size( i+1, component, 0 );
         AdvectedValueRiRi =UF->DOF_value( 
               i+2, j, k, component, advected_level );
     
         thetaRi = fabs( AdvectedValueRiRi - AdvectedValueRi ) > 
      1.e-20 ? ( AdvectedValueRi - AdvectedValueC ) /
      ( AdvectedValueRiRi - AdvectedValueRi ) : 1.e20;
         cRip12 = AdvectedValueRi - ( dxRr / dxR ) 
      * FV_DiscreteField::SuperBee_phi(thetaRi)
            * ( AdvectedValueRiRi - AdvectedValueRi );
         fri = 0.5 * ( ur * ( cRip12 + cLip12 )
          - fabs(ur) * ( cRip12 - cLip12 ) );
       }
           
       // Left (X)
       AdvectorValueToLe = UF->DOF_value(
             i+shift.i-1, j+shift.j, k, 0, advecting_level );
       AdvectorValueBoLe = UF->DOF_value(
             i+shift.i-1, j+shift.j-1, k, 0, advecting_level );
       ul = 0.5 * ( AdvectorValueToLe + AdvectorValueBoLe );
       if ( UF->DOF_color( i-1, j, k, component ) == FV_BC_LEFT
          || UF->DOF_color( i-1, j, k, component ) == FV_BC_BOTTOM_LEFT
          || UF->DOF_color( i-1, j, k, component ) == FV_BC_TOP_LEFT)
       {
         if ( ul > 0. ) fle = ul * AdvectedValueLe;
         else fle = ul * AdvectedValueC;
       }
       else
       {
         xl =UF->get_DOF_coordinate( i+shift.i-1, 0, 0 );
         xL =UF->get_DOF_coordinate( i-1, component, 0 );
         dxl  = xC - xL;
         if ( UF->DOF_color( i, j, k, component ) == FV_BC_RIGHT ) 
     cRim12 = AdvectedValueC;
         else
         {
     xR =UF->get_DOF_coordinate( i+1, component, 0 );
     dxr  = xR - xC;
     dxCl = xC - xl;
     cRim12 = AdvectedValueC
      - ( dxCl / dxr ) 
      * FV_DiscreteField::SuperBee_phi(thetaC)
      * ( AdvectedValueRi - AdvectedValueC );
         }
         dxLl = xl - xL;
               AdvectedValueLeLe =UF->DOF_value( 
               i-2, j, k, component, advected_level );
     
         thetaLe = fabs( AdvectedValueC - AdvectedValueLe ) > 1.e-20 ?
        ( AdvectedValueLe - AdvectedValueLeLe ) /
      ( AdvectedValueC - AdvectedValueLe ) : 1.e20;
         cLim12 = AdvectedValueLe + ( dxLl / dxl ) 
      * FV_DiscreteField::SuperBee_phi(thetaLe)
            * ( AdvectedValueC - AdvectedValueLe );
         fle = 0.5 * ( ul * ( cRim12 + cLim12 )
          - fabs(ul) * ( cRim12 - cLim12 ) );
       }
   
       // Top and Bottom
       // --------------
       if ( UF->DOF_color( i, j, k, component ) == FV_BC_TOP )
       {
         AdvectorValueTo = AdvectorValueC;
         AdvectedValueTo = AdvectedValueC;
       }
       else
       {
         AdvectorValueTo = UF->DOF_value( i, j+1, k, 
            component, advecting_level );      
               AdvectedValueTo =UF->DOF_value( 
             i, j+1, k, component, advected_level );
       }

       if ( UF->DOF_color( i, j, k, component ) == FV_BC_BOTTOM )
       {
         AdvectorValueBo = AdvectorValueC;
         AdvectedValueBo = AdvectedValueC;
       }
       else
       {
               AdvectorValueBo = UF->DOF_value( i, j-1, k, 
            component, advecting_level );
               AdvectedValueBo =UF->DOF_value( 
               i, j-1, k, component, advected_level );
       }

       thetaC = fabs( AdvectedValueTo - AdvectedValueC ) > 1.e-20 ? 
      ( AdvectedValueC - AdvectedValueBo ) /
    ( AdvectedValueTo - AdvectedValueC ) : 1.e20;
       
       // Top (Y)
       if ( UF->DOF_color( i, j, k, component ) == FV_BC_TOP )
         fto = AdvectorValueC * AdvectedValueC;
       else
       {       
               vt = 0.5 * ( AdvectorValueTo + AdvectorValueC );
         if ( UF->DOF_color( i, j+1, k, component ) == FV_BC_TOP )
         {
           if ( vt > 0. ) fto = vt * AdvectedValueC;
           else fto = vt * AdvectedValueTo;
         }
         else
         {
           yt =UF->get_DOF_coordinate( j+shift.j, 0, 1 );
           yT =UF->get_DOF_coordinate( j+1, component, 1 );
           dyCt = yt - yC;
           dyt  = yT - yC;
     cLip12 = AdvectedValueC + ( dyCt / dyt ) 
      * FV_DiscreteField::SuperBee_phi(thetaC)
        * ( AdvectedValueTo - AdvectedValueC );

           dyTt = yT - yt;
           dyT =UF->get_cell_size( j+1, component, 1 );
                 AdvectedValueToTo =UF->DOF_value( 
               i, j+2, k, component, advected_level );
       
           thetaTo = fabs( AdvectedValueToTo - AdvectedValueTo ) > 
      1.e-20 ? ( AdvectedValueTo - AdvectedValueC ) /
          ( AdvectedValueToTo - AdvectedValueTo ) : 1.e20;
           cRip12 = AdvectedValueTo - ( dyTt / dyT ) 
      * FV_DiscreteField::SuperBee_phi(thetaTo)
      * ( AdvectedValueToTo - AdvectedValueTo );
           fto = 0.5 * ( vt * ( cRip12 + cLip12 )
            - fabs(vt) * ( cRip12 - cLip12 ) );
         }     
             }

       // Bottom (Y)
       if ( UF->DOF_color( i, j, k, component ) == FV_BC_BOTTOM )
         fbo = AdvectorValueC * AdvectedValueC;
       else
       {
         vb = 0.5 * ( AdvectorValueBo + AdvectorValueC );
         if ( UF->DOF_color(i,j-1,k,component) == FV_BC_BOTTOM )
         {
           if ( vb > 0. ) fbo = vb * AdvectedValueBo;
           else fbo = vb * AdvectedValueC;
         }
         else
         {
           yb =UF->get_DOF_coordinate( j+shift.j-1, 0, 1 );
           yB =UF->get_DOF_coordinate( j-1, component, 1 );
           dyb  = yC - yB;
         
           dyBb = yb - yB;
                 AdvectedValueBoBo =UF->DOF_value( 
               i, j-2, k, component, advected_level );
       
           thetaBo = fabs( AdvectedValueC - AdvectedValueBo ) > 1.e-20 ?
      ( AdvectedValueBo - AdvectedValueBoBo ) /
      ( AdvectedValueC - AdvectedValueBo ) : 1.e20;
           cLim12 = AdvectedValueBo + ( dyBb / dyb ) 
      * FV_DiscreteField::SuperBee_phi(thetaBo)
      * ( AdvectedValueC - AdvectedValueBo );
     
           if ( UF->DOF_color( i, j, k, component ) == FV_BC_TOP )
       cRim12 = AdvectedValueC;
           else
     {
       yT =UF->get_DOF_coordinate( j+1, component, 1 );
       dyt  = yT - yC;
       dyCb = yC - yb;
             cRim12 = AdvectedValueC - ( dyCb / dyt )
        * FV_DiscreteField::SuperBee_phi(thetaC)
      * ( AdvectedValueTo - AdvectedValueC );
           }
     fbo = 0.5 * ( vb * ( cRim12 + cLim12 )
            - fabs(vb) * ( cRim12 - cLim12 ) );
         }
             }
     }

           flux = ( fto - fbo ) * dxC + ( fri - fle ) * dyC;
   }
   else // DIM = 3
   {
     zC =UF->get_DOF_coordinate( k, component, 2 ) ;
       dzC =UF->get_cell_size( k, component, 2 ) ;      
             AdvectorValueC = UF->DOF_value( i, j, k, component, 
        advecting_level );
             AdvectedValueC =UF->DOF_value( 
             i, j, k, component, advected_level );

       // The First component (u)
       if ( component == 0 )
       {
         // Right and Left
         // --------------
         if ( UF->DOF_color( i, j, k, component ) == FV_BC_RIGHT )
         {
           AdvectorValueRi = AdvectorValueC;
           AdvectedValueRi = AdvectedValueC;
         }
         else
         {
           AdvectorValueRi = UF->DOF_value( i+1, j, k, 
      component, advecting_level );
                 AdvectedValueRi =UF->DOF_value( 
                i+1, j, k, component, advected_level );
         }

         if ( UF->DOF_color( i, j, k, component ) == FV_BC_LEFT )
         {
           AdvectorValueLe = AdvectorValueC;
           AdvectedValueLe = AdvectedValueC;
         }
         else
         {
           AdvectorValueLe = UF->DOF_value( i-1, j, k, 
      component, advecting_level );
                 AdvectedValueLe =UF->DOF_value( 
               i-1, j, k, component, advected_level );
        }

         thetaC = fabs( AdvectedValueRi - AdvectedValueC ) > 1.e-20 ? 
            ( AdvectedValueC - AdvectedValueLe ) /
      ( AdvectedValueRi - AdvectedValueC ) : 1.e20;
       
         // Right (X)
         if ( UF->DOF_color( i, j, k, component ) == FV_BC_RIGHT )
           fri = AdvectorValueC * AdvectedValueC;
         else
         {       
           ur = 0.5 * ( AdvectorValueRi + AdvectorValueC );
           if ( UF->DOF_color( i+1, j, k, component ) == FV_BC_RIGHT )
     {
                   if ( ur > 0. ) fri = ur * AdvectedValueC;
                   else fri = ur * AdvectedValueRi;
     }
           else
           {     
       xr =UF->get_DOF_coordinate( i+shift.i, 1, 0 );
             xR =UF->get_DOF_coordinate( i+1, component, 0 );
             dxCr = xr - xC;
             dxr  = xR - xC;
       cLip12 = AdvectedValueC + ( dxCr / dxr ) 
        * FV_DiscreteField::SuperBee_phi(thetaC)
        * ( AdvectedValueRi - AdvectedValueC );
      
             dxRr = xR - xr;
             dxR =UF->get_cell_size( i+1, component, 0 );
                   AdvectedValueRiRi =UF->DOF_value( 
                 i+2, j, k, component, advected_level );
       
             thetaRi = fabs( AdvectedValueRiRi - AdvectedValueRi ) > 
        1.e-20 ? ( AdvectedValueRi - AdvectedValueC ) /
      ( AdvectedValueRiRi - AdvectedValueRi ) : 1.e20 ;
             cRip12 = AdvectedValueRi - ( dxRr / dxR ) 
        * FV_DiscreteField::SuperBee_phi(thetaRi)
      * ( AdvectedValueRiRi - AdvectedValueRi );
             fri = 0.5 * ( ur * ( cRip12 + cLip12 )
            - fabs(ur) * ( cRip12 - cLip12 ) );
    }
               }

         // Left (X)
         if ( UF->DOF_color( i, j, k, component ) == FV_BC_LEFT )
           fle = AdvectorValueC * AdvectedValueC;
         else
         {       
           ul = 0.5 * ( AdvectorValueLe + AdvectorValueC );
           if ( UF->DOF_color( i-1, j, k, component ) == FV_BC_LEFT )
     {
             if ( ul > 0. ) fle = ul * AdvectedValueLe;
             else fle = ul * AdvectedValueC;
     }
           else
           {
       xl =UF->get_DOF_coordinate( i+shift.i-1, 1, 0 );
       xL =UF->get_DOF_coordinate( i-1, component, 0 );
       dxl  = xC - xL;       
       dxLl = xl - xL;
                   AdvectedValueLeLe =UF->DOF_value( 
               i-2, j, k, component, advected_level );

             thetaLe = fabs( AdvectedValueC - AdvectedValueLe ) > 1.e-20 ?
      ( AdvectedValueLe - AdvectedValueLeLe ) /
      ( AdvectedValueC - AdvectedValueLe ) : 1.e20;
             cLim12 = AdvectedValueLe + ( dxLl / dxl ) 
        * FV_DiscreteField::SuperBee_phi(thetaLe)
      * ( AdvectedValueC - AdvectedValueLe );
     
             if ( UF->DOF_color( i, j, k, component ) == FV_BC_RIGHT ) 
         cRim12 = AdvectedValueC;
             else
       {
               xR =UF->get_DOF_coordinate( i+1, component, 0 );
         dxr  = xR - xC;
         dxCl = xC - xl;
         cRim12 = AdvectedValueC - ( dxCl / dxr )
      * FV_DiscreteField::SuperBee_phi(thetaC)
      * ( AdvectedValueRi - AdvectedValueC );
             }
             fle = 0.5 * ( ul * ( cRim12 + cLim12 )
            - fabs(ul) * ( cRim12 - cLim12 ) );
           }
               }

   
         // Top and Bottom
         // --------------
         if ( UF->DOF_color( i, j, k, component ) == FV_BC_TOP )
         {
           AdvectorValueTo = AdvectorValueC;
           AdvectedValueTo = AdvectedValueC;
         }
         else
         {
           AdvectorValueTo = UF->DOF_value( i, j+1, k, 
      component, advecting_level );      
                 AdvectedValueTo =UF->DOF_value( 
                i, j+1, k, component, advected_level );
               }
         
         if ( UF->DOF_color( i, j, k, component ) == FV_BC_BOTTOM )
         {
           AdvectorValueBo = AdvectorValueC;
           AdvectedValueBo = AdvectedValueC;
         }
         else
               {
           AdvectorValueBo = UF->DOF_value( i, j-1, k, 
      component, advecting_level );
                 AdvectedValueBo =UF->DOF_value( 
                i, j-1, k, component, advected_level );
               }
         
         thetaC = fabs( AdvectedValueTo - AdvectedValueC ) > 1.e-20 ? 
      ( AdvectedValueC - AdvectedValueBo ) /
      ( AdvectedValueTo - AdvectedValueC ) : 1.e20;

         // Top (Y)
         AdvectorValueToLe = UF->DOF_value(
         i+shift.i-1, j+shift.j, k, 1, advecting_level );
         AdvectorValueToRi = UF->DOF_value( 
         i+shift.i, j+shift.j, k, 1, advecting_level );
         vt = 0.5 * ( AdvectorValueToLe + AdvectorValueToRi );
         if ( UF->DOF_color( i, j+1, k, component ) == FV_BC_TOP
              || UF->DOF_color( i, j+1, k, component ) == FV_BC_TOP_LEFT
              || UF->DOF_color( i, j+1, k, component ) 
        == FV_BC_TOP_RIGHT )
         {
           if ( vt > 0. ) fto = vt * AdvectedValueC;
           else fto = vt * AdvectedValueTo;
         }
         else
         {
           yt =UF->get_DOF_coordinate( j+shift.j, 1, 1 );
     yT =UF->get_DOF_coordinate( j+1, component, 1 );
     dyCt = yt - yC;
     dyt  = yT - yC;
     cLip12 = AdvectedValueC + ( dyCt / dyt ) 
        * FV_DiscreteField::SuperBee_phi(thetaC)
      * ( AdvectedValueTo - AdvectedValueC );
           dyTt = yT - yt;
           dyT =UF->get_cell_size( j+1, component, 1 );
                 AdvectedValueToTo =UF->DOF_value( 
                 i, j+2, k, component, advected_level );
     
           thetaTo = fabs( AdvectedValueToTo - AdvectedValueTo ) > 
        1.e-20 ? ( AdvectedValueTo - AdvectedValueC ) /
      ( AdvectedValueToTo - AdvectedValueTo ) : 1.e20;
           cRip12 = AdvectedValueTo - ( dyTt / dyT ) 
        * FV_DiscreteField::SuperBee_phi(thetaTo)
            * ( AdvectedValueToTo - AdvectedValueTo );
     fto = 0.5 * ( vt * ( cRip12 + cLip12 )
          - fabs(vt) * ( cRip12 - cLip12 ) );
         }

         // Bottom (Y)
         AdvectorValueBoLe = UF->DOF_value( 
         i+shift.i-1, j+shift.j-1, k, 1, advecting_level );
         AdvectorValueBoRi = UF->DOF_value( 
         i+shift.i, j+shift.j-1, k, 1, advecting_level );
         vb = 0.5 * ( AdvectorValueBoLe + AdvectorValueBoRi );
         if ( UF->DOF_color( i, j-1, k, component ) == FV_BC_BOTTOM
              || UF->DOF_color( i, j-1, k, component ) 
        == FV_BC_BOTTOM_LEFT
              || UF->DOF_color( i, j-1, k, component ) 
        == FV_BC_BOTTOM_RIGHT )
         {
           if ( vb > 0. ) fbo = vb * AdvectedValueBo;
           else fbo = vb * AdvectedValueC;
         }
         else
         {
           yb =UF->get_DOF_coordinate( j+shift.j-1, 1, 1 );
           yB =UF->get_DOF_coordinate( j-1, component, 1 );
           dyb  = yC - yB;
     if ( UF->DOF_color( i, j, k, component ) == FV_BC_TOP ) 
       cRim12 = AdvectedValueC;
     else
     {
       yT =UF->get_DOF_coordinate( j+1, component, 1 );
       dyt  = yT - yC;
       dyCb = yC - yb;
       cRim12 = AdvectedValueC - ( dyCb / dyt )
          * FV_DiscreteField::SuperBee_phi(thetaC)
          * ( AdvectedValueTo - AdvectedValueC );
           }
     dyBb = yb - yB;
                 AdvectedValueBoBo =UF->DOF_value( 
                 i, j-2, k, component, advected_level );
     
           thetaBo = fabs( AdvectedValueC - AdvectedValueBo ) > 1.e-20 ?
        ( AdvectedValueBo - AdvectedValueBoBo ) /
      ( AdvectedValueC - AdvectedValueBo ) : 1.e20;
           cLim12 = AdvectedValueBo + ( dyBb / dyb )
        * FV_DiscreteField::SuperBee_phi(thetaBo)
            * ( AdvectedValueC - AdvectedValueBo );
           fbo = 0.5 * ( vb * ( cRim12 + cLim12 )
          - fabs(vb) * ( cRim12 - cLim12 ) );
         }

         // Front and Behind
         // ----------------
         if ( UF->DOF_color( i, j, k, component ) == FV_BC_FRONT )
         {
           AdvectorValueFr = AdvectorValueC;
           AdvectedValueFr = AdvectedValueC;
         }
         else
         {
           AdvectorValueFr = UF->DOF_value( i, j, k+1, 
      component, advecting_level );      
                 AdvectedValueFr =UF->DOF_value( 
                i, j, k+1, component, advected_level );
               }
         
         if ( UF->DOF_color( i, j, k, component ) == FV_BC_BEHIND )
         {
           AdvectorValueBe = AdvectorValueC;
           AdvectedValueBe = AdvectedValueC;
         }
         else
               {
           AdvectorValueBe = UF->DOF_value( i, j, k-1, 
      component, advecting_level );
                 AdvectedValueBe =UF->DOF_value( 
                i, j, k-1, component, advected_level );
               }
         
         thetaC = fabs( AdvectedValueFr - AdvectedValueC ) > 1.e-20 ? 
            ( AdvectedValueC - AdvectedValueBe ) /
      ( AdvectedValueFr - AdvectedValueC ) : 1.e20;
         
         // Front (Z)
         AdvectorValueFrLe = UF->DOF_value( 
         i+shift.i-1, j, k+shift.k, 2, advecting_level );
         AdvectorValueFrRi = UF->DOF_value( 
         i+shift.i, j, k+shift.k, 2, advecting_level );
         wf = 0.5 * (AdvectorValueFrLe + AdvectorValueFrRi);
         if ( UF->DOF_color( i, j, k+1, component ) == FV_BC_FRONT
            || UF->DOF_color( i, j, k+1, component ) 
        == FV_BC_FRONT_LEFT 
            || UF->DOF_color( i, j, k+1, component ) 
        == FV_BC_FRONT_RIGHT )
         {
     if ( wf > 0. ) ffr = wf * AdvectedValueC;
           else ffr = wf * AdvectedValueFr;
         }
         else
         {
     zf =UF->get_DOF_coordinate( k+shift.k, 2, 2 );
     zF =UF->get_DOF_coordinate( k+1, component, 2 );
     dzCf = zf - zC;
     dzf  = zF - zC;
     cLip12 = AdvectedValueC + ( dzCf / dzf )
        * FV_DiscreteField::SuperBee_phi(thetaC)
      * ( AdvectedValueFr - AdvectedValueC );
     dzFf = zF - zf;
     dzF =UF->get_cell_size( k+1, component, 2 );
                 AdvectedValueFrFr =UF->DOF_value( 
                 i, j, k+2, component, advected_level );

     thetaFr = fabs( AdvectedValueFrFr - AdvectedValueFr ) > 
        1.e-20 ? ( AdvectedValueFr - AdvectedValueC ) /
      ( AdvectedValueFrFr - AdvectedValueFr ) : 1.e20;
           cRip12 = AdvectedValueFr - ( dzFf / dzF )
      * FV_DiscreteField::SuperBee_phi(thetaFr)
            * ( AdvectedValueFrFr - AdvectedValueFr );
     ffr = 0.5 * ( wf * ( cRip12 + cLip12 )
          - fabs(wf) * ( cRip12 - cLip12 ) );
         }

         // Behind (Z)
         AdvectorValueBeLe = UF->DOF_value( 
         i+shift.i-1, j, k+shift.k-1, 2, advecting_level );
         AdvectorValueBeRi = UF->DOF_value(
         i+shift.i, j, k+shift.k-1, 2, advecting_level );
         wb = 0.5 * ( AdvectorValueBeLe + AdvectorValueBeRi );
         if ( UF->DOF_color( i, j, k-1, component ) == FV_BC_BEHIND
      || UF->DOF_color( i, j, k-1, component ) 
        == FV_BC_BEHIND_LEFT
      || UF->DOF_color( i, j, k-1, component ) 
        == FV_BC_BEHIND_RIGHT )
         {
           if ( wb > 0. ) fbe = wb * AdvectedValueBe;
           else fbe = wb * AdvectedValueC;
         }
         else
         {
     zb =UF->get_DOF_coordinate( k+shift.k-1, 2, 2 );
     zB =UF->get_DOF_coordinate( k-1, component, 2 );
     dzb  = zC - zB;
     if ( UF->DOF_color( i, j, k, component ) == FV_BC_FRONT ) 
       cRim12 = AdvectedValueC;
     else
     {
       zF =UF->get_DOF_coordinate( k+1, component, 2 );
       dzf  = zF - zC;
       dzCb = zC - zb;
       cRim12 = AdvectedValueC - ( dzCb / dzf )
      * FV_DiscreteField::SuperBee_phi(thetaC)
      * ( AdvectedValueFr - AdvectedValueC );
     }
     dzBb = zb - zB;
                 AdvectedValueBeBe =UF->DOF_value( 
      i, j, k-2, component, advected_level );
     
           thetaBe = fabs( AdvectedValueC - AdvectedValueBe ) > 1.e-20 ?
        ( AdvectedValueBe - AdvectedValueBeBe ) /
      ( AdvectedValueC - AdvectedValueBe ) : 1.e20;
           cLim12 = AdvectedValueBe + ( dzBb / dzb )
      * FV_DiscreteField::SuperBee_phi(thetaBe)
            * ( AdvectedValueC - AdvectedValueBe );
     fbe = 0.5 * ( wb * ( cRim12 + cLim12 )
          - fabs(wb) * ( cRim12 - cLim12 ) );
         }
       }
       
       // The Second component (v)
       else if ( component == 1 )
       {
         // Right and Left
         // --------------
         if ( UF->DOF_color( i, j, k, component ) == FV_BC_RIGHT )
         {
           AdvectorValueRi = AdvectorValueC;
           AdvectedValueRi = AdvectedValueC;
         }
         else
         {
           AdvectorValueRi = UF->DOF_value( i+1, j, k, 
      component, advecting_level );      
                 AdvectedValueRi =UF->DOF_value( 
                i+1, j, k, component, advected_level );
               }
         
         if ( UF->DOF_color( i, j, k, component ) == FV_BC_LEFT )
         {
           AdvectorValueLe = AdvectorValueC;
           AdvectedValueLe = AdvectedValueC;
         }
         else
         {
                 AdvectorValueLe = UF->DOF_value( i-1, j, k, 
      component, advecting_level );
                 AdvectedValueLe =UF->DOF_value( 
                i-1, j, k, component, advected_level );
         }

         thetaC = fabs( AdvectedValueRi - AdvectedValueC ) > 1.e-20 ? 
            ( AdvectedValueC - AdvectedValueLe ) /
      ( AdvectedValueRi - AdvectedValueC ) : 1.e20;

         // Right (X)
         AdvectorValueToRi = UF->DOF_value( 
         i+shift.i, j+shift.j, k, 0, advecting_level );
         AdvectorValueBoRi = UF->DOF_value(
         i+shift.i, j+shift.j-1, k, 0, advecting_level );
         ur = 0.5 * ( AdvectorValueToRi + AdvectorValueBoRi );
         if ( UF->DOF_color( i+1, j, k, component ) == FV_BC_RIGHT
            || UF->DOF_color( i+1, j, k, component ) 
        == FV_BC_BOTTOM_RIGHT
            || UF->DOF_color( i+1, j, k, component ) 
        == FV_BC_TOP_RIGHT )
         {
     if ( ur > 0. ) fri = ur * AdvectedValueC;
     else fri = ur * AdvectedValueRi;
         }
         else
         {
     xr =UF->get_DOF_coordinate( i+shift.i, 0, 0 );
     xR =UF->get_DOF_coordinate( i+1, component, 0 );
     dxCr = xr - xC;
     dxr  = xR - xC;
     cLip12 = AdvectedValueC + ( dxCr / dxr ) 
        * FV_DiscreteField::SuperBee_phi(thetaC)
      * ( AdvectedValueRi - AdvectedValueC );
       
     dxRr = xR - xr;
           dxR =UF->get_cell_size( i+1, component, 0 );
                 AdvectedValueRiRi =UF->DOF_value( 
      i+2, j, k, component, advected_level );
     
           thetaRi = fabs( AdvectedValueRiRi - AdvectedValueRi ) > 
        1.e-20 ? ( AdvectedValueRi - AdvectedValueC ) /
      ( AdvectedValueRiRi - AdvectedValueRi ) : 1.e20;
           cRip12 = AdvectedValueRi - ( dxRr / dxR )
      * FV_DiscreteField::SuperBee_phi(thetaRi)
      * ( AdvectedValueRiRi - AdvectedValueRi );
     fri = 0.5 * ( ur * ( cRip12 + cLip12 )
          - fabs(ur) * ( cRip12 - cLip12 ) );
         }
         
         // Left (X)
         AdvectorValueToLe = UF->DOF_value(
         i+shift.i-1, j+shift.j, k, 0, advecting_level );
         AdvectorValueBoLe = UF->DOF_value(
         i+shift.i-1, j+shift.j-1, k, 0, advecting_level );
         ul = 0.5 * (AdvectorValueToLe + AdvectorValueBoLe);
         if ( UF->DOF_color( i-1, j, k, component ) == FV_BC_LEFT
            || UF->DOF_color( i-1, j, k, component ) 
        == FV_BC_BOTTOM_LEFT
            || UF->DOF_color( i-1, j, k, component ) 
        == FV_BC_TOP_LEFT )
         {
           if ( ul > 0. ) fle = ul * AdvectedValueLe;
           else fle = ul * AdvectedValueC;
         }
         else
         {
     xl =UF->get_DOF_coordinate( i+shift.i-1, 0, 0 );
     xL =UF->get_DOF_coordinate( i-1, component, 0 );
     dxl  = xC - xL;
     if ( UF->DOF_color( i, j, k, component ) == FV_BC_RIGHT ) 
       cRim12 = AdvectedValueC;
     else
     {
       xR =UF->get_DOF_coordinate( i+1, component, 0 );
       dxr  = xR - xC;
       dxCl = xC - xl;
       cRim12 = AdvectedValueC - ( dxCl / dxr )
      * FV_DiscreteField::SuperBee_phi(thetaC)
          * ( AdvectedValueRi - AdvectedValueC );
     }
       
     dxLl = xl - xL;
                 AdvectedValueLeLe =UF->DOF_value( 
      i-2, j, k, component, advected_level );

           thetaLe = fabs( AdvectedValueC - AdvectedValueLe ) > 1.e-20 ?
        ( AdvectedValueLe - AdvectedValueLeLe ) /
      ( AdvectedValueC - AdvectedValueLe ) : 1.e20;
           cLim12 = AdvectedValueLe + ( dxLl / dxl )
      * FV_DiscreteField::SuperBee_phi(thetaLe)
      * ( AdvectedValueC - AdvectedValueLe );
     fle = 0.5 * ( ul * ( cRim12 + cLim12 )
          - fabs(ul) * ( cRim12 - cLim12 ) );
         }


         // Top and Bottom
         // --------------
         if ( UF->DOF_color( i, j, k, component ) == FV_BC_TOP )
         {
           AdvectorValueTo = AdvectorValueC;
           AdvectedValueTo = AdvectedValueC;
         }
         else
         {
           AdvectorValueTo = UF->DOF_value( i, j+1, k, 
      component, advecting_level );      
                 AdvectedValueTo =UF->DOF_value( 
                i, j+1, k, component, advected_level );
         }

         if ( UF->DOF_color( i, j, k, component ) == FV_BC_BOTTOM )
         {
           AdvectorValueBo = AdvectorValueC;
           AdvectedValueBo = AdvectedValueC;
         }
         else
         {
                 AdvectorValueBo = UF->DOF_value( i, j-1, k, 
      component, advecting_level );
                 AdvectedValueBo =UF->DOF_value( 
                i, j-1, k, component, advected_level );
               }
         thetaC = fabs( AdvectedValueTo - AdvectedValueC ) > 1.e-20 ? 
            ( AdvectedValueC - AdvectedValueBo ) /
      ( AdvectedValueTo - AdvectedValueC ) : 1.e20;
       
         // Top (Y)
         if ( UF->DOF_color( i, j, k, component ) == FV_BC_TOP )
           fto = AdvectorValueC * AdvectedValueC;
         else
         {       
           vt = 0.5 * ( AdvectorValueTo + AdvectorValueC );
           if ( UF->DOF_color( i, j+1, k, component ) == FV_BC_TOP )
     {
             if ( vt > 0. ) fto = vt * AdvectedValueC;
             else fto = vt * AdvectedValueTo;
     }
           else
           {
       yt =UF->get_DOF_coordinate( j+shift.j, 0, 1 );
             yT =UF->get_DOF_coordinate( j+1, component, 1 );
             dyCt = yt - yC;
             dyt  = yT - yC;
       cLip12 = AdvectedValueC + ( dyCt / dyt )
        * FV_DiscreteField::SuperBee_phi(thetaC)
        * ( AdvectedValueTo - AdvectedValueC );
      
             dyTt = yT - yt;
             dyT =UF->get_cell_size( j+1, component, 1 );
                   AdvectedValueToTo =UF->DOF_value( 
                 i, j+2, k, component, advected_level );
       
             thetaTo = fabs( AdvectedValueToTo - AdvectedValueTo ) > 
        1.e-20 ? ( AdvectedValueTo - AdvectedValueC ) /
           ( AdvectedValueToTo - AdvectedValueTo ) : 1.e20;
             cRip12 = AdvectedValueTo - ( dyTt / dyT )
      * FV_DiscreteField::SuperBee_phi(thetaTo)
      * ( AdvectedValueToTo - AdvectedValueTo );
             fto = 0.5 * ( vt * ( cRip12 + cLip12 )
            - fabs(vt) * ( cRip12 - cLip12 ) );
     }
               }
   
         // Bottom (Y)
         if ( UF->DOF_color( i, j, k, component ) == FV_BC_BOTTOM )
           fbo = AdvectorValueC * AdvectedValueC;
         else
         {
           vb = 0.5 * ( AdvectorValueBo + AdvectorValueC );
     if ( UF->DOF_color(i, j-1, k, component ) == FV_BC_BOTTOM )
     {
       if ( vb > 0. ) fbo = vb * AdvectedValueBo;
             else fbo = vb * AdvectedValueC;
     }
           else
           {
             yb =UF->get_DOF_coordinate( j+shift.j-1, 0, 1 );
             yB =UF->get_DOF_coordinate( j-1, component, 1 );
             dyb  = yC - yB;
       if ( UF->DOF_color( i, j, k, component ) == FV_BC_TOP ) 
         cRim12 = AdvectedValueC;
             else
       {
               yT =UF->get_DOF_coordinate( j+1, component, 1 );
         dyt  = yT - yC;
         dyCb = yC - yb;
               cRim12 = AdvectedValueC - ( dyCb / dyt )
      * FV_DiscreteField::SuperBee_phi(thetaC)
      * ( AdvectedValueTo - AdvectedValueC );
             }
       dyBb = yb - yB;
                   AdvectedValueBoBo =UF->DOF_value( 
                 i, j-2, k, component, advected_level );
       
             thetaBo = fabs( AdvectedValueC - AdvectedValueBo ) > 1.e-20 ?
      ( AdvectedValueBo - AdvectedValueBoBo ) /
      ( AdvectedValueC - AdvectedValueBo ) : 1.e20;
             cLim12 = AdvectedValueBo + ( dyBb / dyb )
      * FV_DiscreteField::SuperBee_phi(thetaBo)
      * ( AdvectedValueC - AdvectedValueBo );
             fbo = 0.5 * ( vb * ( cRim12 + cLim12 )
            - fabs(vb) * ( cRim12 - cLim12 ) );
           }
               }


         // Front and Behind
         // ----------------
         if ( UF->DOF_color( i, j, k, component ) == FV_BC_FRONT )
         {
           AdvectorValueFr = AdvectorValueC;
           AdvectedValueFr = AdvectedValueC;
         }
         else
         {
           AdvectorValueFr = UF->DOF_value( i, j, k+1, 
      component, advecting_level );      
                 AdvectedValueFr =UF->DOF_value( 
                i, j, k+1, component, advected_level );
         }

         if ( UF->DOF_color( i, j, k, component ) == FV_BC_BEHIND )
         {
           AdvectorValueBe = AdvectorValueC;
           AdvectedValueBe = AdvectedValueC;
         }
         else
         {
                 AdvectorValueBe = UF->DOF_value( i, j, k-1, 
      component, advecting_level );
                 AdvectedValueBe =UF->DOF_value( 
                i, j, k-1, component, advected_level );
         }

         thetaC = fabs( AdvectedValueFr - AdvectedValueC ) > 1.e-20 ? 
            ( AdvectedValueC - AdvectedValueBe ) /
      ( AdvectedValueFr - AdvectedValueC ) : 1.e20;
         
         // Front (Z)
         AdvectorValueFrBo = UF->DOF_value( 
         i, j+shift.j-1, k+shift.k, 2, advecting_level );
         AdvectorValueFrTo = UF->DOF_value( 
         i, j+shift.j, k+shift.k, 2, advecting_level );
               wf = 0.5 * ( AdvectorValueFrBo + AdvectorValueFrTo );
         if ( UF->DOF_color( i, j, k+1, component ) == FV_BC_FRONT
            || UF->DOF_color( i, j, k+1, component ) 
        == FV_BC_FRONT_BOTTOM
            || UF->DOF_color( i, j, k+1, component ) 
        == FV_BC_FRONT_TOP )
         {
     if ( wf > 0. ) ffr = wf * AdvectedValueC;
           else ffr = wf * AdvectedValueFr;
         }
         else
         {
     zf =UF->get_DOF_coordinate( k+shift.k, 2, 2 );
     zF =UF->get_DOF_coordinate( k+1, component, 2 );
     dzCf = zf - zC;
     dzf  = zF - zC;
     cLip12 = AdvectedValueC + ( dzCf / dzf )
      * FV_DiscreteField::SuperBee_phi(thetaC)
      * ( AdvectedValueFr - AdvectedValueC );
     dzFf = zF - zf;
     dzF =UF->get_cell_size( k+1, component, 2 );
                 AdvectedValueFrFr =UF->DOF_value( 
      i, j, k+2, component, advected_level );

     thetaFr = fabs( AdvectedValueFrFr - AdvectedValueFr ) > 
        1.e-20 ? ( AdvectedValueFr - AdvectedValueC ) /
      ( AdvectedValueFrFr - AdvectedValueFr ) : 1.e20;
           cRip12 = AdvectedValueFr - ( dzFf / dzF )
      * FV_DiscreteField::SuperBee_phi(thetaFr)
      * ( AdvectedValueFrFr - AdvectedValueFr );
     ffr = 0.5 * ( wf * ( cRip12 + cLip12 )
          - fabs(wf) * ( cRip12 - cLip12 ) );
         }

         // Behind (Z)
         AdvectorValueBeBo = UF->DOF_value( 
         i, j+shift.j-1, k+shift.k-1, 2, advecting_level );
         AdvectorValueBeTo = UF->DOF_value(
         i, j+shift.j, k+shift.k-1, 2, advecting_level );
         wb = 0.5 * ( AdvectorValueBeBo + AdvectorValueBeTo );
         if ( UF->DOF_color( i, j, k-1, component ) == FV_BC_BEHIND
            || UF->DOF_color( i, j, k-1, component ) 
        == FV_BC_BEHIND_BOTTOM
            || UF->DOF_color( i, j, k-1, component ) 
        == FV_BC_BEHIND_TOP )
         {
           if ( wb > 0. ) fbe = wb * AdvectedValueBe;
           else fbe = wb * AdvectedValueC;
         }
         else
         {
     zb =UF->get_DOF_coordinate( k+shift.k-1, 2, 2 );
     zB =UF->get_DOF_coordinate( k-1, component, 2 );
     dzb  = zC - zB;
     if ( UF->DOF_color( i, j, k, component ) == FV_BC_FRONT ) 
       cRim12 = AdvectedValueC;
     else
     {
       zF =UF->get_DOF_coordinate( k+1, component, 2 );
       dzf  = zF - zC;
       dzCb = zC - zb;
       cRim12 = AdvectedValueC - ( dzCb / dzf )
      * FV_DiscreteField::SuperBee_phi(thetaC)
      * ( AdvectedValueFr - AdvectedValueC );
     }
     dzBb = zb - zB;
                 AdvectedValueBeBe =UF->DOF_value( 
                 i, j, k-2, component, advected_level );
     
           thetaBe = fabs( AdvectedValueC - AdvectedValueBe ) > 1.e-20 ?
        ( AdvectedValueBe - AdvectedValueBeBe ) /
      ( AdvectedValueC - AdvectedValueBe ) : 1.e20;
           cLim12 = AdvectedValueBe + ( dzBb / dzb )
      * FV_DiscreteField::SuperBee_phi(thetaBe)
      * ( AdvectedValueC - AdvectedValueBe );
     fbe = 0.5 * ( wb * ( cRim12 + cLim12 )
          - fabs(wb) * ( cRim12 - cLim12 ) );
         }
       }
       
       // The Third component (w)
       else
       {
         // Right and Left
         // --------------
         if ( UF->DOF_color( i, j, k, component ) == FV_BC_RIGHT )
         {
           AdvectorValueRi = AdvectorValueC;
           AdvectedValueRi = AdvectedValueC;
         }
         else
         {
           AdvectorValueRi = UF->DOF_value( i+1, j, k, 
      component, advecting_level );
                 AdvectedValueRi =UF->DOF_value( 
                i+1, j, k, component, advected_level );
         }

         if ( UF->DOF_color( i, j, k, component ) == FV_BC_LEFT )
         {
           AdvectorValueLe = AdvectorValueC;
           AdvectedValueLe = AdvectedValueC;
         }
         else
         {
                 AdvectorValueLe = UF->DOF_value( i-1, j, k, 
      component, advecting_level );
                 AdvectedValueLe =UF->DOF_value( 
                i-1, j, k, component, advected_level );
         }

         thetaC = fabs( AdvectedValueRi - AdvectedValueC ) > 1.e-20 ? 
            ( AdvectedValueC - AdvectedValueLe ) /
      ( AdvectedValueRi - AdvectedValueC ) : 1.e20;

         // Right (X)
         AdvectorValueFrRi = UF->DOF_value(
         i+shift.i, j, k+shift.k, 0, advecting_level );
         AdvectorValueBeRi = UF->DOF_value( 
         i+shift.i, j, k+shift.k-1, 0, advecting_level );
         ur = 0.5 * (AdvectorValueFrRi + AdvectorValueBeRi);
         if ( UF->DOF_color( i+1, j, k, component ) == FV_BC_RIGHT
            || UF->DOF_color( i+1, j, k, component ) 
        == FV_BC_BEHIND_RIGHT
            || UF->DOF_color( i+1, j, k, component ) 
        == FV_BC_FRONT_RIGHT )
         {
     if ( ur > 0. ) fri = ur * AdvectedValueC;
     else fri = ur * AdvectedValueRi;
         }
         else
         {
     xr =UF->get_DOF_coordinate( i+shift.i, 0, 0 );
     xR =UF->get_DOF_coordinate( i+1, component, 0 );
     dxCr = xr - xC;
     dxr  = xR - xC;
     cLip12 = AdvectedValueC + ( dxCr / dxr ) 
      * FV_DiscreteField::SuperBee_phi(thetaC)
      * ( AdvectedValueRi - AdvectedValueC );
       
     dxRr = xR - xr;
           dxR =UF->get_cell_size( i+1, component, 0 );
                 AdvectedValueRiRi =UF->DOF_value( 
      i+2, j, k, component, advected_level );
     
           thetaRi = fabs( AdvectedValueRiRi - AdvectedValueRi ) > 
        1.e-20 ? ( AdvectedValueRi - AdvectedValueC ) /
      ( AdvectedValueRiRi - AdvectedValueRi ) : 1.e20;
           cRip12 = AdvectedValueRi - ( dxRr / dxR )
      * FV_DiscreteField::SuperBee_phi(thetaRi)
      * ( AdvectedValueRiRi - AdvectedValueRi );
     fri = 0.5 * ( ur * ( cRip12 + cLip12 )
          - fabs(ur) * ( cRip12 - cLip12 ) );
         }
         
         // Left (X)
         AdvectorValueFrLe = UF->DOF_value(
         i+shift.i-1, j, k+shift.k, 0, advecting_level );
         AdvectorValueBeLe = UF->DOF_value(
         i+shift.i-1, j, k+shift.k-1, 0, advecting_level );
         ul = 0.5 * ( AdvectorValueFrLe + AdvectorValueBeLe );
         if ( UF->DOF_color( i-1, j, k, component ) == FV_BC_LEFT
            || UF->DOF_color( i-1, j, k, component ) 
        == FV_BC_BEHIND_LEFT
            || UF->DOF_color( i-1, j, k, component ) 
        == FV_BC_FRONT_LEFT )
         {
           if ( ul > 0. ) fle = ul * AdvectedValueLe;
           else fle = ul * AdvectedValueC;
         }
         else
         {
     xl =UF->get_DOF_coordinate( i+shift.i-1, 0, 0 );
     xL =UF->get_DOF_coordinate( i-1, component, 0 );
     dxl  = xC - xL;
     if ( UF->DOF_color( i, j, k, component ) == FV_BC_RIGHT ) 
       cRim12 = AdvectedValueC;
     else
     {
       xR =UF->get_DOF_coordinate( i+1, component, 0 );
       dxr  = xR - xC;
       dxCl = xC - xl;
       cRim12 = AdvectedValueC - ( dxCl / dxr )
      * FV_DiscreteField::SuperBee_phi(thetaC)
      * ( AdvectedValueRi - AdvectedValueC );
     }
       
     dxLl = xl - xL;
                 AdvectedValueLeLe =UF->DOF_value( 
                 i-2, j, k, component, advected_level );

           thetaLe = fabs( AdvectedValueC - AdvectedValueLe ) > 1.e-20 ?
        ( AdvectedValueLe - AdvectedValueLeLe ) /
      ( AdvectedValueC - AdvectedValueLe ) : 1.e20;
           cLim12 = AdvectedValueLe + ( dxLl / dxl )
      * FV_DiscreteField::SuperBee_phi(thetaLe)
      * ( AdvectedValueC - AdvectedValueLe );
     fle = 0.5 * ( ul * ( cRim12 + cLim12 )
          - fabs(ul) * ( cRim12 - cLim12 ) );
         }

         // Top and Bottom
         // --------------
         if ( UF->DOF_color( i, j, k, component ) == FV_BC_TOP )
         {
           AdvectorValueTo = AdvectorValueC;
           AdvectedValueTo = AdvectedValueC;
         }
         else
         {
           AdvectorValueTo = UF->DOF_value( i, j+1, k, 
      component, advecting_level );
                 AdvectedValueTo =UF->DOF_value( 
                i, j+1, k, component, advected_level );
         }

         if ( UF->DOF_color( i, j, k, component ) == FV_BC_BOTTOM )
         {
           AdvectorValueBo = AdvectorValueC;
           AdvectedValueBo = AdvectedValueC;
         }
         else
         {
                 AdvectorValueBo = UF->DOF_value( i, j-1, k, 
      component, advecting_level );
                 AdvectedValueBo =UF->DOF_value( 
                i, j-1, k, component, advected_level );
         }

         thetaC = fabs( AdvectedValueTo - AdvectedValueC ) > 1.e-20 ? 
            ( AdvectedValueC - AdvectedValueBo ) /
      ( AdvectedValueTo - AdvectedValueC ) : 1.e20;

         // Top (Y)
         AdvectorValueBeTo = UF->DOF_value(
         i, j+shift.j, k+shift.k-1, 1, advecting_level );
               AdvectorValueFrTo = UF->DOF_value( 
         i, j+shift.j, k+shift.k, 1, advecting_level );
         vt = 0.5 * ( AdvectorValueBeTo + AdvectorValueFrTo );
         if ( UF->DOF_color( i, j+1, k, component ) == FV_BC_TOP
            || UF->DOF_color( i, j+1, k, component ) 
        == FV_BC_BEHIND_TOP
            || UF->DOF_color( i, j+1, k, component ) 
        == FV_BC_FRONT_TOP )
         {
           if ( vt > 0. ) fto = vt * AdvectedValueC;
           else fto = vt * AdvectedValueTo;
         }
         else
         {
           yt =UF->get_DOF_coordinate( j+shift.j, 1, 1 );
     yT =UF->get_DOF_coordinate( j+1, component, 1 );
     dyCt = yt - yC;
     dyt  = yT - yC;
     cLip12 = AdvectedValueC + ( dyCt / dyt ) 
        * FV_DiscreteField::SuperBee_phi(thetaC)
      * ( AdvectedValueTo - AdvectedValueC );
           dyTt = yT - yt;
           dyT =UF->get_cell_size( j+1, component, 1 );
                 AdvectedValueToTo =UF->DOF_value( 
                 i, j+2, k, component, advected_level );
     
           thetaTo = fabs( AdvectedValueToTo - AdvectedValueTo ) > 
        1.e-20 ? ( AdvectedValueTo - AdvectedValueC ) /
      ( AdvectedValueToTo - AdvectedValueTo ) : 1.e20;
           cRip12 = AdvectedValueTo - ( dyTt / dyT ) 
        * FV_DiscreteField::SuperBee_phi(thetaTo)
            * ( AdvectedValueToTo - AdvectedValueTo );
     fto = 0.5 * ( vt * ( cRip12 + cLip12 )
          - fabs(vt) * ( cRip12 - cLip12 ) );
         }

         // Bottom (Y)
         AdvectorValueBeBo = UF->DOF_value( 
         i, j+shift.j-1, k+shift.k-1, 1, advecting_level );
         AdvectorValueFrBo = UF->DOF_value( 
         i, j+shift.j-1, k+shift.k, 1, advecting_level );
         vb = 0.5 * ( AdvectorValueBeBo + AdvectorValueFrBo );
         if ( UF->DOF_color( i, j-1, k, component ) == FV_BC_BOTTOM
            || UF->DOF_color( i, j-1, k, component ) 
        == FV_BC_BEHIND_BOTTOM
            || UF->DOF_color( i, j-1, k, component ) 
        == FV_BC_FRONT_BOTTOM )
         {
           if ( vb > 0. ) fbo = vb * AdvectedValueBo;
           else fbo = vb * AdvectedValueC;
         }
         else
         {
           yb =UF->get_DOF_coordinate( j+shift.j-1, 1, 1 );
           yB =UF->get_DOF_coordinate( j-1, component, 1 );
           dyb  = yC - yB;
     if ( UF->DOF_color( i, j, k, component ) == FV_BC_TOP ) 
       cRim12 = AdvectedValueC;
     else
     {
       yT =UF->get_DOF_coordinate( j+1, component, 1 );
       dyt  = yT - yC;
       dyCb = yC - yb;
       cRim12 = AdvectedValueC - ( dyCb / dyt )
      * FV_DiscreteField::SuperBee_phi(thetaC)
      * ( AdvectedValueTo - AdvectedValueC );
           }
     dyBb = yb - yB;
                 AdvectedValueBoBo =UF->DOF_value( 
                 i, j-2, k, component, advected_level );
     
           thetaBo = fabs( AdvectedValueC - AdvectedValueBo ) > 1.e-20 ?
        ( AdvectedValueBo - AdvectedValueBoBo ) /
      ( AdvectedValueC - AdvectedValueBo ) : 1.e20;
           cLim12 = AdvectedValueBo + ( dyBb / dyb )
        * FV_DiscreteField::SuperBee_phi(thetaBo)
            * ( AdvectedValueC - AdvectedValueBo );
           fbo = 0.5 * ( vb * ( cRim12 + cLim12 )
          - fabs(vb) * ( cRim12 - cLim12 ) );
         }


         // Front and Behind
         // ----------------
         if ( UF->DOF_color( i, j, k, component ) == FV_BC_FRONT )
         {
           AdvectorValueFr = AdvectorValueC;
           AdvectedValueFr = AdvectedValueC;
         }
         else
         {
           AdvectorValueFr = UF->DOF_value( i, j, k+1, 
      component, advecting_level );      
                 AdvectedValueFr =UF->DOF_value( 
                i, j, k+1, component, advected_level );
         }

         if ( UF->DOF_color( i, j, k, component ) == FV_BC_BEHIND )
         {
           AdvectorValueBe = AdvectorValueC;
           AdvectedValueBe = AdvectedValueC;
         }
         else
         {
                 AdvectorValueBe = UF->DOF_value( i, j, k-1, 
      component, advecting_level );
                 AdvectedValueBe =UF->DOF_value( 
                i, j, k-1, component, advected_level );
         }

         thetaC = fabs( AdvectedValueFr - AdvectedValueC ) > 1.e-20 ? 
            ( AdvectedValueC - AdvectedValueBe ) /
      ( AdvectedValueFr - AdvectedValueC ) : 1.e20;
         
         // Front (Z)
         if ( UF->DOF_color( i, j, k, component ) == FV_BC_FRONT )
           ffr = AdvectorValueC * AdvectedValueC;
         else
         {       
           wf = 0.5 * ( AdvectorValueFr + AdvectorValueC );
           if ( UF->DOF_color( i, j, k+1, component ) == FV_BC_FRONT )
     {
             if ( wf > 0. ) ffr = wf * AdvectedValueC;
             else ffr = wf * AdvectedValueFr;
     }
           else
           {
       zf =UF->get_DOF_coordinate( k+shift.k, 0, 2 );
             zF =UF->get_DOF_coordinate( k+1, component, 2 );
             dzCf = zf - zC;
             dzf  = zF - zC;
       cLip12 = AdvectedValueC + ( dzCf / dzf ) 
        * FV_DiscreteField::SuperBee_phi(thetaC)
        * ( AdvectedValueFr - AdvectedValueC );    

             dzFf = zF - zf;
             dzF =UF->get_cell_size( k+1, component, 2 );
                   AdvectedValueFrFr =UF->DOF_value( 
                 i, j, k+2, component, advected_level );

             thetaFr = fabs( AdvectedValueFrFr - AdvectedValueFr ) > 
        1.e-20 ? ( AdvectedValueFr - AdvectedValueC ) /
      ( AdvectedValueFrFr - AdvectedValueFr ) : 1.e20;
             cRip12 = AdvectedValueFr - ( dzFf / dzF )
      * FV_DiscreteField::SuperBee_phi(thetaFr)
      * ( AdvectedValueFrFr - AdvectedValueFr );
             ffr = 0.5 * ( wf * ( cRip12 + cLip12 )
            - fabs(wf) * ( cRip12 - cLip12 ) );
     } 
               }

         // Behind (Z)
         if ( UF->DOF_color( i, j, k, component ) == FV_BC_BEHIND )
           fbe = AdvectorValueC * AdvectedValueC;
         else
         {       
           wb = 0.5 * (AdvectorValueBe + AdvectorValueC);
           if ( UF->DOF_color(i, j, k-1, component ) == FV_BC_BEHIND )
     {
             if (wb > 0.) fbe = wb * AdvectedValueBe;
             else fbe = wb * AdvectedValueC;
     }
           else
           {
             zb =UF->get_DOF_coordinate( k+shift.k-1, 0, 2 );
             zB =UF->get_DOF_coordinate( k-1, component, 2 );
             dzb  = zC - zB;
       if ( UF->DOF_color( i, j, k, component ) == FV_BC_FRONT ) 
         cRim12 = AdvectedValueC;
             else
       {
               zF =UF->get_DOF_coordinate( k+1, component, 2 );
         dzf  = zF - zC;
         dzCb = zC - zb;
               cRim12 = AdvectedValueC - ( dzCb / dzf )
      * FV_DiscreteField::SuperBee_phi(thetaC)
      * ( AdvectedValueFr - AdvectedValueC );
       }
       dzBb = zb - zB;
                   AdvectedValueBeBe =UF->DOF_value( 
                 i, j, k-2, component, advected_level );
       
             thetaBe = fabs( AdvectedValueC - AdvectedValueBe ) > 1.e-20 ?
      ( AdvectedValueBe - AdvectedValueBeBe ) /
      ( AdvectedValueC - AdvectedValueBe ) : 1.e20;
             cLim12 = AdvectedValueBe + ( dzBb / dzb )
        * FV_DiscreteField::SuperBee_phi(thetaBe)
      * ( AdvectedValueC - AdvectedValueBe );
             fbe = 0.5 * ( wb * ( cRim12 + cLim12 )
            - fabs(wb) * ( cRim12 - cLim12 ) );
           }
               }
       }
       
             flux = ( fto - fbo ) * dxC * dzC
                  + ( fri - fle ) * dyC * dzC
                  + ( ffr - fbe ) * dxC * dyC;
   }

   return ( coef * flux ); 
}




//----------------------------------------------------------------------
void
DDS_NavierStokes:: compute_CFL( FV_TimeIterator const* t_it,
  size_t level ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_DiscreteField_Staggered:: compute_CFL" ); 
   
   double local_cfl = 0., CFL = 0., dl, velo = 0., dt = t_it->time_step() ;
   size_t component;
  
   size_t_vector min_unknown_index(dim,0);
   size_t_vector max_unknown_index(dim,0);

   // X component
   component = 0;

   // Get local max indices
   for (size_t l=0;l<dim;++l)
     min_unknown_index(l) =
      UF->get_min_index_unknown_handled_by_proc( component, l ) ;
   for (size_t l=0;l<dim;++l)
     max_unknown_index(l) =
      UF->get_max_index_unknown_handled_by_proc( component, l ) ;

   for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i)
   {          
     dl = UF->get_cell_size(i,component,0) ;      
     for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j)
     {
       if ( dim == 2 )
       {
         size_t k = 0 ;
         velo = fabs( UF->DOF_value(i,j,k,component,level));
         local_cfl = velo * dt / dl;
         CFL = local_cfl > CFL ? local_cfl : CFL; 
       }
       else
       {
         for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k)
         {
           velo = fabs( UF->DOF_value(i,j,k,component,level));
                 local_cfl = velo * dt / dl;
                 CFL = local_cfl > CFL ? local_cfl : CFL;
         }
       }
     }
   }
    
   // Y component
   component = 1;

   // Get local max indices
   for (size_t l=0;l<dim;++l)
     min_unknown_index(l) =
      UF->get_min_index_unknown_handled_by_proc( component, l ) ;
   for (size_t l=0;l<dim;++l)
     max_unknown_index(l) =
      UF->get_max_index_unknown_handled_by_proc( component, l ) ;

   for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j)
   {          
     dl = UF->get_cell_size(j,component,1) ;      
     for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i)
     {
       if ( dim == 2 )
       {
         size_t k = 0 ;
   velo = fabs( UF->DOF_value(i,j,k,component,level));
         local_cfl = velo * dt / dl;
         CFL = local_cfl > CFL ? local_cfl : CFL; 
       }
       else
       {
         for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k)
   {
     velo = fabs( UF->DOF_value(i,j,k,component,level));
           local_cfl = velo * dt / dl;
           CFL = local_cfl > CFL ? local_cfl : CFL;
   }
       }
     }
   }

   if ( dim == 3 )
   {
     // Z component
     component = 2;

   // Get local max indices
   for (size_t l=0;l<dim;++l)
     min_unknown_index(l) =
      UF->get_min_index_unknown_handled_by_proc( component, l ) ;
   for (size_t l=0;l<dim;++l)
     max_unknown_index(l) =
      UF->get_max_index_unknown_handled_by_proc( component, l ) ;

     for (size_t k=min_unknown_index(2);k<=max_unknown_index(2);++k)
     {          
       dl = UF->get_cell_size(k,component,2) ;       ;
       for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i)
       {
         for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j)
   {
     velo = fabs( UF->DOF_value(i,j,k,component,level));
           local_cfl = velo * dt / dl;
           CFL = local_cfl > CFL ? local_cfl : CFL;
   }
       }
     }
   }
   
   double collective_CFL = MAC_Exec::communicator()->max(CFL);

   if(my_rank == is_master)
      MAC::out()<<"CFL: "<<collective_CFL<<endl;
   // return collective_CFL;

}
