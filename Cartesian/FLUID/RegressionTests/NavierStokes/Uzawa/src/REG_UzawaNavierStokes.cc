#include <REG_UzawaNavierStokes.hh>
#include <FV_DomainAndFields.hh>
#include <FV_DomainBuilder.hh>
#include <FV_DiscreteField.hh>
#include <FV_DiscreteField_Centered.hh>
#include <FV_SystemNumbering.hh>
#include <FV_Mesh.hh>
#include <FV.hh>
#include <FV_TimeIterator.hh>
#include <MAC.hh>
#include <MAC_Root.hh>
#include <MAC_Error.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_ModuleIterator.hh>
#include <MAC_Vector.hh>
#include <MAC_Communicator.hh>
#include <MAC_Exec.hh>
#include <MAC_Application.hh>
#include <MAC_ListIdentity.hh>
#include <intVector.hh>
#include <LA_Vector.hh>
#include <PAC_Misc.hh>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <limits>
#include <cstdlib>


REG_UzawaNavierStokes const* REG_UzawaNavierStokes::PROTOTYPE
			= new REG_UzawaNavierStokes() ;


//---------------------------------------------------------------------------
REG_UzawaNavierStokes:: REG_UzawaNavierStokes( void )
//--------------------------------------------------------------------------
   : FV_OneStepIteration( "REG_UzawaNavierStokes" )
   , PAC_ComputingTime("Solver")
{
   MAC_LABEL( "REG_UzawaNavierStokes:: REG_UzawaNavierStokes" ) ;
}




//---------------------------------------------------------------------------
REG_UzawaNavierStokes*
REG_UzawaNavierStokes:: create_replica( MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "REG_UzawaNavierStokes:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, dom, exp ) ) ;

   REG_UzawaNavierStokes* result =
                        new REG_UzawaNavierStokes( a_owner, dom, exp ) ;

   MAC_CHECK( create_replica_POST( result, a_owner, dom, exp ) ) ;
   return( result ) ;

}




//---------------------------------------------------------------------------
REG_UzawaNavierStokes:: REG_UzawaNavierStokes( MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : FV_OneStepIteration( a_owner, dom, exp )
   , PAC_ComputingTime("Solver")
   , UU ( dom->discrete_field( "velocity" ) )
   , PP ( dom->discrete_field( "pressure" ) )
   , GLOBAL_EQ( 0 )
   , density( exp->double_data( "Density") )
   , viscosity( exp->double_data( "Viscosity") )
   , imposed_CFL( 0.5 )
   , AdvectionScheme( "TVD" )
   , resultsDirectory( "Res" )
   , ViscousTimeAccuracy( 1 )
   , AdvectionTimeAccuracy( 1 )
   , b_PreconditionedWithLapP( false )
   , b_restart( false )
   , restartFileDirectory( MAC::undefined_string )
   , ugunm2_restartFilename_Prefix( "UGUNM2_" )
   , b_pressure_rescaling( false )
{
   MAC_LABEL( "REG_UzawaNavierStokes:: REG_UzawaNavierStokes" ) ;
   MAC_ASSERT( PP->discretization_type() == "centered" ) ;
   MAC_ASSERT( UU->discretization_type() == "staggered" ) ;


   // Call of MAC_Communicator routine to set the rank of each proces and
   // the number of processes during execution of REG_UzawaNavierStokes
   macCOMM = MAC_Exec::communicator();
   my_rank = macCOMM->rank();
   nb_ranks = macCOMM->nb_ranks();
   is_master = 0;


   // Timing routines
   if ( my_rank == is_master )
   {
     CT_set_start();
     SCT_insert_app("Objects_Creation");
     SCT_set_start("Objects_Creation");
   }


   // Is the run a follow up of a previous job 
   b_restart = MAC_Application::is_follow();    


   // Results directory
   if( exp->has_entry( "Results_directory" ) )
     resultsDirectory = exp->string_data( "Results_directory" );


   // Clear results directory in case of a new run
   if( !b_restart ) PAC_Misc::clearAllFiles( resultsDirectory, "Savings", 
   	my_rank ) ;


   // Get space dimension
   dim = UU->primary_grid()->nb_space_dimensions() ;
   if( dim == 1 )
   {
     string error_message="Space dimension should either 2 or 3";
     MAC_Error::object()->raise_bad_data_value(
         exp, "nb_space_dimensions", error_message );
   }


   // Pressure rescaling
   b_pressure_rescaling = PP->all_BCs_nonDirichlet( 0 ) ;
   
   
   // Precondition Uzawa with Lap(p)
   if ( exp->has_entry( "PreconditionedWithLapP" ) ) 
     b_PreconditionedWithLapP = exp->bool_data( "PreconditionedWithLapP" );
   

   // Imposed CFL
   if ( exp->has_entry( "Imposed_CFL" ) )
     imposed_CFL = exp->double_data( "Imposed_CFL" );


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
   	&& UU->primary_grid()->get_security_bandwidth() < 2 )
   {
     string error_message="   >= 2 with TVD scheme";
     MAC_Error::object()->raise_bad_data_value( exp,
        "security_bandwidth", error_message );
   }   


   // Viscous term time accuracy
   if ( exp->has_entry( "ViscousTimeAccuracy" ) )
     ViscousTimeAccuracy = exp->int_data( "ViscousTimeAccuracy" );
   if ( ViscousTimeAccuracy != 1 && ViscousTimeAccuracy != 2 )
   {
     string error_message ="   - 1\n   - 2\n   ";
     MAC_Error::object()->raise_bad_data_value( exp,
         "ViscousTimeAccuracy", error_message );
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


   // NS parameters
   if ( my_rank == is_master )
   {
     MAC::out() << endl << "*** Navier & Stokes algorithm" << endl
     	<< endl;
     MAC::out() << "   Name = Uzawa" << endl;
     MAC::out() << "   Viscous term time accuracy = " <<
     	ViscousTimeAccuracy << endl;
     MAC::out() << "   Advection term time accuracy = " <<
     	AdvectionTimeAccuracy << endl;
     MAC::out() << "   Imposed CFL = " << imposed_CFL << endl;
     MAC::out() << "   Advection spatial scheme = " << AdvectionScheme
     	<< endl;
     MAC::out() << endl << endl;
   }


   // Build the matrix system
   MAC_ModuleExplorer* se =
     exp->create_subexplorer( 0, "REG_UzawaNavierStokesSystem" ) ;
   GLOBAL_EQ = REG_UzawaNavierStokesSystem::create( this, se, UU, PP,
   	ViscousTimeAccuracy, AdvectionTimeAccuracy, b_PreconditionedWithLapP,
	b_pressure_rescaling ) ;
   se->destroy() ;
   

   // Timing routines
   if ( my_rank == is_master )
   {
     SCT_insert_app("Matrix_Assembly&Initialization"); 
     SCT_insert_app("NavierStokes_problem");           
     SCT_get_elapsed_time("Objects_Creation");
   }

}




//---------------------------------------------------------------------------
REG_UzawaNavierStokes:: ~REG_UzawaNavierStokes( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "REG_UzawaNavierStokes:: ~REG_UzawaNavierStokes" ) ;

}




//---------------------------------------------------------------------------
void
REG_UzawaNavierStokes:: do_one_inner_iteration( FV_TimeIterator const* t_it ) 
//---------------------------------------------------------------------------
{
  MAC_LABEL( "REG_UzawaNavierStokes:: do_one_inner_iteration" ) ;
  MAC_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;

  start_total_timer( "REG_UzawaNavierStokes:: do_one_inner_iteration" ) ;
  start_solving_timer() ;
  macCOMM->barrier();

  // Solve the constrained problem momentum + mass conservation
  // by an Uzawa algorithm
  NavierStokes_Uzawa( t_it );

  stop_solving_timer() ;
  stop_total_timer() ;   	
   
}




//---------------------------------------------------------------------------
void
REG_UzawaNavierStokes:: do_before_time_stepping( FV_TimeIterator const* t_it, 
      	std::string const& basename )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "REG_UzawaNavierStokes:: do_before_time_stepping" ) ;
   
   start_total_timer( "REG_UzawaNavierStokes:: do_before_time_stepping" ) ;

   FV_OneStepIteration::do_before_time_stepping( t_it, basename ) ;

   if ( my_rank == is_master ) SCT_set_start("Matrix_Assembly&Initialization");
      
   // Assemble matrices  
   if ( my_rank == is_master ) 
     MAC::out() << "         Assemble matrices & vectors" << endl;   

   // Velocity unsteady matrix
   if ( my_rank == is_master ) 
     MAC::out() << "            Velocity unsteady matrix" << endl;    
   GLOBAL_EQ->assemble_velocity_unsteady_matrix( density / t_it->time_step() );	
	
   // Velocity viscous matrix and rhs
   // Note: we assemble here the total viscous matrix, in case of 2nd order 
   // Crank-Nicholson scheme for the viscous term, half viscosity for the 
   // velocity operator is taken care of at the matrix level in
   // REG_UzawaNavierStokesSystem:: finalize_constant_matrices       
   if ( my_rank == is_master ) 
     MAC::out() << "            Velocity viscous matrix & rhs" << endl;
   GLOBAL_EQ->assemble_velocity_viscous_matrix_rhs( - viscosity );	       
	       	   
   // Velocity divergence matrix and rhs   
   if ( my_rank == is_master ) 
     MAC::out() << "            Velocity divergence matrix & rhs" << endl;
   GLOBAL_EQ->assemble_pdivv_matrix_rhs( -1. ) ; 

   // Pressure laplacian matrix and rhs 
   if ( b_PreconditionedWithLapP )
   {  
     if ( my_rank == is_master ) 
       MAC::out() << "            Pressure laplacian matrix & rhs" << endl;
     GLOBAL_EQ->assemble_pressure_laplacian_matrix_rhs( -1. ) ;	
     GLOBAL_EQ->pressure_laplacian_correction();	 
   }	   
   
   // Synchronize and finalize matrices 
   GLOBAL_EQ->finalize_constant_matrices() ;

   // Initialize velocity & pressure
   GLOBAL_EQ->initialize_velocity() ;
   GLOBAL_EQ->initialize_pressure() ; 

   // Do additional reload
   do_additional_reload( basename ) ;
      
   if ( my_rank == is_master )
     SCT_get_elapsed_time( "Matrix_Assembly&Initialization" ); 
        
   stop_total_timer() ;  
      
}




//---------------------------------------------------------------------------
void
REG_UzawaNavierStokes:: do_after_time_stepping( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "REG_UzawaNavierStokes:: do_after_time_stepping" ) ;  

   // Elapsed time by sub-problems
   if ( my_rank == is_master ) 
   {
     double cputime = CT_get_elapsed_time();
     MAC::out() << endl << "Full problem" << endl;
     write_elapsed_time_smhd( MAC::out(), cputime, "Computing time");     
     SCT_get_summary( MAC::out(), cputime );
   }
     
}




//---------------------------------------------------------------------------
void
REG_UzawaNavierStokes:: do_before_inner_iterations_stage( 
	FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "REG_UzawaNavierStokes:: do_before_inner_iterations_stage" ) ;

   start_total_timer( 
   	"REG_UzawaNavierStokes:: do_before_inner_iterations_stage" ) ;

   FV_OneStepIteration::do_before_inner_iterations_stage( t_it ) ;

   // Perform matrix level operations before each time step 
   GLOBAL_EQ->at_each_time_step( );  

   stop_total_timer() ;
      
}




//---------------------------------------------------------------------------
void
REG_UzawaNavierStokes:: do_after_inner_iterations_stage( 
	FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "REG_UzawaNavierStokes:: do_after_inner_iterations_stage" ) ;
   
   start_total_timer( 
   	"REG_UzawaNavierStokes:: do_after_inner_iterations_stage" ) ;

   FV_OneStepIteration::do_after_inner_iterations_stage( t_it ) ;  

   // Compute velocity change over the time step  
   double velocity_time_change = GLOBAL_EQ->compute_velocity_change()
   	/ t_it->time_step() ;
   if ( my_rank == is_master )
     MAC::out() << "         Velocity change = " << 
     	MAC::doubleToString( ios::scientific, 14, velocity_time_change ) 
	<< endl;

   // Compute Norm(div(u))
   double normdivu = GLOBAL_EQ->compute_velocity_divergence_norm() ;
   if ( my_rank == is_master )   
     MAC::out() << "         Norm div(u) = " << 
     	MAC::doubleToString( ios::scientific, 14, normdivu ) << endl;
   
   stop_total_timer() ;
   
}
  



//---------------------------------------------------------------------------
void
REG_UzawaNavierStokes:: do_additional_savings( FV_TimeIterator const* t_it,
      	int const& cycleNumber )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "REG_UzawaNavierStokes:: do_additional_savings" ) ;

  start_total_timer( "REG_UzawaNavierStokes:: do_additional_savings" ) ;
  
  // Elapsed time by sub-problems
  if ( my_rank == is_master ) 
  {
    double cputime = CT_get_elapsed_time();
    MAC::out() << endl << "Full problem" << endl;
    write_elapsed_time_smhd( MAC::out(), cputime, "Computation time" );
    SCT_get_summary( MAC::out(), cputime );
  } 

  stop_total_timer() ;  

}




//---------------------------------------------------------------------------
void
REG_UzawaNavierStokes:: do_additional_save_for_restart( 
	FV_TimeIterator const* t_it,
      	size_t const& restartCycleNumber, std::string const& basename )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "REG_UzawaNavierStokes:: do_additional_save_for_restart" ) ;

  start_total_timer( "REG_UzawaNavierStokes:: do_additional_save_for_restart" );
  
//   This Block of code is kept commented here for archival purpose
//   Now distributed LA_Vector (of PETSC type) are saved in the standard MAC
//   restart files
//   -----------------------------------------------------------------------
//   // Extract saving directory name and root file name
//   size_t pos = basename.find_last_of( "/" );  
//   if ( restartFileDirectory == MAC::undefined_string )
//   {
//     if ( pos == std::string::npos ) restartFileDirectory = "."; 
//     else restartFileDirectory = basename.substr( 0, pos );
//   }
//   string rootfilename;
//   if ( pos == std::string::npos ) rootfilename = basename;
//   else rootfilename = basename.substr( pos + 1 );
//   pos = rootfilename.find_first_of( "." );
//   if ( pos != std::string::npos ) 
//     rootfilename = rootfilename.substr( 0, pos );
//
//   if ( restartFileDirectory == MAC::undefined_string )
//     restartFileDirectory = MAC::extract_root_directory( basename ) ;
//   string rootfilename = MAC::extract_root_file_name( basename ) ;   
// 
//   // Save vector u.grad(u)^(n-2)
//   ugunm2_restartFilename = restartFileDirectory + "/" 
//   	+ ugunm2_restartFilename_Prefix + rootfilename;
//   if ( my_rank == is_master ) 
//     MAC::out() << "         Save u.grad(u)^(n-2)" << endl;
//   GLOBAL_EQ->do_additional_savings( ugunm2_restartFilename ) ;	
//   --------------------------------------------------------------------------

  stop_total_timer() ;    
    
}




//---------------------------------------------------------------------------
void 
REG_UzawaNavierStokes::NavierStokes_Uzawa( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------	
{   
   MAC_LABEL( "REG_UzawaNavierStokes:: NavierStokes_Uzawa" ) ;

   // Compute CFL
   double computed_CFL = UU->compute_CFL( t_it, 0 );   
   size_t n_advection_subtimesteps = unsigned( computed_CFL / imposed_CFL ) + 1;
   
   if ( my_rank == is_master )
   {		       
     MAC::out() << "CFL : imposed = " << imposed_CFL << " computed = " << 
       		MAC::doubleToString( ios::scientific, 6, computed_CFL ) 
		<< "  Nb of sub time steps = " << 
		n_advection_subtimesteps << endl;        
     SCT_set_start( "NavierStokes_problem" );
     if ( n_advection_subtimesteps != 1 )
       MAC::out() << "Warning : CFL not satisfied" << endl;     
   }

   // Compute velocity advection rhs
   GLOBAL_EQ->assemble_velocity_advection( AdvectionScheme, 0, - density, 0 ) ;

   // Compute velocity advection-diffusion rhs
   GLOBAL_EQ->compute_velocityAdvectionDiffusion_rhs( 
   	b_restart, t_it->iteration_number(), true ); 

   // Uzawa saddle point problem solution
   GLOBAL_EQ->Uzawa_NavierStokes_solver( viscosity, density, 
   	t_it->time_step() ) ;

   // Copy back pressure solution in field
   PP->update_free_DOFs_value( 0, GLOBAL_EQ->get_solution_pressure() ) ;

   // Copy back velocity solution in field
   UU->update_free_DOFs_value( 0, GLOBAL_EQ->get_solution_velocity() ) ; 
      
   if ( my_rank == is_master ) SCT_get_elapsed_time("NavierStokes_problem");
   
}




//---------------------------------------------------------------------------
void 
REG_UzawaNavierStokes:: do_additional_reload( string const& basename )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "REG_UzawaNavierStokes:: do_additional_reload" ) ;

//   This Block of code is kept commented here for archival purpose
//   Now distributed LA_Vector (of PETSC type) are restored from the standard
//   MAC restart files
//   -----------------------------------------------------------------------
//    // Read values of u.grad(u)^(n-2) from restart file in case of restart 
//    // with Adams-Bashforth scheme for advection
//    if ( AdvectionTimeAccuracy == 2 )
//    {
//      // Extract saving directory name and root file name
//      string restartDir = MAC::extract_root_directory( basename ) ;
//      string rootfilename = MAC::extract_root_file_name( basename ) ;     
// 
//      // Reload vector u.grad(u)^(n-2)     
//      if ( my_rank == is_master ) 
//        MAC::out() << "            Reload u.grad(u)^(n-2)" << endl;
//      ugunm2_restartFilename = restartDir + "/" 
//   	+ ugunm2_restartFilename_Prefix + rootfilename;  	 
//      GLOBAL_EQ->initialize_ugradu_Nm2( b_restart, ugunm2_restartFilename ) ;
//    }
//   --------------------------------------------------------------------------

}




//---------------------------------------------------------------------------
void 
REG_UzawaNavierStokes:: add_storable_objects( MAC_ListIdentity* list )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "REG_UzawaNavierStokes:: add_storable_objects" ) ;

   GLOBAL_EQ->add_storable_objects( list ) ;

}
