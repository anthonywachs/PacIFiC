#include <REG_ViscoplasticFISTA.hh>
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


REG_ViscoplasticFISTA const* REG_ViscoplasticFISTA::PROTOTYPE
			= new REG_ViscoplasticFISTA() ;


//---------------------------------------------------------------------------
REG_ViscoplasticFISTA:: REG_ViscoplasticFISTA( void )
//--------------------------------------------------------------------------
   : FV_OneStepIteration( "REG_ViscoplasticFISTA" )
   , ComputingTime("Solver")
{
   MAC_LABEL( "REG_ViscoplasticFISTA:: REG_ViscoplasticFISTA" ) ;
}




//---------------------------------------------------------------------------
REG_ViscoplasticFISTA*
REG_ViscoplasticFISTA:: create_replica( MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "REG_ViscoplasticFISTA:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, dom, exp ) ) ;

   REG_ViscoplasticFISTA* result =
                        new REG_ViscoplasticFISTA( a_owner, dom, exp ) ;

   MAC_CHECK( create_replica_POST( result, a_owner, dom, exp ) ) ;
   return( result ) ;

}




//---------------------------------------------------------------------------
REG_ViscoplasticFISTA:: REG_ViscoplasticFISTA( MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : FV_OneStepIteration( a_owner, dom, exp )
   , ComputingTime("Solver")
   , UU( dom->discrete_field( "velocity" ) )
   , PP( dom->discrete_field( "pressure" ) )
   , DD( dom->discrete_field( "dDStensors" ) )
   , GLOBAL_EQ( 0 )
   , density( exp->double_data( "Density") )
   , viscosity( exp->double_data( "Viscosity") )
   , yield_stress( exp->double_data( "YieldStress") )
   , imposed_CFL( 0.5 )
   , AdvectionScheme( "TVD" )
   , resultsDirectory( "Res" )
   , ViscousTimeAccuracy( 1 )
   , AdvectionTimeAccuracy( 1 )
   , b_PreconditionedWithLapP( false )
   , FISTA_1overL_param( 2. * viscosity )
   , FISTA_tolerance( 1.e-4 )
   , FISTA_maxiter( 1000 )
   , d_level( 0 )
   , D_level( 1 )
   , tau_hat_level( 2 )
   , tau_level( 3 )
   , tau_old_level( 4 )
   , b_restart( false )
   , b_pressure_rescaling( false )
{
   MAC_LABEL( "REG_ViscoplasticFISTA:: REG_ViscoplasticFISTA" ) ;
   MAC_ASSERT( PP->discretization_type() == "centered" ) ;
   MAC_ASSERT( UU->discretization_type() == "staggered" ) ;
   MAC_ASSERT( DD->discretization_type() == "tensor" ) ;   
   
   // Tensor fields stored over the 5 levels of DD
   // level 0 : true strain rate tensor d
   // level 1 : strain rate tensor D(u)
   // level 2 : stress tensor tau_hat
   // level 3 : stress tensor tau 
   // level 4 : stress tensor tau at previous iteration      
   MAC_ASSERT( DD->storage_depth() == 5 );


   // Call of MAC_Communicator routine to set the rank of each proces and
   // the number of processes during execution of REG_ViscoplasticFISTA
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
   // Force to 1


   // Advection term time accuracy
   if ( exp->has_entry( "AdvectionTimeAccuracy" ) )
     AdvectionTimeAccuracy = exp->int_data( "AdvectionTimeAccuracy" );
   if ( AdvectionTimeAccuracy != 1 && AdvectionTimeAccuracy != 2 )
   {
     string error_message ="   - 1\n   - 2\n   ";
     MAC_Error::object()->raise_bad_data_value( exp,
	"AdvectionTimeAccuracy", error_message );
   }


   // ALG2 numerical parameters  
   if ( exp->has_entry( "FISTA_tolerance" ) )
     FISTA_tolerance = exp->double_data( "FISTA_tolerance" ); 
   if ( exp->has_entry( "FISTA_maxiter" ) )
     FISTA_maxiter = exp->int_data( "FISTA_maxiter" );     
     

   // Viscoplastic parameters
   if ( my_rank == is_master )
   {
     MAC::out() << endl << "*** Viscoplastic algorithm" << endl
     	<< endl;
     MAC::out() << "   Name = ALG2" << endl;
     MAC::out() << "   Viscous term time accuracy = " <<
     	ViscousTimeAccuracy << endl;
     MAC::out() << "   Advection term time accuracy = " <<
     	AdvectionTimeAccuracy << endl;
     MAC::out() << "   Imposed CFL = " << imposed_CFL << endl;
     MAC::out() << "   Advection spatial scheme = " << AdvectionScheme
     	<< endl;
     MAC::out() << "   Convergence criterion = " << FISTA_tolerance << endl;
     MAC::out() << "   Maximum number of iterations = " << FISTA_maxiter 
     	<< endl;
     MAC::out() << endl << endl;
   }


   // Build the matrix system
   MAC_ModuleExplorer* se =
     exp->create_subexplorer( 0, "REG_ViscoplasticFISTASystem" ) ;
   GLOBAL_EQ = REG_ViscoplasticFISTASystem::create( this, se, UU, PP, DD,
   	ViscousTimeAccuracy, AdvectionTimeAccuracy, b_PreconditionedWithLapP,
	b_pressure_rescaling ) ;
   se->destroy() ;
   

   // Timing routines
   if ( my_rank == is_master )
   {
     SCT_insert_app("Matrix_Assembly&Initialization"); 
     SCT_insert_app("ViscoplasticFISTA_problem");           
     SCT_get_elapsed_time("Objects_Creation");
   }

}




//---------------------------------------------------------------------------
REG_ViscoplasticFISTA:: ~REG_ViscoplasticFISTA( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "REG_ViscoplasticFISTA:: ~REG_ViscoplasticFISTA" ) ;

}




//---------------------------------------------------------------------------
void
REG_ViscoplasticFISTA:: do_one_inner_iteration( FV_TimeIterator const* t_it ) 
//---------------------------------------------------------------------------
{
  MAC_LABEL( "REG_ViscoplasticFISTA:: do_one_inner_iteration" ) ;
  MAC_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;

  start_total_timer( "REG_ViscoplasticFISTA:: do_one_inner_iteration" ) ;
  start_solving_timer() ;
  macCOMM->barrier();

  // Solve the constrained problem momentum + mass conservation
  // by a FISTA algorithm
  FISTA_solver( t_it );

  stop_solving_timer() ;
  stop_total_timer() ;   	
   
}




//---------------------------------------------------------------------------
void
REG_ViscoplasticFISTA:: do_before_time_stepping( FV_TimeIterator const* t_it, 
      	std::string const& basename )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "REG_ViscoplasticFISTA:: do_before_time_stepping" ) ;
   
   start_total_timer( "REG_ViscoplasticFISTA:: do_before_time_stepping" ) ;

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
   // REG_ViscoplasticFISTASystem:: finalize_constant_matrices       
   if ( my_rank == is_master ) 
     MAC::out() << "            Velocity viscous matrix & rhs" << endl;
   GLOBAL_EQ->assemble_velocity_viscous_matrix_rhs( 
   	- 0.5 * FISTA_1overL_param );	       
	       	   
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

   // Tensoor divergence matrix
   GLOBAL_EQ->assemble_tauGradv_tensor_divergence_matrix( 1. );
      
   // Synchronize and finalize matrices 
   GLOBAL_EQ->finalize_constant_matrices() ;

   // Initialize velocity & pressure
   GLOBAL_EQ->initialize_velocity() ;
   GLOBAL_EQ->initialize_pressure() ; 
   GLOBAL_EQ->initialize_tau_hat( tau_hat_level );

   // Do additional reload
   do_additional_reload( basename ) ;
      
   if ( my_rank == is_master )
     SCT_get_elapsed_time( "Matrix_Assembly&Initialization" ); 
        
   stop_total_timer() ;  
      
}




//---------------------------------------------------------------------------
void
REG_ViscoplasticFISTA:: do_after_time_stepping( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "REG_ViscoplasticFISTA:: do_after_time_stepping" ) ;  

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
REG_ViscoplasticFISTA:: do_before_inner_iterations_stage( 
	FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "REG_ViscoplasticFISTA:: do_before_inner_iterations_stage" ) ;

   start_total_timer( 
   	"REG_ViscoplasticFISTA:: do_before_inner_iterations_stage" ) ;

   FV_OneStepIteration::do_before_inner_iterations_stage( t_it ) ;

   // Perform matrix level operations before each time step 
   GLOBAL_EQ->at_each_time_step( );  

   stop_total_timer() ;
      
}




//---------------------------------------------------------------------------
void
REG_ViscoplasticFISTA:: do_after_inner_iterations_stage( 
	FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "REG_ViscoplasticFISTA:: do_after_inner_iterations_stage" ) ;
   
   start_total_timer( 
   	"REG_ViscoplasticFISTA:: do_after_inner_iterations_stage" ) ;

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
REG_ViscoplasticFISTA:: do_additional_savings( FV_TimeIterator const* t_it,
      	int const& cycleNumber )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "REG_ViscoplasticFISTA:: do_additional_savings" ) ;

  start_total_timer( "REG_ViscoplasticFISTA:: do_additional_savings" ) ;
  
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
REG_ViscoplasticFISTA:: do_additional_save_for_restart( 
	FV_TimeIterator const* t_it,
      	size_t const& restartCycleNumber, std::string const& basename )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "REG_ViscoplasticFISTA:: do_additional_save_for_restart" ) ;

  start_total_timer( "REG_ViscoplasticFISTA:: do_additional_save_for_restart" );

  stop_total_timer() ;    
    
}




//---------------------------------------------------------------------------
void 
REG_ViscoplasticFISTA::FISTA_solver( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------	
{   
   MAC_LABEL( "REG_ViscoplasticFISTA:: FISTA_solver" ) ;
   
   double norm_strainrate = 1., norm_velocity = 1.;
   double fista_t = 1., fista_t_old = 1.;
   size_t iter = 0;

   // Compute CFL
   double computed_CFL = UU->compute_CFL( t_it, 0 );   
   size_t n_advection_subtimesteps = unsigned( computed_CFL / imposed_CFL ) + 1;
   
   if ( my_rank == is_master )
   {		       
     MAC::out() << "CFL : imposed = " << imposed_CFL << " computed = " << 
       		MAC::doubleToString( ios::scientific, 6, computed_CFL ) 
		<< "  Nb of sub time steps = " << 
		n_advection_subtimesteps << endl;        
     SCT_set_start( "ViscoplasticFISTA_problem" );
     if ( n_advection_subtimesteps != 1 )
       MAC::out() << "Warning : CFL not satisfied" << endl;     
   }

   // Compute velocity advection rhs
   GLOBAL_EQ->assemble_velocity_advection( AdvectionScheme, 0, - density, 0 ) ;

   // Compute velocity advection-diffusion rhs
   GLOBAL_EQ->compute_velocityAdvectionDiffusion_rhs( 
   	b_restart, t_it->iteration_number(), true ); 

   // Initialize velocity vector at previous iteration
   GLOBAL_EQ->initialize_previous_iteration_velocity();	

   // ALG fixed point loop
   while( ( norm_strainrate > FISTA_tolerance 
   	|| norm_velocity > FISTA_tolerance )
	&& iter < FISTA_maxiter )
   {
     iter++;

     // Compute d from the constitutive Bingham equation
     update_strain_rate_tensor_d( GLOBAL_EQ->get_tensor_unknown_vector(),
	2. * viscosity, 0., tau_hat_level );
	
     // Compute tau_hat-(1/L).d vector at the matrix level
     // Assumes that strain rate d and tau_hat tensors values
     // were previously transferred from field level to matrix level
     // tau_hat in GLOBAL_EQ::VEC_TauhatMinus1overLd 
     // d in GLOBAL_EQ::VEC_TENSOR
     GLOBAL_EQ->compute_viscoplastic_TauhatMinus1overLd_vector( 
     	FISTA_1overL_param );     	
      
     // Navier-Stokes problem solution
     // The method Uzawa_NavierStokes_solver automatically adds 
     // div(tau_hat-(1/L).d)
     GLOBAL_EQ->Uzawa_NavierStokes_solver( 0.5 * FISTA_1overL_param, 
     	density, t_it->time_step() ) ;

     // Copy back velocity solution in field
     UU->update_free_DOFs_value( 0, GLOBAL_EQ->get_solution_velocity() ) ;    
     
     // Compute strain rate tensor D(u)
     // Values are stored in the field and in the GLOBAL_EQ::VEC_TENSOR vector
     // at the matrix level
     compute_strain_rate_tensor_D( GLOBAL_EQ->get_tensor_unknown_vector() );

     // Compute tau using d, D(u) and tau_hat
     // Values are stored in the field and in the 
     // GLOBAL_EQ::VEC_TauhatMinus1overLd vector at the matrix level     
     update_tau( GLOBAL_EQ->get_TauhatMinus1overLd_vector() );
     
     // Update fista_t parameter
     fista_t = 0.5 * ( 1. + sqrt( 1. + 4. * pow( fista_t_old, 2. ) ) );

     // Compute tau_hat using tau, fista_t and fista_t_old
     // Values are stored in the field and in the 
     // GLOBAL_EQ::VEC_TauhatMinus1overLd vector at the matrix level      
     update_tau_hat( GLOBAL_EQ->get_TauhatMinus1overLd_vector(),
     	( fista_t_old - 1. ) / fista_t );
//      update_tau_hat( GLOBAL_EQ->get_TauhatMinus1overLd_vector(),
//      	MAC::min( ( double(iter) - 1. ) / ( double(iter) - 1. + 4. ),
// 		0.5 ) );	 
     
     // Update fista_t_old
     fista_t_old = fista_t ; 
     
     // Update tau_old = copy tau on tau_old
     DD->copy_DOFs_value( tau_level, tau_old_level );     
                  
     // Compute ||D-d||
     norm_strainrate = compute_norm_Dminusd();
     
     // Compute velocity change and update previous iteration velocity vector
     norm_velocity = GLOBAL_EQ->compute_velocity_change_iteration();

     if ( my_rank == is_master ) 
     {
       MAC::out().precision(8) ;
       MAC::out() << "Viscoplastic Uzawa Iteration  = "
       		<< iter << "   CV = "
		<< norm_strainrate << " " << norm_velocity << endl;  
     }
   }	   

   // Copy back pressure solution in field
   PP->update_free_DOFs_value( 0, GLOBAL_EQ->get_solution_pressure() ) ;
   
   if ( my_rank == is_master )
     SCT_get_elapsed_time( "ViscoplasticFISTA_problem" );
   
}




//---------------------------------------------------------------------------
void 
REG_ViscoplasticFISTA:: do_additional_reload( string const& basename )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "REG_ViscoplasticFISTA:: do_additional_reload" ) ;

}




//---------------------------------------------------------------------------
void 
REG_ViscoplasticFISTA:: add_storable_objects( MAC_ListIdentity* list )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "REG_ViscoplasticFISTA:: add_storable_objects" ) ;

   GLOBAL_EQ->add_storable_objects( list ) ;

}
