#include <VPH_HeatTransferSystem.hh>
#include <LA_Matrix.hh>
#include <LA_Vector.hh>
#include <LA_Scatter.hh>
#include <LA_SeqVector.hh>
#include <LA_SeqMatrix.hh>
#include <LA_Solver.hh>
#include <LA_MatrixIterator.hh>
#include <LA_StorableVectors.hh>
#include <intVector.hh>
#include <intVector.hh>
#include <MAC.hh>
#include <MAC_Exec.hh>
#include <MAC_Error.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_Timer.hh>
#include <MAC_Vector.hh>
#include <MAC_Communicator.hh>
#include <MAC_ListIdentity.hh>
#include <FV_DiscreteField.hh>
#include <FV_DiscreteField_Centered.hh>
#include <FV_SystemNumbering.hh>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>



//----------------------------------------------------------------------
VPH_HeatTransferSystem*
VPH_HeatTransferSystem:: create( MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
        FV_DiscreteField* mac_tt,
        FV_DiscreteField const* mac_uu,
        size_t const& HEAT_Diffusion_TimeAccuracy_,	
	size_t const& HEAT_Advection_TimeAccuracy_ )
//----------------------------------------------------------------------
{
   MAC_LABEL( "VPH_HeatTransferSystem:: create" ) ;
   MAC_CHECK_PRE( exp != 0 ) ;
   MAC_CHECK_PRE( mac_uu != 0 ) ;
   MAC_CHECK_PRE( mac_tt != 0 ) ;     

   VPH_HeatTransferSystem* result = 
         new VPH_HeatTransferSystem( a_owner, exp, mac_tt, mac_uu,
	 	HEAT_Diffusion_TimeAccuracy_, HEAT_Advection_TimeAccuracy_ ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;	  
         
   return( result ) ;
   
}



//----------------------------------------------------------------------
VPH_HeatTransferSystem:: VPH_HeatTransferSystem(
	MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
        FV_DiscreteField* mac_tt,
        FV_DiscreteField const* mac_uu,
        size_t const& HEAT_Diffusion_TimeAccuracy_,	
	size_t const& HEAT_Advection_TimeAccuracy_ )
//----------------------------------------------------------------------
   : MAC_Object( a_owner )
   , TT( mac_tt )
   , UU( mac_uu )
   , TT_NUM( 0 )
   , VECS_Storage( 0 )
   , HEAT_Diffusion_TimeAccuracy( HEAT_Diffusion_TimeAccuracy_ )
   , HEAT_Advection_TimeAccuracy( HEAT_Advection_TimeAccuracy_ )
{
   MAC_LABEL( 
     "VPH_HeatTransferSystem:: VPH_HeatTransferSystem" ) ;

   // Build the matrices & vectors
   build_system( exp ) ;
   re_initialize() ;
   
}




//----------------------------------------------------------------------
void
VPH_HeatTransferSystem:: build_system( MAC_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "VPH_HeatTransferSystem:: build_system" ) ;

   // Temperature unsteady matrix and rhs
   MAT_A_TemperatureUnsteady = LA_Matrix::make( this,
	exp->create_subexplorer( this, "MAT_A_TemperatureUnsteady" )  ) ;
   VEC_rhs_A_Temperature = MAT_A_TemperatureUnsteady->create_vector( this ) ;

   // Temperature unsteady + diffusion matrix and rhs
   MAT_A_TemperatureUnsteadyPlusDiffusion = LA_Matrix::make( this,
	exp->create_subexplorer( this, 
	"MAT_A_TemperatureUnsteadyPlusDiffusion" ) ) ;
   VEC_rhs_A_TemperatureDiffusion = 
	MAT_A_TemperatureUnsteadyPlusDiffusion->create_vector( this ) ;  
	
   // Local vectors
   T_LOC = LA_SeqVector::create( this, 0 ) ;	
	
   // Unknowns vectors
   VEC_T = MAT_A_TemperatureUnsteady->create_vector( this ) ;
   VEC_T_previoustime = MAT_A_TemperatureUnsteady->create_vector( this ) ;
   VEC_T_timechange = MAT_A_TemperatureUnsteady->create_vector( this ) ;
	
   // Temperature advection rhs
   VEC_rhs_TemperatureAdvection = MAT_A_TemperatureUnsteady->create_vector( 
   	this ) ;
   VEC_rhs_TemperatureAdvectionDiffusion = 
   	MAT_A_TemperatureUnsteady->create_vector( this ) ; 
   VEC_rhs_TemperatureAdvection_Nm2 = MAT_A_TemperatureUnsteady->create_vector( 
   	this ) ; 
	
   // Solver
   SOLVER_A_TemperatureUnsteadyPlusDiffusion =
         LA_Solver::make( this,
               exp->create_subexplorer( this, 
	       "SOLVER_A_TemperatureUnsteadyPlusDiffusion" ) ) ;
   SOLVER_A_TemperatureUnsteadyPlusDiffusion
   	->set_initial_guess_nonzero( true );	        
   SOLVER_A_TemperatureUnsteady =
         LA_Solver::make( this,
               exp->create_subexplorer( this,
	       "SOLVER_A_TemperatureUnsteady" ) ) ; 	       

   // Velocity & pressure numbering
   TT_NUM = FV_SystemNumbering::create( this, TT ) ;

   // Vectors storage for reload
   VECS_Storage = LA_StorableVectors::create( this );	
		    
}




//----------------------------------------------------------------------
void
VPH_HeatTransferSystem:: re_initialize( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "VPH_HeatTransferSystem:: re_initialize" ) ;

   size_t t_glob = TT->nb_global_unknowns() ;
   size_t t_loc = TT->nb_local_unknowns() ;
   
   // Temperature unsteady matrix and rhs
   MAT_A_TemperatureUnsteady->re_initialize( t_glob, t_glob ) ;
   VEC_rhs_A_Temperature->re_initialize( t_glob ) ;

   // Temperature diffusion + unsteady matrix and diffusion rhs
   MAT_A_TemperatureUnsteadyPlusDiffusion->re_initialize( t_glob, t_glob );
   VEC_rhs_A_TemperatureDiffusion->re_initialize( t_glob ) ; 
   
   // Local vectors
   T_LOC->re_initialize( t_loc ) ; 
   
   // Unknowns vectors
   VEC_T->re_initialize( t_glob ) ; 
   VEC_T_previoustime->re_initialize( t_glob ) ; 
   VEC_T_timechange->re_initialize( t_glob ) ; 
   
   // Temperature advection rhs
   VEC_rhs_TemperatureAdvection->re_initialize( t_glob ) ; 
   VEC_rhs_TemperatureAdvectionDiffusion->re_initialize( t_glob ) ; 
   VEC_rhs_TemperatureAdvection_Nm2->re_initialize( t_glob ) ; 
   
   // Temperature scatter
   TT_NUM->define_scatter( VEC_T ) ;   
      
}




//----------------------------------------------------------------------
VPH_HeatTransferSystem:: ~VPH_HeatTransferSystem( void )
//----------------------------------------------------------------------
{}




//----------------------------------------------------------------------
void
VPH_HeatTransferSystem::at_each_time_step( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "VPH_HeatTransferSystem:: at_each_time_step" ) ;

   // Store temperature at previous time
   VEC_T->synchronize() ; 
   VEC_T_previoustime->set( VEC_T ) ;
   
}




//----------------------------------------------------------------------
double 
VPH_HeatTransferSystem:: compute_temperature_change( void )
//----------------------------------------------------------------------
{ 	
   MAC_LABEL( "VPH_HeatTransferSystem:: compute_temperature_change" ) ;

   VEC_T->synchronize() ;
   VEC_T_timechange->set( VEC_T ) ;
   VEC_T_timechange->sum( VEC_T_previoustime, -1.0 ) ;

   double norm_T = VEC_T->two_norm() ;
   double time_change = VEC_T_timechange->two_norm() ;
   if ( norm_T > 1e-4 ) time_change /= norm_T;   

   return ( time_change );

}




//----------------------------------------------------------------------
void
VPH_HeatTransferSystem::nullify_temperature_advection_rhs( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( 
   	"VPH_HeatTransferSystem:: nullify_temperature_advection_rhs" ) ;

   // Nullify inertia rhs
   VEC_rhs_TemperatureAdvection->nullify() ;
   
}




//----------------------------------------------------------------------
LA_SeqVector const*
VPH_HeatTransferSystem:: get_solution_temperature( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "VPH_HeatTransferSystem:: get_solution_temperature" ) ;

   TT_NUM->scatter()->get( VEC_T, T_LOC ) ;

   LA_SeqVector const* result = T_LOC ;

   return( result ) ;
   
}




//----------------------------------------------------------------------
void
VPH_HeatTransferSystem::initialize_temperature( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_NavierStokesSystem:: initialize_temperature" ) ;

   TT->extract_unknown_DOFs_value( 0, T_LOC ) ;
   TT_NUM->scatter()->set( T_LOC, VEC_T ) ;
         
}




//----------------------------------------------------------------------
void
VPH_HeatTransferSystem:: assemble_temperature_advection( 
	string const& AdvectionScheme,
      	size_t advecting_level, double const& coef, 
	size_t advected_level  )
//----------------------------------------------------------------------
{
   MAC_LABEL( 
   	"VPH_HeatTransferSystem:: assemble_temperature_advection" ) ;

   if ( AdvectionScheme == "Upwind" )
     TT->assemble_advection_Upwind( UU, advecting_level, coef,
     	advected_level, VEC_rhs_TemperatureAdvection ); 
   else
     TT->assemble_advection_TVD( UU, advecting_level, coef,
     	advected_level, VEC_rhs_TemperatureAdvection );

} 




//----------------------------------------------------------------------
void
VPH_HeatTransferSystem:: finalize_constant_matrices( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( 
   	"VPH_HeatTransferSystem:: finalize_constant_matrices" ) ;

   bool same_pattern = false ;   
   MAC_Communicator const* macCOMM = MAC_Exec::communicator();    
   
   if ( macCOMM->rank() == 0 ) 
     MAC::out() << "            Temperature unsteady + diffusion matrix" 
     	<< endl;

   // Use half thermal conductivity k such that the temperature operator is
   // MAT = rho * Cp / dt - 0.5 * k * L(T)
   // !!! Note that the BC rhs vector VEC_rhs_A_TemperatureDiffusion is the rhs
   // of the total thermal conductivity (not 0.5 * k) operator, i.e.,
   // rhs = l * rhs_minus_L(T) !!!
   if ( HEAT_Diffusion_TimeAccuracy == 2 ) 
     MAT_A_TemperatureUnsteadyPlusDiffusion->scale( 0.5 ) ;
     
   MAT_A_TemperatureUnsteadyPlusDiffusion->add_Mat( 
      	MAT_A_TemperatureUnsteady, 1., same_pattern ) ;    

   SOLVER_A_TemperatureUnsteadyPlusDiffusion->set_matrix( 
   	MAT_A_TemperatureUnsteadyPlusDiffusion ) ;
   
   SOLVER_A_TemperatureUnsteady->set_matrix( MAT_A_TemperatureUnsteady ) ;
      
}




//----------------------------------------------------------------------
void
VPH_HeatTransferSystem:: compute_temperatureAdvectionDiffusion_rhs( 
	bool const& b_restart,
	size_t const& iteration_number, 
	bool const& b_with_advection )
//----------------------------------------------------------------------
{
   MAC_LABEL( "VPH_HeatTransferSystem:: "
   	"compute_temperatureAdvectionDiffusion_rhs" ) ;     
   
   // Unsteady + diffusion term
   // First order scheme
   if( HEAT_Diffusion_TimeAccuracy == 1 )
   {
     // Compute temperature unsteady rhs = (rho*Cp/dt) * T
     MAT_A_TemperatureUnsteady->multiply_vec_then_add( VEC_T,
           VEC_rhs_TemperatureAdvectionDiffusion ) ;

     // Add diffusion BC rhs
     // rhs = (rho*Cp/dt) * T + k * rhs_minus_L(T)   
     VEC_rhs_TemperatureAdvectionDiffusion->sum( VEC_rhs_A_TemperatureDiffusion,
   	1. ) ;     			   
   }
   // Second order semi-implicit Crank-Nicholson scheme
   else    
   {
     // Compute velocity unsteady rhs
     // rhs = 2 * (rho*Cp/dt) * T
     MAT_A_TemperatureUnsteady->multiply_vec_then_add( VEC_T,
   	VEC_rhs_TemperatureAdvectionDiffusion, 2. ) ;
	
     // Compute velocity unsteady + viscous rhs
     // rhs = 2 * (rho*Cp/Dt) * T + k * rhs_minus_L(T)   
     VEC_rhs_TemperatureAdvectionDiffusion->sum( VEC_rhs_A_TemperatureDiffusion,
   	1. ) ;
	
     // Add semi-explicit diffusion part
     // rhs = 2 * (rho*Cp/dt) * T + k * rhs_L(T)
     // 		 - ( rho*Cp/dt - ( k / 2 ) * L ) * T 
     //     = (rho*Cp/dt) * T + k * rhs_minus_L(T) + ( k / 2 ) * L * T     	
     MAT_A_TemperatureUnsteadyPlusDiffusion->multiply_vec_then_add( VEC_T,
   	VEC_rhs_TemperatureAdvectionDiffusion, -1., 1. ) ;        
   }

   // Advection term
   if ( b_with_advection )
   {      
     // First order scheme
     if ( HEAT_Advection_TimeAccuracy == 1 ||
     	( iteration_number == 1 && !b_restart ) )
     {
       // Add u.gradT RHS 
       // rhs += rho * Cp * (u.gradT)^n-1
       VEC_rhs_TemperatureAdvectionDiffusion->sum( 
      	VEC_rhs_TemperatureAdvection, 1. ) ; 
     }
     // Second order explicit Adams-Bashforth scheme 
     else   
     {   
       // Add u.gradT RHS with Adams-Bashforth scheme
       // rhs += 1.5 * rho * Cp * (u.gradu)^n-1 
       //      - 0.5 * rho * Cp * (u.gradu)^(n-2)
       VEC_rhs_TemperatureAdvectionDiffusion->sum( 
     		VEC_rhs_TemperatureAdvection, 1.5 ) ;
       VEC_rhs_TemperatureAdvectionDiffusion->sum( 
     		VEC_rhs_TemperatureAdvection_Nm2, -0.5 ) ; 
     }
      	
     // Store rho * Cp * (u.gradT)^n-1 for next time step, in which this term
     // will be considered as rho * Cp * (u.gradT)^n-2
     if ( HEAT_Advection_TimeAccuracy == 2 ) 
       VEC_rhs_TemperatureAdvection_Nm2->set( VEC_rhs_TemperatureAdvection );	
   }   
   
}




//----------------------------------------------------------------------
void
VPH_HeatTransferSystem:: assemble_temperature_diffusion_matrix_rhs( 
	double const& coef_lap )
//----------------------------------------------------------------------
{
   MAC_LABEL( 
   "VPH_HeatTransferSystem:: assemble_temperature_diffusion_matrix_rhs" ) ;

   TT->assemble_constantcoef_laplacian_matrix( coef_lap,
	MAT_A_TemperatureUnsteadyPlusDiffusion, 
	VEC_rhs_A_TemperatureDiffusion );	
   
}




//----------------------------------------------------------------------
void
VPH_HeatTransferSystem:: assemble_temperature_unsteady_matrix( 
	double const& coef )
//----------------------------------------------------------------------
{
   MAC_LABEL( 
   "VPH_HeatTransferSystem:: assemble_temperature_unsteady_matrix" ) ;

   TT->assemble_mass_matrix( coef, MAT_A_TemperatureUnsteady );
   
}




//----------------------------------------------------------------------
void
VPH_HeatTransferSystem::add_storable_objects(
	MAC_ListIdentity* list )
//----------------------------------------------------------------------
{
   MAC_LABEL( "VPH_HeatTransferSystem:: add_storable_objects" ) ; 

   list->extend( VECS_Storage );
   if ( HEAT_Advection_TimeAccuracy == 2 )
   {
     VEC_rhs_TemperatureAdvection_Nm2->set_name( "UgradT_nm2" );
     VECS_Storage->add_vector_to_store( VEC_rhs_TemperatureAdvection_Nm2 ) ;
   }
	
}




//----------------------------------------------------------------------
void
VPH_HeatTransferSystem::store_ugradT_Nm2( 
	size_t const& n_advection_subtimesteps )
//----------------------------------------------------------------------
{
   MAC_LABEL( "VPH_HeatTransferSystem:: store_ugradT_Nm2" ) ; 

   VEC_rhs_TemperatureAdvection_Nm2->set( VEC_rhs_TemperatureAdvection );
   VEC_rhs_TemperatureAdvection_Nm2->scale( double(n_advection_subtimesteps) );
	
}




//----------------------------------------------------------------------
bool
VPH_HeatTransferSystem::TemperatureDiffusion_solver( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "VPH_HeatTransferSystem:: TemperatureDiffusion_solver" ) ;

   // Copy the advection-diffusion rhs to VEC_rhs_A_Temperature
   VEC_rhs_A_Temperature->set( VEC_rhs_TemperatureAdvectionDiffusion );   

   // Solve diffusion problem
   SOLVER_A_TemperatureUnsteadyPlusDiffusion->solve( VEC_rhs_A_Temperature, 
   	VEC_T ) ;
   MAC_ASSERT( 
   	SOLVER_A_TemperatureUnsteadyPlusDiffusion->solution_is_achieved() ) ;
   if ( MAC_Exec::communicator()->rank() == 0 ) 
     MAC::out() << "Nb iterations = " << 
     	SOLVER_A_TemperatureUnsteadyPlusDiffusion->nb_iterations_achieved() 
	<< endl;

   return( true ) ;

}




//----------------------------------------------------------------------
bool
VPH_HeatTransferSystem:: TemperatureAdvection_solver( void  )
//----------------------------------------------------------------------
{
   MAC_LABEL( "VPH_HeatTransferSystem:: TemperatureAdvection_solver" ) ;

   // Compute temperature unsteady rhs
   MAT_A_TemperatureUnsteady->multiply_vec_then_add(
   	VEC_T, VEC_rhs_A_Temperature ) ;

   // Compute temperature unsteady + advection rhs
   VEC_rhs_A_Temperature->sum( VEC_rhs_TemperatureAdvection );

   // Solve advection problem
   SOLVER_A_TemperatureUnsteady->solve( VEC_rhs_A_Temperature, VEC_T );   
   MAC_ASSERT( SOLVER_A_TemperatureUnsteady->solution_is_achieved() ) ;

   return( true ) ;
   
}
