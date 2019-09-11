#include <REG_ProjectionNavierStokesSystem.hh>
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
REG_ProjectionNavierStokesSystem*
REG_ProjectionNavierStokesSystem:: create( MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
	FV_DiscreteField* mac_uu,
	FV_DiscreteField* mac_pp,
        const size_t &NS_Viscous_TimeAccuracy_,	
	const size_t &NS_Advection_TimeAccuracy_,
	const bool& b_pressure_rescaling_, 
	bool const& b_ExplicitPressureGradient_,
	bool const& b_HighOrderPressureCorrection_ )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ProjectionNavierStokesSystem:: create" ) ;
   MAC_CHECK_PRE( exp != 0 ) ;
   MAC_CHECK_PRE( mac_uu != 0 ) ;
   MAC_CHECK_PRE( mac_pp != 0 ) ;     

   REG_ProjectionNavierStokesSystem* result = 
         new REG_ProjectionNavierStokesSystem( a_owner, exp, mac_uu, mac_pp,
	 	NS_Viscous_TimeAccuracy_, NS_Advection_TimeAccuracy_,
		b_pressure_rescaling_, b_ExplicitPressureGradient_,
		b_HighOrderPressureCorrection_ ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;	  
         
   return( result ) ;
   
}



//----------------------------------------------------------------------
REG_ProjectionNavierStokesSystem:: REG_ProjectionNavierStokesSystem(
	MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
	FV_DiscreteField* mac_uu,
	FV_DiscreteField* mac_pp,
        const size_t &NS_Viscous_TimeAccuracy_,	
	const size_t &NS_Advection_TimeAccuracy_,
	const bool& b_pressure_rescaling_, 
	bool const& b_ExplicitPressureGradient_,
	bool const& b_HighOrderPressureCorrection_ )
//----------------------------------------------------------------------
   : MAC_Object( a_owner )
   , UU( mac_uu )
   , PP( mac_pp )
   , UU_NUM( 0 )
   , PP_NUM( 0 )    
   , MAT_A_VelocityUnsteady( 0 )
   , VEC_rhs_A_Velocity( 0 )
   , MAT_A_VelocityUnsteadyPlusViscous( 0 )
   , VEC_rhs_A_VelocityViscous( 0 )
   , MAT_B_VelocityDivergence( 0 )
   , VEC_rhs_B_VelocityDivergence( 0 )
   , VEC_rhs_Bt_PressureGradient( 0 )
   , MAT_D_PressureLaplacian( 0 )
   , VEC_rhs_D_PressureLaplacian( 0 )
   , U_LOC( 0 )
   , P_LOC( 0 )
   , VEC_U( 0 )
   , VEC_U_previoustime( 0 )
   , VEC_U_timechange( 0 )
   , VEC_P( 0 )
   , VEC_rhs_VelocityAdvection( 0 )
   , VEC_rhs_VelocityAdvectionDiffusion( 0 )
   , VEC_rhs_VelocityAdvection_Nm2( 0 )
   , VEC_rhs_A_UnitaryPeriodicPressureGradient( 0 )
   , VEC_r( 0 )
   , VEC_w( 0 )
   , VEC_mean_pressure( 0 )
   , VEC_unit_pressure( 0 )
   , VEC_inverse_mean_pressure( 0 )
   , SOLVER_A_VelocityUnsteadyPlusViscous( 0 )
   , SOLVER_D_PressureLaplacian( 0 )
   , SOLVER_A_VelocityUnsteady( 0 )
   , VECS_Storage( 0 )
   , NS_Viscous_TimeAccuracy( NS_Viscous_TimeAccuracy_ )
   , NS_Advection_TimeAccuracy( NS_Advection_TimeAccuracy_ )
   , b_pressure_rescaling( b_pressure_rescaling_ )
   , b_ExplicitPressureGradient( b_ExplicitPressureGradient_ )
   , b_HighOrderPressureCorrection( b_HighOrderPressureCorrection_ )
{
   MAC_LABEL( 
     "REG_ProjectionNavierStokesSystem:: REG_ProjectionNavierStokesSystem" ) ;

   // Build the matrices & vectors
   build_system( exp ) ;
   re_initialize() ;
   
}




//----------------------------------------------------------------------
void
REG_ProjectionNavierStokesSystem:: build_system( MAC_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ProjectionNavierStokesSystem:: build_system" ) ;

   // Velocity unsteady matrix and rhs
   MAT_A_VelocityUnsteady = LA_Matrix::make( this,
	exp->create_subexplorer( this, "MAT_A_VelocityUnsteady" )  ) ;
   VEC_rhs_A_Velocity = MAT_A_VelocityUnsteady->create_vector( this ) ;
   
   // Velocity unsteady + viscous matrix and rhs for the Stokes problem
   MAT_A_VelocityUnsteadyPlusViscous = LA_Matrix::make( this,
	exp->create_subexplorer( this, 
	"MAT_A_VelocityUnsteadyPlusViscous" ) ) ;
   VEC_rhs_A_VelocityViscous = 
	MAT_A_VelocityUnsteady->create_vector( this ) ;   

   // Velocity divergence p.Div(v) matrix and rhs
   MAT_B_VelocityDivergence = LA_Matrix::make( this,
         exp->create_subexplorer( this, "MAT_B_VelocityDivergence" ) ) ;
   VEC_rhs_B_VelocityDivergence = MAT_B_VelocityDivergence->create_vector( 
   	this );
	
   // Pressure gradient rhs 
   VEC_rhs_Bt_PressureGradient = MAT_B_VelocityDivergence->create_vector( 
   	this ) ;	

   // Pressure Laplacian matrix and rhs
   MAT_D_PressureLaplacian = LA_Matrix::make( this,
         exp->create_subexplorer( this, "MAT_D_PressureLaplacian" ) ) ;
   VEC_rhs_D_PressureLaplacian = MAT_D_PressureLaplacian->create_vector( 
   	this ) ;
   
   // Local vectors
   U_LOC = LA_SeqVector::create( this, 0 ) ;
   P_LOC = LA_SeqVector::create( this, 0 ) ; 
   
   // Unknowns vectors
   VEC_U = MAT_A_VelocityUnsteady->create_vector( this ) ;
   VEC_U_previoustime = MAT_A_VelocityUnsteady->create_vector( this ) ;
   VEC_U_timechange = MAT_A_VelocityUnsteady->create_vector( this ) ;
   VEC_P = VEC_rhs_B_VelocityDivergence->create_vector( this ) ;

   // Velocity advection rhs
   VEC_rhs_VelocityAdvection = MAT_A_VelocityUnsteady->create_vector( this ) ;
   VEC_rhs_VelocityAdvectionDiffusion = MAT_A_VelocityUnsteady->create_vector( 
   	this ) ; 
   VEC_rhs_VelocityAdvection_Nm2 = MAT_A_VelocityUnsteady->create_vector( 
   	this ) ; 	  

   // Periodic pressure drop rhs
   VEC_rhs_A_UnitaryPeriodicPressureGradient = 
   	MAT_A_VelocityUnsteady->create_vector( this ) ; 	  
   
   // Work vectors
   // Rhs vectors         
   VEC_r = MAT_D_PressureLaplacian->create_vector( this ) ;
   VEC_w = MAT_D_PressureLaplacian->create_vector( this ) ;
   VEC_mean_pressure = MAT_D_PressureLaplacian->create_vector( this ) ;
   VEC_unit_pressure = MAT_D_PressureLaplacian->create_vector( this ) ;
   VEC_inverse_mean_pressure = MAT_D_PressureLaplacian->create_vector( this ) ;
   
   // Solver
   SOLVER_A_VelocityUnsteadyPlusViscous =
         LA_Solver::make( this,
               exp->create_subexplorer( this, 
	       "SOLVER_A_VelocityUnsteadyPlusViscous" ) ) ;
   SOLVER_A_VelocityUnsteadyPlusViscous
   	->set_initial_guess_nonzero( true );	       

   SOLVER_D_PressureLaplacian =
         LA_Solver::make( this,
               exp->create_subexplorer( this,
	       "SOLVER_D_PressureLaplacian" ) ) ; 
   SOLVER_A_VelocityUnsteady =
         LA_Solver::make( this,
               exp->create_subexplorer( this,
	       "SOLVER_A_VelocityUnsteady" ) ) ; 	       

   // Velocity & pressure numbering
   UU_NUM = FV_SystemNumbering::create( this, UU ) ;
   PP_NUM = FV_SystemNumbering::create( this, PP ) ;

   // Vectors storage for reload
   VECS_Storage = LA_StorableVectors::create( this );   
}




//----------------------------------------------------------------------
void
REG_ProjectionNavierStokesSystem:: re_initialize( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ProjectionNavierStokesSystem:: re_initialize" ) ;

   size_t u_glob = UU->nb_global_unknowns() ;
   size_t u_loc = UU->nb_local_unknowns() ;
   
   size_t p_glob = PP->nb_global_unknowns() ;
   size_t p_loc = PP->nb_local_unknowns() ; 
   
   // Velocity unsteady matrix and rhs
   MAT_A_VelocityUnsteady->re_initialize( u_glob, u_glob ) ;
   VEC_rhs_A_Velocity->re_initialize( u_glob ) ;

   // Velocity viscous + unsteady matrix and viscous rhs
   MAT_A_VelocityUnsteadyPlusViscous->re_initialize( u_glob, u_glob );
   VEC_rhs_A_VelocityViscous->re_initialize( u_glob ) ; 
   
   // Velocity divergence p.Div(v) matrix and rhs
   MAT_B_VelocityDivergence->re_initialize( p_glob, u_glob ) ;
   VEC_rhs_B_VelocityDivergence->re_initialize( p_glob ) ;

   // Pressure gradient rhs
   VEC_rhs_Bt_PressureGradient->re_initialize( u_glob ) ;
    
   // Pressure Laplacian preconditioner
   MAT_D_PressureLaplacian->re_initialize( p_glob, p_glob ) ;
   VEC_rhs_D_PressureLaplacian->re_initialize( p_glob ) ;

   // Local vectors
   U_LOC->re_initialize( u_loc ) ;
   P_LOC->re_initialize( p_loc ) ;      

   // Unknowns vectors
   VEC_U->re_initialize( u_glob ) ;
   VEC_U_previoustime->re_initialize( u_glob ) ;
   VEC_U_timechange->re_initialize( u_glob ) ;
   VEC_P->re_initialize( p_glob ) ; 
   
   // Velocity advection rhs
   VEC_rhs_VelocityAdvection->re_initialize( u_glob ) ;
   VEC_rhs_VelocityAdvectionDiffusion->re_initialize( u_glob ) ;
   VEC_rhs_VelocityAdvection_Nm2->re_initialize( u_glob ) ;   

   // Periodic pressure drop rhs
   VEC_rhs_A_UnitaryPeriodicPressureGradient->re_initialize( u_glob ) ;

   // Work vectors
   // Rhs vectors
   VEC_r->re_initialize( p_glob ) ;
   VEC_w->re_initialize( p_glob ) ;
   VEC_mean_pressure->re_initialize( p_glob ) ;
   VEC_unit_pressure->re_initialize( p_glob ) ;
   VEC_inverse_mean_pressure->re_initialize( p_glob ) ;
   
   // Velocity & pressure scatters
   UU_NUM->define_scatter( VEC_U ) ;
   PP_NUM->define_scatter( VEC_P ) ;   
   
}




//----------------------------------------------------------------------
REG_ProjectionNavierStokesSystem:: ~REG_ProjectionNavierStokesSystem( void )
//----------------------------------------------------------------------
{}




//----------------------------------------------------------------------
void
REG_ProjectionNavierStokesSystem::at_each_time_step( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ProjectionNavierStokesSystem:: at_each_time_step" ) ;

   // Store velocity at previous time
   VEC_U->synchronize() ; 
   VEC_U_previoustime->set( VEC_U ) ;
   
}




//----------------------------------------------------------------------
double 
REG_ProjectionNavierStokesSystem:: compute_velocity_change( void )
//----------------------------------------------------------------------
{ 	
   MAC_LABEL( "REG_ProjectionNavierStokesSystem:: compute_velocity_change" ) ;

   VEC_U->synchronize() ;
   VEC_U_timechange->set( VEC_U ) ;
   VEC_U_timechange->sum( VEC_U_previoustime, -1.0 ) ;

   double norm_U = VEC_U->two_norm() ;
   double time_change = VEC_U_timechange->two_norm() ;
   if ( norm_U > 1e-4 ) time_change /= norm_U;   

   return ( time_change );

}




//----------------------------------------------------------------------
double
REG_ProjectionNavierStokesSystem:: compute_velocity_divergence_norm( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( 
     "REG_ProjectionNavierStokesSystem:: compute_velocity_divergence_norm" ) ;

   MAT_B_VelocityDivergence->multiply_vec_then_add( VEC_U, VEC_r, -1.0 ) ;
   VEC_r->sum( VEC_rhs_B_VelocityDivergence ) ;

   double normdivu = VEC_r->two_norm() ;

   return normdivu;
	
}




//----------------------------------------------------------------------
void
REG_ProjectionNavierStokesSystem::nullify_velocity_advection_rhs( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( 
   	"REG_ProjectionNavierStokesSystem:: nullify_velocity_advection_rhs" ) ;

   // Nullify inertia rhs
   VEC_rhs_VelocityAdvection->nullify() ;
   
}




//----------------------------------------------------------------------
LA_SeqVector const*
REG_ProjectionNavierStokesSystem:: get_solution_velocity( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ProjectionNavierStokesSystem:: get_solution_U" ) ;

   UU_NUM->scatter()->get( VEC_U, U_LOC ) ;

   LA_SeqVector const* result = U_LOC ;

   return( result ) ;
   
}




//----------------------------------------------------------------------
LA_SeqVector const*
REG_ProjectionNavierStokesSystem:: get_solution_pressure( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ProjectionNavierStokesSystem:: get_solution_P" ) ;

   PP_NUM->scatter()->get( VEC_P, P_LOC ) ;

   LA_SeqVector const* result = P_LOC ;

   return( result ) ;
   
}




//----------------------------------------------------------------------
void
REG_ProjectionNavierStokesSystem::initialize_velocity( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_NavierStokesSystem:: initialize_velocity" ) ;

   UU->extract_unknown_DOFs_value( 0, U_LOC ) ;
   UU_NUM->scatter()->set( U_LOC, VEC_U ) ;
         
}




//----------------------------------------------------------------------
void
REG_ProjectionNavierStokesSystem::initialize_pressure( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_NavierStokesSystem:: initialize_pressure" ) ;
   
   PP->extract_unknown_DOFs_value( 0, P_LOC ) ;
   PP_NUM->scatter()->set( P_LOC, VEC_P ) ;
         
}




//----------------------------------------------------------------------
void
REG_ProjectionNavierStokesSystem:: finalize_constant_matrices( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( 
   	"REG_ProjectionNavierStokesSystem:: finalize_constant_matrices" ) ;

   bool same_pattern = false ;   
   MAC_Communicator const* macCOMM = MAC_Exec::communicator();    
   
   if ( macCOMM->rank() == 0 ) 
     MAC::out() << "            Velocity unsteady + viscous matrix" 
     	<< endl;

   // Use half viscosity such that the velocity operator is
   // MAT = rho / dt - 0.5 * mu * L(u)
   // !!! Note that the BC rhs vector VEC_rhs_A_VelocityViscous is the rhs of 
   // the total viscosity (not 0.5 * viscosity) operator, i.e.,
   // rhs = mu * rhs_minus_L(u) !!!
   if ( NS_Viscous_TimeAccuracy == 2 ) 
     MAT_A_VelocityUnsteadyPlusViscous->scale( 0.5 ) ;
     
   MAT_A_VelocityUnsteadyPlusViscous->add_Mat( 
      	MAT_A_VelocityUnsteady, 1., same_pattern ) ;    

   SOLVER_A_VelocityUnsteadyPlusViscous->set_matrix( 
   	MAT_A_VelocityUnsteadyPlusViscous ) ;
	
   SOLVER_D_PressureLaplacian->set_matrix( MAT_D_PressureLaplacian ) ;
   
   SOLVER_A_VelocityUnsteady->set_matrix( MAT_A_VelocityUnsteady ) ;
      
}




//----------------------------------------------------------------------
void
REG_ProjectionNavierStokesSystem:: compute_velocityAdvectionDiffusion_rhs( 
	bool const& b_restart,
	size_t const& iteration_number, 
	bool const& b_with_advection,
	double const& dpdl )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ProjectionNavierStokesSystem:: "
   	"compute_velocityAdvectionDiffusion_rhs" ) ;     
   
   // Unsteady + viscous term
   // First order scheme
   if( NS_Viscous_TimeAccuracy == 1 )
   {
     // Compute velocity unsteady rhs = (rho/Dt) * u
     MAT_A_VelocityUnsteady->multiply_vec_then_add( VEC_U,
           VEC_rhs_VelocityAdvectionDiffusion ) ;

     // Add viscous BC rhs
     // rhs = (rho/Dt) * u + mu * rhs_minus_L(u)   
     VEC_rhs_VelocityAdvectionDiffusion->sum( VEC_rhs_A_VelocityViscous, 
   	1. ) ;     			   
   }
   // Second order semi-implicit Crank-Nicholson scheme
   else    
   {
     // Compute velocity unsteady rhs
     // rhs = 2 * (rho/Dt) * u
     MAT_A_VelocityUnsteady->multiply_vec_then_add( VEC_U,
   	VEC_rhs_VelocityAdvectionDiffusion, 2. ) ;
	
     // Compute velocity unsteady + viscous rhs
     // rhs = 2 * (rho/Dt) * u + mu * rhs_minus_L(u)   
     VEC_rhs_VelocityAdvectionDiffusion->sum( VEC_rhs_A_VelocityViscous, 
   	1. ) ;
	
     // Add semi-explicit viscous part
     // rhs = 2 * (rho/Dt) * u + mu * rhs_L(u)
     // 		 - ( rho/Dt - ( mu / 2 ) * L ) * u 
     //     = (rho/Dt) * u + mu * rhs_minus_L(u) + ( mu / 2 ) * L * u     	
     MAT_A_VelocityUnsteadyPlusViscous->multiply_vec_then_add( VEC_U,
   	VEC_rhs_VelocityAdvectionDiffusion, -1., 1. ) ;        
   }

   // Advection term
   if ( b_with_advection )
   {      
     // First order scheme
     if ( NS_Advection_TimeAccuracy == 1 ||
     	( iteration_number == 1 && !b_restart ) )
     {
       // Add u.gradu RHS 
       // rhs += rho * (u.gradu)^n-1
       VEC_rhs_VelocityAdvectionDiffusion->sum( 
      	VEC_rhs_VelocityAdvection, 1. ) ; 
     }
     // Second order explicit Adams-Bashforth scheme 
     else   
     {   
       // Add u.gradu RHS with Adams-Bashforth scheme
       // rhs += 1.5 * rho * (u.gradu)^n-1 - 0.5 * rho * (u.gradu)^(n-2)
       VEC_rhs_VelocityAdvectionDiffusion->sum( 
     		VEC_rhs_VelocityAdvection, 1.5 ) ;
       VEC_rhs_VelocityAdvectionDiffusion->sum( 
     		VEC_rhs_VelocityAdvection_Nm2, -0.5 ) ; 
     }
      	
     // Store rho * (u.gradu)^n-1 for next time step, in which this term
     // will be considered as rho * (u.gradu)^n-2
     if ( NS_Advection_TimeAccuracy == 2 ) 
       VEC_rhs_VelocityAdvection_Nm2->set( VEC_rhs_VelocityAdvection );	
   }   

   // Add periodic pressure gradient source term
   VEC_rhs_VelocityAdvectionDiffusion->sum( 
   	VEC_rhs_A_UnitaryPeriodicPressureGradient, dpdl );
   
}




//----------------------------------------------------------------------
void
REG_ProjectionNavierStokesSystem:: assemble_velocity_viscous_matrix_rhs( 
	double const& coef_lap )
//----------------------------------------------------------------------
{
   MAC_LABEL( 
   "REG_ProjectionNavierStokesSystem:: assemble_velocity_viscous_matrix_rhs" ) ;

   UU->assemble_constantcoef_laplacian_matrix( coef_lap,
	MAT_A_VelocityUnsteadyPlusViscous, VEC_rhs_A_VelocityViscous );	
   
}




//----------------------------------------------------------------------
void
REG_ProjectionNavierStokesSystem:: assemble_velocity_unsteady_matrix( 
	double const& coef )
//----------------------------------------------------------------------
{
   MAC_LABEL( 
   "REG_ProjectionNavierStokesSystem:: assemble_velocity_unsteady_matrix" ) ;

   UU->assemble_mass_matrix( coef, MAT_A_VelocityUnsteady );
   
}




//----------------------------------------------------------------------
void
REG_ProjectionNavierStokesSystem:: assemble_pdivv_matrix_rhs( 
	double const& coef )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ProjectionNavierStokesSystem:: assemble_pdivv_matrix_rhs" ) ;

   PP->assemble_pDivv_matrix( UU, coef, MAT_B_VelocityDivergence,
   	VEC_rhs_B_VelocityDivergence );		
   
}




//----------------------------------------------------------------------
void
REG_ProjectionNavierStokesSystem:: assemble_velocity_advection( 
	string const& AdvectionScheme,
      	size_t advecting_level, double const& coef, 
	size_t advected_level  )
//----------------------------------------------------------------------
{
   MAC_LABEL( 
   	"REG_ProjectionNavierStokesSystem:: assemble_velocity_advection" ) ;

   if ( AdvectionScheme == "Upwind" )
     UU->assemble_advection_Upwind( UU, advecting_level, coef,
     	advected_level, VEC_rhs_VelocityAdvection ); 
   else
     UU->assemble_advection_TVD( UU, advecting_level, coef,
     	advected_level, VEC_rhs_VelocityAdvection );

} 




//----------------------------------------------------------------------
void
REG_ProjectionNavierStokesSystem:: assemble_pressure_laplacian_matrix_rhs( 
	double const& coef_lap )
//----------------------------------------------------------------------
{
   MAC_LABEL( 
 "REG_ProjectionNavierStokesSystem:: assemble_pressure_laplacian_matrix_rhs" ) ;

   PP->assemble_constantcoef_laplacian_matrix( coef_lap, 
   	MAT_D_PressureLaplacian, VEC_rhs_D_PressureLaplacian, 
	b_pressure_rescaling );
   
}




//----------------------------------------------------------------------
void 
REG_ProjectionNavierStokesSystem:: pressure_laplacian_correction( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( 
   	"REG_ProjectionNavierStokesSystem:: pressure_laplacian_correction" ) ;
      
   MAC_Communicator const* macCOMM = MAC_Exec::communicator();  

   if ( b_pressure_rescaling )
   {
     if ( macCOMM->rank() == 0 ) 
       MAC::out() << "            Set unit pressure vector" << endl;
     VEC_unit_pressure->set( 1.0 ) ;
     if ( macCOMM->rank() == 0 ) 
       MAC::out() << "            Assemble mean pressure vector" << endl;     
     PP->assemble_mass_vector( 1., VEC_mean_pressure );
   }
   
   if ( b_HighOrderPressureCorrection )
   {
     if ( macCOMM->rank() == 0 ) 
       MAC::out() << "            Assemble inverse mean pressure vector" 
       		<< endl;     
     PP->assemble_mass_vector( 1., VEC_inverse_mean_pressure, -1. );   
   }      

}




//----------------------------------------------------------------------
void
REG_ProjectionNavierStokesSystem::add_storable_objects(
	MAC_ListIdentity* list ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ProjectionNavierStokesSystem:: add_storable_objects" ) ; 

   list->extend( VECS_Storage );
   if ( NS_Advection_TimeAccuracy == 2 )
   {
     VEC_rhs_VelocityAdvection_Nm2->set_name( "UgradU_nm2" );
     VECS_Storage->add_vector_to_store( VEC_rhs_VelocityAdvection_Nm2 ) ;     
   }
	
}




//----------------------------------------------------------------------
bool
REG_ProjectionNavierStokesSystem::VelocityDiffusion_solver( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ProjectionNavierStokesSystem:: VelocityDiffusion_solver" ) ;

   // Copy the advection-diffusion rhs to VEC_rhs_A_Velocity
   VEC_rhs_A_Velocity->set( VEC_rhs_VelocityAdvectionDiffusion );
   
   // Add explicit pressure gradient to the rhs
   if ( b_ExplicitPressureGradient )
   {
     // Add - grad(p)  
     MAT_B_VelocityDivergence->tr_multiply_vec_then_add( VEC_P, 
   	VEC_rhs_A_Velocity, -1., 1. ) ; 	
     
     // Add grad(p) boundary conditions
     VEC_rhs_A_Velocity->sum( VEC_rhs_Bt_PressureGradient ) ;
   }

   // Solve viscous problem
   SOLVER_A_VelocityUnsteadyPlusViscous->solve( VEC_rhs_A_Velocity, 
   	VEC_U ) ;
   MAC_ASSERT( 
   	SOLVER_A_VelocityUnsteadyPlusViscous->solution_is_achieved() ) ;
   if ( MAC_Exec::communicator()->rank() == 0 ) 
     MAC::out() << "Nb iterations = " << 
     	SOLVER_A_VelocityUnsteadyPlusViscous->nb_iterations_achieved() 
	<< endl;

   return( true ) ;

}




//----------------------------------------------------------------------
double
REG_ProjectionNavierStokesSystem:: VelocityPressure_correction_solver( 
	double const& density, double const& viscosity, 
	double const& timestep )
//----------------------------------------------------------------------
{
   MAC_LABEL( 
   "REG_ProjectionNavierStokesSystem:: VelocityPressure_correction_solver" ) ;

   size_t niter_sublin = 0 ;
   double normdivu = 0. ;
   MAC_Communicator const* macCOMM = MAC_Exec::communicator();
   
   // WARNING: if b_ExplicitPressureGradient, VEC_P contains here the pressure 
   // correction (while it contains the pressure if b_ExplicitPressureGradient 
   // is false)   

   // Second step: pressure/pressure correction Poisson problem
   // ---------------------------------------------------------
   // Compute r = (density/dt)*div(u) 
   MAT_B_VelocityDivergence->multiply_vec_then_add( VEC_U, VEC_r ) ;
   VEC_r->sum( VEC_rhs_B_VelocityDivergence, -1.0 ) ; 

   normdivu = VEC_r->two_norm();
      
   VEC_r->scale( density / timestep );
     
   if ( b_pressure_rescaling ) 
   {
     if ( macCOMM->rank() == 0 ) VEC_r->set_item( 0, 0. ) ; 
     VEC_r->synchronize() ; 
   }

   // Add pressure Laplacian boundary conditions if !b_ExplicitPressureGradient
   // If we solve for pressure, Dirichlet BC on pressure are p=p_0 and we need
   // to add these terms to the rhs
   // If we solve for pressure correction, Dirichlet BC on pressure correction 
   // are pc=0 and there is nothing to add
   if ( !b_ExplicitPressureGradient ) 
     VEC_r->sum( VEC_rhs_D_PressureLaplacian ) ;
   
   // Solve pressure correction or pressure laplacian problem lap(pc)=r
   if ( b_ExplicitPressureGradient ) VEC_P->nullify();

   SOLVER_D_PressureLaplacian->solve( VEC_r, VEC_P ) ;
   MAC_ASSERT( SOLVER_D_PressureLaplacian->solution_is_achieved() ) ; 

   niter_sublin += SOLVER_D_PressureLaplacian->nb_iterations_achieved();

   if ( b_pressure_rescaling )
   {
     double smean = VEC_P->dot( VEC_mean_pressure );
     VEC_P->sum( VEC_unit_pressure, -smean ) ;
   }    


   // Third step: update velocity
   // ---------------------------
   // Compute unsteady rhs = (ro/dt)*u^n-1 
   MAT_A_VelocityUnsteady->multiply_vec_then_add( VEC_U,
	VEC_rhs_A_Velocity ) ;
     
   // Add - grad(p) or - grad(pc)
   // rhs = - grad(p or pc) + (ro/dt)*u^n-1   
   MAT_B_VelocityDivergence->tr_multiply_vec_then_add( VEC_P, 
   	VEC_rhs_A_Velocity, -1., 1. ) ;   

   // Add pressure gradient boundary conditions if !b_ExplicitPressureGradient
   if ( !b_ExplicitPressureGradient )
     VEC_rhs_A_Velocity->sum( VEC_rhs_Bt_PressureGradient ) ;

   // Solve (ro/dt)*u^n = rhs = - grad(pc) + (ro/dt)*u^n-1
   SOLVER_A_VelocityUnsteady->solve( VEC_rhs_A_Velocity, VEC_U );
   MAC_ASSERT( SOLVER_A_VelocityUnsteady->solution_is_achieved() ) ;   

   niter_sublin += SOLVER_A_VelocityUnsteady->nb_iterations_achieved();  

   if ( b_ExplicitPressureGradient && b_HighOrderPressureCorrection )
   {
     // Final step: add second order pressure correction 
     // ------------------------------------------------
     // Compute - lap(pc) = L * pc
     MAT_D_PressureLaplacian->multiply_vec_then_add( VEC_P, VEC_w ) ;
     if ( b_pressure_rescaling )
     {
       if ( macCOMM->rank() == 0 ) VEC_w->set_item( 0, 0. ) ; 
       VEC_w->synchronize() ;    
     }
     
     // We need to divide by the volume of the pressure cell 
     // as L = MAT_D_PressureLaplacian is a Finite Volume matrix and here we
     // perform point wise calculations on the pressure unknowns
     VEC_r->set_as_v_product( VEC_w, VEC_inverse_mean_pressure );   
   
     // Add to pressure correction such that 
     // pc = pc - coef * mu * dt * lap(pc) / density 
     // with coef = 1 for 1st order viscous and 0.5 for 2nd order viscous
     double pc_coef = NS_Viscous_TimeAccuracy == 1 ? 1. : 0.5;
     VEC_P->sum( VEC_r, pc_coef * viscosity * timestep / density );
   }
   
   if ( MAC_Exec::communicator()->rank() == 0 )
     MAC::out() << "Nb iterations sub-linear systems = " << 
     	niter_sublin << endl;
       	
   return( normdivu ) ;
   
}




//----------------------------------------------------------------------
LA_Vector*
REG_ProjectionNavierStokesSystem::get_pressure_DirichletBC_vector(
	void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( 
     "REG_ProjectionNavierStokesSystem:: get_pressure_DirichletBC_vector" ) ; 

   return ( VEC_rhs_Bt_PressureGradient ) ;
	
}




//----------------------------------------------------------------------
LA_Vector*
REG_ProjectionNavierStokesSystem::get_unitary_periodic_pressure_drop_vector(
	void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ProjectionNavierStokesSystem:: "
   	"get_unitary_periodic_pressure_drop_vector" ) ; 

   return ( VEC_rhs_A_UnitaryPeriodicPressureGradient ) ;
	
}




//----------------------------------------------------------------------
void
REG_ProjectionNavierStokesSystem::store_ugradu_Nm2( 
	size_t const& n_advection_subtimesteps )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ProjectionNavierStokesSystem:: store_ugradu_Nm2" ) ; 

   VEC_rhs_VelocityAdvection_Nm2->set( VEC_rhs_VelocityAdvection );
   VEC_rhs_VelocityAdvection_Nm2->scale( double(n_advection_subtimesteps) );
	
}




//----------------------------------------------------------------------
bool
REG_ProjectionNavierStokesSystem:: VelocityAdvection_solver( void  )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ProjectionNavierStokesSystem:: VelocityAdvection_solver" ) ;

   // Compute velocity unsteady rhs
   MAT_A_VelocityUnsteady->multiply_vec_then_add(
   	VEC_U, VEC_rhs_A_Velocity ) ;

   // Compute velocity unsteady + advection rhs
   VEC_rhs_A_Velocity->sum( VEC_rhs_VelocityAdvection );

   // Solve advection problem
   SOLVER_A_VelocityUnsteady->solve( VEC_rhs_A_Velocity, VEC_U );   
   MAC_ASSERT( SOLVER_A_VelocityUnsteady->solution_is_achieved() ) ;

   return( true ) ;
   
}
