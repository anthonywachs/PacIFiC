#include <REG_ViscoplasticALG2System.hh>
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
REG_ViscoplasticALG2System*
REG_ViscoplasticALG2System:: create( MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
	FV_DiscreteField* mac_uu,
	FV_DiscreteField* mac_pp,
	FV_DiscreteField* mac_dd,
        const size_t &NS_Viscous_TimeAccuracy_,	
	const size_t &NS_Advection_TimeAccuracy_,
	const bool& b_PreconditionedWithLapP_, 
	const bool& b_pressure_rescaling_ )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ViscoplasticALG2System:: create" ) ;
   MAC_CHECK_PRE( exp != 0 ) ;
   MAC_CHECK_PRE( mac_uu != 0 ) ;
   MAC_CHECK_PRE( mac_pp != 0 ) ;     

   REG_ViscoplasticALG2System* result = 
         new REG_ViscoplasticALG2System( a_owner, exp, mac_uu, mac_pp,
	 	mac_dd, NS_Viscous_TimeAccuracy_, NS_Advection_TimeAccuracy_,
		b_PreconditionedWithLapP_, b_pressure_rescaling_ ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;	  
         
   return( result ) ;
   
}



//----------------------------------------------------------------------
REG_ViscoplasticALG2System:: REG_ViscoplasticALG2System(
	MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
	FV_DiscreteField* mac_uu,
	FV_DiscreteField* mac_pp,
	FV_DiscreteField* mac_dd,		
        const size_t &NS_Viscous_TimeAccuracy_,	
	const size_t &NS_Advection_TimeAccuracy_,
	const bool& b_PreconditionedWithLapP_, 
	const bool& b_pressure_rescaling_ )
//----------------------------------------------------------------------
   : MAC_Object( a_owner )
   , UU( mac_uu )
   , PP( mac_pp )
   , DD( mac_dd )
   , UU_NUM( 0 )
   , PP_NUM( 0 )
   , DD_NUM( 0 )    
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
   , VEC_q( 0 )
   , VEC_t( 0 )
   , VEC_r( 0 )
   , VEC_w( 0 )
   , VEC_x( 0 )
   , VEC_s( 0 )
   , VEC_z( 0 )
   , VEC_mean_pressure( 0 )
   , VEC_unit_pressure( 0 )
   , DD_LOC( 0 )
   , VEC_TENSOR( 0 )
   , VEC_U_previous_iteration( 0 )
   , VEC_rhs_FullViscoplastic( 0 )
   , MAT_T_LagrangeMultiplierDivergence( 0 )
   , VEC_LambdaMinusrd( 0 )
   , SOLVER_A_VelocityUnsteadyPlusViscous( 0 )
   , SOLVER_D_PressureLaplacian( 0 )
   , SOLVER_Uzawa( 0 )
   , VECS_Storage( 0 )
   , Uzawa_Stokes_precision( 1.e-8  )   
   , Uzawa_Stokes_maxiter( 100 )
   , NS_Viscous_TimeAccuracy( NS_Viscous_TimeAccuracy_ )
   , NS_Advection_TimeAccuracy( NS_Advection_TimeAccuracy_ )
   , b_PreconditionedWithLapP( b_PreconditionedWithLapP_ )
   , b_pressure_rescaling( b_pressure_rescaling_ )
{
   MAC_LABEL( "REG_ViscoplasticALG2System:: REG_ViscoplasticALG2System" ) ;

   // Read the Uzawa/Stokes precision
   if ( exp->has_entry( "Uzawa_Stokes_precision" ) )
     Uzawa_Stokes_precision=exp->double_data( "Uzawa_Stokes_precision" );


   // Read the maximum iterations number allowed for Uzawa/Stokes
   if ( exp->has_entry( "Uzawa_Stokes_maxiter" ) )
     Uzawa_Stokes_maxiter=exp->int_data( "Uzawa_Stokes_maxiter");


   // Build the matrices & vectors
   build_system( exp ) ;
   re_initialize() ;
   
}




//----------------------------------------------------------------------
void
REG_ViscoplasticALG2System:: build_system( MAC_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ViscoplasticALG2System:: build_system" ) ;

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
   
   // Work vectors for Uzawa
   // Unknown vectors
   VEC_t = MAT_A_VelocityUnsteady->create_vector( this ) ;
   VEC_s = MAT_D_PressureLaplacian->create_vector( this ) ;
   // Rhs vectors         
   VEC_q = MAT_A_VelocityUnsteady->create_vector( this ) ;
   VEC_r = MAT_D_PressureLaplacian->create_vector( this ) ;
   VEC_w = MAT_D_PressureLaplacian->create_vector( this ) ;
   VEC_x = MAT_D_PressureLaplacian->create_vector( this ) ;
   VEC_z = MAT_D_PressureLaplacian->create_vector( this ) ;
   VEC_mean_pressure = MAT_D_PressureLaplacian->create_vector( this ) ;
   VEC_unit_pressure = MAT_D_PressureLaplacian->create_vector( this ) ;

   // Specific to viscoplastic ALG2
   DD_LOC = LA_SeqVector::create( this, 0 ) ; 
   VEC_TENSOR = MAT_A_VelocityUnsteady->create_vector(this) ;
   VEC_rhs_FullViscoplastic = MAT_A_VelocityUnsteady
     	->create_vector(this) ;
   MAT_T_LagrangeMultiplierDivergence = LA_Matrix::make( this,
           exp->create_subexplorer( this, 
	   "MAT_T_LagrangeMultiplierDivergence" )  ) ;
   VEC_LambdaMinusrd = MAT_A_VelocityUnsteady->create_vector(this) ;
   VEC_U_previous_iteration = MAT_A_VelocityUnsteady->create_vector(this) ;      
   
   // Solver
   SOLVER_A_VelocityUnsteadyPlusViscous =
         LA_Solver::make( this,
               exp->create_subexplorer( this, 
	       "SOLVER_A_VelocityUnsteadyPlusViscous" ) ) ;
   SOLVER_D_PressureLaplacian =
         LA_Solver::make( this,
               exp->create_subexplorer( this,
	       "SOLVER_D_PressureLaplacian" ) ) ; 

   // Velocity, pressure and tensor numbering
   UU_NUM = FV_SystemNumbering::create( this, UU ) ;
   PP_NUM = FV_SystemNumbering::create( this, PP ) ;
   DD_NUM = FV_SystemNumbering::create( this, DD ) ;   

   // Vectors storage for reload
   VECS_Storage = LA_StorableVectors::create( this );   
}




//----------------------------------------------------------------------
void
REG_ViscoplasticALG2System:: re_initialize( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ViscoplasticALG2System:: re_initialize" ) ;

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

   // Work vectors for Uzawa
   // Unknown vectors
   VEC_t->re_initialize( u_glob ) ;
   VEC_s->re_initialize( p_glob ) ;
   // Rhs vectors
   VEC_q->re_initialize( u_glob ) ;
   VEC_r->re_initialize( p_glob ) ;
   VEC_w->re_initialize( p_glob ) ;
   VEC_x->re_initialize( p_glob ) ;
   VEC_z->re_initialize( p_glob ) ;
   VEC_mean_pressure->re_initialize( p_glob ) ;
   VEC_unit_pressure->re_initialize( p_glob ) ;

   // Specific to viscoplastic ALG2
   DD_LOC->re_initialize( DD->nb_local_unknowns() ) ;
   size_t dd_global = DD->nb_global_unknowns() ;
   VEC_TENSOR->re_initialize( dd_global ) ;
   VEC_rhs_FullViscoplastic->re_initialize( u_glob ) ;
   MAT_T_LagrangeMultiplierDivergence->re_initialize( u_glob, dd_global ) ;
   VEC_LambdaMinusrd->re_initialize( dd_global ) ;
   VEC_U_previous_iteration->re_initialize( u_glob ) ;
   
   // Velocity & pressure scatters
   UU_NUM->define_scatter( VEC_U ) ;
   PP_NUM->define_scatter( VEC_P ) ;   
   DD_NUM->define_scatter( VEC_LambdaMinusrd ) ;

}




//----------------------------------------------------------------------
REG_ViscoplasticALG2System:: ~REG_ViscoplasticALG2System( void )
//----------------------------------------------------------------------
{}




//----------------------------------------------------------------------
double 
REG_ViscoplasticALG2System:: get_Stokes_convergence_criterion( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( 
     "REG_ViscoplasticALG2System:: get_Stokes_convergence_criterion" ) ;
  
   return( Uzawa_Stokes_precision ) ;
   
}




//----------------------------------------------------------------------
void
REG_ViscoplasticALG2System::at_each_time_step( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ViscoplasticALG2System:: at_each_time_step" ) ;

   // Store velocity at previous time
   VEC_U->synchronize() ; 
   VEC_U_previoustime->set( VEC_U ) ;
   
}




//----------------------------------------------------------------------
double 
REG_ViscoplasticALG2System:: compute_velocity_change( void )
//----------------------------------------------------------------------
{ 	
   MAC_LABEL( "REG_ViscoplasticALG2System:: compute_velocity_change" ) ;

   VEC_U->synchronize() ;
   VEC_U_timechange->set( VEC_U ) ;
   VEC_U_timechange->sum( VEC_U_previoustime, -1.0 ) ;

   double norm_U = VEC_U->two_norm() ;
   double time_change = VEC_U_timechange->two_norm() ;
   if ( norm_U>1e-4 ) time_change /= norm_U;   

   return ( time_change );

}




//----------------------------------------------------------------------
double 
REG_ViscoplasticALG2System:: compute_velocity_change_iteration( void )
//----------------------------------------------------------------------
{ 	
   MAC_LABEL( "REG_ViscoplasticALG2System:: compute_velocity_change_iteration" ) ;

   VEC_U_timechange->set( VEC_U ) ;
   VEC_U_timechange->sum( VEC_U_previous_iteration, -1.0 ) ;

   double norm_U = VEC_U->two_norm() ;
   double iter_change = VEC_U_timechange->two_norm() ;
   if ( norm_U > 1e-4 ) iter_change /= norm_U;

   // Update previous iteration velocity vector
   VEC_U_previous_iteration->set( VEC_U ) ;
         	 
   return ( iter_change );

}




//----------------------------------------------------------------------
double
REG_ViscoplasticALG2System:: compute_velocity_divergence_norm( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( 
   	"REG_ViscoplasticALG2System:: compute_velocity_divergence_norm" ) ;

   MAT_B_VelocityDivergence->multiply_vec_then_add( VEC_U, VEC_r, -1.0 ) ;
   VEC_r->sum( VEC_rhs_B_VelocityDivergence ) ;

   double normdivu = VEC_r->two_norm() ;

   return normdivu;
	
}




//----------------------------------------------------------------------
void
REG_ViscoplasticALG2System::nullify_velocity_advection_rhs( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ViscoplasticALG2System:: nullify_velocity_advection_rhs" ) ;

   // Nullify inertia rhs
   VEC_rhs_VelocityAdvection->nullify() ;
   
}




//----------------------------------------------------------------------
bool 
REG_ViscoplasticALG2System:: Uzawa_solver( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ViscoplasticALG2System:: Uzawa_solver" ) ;   

   double nr,nrkm1,nwx,alpha,beta;
   size_t iter=0;
   MAC_Communicator const* macCOMM = MAC_Exec::communicator();

   SOLVER_A_VelocityUnsteadyPlusViscous->set_initial_guess_nonzero( false );

   // Initialisation
   // --------------
   // Compute q=f-Bt.p
   VEC_q->set( VEC_rhs_A_Velocity ) ;
   VEC_P->synchronize() ;
   MAT_B_VelocityDivergence->tr_multiply_vec_then_add( VEC_P, VEC_q, -1.0, 
   	1.0) ;

   // Solution of A.u=q
   SOLVER_Uzawa->solve( VEC_q, VEC_U ) ;
   MAC_ASSERT( SOLVER_Uzawa->solution_is_achieved() ) ;
 
   // Compute r=-B.u+g
   VEC_rhs_B_VelocityDivergence->synchronize() ;
   VEC_r->set( VEC_rhs_B_VelocityDivergence ) ;
   MAT_B_VelocityDivergence->multiply_vec_then_add( VEC_U, VEC_r, -1.0, 1.0 ) ;

   // Compute nr=||r||^2
   nr = MAC::sqr( VEC_r->two_norm() ) ;

   // Compute w=r   
   VEC_w->set( VEC_r ) ;

   // Iterative process
   // -----------------
   while( sqrt(nr) > Uzawa_Stokes_precision && iter < Uzawa_Stokes_maxiter )
   {
      iter++;
      // Compute q=Bt.w
      MAT_B_VelocityDivergence->tr_multiply_vec_then_add( VEC_w, VEC_q ) ;

      // Solution of A.t=q
      SOLVER_Uzawa->solve( VEC_q, VEC_t ) ;
      MAC_ASSERT( SOLVER_Uzawa->solution_is_achieved() ) ;

      // Compute x=B.t
      MAT_B_VelocityDivergence->multiply_vec_then_add( VEC_t, VEC_x, 1.0, 
      	0.0 ) ;

      // Compute nrkm1=r.r
      nrkm1 = MAC::sqr( VEC_r->two_norm() ) ;

      // Compute nws=w.x
      nwx = VEC_w->dot( VEC_x );
      
      // Compute alpha=nrkm1/nwx
      alpha = nrkm1 / nwx ;

      // Update p-=alpha.w, r-=alpha.x, u+=alpha.t and compute nr
      VEC_P->sum( VEC_w , -alpha ) ;
      VEC_r->sum( VEC_x , -alpha ) ;
      VEC_U->sum( VEC_t , alpha ) ;

      nr = MAC::sqr( VEC_r->two_norm() );

      // Compute beta=nr/nrkm1
      beta = nr / nrkm1 ;

      // Update w=r+beta.w
      VEC_w->scale( beta ) ;
      VEC_w->sum( VEC_r ) ;

//       if ( macCOMM->rank() == 0 )
//          MAC::out() << "Residuals = " << 
// 	 	MAC::doubleToString( ios::scientific, 8, sqrt(nr) ) 
// 		<< "  Iteration = " << iter << endl;
   }

   if ( macCOMM->rank() == 0 )
      MAC::out() << "Residuals = " << 
	MAC::doubleToString( ios::scientific, 8, sqrt(nr) ) 
	<< "  Iteration = " << iter << endl;
   
   return( sqrt(nr) < Uzawa_Stokes_precision ) ;
   
}




//----------------------------------------------------------------------
bool 
REG_ViscoplasticALG2System:: PreconditionedStokes_Uzawa_solver( 
	const double &mupr, const double &gamma )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ViscoplasticALG2System:: Preconditioned_Uzawa_solver" ) ;
   
   double nr, nrzkm1, nwx, nrz, alpha, beta, smean, tmp = 0.;
   size_t iter = 0, niter_sublin = 0;
   MAC_Communicator const* macCOMM = MAC_Exec::communicator();
   
   // Initialisation
   // --------------
   // Compute q=f-Bt.p
   VEC_q->set( VEC_rhs_A_Velocity ) ;   
   MAT_B_VelocityDivergence->tr_multiply_vec_then_add( VEC_P, VEC_q, -1.0, 
   	1.0 );

   // Solution of A.u=q
   SOLVER_Uzawa->solve( VEC_q, VEC_U );
   MAC_ASSERT( SOLVER_Uzawa->solution_is_achieved() ) ;

   niter_sublin += SOLVER_Uzawa->nb_iterations_achieved();

   // Compute r=-B.u+g
   VEC_rhs_B_VelocityDivergence->synchronize();
   VEC_r->set( VEC_rhs_B_VelocityDivergence );
   MAT_B_VelocityDivergence->multiply_vec_then_add( VEC_U, VEC_r, -1.0, 1.0 );

   // Compute nr=||r||^2
   nr = MAC::sqr(VEC_r->two_norm());

   // Compute w=r
   VEC_w->set(VEC_r);

   // Preconditioner
   if ( b_pressure_rescaling ) 
   {
     if ( macCOMM->rank() == 0 ) 
     {
       tmp = VEC_r->item( 0 ) ;
       VEC_r->set_item( 0, 0. ) ;
     }
     VEC_r->synchronize() ; 
   } 
   SOLVER_D_PressureLaplacian->solve( VEC_r, VEC_s );
   MAC_ASSERT( SOLVER_D_PressureLaplacian->solution_is_achieved() ) ;

   niter_sublin += SOLVER_D_PressureLaplacian->nb_iterations_achieved();
   if( b_pressure_rescaling )
   {
      if ( macCOMM->rank() == 0 ) VEC_r->set_item( 0, tmp ) ; 
      VEC_r->synchronize() ;
      smean = VEC_s->dot( VEC_mean_pressure );
      VEC_s->sum( VEC_unit_pressure, -smean );
   }

   // Compute z and w
   VEC_z->set( VEC_r );
   VEC_z->scale( mupr / ( mupr + gamma ) );
   VEC_z->sum( VEC_s, gamma / ( mupr + gamma ) );
   VEC_w->set( VEC_z );

   // Iterative process
   // -----------------
   while ((sqrt(nr) > Uzawa_Stokes_precision) && (iter < Uzawa_Stokes_maxiter))
   {
      iter++;
      // Compute q=Bt.w
      MAT_B_VelocityDivergence->tr_multiply_vec_then_add( VEC_w, VEC_q, 1.0, 
      	0.0 );

      // Solution of A.t=q
      SOLVER_Uzawa->solve( VEC_q, VEC_t );
      niter_sublin += SOLVER_Uzawa->nb_iterations_achieved();

      // Compute x=B.t
      MAT_B_VelocityDivergence->multiply_vec_then_add( VEC_t, VEC_x, 1.0, 0.0) ;

      // Compute nrkm1=r.z
      nrzkm1 = VEC_r->dot( VEC_z );

      // Compute nws=w.x
      nwx = VEC_w->dot( VEC_x );

      // Compute alpha=nrkm1/nwx
      alpha = nrzkm1 / nwx ;

      // Update p-=alpha.w, r-=alpha.x and compute nr
      VEC_P->sum( VEC_w , -alpha ) ;
      VEC_r->sum( VEC_x , -alpha ) ;
      VEC_U->sum( VEC_t , alpha ) ;

      nr=MAC::sqr( VEC_r->two_norm() );

      // Preconditioner
      if ( b_pressure_rescaling ) 
      {
        if ( macCOMM->rank() == 0 ) 
	{
	  tmp = VEC_x->item( 0 ) ;
          VEC_x->set_item( 0, 0. ) ; 
	}
	VEC_x->synchronize() ; 
      }  
      SOLVER_D_PressureLaplacian->solve( VEC_x, VEC_s ) ;
      MAC_ASSERT( SOLVER_D_PressureLaplacian->solution_is_achieved() ) ;
      niter_sublin += SOLVER_D_PressureLaplacian->nb_iterations_achieved();
      if ( b_pressure_rescaling )
      {
         if ( macCOMM->rank() == 0 ) VEC_x->set_item( 0, tmp ) ;
	 VEC_x->synchronize() ;  
	 smean=VEC_s->dot( VEC_mean_pressure );
         VEC_s->sum( VEC_unit_pressure, -smean ) ;
      }

      // Compute z
      VEC_z->sum( VEC_x, -alpha * mupr / ( mupr + gamma ) );
      VEC_z->sum( VEC_s, -alpha * gamma / ( mupr + gamma ) );
      nrz = VEC_r->dot( VEC_z );
      
      // Compute beta=nr/nrkm1
      beta = nrz / nrzkm1 ;

      // Update w=r+beta.w
      VEC_w->scale( beta ) ;
      VEC_w->sum( VEC_z );

//       if ( macCOMM->rank() == 0 )
//          MAC::out() << "Residuals = " << 
// 	 	MAC::doubleToString( ios::scientific, 8, sqrt(nr) ) 
// 		<< "  Iteration = " << iter << endl;
   }
   
   if ( b_pressure_rescaling )
   {
      smean=VEC_P->dot( VEC_mean_pressure );
      VEC_P->sum( VEC_unit_pressure, -smean ) ;
   }

   if ( macCOMM->rank() == 0 )
     MAC::out() << "Residuals = " << 
     	MAC::doubleToString( ios::scientific, 8, sqrt(nr) )  << 
	"  Nb iterations = " << iter << endl
   	<< "Nb iterations sub-linear systems = " << niter_sublin << endl;

   return( sqrt(nr) < Uzawa_Stokes_precision );
   
}




//----------------------------------------------------------------------
LA_SeqVector const*
REG_ViscoplasticALG2System:: get_solution_velocity( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ViscoplasticALG2System:: get_solution_U" ) ;

   UU_NUM->scatter()->get( VEC_U, U_LOC ) ;

   LA_SeqVector const* result = U_LOC ;

   return( result ) ;
   
}




//----------------------------------------------------------------------
LA_SeqVector const*
REG_ViscoplasticALG2System:: get_solution_pressure( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ViscoplasticALG2System:: get_solution_P" ) ;

   PP_NUM->scatter()->get( VEC_P, P_LOC ) ;

   LA_SeqVector const* result = P_LOC ;

   return( result ) ;
   
}




//----------------------------------------------------------------------
void
REG_ViscoplasticALG2System::initialize_velocity( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ViscoplasticALG2System:: initialize_velocity" ) ;

   UU->extract_unknown_DOFs_value( 0, U_LOC ) ;
   UU_NUM->scatter()->set( U_LOC, VEC_U ) ;
         
}




//----------------------------------------------------------------------
void
REG_ViscoplasticALG2System::initialize_pressure( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ViscoplasticALG2System:: initialize_pressure" ) ;
   
   PP->extract_unknown_DOFs_value( 0, P_LOC ) ;
   PP_NUM->scatter()->set( P_LOC, VEC_P ) ;
         
}




//----------------------------------------------------------------------
void
REG_ViscoplasticALG2System::initialize_lambda( size_t lambda_level )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ViscoplasticALG2System:: initialize_lambda" ) ;
   
   DD->extract_unknown_DOFs_value( lambda_level, DD_LOC ) ;
   DD_NUM->scatter()->set( DD_LOC, VEC_LambdaMinusrd ) ;
         
}




//----------------------------------------------------------------------
void
REG_ViscoplasticALG2System::initialize_d( size_t d_level )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ViscoplasticALG2System:: initialize_d" ) ;
   
   DD->extract_unknown_DOFs_value( d_level, DD_LOC ) ;
   DD_NUM->scatter()->set( DD_LOC, VEC_TENSOR ) ;
         
}




//----------------------------------------------------------------------
void
REG_ViscoplasticALG2System:: finalize_constant_matrices( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ViscoplasticALG2System:: finalize_constant_matrices" ) ;

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
	
   if ( b_PreconditionedWithLapP ) 
     SOLVER_D_PressureLaplacian->set_matrix( MAT_D_PressureLaplacian ) ;
      
}




//----------------------------------------------------------------------
void
REG_ViscoplasticALG2System:: compute_velocityAdvectionDiffusion_rhs( 
	bool const& b_restart,
	size_t const& iteration_number, 
	bool const& b_with_advection )
//----------------------------------------------------------------------
{
   MAC_LABEL(
     "REG_ViscoplasticALG2System:: compute_velocityAdvectionDiffusion_rhs" ) ;
     
   VEC_U->synchronize() ;
   
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

}




//----------------------------------------------------------------------
void
REG_ViscoplasticALG2System:: Uzawa_NavierStokes_solver( 
	const double &viscosity, const double &density, 
	const double &timestep )
//----------------------------------------------------------------------
{
   MAC_LABEL( 
     "REG_ViscoplasticALG2System:: Uzawa_NavierStokes_solver" ) ;

   // Set velocity rhs
   VEC_rhs_A_Velocity->set( VEC_rhs_VelocityAdvectionDiffusion ) ;

   // Add the div(lambda-r.d) term to velocity rhs
   MAT_T_LagrangeMultiplierDivergence->multiply_vec_then_add(
   	VEC_LambdaMinusrd, VEC_rhs_A_Velocity, 1.0, 1.0 ) ;

   // Set the velocity solver
   SOLVER_A_VelocityUnsteadyPlusViscous
   	->set_initial_guess_nonzero( true );
   SOLVER_Uzawa = SOLVER_A_VelocityUnsteadyPlusViscous ;
   
   // Solve the saddle-point Stokes problem
   if ( b_PreconditionedWithLapP )
     PreconditionedStokes_Uzawa_solver( viscosity, density / timestep );
   else Uzawa_solver();	
   
}




//----------------------------------------------------------------------
void
REG_ViscoplasticALG2System:: assemble_velocity_viscous_matrix_rhs( 
	double const& coef_lap )
//----------------------------------------------------------------------
{
   MAC_LABEL( 
   	"REG_ViscoplasticALG2System:: assemble_velocity_viscous_matrix_rhs" ) ;

   UU->assemble_constantcoef_laplacian_matrix( coef_lap,
	MAT_A_VelocityUnsteadyPlusViscous, VEC_rhs_A_VelocityViscous );	
   
}




//----------------------------------------------------------------------
void
REG_ViscoplasticALG2System:: assemble_velocity_unsteady_matrix( 
	double const& coef )
//----------------------------------------------------------------------
{
   MAC_LABEL( 
   	"REG_ViscoplasticALG2System:: assemble_velocity_unsteady_matrix" ) ;

   UU->assemble_mass_matrix( coef, MAT_A_VelocityUnsteady );
   
}




//----------------------------------------------------------------------
void
REG_ViscoplasticALG2System:: assemble_pdivv_matrix_rhs( double const& coef )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ViscoplasticALG2System:: assemble_pdivv_matrix_rhs" ) ;

   PP->assemble_pDivv_matrix( UU, coef, MAT_B_VelocityDivergence,
   	VEC_rhs_B_VelocityDivergence );		
   
}




//----------------------------------------------------------------------
void
REG_ViscoplasticALG2System:: assemble_velocity_advection( 
	string const& AdvectionScheme,
      	size_t advecting_level, double const& coef, 
	size_t advected_level  )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ViscoplasticALG2System:: assemble_velocity_advection" ) ;

   if ( AdvectionScheme == "Upwind" )
     UU->assemble_advection_Upwind( UU, advecting_level, coef,
     	advected_level, VEC_rhs_VelocityAdvection ); 
   else
     UU->assemble_advection_TVD( UU, advecting_level, coef,
     	advected_level, VEC_rhs_VelocityAdvection );

} 




//----------------------------------------------------------------------
void
REG_ViscoplasticALG2System:: assemble_pressure_laplacian_matrix_rhs( 
	double const& coef_lap )
//----------------------------------------------------------------------
{
   MAC_LABEL( 
     "REG_ViscoplasticALG2System:: assemble_pressure_laplacian_matrix_rhs" ) ;

   PP->assemble_constantcoef_laplacian_matrix( coef_lap, 
   	MAT_D_PressureLaplacian, VEC_rhs_D_PressureLaplacian, 
	b_pressure_rescaling );
   
}




//----------------------------------------------------------------------
void 
REG_ViscoplasticALG2System:: pressure_laplacian_correction( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ViscoplasticALG2System:: pressure_laplacian_correction" ) ;
      
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

}




//----------------------------------------------------------------------
void
REG_ViscoplasticALG2System:: assemble_tauGradv_tensor_divergence_matrix( 
	double const& coef )
//----------------------------------------------------------------------
{
   MAC_LABEL( 
   "REG_ViscoplasticALG2System:: assemble_tauGradv_tensor_divergence_matrix" ) ;

   UU->assemble_tauGradv_tensor_divergence_matrix( DD, coef, 
   	MAT_T_LagrangeMultiplierDivergence );
   
}




//----------------------------------------------------------------------
void
REG_ViscoplasticALG2System::initialize_previous_iteration_velocity( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( 
     "REG_ViscoplasticALG2System:: initialize_previous_iteration_velocity" ) ;
   
   VEC_U_previous_iteration->set( VEC_U ) ;
   
}




//----------------------------------------------------------------------
void
REG_ViscoplasticALG2System::compute_viscoplastic_LambdaMinusrd_vector( 
	double const& ALG2_r )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ViscoplasticALG2System:: "
   	"compute_viscoplastic_LambdaMinusrd_vector" ) ;
   
   VEC_LambdaMinusrd->sum( VEC_TENSOR, - ALG2_r ) ;
   
}




//----------------------------------------------------------------------
LA_Vector*
REG_ViscoplasticALG2System::get_tensor_unknown_vector( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ViscoplasticALG2System:: get_tensor_unknown_vector" ) ;
   
   return ( VEC_TENSOR ) ;
   
}




//----------------------------------------------------------------------
LA_SeqVector const*
REG_ViscoplasticALG2System:: get_solution_strainrate_tensor( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_NavierStokesSystem:: get_solution_strainrate_tensor" ) ;

   DD_NUM->scatter()->get( VEC_TENSOR, DD_LOC ) ;

   LA_SeqVector const* result = DD_LOC ;

   return( result ) ;
   
}




//----------------------------------------------------------------------
LA_Vector*
REG_ViscoplasticALG2System::get_LambdaMinusrd_vector( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ViscoplasticALG2System:: get_LambdaMinusrd_vector" ) ;
   
   return ( VEC_LambdaMinusrd ) ;
   
}




//----------------------------------------------------------------------
LA_SeqVector const*
REG_ViscoplasticALG2System:: get_solution_lambda_tensor( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ViscoplasticALG2System:: get_solution_lambda_tensor" ) ;

   DD_NUM->scatter()->get( VEC_LambdaMinusrd, DD_LOC ) ;

   LA_SeqVector const* result = DD_LOC ;

   return( result ) ;
   
}




//----------------------------------------------------------------------
void
REG_ViscoplasticALG2System::add_storable_objects(
	MAC_ListIdentity* list ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_ViscoplasticALG2System:: add_storable_objects" ) ; 

   list->extend( VECS_Storage );
   if ( NS_Advection_TimeAccuracy == 2 )
   {
     VEC_rhs_VelocityAdvection_Nm2->set_name( "UgradU_nm2" );
     VECS_Storage->add_vector_to_store( VEC_rhs_VelocityAdvection_Nm2 ) ;     
   }
	
}
