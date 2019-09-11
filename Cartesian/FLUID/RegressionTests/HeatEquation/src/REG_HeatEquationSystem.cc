#include <REG_HeatEquationSystem.hh>
#include <LA_Matrix.hh>
#include <LA_Vector.hh>
#include <LA_Scatter.hh>
#include <LA_SeqVector.hh>
#include <LA_SeqMatrix.hh>
#include <LA_Solver.hh>
#include <LA_MatrixIterator.hh>
#include <intVector.hh>
#include <MAC.hh>
#include <MAC_Error.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_Timer.hh>
#include <MAC_Vector.hh>
#include <MAC_Communicator.hh>
#include <MAC_Exec.hh>
#include <FV_DiscreteField.hh>
#include <FV_SystemNumbering.hh>
#include <iostream>
#include <math.h>


//----------------------------------------------------------------------
REG_HeatEquationSystem*
REG_HeatEquationSystem:: create( MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
	FV_DiscreteField* mac_tf )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_HeatEquationSystem:: create" ) ;
   MAC_CHECK_PRE( exp != 0 ) ;
   MAC_CHECK_PRE( mac_tf != 0 ) ;  

   REG_HeatEquationSystem* result = 
         new REG_HeatEquationSystem( a_owner, exp, mac_tf ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;	  
         
   return( result ) ;
   
}




//----------------------------------------------------------------------
REG_HeatEquationSystem:: REG_HeatEquationSystem(
	MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
	FV_DiscreteField* mac_tf )
//----------------------------------------------------------------------
   : MAC_Object( a_owner )
   , TF( mac_tf )
{
   MAC_LABEL( "REG_HeatEquationSystem:: REG_HeatEquationSystem" ) ; 

   // Build the matrices & vectors
   build_system(exp) ;
   re_initialize() ;
 
}




//----------------------------------------------------------------------
void
REG_HeatEquationSystem:: build_system( MAC_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_HeatEquationSystem:: build_system" ) ;

   string key = TF->discretization_type();
   
   // Temperature unsteady
   MAT_A_TemperatureUnsteady = LA_Matrix::make( this,
         exp->create_subexplorer( this, 
	 "MAT_A_TemperatureUnsteady_" + key ) ) ;
   VEC_rhs_A_TemperatureUnsteady = 
     	MAT_A_TemperatureUnsteady->create_vector(this) ;

   // Temperature Laplacian
   MAT_D_TemperatureUnsteadyPlusDiffusion = LA_Matrix::make( this,
         exp->create_subexplorer( this, 
	 "MAT_D_TemperatureDiffusion_" + key ) ) ;
   VEC_rhs_D_TemperatureDiffusionPlusBodyTerm = 
     	MAT_D_TemperatureUnsteadyPlusDiffusion->create_vector(this) ;

   // Unknowns vectors
   VEC_TF = MAT_D_TemperatureUnsteadyPlusDiffusion->create_vector( this ) ;
   VEC_TF_previoustime = 
     MAT_D_TemperatureUnsteadyPlusDiffusion->create_vector( this ) ;
   VEC_TF_timechange = 
     MAT_D_TemperatureUnsteadyPlusDiffusion->create_vector( this ) ;   

   // Local vector
   TF_LOC = LA_SeqVector::create( this, 0 ) ;
 
   // Solvers
   SOLVER_Temperature = LA_Solver::make( this,
               exp->create_subexplorer( this,
	       "SOLVER_Temperature_" + key ) ) ;
   SOLVER_Temperature->set_initial_guess_nonzero( true );	        
   
   // Temperature numbering
   TF_NUM = FV_SystemNumbering::create( this, TF ) ; 
     
}




//----------------------------------------------------------------------
void
REG_HeatEquationSystem:: re_initialize( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_HeatEquationSystem:: re_initialize" ) ;

   size_t tf_glob = TF->nb_global_unknowns() ;
   size_t tf_loc = TF->nb_local_unknowns() ;
   
   // Temperature unsteady
   MAT_A_TemperatureUnsteady->re_initialize( tf_glob, tf_glob ) ;
   VEC_rhs_A_TemperatureUnsteady->re_initialize( tf_glob ) ;      

   // Temperature Laplacian
   MAT_D_TemperatureUnsteadyPlusDiffusion->re_initialize( tf_glob, tf_glob ) ;
   VEC_rhs_D_TemperatureDiffusionPlusBodyTerm->re_initialize( tf_glob ) ;  

   // Unknowns vectors
   VEC_TF->re_initialize( tf_glob ) ;
   VEC_TF_previoustime->re_initialize( tf_glob ) ;
   VEC_TF_timechange->re_initialize( tf_glob ) ;      
   
   // Local vector
   TF_LOC->re_initialize( tf_loc ) ;
   
   // Temperature numbering
   TF_NUM->define_scatter( VEC_TF ) ;
      
}




//----------------------------------------------------------------------
REG_HeatEquationSystem:: ~REG_HeatEquationSystem( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_HeatEquationSystem:: ~REG_HeatEquationSystem" ) ;
}




//----------------------------------------------------------------------
void
REG_HeatEquationSystem::initialize_temperature( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_HeatEquationSystem:: initialize_temperature" ) ;

   TF->extract_unknown_DOFs_value( 0, TF_LOC ) ;
   TF_NUM->scatter()->set( TF_LOC, VEC_TF ) ;
      
}




//----------------------------------------------------------------------
LA_SeqVector const*
REG_HeatEquationSystem:: get_solution_temperature( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_HeatEquationSystem:: get_solution_temperature" ) ;

   TF_NUM->scatter()->get( VEC_TF, TF_LOC ) ;

   LA_SeqVector const* result = TF_LOC ;

   return( result ) ;
   
}




//----------------------------------------------------------------------
void
REG_HeatEquationSystem:: finalize_constant_matrices( 
	size_t DiffusionTimeAccuracy )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_HeatEquationSystem:: finalize_constant_matrices" ) ;

   bool same_pattern = false ;   
   MAC_Communicator const* macCOMM = MAC_Exec::communicator();
   
   if ( macCOMM->rank() == 0 ) 
     cout << "Unsteady + laplacian matrix" << endl;  

   // Use half diffusion coef 1/Pe such that the temperature operator is
   // MAT = 1 / dt - 0.5 * ( 1 / Pe ) * L(T)
   // !!! Note that the BC rhs vector VEC_rhs_D_TemperatureDiffusionPlusBodyTerm
   // contains the rhs of the total diffusion coef (not 0.5 * diffusion coef) 
   // operator, i.e., rhs = ( 1 / Pe ) * rhs_minus_L(T) !!!
   if ( DiffusionTimeAccuracy == 2 ) 
     MAT_D_TemperatureUnsteadyPlusDiffusion->scale( 0.5 ) ;   

   MAT_D_TemperatureUnsteadyPlusDiffusion->add_Mat( 
      	MAT_A_TemperatureUnsteady, 1., same_pattern ) ;    

   SOLVER_Temperature->set_matrix( MAT_D_TemperatureUnsteadyPlusDiffusion );
      
}




//----------------------------------------------------------------------
bool
REG_HeatEquationSystem:: HeatEquation_solver( size_t DiffusionTimeAccuracy,
	size_t iter_number )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_HeatEquationSystem:: HeatEquation_solver" ) ;

   // Synchronize matrices & vectors
   VEC_TF->synchronize() ;

   // First order scheme
   if( DiffusionTimeAccuracy == 1 )
   {
     // Compute temperature unsteady rhs = (1/Dt) * T
     MAT_A_TemperatureUnsteady->multiply_vec_then_add( VEC_TF,
   	VEC_rhs_A_TemperatureUnsteady ) ;

     // Add temperature BC rhs
     // rhs = (1/Dt) * u + ( 1 / Pe ) * rhs_minus_L(T) + body term   
     VEC_rhs_A_TemperatureUnsteady->sum( 
   	VEC_rhs_D_TemperatureDiffusionPlusBodyTerm, 1. ) ;
   }
   // Second order semi-implicit Crank-Nicholson scheme
   else    
   {
     // Compute temperature unsteady rhs
     // rhs = 2 * (1/Dt) * T
     MAT_A_TemperatureUnsteady->multiply_vec_then_add( VEC_TF,
   	VEC_rhs_A_TemperatureUnsteady, 2. ) ;
	
     // Compute velocity unsteady + viscous rhs
     // rhs = 2 * (1/Dt) * T + ( 1 / Pe ) * rhs_minus_L(T) + body term     
     VEC_rhs_A_TemperatureUnsteady->sum( 
   	VEC_rhs_D_TemperatureDiffusionPlusBodyTerm, 
	1. ); 
	
     // Add semi-explicit diffusion part
     // rhs = 2 * (1/Dt) * T + ( 1 / Pe ) * rhs_minus_L(T) + body term 
     // 		 - ( 1/Dt - ( ( 1 / Pe ) / 2 ) * L ) * T 
     //     = (1/Dt) * T + ( 1 / Pe ) * rhs_minus_L(T) + body term  
     // 		 + ( ( 1 / Pe ) / 2 ) * L * T     	
     MAT_D_TemperatureUnsteadyPlusDiffusion->multiply_vec_then_add( VEC_TF,
   	VEC_rhs_A_TemperatureUnsteady, -1., 1. ) ;        
   }

// 
//    // Compute temperature unsteady rhs   
//    MAT_A_TemperatureUnsteady->multiply_vec_then_add( VEC_TF,
//    	VEC_rhs_A_TemperatureUnsteady ) ;
// 
//    // Add temperature BC rhs
//    VEC_rhs_A_TemperatureUnsteady->sum( 
//    	VEC_rhs_D_TemperatureDiffusionPlusBodyTerm ) ;   

   // Solve unsteady laplacian problem
   SOLVER_Temperature->solve( VEC_rhs_A_TemperatureUnsteady, 
   	VEC_TF );   
   MAC_ASSERT( SOLVER_Temperature->solution_is_achieved() ) ;

   return( true ) ;
   
}




//----------------------------------------------------------------------
void
REG_HeatEquationSystem::at_each_time_step( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_HeatEquationSystem:: at_each_time_step" ) ;

   // Store temperature at previous time
   VEC_TF->synchronize() ;   
   VEC_TF_previoustime->set( VEC_TF ) ;
   
}




//----------------------------------------------------------------------
double 
REG_HeatEquationSystem:: compute_temperature_change( void )
//----------------------------------------------------------------------
{ 	
   MAC_LABEL( "REG_HeatEquationSystem:: compute_temperature_change" ) ;

   VEC_TF->synchronize() ;
   VEC_TF_timechange->set( VEC_TF ) ;
   VEC_TF_timechange->sum( VEC_TF_previoustime, -1.0 ) ;

   double norm_TF = VEC_TF->two_norm() ;
   double time_change = VEC_TF_timechange->two_norm() ;
   if ( norm_TF > 1e-4 ) time_change /= norm_TF;   

   return ( time_change ) ;

}



//----------------------------------------------------------------------
void
REG_HeatEquationSystem:: assemble_temperature_unsteady_matrix( 
	double const& coef )
//----------------------------------------------------------------------
{
   MAC_LABEL( 
   "REG_ProjectionNavierStokesSystem:: assemble_temperature_unsteady_matrix" ) ;

   TF->assemble_mass_matrix( coef, MAT_A_TemperatureUnsteady );
   
}




//----------------------------------------------------------------------
void
REG_HeatEquationSystem:: assemble_temperature_diffusion_matrix_rhs( 
	double const& coef_lap )
//----------------------------------------------------------------------
{
   MAC_LABEL( 
   "REG_HeatEquationSystem:: assemble_temperature_diffusion_matrix_rhs" ) ;

   TF->assemble_constantcoef_laplacian_matrix( coef_lap,
	MAT_D_TemperatureUnsteadyPlusDiffusion, 
	VEC_rhs_D_TemperatureDiffusionPlusBodyTerm );	
   
}




//----------------------------------------------------------------------
LA_Vector*
REG_HeatEquationSystem::get_diffrhs_plus_bodyterm_vector( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "REG_HeatEquationSystem:: get_diffrhs_plus_bodyterm_vector" ) ; 

   return ( VEC_rhs_D_TemperatureDiffusionPlusBodyTerm ) ;
	
}
