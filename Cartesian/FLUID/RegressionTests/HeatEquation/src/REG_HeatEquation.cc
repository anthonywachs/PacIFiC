#include <REG_HeatEquation.hh>
#include <FV_DomainAndFields.hh>
#include <FV_DiscreteField.hh>
#include <REG_HeatEquationSystem.hh>
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
#include <PAC_Misc.hh>
#include <time.h>
#include <sys/time.h>
#include <math.h>


REG_HeatEquation const* REG_HeatEquation::PROTOTYPE 
                                                 = new REG_HeatEquation() ;


//---------------------------------------------------------------------------
REG_HeatEquation:: REG_HeatEquation( void )
//--------------------------------------------------------------------------
   : FV_OneStepIteration( "REG_HeatEquation" )
   , PAC_ComputingTime("Solver")      
{
   MAC_LABEL( "REG_HeatEquation:: REG_HeatEquation" ) ;

}




//---------------------------------------------------------------------------
REG_HeatEquation*
REG_HeatEquation:: create_replica( MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "REG_HeatEquation:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, dom, exp ) ) ;

   REG_HeatEquation* result = 
                        new REG_HeatEquation( a_owner, dom, exp ) ;

   MAC_CHECK( create_replica_POST( result, a_owner, dom, exp ) ) ;
   return( result ) ;
   
}




//---------------------------------------------------------------------------
REG_HeatEquation:: REG_HeatEquation( MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : FV_OneStepIteration( a_owner, dom, exp )
   , PAC_ComputingTime("Solver") 
   , TF ( dom->discrete_field( "temperature" ) ) 
   , TF_STAG ( dom->discrete_field( "temperature_staggered" ) ) 
   , TF_VERTEX ( dom->discrete_field( "temperature_vertex" ) )    
   , TF_ERROR( 0 )
   , TF_STAG_ERROR( 0 )
   , TF_VERTEX_ERROR( 0 )         
   , GLOBAL_EQ( 0 ) 
   , GLOBAL_EQ_STAG( 0 ) 
   , GLOBAL_EQ_VERTEX( 0 )      
   , peclet( 1. ) 
   , b_bodyterm( false )
   , DiffusionTimeAccuracy( 1 ) 
{
   MAC_LABEL( "REG_HeatEquation:: REG_HeatEquation" ) ;
   MAC_ASSERT( TF->discretization_type() == "centered" ) ;
   MAC_ASSERT( TF_STAG->discretization_type() == "staggered" ) ;
   MAC_ASSERT( TF_VERTEX->discretization_type() == "vertex" ) ;    
   MAC_ASSERT( TF->storage_depth() == 1 ) ;
   MAC_ASSERT( TF_STAG->storage_depth() == 1 ) ; 
   MAC_ASSERT( TF_VERTEX->storage_depth() == 1 ) ;         

   // Call of MAC_Communicator routine to set the rank of each proces and  
   // the number of processes during execution 	   
   pelCOMM = MAC_Exec::communicator();	   		
   my_rank = pelCOMM->rank();
   nb_procs = pelCOMM->nb_ranks();
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


   // Clear results directory in case of a new run
   if( !b_restart ) PAC_Misc::clearAllFiles( "Res", "Savings", my_rank ) ;
   
      
   // Get space dimension
   dim = TF->primary_grid()->nb_space_dimensions() ;   
   if ( dim == 1 )
   {
     string error_message="Space dimension should either 2 or 3";
     MAC_Error::object()->raise_bad_data_value(exp, 
	"nb_space_dimensions",
	error_message );     
   } 


   // Read Peclet number
   if ( exp->has_entry( "Peclet" ) )
   {
     peclet = exp->double_data( "Peclet" ) ; 
     exp->test_data( "Peclet", "Peclet>0." ) ; 
   }    


   // Read with or without body term
   if ( exp->has_entry( "BodyTerm" ) )
     b_bodyterm = exp->bool_data( "BodyTerm" ) ; 


   // Diffusion term time accuracy
   if ( exp->has_entry( "DiffusionTimeAccuracy" ) )
     DiffusionTimeAccuracy = exp->int_data( "DiffusionTimeAccuracy" );
   if ( DiffusionTimeAccuracy != 1 && DiffusionTimeAccuracy != 2 )
   {
     string error_message ="   - 1\n   - 2\n   ";
     MAC_Error::object()->raise_bad_data_value( exp,
         "DiffusionTimeAccuracy", error_message );
   }


   // Build the matrix system
   MAC_ModuleExplorer* se = 
	exp->create_subexplorer( 0, 
	"REG_HeatEquationSystem_"+TF->discretization_type() ) ; 
   GLOBAL_EQ = REG_HeatEquationSystem::create( this, se, TF ) ; 
   se->destroy() ;
   
   MAC_ModuleExplorer* se_staggered = 
	exp->create_subexplorer( 0, 
	"REG_HeatEquationSystem_"+TF_STAG->discretization_type() ) ; 
   GLOBAL_EQ_STAG = REG_HeatEquationSystem::create( this, se_staggered, 
   	TF_STAG ) ; 
   se_staggered->destroy() ; 
   
   MAC_ModuleExplorer* se_vertex = 
	exp->create_subexplorer( 0, 
	"REG_HeatEquationSystem_"+TF_VERTEX->discretization_type() ) ; 
   GLOBAL_EQ_VERTEX = REG_HeatEquationSystem::create( this, se_staggered, 
   	TF_VERTEX ) ; 
   se_vertex->destroy() ;
         

   // Duplicates field for error calculation in case of sinusoidal
   // solution with body term 3.pi^2.sin(pi.x).sin(pi.y).sin(pi.z)
   // on a [0:1]x[0:1]x[0:1] domain
   if ( b_bodyterm )
   {
     const_cast<FV_DomainAndFields*>(dom)->duplicate_field( 
   	"temperature", "tf_error" ) ;
     const_cast<FV_DomainAndFields*>(dom)->duplicate_field( 
   	"temperature_staggered", "tfstag_error" ) ;
     const_cast<FV_DomainAndFields*>(dom)->duplicate_field( 
   	"temperature_vertex", "tfvertex_error" ) ;

     TF_ERROR = dom->discrete_field( "tf_error" ) ;
     TF_ERROR->set_BC_values_modif_status( true ) ;
     TF_STAG_ERROR = dom->discrete_field( "tfstag_error" ) ;
     TF_STAG_ERROR->set_BC_values_modif_status( true ) ;     
     TF_VERTEX_ERROR = dom->discrete_field( "tfvertex_error" ) ;
     TF_VERTEX_ERROR->set_BC_values_modif_status( true ) ;                    
   }	
   

   // Timing routines
   if ( my_rank == is_master )
   {
     SCT_insert_app("Matrix_Assembly&Initialization"); 
     SCT_insert_app("Matrix_Solution");           
     SCT_get_elapsed_time("Objects_Creation");
   }
}




//---------------------------------------------------------------------------
REG_HeatEquation:: ~REG_HeatEquation( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "REG_HeatEquation:: ~REG_HeatEquation" ) ;

}




//---------------------------------------------------------------------------
void
REG_HeatEquation:: do_one_inner_iteration( FV_TimeIterator const* t_it ) 
//---------------------------------------------------------------------------
{
   MAC_LABEL( "REG_HeatEquation:: do_one_inner_iteration" ) ;
   MAC_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;

   start_total_timer( "REG_HeatEquation:: do_one_inner_iteration" ) ;
   start_solving_timer() ;

   if ( my_rank == is_master ) SCT_set_start("Matrix_Solution");

   // Solve heat equation
   GLOBAL_EQ->HeatEquation_solver( DiffusionTimeAccuracy,
   	t_it->iteration_number() );
   GLOBAL_EQ_STAG->HeatEquation_solver( DiffusionTimeAccuracy,
   	t_it->iteration_number() ); 
   GLOBAL_EQ_VERTEX->HeatEquation_solver( DiffusionTimeAccuracy,
   	t_it->iteration_number() );      

   if ( my_rank == is_master ) SCT_get_elapsed_time( "Matrix_Solution" );
   TF->update_free_DOFs_value( 0, GLOBAL_EQ->get_solution_temperature() ) ;
   TF_STAG->update_free_DOFs_value( 0, 
   	GLOBAL_EQ_STAG->get_solution_temperature() ) ;
   TF_VERTEX->update_free_DOFs_value( 0, 
   	GLOBAL_EQ_VERTEX->get_solution_temperature() ) ;
   
   stop_solving_timer() ;
   stop_total_timer() ;  
   
}




//---------------------------------------------------------------------------
void
REG_HeatEquation:: do_before_time_stepping( FV_TimeIterator const* t_it, 
      	std::string const& basename )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "REG_HeatEquation:: do_before_time_stepping" ) ;
   
   start_total_timer( "REG_HeatEquation:: do_before_time_stepping" ) ;

   if ( my_rank == is_master ) SCT_set_start("Matrix_Assembly&Initialization");
    
   FV_OneStepIteration::do_before_time_stepping( t_it, basename ) ;

   // Diffusion matrix and rhs
   // Note: we assemble here the total diffusion matrix, in case of 2nd order 
   // Crank-Nicholson scheme for the diffusion term, half diffusion coef for the
   // velocity operator is taken care of at the matrix level in
   // REG_HeatEquationSystem:: finalize_constant_matrices
   
   // Assemble constant matrices and operators
   // Centered
   GLOBAL_EQ->assemble_temperature_unsteady_matrix( 1. / t_it->time_step() ) ;
   GLOBAL_EQ->assemble_temperature_diffusion_matrix_rhs( - 1. / peclet );
   if ( b_bodyterm ) assemble_temperature_bodyterm_rhs( TF, 
   	GLOBAL_EQ->get_diffrhs_plus_bodyterm_vector() ) ;      
   GLOBAL_EQ->finalize_constant_matrices( DiffusionTimeAccuracy ) ;

   // Staggered
   GLOBAL_EQ_STAG->assemble_temperature_unsteady_matrix( 
   	1. / t_it->time_step() ) ;
   GLOBAL_EQ_STAG->assemble_temperature_diffusion_matrix_rhs( - 1. / peclet );
   if ( b_bodyterm ) assemble_temperature_bodyterm_rhs( TF_STAG, 
   	GLOBAL_EQ_STAG->get_diffrhs_plus_bodyterm_vector() ) ;        
   GLOBAL_EQ_STAG->finalize_constant_matrices( DiffusionTimeAccuracy ) ;
   
   // Vertex
   GLOBAL_EQ_VERTEX->assemble_temperature_unsteady_matrix( 
   	1. / t_it->time_step() ) ;
   GLOBAL_EQ_VERTEX->assemble_temperature_diffusion_matrix_rhs( - 1. / peclet );
   if ( b_bodyterm ) assemble_temperature_bodyterm_rhs( TF_VERTEX, 
   	GLOBAL_EQ_VERTEX->get_diffrhs_plus_bodyterm_vector() ) ;        
   GLOBAL_EQ_VERTEX->finalize_constant_matrices( DiffusionTimeAccuracy ) ;

   // Initialize temperature vector at the matrix level
   GLOBAL_EQ->initialize_temperature();
   GLOBAL_EQ_STAG->initialize_temperature();
   GLOBAL_EQ_VERTEX->initialize_temperature();   

   if ( my_rank == is_master )
     SCT_get_elapsed_time( "Matrix_Assembly&Initialization" );

   stop_total_timer() ;  

}




//---------------------------------------------------------------------------
void
REG_HeatEquation:: do_after_time_stepping( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "REG_HeatEquation:: do_after_time_stepping" ) ;  

   // Elapsed time by sub-problems
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
REG_HeatEquation:: do_before_inner_iterations_stage( 
	FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "REG_HeatEquation:: do_before_inner_iterations_stage" ) ;
   
   start_total_timer( "REG_HeatEquation:: do_before_inner_iterations_stage" ) ;

   FV_OneStepIteration::do_before_inner_iterations_stage( t_it ) ;

   // Perform matrix level operations before each time step 
   GLOBAL_EQ->at_each_time_step( );
   GLOBAL_EQ_STAG->at_each_time_step( );
   GLOBAL_EQ_VERTEX->at_each_time_step( );     

   stop_total_timer() ;
   
}




//---------------------------------------------------------------------------
void
REG_HeatEquation:: do_after_inner_iterations_stage( 
	FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "REG_HeatEquation:: do_after_inner_iterations_stage" ) ;
   
   start_total_timer( "REG_HeatEquation:: do_after_inner_iterations_stage" ) ;
   
   FV_OneStepIteration::do_after_inner_iterations_stage( t_it ) ;

   // Compute temperature change over the time step  
   double temperature_time_change = GLOBAL_EQ->compute_temperature_change()
   	/ t_it->time_step() ;
   if ( my_rank == is_master )
     cout << "Centered Temperature change = " << 
     	MAC::doubleToString( ios::scientific, 5, temperature_time_change ) 
	<< endl;     

   temperature_time_change = GLOBAL_EQ_STAG->compute_temperature_change()
   	/ t_it->time_step() ;
   if ( my_rank == is_master )
     cout << "Staggered Temperature change = " << 
     	MAC::doubleToString( ios::scientific, 5, temperature_time_change ) 
	<< endl;

   temperature_time_change = GLOBAL_EQ_VERTEX->compute_temperature_change()
   	/ t_it->time_step() ;
   if ( my_rank == is_master )
     cout << "Vertex Temperature change = " << 
     	MAC::doubleToString( ios::scientific, 5, temperature_time_change ) 
	<< endl;

   stop_total_timer() ;
   
}
  



//---------------------------------------------------------------------------
void
REG_HeatEquation:: do_additional_savings( FV_TimeIterator const* t_it,
      	int const& cycleNumber )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "REG_HeatEquation:: do_additional_savings" ) ;

   start_total_timer( "REG_HeatEquation:: do_additional_savings" ) ;

   if ( b_bodyterm )
   {
     error_with_analytical_solution( TF, TF_ERROR );
     error_with_analytical_solution( TF_STAG, TF_STAG_ERROR ); 
     error_with_analytical_solution( TF_VERTEX, TF_VERTEX_ERROR );      
   }
   
   compute_L2Norm_solution ( TF, 0 );   
   for (size_t comp=0;comp<TF_STAG->nb_components();++comp)
     compute_L2Norm_solution ( TF_STAG, comp );
   compute_L2Norm_solution ( TF_VERTEX, 0 );      

   stop_total_timer() ;
    
}




//---------------------------------------------------------------------------
void
REG_HeatEquation:: assemble_temperature_bodyterm_rhs ( 
	FV_DiscreteField const* FF,
	LA_Vector* VEC_rhs )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( 
   "REG_HeatEquation:: assemble_temperature_bodyterm_rhs" ) ;

   if ( my_rank == is_master ) cout << "Temperature body term rhs "
   	<< endl;

   // Parameters
   size_t nb_comps = FF->nb_components() ;
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
REG_HeatEquation:: error_with_analytical_solution (
	FV_DiscreteField const* FF,
	FV_DiscreteField* FF_ERROR )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( 
   "REG_HeatEquation:: error_with_analytical_solution" ) ;

   // Parameters
   size_t nb_comps = FF->nb_components() ;
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
     error_L2 = MAC::sqrt( error_L2 ) ;
     if ( my_rank == 0 )
       cout << "L2 Error, field " << FF->name() << ", component " << comp
     	<< " = " << error_L2 << endl;      	       	   
   }

}




//---------------------------------------------------------------------------
double
REG_HeatEquation:: compute_L2Norm_solution (
	FV_DiscreteField const* FF, size_t const& comp )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( 
   "REG_HeatEquation:: compute_L2Norm_solution" ) ;
   
   double L2Norm = 0. ;
   size_t k = 0 ;   

   // Get nb of local dof
   size_t_vector local_dof_number( dim, 0 );
   for (size_t l=0;l<dim;++l) 
     local_dof_number(l) = FF->get_local_nb_dof( comp, l ) ;
	
   // Compute L2 norm
   for (size_t i=0;i<local_dof_number(0);++i)     
     for (size_t j=0;j<local_dof_number(1);++j) 
       if ( dim == 2 )
       {
	 if ( FF->DOF_is_unknown_handled_by_proc( i, j, k, comp ) )
	   L2Norm += MAC::sqr( FF->DOF_value( i, j, k, comp, 0 ) )
	     	* FF->get_cell_measure( i, j, k, comp ) ;
       }	 
       else
       {
	 for (k=0;k<local_dof_number(2);++k)
	   if ( FF->DOF_is_unknown_handled_by_proc( i, j, k, comp ) )
	     L2Norm += MAC::sqr( FF->DOF_value( i, j, k, comp, 0 ) )
	     	* FF->get_cell_measure( i, j, k, comp ) ;
       }	     
     
   L2Norm = pelCOMM->sum( L2Norm ) ;
   L2Norm = MAC::sqrt( L2Norm ) ;
   if ( my_rank == 0 )
     cout << "L2 norm of field " << FF->name() << ", component " << comp
     	<< " = " << MAC::doubleToString( ios::scientific, 10, L2Norm ) << endl;
   
   return( L2Norm ); 
   
}
