#include <FV_StreamFunction.hh>
#include <FV_DomainAndFields.hh>
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
#include <MAC_Vector.hh>
#include <MAC_Communicator.hh>
#include <MAC_Exec.hh>
#include <LA_Matrix.hh>
#include <LA_Vector.hh>
#include <LA_Scatter.hh>
#include <LA_SeqVector.hh>
#include <LA_SeqMatrix.hh>
#include <LA_Solver.hh>
#include <size_t_vector.hh>
#include <math.h>
using std::cout ; 
using std::endl ;
using std::ostringstream ; 

FV_StreamFunction const* FV_StreamFunction::PROTOTYPE = 
	new FV_StreamFunction() ;

struct FV_StreamFunction_ERROR
{
   static void n1( std::string const& name ) ;   
} ;


//---------------------------------------------------------------------------
FV_StreamFunction:: FV_StreamFunction( void )
//--------------------------------------------------------------------------
   : FV_OneStepIteration( "FV_StreamFunction" ) 
{
   MAC_LABEL( "FV_StreamFunction:: FV_StreamFunction" ) ;
}




//---------------------------------------------------------------------------
FV_StreamFunction*
FV_StreamFunction:: create_replica( MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer* exp ) const
//---------------------------------------------------------------------------
{
   MAC_LABEL( "FV_StreamFunction:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, dom, exp ) ) ;

   FV_StreamFunction* result = 
                        new FV_StreamFunction( a_owner, dom, exp ) ;

   MAC_CHECK( create_replica_POST( result, a_owner, dom, exp ) ) ;
   return( result ) ;
   
}




//---------------------------------------------------------------------------
FV_StreamFunction:: FV_StreamFunction( MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer const* exp )
//---------------------------------------------------------------------------
   : FV_OneStepIteration( a_owner, dom, exp ) 
   , UU( 0 ) 
   , SF( 0 )
   , velocity_level( 0 ) 
   , velocity_name( "velocity" )
   , streamfunction_name( "streamfunction" )       
{
   MAC_LABEL( "FV_StreamFunction:: FV_StreamFunction" ) ;

   // Fields
   if ( exp->has_entry( "velocity" ) )
     velocity_name = exp->string_data( "velocity" ) ;
     
   if ( exp->has_entry( "stream_function" ) )
     streamfunction_name = exp->string_data( "stream_function" ) ;     

   if ( !dom->has_discrete_field( velocity_name ) )
     FV_StreamFunction_ERROR::n1( velocity_name ) ;
     
   if ( !dom->has_discrete_field( streamfunction_name ) )
     FV_StreamFunction_ERROR::n1( streamfunction_name ) ;     

   UU = dom->discrete_field( velocity_name ) ;
   SF = dom->discrete_field( streamfunction_name ) ;

   MAC_ASSERT( UU->primary_grid()->nb_space_dimensions() == 2 ) ;
   MAC_ASSERT( SF->discretization_type() == "vertex" ) ;
   MAC_ASSERT( SF->nb_components() == 1 ) ;
   MAC_ASSERT( SF->storage_depth() == 1 ) ;      
   MAC_ASSERT( UU->discretization_type() == "staggered" ) ; 
   MAC_ASSERT( UU->nb_components() == 2 ) ;             
      
   // MPI
   macCOMM = MAC_Exec::communicator();	   		
   my_rank = macCOMM->rank();
   nb_procs = macCOMM->nb_ranks();
   is_master = 0;  
   
   // Matrix system
   MAC_ModuleExplorer* se = 
	exp->create_subexplorer( 0, "FV_StreamFunctionSystem" ) ;    
   MAT_SF_Laplacian = LA_Matrix::make( this,
         se->create_subexplorer( this, 
	 "MAT_StreamFunction" ) ) ;
   VEC_SF_rhs = MAT_SF_Laplacian->create_vector( this ) ;
   VEC_SF = MAT_SF_Laplacian->create_vector( this ) ;
   SF_LOC = LA_SeqVector::create( this, 0 ) ;
   SF_NUM = FV_SystemNumbering::create( this, SF ) ;
   SOLVER_SF = LA_Solver::make( this,
               se->create_subexplorer( this,
	       "SOLVER_StreamFunction" ) ) ;
   se->destroy() ;
   se = 0 ;
   
   size_t sf_glob = SF->nb_global_unknowns() ;
   size_t sf_loc = SF->nb_local_unknowns() ; 
   
   MAT_SF_Laplacian->re_initialize( sf_glob, sf_glob ) ;
   VEC_SF_rhs->re_initialize( sf_glob ) ; 
   VEC_SF->re_initialize( sf_glob ) ;
   SF_LOC->re_initialize( sf_loc ) ; 
   SF_NUM->define_scatter( VEC_SF ) ;        
}




//---------------------------------------------------------------------------
FV_StreamFunction:: ~FV_StreamFunction( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "FV_StreamFunction:: ~FV_StreamFunction" ) ;

}




//---------------------------------------------------------------------------
void
FV_StreamFunction:: do_one_inner_iteration( FV_TimeIterator const* t_it ) 
//---------------------------------------------------------------------------
{
   MAC_LABEL( "FV_StreamFunction:: do_one_inner_iteration" ) ;
   MAC_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;	      
}




//---------------------------------------------------------------------------
void
FV_StreamFunction:: do_before_time_stepping( FV_TimeIterator const* t_it, 
      	std::string const& basename )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "FV_StreamFunction:: do_before_time_stepping" ) ;
   
   start_total_timer( "FV_StreamFunction:: do_before_time_stepping" ) ;

   // Assemble stream function laplacian matrix
   assemble_StreamFunctionLaplacian_matrix() ; 

   stop_total_timer() ;  
}




//---------------------------------------------------------------------------
void
FV_StreamFunction:: do_after_time_stepping( void )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "FV_StreamFunction:: do_after_time_stepping" ) ;  
}




//---------------------------------------------------------------------------
void
FV_StreamFunction:: do_before_inner_iterations_stage( 
	FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "FV_StreamFunction:: do_before_inner_iterations_stage" ) ;  
}




//---------------------------------------------------------------------------
void
FV_StreamFunction:: do_after_inner_iterations_stage( 
	FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "FV_StreamFunction:: do_after_inner_iterations_stage" ) ;  
}
  



//---------------------------------------------------------------------------
void
FV_StreamFunction:: do_additional_savings( FV_TimeIterator const* t_it,
      	int const& cycleNumber )
//---------------------------------------------------------------------------
{
   MAC_LABEL( "FV_StreamFunction:: do_additional_savings" ) ;

   start_total_timer( "FV_StreamFunction:: do_additional_savings" ) ;

   if ( my_rank == is_master ) 
     FV::out() << "         Solve StreamFunction problem" << endl;

   // Compute rhs
   assemble_StreamFunctionRHS() ;

   // Solve linear system
   SOLVER_SF->solve( VEC_SF_rhs, VEC_SF ) ;     

   // Copy back to field
   SF_NUM->scatter()->get( VEC_SF, SF_LOC ) ; 
   SF->update_free_DOFs_value( 0, SF_LOC ) ;      

   stop_total_timer() ;    
}




//---------------------------------------------------------------------------
void
FV_StreamFunction:: assemble_StreamFunctionLaplacian_matrix ( void )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( 
   "FV_Stokes:: assemble_StreamFunctionLaplacian_matrix" ) ;

   if ( my_rank == is_master ) 
     FV::out() << "         Assemble stream function laplacian matrix" << endl;

   // Parameters
   size_t dim = 2 ;
   double xC, xR, xL, yC, yT, yB, dxr, dxl, dyt, dyb, dxC, dyC ;
   double ac, arx, alx, aty, aby ;
   size_t center_pos_in_matrix = 0 ;
   
   // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l) min_unknown_index(l) = 
       	SF->get_min_index_unknown_handled_by_proc( 0, l ) ;
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l) max_unknown_index(l) = 
       	SF->get_max_index_unknown_handled_by_proc( 0, l ) ;
	
   // Perform assembling
   for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i)
   {          
     xC = SF->get_DOF_coordinate( i, 0, 0 ) ;
     xR = SF->get_DOF_coordinate_Assembling( i+1, 0, 0 ) ;
     xL = SF->get_DOF_coordinate_Assembling( i-1, 0, 0 ) ;
     dxr = xR - xC ; 
     dxl = xC - xL ; 
     dxC = SF->get_cell_size( i, 0, 0 ) ;      
     for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j)
     {
       yC = SF->get_DOF_coordinate( j, 0, 1 ) ; 
       yT = SF->get_DOF_coordinate_Assembling( j+1, 0, 1 ) ;
       yB = SF->get_DOF_coordinate_Assembling( j-1, 0, 1 ) ;
       dyt = yT - yC ;
       dyb = yC - yB ;
       dyC = SF->get_cell_size( j, 0, 1 ) ; 
 
       center_pos_in_matrix = SF->DOF_global_number( i, j, 0, 0 );
	   
       // Right (X)
       arx = - dyC / dxr ;
       arx = one_DOF_StreamFunctionLaplacian( i+1, j, center_pos_in_matrix, 
       		arx ) ;
	 
       // Left (X)
       alx = - dyC / dxl ;
       alx = one_DOF_StreamFunctionLaplacian( i-1, j, center_pos_in_matrix, 
       		alx ) ;
	   
       // Top (Y)
       aty = - dxC / dyt ;
       aty = one_DOF_StreamFunctionLaplacian( i, j+1, center_pos_in_matrix, 
       		aty ) ;
	   
       // Bottom (Y)
       aby = - dxC / dyb ;
       aby = one_DOF_StreamFunctionLaplacian( i, j-1, center_pos_in_matrix, 
       		aby ) ;
	   
       // Center
       ac = - arx - alx - aty - aby ;
       MAT_SF_Laplacian->add_to_item( center_pos_in_matrix, 
       		center_pos_in_matrix, ac ) ;

     }      	       	   
   }

   // Synchronize
   MAT_SF_Laplacian->synchronize() ;     
   SOLVER_SF->set_matrix( MAT_SF_Laplacian ) ;
}




//---------------------------------------------------------------------------
double
FV_StreamFunction:: one_DOF_StreamFunctionLaplacian( 
	int i, int j, 
	size_t center_pos_in_matrix, 
	double ai )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "FV_Stokes:: one_DOF_StreamFunctionLaplacian" ) ;

   size_t pos_in_matrix = SF->DOF_global_number( i, j, 0, 0 );
   if ( SF->DOF_is_unknown( i, j, 0, 0 ) )
     MAT_SF_Laplacian->add_to_item( center_pos_in_matrix, pos_in_matrix, ai ) ;
   else if ( SF->DOF_has_imposed_Dirichlet_value( i, j, 0, 0 ) )
   {
     double dirichlet_value = SF->DOF_value( i, j, 0, 0, 0 ) ;
     pair< size_t, double > ppp( center_pos_in_matrix, - ai * dirichlet_value );
     SF_DirichletBCValues.push_back( ppp ) ;
   }
   else ai = 0. ;
 
   return ( ai ) ;   
}




//---------------------------------------------------------------------------
void 
FV_StreamFunction:: assemble_StreamFunctionRHS( void )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "FV_Stokes:: assemble_StreamFunctionRHS" ) ;

   // Parameters
   size_t dim = 2 ;
   double dxC, dyC, xR, xL, yT, yB, vorticity_source_term = 0. ;
   size_t center_pos_in_matrix = 0 ;
   FV_SHIFT_TRIPLET shift = SF->shift_vertexToStaggered() ;
   
   // Set rhs to zero
   VEC_SF_rhs->nullify() ;
   
   // Add Dirichlet BC values
   list< pair< size_t, double > >::iterator il;
   for (il=SF_DirichletBCValues.begin();il!=SF_DirichletBCValues.end();il++)
     VEC_SF_rhs->add_to_item( il->first, il->second ) ;
     
   // Add vorticity source term
   // Get local min and max indices
   size_t_vector min_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l) min_unknown_index(l) = 
       	SF->get_min_index_unknown_handled_by_proc( 0, l ) ;
   size_t_vector max_unknown_index(dim,0);
   for (size_t l=0;l<dim;++l) max_unknown_index(l) = 
       	SF->get_max_index_unknown_handled_by_proc( 0, l ) ;
	
   // Perform assembling
   for (size_t i=min_unknown_index(0);i<=max_unknown_index(0);++i)
   {          
     dxC = SF->get_cell_size( i, 0, 0 ) ; 
     xR = UU->get_DOF_coordinate_Assembling( i+shift.i, 1, 0 ) ;
     xL = UU->get_DOF_coordinate_Assembling( i+shift.i-1, 1, 0 ) ;     
     for (size_t j=min_unknown_index(1);j<=max_unknown_index(1);++j)
     {
       dyC = SF->get_cell_size( j, 0, 1 ) ;  
       center_pos_in_matrix = SF->DOF_global_number( i, j, 0, 0 );
       yT = UU->get_DOF_coordinate_Assembling( j+shift.j, 0, 1 ) ;
       yB = UU->get_DOF_coordinate_Assembling( j+shift.j-1, 0, 1 ) ;
       
       vorticity_source_term = 
       	( UU->DOF_value( i+shift.i, j, 0, 1, velocity_level ) 
       	- UU->DOF_value( i+shift.i-1, j, 0, 1, velocity_level ) ) 
	/ ( xR - xL )
	- ( UU->DOF_value( i, j+shift.j, 0, 0, velocity_level ) 
       	- UU->DOF_value( i, j+shift.j-1, 0, 0, velocity_level ) ) 
	/ ( yT - yB ) ;

       VEC_SF_rhs->add_to_item( center_pos_in_matrix, 
       		vorticity_source_term * dxC * dyC ) ;
     }      	       	   
   } 
   
   VEC_SF_rhs->synchronize() ;
}	




//----------------------------------------------------------------------------
void
FV_StreamFunction:: do_more_post_processing(
        FV_DomainAndFields* dom,
        MAC_ModuleExplorer const* exp )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_StreamFunction:: do_more_post_processing" ) ;

   if ( my_rank == is_master ) 
     FV::out() << "         Nothing to do in FV_StreamFunction::"
     	"do_more_post_processing " << endl << endl;
}




//internal--------------------------------------------------------------
void
FV_StreamFunction_ERROR:: n1( std::string const& name )
//internal--------------------------------------------------------------
{
   ostringstream mesg ;
   mesg << "Field \"" << name << "\" is missing in MODULE interior_fields "
   	<< "to build FV_StreamFunction" << endl ;
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}
