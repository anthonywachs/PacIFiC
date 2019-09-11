#include <LA_Solver.hh>

#include <MAC.hh>
#include <MAC_Communicator.hh>
#include <MAC_DistributedPartition.hh>
#include <MAC_Error.hh>
#include <MAC_Exec.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_ObjectRegister.hh>
#include <MAC_Root.hh>
#include <MAC_System.hh>
#include <MAC_Timer.hh>

//#include <LA_PreconditionedSolver.hh>
#include <LA_SeqMatrix.hh>
#include <LA_SeqVector.hh>

#include <iostream>
#include <sstream>

using std::endl ;
using std::string ;
using std::ostringstream ;

//----------------------------------------------------------------------
bool
LA_Solver:: is_registered( std::string const solver_name )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Solver:: is_registered" ) ;
   MAC_CHECK_PRE( !solver_name.empty() ) ;
   return( plugins_map()->has( solver_name ) ) ;
}

//----------------------------------------------------------------------
LA_Solver*
LA_Solver:: make( MAC_Object* a_owner,
                  MAC_ModuleExplorer const* exp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Solver:: make" ) ;
   MAC_CHECK_PRE( exp != 0 ) ;

   std::string name = exp->string_data( "concrete_name" ) ;
   LA_Solver const* proto =
      static_cast<LA_Solver const*>(
         plugins_map()->item( name ) ) ;
   MAC_ASSERT( proto->is_a_prototype() ) ;

   LA_Solver* result = proto->create_replica( a_owner, exp ) ;
   result->SOLVER_NAME = proto->SOLVER_NAME ;
   if( exp->has_entry( "name" ) )
      result->SOLVER_NAME = exp->string_data( "name" ) ;
   if( exp->has_entry( "verbose" ) )
      result->VERBOSE = exp->bool_data( "verbose" ) ;
   if( exp->has_entry( "save_matrix_at_iteration" ) )
      result->SAVE_MATRIX_ITER_NB = exp->int_data( "save_matrix_at_iteration" ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->stop_on_error() ) ;
   MAC_CHECK_POST( result->zero_initial_guess() ) ;
   MAC_CHECK_POST( !result->matrix_is_set() ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
LA_Solver:: LA_Solver( MAC_Object* a_owner )
//----------------------------------------------------------------------
   : MAC_Object( a_owner )
   , IS_PROTO( false )
   , SOLVER_NAME( "" )
   , ITERATIVE( false )
   , ZERO_INIT( true )
   , STOP( true )
   , VERBOSE( false )
   , SIZE( MAC::bad_index() )
   , NB_ITER( MAC::bad_index() )
   , SOL_ACHIEVED( false )
   , MATRIX( 0 )
   , SAVE_MATRIX_ITER_NB( -1 )
{
   MAC_LABEL( "LA_Solver:: LA_Solver" ) ;
}

//----------------------------------------------------------------------
LA_Solver:: LA_Solver( MAC_Object* a_owner, LA_Solver const* other )
//----------------------------------------------------------------------
   : MAC_Object( a_owner )
   , IS_PROTO( false )
   , SOLVER_NAME( "copy of "+other->SOLVER_NAME )
   , ITERATIVE( other->ITERATIVE )
   , ZERO_INIT( other->ZERO_INIT )
   , STOP( other->STOP )
   , VERBOSE( other->VERBOSE )
   , SIZE( MAC::bad_index() )
   , NB_ITER( MAC::bad_index() )
   , SOL_ACHIEVED( false )
   , MATRIX( 0 )
   , SAVE_MATRIX_ITER_NB( other->SAVE_MATRIX_ITER_NB )
{
   MAC_LABEL( "LA_Solver:: LA_Solver" ) ;
}

//----------------------------------------------------------------------
LA_Solver:: LA_Solver( std::string const& a_name )
//----------------------------------------------------------------------
   : MAC_Object( plugins_map() )
   , IS_PROTO( true )
   , SOLVER_NAME( a_name )
   , ITERATIVE( false )
   , ZERO_INIT( true )
   , STOP( true )
   , VERBOSE( false )
   , SIZE( MAC::bad_index() )
   , NB_ITER( MAC::bad_index() )
   , SOL_ACHIEVED( false )
   , MATRIX( 0 )
   , SAVE_MATRIX_ITER_NB( -1 )
{
   MAC_LABEL( "LA_Solver:: LA_Solver" ) ;
   plugins_map()->register_item( a_name, this ) ;
}

//----------------------------------------------------------------------
LA_Solver:: ~LA_Solver( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Solver:: ~LA_Solver" ) ;
   //MAC_CHECK( is_a_prototype() || ! matrix_is_set() ) ;
}

//----------------------------------------------------------------------
void
LA_Solver:: set_stop_on_error( bool stop )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Solver:: set_stop_on_error" ) ;
   MAC_CHECK_PRE( !is_a_prototype() ) ;

   STOP = stop ;

   MAC_CHECK_POST( stop_on_error() == stop ) ;
}

//----------------------------------------------------------------------
bool
LA_Solver:: stop_on_error( void ) const
//----------------------------------------------------------------------
{
   return( STOP ) ;
}

//----------------------------------------------------------------------
bool
LA_Solver:: matrix_is_set( void ) const
//----------------------------------------------------------------------
{
   return( MATRIX!=0 ) ;
}


//----------------------------------------------------------------------
LA_Matrix const*
LA_Solver:: matrix( void ) const
//----------------------------------------------------------------------
{
   MAC_CHECK_PRE( matrix_is_set() ) ;
   return( MATRIX ) ;
}

//----------------------------------------------------------------------
LA_Implementation const*
LA_Solver:: matrix_implementation( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Solver:: matrix_implementation" ) ;
   MAC_CHECK_PRE( matrix_is_set() ) ;
   return( MATRIX->implementation() ) ;
}

//----------------------------------------------------------------------
bool
LA_Solver:: is_verbose( void ) const
//----------------------------------------------------------------------
{
   return( VERBOSE ) ;
}

//----------------------------------------------------------------------
size_t
LA_Solver:: size( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Solver:: size" ) ;
   MAC_CHECK_PRE( matrix_is_set() ) ;

   return( SIZE ) ;
}

//----------------------------------------------------------------------
bool
LA_Solver:: solution_is_achieved( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Solver:: solution_is_achieved" ) ;

   return( SOL_ACHIEVED ) ;
}

//----------------------------------------------------------------------
size_t
LA_Solver:: nb_iterations_achieved( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Solver:: nb_iterations_achieved" ) ;
   MAC_CHECK_PRE( solution_is_achieved() ) ;

   size_t result = NB_ITER ;

   MAC_CHECK_POST( IMPLIES( !is_iterative(), result==1 ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
bool
LA_Solver:: is_iterative( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Solver:: is_iterative" ) ;
   return( ITERATIVE ) ;
}

//----------------------------------------------------------------------
void
LA_Solver:: set_iterative( bool iterative )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Solver:: set_iterative" ) ;
   ITERATIVE = iterative ;
   MAC_CHECK_POST( is_iterative() == iterative ) ;
}

//----------------------------------------------------------------------
void
LA_Solver:: set_initial_guess_nonzero( bool flg )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Solver:: set_initial_guess_nonzero" ) ;

   if( is_iterative() ) ZERO_INIT = !flg ;

   MAC_CHECK_POST( IMPLIES( is_iterative(),
                            zero_initial_guess() == !flg ) ) ;
}

//----------------------------------------------------------------------
bool
LA_Solver:: zero_initial_guess( void ) const
//----------------------------------------------------------------------
{
   return( ZERO_INIT ) ;
}

//----------------------------------------------------------------------
void
LA_Solver:: set_matrix( LA_Matrix const* mat )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Solver:: set_matrix" ) ;
   MAC_CHECK_PRE( mat != 0 ) ;
   MAC_CHECK_PRE( mat->nb_rows() > 0 ) ;
   MAC_CHECK_PRE( mat->nb_rows() == mat->nb_cols() ) ;
   MAC_CHECK_PRE( !matrix_is_set() ) ;
   MAC_CHECK_PRE( mat->is_synchronized() ) ;

   if( is_verbose() )
   {
      MAC::out() << SOLVER_NAME
                 << ": setting matrix of size "
                 << mat->nb_rows() << std::endl ;
   }

   if( SAVE_MATRIX_ITER_NB==0 )
   {
      std::ostringstream os ;
      static size_t save_index = 0 ;

      os << SOLVER_NAME << save_index++ << ".mtx" ;
      MAC::out() << "Write matrix to file " << os.str() << std::endl ;
      mat->writeMM( os.str() ) ;
   }
   SAVE_MATRIX_ITER_NB-- ;

   SOL_ACHIEVED = false ;
   bool ok = false ;
   set_matrix_self( mat, ok, false ) ;
   if( ok )
   {
      MATRIX = mat ;
      SIZE = mat->nb_rows() ;
   }
   if( !ok && stop_on_error() )
   {
      MAC_Error::object()->raise_plain(
         "*** LA_Solver error\n"
         "    set_matrix fails" ) ;
   }

   if( is_verbose() )
   {
      MAC::out() << SOLVER_NAME << ": usage memory after matrix setting (Mo) : "
                 << MAC_System::used_memory() /1024 / 1024 << std::endl ;
   }
   MAC_CHECK_POST( !solution_is_achieved() ) ;
   MAC_CHECK_POST( IMPLIES( stop_on_error(),
                            matrix_is_set() ) ) ;
   MAC_CHECK_POST( IMPLIES( matrix_is_set(),
                            matrix_implementation() == mat->implementation() ) ) ;
}

//----------------------------------------------------------------------
void
LA_Solver:: reset_matrix( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Solver:: reset_matrix" ) ;
   MAC_CHECK_PRE( matrix_is_set() ) ;
   MAC_CHECK_PRE( size()==matrix()->nb_rows() ) ; // pas suffisant, mais c'est deja cela

   if( is_verbose() )
   {
      MAC::out() << SOLVER_NAME
                 << ": resetting matrix of size "
                 << MATRIX->nb_rows() << std::endl ;
   }

   if( SAVE_MATRIX_ITER_NB==0 )
   {
      std::ostringstream os ;
      os << SOLVER_NAME << hash_code() << ".mtx" ;
      MAC::out() << "Write matrix to file " << os.str() << std::endl ;
      MATRIX->writeMM( os.str() ) ;
   }
   SAVE_MATRIX_ITER_NB-- ;

   SOL_ACHIEVED = false ;
   bool ok = false ;
   set_matrix_self( MATRIX, ok, true ) ;
   if( !ok && stop_on_error() )
   {
      MAC_Error::object()->raise_plain(
         "*** LA_Solver error\n"
         "    set_matrix fails" ) ;
   }

   if( is_verbose() )
   {
      MAC::out() << SOLVER_NAME << ": usage memory after matrix setting (Mo) : "
                 << MAC_System::used_memory() /1024 / 1024 << std::endl ;
   }
   MAC_CHECK_POST( !solution_is_achieved() ) ;
   MAC_CHECK_POST( IMPLIES( stop_on_error(),
                            matrix_is_set() ) ) ;
}

//----------------------------------------------------------------------
void
LA_Solver:: unset_matrix( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Solver:: unset_matrix" ) ;
   MAC_CHECK_PRE( matrix_is_set() ) ;

   unset_matrix_self() ;
   MATRIX=0 ;

   MAC_CHECK_POST( !matrix_is_set() ) ;
}

//----------------------------------------------------------------------
void
LA_Solver:: solve( LA_Vector const* b, LA_Vector* x )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Solver:: solve" ) ;
   MAC_CHECK_PRE( matrix_is_set() ) ;
   MAC_CHECK_PRE( b != 0 ) ;
   MAC_CHECK_PRE( b->nb_rows() == size() ) ;
   MAC_CHECK_PRE( b->implementation() == matrix_implementation() ) ;
   MAC_CHECK_PRE( b->is_synchronized() ) ;
   MAC_CHECK_PRE(
       b->row_distribution()->is_compatible( matrix()->row_distribution() ) ) ;
   MAC_CHECK_PRE( x != 0 ) ;
   MAC_CHECK_PRE( x->nb_rows() == size() ) ;
   MAC_CHECK_PRE( x->implementation() == matrix_implementation() ) ;
   MAC_CHECK_PRE(
       x->row_distribution()->is_compatible( matrix()->row_distribution() ) ) ;

   if( is_verbose() )
      MAC::out() << SOLVER_NAME << " : solving system" <<std::endl ;

   solve_self( b, x, NB_ITER, SOL_ACHIEVED ) ;
   x->synchronize() ;

   if( !SOL_ACHIEVED && stop_on_error() )
   {
      MAC_Error::object()->raise_plain( "LA_Solver:: solve fail" ) ;
   }
   MAC_CHECK_POST( IMPLIES( stop_on_error(), solution_is_achieved() ) ) ;
   MAC_CHECK_POST( x->is_synchronized() ) ;
}

//----------------------------------------------------------------------
void
LA_Solver:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Solver:: print" ) ;
   std::string const s( indent_width, ' ' ) ;
   os << s << SOLVER_NAME ;
}

//----------------------------------------------------------------------
bool
LA_Solver:: is_a_prototype( void ) const
//----------------------------------------------------------------------
{
   return( IS_PROTO ) ;
}

//----------------------------------------------------------------------
void
LA_Solver:: raise_fatal_error_if_not_sequential( void )
//----------------------------------------------------------------------
{
   if( MAC_Exec::communicator()->nb_ranks() > 1 )
   {
      ostringstream mesg ;
      mesg << "*** " << type_name() << " error:" << std::endl ;
      mesg << "    parallel execution of multiple processes" << std::endl ;
      mesg << "    is not supported" ;
      MAC_Error::object()->raise_plain( mesg.str() ) ;
   }
}

//----------------------------------------------------------------------
bool
LA_Solver:: create_clone_POST( LA_Solver const* result,
                               MAC_Object* a_owner ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( MAC_Object::create_clone_POST( result, a_owner ) ) ;
   MAC_ASSERT( !result->is_a_prototype() ) ;
   MAC_ASSERT( result->is_iterative() == is_iterative() ) ;
   MAC_ASSERT( result->stop_on_error() == stop_on_error() ) ;
   MAC_ASSERT( result->is_verbose() == is_verbose() ) ;
   MAC_ASSERT( result->zero_initial_guess() == zero_initial_guess() ) ;
   MAC_ASSERT( !result->matrix_is_set() ) ;
   MAC_ASSERT( !result->solution_is_achieved() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Solver:: create_replica_PRE( MAC_Object* a_owner,
                                MAC_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( exp != 0 ) ;
   MAC_ASSERT( is_a_prototype() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Solver:: create_replica_POST( LA_Solver const* result,
                                 MAC_Object* a_owner,
                                 MAC_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( result != 0 ) ;
   MAC_ASSERT( !result->is_a_prototype() ) ;
   MAC_ASSERT( result->owner()==a_owner ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Solver:: unset_matrix_self_PRE( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( !is_a_prototype() ) ;
   MAC_ASSERT( matrix_is_set() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Solver:: set_matrix_self_PRE(
                         LA_Matrix const* mat, bool same_pattern ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( !is_a_prototype() ) ;
   MAC_ASSERT( EQUIVALENT( same_pattern, matrix_is_set() ) ) ;
   MAC_ASSERT( mat != 0 ) ;
   MAC_ASSERT( mat->nb_rows() > 0 ) ;
   MAC_ASSERT( mat->nb_rows() == mat->nb_cols() ) ;
   MAC_ASSERT( mat->is_synchronized() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Solver:: set_matrix_self_POST( LA_Matrix const* mat, bool ok ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Solver:: solve_self_PRE( LA_Vector const* b,
                            LA_Vector const* x ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( !is_a_prototype() ) ;
   MAC_ASSERT( matrix_is_set() ) ;
   MAC_ASSERT( b != 0 ) ;
   MAC_ASSERT( b->implementation() == matrix_implementation() ) ;
   MAC_ASSERT( b->nb_rows() == size() ) ;
   MAC_ASSERT( x != 0 ) ;
   MAC_ASSERT( x->implementation() == matrix_implementation() ) ;
   MAC_ASSERT( x->nb_rows() == size() ) ;
   return( true ) ;
}

//----------------------------------------------------------------------
bool
LA_Solver:: solve_self_POST( LA_Vector const* b,
                             LA_Vector const* x,
                             size_t nb_iter,
                             bool ok ) const
//----------------------------------------------------------------------
{
   return( true ) ;
}

//----------------------------------------------------------------------
MAC_ObjectRegister*
LA_Solver:: plugins_map( void )
//----------------------------------------------------------------------
{
   static MAC_ObjectRegister* result =
      MAC_ObjectRegister::create( MAC_Root::object(),
                                  "LA_Solver descendant" ) ;
   return( result ) ;
}

