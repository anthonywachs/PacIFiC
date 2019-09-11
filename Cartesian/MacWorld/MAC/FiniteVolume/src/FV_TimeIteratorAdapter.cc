#include <FV_TimeIteratorAdapter.hh>

#include <FV_DomainAndFields.hh>

#include <MAC.hh>
#include <MAC_Communicator.hh>
#include <MAC_Error.hh>
#include <MAC_Exec.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_ObjectRegister.hh>
#include <MAC_ObjectReader.hh>
#include <MAC_ObjectWriter.hh>
#include <MAC_Root.hh>
#include <MAC_assertions.hh>

#include <FV_TimeIterator.hh>

FV_TimeIteratorAdapter const* FV_TimeIteratorAdapter::DEFAULT_PROTOTYPE = 0 ;

//-------------------------------------------------------------------------
FV_TimeIteratorAdapter*
FV_TimeIteratorAdapter:: make( FV_TimeIterator* t_it,
                         FV_DomainAndFields const* dom,
                         MAC_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIteratorAdapter:: make(single domain)" ) ;
   MAC_CHECK_PRE( t_it != 0 ) ;
   MAC_CHECK_PRE( !t_it->is_started() ) ;
   MAC_CHECK_PRE( dom != 0 ) ;
   MAC_CHECK_PRE( exp != 0 ) ;

   std::string const& name = exp->string_data( "concrete_name" ) ;
   FV_TimeIteratorAdapter const* proto =
      static_cast<FV_TimeIteratorAdapter const*>(
                                            plugins_map()->item( name ) ) ;
      
   FV_TimeIteratorAdapter* result =
                            proto->create_replica( t_it, dom, exp ) ;
   
   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == t_it ) ;
   MAC_CHECK_POST( result->time_iterator() == t_it ) ;
   return( result ) ;
}




//-------------------------------------------------------------------------
FV_TimeIteratorAdapter*
FV_TimeIteratorAdapter:: make( FV_TimeIterator* t_it,
                         MAC_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIteratorAdapter:: make(multi domain)" ) ;
   MAC_CHECK_PRE( t_it != 0 ) ;
   MAC_CHECK_PRE( !t_it->is_started() ) ;
   MAC_CHECK_PRE( exp != 0 ) ;

   std::string const& name = exp->string_data( "concrete_name" ) ;
   FV_TimeIteratorAdapter const* proto =
      static_cast<FV_TimeIteratorAdapter const*>(
                                            plugins_map()->item( name ) ) ;
      
   FV_TimeIteratorAdapter* result =
                          proto->create_replica( t_it, exp ) ;
   
   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == t_it ) ;
   MAC_CHECK_POST( result->time_iterator() == t_it ) ;
   return( result ) ;
}




//-------------------------------------------------------------------------
FV_TimeIteratorAdapter*
FV_TimeIteratorAdapter:: make_default( FV_TimeIterator* t_it,
                                 FV_DomainAndFields const* dom )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIteratorAdapter:: make_default(single domain)" ) ;
   MAC_CHECK_PRE( t_it != 0 ) ;
   MAC_CHECK_PRE( !t_it->is_started() ) ;
   MAC_CHECK_PRE( dom != 0 ) ;

   FV_TimeIteratorAdapter* result = 0 ;
   if( DEFAULT_PROTOTYPE != 0 )
   {
      result = DEFAULT_PROTOTYPE->create_replica( t_it, dom ) ;
   }
   else
   {
      result = new FV_TimeIteratorAdapter( t_it ) ;
   }
   
   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == t_it ) ;
   MAC_CHECK_POST( result->time_iterator() == t_it ) ;
   return( result ) ;
}




//-------------------------------------------------------------------------
FV_TimeIteratorAdapter*
FV_TimeIteratorAdapter:: make_default( FV_TimeIterator* t_it )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIteratorAdapter:: make_default(multi domain)" ) ;
   MAC_CHECK_PRE( t_it != 0 ) ;
   MAC_CHECK_PRE( !t_it->is_started() ) ;

   FV_TimeIteratorAdapter* result = 0 ;
   if( DEFAULT_PROTOTYPE != 0 )
   {
      result = DEFAULT_PROTOTYPE->create_replica( t_it ) ;
   }
   else
   {
      result = new FV_TimeIteratorAdapter( t_it ) ;
   }
   
   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == t_it ) ;
   MAC_CHECK_POST( result->time_iterator() == t_it ) ;
   return( result ) ;
}




//-------------------------------------------------------------------------
FV_TimeIteratorAdapter:: FV_TimeIteratorAdapter( FV_TimeIterator* t_it )
//-------------------------------------------------------------------------
   : MAC_Object( t_it )
   , IS_PROTO( false )
   , T_IT( t_it )
   , DT( t_it->time_step() )
   , FINISHED( false )
{
   MAC_LABEL( "FV_TimeIteratorAdapter:: FV_TimeIteratorAdapter" ) ;
   MAC_CHECK_POST( !is_a_prototype() ) ;
   MAC_CHECK_POST( owner() == t_it ) ;
   MAC_CHECK_POST( time_iterator() == t_it ) ;
}




//-------------------------------------------------------------------------
FV_TimeIteratorAdapter:: FV_TimeIteratorAdapter(
                                       std::string const& a_concrete_name )
//-------------------------------------------------------------------------
   : MAC_Object( plugins_map() )
   , IS_PROTO( true )
   , T_IT( 0 )
   , DT( MAC::bad_double() )
   , FINISHED( false )
{
   MAC_LABEL( "FV_TimeIteratorAdapter:: FV_TimeIteratorAdapter" ) ;

   plugins_map()->register_item( a_concrete_name, this ) ;

   MAC_CHECK_POST( is_a_prototype() ) ;
}




//-------------------------------------------------------------------------
FV_TimeIteratorAdapter:: FV_TimeIteratorAdapter( void )
//-------------------------------------------------------------------------
   : MAC_Object( plugins_map() )
   , IS_PROTO( true )
   , T_IT( 0 )
   , DT( MAC::bad_double() )
   , FINISHED( false )
{
   MAC_LABEL( "FV_TimeIteratorAdapter:: FV_TimeIteratorAdapter" ) ;

   static bool first = true ;
   if( !first )
   {
      MAC_Error::object()->raise_internal(
         "Several default FV_TimeIteratorAdapter objects are referenced" ) ;
   }
   first = false ;
   DEFAULT_PROTOTYPE = this ;
   
   MAC_CHECK_POST( is_a_prototype() ) ;
}




//-------------------------------------------------------------------------
FV_TimeIteratorAdapter:: ~FV_TimeIteratorAdapter( void )
//-------------------------------------------------------------------------
{
   if( this == DEFAULT_PROTOTYPE )
   {
      DEFAULT_PROTOTYPE = 0 ;
   }
}




//-------------------------------------------------------------------------
FV_TimeIterator const*
FV_TimeIteratorAdapter:: time_iterator( void ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIteratorAdapter:: time_iterator" ) ;
   FV_TimeIterator const* result = T_IT ;
   MAC_CHECK_POST( result != 0 ) ;
   return( result ) ;
}




//-------------------------------------------------------------------------
void
FV_TimeIteratorAdapter:: initialize_time_step( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIteratorAdapter:: initialize_time_step" ) ;
   MAC_CHECK_PRE( time_iterator()->is_started() ) ;
   MAC_CHECK_PRE( !time_iterator()->is_finished() ) ; 
   DT = T_IT->time_step() ;
   initialize_inner() ;
}




//-------------------------------------------------------------------------
void
FV_TimeIteratorAdapter:: adapt_time_iterator( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIteratorAdapter:: adapt_time_iterator" ) ;
   MAC_CHECK_PRE( time_iterator()->is_started() ) ;

   MAC_Communicator const* com = MAC_Exec::communicator() ;
   FINISHED = com->boolean_and( FINISHED ) ;
   
   if( FINISHED )
   {
      T_IT->finish_iterations() ;
   }
   else
   {
      bool restart = false ;
      double next_dt = DT ;
      define_parameters_for_next_iteration( FINISHED, restart, next_dt ) ;
      if( restart )
      {
         restart_iteration_with_new_time_step( next_dt ) ;
      }
      else if( !T_IT->time_step_is_fixed() )
      {
         DT = next_dt ;
         T_IT->set_next_time_step( DT ) ;
      }
      FINISHED = com->boolean_and( FINISHED ) ;
      if( FINISHED )
      {
         T_IT->finish_iterations() ;
      }
   }
   
   MAC_CHECK_POST(
      EQUIVALENT( iterations_are_finished(),
                  time_iterator()->is_finished() ) ) ;
   MAC_CHECK_POST(
      IMPLIES( !iterations_are_finished() &&
                             !time_iterator()->time_step_is_fixed(),
               time_iterator()->next_time_step() == next_time_step() ) ) ;
}




//-------------------------------------------------------------------------
void
FV_TimeIteratorAdapter:: propose_next_time_step( double dt )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIteratorAdapter:: propose_next_time_step" ) ;
   MAC_CHECK_PRE( time_iterator()->is_started() ) ;
   MAC_CHECK_PRE( !time_iterator()->is_finished() ) ;
   MAC_CHECK_PRE( dt > 0. ) ;
   DT = next_time_step_from_proposed_one( dt ) ;
}




//-------------------------------------------------------------------------
double
FV_TimeIteratorAdapter:: next_time_step( void ) const 
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIteratorAdapter:: next_time_step" ) ;
   double result = DT ;
   MAC_CHECK_POST( result > 0. ) ;
   return( result ) ;
}




//-------------------------------------------------------------------------
void
FV_TimeIteratorAdapter:: propose_to_finish_iterations( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIteratorAdapter:: propose_to_finish_iterations" ) ;
   MAC_CHECK_PRE( time_iterator()->is_started() ) ;
   MAC_CHECK_PRE( !time_iterator()->is_finished() ) ;
   
   FINISHED = true ;
   
   MAC_CHECK_POST( iterations_are_finished() ) ;
}




//-------------------------------------------------------------------------
bool
FV_TimeIteratorAdapter:: iterations_are_finished( void ) const
//-------------------------------------------------------------------------
{
   return( FINISHED ) ;
}




//-------------------------------------------------------------------------
void
FV_TimeIteratorAdapter:: set_time_iteration_failed( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIteratorAdapter:: set_time_iteration_failed" ) ;
   MAC_CHECK_PRE( set_time_iteration_failed_PRE() ) ;

   MAC_Error::object()->raise_plain(
      "*** Time step iteration error: iteration failed" ) ;
   
   MAC_CHECK_POST( set_time_iteration_failed_POST() ) ;
}




//----------------------------------------------------------------------
void
FV_TimeIteratorAdapter:: save_state( MAC_ObjectWriter* writer ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIteratorAdapter:: save_state" ) ;
   MAC_CHECK_PRE( save_state_PRE( writer ) ) ;

   T_IT->save_state( writer ) ;

   MAC_CHECK_POST( save_state_POST( writer ) ) ;
}




//----------------------------------------------------------------------
void
FV_TimeIteratorAdapter:: restore_state( MAC_ObjectReader* reader )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIteratorAdapter:: restore_state" ) ;
   MAC_CHECK_PRE( restore_state_PRE( reader ) ) ;

   T_IT->restore_state( reader ) ;
   
   MAC_CHECK_POST( restore_state_POST( reader ) ) ;
}




//----------------------------------------------------------------------
void
FV_TimeIteratorAdapter:: initialize_inner( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIteratorAdapter:: initialize_inner" ) ;
   MAC_CHECK_PRE( initialize_inner_PRE() ) ;
   MAC_CHECK_POST( initialize_inner_POST() ) ;
}




//-------------------------------------------------------------------------
void
FV_TimeIteratorAdapter:: define_parameters_for_next_iteration(
                           bool& finished, bool& restart, double& next_dt )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIteratorAdapter:: define_parameters_for_next_iteration" ) ;
   MAC_CHECK_PRE( define_parameters_for_next_iteration_PRE( finished, restart, next_dt ) ) ;
   MAC_CHECK_POST( define_parameters_for_next_iteration_POST( finished, restart, next_dt ) ) ;
}




//----------------------------------------------------------------------
double
FV_TimeIteratorAdapter:: next_time_step_from_proposed_one( double dt ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIteratorAdapter:: next_time_step_from_proposed_one" ) ;
   MAC_CHECK_PRE( next_time_step_from_proposed_one_PRE( dt ) ) ;
   double result = dt ;
   MAC_CHECK_POST( next_time_step_from_proposed_one_POST( result, dt ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
void
FV_TimeIteratorAdapter:: restart_iteration_with_new_time_step( double dt )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIteratorAdapter:: restart_iteration_with_new_time_step" ) ;
   MAC_CHECK_PRE( time_iterator()->is_started() ) ;
   MAC_CHECK_PRE( !time_iterator()->is_finished() ) ;
   MAC_CHECK_PRE( !time_iterator()->time_step_is_fixed() ) ;
   MAC_CHECK_PRE( !time_iterator()->just_went_back() ) ;
   MAC_CHECK_PRE( dt > 0. ) ;

   DT = dt ;
   T_IT->go_back() ;
   T_IT->set_next_time_step( DT ) ;

   MAC_CHECK_POST( time_iterator()->just_went_back() ) ;
   MAC_CHECK_POST( time_iterator()->next_time_step() == dt ) ;
   MAC_CHECK_POST( next_time_step() == dt ) ;
}




//-------------------------------------------------------------------------
FV_TimeIteratorAdapter*
FV_TimeIteratorAdapter:: create_replica( FV_TimeIterator* t_it,
                                   FV_DomainAndFields const* dom,
                                   MAC_ModuleExplorer const* exp ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIteratorAdapter:: create_replica" ) ;
   MAC_CHECK_PRE( create_replica_PRE( t_it, dom, exp ) ) ;

   MAC_Error::object()->raise_internal(
                                 "The function has not been overloaded" ) ;

   FV_TimeIteratorAdapter* result = 0 ;
   
   MAC_CHECK_POST( create_replica_POST( result, t_it, dom, exp ) ) ;
   return( result ) ;
}




//-------------------------------------------------------------------------
FV_TimeIteratorAdapter*
FV_TimeIteratorAdapter:: create_replica( FV_TimeIterator* t_it,
                                   MAC_ModuleExplorer const* exp ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIteratorAdapter:: create_replica" ) ;
   MAC_CHECK_PRE( create_replica_PRE( t_it, exp ) ) ;

   MAC_Error::object()->raise_internal(
                                 "The function has not been overloaded" ) ;

   FV_TimeIteratorAdapter* result = 0 ;
   
   MAC_CHECK_POST( create_replica_POST( result, t_it, exp ) ) ;
   return( result ) ;
}




//-------------------------------------------------------------------------
FV_TimeIteratorAdapter*
FV_TimeIteratorAdapter:: create_replica( FV_TimeIterator* t_it,
                                   FV_DomainAndFields const* dom ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIteratorAdapter:: create_replica" ) ;
   MAC_CHECK_PRE( create_replica_PRE( t_it, dom ) ) ;

   MAC_Error::object()->raise_internal(
                                 "The function has not been overloaded" ) ;

   FV_TimeIteratorAdapter* result = 0 ;
   
   MAC_CHECK_POST( create_replica_POST( result, t_it, dom ) ) ;
   return( result ) ;
}




//-------------------------------------------------------------------------
FV_TimeIteratorAdapter*
FV_TimeIteratorAdapter:: create_replica( FV_TimeIterator* t_it ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIteratorAdapter:: create_replica" ) ;
   MAC_CHECK_PRE( create_replica_PRE( t_it ) ) ;

   MAC_Error::object()->raise_internal(
                                 "The function has not been overloaded" ) ;

   FV_TimeIteratorAdapter* result = 0 ;
   
   MAC_CHECK_POST( create_replica_POST( result, t_it ) ) ;
   return( result ) ;
}




//-------------------------------------------------------------------------
bool
FV_TimeIteratorAdapter:: is_a_prototype( void ) const
//-------------------------------------------------------------------------
{
   return( IS_PROTO ) ;
}




//-------------------------------------------------------------------------
bool
FV_TimeIteratorAdapter:: invariant( void ) const
//-------------------------------------------------------------------------
{
   MAC_ASSERT( MAC_Object::invariant() ) ;
   return( true ) ;
}




//--------------------------------------------------------------------------
bool
FV_TimeIteratorAdapter:: set_time_iteration_failed_PRE( void ) const
//--------------------------------------------------------------------------
{
   MAC_ASSERT( time_iterator()->is_started() ) ;
   MAC_ASSERT( !time_iterator()->is_finished() ) ;
   return( true ) ;
}




//--------------------------------------------------------------------------
bool
FV_TimeIteratorAdapter:: set_time_iteration_failed_POST( void ) const
//--------------------------------------------------------------------------
{
   return( true ) ;
}




//--------------------------------------------------------------------------
bool
FV_TimeIteratorAdapter:: initialize_inner_PRE( void ) const 
//--------------------------------------------------------------------------
{
   MAC_ASSERT( time_iterator()->is_started() ) ;
   MAC_ASSERT( !time_iterator()->is_finished() ) ;
   MAC_ASSERT( next_time_step() == time_iterator()->time_step() ) ;
   return( true ) ;
}




//--------------------------------------------------------------------------
bool
FV_TimeIteratorAdapter:: initialize_inner_POST( void ) const 
//--------------------------------------------------------------------------
{
   return( true ) ;
}




//--------------------------------------------------------------------------
bool
FV_TimeIteratorAdapter:: define_parameters_for_next_iteration_PRE( 
                         bool finished, bool restart, double next_dt ) const
//--------------------------------------------------------------------------
{
   MAC_ASSERT( time_iterator()->is_started() ) ;
   MAC_ASSERT( !time_iterator()->is_finished() ) ;
   MAC_ASSERT( finished == false ) ;
   MAC_ASSERT( restart == false ) ;
   MAC_ASSERT( next_dt == next_time_step() ) ;
   return( true ) ;
}




//--------------------------------------------------------------------------
bool
FV_TimeIteratorAdapter:: define_parameters_for_next_iteration_POST(
                         bool finished, bool restart, double next_dt ) const
//--------------------------------------------------------------------------
{
   MAC_ASSERT( IMPLIES( finished, !restart ) ) ;
   MAC_ASSERT( IMPLIES( restart, !finished ) ) ;
   MAC_ASSERT( IMPLIES( time_iterator()->time_step_is_fixed(), !restart ) ) ;
   return( true ) ;
}




//--------------------------------------------------------------------------
bool
FV_TimeIteratorAdapter:: next_time_step_from_proposed_one_PRE(
                                                           double dt ) const
//--------------------------------------------------------------------------
{
   MAC_ASSERT( time_iterator()->is_started() ) ;
   MAC_ASSERT( !time_iterator()->is_finished() ) ;
   MAC_ASSERT( dt > 0. ) ;
   return( true ) ;
}




//--------------------------------------------------------------------------
bool
FV_TimeIteratorAdapter:: next_time_step_from_proposed_one_POST(
                                            double result, double dt ) const
//--------------------------------------------------------------------------
{
   MAC_ASSERT( result > 0. ) ;
   return( true ) ;
}




//--------------------------------------------------------------------------
bool
FV_TimeIteratorAdapter:: create_replica_PRE(
                                       FV_TimeIterator const* t_it,
                                       FV_DomainAndFields const* dom,
                                       MAC_ModuleExplorer const* exp ) const
//--------------------------------------------------------------------------
{
   MAC_ASSERT( is_a_prototype() ) ;
   MAC_ASSERT( t_it != 0 ) ;
   MAC_ASSERT( !t_it->is_started() ) ;
   MAC_ASSERT( dom != 0 ) ;
   MAC_ASSERT( exp != 0 ) ;
   return( true ) ;
}




//--------------------------------------------------------------------------
bool
FV_TimeIteratorAdapter:: create_replica_POST(
                                        FV_TimeIteratorAdapter const* result,
                                        FV_TimeIterator const* t_it,
                                        FV_DomainAndFields const* dom,
                                        MAC_ModuleExplorer const* exp ) const
//--------------------------------------------------------------------------
{
   MAC_ASSERT( result != 0 ) ;
   MAC_ASSERT( result->owner() == t_it ) ;
   MAC_ASSERT( !result->is_a_prototype() ) ;
   return( true ) ;
}




//--------------------------------------------------------------------------
bool
FV_TimeIteratorAdapter:: create_replica_PRE(
                                       FV_TimeIterator const* t_it,
                                       MAC_ModuleExplorer const* exp ) const
//--------------------------------------------------------------------------
{
   MAC_ASSERT( is_a_prototype() ) ;
   MAC_ASSERT( t_it != 0 ) ;
   MAC_ASSERT( !t_it->is_started() ) ;
   MAC_ASSERT( exp != 0 ) ;
   return( true ) ;
}




//--------------------------------------------------------------------------
bool
FV_TimeIteratorAdapter:: create_replica_POST(
                                        FV_TimeIteratorAdapter const* result,
                                        FV_TimeIterator const* t_it,
                                        MAC_ModuleExplorer const* exp ) const
//--------------------------------------------------------------------------
{
   MAC_ASSERT( result != 0 ) ;
   MAC_ASSERT( result->owner() == t_it ) ;
   MAC_ASSERT( !result->is_a_prototype() ) ;
   return( true ) ;
}




//--------------------------------------------------------------------------
bool
FV_TimeIteratorAdapter:: create_replica_PRE(
                                       FV_TimeIterator const* t_it,
                                       FV_DomainAndFields const* dom ) const
//--------------------------------------------------------------------------
{
   MAC_ASSERT( is_a_prototype() ) ;
   MAC_ASSERT( t_it != 0 ) ;
   MAC_ASSERT( !t_it->is_started() ) ;
   MAC_ASSERT( dom != 0 ) ;
   return( true ) ;
}




//--------------------------------------------------------------------------
bool
FV_TimeIteratorAdapter:: create_replica_POST(
                                        FV_TimeIteratorAdapter const* result,
                                        FV_TimeIterator const* t_it,
                                        FV_DomainAndFields const* dom ) const
//--------------------------------------------------------------------------
{
   MAC_ASSERT( result != 0 ) ;
   MAC_ASSERT( result->owner() == t_it ) ;
   MAC_ASSERT( !result->is_a_prototype() ) ;
   return( true ) ;
}




//--------------------------------------------------------------------------
bool
FV_TimeIteratorAdapter:: create_replica_PRE(
                                       FV_TimeIterator const* t_it ) const
//--------------------------------------------------------------------------
{
   MAC_ASSERT( is_a_prototype() ) ;
   MAC_ASSERT( t_it != 0 ) ;
   MAC_ASSERT( !t_it->is_started() ) ;
   return( true ) ;
}




//--------------------------------------------------------------------------
bool
FV_TimeIteratorAdapter:: create_replica_POST(
                                        FV_TimeIteratorAdapter const* result,
                                        FV_TimeIterator const* t_it ) const
//--------------------------------------------------------------------------
{
   MAC_ASSERT( result != 0 ) ;
   MAC_ASSERT( result->owner() == t_it ) ;
   MAC_ASSERT( !result->is_a_prototype() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
MAC_ObjectRegister*
FV_TimeIteratorAdapter:: plugins_map( void )
//----------------------------------------------------------------------
{
   static MAC_ObjectRegister* result =
      MAC_ObjectRegister::create( MAC_Root::object(),
                                  "FV_TimeIteratorAdapter descendant" ) ;
   return( result ) ;
}
