#include <FV_OneStepIteration.hh>
#include <FV.hh>
#include <MAC.hh>
#include <MAC_Communicator.hh>
#include <MAC_Error.hh>
#include <MAC_Exec.hh>
#include <MAC_ListIdentity.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_ObjectRegister.hh>
#include <MAC_Root.hh>
#include <MAC_Timer.hh>
#include <MAC_assertions.hh>
#include <FV_TimeIteratorAdapter.hh>
#include <FV_TimeIterator.hh>
#include <iomanip>
#include <iostream>
#include <sstream>

using std::ios_base ; 
using std::setprecision ; 
using std::setw ;
using std::string ;
using std::endl ;
using std::ostringstream ;
using std::set ;

std::string FV_OneStepIteration::INDENT ;

std::map< std::string, MAC_Timer* > FV_OneStepIteration::ASS_TIMER ;
std::map< std::string, MAC_Timer* > FV_OneStepIteration::SOL_TIMER ;
std::map< std::string, MAC_Timer* > FV_OneStepIteration::TOT_TIMER ;

//----------------------------------------------------------------------
FV_OneStepIteration*
FV_OneStepIteration:: make( MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer* exp )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: make" ) ;
   MAC_CHECK_PRE( exp != 0 ) ;
   
   string name = exp->string_data( "concrete_name" ) ;
   FV_OneStepIteration const* proto =
      static_cast<FV_OneStepIteration const*>(
                                        plugins_map()->item( name ) ) ;
   MAC_ASSERT( proto->is_a_prototype() ) ;

   FV_OneStepIteration* result = proto->create_replica( a_owner, dom, exp ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
//??????????????????? beaucoup d'autres choses
   return( result ) ;
}




//----------------------------------------------------------------------------
FV_OneStepIteration:: FV_OneStepIteration( MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer const* exp )
//----------------------------------------------------------------------------
   : MAC_Object( a_owner )
   , IS_PROTO( false )
   , FAILURE( false )
   , VERBOSE_LEVEL( exp->has_entry( "verbose_level" ) ?
                                    exp->int_data( "verbose_level" ) : 2 )
   , DOMAIN( dom )
   , RANK( MAC_Exec::communicator()->rank() )				    
{
   MAC_LABEL( "FV_OneStepIteration:: FV_OneStepIteration" ) ;
}




//----------------------------------------------------------------------------
FV_OneStepIteration:: FV_OneStepIteration( std::string const& name )
//----------------------------------------------------------------------------
   : MAC_Object( plugins_map() )
   , IS_PROTO( true )
   , FAILURE( false )
   , VERBOSE_LEVEL( 0 )
   , DOMAIN( 0 )
   , RANK( 0 )	   
{
   MAC_LABEL( "FV_OneStepIteration:: FV_OneStepIteration" ) ;

   plugins_map()->register_item( name, this ) ;

   MAC_CHECK_POST( is_a_prototype() ) ;
}




//----------------------------------------------------------------------------
FV_OneStepIteration:: ~FV_OneStepIteration( void )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: ~FV_OneStepIteration" ) ;
}




//----------------------------------------------------------------------------
void
FV_OneStepIteration:: do_before_time_stepping( FV_TimeIterator const* t_it, 
      	std::string const& basename )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: do_before_time_stepping" ) ;
   MAC_CHECK_PRE( do_before_time_stepping_PRE( t_it ) ) ;
}




//----------------------------------------------------------------------------
void
FV_OneStepIteration:: do_before_inner_iterations_stage(
                                               FV_TimeIterator const* t_it )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: do_before_inner_iterations_stage" ) ;
   MAC_CHECK_PRE( do_before_inner_iterations_stage_PRE( t_it ) ) ;
   FAILURE = false ;
   MAC_CHECK_POST( do_before_inner_iterations_stage_POST( t_it ) ) ;
}




//----------------------------------------------------------------------------
void
FV_OneStepIteration:: do_inner_iterations_stage(
                                               FV_TimeIterator const* t_it )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: do_inner_iterations_stage" ) ;
   MAC_CHECK_PRE( do_inner_iterations_stage_PRE( t_it ) ) ;
   do_one_inner_iteration( t_it ) ;
}




//----------------------------------------------------------------------------
bool
FV_OneStepIteration:: inner_iterations_stage_failed( void ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: inner_iterations_stage_failed" ) ;
   return( FAILURE ) ;
}




//----------------------------------------------------------------------------
bool
FV_OneStepIteration:: inner_iterations_are_completed(
                                        FV_TimeIterator const* t_it ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: inner_iterations_are_completed" ) ;
   MAC_CHECK_PRE( inner_iterations_are_completed_PRE( t_it ) ) ;
   return( true ) ;
}




//----------------------------------------------------------------------------
void
FV_OneStepIteration:: do_after_inner_iterations_stage(
                                               FV_TimeIterator const* t_it )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: do_after_inner_iterations_stage" ) ;
   MAC_CHECK_PRE( do_after_inner_iterations_stage_PRE( t_it ) ) ;
}




//----------------------------------------------------------------------------
void
FV_OneStepIteration:: do_after_time_adaptation( FV_TimeIterator const* t_it )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: do_after_time_adaptation" ) ;
   MAC_CHECK_PRE( do_after_time_adaptation_PRE( t_it ) ) ;
}




//----------------------------------------------------------------------------
void
FV_OneStepIteration:: do_after_time_stepping( void )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: do_after_time_stepping" ) ;
}




//----------------------------------------------------------------------------
void
FV_OneStepIteration:: reset_standard_times( void )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: reset_standard_times" ) ;

   std::map< std::string, MAC_Timer* >::const_iterator it ;
   for( it = ASS_TIMER.begin() ; it != ASS_TIMER.end() ; ++it )
   {
      MAC_Timer* tt = it->second ;
      tt->reset() ;
   }
   for( it = SOL_TIMER.begin() ; it != SOL_TIMER.end() ; ++it )
   {
      MAC_Timer* tt = it->second ;
      tt->reset() ;
   }
   for( it = TOT_TIMER.begin() ; it != TOT_TIMER.end() ; ++it )
   {
      MAC_Timer* tt = it->second ;
      tt->reset() ;
   }
}




//----------------------------------------------------------------------------
void
FV_OneStepIteration:: print_standard_times( std::ostream& os,
                                            size_t indent_width )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: print_standard_times" ) ;

   std::string space( indent_width, ' ' ) ;

   size_t l_name = 0 ;
   size_t l_ass = 0 ;
   size_t l_sol = 0 ;
   size_t l_tot = 0 ;
   std::map< std::string, MAC_Timer* >::const_iterator it ;
   std::set< std::string > names ;
   for( it = ASS_TIMER.begin() ; it != ASS_TIMER.end() ; ++it )
   {
      std::string const& nn = it->first ;
      names.insert( nn ) ;
      l_name = MAC::max( l_name, (size_t)nn.length() ) ;
      MAC_Timer* tt = it->second ;
      std::ostringstream time ;
      tt->print( time, 0 ) ;
      l_ass = MAC::max( l_ass, (size_t)time.str().length() ) ;
   }
   for( it = SOL_TIMER.begin() ; it != SOL_TIMER.end() ; ++it )
   {
      std::string const& nn = it->first ;
      names.insert( nn ) ;
      l_name = MAC::max( l_name, nn.length() ) ;
      MAC_Timer* tt = it->second ;
      std::ostringstream time ;
      tt->print( time, 0 ) ;
      l_sol = MAC::max( l_sol, (size_t)time.str().length() ) ;
   }
   for( it = TOT_TIMER.begin() ; it != TOT_TIMER.end() ; ++it )
   {
      std::string const& nn = it->first ;
      names.insert( nn ) ;
      l_name = MAC::max( l_name, nn.length() ) ;
      MAC_Timer* tt = it->second ;
      std::ostringstream time ;
      tt->print( time, 0 ) ;
      l_tot = MAC::max( l_tot, time.str().length() ) ;
   }

   l_ass = MAC::max( l_ass, (size_t)10 ) ; // 10 = length of "Assembling"
   l_sol = MAC::max( l_sol, (size_t) 7 ) ; //  7 = length of "Solving"
   l_tot = MAC::max( l_tot, (size_t) 5 ) ; //  5 = length of "Total"
   os << space
      << setw( l_ass+3 ) << "Assembling" << "  "
      << setw( l_sol+2 ) << "Solving" << "  "
      << setw( l_tot+2 ) << "Total" << std::endl ;
   std::set< std::string >::const_iterator itn = names.begin() ;
   for( ; itn != names.end() ; ++itn )
   {
      std::string const& nn = *itn ;

      std::ostringstream time_ass ;
      it = ASS_TIMER.find( nn ) ;
      if( it != ASS_TIMER.end() )
      {
         it->second->print( time_ass, 0 ) ;
      }
      std::ostringstream time_sol ;
      it = SOL_TIMER.find( nn ) ;
      if( it != SOL_TIMER.end() )
      {
         it->second->print( time_sol, 0 ) ;
      }
      std::ostringstream time_tot ;
      it = TOT_TIMER.find( nn ) ;
      if( it != TOT_TIMER.end() )
      {
         it->second->print( time_tot, 0 ) ;
      }

      os << space
         << setw( l_ass+3 ) << time_ass.str() << "  "
         << setw( l_sol+2 ) << time_sol.str() << "  "
         << setw( l_tot+2 ) << time_tot.str() << "  "
         << nn << std::endl ;
   }
}




//----------------------------------------------------------------------------
void
FV_OneStepIteration:: print_additional_times( std::ostream& os,
                                              size_t indent_width ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: print_additional_times" ) ;
}




//----------------------------------------------------------------------------
void
FV_OneStepIteration:: adapt_time_iterator( FV_TimeIteratorAdapter* t_adapter )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: adapt_time_iterator" ) ;
   MAC_CHECK_PRE( adapt_time_iterator_PRE( t_adapter ) ) ;
}




//----------------------------------------------------------------------------
void
FV_OneStepIteration:: notify_inner_iterations_stage_failure( void )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: notify_inner_iterations_stage_failure" ) ;
   MAC_CHECK_PRE( notify_inner_iterations_stage_failure_PRE() ) ;
   FAILURE = true ;
   MAC_CHECK_POST( notify_inner_iterations_stage_failure_POST() ) ;
}




//----------------------------------------------------------------------------
void
FV_OneStepIteration:: reset_after_inner_iterations_stage_failure( void )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration::"
              "reset_after_inner_iterations_stage_failure" ) ;
   MAC_CHECK_PRE( reset_after_inner_iterations_stage_failure_PRE() ) ;
   FAILURE = false ;
   MAC_CHECK_POST( reset_after_inner_iterations_stage_failure_POST() ) ;
}




//----------------------------------------------------------------------------
void
FV_OneStepIteration:: do_additional_savings( FV_TimeIterator const* t_it,
	int const& cycleNumber )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: do_additional_savings" ) ;
   MAC_CHECK_PRE( do_additional_savings_PRE( t_it ) ) ;

}




//----------------------------------------------------------------------------
void
FV_OneStepIteration:: do_additional_save_for_restart( 
	FV_TimeIterator const* t_it, size_t const& restartCycleNumber,
	std::string const& basename )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: do_additional_savings" ) ;
   MAC_CHECK_PRE( do_additional_savings_PRE( t_it ) ) ;

}




//---------------------------------------------------------------------------- 
void
FV_OneStepIteration:: register_storable_objects( MAC_ListIdentity* list )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: register_storable_objects" ) ;

   add_storable_objects( list ) ;
}




//----------------------------------------------------------------------------
void
FV_OneStepIteration:: add_storable_objects( MAC_ListIdentity* list )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: add_storable_objects" ) ;

   MAC_CHECK_PRE( add_storable_objects_PRE( list ) ) ;
}




//----------------------------------------------------------------------------
void
FV_OneStepIteration:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: print" ) ;

   std::string const s( indent_width, ' ' ) ;
   os << s << type_name() << std::endl << std::endl;
}




//----------------------------------------------------------------------------
std::string const&
FV_OneStepIteration:: indent( void )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: indent" ) ;

   return( INDENT ) ;
}




//----------------------------------------------------------------------------
void
FV_OneStepIteration:: increase_indent( void )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: increase_indent" ) ;

   INDENT += "   " ;
}




//----------------------------------------------------------------------------
void
FV_OneStepIteration:: decrease_indent( void )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: decrease_indent" ) ;

   if( INDENT.size() < 3 )
   {
      MAC_Error::object()->raise_plain( "attempt to decrease indentation "
                                        "below the zero limit" ) ;
   }
   INDENT.erase( INDENT.length()-3 ) ;
}




//----------------------------------------------------------------------------
bool
FV_OneStepIteration:: is_a_prototype( void ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: is_a_prototype" ) ;

   return( IS_PROTO ) ;
}




//----------------------------------------------------------------------------
MAC_Communicator const*
FV_OneStepIteration:: communicator( void )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: communicator" ) ;

   static MAC_Communicator const* result = MAC_Exec::communicator() ;

   MAC_CHECK_POST( result == MAC_Exec::communicator() ) ;
   return result ;
}




//----------------------------------------------------------------------------
size_t
FV_OneStepIteration:: verbose_level( void ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: verbose_level" ) ;

   return( VERBOSE_LEVEL ) ;
}




//----------------------------------------------------------------------------
void
FV_OneStepIteration:: start_assembling_timer( bool silent ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: start_assembling_timer" ) ;

   bool do_display = ( !silent && verbose_level() >= 2 ) ;
   start_internal_timer( ASS_TIMER, "assemble... ", do_display ) ;
}




//----------------------------------------------------------------------------
void
FV_OneStepIteration:: stop_assembling_timer( bool silent ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: stop_assembling_timer" ) ;

   bool do_display = ( !silent && verbose_level() >= 2 ) ;
   stop_internal_timer( ASS_TIMER, do_display ) ;
}




//----------------------------------------------------------------------------
void
FV_OneStepIteration:: start_solving_timer( bool silent ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: start_solving_timer" ) ;

   bool do_display = ( !silent && verbose_level() >= 2 ) ;
   start_internal_timer( SOL_TIMER, "solve...", do_display ) ;
}




//----------------------------------------------------------------------------
void
FV_OneStepIteration:: stop_solving_timer( bool silent ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: stop_solving_timer" ) ;

   bool do_display = ( !silent && verbose_level() >= 2 ) ;
   stop_internal_timer( SOL_TIMER, do_display ) ;
}




//----------------------------------------------------------------------------
void
FV_OneStepIteration:: start_total_timer( std::string const& mesg,
                                         bool silent ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: start_total_timer" ) ;

   bool do_display = ( !silent && verbose_level() >= 1 ) ;
   start_internal_timer( TOT_TIMER, mesg+"...", do_display ) ;
}




//----------------------------------------------------------------------------
void
FV_OneStepIteration:: stop_total_timer( bool silent ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: stop_total_timer" ) ;

   bool do_display = ( !silent && verbose_level() >= 1 ) ;
   stop_internal_timer( TOT_TIMER, do_display ) ;
}




//----------------------------------------------------------------------
bool
FV_OneStepIteration:: do_before_time_stepping_PRE(
                                         FV_TimeIterator const* t_it ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( t_it != 0 ) ;
   MAC_ASSERT( !t_it->is_started() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
FV_OneStepIteration:: do_before_inner_iterations_stage_PRE(
                                         FV_TimeIterator const* t_it ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( t_it != 0 ) ;
   MAC_ASSERT( t_it->is_started() ) ;
   MAC_ASSERT( !t_it->is_finished() ) ;
   MAC_ASSERT( !inner_iterations_stage_failed() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
FV_OneStepIteration:: do_before_inner_iterations_stage_POST(
                                         FV_TimeIterator const* t_it ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( !inner_iterations_stage_failed() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
FV_OneStepIteration:: do_inner_iterations_stage_PRE(
                                         FV_TimeIterator const* t_it ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( t_it != 0 ) ;
   MAC_ASSERT( t_it->is_started() ) ;
   MAC_ASSERT( !t_it->is_finished() ) ;
   MAC_ASSERT( !inner_iterations_stage_failed() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
FV_OneStepIteration:: do_one_inner_iteration_PRE(
                                         FV_TimeIterator const* t_it ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( t_it != 0 ) ;
   MAC_ASSERT( t_it->is_started() ) ;
   MAC_ASSERT( !t_it->is_finished() ) ;
   MAC_ASSERT( !inner_iterations_stage_failed() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
FV_OneStepIteration:: inner_iterations_are_completed_PRE(
                                         FV_TimeIterator const* t_it ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( t_it != 0 ) ;
   MAC_ASSERT( t_it->is_started() ) ;
   MAC_ASSERT( !t_it->is_finished() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
FV_OneStepIteration:: do_after_inner_iterations_stage_PRE(
                                         FV_TimeIterator const* t_it ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( t_it != 0 ) ;
   MAC_ASSERT( t_it->is_started() ) ;
   MAC_ASSERT( !t_it->is_finished() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
FV_OneStepIteration:: do_after_time_adaptation_PRE(
                                         FV_TimeIterator const* t_it ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( t_it != 0 ) ;
   MAC_ASSERT( t_it->is_started() ) ;
   MAC_ASSERT( !t_it->just_went_back() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
FV_OneStepIteration:: do_additional_savings_PRE(
                                         FV_TimeIterator const* t_it ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( t_it != 0 ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
FV_OneStepIteration:: notify_inner_iterations_stage_failure_PRE( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( !inner_iterations_stage_failed() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
FV_OneStepIteration:: notify_inner_iterations_stage_failure_POST( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( inner_iterations_stage_failed() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
FV_OneStepIteration::reset_after_inner_iterations_stage_failure_PRE( 
    void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( inner_iterations_stage_failed() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
FV_OneStepIteration::reset_after_inner_iterations_stage_failure_POST(
    void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( !inner_iterations_stage_failed() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
FV_OneStepIteration:: adapt_time_iterator_PRE(
                         FV_TimeIteratorAdapter const* t_adapter ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( t_adapter != 0 ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
FV_OneStepIteration:: add_storable_objects_PRE( MAC_ListIdentity* list ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( list != 0 ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
FV_OneStepIteration:: create_replica_PRE(
		MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( is_a_prototype() ) ;
   MAC_ASSERT( exp != 0 ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
bool
FV_OneStepIteration:: create_replica_POST(
		FV_OneStepIteration const* result,
		MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( result != 0 ) ;
   MAC_ASSERT( result->owner() == a_owner ) ;
   MAC_ASSERT( !result->is_a_prototype() ) ;
   return( true ) ;
}




//----------------------------------------------------------------------
MAC_ObjectRegister*
FV_OneStepIteration:: plugins_map( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: plugins_map" ) ;

   static MAC_ObjectRegister* result =
        MAC_ObjectRegister::create( MAC_Root::object(),
                                    "FV_OneStepIteration descendant" ) ;
   return( result ) ;
}




//----------------------------------------------------------------------------
void
FV_OneStepIteration:: start_internal_timer(
                            std::map< std::string, MAC_Timer* >& timers,
                            std::string const& display,
                            bool do_display ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: start_internal_timers" ) ;

   if( do_display && RANK == 0 )
   {
      increase_indent() ;
      FV::out() << indent() << display << endl ;
   }

   MAC_Timer* tt = 0 ;
   std::string const& nn = type_name() ;
   std::map< std::string, MAC_Timer* >::const_iterator it = timers.find( nn ) ;
   if( it == timers.end() )
   {
      tt = MAC_Timer::create( MAC_Root::object() ) ;
      timers[ nn ] = tt ;
   }
   else
   {
      tt = it->second ;
   }
   tt->start() ;
}




//----------------------------------------------------------------------------
void
FV_OneStepIteration:: stop_internal_timer(
                           std::map< std::string, MAC_Timer* > const& timers,
                           bool do_display ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: stop_internal_timers" ) ;

   if( do_display && RANK == 0 )
   {
      decrease_indent() ;
   }

   std::map< std::string, MAC_Timer* >::const_iterator it =
                                              timers.find( type_name() ) ;
   MAC_ASSERT( it != timers.end() ) ;

   MAC_Timer* tt = it->second ;
   tt->stop() ;
}




//----------------------------------------------------------------------------
void
FV_OneStepIteration:: do_more_post_processing( 
        FV_DomainAndFields* dom,
        MAC_ModuleExplorer const* exp )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_OneStepIteration:: do_more_post_processing" ) ;

}
