#include <FV_SplitSystem.hh>
#include <MAC_Error.hh>
#include <MAC_List.hh>
#include <MAC_ListIterator.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_assertions.hh>
#include <iostream>
using std::endl ;
using std::string ;

FV_SplitSystem const* FV_SplitSystem::PROTOTYPE = new FV_SplitSystem() ;


//-------------------------------------------------------------------------
FV_SplitSystem:: FV_SplitSystem( void )
//-------------------------------------------------------------------------
   : FV_OneStepIteration( "FV_SplitSystem" )
{
}




//-------------------------------------------------------------------------
FV_SplitSystem*
FV_SplitSystem:: create_replica( MAC_Object* a_owner,
	FV_DomainAndFields const* dom,
	MAC_ModuleExplorer* exp ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_SplitSystem:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, dom, exp ) ) ;

   FV_SplitSystem* result = new FV_SplitSystem( a_owner, dom, exp ) ;

   MAC_CHECK( create_replica_POST( result, a_owner, dom, exp ) ) ;
   return( result ) ;
}




//-------------------------------------------------------------------------
FV_SplitSystem:: FV_SplitSystem( MAC_Object* a_owner,
	FV_DomainAndFields const* dom,
	MAC_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
   : FV_OneStepIteration( a_owner, dom, exp )
   , CMPS( MAC_List::create( this ) )
   , IT( 0 )
{
   MAC_ModuleExplorer* e = exp->create_subexplorer( 0, 
                                             "list_of_FV_OneStepIteration" ) ;
   e->start_module_iterator() ;
   for( ; e->is_valid_module() ; e->go_next_module() )
   {
      MAC_ModuleExplorer* ee = e->create_subexplorer( 0 ) ; 
      FV_OneStepIteration* cmp = 
                           FV_OneStepIteration::make( CMPS, dom, ee ) ;
      CMPS->append( cmp ) ;
      ee->destroy() ; ee = 0 ;
   }
   e->destroy() ; e = 0 ;

   IT = MAC_ListIterator::create( this, CMPS ) ;
}



//-------------------------------------------------------------------------
FV_SplitSystem:: ~FV_SplitSystem( void ) 
//-------------------------------------------------------------------------
{
}




//-------------------------------------------------------------------------
void
FV_SplitSystem:: do_before_time_stepping( FV_TimeIterator const* t_it, 
      	std::string const& basename )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_SplitSystem:: do_before_time_stepping" ) ;
   MAC_CHECK_PRE( do_before_time_stepping_PRE( t_it ) ) ;

   start_total_timer( "FV_SplitSystem:: do_before_time_stepping" ) ;
   
   for( IT->start() ; IT->is_valid() ; IT->go_next() )
   {
      FV_OneStepIteration* pb = 
                           static_cast<FV_OneStepIteration*>( IT->item() ) ;
      pb->do_before_time_stepping( t_it, basename ) ;
   }

   stop_total_timer() ;
}




//-------------------------------------------------------------------------
void
FV_SplitSystem:: do_before_inner_iterations_stage( 
                                             FV_TimeIterator const* t_it )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_SplitSystem:: do_before_inner_iterations_stage" ) ;
   MAC_CHECK_PRE( do_before_inner_iterations_stage_PRE( t_it ) ) ;

   start_total_timer( "FV_SplitSystem:: do_before_inner_iterations_stage" ) ;

   for( IT->start() ; IT->is_valid() ; IT->go_next() )
   {
      FV_OneStepIteration* pb = 
                           static_cast<FV_OneStepIteration*>( IT->item() ) ;
      pb->do_before_inner_iterations_stage( t_it ) ;
   }

   stop_total_timer() ;

   MAC_CHECK_POST( do_before_inner_iterations_stage_POST( t_it ) ) ;
}




//-------------------------------------------------------------------------
void
FV_SplitSystem:: do_inner_iterations_stage( FV_TimeIterator const* t_it )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_SplitSystem:: do_inner_iterations_stage" ) ;
   MAC_CHECK_PRE( do_inner_iterations_stage_PRE( t_it ) ) ;

   start_total_timer( "FV_SplitSystem:: do_inner_iterations_stage" ) ;

   IT->start() ;
   for( ; !inner_iterations_stage_failed() && IT->is_valid() ; IT->go_next() )
   {
      FV_OneStepIteration* pb = 
                           static_cast<FV_OneStepIteration*>( IT->item() ) ;
      pb->do_inner_iterations_stage( t_it ) ;
      if( pb->inner_iterations_stage_failed() )
      {
         notify_inner_iterations_stage_failure() ;
      }
   }

   stop_total_timer() ;
}




//-------------------------------------------------------------------------
void
FV_SplitSystem:: do_one_inner_iteration( FV_TimeIterator const* t_it )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_SplitSystem:: do_one_inner_iteration" ) ;
   MAC_CHECK_PRE( do_one_inner_iteration_PRE( t_it ) ) ;

   MAC_Error::object()->raise_plain( "should not be called" ) ;
}




//-------------------------------------------------------------------------
bool
FV_SplitSystem:: inner_iterations_are_completed( 
                                         FV_TimeIterator const* t_it ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_SplitSystem:: inner_iterations_are_completed" ) ;
   MAC_CHECK_PRE( inner_iterations_are_completed_PRE( t_it ) ) ;

   return( true ) ;
}




//-------------------------------------------------------------------------
void
FV_SplitSystem:: do_after_inner_iterations_stage( FV_TimeIterator const* t_it )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_SplitSystem:: do_after_inner_iterations_stage" ) ;
   MAC_CHECK_PRE( do_after_inner_iterations_stage_PRE( t_it ) ) ;

   start_total_timer( "FV_SplitSystem:: do_after_inner_iterations_stage" ) ;
   
   for( IT->start() ; IT->is_valid() ; IT->go_next() )
   {
      FV_OneStepIteration* pb = 
                           static_cast<FV_OneStepIteration*>( IT->item() ) ;
      pb->do_after_inner_iterations_stage( t_it ) ;
   }
   
   stop_total_timer() ;
}




//-------------------------------------------------------------------------
void
FV_SplitSystem:: do_after_time_adaptation( FV_TimeIterator const* t_it )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_SplitSystem:: do_after_time_adaptation" ) ;
   MAC_CHECK_PRE( do_after_time_adaptation_PRE( t_it ) ) ;

   start_total_timer( "FV_SplitSystem:: do_after_time_adaptation" ) ;
   for( IT->start() ; IT->is_valid() ; IT->go_next() )
   {
      FV_OneStepIteration* pb =  
                           static_cast<FV_OneStepIteration*>( IT->item() ) ;
      pb->do_after_time_adaptation( t_it ) ;
   }
   
   stop_total_timer() ;
}




//-------------------------------------------------------------------------
void
FV_SplitSystem:: do_additional_savings( FV_TimeIterator const* t_it,
      	int const& cycleNumber  )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_SplitSystem:: do_additional_savings" ) ;
   MAC_CHECK_PRE( do_additional_savings_PRE( t_it ) ) ;

   start_total_timer( "FV_SplitSystem:: do_additional_savings" ) ;

   for( IT->start() ; IT->is_valid() ; IT->go_next() )
   {
      FV_OneStepIteration* pb = 
                           static_cast<FV_OneStepIteration*>( IT->item() ) ;
      pb->do_additional_savings( t_it, cycleNumber ) ;
   }

   stop_total_timer() ;
}




//-------------------------------------------------------------------------
void
FV_SplitSystem:: do_additional_save_for_restart( FV_TimeIterator const* t_it,
      	size_t const& restartCycleNumber, std::string const& basename )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_SplitSystem:: do_additional_save_for_restart" ) ;
   MAC_CHECK_PRE( do_additional_savings_PRE( t_it ) ) ;

   start_total_timer( "FV_SplitSystem:: do_additional_save_for_restart" ) ;

   for( IT->start() ; IT->is_valid() ; IT->go_next() )
   {
      FV_OneStepIteration* pb = 
                           static_cast<FV_OneStepIteration*>( IT->item() ) ;
      pb->do_additional_save_for_restart( t_it, restartCycleNumber, basename ) ;
   }

   stop_total_timer() ;
}




//-------------------------------------------------------------------------
void
FV_SplitSystem:: do_after_time_stepping( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_SplitSystem:: do_after_time_stepping" ) ;

   start_total_timer( "FV_SplitSystem:: do_after_time_stepping" ) ;

   for( IT->start() ; IT->is_valid() ; IT->go_next() )
   {
      FV_OneStepIteration* pb = 
                           static_cast<FV_OneStepIteration*>( IT->item() ) ;
      pb->do_after_time_stepping() ;
   }

   stop_total_timer() ;
}




//----------------------------------------------------------------------------
void
FV_SplitSystem:: print_additional_times( std::ostream& os,
                                         size_t indent_width ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_SplitSystem:: print_additional_times" ) ;

   for( IT->start() ; IT->is_valid() ; IT->go_next() )
   {
      FV_OneStepIteration* pb = 
                           static_cast<FV_OneStepIteration*>( IT->item() ) ;
      pb->print_additional_times( os, indent_width ) ;
   }
}




//----------------------------------------------------------------------------
void
FV_SplitSystem:: notify_inner_iterations_stage_failure( void )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_SplitSystem:: notify_inner_iterations_stage_failure" ) ;
   MAC_CHECK_PRE( notify_inner_iterations_stage_failure_PRE() ) ;

   FV_OneStepIteration::notify_inner_iterations_stage_failure() ;
   MAC_ListIterator* it = MAC_ListIterator::create( 0, CMPS ) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      FV_OneStepIteration* pb = 
                           static_cast<FV_OneStepIteration*>( it->item() ) ;
      if( !pb->inner_iterations_stage_failed() )
      {
         pb->notify_inner_iterations_stage_failure() ;
      }
   }
   it->destroy() ; it = 0 ;
   
   MAC_CHECK_POST( notify_inner_iterations_stage_failure_POST() ) ;
}




//----------------------------------------------------------------------------
void
FV_SplitSystem:: reset_after_inner_iterations_stage_failure( void )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_SplitSystem:: reset_after_inner_iterations_stage_failure" ) ;
   MAC_CHECK( reset_after_inner_iterations_stage_failure_PRE() ) ;

   FV_OneStepIteration::reset_after_inner_iterations_stage_failure() ;
   for( IT->start() ; IT->is_valid() ; IT->go_next() )
   {
      FV_OneStepIteration* pb = 
                           static_cast<FV_OneStepIteration*>( IT->item() ) ;
      pb->reset_after_inner_iterations_stage_failure() ;
   }
   
   MAC_CHECK_POST( reset_after_inner_iterations_stage_failure_POST() ) ;
}




//----------------------------------------------------------------------------
void
FV_SplitSystem:: adapt_time_iterator( FV_TimeIteratorAdapter* t_adapter )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "FV_SplitSystem:: adapt_time_iterator" ) ;
   MAC_CHECK_PRE( adapt_time_iterator_PRE( t_adapter ) ) ;

   for( IT->start() ; IT->is_valid() ; IT->go_next() )
   {
      FV_OneStepIteration* pb = 
                           static_cast<FV_OneStepIteration*>( IT->item() ) ;
      pb->adapt_time_iterator( t_adapter ) ;
   }
}




//-------------------------------------------------------------------------
void
FV_SplitSystem:: add_storable_objects( MAC_ListIdentity* list )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_SplitSystem:: add_storable_objects" ) ;

   for( IT->start() ; IT->is_valid() ; IT->go_next() )
   {
      FV_OneStepIteration* pb = 
                           static_cast<FV_OneStepIteration*>( IT->item() ) ;
      pb->add_storable_objects( list ) ;
   }
}




//-------------------------------------------------------------------------
void
FV_SplitSystem:: print( std::ostream& os, size_t indent_width ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_SplitSystem:: print" ) ;

   os << "*** Applications hierarchy" << std::endl << std::endl;

   FV_OneStepIteration:: print( os, indent_width ) ;
   
   for( IT->start() ; IT->is_valid() ; IT->go_next() )
   {
      FV_OneStepIteration* pb = 
                           static_cast<FV_OneStepIteration*>( IT->item() ) ;
      pb->print( os, indent_width+3 ) ;
   }
}




//-------------------------------------------------------------------------
void
FV_SplitSystem:: do_more_post_processing( 
        FV_DomainAndFields* dom,
        MAC_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "FV_SplitSystem:: do_more_post_processing" ) ;

   for( IT->start() ; IT->is_valid() ; IT->go_next() )
   {
      FV_OneStepIteration* pb = 
                           static_cast<FV_OneStepIteration*>( IT->item() ) ;
      pb->do_more_post_processing( dom, exp ) ;
   }
}
