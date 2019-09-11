#include <FV_StepByStepProgression.hh>
#include <FV_DomainAndFields.hh>
#include <FV.hh>
#include <MAC.hh>
#include <MAC_Communicator.hh>
#include <MAC_Communicator.hh>
#include <MAC_Exec.hh>
#include <MAC_Error.hh>
#include <MAC_ListIdentity.hh>
#include <MAC_MemoryTracer.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_ApplicationRestorer.hh>
#include <MAC_System.hh>
#include <MAC_Root.hh>
#include <MAC_Timer.hh>
#include <FV_OneStepIteration.hh>
#include <FV_TimeIteratorAdapter.hh>
#include <FV_TimeIterator.hh>
#include <intVector.hh>
#include <doubleVector.hh>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <sys/types.h>
#include <unistd.h>

using std::cout ;
using std::endl ;
using std::ostringstream ;
using std::string ;

size_t FV_StepByStepProgression::graphics_level = 0 ;

FV_StepByStepProgression const* 
FV_StepByStepProgression:: PROTOTYPE = new FV_StepByStepProgression() ;


//----------------------------------------------------------------------
FV_StepByStepProgression:: FV_StepByStepProgression( void )
//----------------------------------------------------------------------
   : MAC_Application( "FV_StepByStepProgression" )
   , TIME_IT( 0 )
   , ONE_IT( 0 )
   , GRAPHICS_TIMES( 0 )
   , GRAPHICS_NEXT_TIME( MAC::bad_double() )
   , SAVER_TIMES( 0 )
   , SAVER_NEXT_TIME( MAC::bad_double() )
   , overall( 0 )
   , POST_TIMER( 0 )
   , SAVE_TIMER( 0 )
   , SAVEFG( true )
   , VERBOSE( true )
   , POST_CYCLE_NUMBER( 1 )
   , RANK( 0 )
{
   MAC_LABEL( "FV_StepByStepProgression:: FV_StepByStepProgression" ) ;
}




//----------------------------------------------------------------------
FV_StepByStepProgression*
FV_StepByStepProgression:: create( MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
	double const& initial_time )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_StepByStepProgression:: create" ) ;
   MAC_CHECK_PRE( exp != 0 ) ;

   FV_StepByStepProgression* result = 
	new FV_StepByStepProgression( a_owner, exp, initial_time ) ;
   result->register_storable_objects() ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
FV_StepByStepProgression*
FV_StepByStepProgression:: create_replica( 
	MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
	double const& initial_time ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_StepByStepProgression:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   FV_StepByStepProgression* result = 
	new FV_StepByStepProgression( a_owner, exp, initial_time ) ;

   MAC_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
FV_StepByStepProgression:: FV_StepByStepProgression( 
	MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
	double const& initial_time )
//----------------------------------------------------------------------
   : MAC_Application( a_owner, exp )
   , TIME_IT( 0 )
   , ONE_IT( 0 )
   , GRAPHICS_TIMES( 0 )
   , GRAPHICS_NEXT_TIME( MAC::max_double() )
   , SAVER_TIMES( 0 )
   , SAVER_NEXT_TIME( MAC::max_double() )
   , overall( MAC_Timer::create( this ) )
   , POST_TIMER( MAC_Timer::create( this ) )
   , SAVE_TIMER( MAC_Timer::create( this ) )
   , SAVEFG( true )
   , LAST_MEM_IT( MAC::bad_index() )
   , VERBOSE( true )
   , POST_CYCLE_NUMBER( 1 )
   , RANK( MAC_Exec::communicator()->rank() )   
{
   MAC_LABEL( "FV_StepByStepProgression:: FV_StepByStepProgression" ) ;
   
   if ( exp->has_entry( "verbose" ) ) VERBOSE = exp->bool_data( "verbose" ) ;
   
   if ( exp->has_module( "memory_trace" ) )
   {
      MAC_MemoryTracer::object()->enable_memory_trace() ;
      MAC_ModuleExplorer* ee =
                       exp->create_subexplorer( 0, "memory_trace" ) ;
      std::string const& type = ee->string_data( "type" ) ;
      ee->test_data_in( "type", "trace_all_iterations,trace_first_iterations" );
      if( type == "trace_first_iterations" )
      {
         LAST_MEM_IT = ee->int_data( "last_iteration_checked" ) ;
      }
      ee->destroy() ; ee = 0 ;
   }
   
   if ( VERBOSE ) display_memory_usage_proc_by_proc( FV::out(), 4 ) ;

   if ( MAC_MemoryTracer::object()->memory_trace_enabled() )
   {
      MAC_MemoryTracer::object()->trace(
                   MAC_MemoryTracer::object()->message( "Start program" ) ) ;
   }

   // FV_DomainAndFields:
   if ( exp->has_module( "FV_DomainAndFields" ) )
   {
      MAC_ModuleExplorer* ee =
                       exp->create_subexplorer( 0, "FV_DomainAndFields" ) ;
      MAC_MemoryTracer::object()->start_event( "Building FV_DomainAndFields" ) ;
      DOM = FV_DomainAndFields::create( this, ee, MAC_Exec::communicator() ) ;
      MAC_MemoryTracer::object()->stop_event() ;
      ee->destroy() ; ee = 0 ;
   }
      
   // FV_OneStepIteration:
   if ( exp->has_module( "FV_OneStepIteration" ) )
   {
      if( VERBOSE && RANK == 0 ) FV::out() << endl 
      	<< "*** building FV_OneStepIteration" << endl << endl;
      MAC_ModuleExplorer* ee =
                         exp->create_subexplorer( 0, "FV_OneStepIteration" ) ;
      MAC_MemoryTracer::object()->start_event( "Building FV_OneStepIteration" );
      ONE_IT = FV_OneStepIteration::make( this, DOM, ee ) ;
      MAC_MemoryTracer::object()->stop_event() ;
      ee->destroy() ; ee = 0 ;
      FV_OneStepIteration::reset_standard_times() ;
   }

   // FV_TimeIterator:
   if ( exp->has_module( "FV_TimeIterator" ) )
   {
      MAC_ModuleExplorer* ee = exp->create_subexplorer( 0, "FV_TimeIterator" ) ;
      TIME_IT = FV_TimeIterator::create( this, ee, initial_time ) ;
      ee->destroy() ; ee = 0 ;
   }

   if ( exp->has_entry( "graphics_output_times" ) )
   {
      GRAPHICS_TIMES = exp->doubleVector_data( "graphics_output_times" ) ;
   }
   else if ( exp->has_entry( "number_graphics_output_times" ) )
   {
      size_t ngot = exp->int_data( "number_graphics_output_times" );
      GRAPHICS_TIMES.re_initialize( ngot, TIME_IT->initial_time(),
      	TIME_IT->final_time() ) ;
   }   

   if ( exp->has_entry( "state_saving_times" ) )
   {
      SAVER_TIMES = exp->doubleVector_data( "state_saving_times" ) ;
   }
   else if ( exp->has_entry( "number_state_saving_times" ) )
   {
      size_t nsst = exp->int_data( "number_state_saving_times" );
      SAVER_TIMES.re_initialize( nsst, TIME_IT->initial_time(),
      	TIME_IT->final_time() ) ;
   }
   
   if ( exp->has_entry( "save_grid_and_fields_for_postprocessing" ) )
   {
      SAVEFG = exp->bool_data( "save_grid_and_fields_for_postprocessing" ) ;
   }
   
   if ( exp->has_entry( "POST_INITIAL_CYCLE_NUMBER" ) )
   {
     POST_CYCLE_NUMBER = exp->int_data( "POST_INITIAL_CYCLE_NUMBER" ) ;
   }

   if ( VERBOSE && RANK == 0 )
      ONE_IT->print( FV::out(), 3 ) ;
   MAC_Exec::communicator()->barrier();   
   
   overall->start() ;
}




//----------------------------------------------------------------------
FV_StepByStepProgression:: ~FV_StepByStepProgression( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_StepByStepProgression:: ~FV_StepByStepProgression" ) ;
}




//----------------------------------------------------------------------
FV_TimeIterator const* 
FV_StepByStepProgression:: time_iterator( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_StepByStepProgression:: time_iterator" ) ;

   FV_TimeIterator const* result = TIME_IT ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->is_under_ownership_of( this ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
unsigned long long int
FV_StepByStepProgression:: display_memory_usage(
	std::ostream& os, size_t indent_width )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_StepByStepProgression:: display_memory_usage" ) ;
   MAC_CHECK_PRE( os ) ;

   std::string const s( indent_width, ' ' ) ;
   os << s << "Memory usage: " ;
   unsigned long long int result = MAC::used_memory() ;
   MAC::display_memory( os, result ) ;
   os << endl
      << s << "Number of objects: "
      << MAC_Object::GetNumberOf_MAC_objects()
      << endl ;

   return( result ) ;            
}




//----------------------------------------------------------------------
void
FV_StepByStepProgression:: display_memory_usage_proc_by_proc(
	std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( 
     "FV_StepByStepProgression:: display_memory_usage_proc_by_proc" ) ;
   MAC_CHECK_PRE( os ) ;
   
   MAC_Communicator const* macCOMM = MAC_Exec::communicator();
   size_t nb_ranks = macCOMM->nb_ranks() ;
   unsigned long long int local_mem = 0 ;
   std::string space( indent_width, ' ' ) ;
   
   macCOMM->barrier(); 
   if ( RANK == 0 )
     FV::out() << endl << space << "++++++ MEMORY USAGE ++++++" 
     	<< endl << endl ;     
   macCOMM->barrier(); 
         
   std::streamsize p = os.precision() ;
   os << std::setprecision( 3 );
   os.setf( std::ios::fixed, std::ios::floatfield );  

   for (size_t i=0;i<nb_ranks;++i)
   {
     if ( i == RANK )
     {       
       os << space << "Rank " << RANK << endl ;
       local_mem = display_memory_usage( os, indent_width + 3 ) ;       
     }
     macCOMM->barrier();
   }
   
   unsigned long long int allprocs_mem = macCOMM->sum( local_mem );
   if ( RANK == 0 ) 
   {
     os << endl << space << "Total memory on all procs" << endl 
     	<< space << "   Memory usage: "; 
     MAC::display_memory( os, allprocs_mem );     
     os << endl;
   } 
   
   os << std::setprecision( p );
   os.unsetf( std::ios::floatfield );       
}




//----------------------------------------------------------------------
void
FV_StepByStepProgression:: add_storable_objects( MAC_ListIdentity* list )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_StepByStepProgression:: add_storable_objects" ) ;

   if( DOM != 0 )
   {
      list->extend( DOM ) ;
   }
   ONE_IT->register_storable_objects( list ) ;
}




//----------------------------------------------------------------------
void
FV_StepByStepProgression:: run( std::string const& inputRestartFileName )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_StepByStepProgression:: run" ) ; 

   check_times( "graphics_output_times", GRAPHICS_TIMES ) ;
   check_times( "state_saving_times", SAVER_TIMES ) ;
   
   swap_restart_file_names( inputRestartFileName );
   
   GRAPHICS_NEXT_TIME = TIME_IT->next_time_in_table( GRAPHICS_TIMES ) ;
   SAVER_NEXT_TIME = TIME_IT->next_time_in_table( SAVER_TIMES ) ;
   
   if ( RANK == 0 )
     FV::out() << endl << "++++++ TIME STEPPING " 
              << " *** INITIAL TIME = " << TIME_IT->initial_time()
              << " *** FINAL TIME = " << TIME_IT->final_time()
              << " ++++++" << endl << endl ;

   MAC_MemoryTracer::object()->start_event( "Time step initialization" ) ;
   ONE_IT->do_before_time_stepping( TIME_IT, inputRestartFileName ) ;
   MAC_MemoryTracer::object()->stop_event() ;

   save_for_post_processing( true ) ;

   for( TIME_IT->start() ; !TIME_IT->is_finished() ; TIME_IT->go_next_time() )
   {
      if( VERBOSE && RANK == 0 )
      {
         std::streamsize p = FV::out().precision() ;
         FV::out() << std::setprecision( 12 ) << endl ;;
         FV::out() << "++++++ ITERATION = " << TIME_IT->iteration_number()
                    << " *** TIME = " << TIME_IT->time()
                    << " *** TIME STEP = " << TIME_IT->time_step()
                    << " ++++++" << endl << endl ;
         FV::out() << std::setprecision( p ) ;
      }

      std::ostringstream mem_label ;
      mem_label << "Time iteration " << TIME_IT->iteration_number() ;
      MAC_MemoryTracer::object()->start_event( mem_label.str() ) ;

      ONE_IT->do_before_inner_iterations_stage( TIME_IT ) ;
      ONE_IT->do_inner_iterations_stage( TIME_IT ) ;
      if( ONE_IT->inner_iterations_stage_failed() )
      {
         ONE_IT->reset_after_inner_iterations_stage_failure() ;
      }
      else
      {
         ONE_IT->do_after_inner_iterations_stage( TIME_IT ) ;
         if( !TIME_IT->just_went_back() )
         {
            ONE_IT->do_after_time_adaptation( TIME_IT ) ;
            save_for_post_processing( false ) ;
            save_for_restart( false ) ;	    
         }
      }

      MAC_MemoryTracer::object()->stop_event() ;
      if( MAC_MemoryTracer::object()->memory_trace_enabled() &&
          TIME_IT->iteration_number() == LAST_MEM_IT )
      {
         MAC_MemoryTracer::object()->disable_memory_trace() ;
         FV::out() << std::endl ;
      }
   }

   if ( VERBOSE && RANK == 0 )
   {
      FV::out() << endl ;
      if( TIME_IT->finished_reason() == FV_TimeIterator::UserCheckpoint )
         FV::out() << "++++++ TIME STEPPING ENDED BY USER " ;
      else
         FV::out() << "++++++ TIME STEPPING COMPLETED " ;
      FV::out() << " ++++++" << endl << endl ;
   }

   ONE_IT->do_after_time_stepping() ;
   
   // Remarks independent of the reason the code stops (end of time marching or
   // stopped by the user)
   // 1) post processing files are not written by default at last time, this is
   // the user's responsability i.e.how "graphics_output_times" was defined in 
   // the input file
   // 2) restart files are not written by default at last time unless the
   // keyword "force_write_at_last_time" is set to true in the MAC_ObjectWriter
   // module of the input file 
   
   save_for_restart( force_write_restart_files_at_last_time() ) ;
   
   overall->stop() ;
   
   if ( !VERBOSE && RANK == 0 ) 
   {
      FV::out() << endl << "++++++ " << TIME_IT->iteration_number() 
                 << " ITERATIONS ++++++" << endl ;
   }
   
   display_timers_proc_by_proc( 0 ) ;
   if( RANK == 0 ) FV::out() << endl ;
}





//----------------------------------------------------------------------
void
FV_StepByStepProgression:: do_more_post_processing( 
	MAC_ModuleExplorer const* exp  )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_StepByStepProgression:: do_more_post_processing" ) ;
   
   ONE_IT->do_more_post_processing( DOM, exp );
   
   if ( VERBOSE ) display_memory_usage_proc_by_proc( FV::out(), 4 ) ; 
}




//----------------------------------------------------------------------
void
FV_StepByStepProgression:: save_for_post_processing( bool const& force )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_StepByStepProgression:: save_for_post_processing" ) ;

   POST_TIMER->start() ;
   
   static size_t last_iter_saved = MAC::bad_index() ;
   size_t const iter =
      TIME_IT->is_started() ? TIME_IT->iteration_number() : 0 ;
   bool const save_iteration =
   	FV_TimeIterator::greater_or_equal( TIME_IT->time(),
                                                  GRAPHICS_NEXT_TIME ) ;

   if ( ( save_iteration || force ) && iter!=last_iter_saved )
   {
      do_save_for_post( DOM->post_processing_writer() ) ;
      if ( RANK == 0 ) display_timers( 3 ) ;
      if ( VERBOSE )
      {
        if ( RANK == 0 ) FV::out() << std::endl ;
        display_memory_usage_proc_by_proc( FV::out(), 3 ) ;
      }
      last_iter_saved = iter ;
      GRAPHICS_NEXT_TIME = TIME_IT->next_time_in_table( GRAPHICS_TIMES ) ;
   }
   
   POST_TIMER->stop() ;
}




//----------------------------------------------------------------------
void
FV_StepByStepProgression:: display_timers( size_t indent_width ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_StepByStepProgression:: display_timers" ) ;
   
   std::string space( indent_width, ' ' ) ;
   FV::out() << endl << space << "++++++ TIMERS ON RANK " << RANK 
   	<< " ++++++" << endl << endl ;
   FV::out() << space << "total time : " ;
   overall->print( FV::out(), 0 ) ; FV::out() << endl ;
   FV::out() << space << "   including \"save_for_post_processing\" : " ;
   POST_TIMER->print( FV::out(), 0 ) ; FV::out() << endl ;
   FV::out() << space << "   including \"save_for_restart\" : " ;
   SAVE_TIMER->print( FV::out(), 0 ) ; FV::out() << endl ;
   FV_OneStepIteration::print_standard_times( FV::out(), indent_width ) ;
   ONE_IT->print_additional_times( FV::out(), indent_width ) ;
}




//----------------------------------------------------------------------
void
FV_StepByStepProgression:: display_timers_proc_by_proc(
		size_t indent_width ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( 
     "FV_StepByStepProgression:: display_timers_proc_by_proc" ) ;
   
   size_t nb_ranks = MAC_Exec::communicator()->nb_ranks() ;
        
   for (size_t i=0;i<nb_ranks;++i)
   {
     if ( i == RANK ) display_timers( indent_width + 3 ) ;       
     MAC_Exec::communicator()->barrier();
   }        
}




//----------------------------------------------------------------------
void
FV_StepByStepProgression:: do_save_for_post(
		FV_PostProcessingWriter* mpp_writer )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_StepByStepProgression:: do_save_for_post" ) ;

   static size_t post_counter = 0 ;

   // In case of restart (MAC_ApplicationRestorer), automatically get the last 
   // cycle number in the results directory  
   // Otherwise, clear results files in the results directory 
   if ( post_counter == 0 )
   { 
     if ( is_follow() )
     {
       POST_CYCLE_NUMBER = mpp_writer->getPreviousCycleNumber() ;
       mpp_writer->readTimeFile( TIME_IT, POST_CYCLE_NUMBER ) ;   
     }
     else
       mpp_writer->clearResultFiles() ;
   }
   // Synchronize all procs here, otherwise non-master processes may create 
   // their result files before the master has erased old results files 
   // and hence new created result files by non-master processess are deleted
   MAC_Exec::communicator()->barrier();
   
   if ( RANK == 0 )
   {
     std::streamsize p = FV::out().precision() ;
     FV::out() << std::setprecision(12) ;
     FV::out() << endl
              << "   +++ SAVE FOR POSTPROCESSING"
              << " *** CYCLE = " << POST_CYCLE_NUMBER
              << " *** TIME = "  << TIME_IT->time()
              << " ++++++" << endl << endl ;
     FV::out() << std::setprecision(p) ;
   }   
    
   // User defined savings
   ONE_IT->do_additional_savings( TIME_IT, POST_CYCLE_NUMBER ) ;
   
   // PostProcessingWriter performs the actual writing of the data
   if ( SAVEFG ) mpp_writer->write_cycle( TIME_IT, POST_CYCLE_NUMBER );
   
   ++POST_CYCLE_NUMBER;
   ++post_counter;
}




//----------------------------------------------------------------------
void
FV_StepByStepProgression:: check_times( std::string const& list_name,
                                        doubleVector const& dates ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_StepByStepProgression:: check_times" ) ;

   if( !TIME_IT->table_of_times_is_valid( dates ) )
   {
      TIME_IT->raise_invalid_table_of_times( "FV_StepByStepProgression",
                                             list_name,
                                             dates ) ;
   }
}




//----------------------------------------------------------------------
void
FV_StepByStepProgression:: save_for_restart( bool const& force )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_StepByStepProgression:: save_for_restart" ) ;

   SAVE_TIMER->start() ;
   
   static size_t last_iter_saved = MAC::bad_index() ;
   size_t const iter =
      TIME_IT->is_started() ? TIME_IT->iteration_number() : 0 ;
   bool save_iteration =
      FV_TimeIterator::greater_or_equal( TIME_IT->time(), SAVER_NEXT_TIME ) ;
   
   // The SAVER, i.e., MAC_ObjectWriter, is always implemented by
   // FV_StepByStepProgression and no other MAC_Application
   if ( ( save_iteration || force ) && iter!=last_iter_saved )
   {
      write_storable_objects( TIME_IT->time() ) ;
      ONE_IT->do_additional_save_for_restart( TIME_IT, restart_cycle_number(),
      	output_restart_file_name() );
      SAVER_NEXT_TIME = TIME_IT->next_time_in_table( SAVER_TIMES ) ;
      last_iter_saved = iter ;
   }
   
   SAVE_TIMER->stop() ;
}




//----------------------------------------------------------------------
bool
FV_StepByStepProgression:: invariant( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_StepByStepProgression:: invariant" ) ;
   MAC_ASSERT( MAC_Application::invariant() ) ;

   return( true ) ;
}
