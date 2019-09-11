#include <FV_TimeIterator.hh>

#include <MAC.hh>
#include <MAC_Bool.hh>
#include <MAC_Communicator.hh>
#include <MAC_ContextSimple.hh>
#include <MAC_DataWithContext.hh>
#include <MAC_Double.hh>
#include <MAC_DoubleVector.hh>
#include <MAC_Error.hh>
#include <MAC_Exec.hh>
#include <MAC_Int.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_ObjectReader.hh>
#include <MAC_ObjectWriter.hh>
#include <MAC_System.hh>
#include <MAC_Variable.hh>
#include <MAC_assertions.hh>

#include <fstream>
#include <iostream>
#include <sstream>

using std::endl ;

//----------------------------------------------------------------------
FV_TimeIterator*
FV_TimeIterator:: create( MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
	double const& initial_time )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIterator:: create" ) ;
   MAC_CHECK_PRE( exp != 0 ) ;
   
   FV_TimeIterator* result =  new FV_TimeIterator( a_owner,  exp, 
   	initial_time ) ;
   
   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ; 
   MAC_CHECK_POST( 
      result->initial_time() == exp->double_data( "time_start" ) ) ;
   MAC_CHECK_POST( 
      result->final_time() == exp->double_data( "time_end" ) ) ;
   MAC_CHECK_POST(
      IMPLIES( !result->time_step_is_fixed(),
               result->time_step() == exp->double_data( "time_step" ) ) ) ;
   MAC_CHECK_POST(
      result->initial_iteration_number() == 1 ) ;
   MAC_CHECK_POST(
      IMPLIES( exp->has_entry( "storage_depth" ),
               result->storage_depth() == (size_t) exp->int_data( 
	       "storage_depth" ) ) ) ;
   MAC_CHECK_POST(
      IMPLIES( !exp->has_entry( "storage_depth" ),
	           result->storage_depth() == 5 ) ) ;
   MAC_CHECK_POST( 
      result->time() == result->initial_time() ) ;
   MAC_CHECK_POST( !result->is_started() ) ;
   MAC_CHECK_POST( result->time() == result->initial_time() ) ;
   
   return( result ) ;
}




//----------------------------------------------------------------------
FV_TimeIterator:: FV_TimeIterator( MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
	double const& initial_time )
//----------------------------------------------------------------------
   : MAC_Object( a_owner )
   , STARTED( false )
   , FINISHED_REASON( NotFinished )
   , ITER_INIT( 1 )
   , ITER( MAC::bad_index() )
   , T_INIT( initial_time )
   , T_START( initial_time )
   , T_END( 0 )
   , T( initial_time )
   , VARYING_DT( 0 )
   , TIME( 0 )
   , RESTART_IT( false )
   , BACK_AND_FORTH_IT( false )
   , DT( doubleVector( 5 ) )
   , INI_DT( MAC::bad_double() )
   , NEXT_DT( MAC::bad_double() )
{
   MAC_LABEL( "FV_TimeIterator:: FV_TimeIterator" ) ;

   if ( exp->has_entry( "time_start" ) && initial_time == 0. )
   {
     T_INIT = exp->double_data( "time_start" ) ;
     T_START = exp->double_data( "time_start" ) ;
     T = exp->double_data( "time_start" ) ;     
   }
   
   if ( exp->has_entry( "time_interval" ) )
     T_END = T_INIT + exp->double_data( "time_interval" ) ;
   else if ( exp->has_entry( "time_end" ) )
     T_END = exp->double_data( "time_end" ) ;
   else
     MAC_Error::object()->raise_missing_keyword( exp,
     	"either time_end or time_interval" );

   if ( T_END < T_START )
   {
      MAC_Error::object()->raise_bad_data_value(
         exp, "time_end",
         " a value greater than \"time_start\" is expected" ) ;
   }
   
   if ( exp->has_entry( "storage_depth" ) )
   {
      int s = exp->int_data( "storage_depth" ) ;
      if( s<=0 )
      {
         MAC_Error::object()->raise_bad_data_value(
            exp, "storage_depth",
            "   a positive value is expected" ) ;
      }
      DT.re_initialize( (size_t) s ) ;
   }
   
   MAC_Data* cste_dt = exp->abstract_data( 0, "time_step" ) ;
   if ( cste_dt->value_can_be_evaluated(0) )
   {
      double dt_ini = exp->double_data( "time_step" ) ;
      DT.set( dt_ini ) ;
   }
   else
   {
      MAC_ContextSimple* ctx = MAC_ContextSimple::create( this ) ;
      TIME = MAC_Double::create( ctx, T_INIT ) ;
      ctx->extend( MAC_Variable::object( "DS_T" ), TIME ) ;
      VARYING_DT = exp->abstract_data( this, "time_step", ctx ) ;
      if( !VARYING_DT->value_can_be_evaluated(0) )
      {
         MAC_Error::object()->raise_not_evaluable(
                    exp, "time_step", VARYING_DT->undefined_variables(0) ) ;
      }
      if( VARYING_DT->data_type() != MAC_Data::Double )
      {
         MAC_Error::object()->raise_bad_data_type(
                        exp, "time_step", MAC_Data::Double ) ;
      }
      DT.set( VARYING_DT->to_double() ) ;
   }
   cste_dt->destroy() ; cste_dt = 0 ;
   
   INI_DT = DT(0) ;
}




//----------------------------------------------------------------------
FV_TimeIterator:: ~FV_TimeIterator( void )
//----------------------------------------------------------------------
{
   MAC_Communicator const* com = MAC_Exec::communicator() ;
   com->barrier() ;
   if( com->rank() == 0 )
   {
      MAC_System::erase( CHECKPOINT_FILE ) ;
   }
}




//----------------------------------------------------------------------
size_t
FV_TimeIterator:: initial_iteration_number( void ) const
//----------------------------------------------------------------------
{
   return( ITER_INIT ) ;
}




//----------------------------------------------------------------------
double
FV_TimeIterator:: initial_time( void ) const
//----------------------------------------------------------------------
{
   return( T_INIT ) ;
}




//----------------------------------------------------------------------
double
FV_TimeIterator:: final_time( void ) const
//----------------------------------------------------------------------
{
   return( T_END ) ;
}




//----------------------------------------------------------------------
size_t
FV_TimeIterator:: storage_depth( void ) const
//----------------------------------------------------------------------
{
   return( DT.size() ) ;
}




//----------------------------------------------------------------------
void
FV_TimeIterator:: start( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIterator:: start" ) ;
   MAC_CHECK_PRE( !is_started() || is_finished() ) ;
   
   STARTED = true ;
   FINISHED_REASON = NotFinished ;
   ITER = ITER_INIT ;
   T = T_START ;
   set_time_step() ;
   start_checkpointing() ;

   MAC_CHECK_POST( is_started() ) ;
   MAC_CHECK_POST( !is_finished() ) ;
   MAC_CHECK_POST( iteration_number() == initial_iteration_number() ) ;
   MAC_CHECK_POST( FORMAL( time() == initial_time()+time_step() ) ) ;
   MAC_CHECK_POST( !just_went_back() ) ;
   MAC_CHECK_POST( !just_went_back_and_forth() ) ;
}




//----------------------------------------------------------------------
void
FV_TimeIterator:: go_next_time( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIterator:: go_next_time" ) ;
   MAC_CHECK_PRE( is_started() ) ;
   MAC_SAVEOLD( double, time, time() ) ;
   MAC_SAVEOLD( double, time_step, time_step() ) ;
   MAC_SAVEOLD( size_t, iteration_number, iteration_number() ) ;
   MAC_SAVEOLD( bool, just_went_back, just_went_back() ) ;
   
   test_checkpoint() ;
   
   // Condition to stop exactly at the prescribed time: 
   // 1) works well for a constant dt 
   // 2) and should work unless the time step abruptly changed by more 
   // than a factor of 1000
   if( ( FINISHED_REASON == NotFinished ) 
   	&& ( ( MAC::abs( T-T_END ) < 1.e-3 * DT(1) ) || T>T_END ) )
   {
      FINISHED_REASON = FinalTimeReached ;
   }
   
   if( FINISHED_REASON == NotFinished )
   {
      ITER++ ;
      set_time_step() ;
   }
   
   MAC_CHECK_POST( IMPLIES( is_finished(),
                            iteration_number() == OLD(iteration_number) ) ) ;
   MAC_CHECK_POST( IMPLIES( !is_finished(),
                            iteration_number() == OLD(iteration_number)+1 ) ) ;
   MAC_CHECK_POST( IMPLIES( is_finished(), time() == OLD(time) ) ) ;
   MAC_CHECK_POST( IMPLIES( is_finished(),
                            time_step() == OLD(time_step) ) ) ;
   MAC_CHECK_POST( IMPLIES( !is_finished(),
                            FORMAL( time() == OLD(time) + time_step() ) ) ) ;
   MAC_CHECK_POST( IMPLIES( !is_finished(),
                            !just_went_back() ) ) ;
   MAC_CHECK_POST( IMPLIES( !is_finished(),
		just_went_back_and_forth() == OLD(just_went_back) ) ) ;
   
}




//----------------------------------------------------------------------
bool
FV_TimeIterator:: is_started( void ) const
//----------------------------------------------------------------------
{
   return( STARTED ) ;
}




//----------------------------------------------------------------------
bool
FV_TimeIterator:: is_finished( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIterator:: is_finished" ) ;
   MAC_CHECK_PRE( is_started() ) ;
   
   bool result = ( FINISHED_REASON != NotFinished ) ;
   
   MAC_CHECK_POST( EQUIVALENT( result, ( FINISHED_REASON != NotFinished ) ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
double
FV_TimeIterator:: time( void ) const
//----------------------------------------------------------------------
{
   return( T ) ;
}




//----------------------------------------------------------------------
double
FV_TimeIterator:: time_step( size_t level ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIterator:: time_step" ) ;
   MAC_CHECK_PRE( level<storage_depth() ) ;
   return( DT(level) ) ;
}




//----------------------------------------------------------------------
size_t
FV_TimeIterator:: iteration_number( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIterator:: iteration_number" ) ;
   MAC_CHECK_PRE( is_started() ) ;
   return( ITER ) ;
}




//----------------------------------------------------------------------
bool
FV_TimeIterator:: just_went_back_and_forth( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIterator:: just_went_back_and_forth" ) ;
   MAC_CHECK_PRE( is_started() ) ;
   return( BACK_AND_FORTH_IT ) ;
}




//----------------------------------------------------------------------
FV_TimeIterator::FinishedReason
FV_TimeIterator:: finished_reason( void ) const
//----------------------------------------------------------------------
{
   return( FINISHED_REASON ) ;
}




//----------------------------------------------------------------------
void
FV_TimeIterator:: finish_iterations( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIterator:: finish_iterations" ) ;
   MAC_CHECK_PRE( is_started() ) ;
   MAC_CHECK_PRE( !is_finished() ) ;
   
   FINISHED_REASON = FinishIterationsCalled ;
   
   MAC_CHECK_POST( is_finished() ) ;
   MAC_CHECK_POST( finished_reason() == FinishIterationsCalled ) ;
}




//----------------------------------------------------------------------
bool
FV_TimeIterator:: time_step_is_fixed( void ) const
//----------------------------------------------------------------------
{
   return( VARYING_DT != 0 ) ;
}




//----------------------------------------------------------------------
void
FV_TimeIterator:: set_next_time_step( double dt )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIterator:: set_next_time_step" ) ;
   MAC_CHECK_PRE( is_started() ) ;
   MAC_CHECK_PRE( !is_finished() ) ;
   MAC_CHECK_PRE( !time_step_is_fixed() ) ;
   MAC_CHECK_PRE( dt > 0. ) ;

   NEXT_DT = ( dt == DT(0) ? MAC::bad_double() : dt ) ;

   MAC_CHECK_POST( next_time_step() == dt ) ;
}




//----------------------------------------------------------------------
double
FV_TimeIterator:: next_time_step( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIterator:: next_time_step" ) ;
   double result = ( NEXT_DT == MAC::bad_double() ? DT(0) : NEXT_DT ) ;
   MAC_CHECK_POST( result > 0. ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
void
FV_TimeIterator:: go_back( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIterator:: go_back" ) ;
   MAC_CHECK_PRE( is_started() ) ;
   MAC_CHECK_PRE( !just_went_back() ) ;
   MAC_CHECK_PRE( !is_finished() ) ;
   MAC_CHECK_PRE( !time_step_is_fixed() ) ;
   MAC_SAVEOLD( double, time_step, time_step() ) ;
   
   T = T-DT(0) ;
   RESTART_IT = true ;
 
   MAC_CHECK_POST( !is_finished() ) ;
   MAC_CHECK_POST( time_step() == OLD(time_step) ) ;
   MAC_CHECK_POST( FORMAL( time() == OLD(time) - time_step() ) ) ;
   MAC_CHECK_POST( just_went_back() ) ;
}




//----------------------------------------------------------------------
bool
FV_TimeIterator:: just_went_back( void ) const
//----------------------------------------------------------------------
{
   return( RESTART_IT ) ;
}




//----------------------------------------------------------------------
bool
FV_TimeIterator:: greater_or_equal( double time1, double time2 )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIterator:: greater_or_equal" ) ;

   double const epsilon = ( time2>0. ? 1.E-8 : -1.E-8 ) ;
   bool const result = ( time1 >= time2*( 1.0 - epsilon ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
bool
FV_TimeIterator:: table_of_times_is_valid( doubleVector const& times ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIterator:: table_of_times_is_valid" ) ;

   bool result = true ;
   size_t const n = times.size() ;
   if( n>0 )
   {
      result = greater_or_equal( times(0), initial_time() )
            && greater_or_equal( final_time(), times(n-1) ) ;
      for( size_t i=0 ; result && i<n-1 ; ++i )
      {
         result = ( times(i+1)>times(i) ) ; // increasing values
      }
   }
   return( result ) ;
}




//----------------------------------------------------------------------
void
FV_TimeIterator:: raise_invalid_table_of_times( 
                                         std::string const& class_name,
                                         std::string const& keyword,
                                         doubleVector const& times ) const 
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIterator:: raise_invalid_table_of_times" ) ;
   MAC_CHECK_PRE( !table_of_times_is_valid( times ) ) ;

   std::ostringstream mesg ;
   mesg << "*** " << class_name << " error:" << endl ;
   mesg << "    invalid data of keyword \""<< keyword <<"\"" << endl ;
   mesg << "    the time sequence must " << endl ;
   mesg << "      - start after the initial time (which is: " 
        << initial_time() << ")" << endl ;
   mesg << "      - stop before the final time (which is: " 
        << final_time() << ")" << endl ;
   mesg << "      - be sorted by increasing values" ;
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}




//----------------------------------------------------------------------
double
FV_TimeIterator:: next_time_in_table( doubleVector const& times ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIterator:: next_time_in_table" ) ;
   MAC_CHECK_PRE( table_of_times_is_valid( times ) ) ;
   
   double result = MAC::max_double() ;
   for( size_t i=0 ; i<times.size() ; i++ )
   {
      if( !greater_or_equal( time(), times(i) ) )
      {
         result = times(i) ;
         break ;
      }
   }
   
   MAC_CHECK_POST( !greater_or_equal( time(), result ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
void
FV_TimeIterator:: save_state( MAC_ObjectWriter* writer ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIterator:: save_state" ) ;
   MAC_CHECK_PRE( save_state_PRE( writer ) ) ;
   
   writer->start_new_object( "FV_TimeIterator" ) ;

   // Saving current time
   writer->add_entry( "t", MAC_Double::create( 0, T ) ) ;
   
   // Saving current iter
   writer->add_entry( "iter", MAC_Int::create( 0, ITER ) ) ;
   
   // Initial time step
   writer->add_entry( "initial_dt", MAC_Double::create( 0, INI_DT ) ) ;

   // Time step has been modified
   if( NEXT_DT != MAC::bad_double() )
   {
      writer->add_entry( "new_dt", MAC_Double::create( 0, NEXT_DT ) ) ;
   }
   
   // Restart
   if( RESTART_IT )
   {
     writer->add_entry( "just_went_back", MAC_Bool::create( 0, true ) ) ;
   }
   
   // Restart
   if( BACK_AND_FORTH_IT )
   {
     writer->add_entry( "just_went_back_and_forth", 
     	MAC_Bool::create( 0, true ) ) ;
   }
   
   // Saving current time step
   writer->add_entry( "dt", MAC_DoubleVector::create( 0, DT ) ) ;
   
   writer->finalize_object() ;

   MAC_CHECK_POST( save_state_POST( writer ) ) ;
}




//----------------------------------------------------------------------
void
FV_TimeIterator:: restore_state( MAC_ObjectReader* reader )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIterator:: restore_state" ) ;
   MAC_CHECK_PRE( restore_state_PRE( reader ) ) ;

   reader->start_object_retrieval( "FV_TimeIterator" ) ;
   
   T = reader->data_of_entry( "t" )->to_double() ;
   T_START = T ;
   ITER_INIT = reader->data_of_entry( "iter" )->to_int()+1 ;
   NEXT_DT = MAC::bad_double() ;
   if( !VARYING_DT )
   {
      double ini_dt = reader->data_of_entry( "initial_dt" )->to_double() ;
      if( ini_dt != INI_DT ) // time step is changed in the data deck
      {
         NEXT_DT = INI_DT ;
      }
      else if( reader->has_entry( "new_dt" ) )
      {
         NEXT_DT = reader->data_of_entry( "new_dt" )->to_double() ;
      }
   }
   RESTART_IT = reader->has_entry( "just_went_back" ) ;
   BACK_AND_FORTH_IT = reader->has_entry( "just_went_back_and_forth" ) ;
   DT = reader->data_of_entry( "dt" )->to_double_vector() ;
   ITER = MAC::bad_index() ;
   STARTED = false ;
   FINISHED_REASON = NotFinished ;
   reader->end_object_retrieval() ;
   MAC_CHECK_POST( restore_state_POST( reader ) ) ;
}




//----------------------------------------------------------------------
void
FV_TimeIterator:: set_time_step( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIterator:: set_time_step" ) ;

   if( !RESTART_IT )
   {
      for( size_t i=DT.size()-1 ; i>0 ; --i )
      {
         DT( i ) = DT( i-1 ) ;
      }
   }
   if( NEXT_DT != MAC::bad_double() )
   {
      DT(0) = NEXT_DT ;
      NEXT_DT =  MAC::bad_double() ;
   }
   else if( VARYING_DT != 0 )
   {
      TIME->set( T ) ;
      DT(0) = VARYING_DT->to_double() ;
   }
   T += DT(0) ;
   BACK_AND_FORTH_IT = RESTART_IT ;
   RESTART_IT = false ;
}




//------------------------------------------------------------------------
void
FV_TimeIterator:: start_checkpointing( void )
//------------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIterator:: start_checkpointing" ) ;
   
   MAC_Communicator const* com = MAC_Exec::communicator() ;
   size_t rank = com->rank() ;
   int pid ;
   if( rank == 0 )
   {
      pid = MAC_System::process_id() ;
   }
   com->broadcast( pid, 0 ) ;
   std::ostringstream nstr ;
   nstr << "remove_to_stop_" << pid ;
   CHECKPOINT_FILE = nstr.str() ;
   if( rank == 0 )
   {
      std::ofstream file( CHECKPOINT_FILE.c_str(), std::ios::trunc ) ;
      file.close() ;
   }
}




//------------------------------------------------------------------------
void
FV_TimeIterator:: test_checkpoint( void )
//------------------------------------------------------------------------
{
   MAC_LABEL( "FV_TimeIterator:: test_checkpoint" ) ;
   
   std::ifstream file( CHECKPOINT_FILE.c_str() ) ;
   if( !file ) FINISHED_REASON = UserCheckpoint ;
}
