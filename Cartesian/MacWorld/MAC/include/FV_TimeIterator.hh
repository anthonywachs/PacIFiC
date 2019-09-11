#ifndef FV_TIME_ITERATOR_HH
#define FV_TIME_ITERATOR_HH

#include <MAC_Object.hh>

#include <doubleArray2D.hh>
#include <doubleVector.hh>

class MAC_Data ;
class MAC_Double ;
class MAC_ModuleExplorer ;

/*
Time steppers.

PUBLISHED
*/

class FV_TimeIterator : public MAC_Object
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static FV_TimeIterator* create( MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
	double const& initial_time = 0. ) ;

   //-- Characteristics

      // initial iteration number
      size_t initial_iteration_number( void ) const ;

      // lower bound of the time interval
      double initial_time( void ) const ;

      // upper bound of the time interval
      double final_time( void ) const ;

      // time step storage depth
      size_t storage_depth( void ) const ;

   //-- Iteration

      // Start time stepping.
      void start( void ) ;

      // Advance one time step.
      void go_next_time( void ) ;

      // Is time stepping started ?
      bool is_started( void ) const ;

      // Is time stepping finished ?
      bool is_finished( void ) const ;

   //-- Access

      // current time
      double time( void ) const ;

      // time step of the past `level'-th iteration (eg current time step if 
      // `level' is zero, time step of the previous iteration if 
      // `level' is 1)
      double time_step( size_t level = 0 ) const ;

      // current iteration number
      size_t iteration_number( void ) const ;

      bool just_went_back_and_forth( void ) const ;
      
      enum FinishedReason
      {
         UserCheckpoint         = 1,
         FinishIterationsCalled = 2,
         FinalTimeReached       = 3,
         
         NotFinished = 0
      } ;
      
      FinishedReason finished_reason( void ) const ;

   //-- Adaptation(60.)
      
      void finish_iterations( void ) ;

      bool time_step_is_fixed( void ) const ;
      
      void set_next_time_step( double dt ) ;

      double next_time_step( void ) const ;
      
      void go_back( void ) ;

      bool just_went_back( void ) const ;

    //-- Utilities(70.)

      // Is `time1' greater than or equal to `time2' ?
      static bool greater_or_equal( double time1, double time2 ) ;

      /* 
      Is `times' a valid table of times, in the sense that:
          - items are sorted by increasing values,
          - all items are greater than (or equal to) `::initial_time()' and 
            lower than (or equal to) `::final_time()' ?
      */
      bool table_of_times_is_valid( doubleVector const& times ) const ;

      // Raise a fatal error due to the invalid table `times', this table
      // being the data of keyword `keyword' in a hierarchical data structure
      // read by the class of name `class_name'.
      void raise_invalid_table_of_times( std::string const& class_name,
                                         std::string const& keyword,
                                         doubleVector const& times ) const ;

      // first item of `times' that is `::greater_or_equal()' to `::time'()
      double next_time_in_table( doubleVector const& times ) const ;

    //-- Persistence

      virtual void save_state( MAC_ObjectWriter* writer ) const ;

      virtual void restore_state( MAC_ObjectReader* reader ) ;
      
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      FV_TimeIterator( void ) ;
     ~FV_TimeIterator( void ) ;
      FV_TimeIterator( FV_TimeIterator const& other ) ;
      FV_TimeIterator& operator=( FV_TimeIterator const& other ) ;

      FV_TimeIterator( MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
	double const& initial_time ) ;

   //-- Internals
      
      void set_time_step( void ) ;
      
      void start_checkpointing( void ) ;
      
      void test_checkpoint( void ) ;

   //-- Attributes
      
      bool STARTED ;
      FinishedReason FINISHED_REASON ;

      size_t ITER_INIT ;
      size_t ITER ;

      double T_INIT ;
      double T_START ;
      double T_END  ;
      double T ;

      MAC_Data* VARYING_DT ;
      MAC_Double* TIME ;
      bool RESTART_IT ;
      bool BACK_AND_FORTH_IT ;
      doubleVector  DT ;
      double INI_DT ;
      double NEXT_DT ;
      
      std::string CHECKPOINT_FILE ;
      
       
} ;

#endif
