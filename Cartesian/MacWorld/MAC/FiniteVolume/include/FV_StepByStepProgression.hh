#ifndef FV_STEP_BY_STEP_PROGRESSION_HH
#define FV_STEP_BY_STEP_PROGRESSION_HH

#include <MAC_Application.hh>
#include <MAC.hh>

#include <doubleVector.hh>

class MAC_ModuleExplorer ;
class MAC_Timer ;
class FV_OneStepIteration ;
class FV_TimeIteratorAdapter ;
class FV_TimeIterator ;
class FV_DomainAndFields ;
class FV_PostProcessingWriter ;



class FV_StepByStepProgression : public MAC_Application
{
    public: //----------------------------------------------------------

    //-- Instance delivery and initialization

      static FV_StepByStepProgression* create( 
	MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
	double const& initial_time = 0. ) ;

    //-- Instance characteristics

      FV_TimeIterator const* time_iterator( void ) const ;

      virtual void run( std::string const& inputRestartFileName = 
      	MAC::undefined_string ) ;
      
      /**
        @brief Perform additional post-processing for the special purpose
	of the application FV_MorePostProcessing that calls this method
        @param exp Pointer to the module explorer
        @remarks Called by FV_MorePostProcessing::run
        @par Content 
            - Calls FV_OneStepIteration::do_more_post_processing
            - If verbose mode, calls 
                FV_StepByStepProgression::display_memory_usage_proc_by_proc
      */
      void do_more_post_processing( MAC_ModuleExplorer const* exp ) ;              

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~FV_StepByStepProgression( void ) ;
      FV_StepByStepProgression( FV_StepByStepProgression const& other ) ;
      FV_StepByStepProgression& operator=( 
                                FV_StepByStepProgression const& other ) ;

      FV_StepByStepProgression( MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
	double const& initial_time ) ;

   //-- Plug in

      FV_StepByStepProgression( void ) ;

      virtual FV_StepByStepProgression* create_replica( 
	MAC_Object* a_owner,
	MAC_ModuleExplorer const* exp,
	double const& initial_time = 0. ) const ;

   //-- Internals

      void save_for_restart( bool const& force ) ;

      void save_for_post_processing( bool const& force ) ;
      
      void display_timers( size_t indent_width ) const ;
      
      void display_timers_proc_by_proc( size_t indent_width ) const ;      

      void do_save_for_post(
      		FV_PostProcessingWriter* mpp_writer ) ;
      
      void check_times( std::string const& list_name,
                        doubleVector const& dates ) const ;
      
      static unsigned long long int display_memory_usage(
		std::ostream& os, size_t indent_width ) ;
				  
      void display_memory_usage_proc_by_proc(
		std::ostream& os, size_t indent_width ) const ;

   //-- Persistence
      
      virtual void add_storable_objects( MAC_ListIdentity* list ) ;

   //-- Preconditions, Postconditions, Invariant      

      virtual bool invariant( void ) const ;

   //-- Class attributes

      static FV_StepByStepProgression const* PROTOTYPE ;

      static size_t graphics_level ;

   //-- Attributes

      FV_TimeIterator* TIME_IT ;
      
      FV_OneStepIteration* ONE_IT ;
      FV_DomainAndFields* DOM ;

      doubleVector GRAPHICS_TIMES ;
      double GRAPHICS_NEXT_TIME ;

      doubleVector SAVER_TIMES ;
      double SAVER_NEXT_TIME ;

      MAC_Timer* overall ;
      MAC_Timer* POST_TIMER ;
      MAC_Timer* SAVE_TIMER ;
      bool SAVEFG ;
      size_t LAST_MEM_IT ;

      bool VERBOSE ;
      
      size_t POST_CYCLE_NUMBER;
      
      size_t RANK ;
} ;


#endif

