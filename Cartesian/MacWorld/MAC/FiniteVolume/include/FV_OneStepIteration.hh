#ifndef FV_ONE_STEP_ITERATION_HH
#define FV_ONE_STEP_ITERATION_HH

#include <MAC_Object.hh>
#include <map>
#include <set>

class MAC_Communicator ;
class MAC_ModuleExplorer ;
class MAC_ObjectRegister ;
class MAC_ObjectWriter ;
class MAC_Timer ;
class FV_TimeIteratorAdapter ;
class FV_TimeIterator ;
class FV_DomainAndFields ;


class FV_OneStepIteration : public MAC_Object
{
   public: //-----------------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return an instance according to the data attainable by
      // `exp', devoted to perform computation associated to a system
      // of partial differential equations whose numerical solution requires
      // objects delivered by `dom'.
      static FV_OneStepIteration* make( MAC_Object* a_owner,
      			FV_DomainAndFields const* dom,
			MAC_ModuleExplorer* exp ) ;

   //-- Substeps of the step by step progression

      // Before starting time marching in `FV_StepByStepProgression::run',
      // perform initial computations.
      // IMPLEMENTATION : do nothing.
      virtual void do_before_time_stepping( FV_TimeIterator const* t_it, 
      	std::string const& basename ) ;

      // Within a step of the time marching procedure in
      // `FV_StepByStepProgression::run', perform an initial stage
      // just before start of the inner iterations.
      // IMPLEMENTATION : do nothing.
      virtual void do_before_inner_iterations_stage(
                                           FV_TimeIterator const* t_it ) ;

      virtual void do_inner_iterations_stage( FV_TimeIterator const* t_it ) ;

      bool inner_iterations_stage_failed( void ) const ;

      // Within a step of the time marching procedure in
      // `FV_StepByStepProgression::run', perform an inner iteration stage.
      virtual void do_one_inner_iteration( FV_TimeIterator const* t_it ) = 0;

      // Are inner iterations performed in `::do_one_inner_iteration'
      // completed ?
      // IMPLEMENTATION : true
      virtual bool inner_iterations_are_completed(
                                     FV_TimeIterator const* t_it ) const ;

      // Within a step of the time marching procedure in
      // `FV_StepByStepProgression::run', perform a final stage just
      // after completion of the inner iterations.
      // IMPLEMENTATION : do nothing.
      virtual void do_after_inner_iterations_stage(
                                     FV_TimeIterator const* t_it ) ;

      virtual void do_after_time_adaptation(
                                     FV_TimeIterator const* t_it ) ;

      // After completion of the time marching procedure in
      // `FV_StepByStepProgression::run', perform final computations.
      // IMPLEMENTATION : do nothing.
      virtual void do_after_time_stepping( void ) ;

   //-- Elapsed times

      static void reset_standard_times( void ) ;

      static void print_standard_times( std::ostream& os,
                                        size_t indent_width ) ;

      virtual void print_additional_times( std::ostream& os,
                                           size_t indent_width ) const ;

   //-- Time iterator modification

      virtual void adapt_time_iterator( FV_TimeIteratorAdapter* t_adapter ) ;

      virtual void notify_inner_iterations_stage_failure( void ) ;

      virtual void reset_after_inner_iterations_stage_failure( void ) ;

   //-- Savings for post-processing and restart

      // Use `rs' to save data for post processing other than time and fields 
      // at the current time step (attainable through `t_it') in the time 
      // marching procedure of `FV_StepByStepProgression::run'.
      // IMPLEMENTATION : do nothing.
      virtual void do_additional_savings( FV_TimeIterator const* t_it,
      	int const& cycleNumber ) ;
	
      // Save other data than FV fields for restart 
      // at the current time step (attainable through `t_it') in the time 
      // marching procedure of `FV_StepByStepProgression::run'.
      // IMPLEMENTATION : do nothing.
      virtual void do_additional_save_for_restart( FV_TimeIterator const* t_it,
      	size_t const& restartCycleNumber, std::string const& basename ) ;	

   //-- Persistence

      // Extend `list' so that it contains all objects required by the
      // storage and retrieval mechanisms concerning the `FV_OneStepIteration::'
      // base class subobject, then call `::add_storable_objects'.
      void register_storable_objects( MAC_ListIdentity* list ) ;

      // Extend `list' so that it contains all objects required by the
      // storage and retrieval mechanisms that are not part of the
      // `FV_OneStepIteration::' base class subobject.
      // IMPLEMENTATION : do nothing, i.e. leave `list' unchanged.
      virtual void add_storable_objects( MAC_ListIdentity* list ) ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;

      static std::string const& indent( void ) ;

      static void increase_indent( void ) ;

      static void decrease_indent( void ) ;
      
   //-- Additional post-processing
      /**
        @brief Perform additional post-processing for the special purpose
	of the application FV_MorePostProcessing 
        @param dom Pointer to all fields
        @param exp Pointer to module explorer
        @remarks Called by FV_StepByStepProgression::do_more_post_processing
        @remarks Called by FV_SplitSystem::do_more_post_processing
        @par Content 
          EMPTY VIRTUAL METHOD (see in children's method)
      */
      virtual void do_more_post_processing( FV_DomainAndFields* dom,
      		MAC_ModuleExplorer const* exp ) ;

   protected: //--------------------------------------------------------------

      virtual ~FV_OneStepIteration( void ) ;

      // Construction of an instance whose owner is `a_owner'.
      // On exit, `self' still refers to the object identified by `dom'
      FV_OneStepIteration( MAC_Object* a_owner,
 		FV_DomainAndFields const* dom,     
		MAC_ModuleExplorer const* exp ) ;

   //-- Plug in

      // for prototype registration only
      FV_OneStepIteration( std::string const& name ) ;

      virtual FV_OneStepIteration* create_replica(
		MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer* exp ) const = 0;

      bool is_a_prototype( void ) const ;

   //-- Distributed processing

      static MAC_Communicator const* communicator( void ) ;

   //-- Timer and screen output

      size_t verbose_level( void ) const ;

      void start_assembling_timer( bool silent = false ) const ;
      void stop_assembling_timer( bool silent = false ) const ;

      void start_solving_timer( bool silent = false ) const ;
      void stop_solving_timer( bool silent = false ) const ;

      void start_total_timer( std::string const& mesg,
                              bool silent = false ) const ;
      void stop_total_timer( bool silent = false ) const ;


   //-- Preconditions, Postconditions, Invariant

      bool do_before_time_stepping_PRE( FV_TimeIterator const* t_it ) const ;

      bool do_before_inner_iterations_stage_PRE(
                                       FV_TimeIterator const* t_it ) const ;
      bool do_before_inner_iterations_stage_POST(
                                       FV_TimeIterator const* t_it ) const ;

      bool do_inner_iterations_stage_PRE(
                                       FV_TimeIterator const* t_it ) const ;

      bool do_one_inner_iteration_PRE( FV_TimeIterator const* t_it ) const ;

      bool inner_iterations_are_completed_PRE(
                                       FV_TimeIterator const* t_it ) const ;

      bool do_after_inner_iterations_stage_PRE(
                                       FV_TimeIterator const* t_it ) const ;

      bool do_after_time_adaptation_PRE(
                                       FV_TimeIterator const* t_it ) const ;

      bool notify_inner_iterations_stage_failure_PRE( void ) const ;
      bool notify_inner_iterations_stage_failure_POST( void ) const ;

      bool reset_after_inner_iterations_stage_failure_PRE( void ) const ;
      bool reset_after_inner_iterations_stage_failure_POST( void ) const ;

      bool adapt_time_iterator_PRE(
                           FV_TimeIteratorAdapter const* t_adapter ) const ;

      bool do_additional_savings_PRE( FV_TimeIterator const* t_it ) const ;

      bool add_storable_objects_PRE( MAC_ListIdentity* list ) const ;

      bool create_replica_PRE( MAC_Object* a_owner,
      		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer const* exp ) const ;
      bool create_replica_POST( FV_OneStepIteration const* result,
		MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer const* exp ) const ;


   private: //----------------------------------------------------------------

      FV_OneStepIteration( void ) ;
      FV_OneStepIteration( FV_OneStepIteration const& other ) ;
      FV_OneStepIteration& operator=( FV_OneStepIteration const& other ) ;

      static MAC_ObjectRegister* plugins_map( void ) ;

      void start_internal_timer( std::map< std::string, MAC_Timer* >& timers,
                                 std::string const& display,
                                 bool do_display ) const ;
      void stop_internal_timer(
                         std::map< std::string, MAC_Timer* > const& timers,
                         bool do_display ) const ;

   //-- Class attributes

      static std::string INDENT ;
      static std::map< std::string, MAC_Timer* > ASS_TIMER ;
      static std::map< std::string, MAC_Timer* > SOL_TIMER ;
      static std::map< std::string, MAC_Timer* > TOT_TIMER ;

   //-- Attributes

      bool IS_PROTO ;
      bool FAILURE ;
      size_t VERBOSE_LEVEL ;
      FV_DomainAndFields const* DOMAIN ;
      size_t RANK ;

} ;

#endif
