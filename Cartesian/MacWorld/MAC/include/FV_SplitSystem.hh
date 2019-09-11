#ifndef FV_SPLIT_SYSTEM
#define FV_SPLIT_SYSTEM

#include <FV_OneStepIteration.hh>

class MAC_List ;
class MAC_ListIterator ;

/*
PUBLISHED
*/

class FV_SplitSystem : public FV_OneStepIteration
{
   public: //------------------------------------------------------------

   //-- Substeps of the step by step progression

      // IMPLEMENTATION : Call `FE_OneStepIteration::do_before_time_stepping'
      // for all handled `FE_OneStepIteration::' instances.
      virtual void do_before_time_stepping( FV_TimeIterator const* t_it, 
      	std::string const& basename ) ;

      // IMPLEMENTATION : 
      // Call `FE_OneStepIteration::do_before_inner_iterations_stage'
      // for all handled `FE_OneStepIteration::' instances.
      virtual void do_before_inner_iterations_stage( 
                                           FV_TimeIterator const* t_it ) ;

      // IMPLEMENTATION : 
      // Call `FE_OneStepIteration::do_inner_iterations_stage'
      // for all handled `FE_OneStepIteration::' instances.
      virtual void do_inner_iterations_stage( FV_TimeIterator const* t_it ) ;

      // IMPLEMENTATION : raise a fatal error.
      virtual void do_one_inner_iteration( FV_TimeIterator const* t_it ) ;

      // IMPLEMENTATION : 
      // true if `FE_OneStepIteration::inner_iterations_are_completed' is true
      // for all handled `FE_OneStepIteration::' instances, false otherwise
      virtual bool inner_iterations_are_completed( 
                                      FV_TimeIterator const* t_it ) const ;

      // IMPLEMENTATION : 
      // Call `FE_OneStepIteration::do_after_inner_iterations_stage'
      // for all handled `FE_OneStepIteration::' instances.
      virtual void do_after_inner_iterations_stage(  
                                      FV_TimeIterator const* t_it ) ;

      // IMPLEMENTATION : 
      // Call `FE_OneStepIteration::do_after_time_adaptation'
      // for all handled `FE_OneStepIteration::' instances.
      virtual void do_after_time_adaptation(  
                                      FV_TimeIterator const* t_it ) ;

      // IMPLEMENTATION : 
      // Call `FE_OneStepIteration::do_after_time_stepping'
      // for all handled `FE_OneStepIteration::' instances.
      virtual void do_after_time_stepping( void ) ;

   //-- Elapsed times

      virtual void print_additional_times( std::ostream& os,
                                           size_t indent_width ) const ;

   //-- Time iterator modification

      virtual void adapt_time_iterator( FV_TimeIteratorAdapter* t_adapter ) ;

      virtual void notify_inner_iterations_stage_failure( void ) ;
      
      virtual void reset_after_inner_iterations_stage_failure( void ) ;

   //-- Savings for post-processing

      // IMPLEMENTATION : 
      // Call `FE_OneStepIteration::do_additional_savings'
      // for all handled `FE_OneStepIteration::' instances.
      virtual void do_additional_savings( FV_TimeIterator const* t_it,
      	int const& cycleNumber ) ;
	
      // Save other data than FV fields for restart 
      // at the current time step (attainable through `t_it') in the time 
      // marching procedure of `FV_StepByStepProgression::run'.
      // IMPLEMENTATION : do nothing.
      virtual void do_additional_save_for_restart( FV_TimeIterator const* t_it,
      	size_t const& restartCycleNumber, std::string const& basename ) ;
      
    //-- Persistence
      
      // IMPLEMENTATION : 
      // Call `FE_OneStepIteration::add_storable_objects'
      // for all handled `FE_OneStepIteration::' instances.
      virtual void add_storable_objects( MAC_ListIdentity* list ) ;
           
   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   //-- Additional post-processing
      
      /**
        @brief Perform additional post-processing for the special purpose
	of the application FV_MorePostProcessing 
        @param dom Pointer to all fields
        @param exp Pointer to the module explorer
        @remarks Called through it's mother class
          FV_OneStepIteration::do_more_post_processing (call, call-back !!!)
        @par Content 
          For each item of the MAC_List (i.e. each MAC_ListIterator),
            calls FV_OneStepIteration::do_more_post_processing
      */
      virtual void do_more_post_processing( FV_DomainAndFields* dom,
      		MAC_ModuleExplorer const* exp ) ;
            
   protected: //---------------------------------------------------------

   private: //-----------------------------------------------------------

     ~FV_SplitSystem( void ) ;
      FV_SplitSystem( FV_SplitSystem const& other ) ;
      FV_SplitSystem& operator=( FV_SplitSystem const& other ) ;

      FV_SplitSystem( MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer const* exp ) ;

   //-- Plug in

      FV_SplitSystem( void ) ;

      virtual FV_SplitSystem* create_replica( 
		MAC_Object* a_owner,
		FV_DomainAndFields const* dom,
		MAC_ModuleExplorer* exp ) const ;

   //-- Class attributes

      static FV_SplitSystem const* PROTOTYPE ;

   //-- Attributes

      MAC_List* CMPS ;
      MAC_ListIterator* IT ;
} ;

#endif
