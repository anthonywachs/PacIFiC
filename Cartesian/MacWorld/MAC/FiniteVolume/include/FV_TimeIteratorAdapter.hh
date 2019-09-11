#ifndef FE_TIME_ADAPTATOR_HH
#define FE_TIME_ADAPTATOR_HH

#include <MAC_Object.hh>

class FV_DomainAndFields ;

class PDE_SetOfDiscreteFields ;
class MAC_ModuleExplorer ;
class MAC_ObjectRegister ;

class FV_TimeIterator ;

/*
Servers for adaptation of the time iterations managed by `FV_TimeIterator::'
objects.

Available adaptations are:
   - modification of the time step:
       * several external objects can propose a new time step calling
         the function `::propose_next_time_step()'
       * the resulting new time step is set to `::time_iterator()' by calling
         the function `::adapt_time_iterator()'

   - normal end of the time iterations (for example when steady state is
     reached before the final time of `::time_iterator()') by calling
     the function `::propose_to_finish_iterations()' and then applying
     this change to `::time_iterator() calling `::adapt_time_iterator()'

   - failure of the time iterations (for example after no convergence of an
     inner solver) by calling `::set_time_iteration_failed'

Remarks:
   - A pluggable factory allows the user to define adapters that
     are dynamically choosen according to the data deck.
   - A default behavior is proposed.

FRAMEWORK INSTANTIATION

   1. Derive a concrete subclass, say MyTimeIteratorAdapter.
   2. Choose a name for MyTimeIteratorAdapter, say "MyTimeIteratorAdapter".
   3. Implement a private destructor.
   4. Declare all constructors private.
   5. Define the prototype to be registered :
      5.1 Implement a default constructor that initializes the
          `FV_TimeIteratorAdapter::' subobject by calling
               `FV_TimeIteratorAdapter( std::string const& )'
          with "MyTimeIteratorAdapter" as argument.
      5.2 Define and initialize a static instance by calling the default
          constructor.
   6. Implement a private constructor that initializes the 
      `FV_TimeIteratorAdapter::' by calling
                     `FV_TimeIteratorAdapter( FV_TimeIterator* )'
   7. Implement the `::create_replica' methods that allocate an object
      of type `MyTimeIteratorAdapter' initialized using the private constructor
      described above, and subsequently return a pointer to that object.
   8. Implement relevant virtual functions.

PUBLISHED*/

class FV_TimeIteratorAdapter : public MAC_Object
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return the concrete instance of `FV_TimeIteratorAdapter::'
      // linked to `t_it' and of name `exp'->string_data( "concrete_name" )
      // for single domain applications.
      static FV_TimeIteratorAdapter* make(
                                FV_TimeIterator* t_it,
                                FV_DomainAndFields const* dom,
                                MAC_ModuleExplorer const* exp ) ;
      
      // Create and return the concrete instance of `FV_TimeIteratorAdapter::'
      // linked to `t_it' and of name `exp'->string_data( "concrete_name" )
      // for multi domains applications.
      static FV_TimeIteratorAdapter* make(
                                FV_TimeIterator* t_it,
                                MAC_ModuleExplorer const* exp ) ;

      // Create and return the default `FV_TimeIteratorAdapter::'
      // linked to `t_it' for single domain applications.
      static FV_TimeIteratorAdapter* make_default(
                                FV_TimeIterator* t_it,
                                FV_DomainAndFields const* dom ) ;
      
      // Create and return the default `FV_TimeIteratorAdapter::'
      // linked to `t_it' for multi domains applications.
      static FV_TimeIteratorAdapter* make_default(
                                FV_TimeIterator* t_it ) ;

   //-- Associated FV_TimeIterator object(1.0)

      FV_TimeIterator const* time_iterator( void ) const ;
      
   //-- Adaptation of the associated FV_TimeIterator object(2.0)

      // Initialize for a new time step.
      void initialize_time_step( void ) ;

      // Adapt new time iterator.
      void adapt_time_iterator( void ) ;
      
      // IMPLEMENTATION : raise error and terminate program.
      virtual void set_time_iteration_failed( void ) ;

   //-- Proposed modifications to the associated FV_TimeIterator object(2.2)

      // Propose `dt' for the next time step of `::time_iterator()' 
      // (`::time_iterator()' is not modified and the proposed modification 
      // will be effective only after calling `::adapt_time_iterator()').
      void propose_next_time_step( double dt ) ;
      
      // value currently proposed for the next time step of `::time_iterator()'
      double next_time_step( void ) const ;

      // Propose to terminate current time iterations
      // (`::time_iterator()' is not modified and the proposed modification 
      // will be effective only after calling `::adapt_time_iterator()').
      void propose_to_finish_iterations( void ) ;
      
      // Are the current time iterations terminated ?
      bool iterations_are_finished( void ) const ;
      
   //-- Persistence

      virtual void save_state( MAC_ObjectWriter* writer ) const ;

      virtual void restore_state( MAC_ObjectReader* reader ) ;
      
   protected: //--------------------------------------------------------
      
   //-- Plug in

      virtual ~FV_TimeIteratorAdapter( void ) ;

      FV_TimeIteratorAdapter( std::string const& a_concrete_name ) ;

      FV_TimeIteratorAdapter( FV_TimeIterator* t_it ) ;
      
      virtual FV_TimeIteratorAdapter* create_replica( 
                                     FV_TimeIterator* t_it,
                                     FV_DomainAndFields const* dom,
                                     MAC_ModuleExplorer const* exp ) const ;
      virtual FV_TimeIteratorAdapter* create_replica( 
                                     FV_TimeIterator* t_it,
                                     MAC_ModuleExplorer const* exp ) const ;
      
      bool is_a_prototype( void ) const ;

   //-- Default adaptator

      FV_TimeIteratorAdapter( void ) ;
      virtual FV_TimeIteratorAdapter* create_replica(
                                     FV_TimeIterator* t_it,
                                     FV_DomainAndFields const* dom ) const ;
      virtual FV_TimeIteratorAdapter* create_replica(
                                     FV_TimeIterator* t_it ) const ;
      
   //-- Adaptation of the associated FV_TimeIterator object

      // Inner initialization called by `::initialize_time_step'().
      // IMPLEMENTATION : do nothing.
      virtual void initialize_inner( void ) ;
      
      /*
      Indicate the desired parameters for the next time iteration.
      On exit, the meaning of the arguments is the following :
         - `restart'=true: restart current iteration with 
                          `next_dt' as time step
         - `finished'=true: time iterations are finished
         - `finished'=false: perform next iteration with time step `next_dt'.
      IMPLEMENTATION : arguments are unchanged.
      */
      virtual void define_parameters_for_next_iteration(
                           bool& finished, bool& restart, double& next_dt ) ;
      
      // new time step from the one proposed by the client
      // IMPLEMENTATION : `dt' is returned without modification
      virtual double next_time_step_from_proposed_one( double dt ) const ;

      // Restart time iteration with `dt' as new time step :
      // should be called by `::set_time_iteration_failed'().
      void restart_iteration_with_new_time_step( double dt ) ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

      virtual bool set_time_iteration_failed_PRE( void ) const ;
      virtual bool set_time_iteration_failed_POST( void ) const ;

      virtual bool initialize_inner_PRE( void ) const ;
      virtual bool initialize_inner_POST( void ) const ;
      
      virtual bool define_parameters_for_next_iteration_PRE( 
                       bool finished, bool restart, double next_dt ) const ;
      virtual bool define_parameters_for_next_iteration_POST(
                       bool finished, bool restart, double next_dt ) const ;
      
      virtual bool next_time_step_from_proposed_one_PRE( double dt ) const ;
      virtual bool next_time_step_from_proposed_one_POST(
                                          double result, double dt ) const ;
      
      virtual bool create_replica_PRE(
                               FV_TimeIterator const* t_it,
                               FV_DomainAndFields const* dom,
                               MAC_ModuleExplorer const* exp ) const ;
      virtual bool create_replica_POST(
                               FV_TimeIteratorAdapter const* result,
                               FV_TimeIterator const* t_it,
                               FV_DomainAndFields const* dom,
                               MAC_ModuleExplorer const* exp ) const ;
      
      virtual bool create_replica_PRE(
                               FV_TimeIterator const* t_it,
                               MAC_ModuleExplorer const* exp ) const ;
      virtual bool create_replica_POST(
                               FV_TimeIteratorAdapter const* result,
                               FV_TimeIterator const* t_it,
                               MAC_ModuleExplorer const* exp ) const ;
      
      virtual bool create_replica_PRE(
                               FV_TimeIterator const* t_it,
                               FV_DomainAndFields const* dom ) const ;
      virtual bool create_replica_POST(
                               FV_TimeIteratorAdapter const* result,
                               FV_TimeIterator const* t_it,
                               FV_DomainAndFields const* dom ) const ;
      
      virtual bool create_replica_PRE(
                               FV_TimeIterator const* t_it ) const ;
      virtual bool create_replica_POST(
                               FV_TimeIteratorAdapter const* result,
                               FV_TimeIterator const* t_it ) const ;

      
   private: //----------------------------------------------------------
      
      FV_TimeIteratorAdapter( FV_TimeIteratorAdapter const& other ) ;
      FV_TimeIteratorAdapter& operator=(
                              FV_TimeIteratorAdapter const& other ) ;

      static MAC_ObjectRegister* plugins_map( void ) ;

   //-- Static attributes

      static FV_TimeIteratorAdapter const* DEFAULT_PROTOTYPE ;

   //-- Attributes

      bool const IS_PROTO ;
      FV_TimeIterator* const T_IT ;
      double DT ;
      bool FINISHED ;
} ;

#endif
