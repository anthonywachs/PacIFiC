#ifndef MAC_TIMER_HH
#define MAC_TIMER_HH

#include <MAC_Object.hh>

#include <ctime>

/*
Timers for measuring and reporting the elapsed time passed
between start and stop events.

PUBLISHED
*/

class MAC_Timer : public MAC_Object
{
   public: //----------------------------------------------------------

   //-- Instance delivery and initialization

      static MAC_Timer* create( MAC_Object * a_owner ) ;
      
   //-- Commands(1.)
      
      // Start or restart self.
      void start( void ) ;
      
      // Stop self.
      void stop( void ) ;

      // Stop self and nullify cumulative time.
      void reset( void ) ;

   //-- Access
      
      // cumulative user time spent between 
      // the `::start' and the `::stop' method calls
      double time( void ) const ;

      // cumulative wall clock time spent between 
      // the `::start' and the `::stop' method calls
      double elapsed_time( void ) const ;

      // Is self running ?
      bool is_running( void ) const ;

   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
      static void print_time( double a_time, 
                              std::ostream& os, size_t indent_width ) ;

   protected: //------------------------------------------------------

   private: //--------------------------------------------------------

      MAC_Timer( void ) ;
     ~MAC_Timer( void ) ;
      MAC_Timer( MAC_Timer const& other ) ;
      MAC_Timer& operator=( MAC_Timer const& other ) ;

      MAC_Timer( MAC_Object* a_owner ) ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;
      
   //-- Attributes

      double CUMUL_TIME ;
      double CURRENT_TIME ;
      bool RUNNING ;
      double CUMUL_ELAPSED_TIME ;
      double CURRENT_ELAPSED_TIME ;
      
} ;


#endif

