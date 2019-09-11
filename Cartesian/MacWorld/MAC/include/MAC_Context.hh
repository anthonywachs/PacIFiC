#ifndef MAC_CONTEXT_HH
#define MAC_CONTEXT_HH

#include <MAC_Object.hh>

class MAC_Data ;
class MAC_List ;
class MAC_Variable ;

/*
Contexts in which to interpret expressions.

Contexts are mappings from variables to values.
   - variables are non NULL `MAC_Variables::' objects.
   - values that are `MAC_Data::' objects that are possibly NULL (in which
     case the value is said to be cleared).
     
One context can be said to inherit (`::inherit') from another one : it
means that it shares variables from father context with its owns.
Moreover, one context can be extended (`::extend') from another one : it
means that all variables from source one will be copied in new one.
*/

class MAC_Context : public MAC_Object
{

   public: //-------------------------------------------------------

   //-- Instance delivery and initialization

      virtual MAC_Context* create_clone( MAC_Object* a_owner ) const = 0 ;

   //-- Access

      // variable number
      virtual size_t nb_variables( void ) const = 0 ;

      // `i'th variable
      virtual MAC_Variable const* variable( size_t i ) const = 0 ;

      // Is `var' included as a variable ?
      virtual bool has_variable( MAC_Variable const* var ) const = 0 ;

      // data identified by `var'
      virtual MAC_Data* value( MAC_Variable const* var ) const = 0 ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;

      void print( std::ostream& os, size_t indent_width,
                  MAC_Context const* ctx ) const ;
      
  //-- Observer pattern : subject
      
      void attach_observer( MAC_Context* observer ) const ;

      void detach_observer( MAC_Context* observer ) const ;

   protected: //-------------------------------------------------------

      virtual ~MAC_Context( void ) ;

      MAC_Context( MAC_Object* a_owner ) ;

   //-- Observer pattern : subject

      void update_observers( void ) const ;

      void notify_observers_of_my_destruction( void ) const ;

   //-- Instance delivery and initialization

      virtual void update( void ) = 0 ;

      virtual void update_for_destruction_of( 
                                            MAC_Context const* subject ) = 0 ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool variable_PRE( size_t i  ) const ;
      virtual bool variable_POST( MAC_Variable const* result ) const ;
      virtual bool has_variable_PRE( MAC_Variable const* var ) const ;
      virtual bool value_PRE( MAC_Variable const* var ) const ;
      virtual bool value_POST( MAC_Data* result,
                               MAC_Variable const* var ) const ;
      
   private: //---------------------------------------------------------

      MAC_Context( void ) ;
      MAC_Context( MAC_Context const& other ) ;
      MAC_Context& operator=( MAC_Context const& other ) ;

   //-- Attributes

      MAC_List* OBSERVERS ;
};

#endif
