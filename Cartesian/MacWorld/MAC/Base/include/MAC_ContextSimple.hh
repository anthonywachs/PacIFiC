#ifndef MAC_CONTEXT_SIMPLE_HH
#define MAC_CONTEXT_SIMPLE_HH

#include <MAC_Context.hh>

class MAC_Vector ;
class MAC_Variable ;

class MAC_ContextSimple : public MAC_Context
{

   public: //-------------------------------------------------------

   //-- Instance delivery and initialization

      // Create, initialize and return an instance.
      static MAC_ContextSimple* create( MAC_Object* a_owner ) ;

      // Provide new context pointing on `self' values and upper context.
      virtual MAC_ContextSimple* create_clone( MAC_Object* a_owner ) const ;

   //-- Access

      // variable number
      virtual size_t nb_variables( void ) const ;

      // `i'th variable
      virtual MAC_Variable const* variable( size_t i ) const ;

      // Is `var' included as a variable ?
      virtual bool has_variable( MAC_Variable const* var ) const ;

      // data identified by `var' (possibly null)
      virtual MAC_Data* value( MAC_Variable const* var ) const ;

   //-- Modifier
      
      // Does the evaluation of `a_value' requires `var' ?
      static bool has_circular_definition( MAC_Variable const* var,
                                           MAC_Data const* a_value ) ;
      
      // Ensure that `var' is included as a variable with associated `a_value'.
      void extend( MAC_Variable const* var, MAC_Data const* a_value ) ;

      // Extend `self' with all variable contained in `other'.
      // All value are cloned to `self', upper context is not took in account.
      void extend( MAC_Context const* other ) ;

      // Does `a_value' the value associated to `var'.
      void set_value_of( MAC_Variable const* var,
                         MAC_Data const* a_value ) ;

   protected: //-------------------------------------------------------
      
   private: //---------------------------------------------------------

      MAC_ContextSimple( void ) ;
     ~MAC_ContextSimple( void ) ;
      MAC_ContextSimple( MAC_ContextSimple const& other ) ;
      MAC_ContextSimple& operator=( MAC_ContextSimple const& other ) ;

      MAC_ContextSimple( MAC_Object* a_owner ) ;

   //-- Instance delivery and initialization

      virtual void update( void ) ;

      virtual void update_for_destruction_of( MAC_Context const* subject ) ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

   //-- Attributes

      MAC_Vector* VALUES ;
};

#endif
