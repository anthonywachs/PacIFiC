#ifndef MAC_CONTEXT_PAIR_HH
#define MAC_CONTEXT_PAIR_HH

#include <MAC_Context.hh>

class MAC_Data ;
class MAC_Vector ;
class MAC_Variable ;

class MAC_ContextPair : public MAC_Context
{

   public: //-------------------------------------------------------

   //-- Instance delivery and initialization

      // Create, initialize and return an instance (the result is the merge
      // of `first' and `second' context, `second' having priority).
      static MAC_ContextPair* create( MAC_Object* a_owner,
                                      MAC_Context const* first,
                                      MAC_Context const* second) ;
      
      // Reinitialize the internal state, as if the `::create' method was
      // just completed (`self' is the merge of `first' and `second' context,
      // `second' having priority).
      void re_initialize( MAC_Context const* first,
                          MAC_Context const* second ) ;

      // Provide new context pointing on `self' values and upper context.
      virtual MAC_ContextPair* create_clone( MAC_Object* a_owner ) const ;

   //-- Access

      // variable number
      virtual size_t nb_variables( void ) const ;

      // `i'-th variable
      virtual MAC_Variable const* variable( size_t i ) const ;

      // Is `var' included as a variable ?
      virtual bool has_variable( MAC_Variable const* var ) const ;

      // data identified by `var' (possibly null)
      virtual MAC_Data* value( MAC_Variable const* var ) const ;

   protected: //-------------------------------------------------------

   private: //---------------------------------------------------------

      MAC_ContextPair( void ) ;
     ~MAC_ContextPair( void ) ;
      MAC_ContextPair( MAC_ContextPair const& other ) ;
      MAC_ContextPair& operator=( MAC_ContextPair const& other ) ;

      MAC_ContextPair( MAC_Object* a_owner,
                       MAC_Context const* first,
                       MAC_Context const* second ) ;

   //-- Instance delivery and initialization

      virtual void update( void ) ;

      virtual void update_for_destruction_of( MAC_Context const* subject ) ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

   //-- Attributes
      
      MAC_Context const* CT1 ;
      MAC_Context const* CT2 ;
      MAC_Vector* VALUES ;
};

#endif
