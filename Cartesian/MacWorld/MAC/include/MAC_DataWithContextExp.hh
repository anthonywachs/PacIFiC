#ifndef MAC_DATA_WITH_CONTEXT_EXP_HH
#define MAC_DATA_WITH_CONTEXT_EXP_HH

#include <MAC_TransferExp.hh>

/*

Expressions associated to a set of variables.

---
name     : data_with_context
argument : expression, [ variable name (String), variable value ]
type     : first argument type

The returned value is the evaluation of the expression associated to a
context defined with a set of pairs (variable name, variable value).

Example:
  
  MODULE titi
     $DS_X = 3.
     toto = data_with_context( 3.*$DS_X*$DS_ALPHA, "DS_ALPHA", 2. )
  END MODULE titi

PUBLISHED
*/

class MAC_DataWithContextExp : public MAC_TransferExp
{
   public: //-------------------------------------------------------

   //-- Instance delivery and initialization
      
      static MAC_DataWithContextExp* create(
                      MAC_Object* a_owner,
                      MAC_Data const* data, MAC_Context const* ct ) ;
      
   //-- Context
      
      virtual void declare( MAC_List* lst ) const ;
      
      virtual bool context_has_required_variables(
                                      MAC_Context const* ct ) const ;

      
   //-- Type
      
      virtual MAC_Data::Type data_type( void ) const ;
      
   //-- Value
      
      virtual bool value_can_be_evaluated( MAC_Context const* ct ) const ;

      virtual stringVector const& undefined_variables(
                                           MAC_Context const* ct ) const ;
 
   //-- Formal calculus

      virtual bool is_constant( void ) const ;

   protected: //-----------------------------------------------------
      
   private: //-------------------------------------------------------
      
      MAC_DataWithContextExp( void ) ;
     ~MAC_DataWithContextExp( void ) ;
      MAC_DataWithContextExp( MAC_DataWithContextExp const& other ) ;
      MAC_DataWithContextExp& operator=(
                              MAC_DataWithContextExp const& other ) ;
      
      MAC_DataWithContextExp( MAC_Object* a_owner,
                              std::string const& a_name,
                              MAC_Sequence const* argument_list ) ;

   //-- Plug in

      MAC_DataWithContextExp( std::string const& a_name ) ;

      virtual MAC_DataWithContextExp* create_replica( 
                      MAC_Object* a_owner,
                      MAC_Sequence const* argument_list ) const ;

      
   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments( MAC_Sequence const* some_arguments ) const ;
      
   //-- Transfer implementation
      
      virtual MAC_Data const* data( MAC_Context const* ct ) const ;

   //-- Class attributes

      static MAC_DataWithContextExp const* PROTO ;
      
   //-- Attributes
      
      MAC_Data const* DATA ;
} ;

#endif
