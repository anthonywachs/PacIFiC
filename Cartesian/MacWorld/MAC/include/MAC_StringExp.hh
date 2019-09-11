#ifndef MAC_STRING_EXP_HH
#define MAC_STRING_EXP_HH

#include <MAC_Expression.hh>

/*
Expressions involving strings.

---
name      : empty
arguments : String
type      : Bool

---
name      : to_string
arguments : Int or Double
type      : same as the argument type

___

PUBLISHED
*/

class MAC_StringExp : public MAC_Expression
{
   public: //--------------------------------------------------------

   //-- Type
      
      virtual MAC_Data::Type data_type( void ) const ;
      
  //-- Value
      
      virtual bool to_bool( MAC_Context const* ct = 0 ) const ;

      virtual std::string const& to_string( MAC_Context const* ct = 0 ) const ;
      
   protected: //-----------------------------------------------------
      
   private: //-------------------------------------------------------

      MAC_StringExp( void ) ;
     ~MAC_StringExp( void ) ;
      MAC_StringExp( MAC_StringExp const& other ) ;
      MAC_StringExp& operator=( MAC_StringExp const& other ) ;

      MAC_StringExp( MAC_Object* a_owner,
                     std::string const& a_name,
                     MAC_Sequence const* argument_list ) ;

   //-- Plug in

      MAC_StringExp( std::string const& a_name ) ;

      virtual MAC_StringExp* create_replica( 
                                   MAC_Object * a_owner,
                                   MAC_Sequence const* argument_list ) const ;

   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments( MAC_Sequence const* some_arguments ) const ;
      
   //-- Class attributes
            
      static MAC_StringExp const* PROTOTYPE_EMPTY ;
      static MAC_StringExp const* PROTOTYPE_TO_STRING ;
      
   //-- Attributes      

      mutable std::string STR_RESULT ;
} ;

#endif
