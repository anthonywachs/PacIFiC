#ifndef MAC_BOOLEAN_EXP_HH
#define MAC_BOOLEAN_EXP_HH

#include <MAC_Expression.hh>

/*
Common boolean operators : {AND,OR,NOT}.
Arguments must be boolean.

PUBLISHED
*/

class MAC_BooleanExp : public MAC_Expression
{
   public: //---------------------------------------------------------

   //-- Type
      
      virtual MAC_Data::Type data_type( void ) const ;
      
   //-- Value
      
      virtual bool to_bool( MAC_Context const* ct ) const ;

   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;

   protected: //-------------------------------------------------------
      
   private: //---------------------------------------------------------

      MAC_BooleanExp( void ) ; 
     ~MAC_BooleanExp( void ) ;
      MAC_BooleanExp( MAC_BooleanExp const& other ) ;
      MAC_BooleanExp& operator=( MAC_BooleanExp const& other ) ;

      enum BoolExp{ OR, AND, NOT } ;

      MAC_BooleanExp( MAC_Object* a_owner,
                      BoolExp exp_id,
                      std::string const& a_name,
                      MAC_Sequence const* argument_list ) ;
      
   //-- Plug in

      MAC_BooleanExp( BoolExp exp_id, std::string const& a_name ) ;
      
      virtual MAC_BooleanExp* create_replica( 
                                   MAC_Object * a_owner,
                                   MAC_Sequence const* argument_list ) const ;
      
   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments( MAC_Sequence const* some_arguments ) const ;
      
   //-- Class attributes

      static MAC_BooleanExp const* PROTOTYPE_or ;
      static MAC_BooleanExp const* PROTOTYPE_and ;
      static MAC_BooleanExp const* PROTOTYPE_not ;

   //-- Attributes

      BoolExp const OP ;
      MAC_Data const* const ARG0 ;
      MAC_Data const* const ARG1 ;
};

#endif
