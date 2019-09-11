#ifndef MAC_CONDITIONAL_EXP_HH
#define MAC_CONDITIONAL_EXP_HH

#include <MAC_TransferExp.hh>

/*
Alternative expression :
   ( test1 ? true_alternative1 : test2 ? true_alternative2 : ... false_alternative ).
   If test1 is true, returns true_alternative1 else if test2 is true , returns true_alternative2 else .... else returns false_alternative.

PUBLISHED
*/

class MAC_ConditionalExp : public MAC_TransferExp
{
   public: //----------------------------------------------------------

   //-- Type
      
      virtual MAC_Data::Type data_type( void ) const ;
      
   //-- Formal calculus

      virtual MAC_Data* create_derivative( MAC_Object* a_owner,
                                           MAC_Variable const* var,
                                           MAC_Context const* ct  ) const ;

   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //-------------------------------------------------------
      
   private: //---------------------------------------------------------

     ~MAC_ConditionalExp( void ) ;
      MAC_ConditionalExp( MAC_ConditionalExp const& other ) ;
      MAC_ConditionalExp& operator=( MAC_ConditionalExp const& other ) ;

      MAC_ConditionalExp( MAC_Object* a_owner,
                          MAC_Sequence const* argument_list ) ;

      MAC_Data const* alternative_result( void ) const ;
      
   //-- Plug in
      
      MAC_ConditionalExp( void ) ;

      virtual MAC_ConditionalExp* create_replica( 
                                    MAC_Object * a_owner,
                                    MAC_Sequence const* argument_list ) const ;
  //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments( MAC_Sequence const* some_arguments ) const ;

   //-- Transfer implementation
      
      MAC_Data const* data( MAC_Context const* ct ) const ;      
      
   //-- Class attributes
      
      static MAC_ConditionalExp const* PROTOTYPE ;

   //-- Attributes

      MAC_Data const* const DEFAULT ;
      
} ;

#endif
