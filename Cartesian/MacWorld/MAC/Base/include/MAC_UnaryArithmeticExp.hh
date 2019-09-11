#ifndef MAC_UNARY_ARITHMETIC_EXP_HH
#define MAC_UNARY_ARITHMETIC_EXP_HH

#include <MAC_Expression.hh>

/* Unary arithmetic functions -.
   
PUBLISHED
*/

class MAC_UnaryArithmeticExp : public MAC_Expression
{
   public: //-------------------------------------------------------

   //-- Type
      
      virtual MAC_Data::Type data_type( void ) const ;
      
   //-- Value
      
      virtual double to_double( MAC_Context const* ct = 0 ) const ;

      virtual int to_int( MAC_Context const* ct = 0 ) const ;
      
   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   //-- Formal calculus

      virtual MAC_Data* create_derivative( MAC_Object* a_owner,
                                           MAC_Variable const* var,
                                           MAC_Context const* ct ) const ;

   protected: //-------------------------------------------------------
      
   private: //---------------------------------------------------------

      MAC_UnaryArithmeticExp( void ) ;
     ~MAC_UnaryArithmeticExp( void ) ;
      MAC_UnaryArithmeticExp( MAC_UnaryArithmeticExp const& other ) ;
      MAC_UnaryArithmeticExp& operator=( 
                              MAC_UnaryArithmeticExp const& other ) ;

      enum UnaryOperator { Minus } ;
      
      MAC_UnaryArithmeticExp( MAC_Object* a_owner,
                        std::string const& a_name,
                        MAC_Sequence const* argument_list,
                        UnaryOperator a_op ) ;

      MAC_Data const* alternative_result( void ) const ;
      
   //-- Plug in

      MAC_UnaryArithmeticExp( std::string const& a_name,
                              UnaryOperator a_op ) ;

      virtual MAC_UnaryArithmeticExp* create_replica(
                                   MAC_Object * a_owner,
                                   MAC_Sequence const* argument_list ) const ;
      
   //-- Characteristics

      virtual std::string const& usage( void ) const ;
      
      virtual bool valid_arguments( MAC_Sequence const* some_arguments ) const ;
      
   //-- Class attributes
      
      static MAC_UnaryArithmeticExp const* PROTOTYPE_Minus ;

   //-- Attribute

      UnaryOperator const OP ;
      MAC_Data const* const FIRST ;
};

#endif
