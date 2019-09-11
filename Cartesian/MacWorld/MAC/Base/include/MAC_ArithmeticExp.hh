#ifndef MAC_ARITHMETIC_EXP_HH
#define MAC_ARITHMETIC_EXP_HH

#include <MAC_Expression.hh>

/*

Common algebraic binary operators : {+,-,*,/,modulo}.
Arguments must be integer of real.

PUBLISHED
*/

class MAC_ArithmeticExp : public MAC_Expression
{
   public: //-------------------------------------------------------

   //-- Type
      
      virtual MAC_Data::Type data_type( void ) const ;
      
   //-- Value
      
      virtual double to_double( MAC_Context const* ct ) const ;

      virtual int to_int( MAC_Context const* ct ) const ;

      virtual std::string const& to_string( MAC_Context const* ct ) const ;

   //-- Formal calculus

      virtual MAC_Data* create_derivative( MAC_Object* a_owner,
                                           MAC_Variable const* var,
                                           MAC_Context const* ct ) const ;

   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;
         
   protected: //-----------------------------------------------------
      
   private: //-------------------------------------------------------

      MAC_ArithmeticExp( void ) ;
     ~MAC_ArithmeticExp( void ) ;
      MAC_ArithmeticExp( MAC_ArithmeticExp const& other ) ;
      MAC_ArithmeticExp& operator=( MAC_ArithmeticExp const& other ) ;

      enum AlgebraicOperator { M,L,T,D,MOD } ;
      
      MAC_ArithmeticExp( MAC_Object* a_owner,
                         std::string const& a_name,
                         MAC_Sequence const* argument_list,
                         AlgebraicOperator a_op ) ;

      MAC_Data const* alternative_result( void ) const ;
      
   //-- Plug in

      MAC_ArithmeticExp( std::string const& a_name,
                         AlgebraicOperator a_op ) ;

      virtual MAC_ArithmeticExp* create_replica( 
                                  MAC_Object * a_owner,
                                  MAC_Sequence const* argument_list ) const ;
      
   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments(
                                MAC_Sequence const* some_arguments ) const ;
      
   //-- Formal calculus

      virtual MAC_Data* create_operator_simplification( MAC_Object* a_owner ) ;

   //-- Class attributes
            
      static MAC_ArithmeticExp const* PROTOTYPE_M ;
      static MAC_ArithmeticExp const* PROTOTYPE_L ;
      static MAC_ArithmeticExp const* PROTOTYPE_T ;
      static MAC_ArithmeticExp const* PROTOTYPE_D ;
      static MAC_ArithmeticExp const* PROTOTYPE_MOD ;

   //-- Attributes
      
      AlgebraicOperator const OP ;
      MAC_Data const* const ARG0 ;
      MAC_Data const* const ARG1 ;

      mutable std::string RESULT_STR ;
      
} ;

#endif
