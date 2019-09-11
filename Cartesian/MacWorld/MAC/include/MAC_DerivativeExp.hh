#ifndef MAC_DERIVATIVE_EXP_HH
#define MAC_DERIVATIVE_EXP_HH

#include <MAC_Expression.hh>

/* Differential operators d , dnum.

PUBLISHED
*/

class MAC_DerivativeExp : public MAC_Expression
{
   public: //-------------------------------------------------------
      
   //-- Instance delivery and initialization

      virtual MAC_DerivativeExp* create_clone( MAC_Object* a_owner ) const ;

   //-- Type
      
      virtual MAC_Data::Type data_type( void ) const ;
      
   //-- Value
      
      virtual double to_double( MAC_Context const* ct = 0 ) const ;

      virtual doubleVector const& to_double_vector( 
                                MAC_Context const* ct = 0 ) const ;
      
   //-- Formal calculus

      virtual MAC_Data* create_derivative( MAC_Object* a_owner,
                                           MAC_Variable const* var,
                                           MAC_Context const* ct ) const ;

   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //-------------------------------------------------------
      
   private: //---------------------------------------------------------

      MAC_DerivativeExp( void ) ;
     ~MAC_DerivativeExp( void ) ;
      MAC_DerivativeExp( MAC_DerivativeExp const& other ) ;
      MAC_DerivativeExp& operator=( MAC_DerivativeExp const& other ) ;
      
   //-- Plug in

      enum OP_TYPE { d, dnum } ;
      
      MAC_DerivativeExp( std::string const& a_name, OP_TYPE a_op ) ;

      virtual MAC_DerivativeExp* create_replica(
                                   MAC_Object* a_owner,
                                   MAC_Sequence const* argument_list ) const ;

      MAC_DerivativeExp( MAC_Object* a_owner,
                         std::string const& a_name,
                         OP_TYPE a_op,
                         MAC_Sequence const* argument_list ) ;
      
   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments( MAC_Sequence const* some_arguments ) const ;
      
   //-- Formal calculus

      virtual MAC_Data* create_non_const_simplification( 
                                                 MAC_Object* a_owner ) const ;

      MAC_Data const* derivative( MAC_Context const* ct ) const ;
      
   //-- Class attributes
      
      static MAC_DerivativeExp const* PROTOTYPE_d ;
      static MAC_DerivativeExp const* PROTOTYPE_dnum ;
      
   //-- Attribute

      OP_TYPE const OP ;
      MAC_Data const* const EXP ;
      MAC_Variable const* const VAR ;
      mutable MAC_Data const* DERIVATIVE ;
} ;

#endif
