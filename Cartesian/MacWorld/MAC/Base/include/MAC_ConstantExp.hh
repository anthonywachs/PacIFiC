#ifndef MAC_CONSTANT_EXP_HH
#define MAC_CONSTANT_EXP_HH

#include <MAC_Expression.hh>

/* Constant expressions.
   There are : pi(), e(), euler().

PUBLISHED
*/

class MAC_ConstantExp : public MAC_Expression
{
   public: //-------------------------------------------------------

   //-- Type
      
      virtual MAC_Data::Type data_type( void ) const ;
      
   //-- Value
      
      virtual double to_double( MAC_Context const* ct = 0 ) const ;
      
   //-- Formal calculus

      virtual MAC_Data* create_derivative( MAC_Object* a_owner,
                                           MAC_Variable const* var,
                                           MAC_Context const* ct ) const ;

   protected: //-----------------------------------------------------     
      
   private: //-------------------------------------------------------

      MAC_ConstantExp( void ) ;
     ~MAC_ConstantExp( void ) ;
      MAC_ConstantExp( MAC_ConstantExp const& other ) ;
      MAC_ConstantExp& operator=( MAC_ConstantExp const& other ) ;
      
      MAC_ConstantExp( MAC_Object* a_owner,
                       std::string const& a_name,
                       double a_val,
                       MAC_Sequence const* argument_list ) ;
      
   //-- Plug in

      MAC_ConstantExp( std::string const& a_name, double a_val ) ;
      
      virtual MAC_ConstantExp* create_replica( 
                                 MAC_Object * a_owner,
                                 MAC_Sequence const* argument_list ) const ;
      
   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments( MAC_Sequence const* some_arguments ) const ;
      
   //-- Class attributes
            
      static MAC_ConstantExp const* PROTOTYPE_pi ;
      static MAC_ConstantExp const* PROTOTYPE_e ;
      static MAC_ConstantExp const* PROTOTYPE_eu ;

   //-- Attributes      

      double const VAL ;
} ;

#endif
