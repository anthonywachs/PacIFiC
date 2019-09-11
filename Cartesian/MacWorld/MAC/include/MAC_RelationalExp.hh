#ifndef MAC_RELATIONAL_EXP_HH
#define MAC_RELATIONAL_EXP_HH

#include <MAC_Expression.hh>

/* Common binary comparison operators : { <=, >=, <, >, =, != }.
   Arguments must be integer or real.

PUBLISHED
*/

class MAC_RelationalExp : public MAC_Expression
{
   public: //-------------------------------------------------------

   //-- Type
      
      virtual MAC_Data::Type data_type( void ) const ;
      
  //-- Value
      
      virtual bool to_bool( MAC_Context const* ct ) const ;
      
   //-- Input - Output
      
      virtual void print( std::ostream& os, size_t indent_width ) const ;

   protected: //-------------------------------------------------------
      
   private: //-------------------------------------------------------
      
      MAC_RelationalExp( void ) ;
     ~MAC_RelationalExp( void ) ;     
      MAC_RelationalExp( MAC_RelationalExp const& other ) ;
      MAC_RelationalExp& operator=( MAC_RelationalExp const& other ) ;
      
      enum ComparisonOperator { LT, GT, EQ, LE, GE, NEQ } ;
         
      MAC_RelationalExp( MAC_Object* a_owner,
                        std::string const& a_name,
                        MAC_Sequence const* argument_list,
                        ComparisonOperator a_op ) ;
      
   //-- Plug in

      MAC_RelationalExp( std::string const& a_name,
                        ComparisonOperator a_op ) ;

      virtual MAC_RelationalExp* create_replica( 
                                    MAC_Object * a_owner,
                                    MAC_Sequence const* argument_list ) const ;
      
   //-- Characteristics

      virtual std::string const& usage( void ) const ;
      
      virtual bool valid_arguments( MAC_Sequence const* some_arguments ) const ;


   //-- Class attributes
            
      static MAC_RelationalExp const* PROTOTYPE_LE ;
      static MAC_RelationalExp const* PROTOTYPE_GE ;
      static MAC_RelationalExp const* PROTOTYPE_LT ;
      static MAC_RelationalExp const* PROTOTYPE_GT ;
      static MAC_RelationalExp const* PROTOTYPE_EQ ;
      static MAC_RelationalExp const* PROTOTYPE_NEQ ;

   //-- Attributes      

      ComparisonOperator const OP ;
      MAC_Data const* const ARG0 ;
      MAC_Data const* const ARG1 ;
} ;

#endif
