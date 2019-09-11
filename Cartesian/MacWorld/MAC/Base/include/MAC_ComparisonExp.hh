#ifndef MAC_ComparisonExp_HH
#define MAC_ComparisonExp_HH

#include <MAC_Expression.hh>

/* Comparison : min and  max

PUBLISHED
*/

class MAC_ComparisonExp : public MAC_Expression
{
   public: //-------------------------------------------------------

   //-- Type
      
      virtual MAC_Data::Type data_type( void ) const ;
      
   //-- Value
      
      virtual double to_double( MAC_Context const* ct ) const ;
      
   protected: //-------------------------------------------------------
            
   private: //-------------------------------------------------------
      
      MAC_ComparisonExp( void ) ;
     ~MAC_ComparisonExp( void ) ;
      MAC_ComparisonExp( MAC_ComparisonExp const& other ) ;
      MAC_ComparisonExp& operator=( MAC_ComparisonExp const& other ) ;

      enum Function { min_op, max_op } ;

      MAC_ComparisonExp( MAC_Object* a_owner,
			 std::string const& a_name,
			 MAC_Sequence const* argument_list,
			 Function a_op ) ;
   
   //-- Plug in

      MAC_ComparisonExp( std::string const& a_name,
			 Function a_op ) ;

      virtual MAC_ComparisonExp* create_replica(
                                    MAC_Object * a_owner,
                                    MAC_Sequence const* argument_list ) const ;
      
   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments( MAC_Sequence const* some_arguments ) const ;
      
   //-- Class attributes
      
      static MAC_ComparisonExp const* PROTOTYPE_min ;
      static MAC_ComparisonExp const* PROTOTYPE_max ;
      
   //-- Attributes

      Function const OP ;
      MAC_Data const* const FIRST ;
      MAC_Data const* const SECOND ;
} ;

#endif
