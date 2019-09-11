#ifndef MAC_SORT_EXP_HH
#define MAC_SORT_EXP_HH

#include <MAC_Expression.hh>
#include <doubleVector.hh>

/* Operator to sort vector.
   Syntax :
     sort( DoubleVector , option )
   With option a string in : "<", ">".
   Duplicate items are removed from sorted list.

PUBLISHED
*/

class MAC_SortExp : public MAC_Expression
{
   public: //----------------------------------------------------------
      
   //-- Type
      
      virtual MAC_Data::Type data_type( void ) const ;
      
   //-- Value
      
      virtual doubleVector const& to_double_vector( 
                                          MAC_Context const* ct = 0 ) const ;
      
   protected: //-------------------------------------------------------
      
   private: //---------------------------------------------------------
      
     ~MAC_SortExp( void ) ;    
      MAC_SortExp( MAC_SortExp const& other ) ;
      MAC_SortExp& operator=( MAC_SortExp const& other ) ;

      MAC_SortExp( MAC_Object* a_owner,
                      MAC_Sequence const* argument_list ) ;

   //-- Plug in

      MAC_SortExp( void ) ;

      virtual MAC_SortExp* create_replica( 
                                   MAC_Object * a_owner,
                                   MAC_Sequence const* argument_list ) const ;
     
   //-- Characteristics
      
      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments( MAC_Sequence const* some_arguments ) const ;
      
   //-- Class attributes
            
      static MAC_SortExp const* PROTOTYPE ;

   //-- Attributes      

      MAC_Data const* const VECTOR ;
      bool  const GROWING ;
      mutable doubleVector RESULT ;
} ;

#endif
