#ifndef MAC_DOUBLE_COMPARATOR_EXACT_HH
#define MAC_DOUBLE_COMPARATOR_EXACT_HH

#include <MAC_DoubleComparator.hh>

/*
Server using exact double comparaison.
PUBLISHED
*/

class MAC_DoubleComparatorExact : public MAC_DoubleComparator
{
   public: //-----------------------------------------------------------------

   //-- Instance delivery and initialization

      static MAC_DoubleComparator const* object( void ) ;

   //-- Comparison

      virtual int three_way_comparison( double x, double y ) const ;
      
   protected: //--------------------------------------------------------------
      
   private: //----------------------------------------------------------------
      
      MAC_DoubleComparatorExact( MAC_Object* a_owner ) ;
      
      MAC_DoubleComparatorExact( MAC_DoubleComparatorExact const& other ) ;
      MAC_DoubleComparatorExact& operator=(
                                  MAC_DoubleComparatorExact const& other ) ;

   //-- Plug in

      MAC_DoubleComparatorExact( void ) ;
     ~MAC_DoubleComparatorExact( void ) ;

      virtual MAC_DoubleComparator const* create_replica(
                            MAC_Object* a_owner,
                            MAC_ModuleExplorer const* exp ) const ;

   //-- Class attributes

      static MAC_DoubleComparator const* PROTOTYPE ;
      
} ;

#endif
