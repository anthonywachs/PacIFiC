#ifndef MAC_DOUBLE_COMPARATOR_FLOAT_HH
#define MAC_DOUBLE_COMPARATOR_FLOAT_HH

#include <MAC_DoubleComparator.hh>

/*
Server using foat comparaison in order to compare two no zero double values,
a constant representing the lower bound under which two double values
are undistinguishable from 0 being given.

PUBLISHED
*/

class MAC_DoubleComparatorFloat : public MAC_DoubleComparator
{
   public: //-----------------------------------------------------------------

   //-- Instance delivery and initialization

      static MAC_DoubleComparator const* create( MAC_Object* a_owner,
                                                 double a_dbl_min ) ;

   //-- Comparison

      virtual int three_way_comparison( double x, double y ) const ;
      
   protected: //--------------------------------------------------------------
      
   private: //----------------------------------------------------------------
      
      MAC_DoubleComparatorFloat( MAC_Object* a_owner, double a_dbl_min ) ;
      
      MAC_DoubleComparatorFloat( MAC_DoubleComparatorFloat const& other ) ;
      MAC_DoubleComparatorFloat& operator=(
                                  MAC_DoubleComparatorFloat const& other ) ;

   //-- Plug in

      MAC_DoubleComparatorFloat( void ) ;
     ~MAC_DoubleComparatorFloat( void ) ;

      virtual MAC_DoubleComparator const* create_replica(
                            MAC_Object* a_owner,
                            MAC_ModuleExplorer const* exp ) const ;

   //-- Class attributes

      static MAC_DoubleComparator const* PROTOTYPE ;
 
   //-- Attributes

      double const EPSILON ;
      
} ;

#endif
