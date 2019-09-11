#include <MAC_DoubleComparatorFloat.hh>

#include <MAC_assertions.hh>
#include <MAC.hh>
#include <MAC_ModuleExplorer.hh>

MAC_DoubleComparator const*
MAC_DoubleComparatorFloat:: PROTOTYPE = new MAC_DoubleComparatorFloat() ;

//----------------------------------------------------------------------
MAC_DoubleComparator const*
MAC_DoubleComparatorFloat:: create( MAC_Object* a_owner,
                                    double a_dbl_min )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DoubleComparatorFloat:: create" ) ;
   MAC_CHECK_PRE( a_dbl_min >= 0. ) ;

   MAC_DoubleComparator const* result =
                  new MAC_DoubleComparatorFloat( a_owner, a_dbl_min ) ;
   
   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_DoubleComparator const*
MAC_DoubleComparatorFloat:: create_replica(
              MAC_Object* a_owner, MAC_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DoubleComparatorFloat:: create_replica" ) ;
   MAC_CHECK_PRE( create_replica_PRE( a_owner, exp ) ) ;

   double a_dbl_min = exp->double_data( "dbl_minimum" ) ;
   exp->test_data( "dbl_minimum", "dbl_minimum>=0." ) ;
   
   MAC_DoubleComparator const* result =
                   new MAC_DoubleComparatorFloat( a_owner, a_dbl_min ) ;
   
   MAC_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_DoubleComparatorFloat:: MAC_DoubleComparatorFloat(
                                MAC_Object* a_owner, double a_dbl_min  )
//----------------------------------------------------------------------
   : MAC_DoubleComparator( a_owner )
   , EPSILON( a_dbl_min )
{
}

//----------------------------------------------------------------------
MAC_DoubleComparatorFloat:: MAC_DoubleComparatorFloat( void )
//----------------------------------------------------------------------
   : MAC_DoubleComparator( "MAC_DoubleComparatorFloat" )
   , EPSILON( MAC::bad_double() )
{
}

//----------------------------------------------------------------------
MAC_DoubleComparatorFloat:: ~MAC_DoubleComparatorFloat( void )
//----------------------------------------------------------------------
{
   if( is_a_prototype() ) PROTOTYPE = 0 ;
}

//----------------------------------------------------------------------
int
MAC_DoubleComparatorFloat:: three_way_comparison( double x,
                                                  double y ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DoubleComparatorFloat:: three_way_comparison" ) ;

   int result = 0 ;
   float const xx = ( MAC::abs(x)>EPSILON ? (float) x : 0. ) ;
   float const yy = ( MAC::abs(y)>EPSILON ? (float) y : 0. ) ;
   float diff =  xx - yy ;
   if( diff > 0. )
   {
      result = 1 ;
   }
   else if( diff < 0. )
   {
      result = -1 ;
   }
   return( result ) ;
}

