#include <MAC_DoubleComparatorExact.hh>

#include <MAC_assertions.hh>

MAC_DoubleComparator const*
MAC_DoubleComparatorExact:: PROTOTYPE = new MAC_DoubleComparatorExact() ;

//----------------------------------------------------------------------
MAC_DoubleComparator const*
MAC_DoubleComparatorExact:: object( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DoubleComparatorExact:: object" ) ;

   MAC_DoubleComparator const* result = PROTOTYPE ;
   
   MAC_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_DoubleComparator const*
MAC_DoubleComparatorExact:: create_replica(
              MAC_Object* a_owner, MAC_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_DoubleComparatorExact:: create_replica" ) ;
   MAC_CHECK_PRE( create_replica_PRE( a_owner, exp ) ) ;

   MAC_DoubleComparator const* result =
                              new MAC_DoubleComparatorExact( a_owner ) ;
   
   MAC_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_DoubleComparatorExact:: MAC_DoubleComparatorExact( MAC_Object* a_owner )
//----------------------------------------------------------------------
   : MAC_DoubleComparator( a_owner )
{
}

//----------------------------------------------------------------------
MAC_DoubleComparatorExact:: MAC_DoubleComparatorExact( void )
//----------------------------------------------------------------------
   : MAC_DoubleComparator( "MAC_DoubleComparatorExact" )
{
}

//----------------------------------------------------------------------
MAC_DoubleComparatorExact:: ~MAC_DoubleComparatorExact( void )
//----------------------------------------------------------------------
{
   if( is_a_prototype() ) PROTOTYPE = 0 ;
}

//----------------------------------------------------------------------
int
MAC_DoubleComparatorExact:: three_way_comparison( double x,
                                                  double y ) const
//----------------------------------------------------------------------
{
   return( ( x == y ? 0 : ( x < y ? -1 : 1 ) ) ) ;
}

