#include <MAC_ListIdentity.hh>

#include <iostream>

#include <MAC_assertions.hh>
#include <MAC_ListIterator.hh>

using std::ostream ;
using std::endl ;



//----------------------------------------------------------------------
MAC_ListIdentity*
MAC_ListIdentity:: create( MAC_Object* a_owner )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ListIdentity:: create" ) ;
   MAC_ListIdentity* result = new MAC_ListIdentity( a_owner ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->index_limit() == 0 ) ;
   return( result ) ;
}



//-------------------------------------------------------------------------
MAC_ListIdentity:: MAC_ListIdentity( MAC_Object* a_owner )
//-------------------------------------------------------------------------
   : MAC_List( a_owner )
{
   MAC_LABEL( "MAC_ListIdentity:: MAC_ListIdentity" ) ;
}



//-------------------------------------------------------------------------
MAC_ListIdentity:: ~MAC_ListIdentity( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ListIdentity:: ~MAC_ListIdentity" ) ;
}



//-------------------------------------------------------------------------
MAC_ListIdentity*
MAC_ListIdentity:: create_clone( MAC_Object* a_owner ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ListIdentity:: create_clone" ) ;
   MAC_ListIdentity* result = create( a_owner ) ;
   
   MAC_Iterator* it = create_iterator( 0 ) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      result->append( it->item() ) ;
   }
   it->destroy() ;

   MAC_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}



//-------------------------------------------------------------------------
bool
MAC_ListIdentity:: matching_items( MAC_Object const* object1,
                                   MAC_Object const* object2 ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ListIdentity:: matching_items" ) ;
   MAC_CHECK_PRE( matching_items_PRE( object1, object2 ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   return( object1 == object2 ) ;
}



