#include <MAC_ListIterator.hh>

#include <MAC_List.hh>
#include <MAC_ListItem.hh>
#include <MAC_assertions.hh>


//----------------------------------------------------------------------
MAC_ListIterator*
MAC_ListIterator:: create( MAC_Object* a_owner )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ListIterator:: create" ) ;
   
   MAC_ListIterator* result =  new MAC_ListIterator( a_owner ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( !result->is_valid() ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
MAC_ListIterator:: MAC_ListIterator( MAC_Object* a_owner )
//----------------------------------------------------------------------
   : MAC_Iterator( a_owner )
   , list( 0 )
   , where( 0 )
{
   MAC_LABEL( "MAC_ListIterator:: MAC_ListIterator" ) ;
   MAC_CHECK_INV( invariant() ) ;
}




//----------------------------------------------------------------------
MAC_ListIterator*
MAC_ListIterator:: create( MAC_Object* a_owner, MAC_List const* a_list )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ListIterator:: create" ) ;
   MAC_CHECK_PRE( a_list != 0 ) ;
   
   MAC_ListIterator* result =  new MAC_ListIterator( a_owner, a_list ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( EQUIVALENT( a_list->count()!=0 , result->is_valid() ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
MAC_ListIterator:: MAC_ListIterator( MAC_Object* a_owner,
                                     MAC_List const* aList )
//----------------------------------------------------------------------
   : MAC_Iterator( a_owner, aList )
   , list( aList )
   , where( 0 )
{
   MAC_LABEL( "MAC_ListIterator:: MAC_ListIterator" ) ;
   where = list->theList ;

   MAC_CHECK_INV( invariant() ) ;
}




//----------------------------------------------------------------------
MAC_ListIterator:: ~MAC_ListIterator( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ListIterator:: ~MAC_ListIterator" ) ;
   MAC_CHECK_INV( invariant() ) ;
}




//----------------------------------------------------------------------
void
MAC_ListIterator:: re_initialize( MAC_List const* a_list )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ListIterator:: re_initialize" ) ;
   MAC_CHECK_PRE( a_list != 0 ) ;
   
   base_initialize( a_list ) ;
   
   list = a_list ;
   where = list->theList ;
      
   MAC_CHECK_INV( invariant() ) ;
}




//----------------------------------------------------------------------
bool
MAC_ListIterator:: is_valid( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ListIterator:: is_valid" ) ;
   MAC_CHECK_PRE( is_valid_PRE() ) ;
   MAC_CHECK_INV( invariant() ) ;

   return( where != 0 ) ;
}




//----------------------------------------------------------------------
void
MAC_ListIterator:: start( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ListIterator:: start" ) ;
   MAC_CHECK_INV( invariant() ) ;

   where = list->theList ;
   record_container_state_id() ;
   
   MAC_CHECK_POST( start_POST() ) ;
   MAC_CHECK_INV( invariant() ) ;
}




//----------------------------------------------------------------------
void
MAC_ListIterator:: go_next( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ListIterator:: go_next" ) ;
   MAC_CHECK_PRE( go_next_PRE() ) ;
   
   MAC_CHECK_INV( invariant() ) ;

   where = where->next() ;

   MAC_CHECK_INV( invariant() ) ;
}




//----------------------------------------------------------------------
MAC_Object*
MAC_ListIterator:: item( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ListIterator:: item" ) ;
   MAC_CHECK_PRE( item_PRE() ) ;
   MAC_CHECK_INV( invariant() ) ; 

   MAC_Object* resu = where->val() ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( item_POST( resu ) ) ;

   return( resu ) ;
}




//---------------------------------------------------------------------
bool
MAC_ListIterator:: invariant( void ) const
//---------------------------------------------------------------------
{
   MAC_ASSERT( MAC_Iterator::invariant() ) ;

   MAC_ASSERT( !is_valid() || where != 0 ) ;

   return( true ) ;
}
