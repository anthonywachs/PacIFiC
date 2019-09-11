#include <MAC_MapIterator.hh>

#include <MAC_Map.hh>
#include <MAC_KeyItemPair.hh>
#include <MAC_HashTableSet.hh>
#include <MAC_ListIterator.hh>
#include <MAC_assertions.hh>

//----------------------------------------------------------------------
MAC_MapIterator*
MAC_MapIterator:: create( MAC_Object* a_owner,
                          MAC_Map const* map  )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_MapIterator:: create" ) ;

   MAC_MapIterator* result =  new MAC_MapIterator( a_owner, map ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( EQUIVALENT( map->count()!=0 , result->is_valid() ) ) ;
   return( result ) ;
}



//----------------------------------------------------------------------
MAC_MapIterator:: MAC_MapIterator( MAC_Object* a_owner,
                                   MAC_Map const* hTable )
//----------------------------------------------------------------------
   : MAC_Iterator( a_owner, hTable->hList ), 
     table( hTable ), 
     where( 0 )
{
   MAC_LABEL( "MAC_MapIterator:: MAC_MapIterator" ) ;
   where = hTable->hList->create_iterator( this ) ;
}



//----------------------------------------------------------------------
MAC_MapIterator:: ~MAC_MapIterator( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_MapIterator:: ~MAC_MapIterator" ) ;
}



//----------------------------------------------------------------------
bool
MAC_MapIterator:: is_valid( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_MapIterator:: is_valid" ) ;
   MAC_CHECK_PRE( is_valid_PRE() ) ;
   
   return( where->is_valid() ) ;
}



//----------------------------------------------------------------------
void
MAC_MapIterator:: start( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_MapIterator:: start" ) ;
   
   where->start() ;
   record_container_state_id() ;
   
   MAC_CHECK_POST( start_POST() ) ;
}



//----------------------------------------------------------------------
void
MAC_MapIterator:: go_next( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_MapIterator:: go_next" ) ;
   MAC_CHECK_PRE( go_next_PRE() ) ;
   
   where->go_next() ;
}



//----------------------------------------------------------------------
MAC_Object*
MAC_MapIterator:: key( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_MapIterator:: key" ) ;
   MAC_CHECK_PRE( item_PRE() ) ;

   MAC_KeyItemPair* cont = static_cast<MAC_KeyItemPair *>( where->item() ) ;
   MAC_Object* result = cont->key() ;
   
   MAC_CHECK_POST( result != 0 ) ;
   return( result ) ;
}



//----------------------------------------------------------------------
MAC_Object*
MAC_MapIterator:: item( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_MapIterator:: item" ) ;
   MAC_CHECK_PRE( item_PRE() ) ;

   MAC_KeyItemPair* cont = static_cast< MAC_KeyItemPair* >( where->item() ) ;
   MAC_Object* result = cont->item() ;

   MAC_CHECK_POST( item_POST( result ) ) ;   
   return( result ) ;
}



