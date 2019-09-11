#include <MAC_Map.hh>

#include <MAC_Error.hh>
#include <MAC_KeyItemPair.hh>
#include <MAC_MapIterator.hh>
#include <MAC_HashTableSet.hh>
#include <MAC_assertions.hh>

//-------------------------------------------------------------------------
MAC_Map*
MAC_Map:: create( MAC_Object* a_owner, size_t size ) 
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Map:: create" ) ;
   MAC_Map* result = new MAC_Map( a_owner, size ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->nb_buckets() == size ) ;
   MAC_CHECK_POST( result->count() == 0 ) ;
   return( result ) ;
}



//-------------------------------------------------------------------------
MAC_Map:: MAC_Map( MAC_Object* a_owner, size_t aSize ) 
//-------------------------------------------------------------------------
   : MAC_Container( a_owner ), 
     hList( 0 )
{
   MAC_LABEL( "MAC_Map:: MAC_Map" ) ;
   hList = MAC_HashTableSet::create( this, aSize ) ;

   MAC_CHECK_INV( invariant() ) ;
}



//-------------------------------------------------------------------------
MAC_Map:: ~MAC_Map( void ) 
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Map:: ~MAC_Map" ) ;
   MAC_CHECK_INV( invariant() ) ;
}



//-------------------------------------------------------------------------
MAC_Map*
MAC_Map:: create_clone( MAC_Object* a_owner ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Map:: create_clone" ) ;
   return( new MAC_Map( a_owner, this ) ) ;
}



//-------------------------------------------------------------------------
MAC_Map:: MAC_Map( MAC_Object* a_owner, MAC_Map const* other )
//-------------------------------------------------------------------------
   : MAC_Container( a_owner ),
     hList( 0 )
{
   MAC_LABEL( "MAC_Map:: MAC_Map" ) ;
   hList = other->hList->create_clone( this ) ;

   MAC_CHECK_INV( invariant() ) ;
}



//-------------------------------------------------------------------------
size_t
MAC_Map:: nb_buckets( void ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Map:: nb_buckets" ) ;
   return( hList->nb_buckets() ) ;
}



//-------------------------------------------------------------------------
void
MAC_Map:: set_item_at( MAC_Object* key, MAC_Object* a_item  )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Map:: set_item_at" ) ;
   MAC_CHECK_PRE( key != 0 ) ;
   MAC_CHECK_PRE( a_item != 0 ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_Object* object = hList->item( key ) ;
   MAC_KeyItemPair* cont = 0 ;
   MAC_Object* old = 0 ;
   
   if( object != 0 )
   {
      cont = static_cast<MAC_KeyItemPair*>( object ) ;
      old = cont->item() ;
      cont->set_item( a_item ) ;
   }
   else
   {
      cont = MAC_KeyItemPair::create( hList, key, a_item ) ;
      hList->extend( cont ) ;
      old = a_item ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( item_at( key ) == a_item ) ;
}



//-------------------------------------------------------------------------
size_t
MAC_Map:: count( void ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Map:: count" ) ;
   MAC_CHECK_INV( invariant() ) ;

   size_t resu = hList->count() ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( count_POST( resu ) ) ;

   return( resu ) ;
}



//-------------------------------------------------------------------------
bool
MAC_Map:: has_key( MAC_Object const* key ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Map:: has_key" ) ;
   MAC_CHECK_PRE( key != 0 ) ;
   MAC_CHECK_INV( invariant() ) ;

   return( item_at( key ) != 0 ) ;
}



//-------------------------------------------------------------------------
MAC_Object*
MAC_Map:: item( MAC_Object const* object ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Map:: item" ) ;
   MAC_Error::object()->raise_not_implemented( this, "item" ) ;

   return( 0 ) ;
}



//-------------------------------------------------------------------------
MAC_Object* 
MAC_Map:: item_at( MAC_Object const* key ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Map:: item_at" ) ;
   MAC_CHECK_PRE( key!=0 ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_Object* result = 0 ;

   MAC_Object* pair = hList->item( key ) ;
   if( pair != 0 )
   {
      result = static_cast<MAC_KeyItemPair*>( pair )->item() ;
   }
   MAC_CHECK( pair==0 || result!=0 ) ;

   MAC_CHECK_POST( result==0 || has_key( key ) ) ; //?????????
   return( result ) ;
}



//-------------------------------------------------------------------------
MAC_MapIterator*
MAC_Map:: create_iterator( MAC_Object* a_owner ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Map:: create_iterator" ) ;
   return( MAC_MapIterator::create( a_owner, this ) ) ;
}



//-------------------------------------------------------------------------
void
MAC_Map:: remove_at( MAC_Object const* key )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Map:: remove_at" ) ;
   MAC_CHECK_PRE( key!=0 ) ;
   MAC_CHECK_PRE( has_key( key ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   hList->remove( key ) ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( !has_key( key ) ) ; 
}



//-------------------------------------------------------------------------
void
MAC_Map:: clear( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Map:: clear" ) ;
   hList->clear() ;

   MAC_CHECK_POST( count() == 0 ) ;
}



//-------------------------------------------------------------------------
bool
MAC_Map:: invariant( void ) const
//-------------------------------------------------------------------------
{
   MAC_ASSERT( hList!=0 ) ;

   return( true ) ;
} 
