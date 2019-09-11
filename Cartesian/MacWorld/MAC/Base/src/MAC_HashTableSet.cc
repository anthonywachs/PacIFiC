#include <MAC_HashTableSet.hh>

#include <MAC_Error.hh>
#include <MAC_List.hh>
#include <MAC_Vector.hh>
#include <MAC_HashTableSetIterator.hh>
#include <MAC_assertions.hh>

//-------------------------------------------------------------------------
MAC_HashTableSet*
MAC_HashTableSet:: create( MAC_Object* a_owner, size_t size ) 
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_HashTableSet:: create" ) ;
   MAC_HashTableSet* result =  new MAC_HashTableSet( a_owner, size ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->nb_buckets() == size ) ;
   MAC_CHECK_POST( result->count() == 0 ) ;
   return( result ) ;
}



//-------------------------------------------------------------------------
MAC_HashTableSet:: MAC_HashTableSet( MAC_Object* a_owner, size_t size ) 
//-------------------------------------------------------------------------
   : MAC_Set( a_owner ), 
     the_size( size )
{
   MAC_LABEL( "MAC_HashTableSet:: MAC_HashTableSet" ) ;
   hTable = MAC_Vector::create( this, the_size ) ;
   for( size_t i = 0 ; i<the_size ; i++ )
   {
      hTable->set_at( i, MAC_List::create( this ) ) ;
   }

   MAC_CHECK_INV( invariant() ) ;
}



//-------------------------------------------------------------------------
MAC_HashTableSet:: ~MAC_HashTableSet( void ) 
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_HashTableSet:: ~MAC_HashTableSet" ) ;
   MAC_CHECK_INV( invariant() ) ;
}



//-------------------------------------------------------------------------
MAC_HashTableSet*
MAC_HashTableSet:: create_clone( MAC_Object* a_owner ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_HashTableSet:: create_clone" ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_HashTableSet* result = MAC_HashTableSet::create( a_owner, 1 ) ;
   result->hTable->re_initialize( the_size ) ;
   result->the_size = the_size ;
   for( size_t i = 0 ; i<the_size ; i++ )
   {
      result->hTable->set_at( i, getList(i)->create_clone( result ) ) ;     
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( create_clone_POST( result, a_owner ) ) ;

   return result ; 
}



//-------------------------------------------------------------------------
size_t
MAC_HashTableSet:: nb_buckets( void ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_HashTableSet:: nb_buckets( void )" ) ;
   return( the_size ) ;
}



//-------------------------------------------------------------------------
void
MAC_HashTableSet:: extend( MAC_Object* object )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_HashTableSet:: extend" ) ;
   MAC_CHECK_PRE( extend_PRE( object ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, count, count() ) ;
   MAC_SAVEOLD( bool, has, has(object) ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;
   
   size_t entry = object->hash_code()%the_size ;
   getList(entry)->extend( object ) ;
   report_state_change() ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( extend_POST( OLD(has), OLD(count), object, OLD(state_id) ) ) ;
}



//-------------------------------------------------------------------------
size_t
MAC_HashTableSet:: count( void ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_HashTableSet:: count" ) ;
   MAC_CHECK_INV( invariant() ) ;

   size_t ret = 0 ;
   for( size_t i = 0 ; i<the_size ; i++ )
   {
      ret += getList(i)->count() ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK( count_POST( ret ) ) ;
   return( ret ) ;
}



//-------------------------------------------------------------------------
MAC_Object*
MAC_HashTableSet:: item( MAC_Object const* object ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_HashTableSet:: item" ) ;
   MAC_CHECK_PRE( item_PRE( object ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   size_t entry = object->hash_code()%the_size ;
   MAC_Object* ret = getList(entry)->item( object ) ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( item_POST( ret, object ) ) ;
   return( ret ) ;
}



//-------------------------------------------------------------------------
MAC_HashTableSetIterator*
MAC_HashTableSet:: create_iterator( MAC_Object* a_owner ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_HashTableSet:: create_iterator" ) ;
   MAC_HashTableSetIterator* result = 
                          MAC_HashTableSetIterator::create( a_owner, this ) ;

   MAC_CHECK_POST( create_iterator_POST( result, a_owner ) ) ;
   return( result ) ;
}



//-------------------------------------------------------------------------
void
MAC_HashTableSet:: remove( MAC_Object const* object )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_HashTableSet:: remove" ) ;
   MAC_CHECK_PRE( remove_PRE( object ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, count, count() ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;
   
   size_t entry = object->hash_code()%the_size ;
   getList(entry)->remove( object ) ;
   report_state_change() ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( remove_POST( OLD(count), object, OLD(state_id) ) ) ;
}



//-------------------------------------------------------------------------
void
MAC_HashTableSet:: destroy_item_and_remove( MAC_Object const* object )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_HashTableSet:: destroy_item_and_remove" ) ;
   MAC_Error::object()->raise_not_implemented( this, 
                                               "destroy_item_and_remove" ) ;
}



//-------------------------------------------------------------------------
void
MAC_HashTableSet:: clear( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_HashTableSet:: clear" ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;
   
   for( size_t i = 0 ; i<the_size ; i++ )
   {
       getList(i)->clear() ;
   }
   report_state_change() ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( clear_POST( OLD( state_id ) ) ) ;
}



//-------------------------------------------------------------------------
void
MAC_HashTableSet:: destroy_items_and_clear( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_HashTableSet:: destroy_items_and_clear" ) ;
   MAC_Error::object()->raise_not_implemented( this, 
                                               "destroy_items_and_clear" ) ;
}



//-------------------------------------------------------------------------
MAC_List *
MAC_HashTableSet:: getList( size_t entry )
//-------------------------------------------------------------------------
{
   MAC_List * ret = dynamic_cast<MAC_List*>( hTable->at( entry ) ) ;
   MAC_CHECK( ret!=0 ) ;
   return ret ;
}



//-------------------------------------------------------------------------
MAC_List const*
MAC_HashTableSet:: getList( size_t entry ) const
//-------------------------------------------------------------------------
{
   MAC_List const* ret = dynamic_cast<MAC_List const*>( hTable->at( entry ) ) ;
   MAC_CHECK( ret!=0 ) ;
   return ret ;
}



//-------------------------------------------------------------------------
bool
MAC_HashTableSet:: invariant( void ) const
//-------------------------------------------------------------------------
{
   MAC_ASSERT( the_size>0 ) ;
   size_t cpt = 0 ;
   
   for( size_t i = 0 ; i<the_size ; i++ )
   {
      MAC_List const* pList = getList( i ) ;
      MAC_ASSERT( pList!=0 ) ;
      cpt += pList->count() ;
   }
   MAC_ASSERT( cpt==count() ) ;
   return( true ) ;
} 
