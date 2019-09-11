#include <MAC_HashTableSetIterator.hh>

#include <MAC_HashTableSet.hh>
#include <MAC_ListIterator.hh>
#include <MAC_Vector.hh>
#include <MAC_assertions.hh>

//----------------------------------------------------------------------
MAC_HashTableSetIterator*
MAC_HashTableSetIterator:: create( MAC_Object* a_owner,
                                   MAC_HashTableSet const* a_table  )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_HashTableSetIterator:: create" ) ;

   MAC_HashTableSetIterator* result = 
                               new MAC_HashTableSetIterator( a_owner, a_table ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( EQUIVALENT( a_table->count()!=0 , result->is_valid() ) ) ;
   return( result ) ;
}



//----------------------------------------------------------------------
MAC_HashTableSetIterator:: MAC_HashTableSetIterator( MAC_Object* a_owner,
                                             MAC_HashTableSet const* hTable )
//----------------------------------------------------------------------
   : MAC_Iterator( a_owner, hTable ), 
     table( hTable ), 
     where( 0 ),
     init( true ),
     bucket( 0 )
{
   MAC_LABEL( "MAC_HashTableSetIterator:: MAC_HashTableSetIterator" ) ;
   start() ;
}



//----------------------------------------------------------------------
MAC_HashTableSetIterator:: ~MAC_HashTableSetIterator( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_HashTableSetIterator:: ~MAC_HashTableSetIterator" ) ;
   MAC_CHECK_INV( invariant() ) ;

   if( where != 0 )
   {
      where->destroy() ;
      where = 0 ;
   }
}



//----------------------------------------------------------------------
bool
MAC_HashTableSetIterator:: is_valid( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_HashTableSetIterator:: is_valid" ) ;
   MAC_CHECK_PRE( is_valid_PRE() ) ;
   return( where!=0 ) ;
}



//----------------------------------------------------------------------
void
MAC_HashTableSetIterator:: start( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_HashTableSetIterator:: start" ) ;
   if( where!=0 )
   {
      where->destroy() ;
   }
   where = 0 ;
   
   for( bucket = 0 ; bucket<table->the_size ; bucket++ )
   {
      MAC_List const* list = table->getList( bucket ) ;
      if( list->count() > 0 )
      {
         where = list->create_iterator( 0 ) ;
         MAC_CHECK( where->is_valid() ) ;
         break ;
      }
   }
   record_container_state_id() ;
   MAC_CHECK_POST( start_POST() ) ;
}



//----------------------------------------------------------------------
void
MAC_HashTableSetIterator:: go_next( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_HashTableSetIterator:: go_next" ) ;
   MAC_CHECK_PRE( go_next_PRE() ) ;
   
   if( where!=0 )
   {
      where->go_next() ;
      if( !where->is_valid() )
      {
         where->destroy() ;
         where = 0 ;
         for( bucket++ ; bucket<table->the_size ; bucket++ )
         {
            MAC_List const* list = table->getList( bucket ) ;
            if( list->count() > 0 )
            {
               where = list->create_iterator( 0 ) ;
               MAC_CHECK( where->is_valid() ) ;
               break ;
            }
         }
      }
   }
}



//----------------------------------------------------------------------
MAC_Object*
MAC_HashTableSetIterator:: item( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_HashTableSetIterator:: item" ) ;
   MAC_CHECK_PRE( item_PRE() ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   MAC_Object* resu = where->item() ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( item_POST( resu ) ) ;   

   return( resu ) ;
}



