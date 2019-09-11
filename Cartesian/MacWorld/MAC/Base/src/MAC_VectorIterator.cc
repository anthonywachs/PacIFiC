#include <MAC_VectorIterator.hh>

#include <MAC_Vector.hh>
#include <MAC_assertions.hh>

//----------------------------------------------------------------------
MAC_VectorIterator*
MAC_VectorIterator:: create( MAC_Object* a_owner,
                             MAC_Vector const* vector )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VectorIterator:: create" ) ;
   MAC_CHECK_PRE( vector!=0 ) ;
   
   MAC_VectorIterator* result =  new MAC_VectorIterator( a_owner, vector ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( EQUIVALENT( vector->count()!=0 , result->is_valid() ) ) ;
   return( result ) ;
}



//----------------------------------------------------------------------
MAC_VectorIterator:: MAC_VectorIterator( MAC_Object* a_owner,
                                         MAC_Vector const* vector )
//----------------------------------------------------------------------
   : MAC_Iterator( a_owner, vector ), 
     vec( vector ),
     valid( false ),
     current_i( static_cast<size_t>(~0) ),
     current_item( 0 )
{
   MAC_LABEL( "MAC_VectorIterator:: MAC_VectorIterator" ) ;
   start() ;
}



//----------------------------------------------------------------------
MAC_VectorIterator:: ~MAC_VectorIterator( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VectorIterator:: ~MAC_VectorIterator" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
void
MAC_VectorIterator:: go_i_th( size_t i )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VectorIterator:: go_i_th" ) ;
   MAC_CHECK_INV( invariant() ) ;

   current_i = i ;
   current_item = vec->at( current_i ) ;
   if( current_item != 0 )
   { 
      valid = true ;
   }
   else
   {
      valid = false ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( IMPLIES( is_valid(), index_of_item()==i ) ) ;
}

//----------------------------------------------------------------------
void
MAC_VectorIterator:: start( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VectorIterator:: start" ) ;
   MAC_CHECK_INV( invariant() ) ;

   valid = false ;
   current_item = 0 ;

   for( current_i=0 ; current_i<vec->index_limit() ; current_i++ )
   {
      current_item = vec->at( current_i ) ;
      if( current_item != 0 )
      {
         valid = true ;
         break ;
      }
   }
   record_container_state_id() ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( start_POST() ) ;
}



//----------------------------------------------------------------------
bool
MAC_VectorIterator:: is_valid( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VectorIterator:: is_valid" ) ;
   MAC_CHECK_PRE( is_valid_PRE() ) ;
   MAC_CHECK_INV( invariant() ) ;

   return( valid ) ;
}



//----------------------------------------------------------------------
void
MAC_VectorIterator:: go_next( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VectorIterator:: go_next" ) ;
   MAC_CHECK_PRE( go_next_PRE() ) ;
   MAC_CHECK_INV( invariant() ) ;

   valid = false ;
   current_i++ ;
   current_item = 0 ;
   for( ; current_i<vec->index_limit() ; current_i++ )
   {
      current_item = vec->at( current_i ) ;
      if( current_item != 0 )
      {
         valid = true ;
         break ;
      }
   }

   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
size_t
MAC_VectorIterator:: index_of_item( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VectorIterator:: index_of_item" ) ;
   MAC_CHECK_PRE( is_valid() ) ;
   MAC_CHECK_INV( invariant() ) ;

   return( current_i ) ;
}

//----------------------------------------------------------------------
MAC_Object*
MAC_VectorIterator:: item( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_VectorIterator:: item" ) ;
   MAC_CHECK_PRE( item_PRE() ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_CHECK_POST( item_POST( current_item ) ) ;

   return( current_item ) ;
}



//----------------------------------------------------------------------
bool
MAC_VectorIterator:: invariant( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( MAC_Iterator::invariant() ) ;

   MAC_ASSERT( IMPLIES( !container_has_been_modified() && is_valid(),
                        ( current_item != 0 &&
                          vec->at( current_i ) == current_item ) ) ) ;

   return( true ) ;
}



