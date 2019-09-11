#include <MAC_Vector.hh>

#include <iostream>

#include <MAC_System.hh>
#include <MAC_VectorIterator.hh>
#include <MAC_assertions.hh>

//----------------------------------------------------------------------
MAC_Vector*
MAC_Vector:: create( MAC_Object* a_owner, size_t size )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Vector:: create" ) ;
   MAC_Vector* result =  new MAC_Vector( a_owner, size ) ;

   MAC_CHECK_POST( result != 0  ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->index_limit() == size ) ;
   MAC_CHECK_POST( FORALL( ( size_t i=0 ; i<result->index_limit() ; ++i ),
                            result->at(i) == 0 ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_Vector:: MAC_Vector( MAC_Object* a_owner, size_t size )
//----------------------------------------------------------------------
   : MAC_Sequence( a_owner )
   , VECTOR( 0 )
   , LENGTH( 0 )
   , NB_ENTRIES( 0 )
   , CAPACITY( 0 )
{
   MAC_LABEL( "MAC_Vector:: MAC_Vector" ) ;
   re_initialize( size ) ;

   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_Vector:: ~MAC_Vector( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Vector:: ~MAC_Vector" ) ;
   MAC_CHECK_INV( invariant() ) ;

   clear() ;
}

//----------------------------------------------------------------------
MAC_Vector* 
MAC_Vector:: create_clone( MAC_Object* a_owner ) const 
//----------------------------------------------------------------------
//????????? on devrait faire appel au copy constructeur ?????????????
//????????? trop trop lent
{
   MAC_LABEL( "MAC_Vector:: create_clone" ) ;
   MAC_Vector* result = new MAC_Vector( a_owner, 0 ) ;
   result->copy( this ) ;

   MAC_CHECK_POST( create_clone_POST( result, a_owner ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
void
MAC_Vector:: re_initialize( size_t size ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Vector:: re_initialize" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;

   if( LENGTH!=size )
   {
      clear() ;
      LENGTH = size ;
      CAPACITY = size ;
      if( LENGTH>0 )
      {
         VECTOR = new MAC_Object*[ LENGTH ] ;
      }
      else
      {
         VECTOR=0 ;
      }
   }
   nullify() ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( state_id() != OLD( state_id) ) ;
   MAC_CHECK_POST( index_limit() == size ) ;
   MAC_CHECK_POST( count()==0 ) ;
   MAC_CHECK_POST( FORALL( ( size_t i=0 ; i<index_limit() ; ++i ),
                            at(i) == 0 ) ) ;
   
}

//----------------------------------------------------------------------
void 
MAC_Vector:: copy( MAC_Vector const* other ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Vector:: copy" ) ;
   MAC_CHECK_PRE( other != 0 ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;

   re_initialize( other->LENGTH ) ;
   NB_ENTRIES = other->count() ;
   for( size_t i=0 ; i<LENGTH ; i++ )
      VECTOR[i] = other->VECTOR[i] ;
   report_state_change() ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( count() == other->count() ) ;
   MAC_CHECK_POST( index_limit() == other->index_limit() ) ;
   MAC_CHECK_POST( FORALL( ( size_t i=0 ; i<index_limit() ; ++i ),
                           at(i) == other->at(i) ) ) ;
   MAC_CHECK_POST( state_id() != OLD( state_id ) ) ;
}

//----------------------------------------------------------------------
void 
MAC_Vector:: resize( size_t size )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Vector:: resize" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, index_limit, index_limit() ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;
   
   if( size > CAPACITY ) 
   {
      CAPACITY = MAC_System::new_block_size( CAPACITY, size ) ;
      MAC_Object** oldVector = VECTOR ;
      VECTOR = new MAC_Object*[ CAPACITY ] ;
      for( size_t i=LENGTH ; i<CAPACITY ; i++ )
         VECTOR[i] = 0 ;
      if( oldVector != 0 )
      {
         for( size_t i=0 ; i<LENGTH ; i++ )
            VECTOR[i] = oldVector[i] ;
        delete[] oldVector ;
      }
   }
   else if( size < LENGTH )
   {
      for( size_t i=size ; i<LENGTH  ; ++i )
      {
         if( VECTOR[ i ]!=0 )
	 {
	    NB_ENTRIES-- ;
	 }
         VECTOR[ i ] = 0 ;
      }
   }
   LENGTH = size ;
   report_state_change() ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( index_limit() == size ) ;
   MAC_CHECK_POST( 
      FORALL( ( size_t ii=OLD(index_limit) ; ii<index_limit() ; ++ii ),
              at(ii) == 0 ) ) ;
   MAC_CHECK_POST( state_id() != OLD( state_id ) ) ;
}

//----------------------------------------------------------------------
size_t
MAC_Vector:: index_limit( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Vector:: index_limit" ) ;
   MAC_CHECK_INV( invariant() ) ;

   return( LENGTH ) ;
}

//----------------------------------------------------------------------
void
MAC_Vector:: nullify( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Vector:: nullify" ) ;
   MAC_SAVEOLD( size_t, index_limit, index_limit() ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;
   MAC_CHECK_INV( invariant() ) ;

   if( LENGTH != 0 )
   {
      for( size_t i=0 ; i<LENGTH ; i++ )
         VECTOR[i] = 0 ;
   }
   NB_ENTRIES = 0 ;
   report_state_change() ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( state_id() != OLD( state_id ) ) ;
   MAC_CHECK_POST( index_limit() == OLD( index_limit ) ) ;
   MAC_CHECK_POST( count()==0 ) ;
   MAC_CHECK_POST( FORALL( ( size_t i=0 ; i<index_limit() ; ++i ),
                           at(i) == 0 ) ) ;
}

//----------------------------------------------------------------------
void
MAC_Vector:: append( MAC_Object* object )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Vector:: append" ) ;
   MAC_CHECK_PRE( append_PRE( object ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, index_limit, index_limit() ) ;
   MAC_SAVEOLD( size_t, count, count() ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;

   size_t oldSize = LENGTH ;
   resize( oldSize+1 ) ;
   set_at( oldSize, object ) ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( append_POST( OLD( index_limit ),
                                OLD( count ),
                                object,
                                OLD( state_id ) ) ) ;
}

//----------------------------------------------------------------------
void
MAC_Vector:: prepend( MAC_Object* object )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Vector:: prepend" ) ;
   MAC_CHECK_PRE( prepend_PRE( object ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, index_limit, index_limit() ) ;
   MAC_SAVEOLD( size_t, count, count() ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;

   if( LENGTH == 0 )
   {
      append( object ) ;
   }
   else
   {
      insert_at( 0, object ) ;
   }
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( prepend_POST( OLD( index_limit ),
                                 OLD( count ),
                                 object,
                                 OLD( state_id ) ) ) ;
}

//----------------------------------------------------------------------
void
MAC_Vector:: set_at( size_t i, MAC_Object* object )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Vector:: set_at" ) ;
   MAC_CHECK_PRE( ( object==0 && i < index_limit() ) || set_at_PRE( i, object ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, index_limit, index_limit() ) ;
   MAC_SAVEOLD( size_t, count, count() ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;
   
   MAC_Object* ith = VECTOR[ i ] ;
   VECTOR[ i ] = object ;
   if( ith==0 )
   {
      NB_ENTRIES++ ;
   }
   if( object==0 )
   {
      NB_ENTRIES-- ;
   }
   report_state_change() ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( set_at_POST( OLD( index_limit ),
                                OLD( count ),
                                i,
                                object,
                                OLD( state_id ) ) ) ;
}

//----------------------------------------------------------------------
void
MAC_Vector:: insert_at( size_t i, MAC_Object* object )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Vector:: insert_at" ) ;
   MAC_CHECK_PRE( insert_at_PRE( i, object ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, index_limit, index_limit() ) ;
   MAC_SAVEOLD( size_t, count, count() ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;

   size_t size = LENGTH+1 ;
   if( size > CAPACITY )
   {
      CAPACITY = MAC_System::new_block_size( CAPACITY, size ) ;
      MAC_Object** oldVector = VECTOR ;
      VECTOR = new MAC_Object*[ CAPACITY ] ;
      for( size_t j = 0 ; j < i ; ++j )
      {
         VECTOR[ j ] = oldVector[ j ] ;
      }
      for( size_t j = i ; j < LENGTH ; ++j )
      {
         VECTOR[ j+1 ] = oldVector[ j ] ;
      }
      delete [] oldVector ;
   }
   else
   {
      for( size_t j = LENGTH ; j > i ; --j )
      {
         VECTOR[ j ] = VECTOR[ j-1 ] ;
      }
   }
   VECTOR[ i ] = object ;
   NB_ENTRIES++ ;   
   LENGTH++ ;
   report_state_change() ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( insert_at_POST( OLD( index_limit ),
                                   OLD( count ),
                                   i,
                                   object,
                                   OLD( state_id ) ) ) ;
}

//----------------------------------------------------------------------
size_t
MAC_Vector:: count( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Vector:: count" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( count_POST( NB_ENTRIES ) ) ;
   return( NB_ENTRIES ) ;
}

//----------------------------------------------------------------------
MAC_Object*
MAC_Vector:: item( MAC_Object const* object ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Vector:: item" ) ;
   MAC_CHECK_PRE( item_PRE( object ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_Object* result = 0 ;
   size_t i = 0 ;

   for( i = 0 ; i < LENGTH ; i++ )
   {
      if(  VECTOR[ i ]!=0 && VECTOR[ i ]->is_equal( object ) )
      {
         result = VECTOR[ i ] ;
         break ;
      }
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( item_POST( result, object ) ) ;

   return( result ) ;
}

//----------------------------------------------------------------------
MAC_Object*
MAC_Vector:: at( size_t i ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Vector:: at" ) ;
   MAC_CHECK_PRE( at_PRE( i ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_Object* resu = VECTOR[ i ] ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( at_POST( resu, i ) ) ;

   return( resu ) ;
}

//----------------------------------------------------------------------
size_t
MAC_Vector:: index_of( MAC_Object const* object ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Vector:: index_of" ) ;
   MAC_CHECK_PRE( index_of_PRE( object ) ) ;
   MAC_CHECK_INV( invariant() ) ;

   size_t ret = badIndex ;
   
   for( size_t i = 0 ; i < LENGTH ; i++ )
   {
      if( VECTOR[ i ]!=0 && matching_items( VECTOR[ i ], object ) )
      {
         ret = i ;
         break ;
      }
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( index_of_POST( ret, object ) ) ;

   return( ret ) ;
}

//----------------------------------------------------------------------
MAC_VectorIterator*
MAC_Vector:: create_iterator( MAC_Object* a_owner ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Vector:: create_iterator" ) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_VectorIterator* result = MAC_VectorIterator::create( a_owner, this ) ;

   MAC_CHECK_POST( create_iterator_POST( result, a_owner ) ) ;
   return( result ) ;
}

//---------------------------------------------------------------------
void
MAC_Vector:: remove_at( size_t i )
//---------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Vector:: remove_at" ) ;
   MAC_CHECK_PRE( remove_at_PRE( i ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, index_limit, index_limit() ) ;
   MAC_SAVEOLD( size_t, count, count() ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;
   
   if( VECTOR[ i ] != 0 )
   {
      NB_ENTRIES-- ;
   }
   for( size_t j = i ; j < LENGTH - 1 ; ++j )
   {
      VECTOR[ j ] = VECTOR[ j + 1 ] ;
   }
   VECTOR[ LENGTH - 1 ] = 0 ;
  
   report_state_change() ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( remove_at_POST( OLD( index_limit ),
                                   OLD( count ),
                                   OLD( state_id ) ) ) ;
}

//---------------------------------------------------------------------
void
MAC_Vector:: remove_section( size_t iFirst, size_t length )
//---------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Vector:: remove_section" ) ;
   MAC_CHECK_PRE( remove_section_PRE( iFirst, length ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, index_limit, index_limit() ) ;
   MAC_SAVEOLD( size_t, count, count() ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;

   for( size_t j = iFirst ; j < iFirst+length ; j++ )
   {
      if( VECTOR[ j ]!=0 )
      {
	 NB_ENTRIES-- ;
         VECTOR[ j ] = 0 ;
      }
   }
   for( size_t j = iFirst ; j < LENGTH - length ; j++ )
   {
      VECTOR[ j ] = VECTOR[ j + length ] ;
      VECTOR[ j + length ] = 0 ;
   }
   report_state_change() ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( remove_section_POST( OLD( index_limit ), 
                                        OLD( count ), 
                                        length,
                                        OLD( state_id ) ) ) ;
}

//---------------------------------------------------------------------
void
MAC_Vector:: destroy_items_and_remove_section( size_t iFirst, size_t length )
//---------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Vector:: destroy_items_and_remove_section" ) ;
   MAC_CHECK_PRE( destroy_items_and_remove_section_PRE( iFirst, length ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, index_limit, index_limit() ) ;
   MAC_SAVEOLD( size_t, count, count() ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;
   
   for( size_t j = iFirst ; j < iFirst+length ; j++ )
   {
      if( VECTOR[ j ]!=0 )
      {
	 NB_ENTRIES-- ;
         VECTOR[ j ]->destroy() ;
         VECTOR[ j ] = 0 ;
      }
   }
   for( size_t j = iFirst ; j < LENGTH - length ; j++ )
   {
      VECTOR[ j ] = VECTOR[ j + length ] ;
      VECTOR[ j + length ] = 0 ;
   }
   report_state_change() ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( remove_section_POST( OLD( index_limit ), 
                                        OLD( count ), 
                                        length,
                                        OLD( state_id ) ) ) ;
}

//----------------------------------------------------------------------
void
MAC_Vector:: clear( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Vector:: clear" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;

   if( VECTOR!=0 )
   {
      delete [] VECTOR ;
      VECTOR = 0 ;
      LENGTH = 0 ;
      CAPACITY = 0 ;
      NB_ENTRIES = 0 ;
   }
   report_state_change() ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( clear_POST( OLD( state_id ) ) ) ;
}

//----------------------------------------------------------------------
bool
MAC_Vector:: invariant( void ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( MAC_Sequence::invariant() ) ;

   MAC_ASSERT( count() <= index_limit() ) ;
   MAC_ASSERT( LENGTH==0 || VECTOR!=0 ) ;

   return( true ) ;
}

//----------------------------------------------------------------------
bool
MAC_Vector:: remove_at_POST( size_t old_index_limit,
                             size_t old_count,
                             size_t old_state_id ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( MAC_Sequence::remove_at_POST( old_index_limit,
                                             old_count,
                                             old_state_id  ) ) ;

   MAC_ASSERT( index_limit() == old_index_limit ) ;

   return( true ) ;
}

//----------------------------------------------------------------------
bool
MAC_Vector:: remove_section_POST( size_t old_index_limit,
                                  size_t old_count,
                                  size_t length,
                                  size_t old_state_id  ) const
//----------------------------------------------------------------------
{
   MAC_ASSERT( MAC_Sequence::remove_section_POST( old_index_limit, 
                                                  old_count, 
                                                  length,
                                                  old_state_id ) ) ;

   MAC_ASSERT( index_limit() == old_index_limit ) ;

   return( true ) ;
}

