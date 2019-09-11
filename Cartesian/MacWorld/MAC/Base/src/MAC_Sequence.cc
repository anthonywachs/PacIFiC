#include <MAC_Sequence.hh>

#include <MAC_BalancedBinaryTree.hh>
#include <MAC_Error.hh>
#include <MAC_Iterator.hh>
#include <MAC_assertions.hh>

using std::string ;

//----------------------------------------------------------------------
size_t const MAC_Sequence:: badIndex = static_cast<size_t>(~0) ;
//----------------------------------------------------------------------



//-----------------------------------------------------------------------------
MAC_Sequence:: MAC_Sequence( MAC_Object* a_owner )
//-----------------------------------------------------------------------------
   : MAC_Collection( a_owner )
{
   MAC_LABEL( "MAC_Sequence:: MAC_Sequence" ) ;
}



//-----------------------------------------------------------------------------
MAC_Sequence:: ~MAC_Sequence( void )
//-----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Sequence:: ~MAC_Sequence" ) ;
}



//-------------------------------------------------------------------------
bool
MAC_Sequence:: comparable( MAC_Object const* other ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Sequence:: comparable" ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   // these preconditions are less restrictive than MAC_Object::
   return MAC_Object::comparable( other ) ||
      dynamic_cast<MAC_Sequence const*>( other ) != 0  ;
}


//-------------------------------------------------------------------------
bool
MAC_Sequence:: is_equal( MAC_Object const* other ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Sequence:: is_equal" ) ;
   
   // these preconditions are less restrictive than MAC_Object::is_equal_PRE()
   MAC_CHECK_PRE( is_equal_PRE( other ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   bool result = ( three_way_comparison( other ) == 0 ) ;

   MAC_CHECK_POST( is_equal_POST( result, other ) ) ;
   return( result ) ;
}



//-------------------------------------------------------------------------
int
MAC_Sequence:: three_way_comparison( MAC_Object const* other ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Sequence:: three_way_comparison" ) ;
   
   MAC_CHECK_PRE( three_way_comparison_PRE( other )) ;
   MAC_CHECK_INV( invariant() ) ;

   MAC_Sequence const* otherList = static_cast<MAC_Sequence const*>( other ) ;

   MAC_Iterator* itSelf = create_iterator( 0 ) ;
   MAC_Iterator* itOther = otherList->create_iterator( 0 ) ;
   int result = 0 ;
   
   while( itSelf->is_valid() && itOther->is_valid() )
   {
      int pRet = itSelf->item()->three_way_comparison( itOther->item() ) ;
      if( pRet!=0 )
      {
         result = pRet ;
         break ;
      }
      itSelf->go_next() ;
      itOther->go_next() ;
   }
   
   if( result==0 )
   {
      if( itSelf->is_valid() && !itOther->is_valid() )
      {
         result = 1 ;
      }
      else if( !itSelf->is_valid() && itOther->is_valid() )
      {
         result = -1 ;
      }
   }

   itSelf->destroy() ;
   itOther->destroy() ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( three_way_comparison_POST( result, other ) ) ;
   return( result ) ;
}



//-------------------------------------------------------------------------
size_t
MAC_Sequence:: hash_code( void ) const 
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Sequence:: hash_code" ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   size_t result = 0 ;
   MAC_Iterator* itSelf = create_iterator( 0 ) ;
   for( ; itSelf->is_valid() ; itSelf->go_next() )
   {
      result += itSelf->item()->hash_code() ;
   }
   itSelf->destroy() ;
   return( result ) ;
}



//-------------------------------------------------------------------------
void
MAC_Sequence:: remove( MAC_Object const* object )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Sequence:: remove" ) ;
   
   MAC_CHECK_PRE( remove_PRE( object ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, count, count() ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;

   size_t id = index_of( object ) ;
   MAC_ASSERT( id != badIndex ) ;
   remove_at( id ) ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( remove_POST( OLD( count ), object, OLD( state_id ) ) ) ;
}



//-------------------------------------------------------------------------
void
MAC_Sequence:: sort( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Sequence:: sort" ) ;
   MAC_SAVEOLD( size_t, count, count() ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;

   if( count()>0 )
   {
      MAC_BalancedBinaryTree* bin=MAC_BalancedBinaryTree::create(0) ;
      MAC_Iterator* it=create_iterator(0) ;
      for(it->start();it->is_valid();it->go_next())
      {
         bin->extend( it->item() ) ;
      }
      clear() ;
      it->destroy() ;
      it=bin->create_iterator(0) ;
      for(it->start();it->is_valid();it->go_next())
      {
         append( it->item() ) ;
      }
      bin->destroy() ;
      it->destroy() ;
   }
   
   report_state_change() ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( count()==OLD(count) ) ;
   MAC_CHECK_POST( state_id()!=OLD(state_id) ) ;
   MAC_CHECK_POST(
      FORALL( ( size_t i=1 ; i<count() ; ++i ),
              at(i-1)->three_way_comparison(at(i))<=0 ) ) ;
}



//-------------------------------------------------------------------------
void
MAC_Sequence:: extend( MAC_Object* object )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Sequence:: extend" ) ;
   
   MAC_CHECK_PRE( extend_PRE( object ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( bool, has, has( object ) ) ;
   MAC_SAVEOLD( size_t, count, count() ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;

   if( !has( object ) )
   {
      append( object ) ;
   }
   report_state_change() ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( extend_POST( OLD(has), OLD(count), object, OLD(state_id) ) ) ;
}



//-------------------------------------------------------------------------
void
MAC_Sequence:: destroy_item_and_remove( MAC_Object const* object )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Sequence:: destroy_item_and_remove" ) ;
   
   MAC_CHECK_PRE( object!=0 ) ;
   MAC_CHECK_PRE( has( object ) ) ;
   MAC_CHECK_PRE( object->owner()==0 ) ;
   
   //???????????????????? temps calcul
   MAC_Object* obj = item( object ) ;
   remove( object ) ;
   obj->destroy() ;
}



//---------------------------------------------------------------------
void
MAC_Sequence:: destroy_item_and_remove_at( size_t i )
//---------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Sequence:: destroy_item_and_remove_at" ) ;
   
   MAC_CHECK_PRE( i<index_limit() ) ;
   MAC_CHECK_PRE( at( i )!=0 ) ;
   MAC_CHECK_PRE( at( i )->owner()==0 ) ;
   
   //?????????????????? temps calcul
   MAC_Object* obj = at( i ) ;
   MAC_CHECK( obj!=0 ) ;
   obj->destroy() ;
   remove_at( i ) ;

   MAC_CHECK_POST( at(i)==0 ) ;
}



//----------------------------------------------------------------------
void
MAC_Sequence:: destroy_items_and_clear( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Sequence:: destroy_items_and_clear" ) ;
   
   MAC_CHECK_PRE( destroy_items_and_remove_section_PRE( 0, index_limit() ) ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;

   destroy_items_and_remove_section( 0, index_limit() ) ;
   clear() ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( clear_POST( OLD( state_id ) ) ) ;
}



//----------------------------------------------------------------------
bool
MAC_Sequence:: invariant( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Sequence:: invariant" ) ;
   
   MAC_ASSERT( MAC_Collection::invariant() ) ;
   MAC_ASSERT( count() <= index_limit() ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
MAC_Sequence:: create_clone_POST( MAC_Sequence* result,
                                  MAC_Object* a_owner ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Sequence:: create_clone_POST" ) ;
   
   MAC_ASSERT( MAC_Container::create_clone_POST( result, a_owner ) ) ;
   MAC_ASSERT( result->index_limit() == index_limit() ) ;
   MAC_ASSERT( FORALL( ( size_t i=0 ; i<index_limit() ; ++i ),
                       result->at(i) == at(i) ) ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
MAC_Sequence:: append_PRE( MAC_Object const* object ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Sequence:: append_PRE" ) ;
   
   MAC_ASSERT( object != 0 ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
MAC_Sequence:: append_POST( size_t old_index_limit,
                            size_t old_count,
                            MAC_Object const* object,
                            size_t old_state_id ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Sequence:: append_POST" ) ;

//   MAC_ASSERT( has( object ) ) ;
   MAC_ASSERT( count() == old_count+1 ) ;
   MAC_ASSERT( index_limit() == old_index_limit+1 ) ;
   MAC_ASSERT( at(index_limit()-1) == object ) ;
   MAC_ASSERT( old_state_id != state_id() ) ;
   
   return( true ) ;
}



//----------------------------------------------------------------------
bool
MAC_Sequence:: prepend_PRE( MAC_Object const* object ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Sequence:: prepend_PRE" ) ;
   
   MAC_ASSERT( object!=0 ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
MAC_Sequence:: prepend_POST( size_t old_index_limit,
                             size_t old_count,
                             MAC_Object const* object,
                             size_t old_state_id ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Sequence:: prepend_POST" ) ;

//   MAC_ASSERT( has( object ) ) ;
   MAC_ASSERT( count() == old_count+1 ) ;
   MAC_ASSERT( index_limit() == old_index_limit+1 ) ;
   MAC_ASSERT( at(0) == object ) ;
   MAC_ASSERT( old_state_id != state_id() ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
MAC_Sequence:: set_at_PRE( size_t i,
                           MAC_Object const* object ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Sequence:: set_at_PRE" ) ;
   
   MAC_ASSERT( i < index_limit() ) ;
   MAC_ASSERT( object != 0 ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
MAC_Sequence:: set_at_POST( size_t old_index_limit,
                            size_t old_count,
                            size_t i,
                            MAC_Object const* object,
                            size_t old_state_id  ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Sequence:: set_at_POST" ) ;

//   MAC_ASSERT( has( object ) ) ;
   MAC_ASSERT( old_index_limit == index_limit() ) ;
   MAC_ASSERT( old_count <= count() ) ;
   MAC_ASSERT( at(i) == object ) ;
   MAC_ASSERT( old_state_id != state_id() ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
MAC_Sequence:: insert_at_PRE( size_t i,
                              MAC_Object const* object ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Sequence:: insert_at_PRE" ) ;
   
   MAC_ASSERT( i<index_limit() ) ;
   MAC_ASSERT( object!=0 ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
MAC_Sequence:: insert_at_POST( size_t old_index_limit,
                               size_t old_count,
                               size_t i,
                               MAC_Object const* object,
                               size_t old_state_id  ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Sequence:: insert_at_POST" ) ;

//   MAC_ASSERT( has( object ) ) ;
   MAC_ASSERT( index_limit() == old_index_limit+1 ) ;
   MAC_ASSERT( count() == old_count+1 ) ;
   MAC_ASSERT( at(i) == object ) ;
   MAC_ASSERT( old_state_id != state_id() ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
MAC_Sequence:: remove_at_PRE( size_t i ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Sequence:: remove_at_PRE" ) ;
   
   MAC_ASSERT( i < index_limit() ) ;
   MAC_ASSERT( at(i)!=0 ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
MAC_Sequence:: remove_at_POST( size_t old_index_limit,
                               size_t old_count,
                               size_t old_state_id  ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Sequence:: remove_at_POST" ) ;
   
   MAC_ASSERT( count() <= old_count-1 ) ;
   MAC_ASSERT( index_limit() <= old_index_limit ) ;
   MAC_ASSERT( old_state_id != state_id() ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
MAC_Sequence:: remove_section_PRE( size_t iFirst,
                                   size_t length ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Sequence:: remove_section_PRE" ) ;
   
   MAC_ASSERT( iFirst+length <= index_limit() ) ;
   MAC_ASSERT( length!=0 ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
MAC_Sequence:: remove_section_POST( size_t old_index_limit,
                                    size_t old_count,
                                    size_t length,
                                    size_t old_state_id ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Sequence:: remove_section_POST" ) ;
   
   MAC_ASSERT( index_limit() <= old_index_limit ) ;
   MAC_ASSERT( count() <= old_count ) ;
   MAC_ASSERT( old_state_id != state_id() ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
MAC_Sequence:: at_PRE( size_t i ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Sequence:: at_PRE" ) ;
   
   MAC_ASSERT( i < index_limit() ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
MAC_Sequence:: at_POST( MAC_Object const* resu, size_t i ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Sequence:: at_POST" ) ;
   
   return( true ) ;
}



//----------------------------------------------------------------------
bool
MAC_Sequence:: index_of_PRE( MAC_Object const* object ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Sequence:: index_of_PRE" ) ;
   
   MAC_ASSERT( object != 0 ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
MAC_Sequence:: index_of_POST( size_t result, MAC_Object const* object ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Sequence:: index_of_POST" ) ;

   MAC_ASSERT( EQUIVALENT( !has(object), result==badIndex ) ) ;
   MAC_ASSERT( result==badIndex || ( result < index_limit() ) ) ;

   return( true ) ;
}



//----------------------------------------------------------------------
bool
MAC_Sequence:: destroy_items_and_remove_section_PRE( size_t iFirst,
                                                     size_t length ) const 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Sequence:: destroy_items_and_remove_section_PRE" ) ;
   
   MAC_ASSERT( iFirst+length<=index_limit() ) ;
   MAC_ASSERT( length!=0 ) ;
   MAC_ASSERT( FORALL( ( size_t i=iFirst ; i<iFirst+length ; ++i ),
                       IMPLIES( at(i)!=0, at(i)->owner()==0 ) ) ) ;
   return( true ) ;
}



//----------------------------------------------------------------------
bool
MAC_Sequence:: clear_POST( size_t old_state_id ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Sequence:: clear_POST" ) ;
   MAC_ASSERT( MAC_Collection::clear_POST( old_state_id ) ) ;
   MAC_ASSERT( count()==0 ) ;
   return( true ) ;
}
