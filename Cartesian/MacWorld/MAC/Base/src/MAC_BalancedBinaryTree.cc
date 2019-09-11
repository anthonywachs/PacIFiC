#include <MAC_BalancedBinaryTree.hh>

#include <MAC_BalancedBinaryTreeNode.hh>
#include <MAC_Error.hh>
#include <MAC_assertions.hh>


//-------------------------------------------------------------------------
MAC_BalancedBinaryTree*
MAC_BalancedBinaryTree:: create( MAC_Object* a_owner ) 
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BalancedBinaryTree:: create" ) ;
   MAC_BalancedBinaryTree* result = new MAC_BalancedBinaryTree( a_owner ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( result->count() == 0 ) ;
   return( result ) ;
}



//-------------------------------------------------------------------------
MAC_BalancedBinaryTree:: MAC_BalancedBinaryTree( MAC_Object* a_owner )
//-------------------------------------------------------------------------
   : MAC_Set( a_owner ), 
     top( 0 ), 
     nbEntries( 0 )
{
   MAC_LABEL( "MAC_BalancedBinaryTree:: MAC_BalancedBinaryTree" ) ;
   MAC_CHECK_INV( invariant() ) ;
}



//-------------------------------------------------------------------------
MAC_BalancedBinaryTree:: ~MAC_BalancedBinaryTree( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BalancedBinaryTree:: ~MAC_BalancedBinaryTree" ) ;
   MAC_CHECK_INV( invariant() ) ;

   if( top!=0 )
   {
      top->destroy() ;
   }
}



//-------------------------------------------------------------------------
MAC_BalancedBinaryTree*
MAC_BalancedBinaryTree:: create_clone( MAC_Object* a_owner ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BalancedBinaryTree:: create_clone" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_BalancedBinaryTree* result =  MAC_BalancedBinaryTree::create( a_owner ) ;
   MAC_Iterator* it = create_iterator( 0 ) ;
   for( ; it->is_valid() ; it->go_next() )
   {         
      result->extend( it->item() ) ;
   }
   it->destroy() ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( create_clone_POST( result, a_owner ) ) ;

   return( result ) ; 
}



//-------------------------------------------------------------------------
void
MAC_BalancedBinaryTree:: extend( MAC_Object* object )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BalancedBinaryTree:: extend" ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;
   MAC_CHECK_PRE( extend_PRE( object ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, count, count() ) ;
   MAC_SAVEOLD( bool, has, has( object ) ) ;

   MAC_Object * ret ;   
   
   if( top==0 )
   {
      top = new MAC_BalancedBinaryTreeNode( object ) ; //?????????????????
      ret = object ;
   }
   else
   {
      ret = top->addChild( object, top ) ;
   }
   if( ret==object )
   {
       nbEntries++ ;
   }
   report_state_change() ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( extend_POST( OLD(has), OLD(count), object, OLD(state_id) ) ) ;
//   return( ret ) ;
}



//-------------------------------------------------------------------------
size_t
MAC_BalancedBinaryTree:: count( void ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BalancedBinaryTree:: count" ) ;
   MAC_CHECK_POST( count_POST( nbEntries ) ) ;
   return( nbEntries ) ;
}



//-------------------------------------------------------------------------
MAC_Object*
MAC_BalancedBinaryTree:: item( MAC_Object const* object ) const 
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BalancedBinaryTree:: item" ) ;
   MAC_CHECK_PRE( item_PRE( object ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   MAC_Object* ret = 0 ;
   if( top!=0 )
   {
      ret = top->search( object ) ;
   }

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( item_POST( ret, object ) ) ;

   return( ret ) ;
}



//-------------------------------------------------------------------------
MAC_BalancedBinaryTreeIterator*
MAC_BalancedBinaryTree:: create_iterator( MAC_Object* a_owner ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BalancedBinaryTree:: create_iterator" ) ;
   MAC_BalancedBinaryTreeIterator* result = 
                     MAC_BalancedBinaryTreeIterator::create( a_owner, this ) ;

   MAC_CHECK_POST( create_iterator_POST( result, a_owner ) ) ;
   return( result ) ;
}



//-------------------------------------------------------------------------
void
MAC_BalancedBinaryTree:: remove( MAC_Object const* object ) 
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BalancedBinaryTree:: remove" ) ;
   MAC_CHECK_PRE( remove_PRE( object ) ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, count, count()) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;
   
   MAC_BalancedBinaryTreeNode* node = top->remove( object, top ) ;
   nbEntries-- ;
   node->destroy() ;
   report_state_change() ;
   
   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( remove_POST( OLD( count ), object, OLD(state_id) ) ) ;
}



//-------------------------------------------------------------------------
void
MAC_BalancedBinaryTree:: destroy_item_and_remove( MAC_Object const* object )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BalancedBinaryTree:: destroy_item_and_remove" ) ;
   MAC_Error::object()->raise_not_implemented( this, 
                                               "destroy_item_and_remove" ) ;
}



//-------------------------------------------------------------------------
void
MAC_BalancedBinaryTree:: clear( void )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BalancedBinaryTree:: clear" ) ;
   MAC_CHECK_INV( invariant() ) ;
   MAC_SAVEOLD( size_t, state_id, state_id() ) ;

   if( top != 0 ) top->destroy() ;
   top = 0 ;
   nbEntries = 0 ;
   report_state_change() ;

   MAC_CHECK_INV( invariant() ) ;
   MAC_CHECK_POST( clear_POST( OLD(state_id) ) ) ;   
}



//-------------------------------------------------------------------------
void
MAC_BalancedBinaryTree:: destroy_items_and_clear( void ) 
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BalancedBinaryTree:: destroy_items_and_clear" ) ;
   MAC_Error::object()->raise_not_implemented( this, 
                                               "destroy_items_and_clear" ) ;
}



//-------------------------------------------------------------------------
bool
MAC_BalancedBinaryTree:: invariant( void ) const
//-------------------------------------------------------------------------
{
   MAC_ASSERT( MAC_Collection::invariant() ) ; //??????????????? MAC_Set
   size_t cpt = 0 ;
   
   if( top!=0 )
   {
      cpt = top->getNbElem() ;
   }
   
   MAC_ASSERT( cpt==nbEntries ) ;
   return( true ) ;
}



