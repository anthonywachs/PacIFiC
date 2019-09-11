#include <MAC_BalancedBinaryTreeIterator.hh>

#include <MAC_BalancedBinaryTree.hh>
#include <MAC_BalancedBinaryTreeNode.hh>
#include <MAC_List.hh>
#include <MAC_assertions.hh>

//----------------------------------------------------------------------
MAC_BalancedBinaryTreeIterator*
MAC_BalancedBinaryTreeIterator:: create( MAC_Object* a_owner,
                                         MAC_BalancedBinaryTree const* a_tree )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BalancedBinaryTreeIterator:: create" ) ;
   MAC_BalancedBinaryTreeIterator* result = 
                         new MAC_BalancedBinaryTreeIterator( a_owner, a_tree ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   MAC_CHECK_POST( EQUIVALENT( a_tree->count()!=0 , result->is_valid() ) ) ;
   return( result ) ;
}



//----------------------------------------------------------------------
MAC_BalancedBinaryTreeIterator:: MAC_BalancedBinaryTreeIterator(
                                          MAC_Object* a_owner,
                                          MAC_BalancedBinaryTree const* aTree )
//----------------------------------------------------------------------
   : MAC_Iterator( a_owner, aTree ), 
     tree( aTree ),
     valid( false )
{
   MAC_LABEL( "MAC_BalancedBinaryTreeIterator:: MAC_BalancedBinaryTreeIterator" ) ;
   pathList = MAC_List::create( this ) ;
   reset() ;
   valid = next() ;
}



//----------------------------------------------------------------------
MAC_BalancedBinaryTreeIterator:: ~MAC_BalancedBinaryTreeIterator( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BalancedBinaryTreeIterator:: ~MAC_BalancedBinaryTreeIterator" ) ;
}



//----------------------------------------------------------------------
bool
MAC_BalancedBinaryTreeIterator:: is_valid( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BalancedBinaryTreeIterator:: is_valid" ) ;
   MAC_CHECK_PRE( is_valid_PRE() ) ;
   return( valid ) ;
}



//----------------------------------------------------------------------
void
MAC_BalancedBinaryTreeIterator:: start( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BalancedBinaryTreeIterator:: start" ) ;
   reset() ;
   valid = next() ;
   record_container_state_id() ;
   MAC_CHECK_POST( start_POST() ) ;
}



//----------------------------------------------------------------------
void
MAC_BalancedBinaryTreeIterator:: go_next( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BalancedBinaryTreeIterator:: go_next" ) ;
   MAC_CHECK_PRE( go_next_PRE() ) ;
   valid = next() ;
}



//----------------------------------------------------------------------
bool
MAC_BalancedBinaryTreeIterator:: next( void ) 
//----------------------------------------------------------------------
{
   bool ret = false ;
   if( where )
   {
      MAC_BalancedBinaryTreeNode * less = where->less ;
      MAC_BalancedBinaryTreeNode * more = where->more ;
      
      if( state==LEFT_TO_DO && less!=0 )
      {
         pathList->append( where ) ;
         where = less ;
         ret = next() ;
      }
      else if( ( state==LEFT_TO_DO && less==0 )
               || state == CURRENT_TO_DO )
      {
         ret = true ;
         state = RIGHT_TO_DO ;
      }
      else if( state == RIGHT_TO_DO && more!=0 )
      {
         state = LEFT_TO_DO ;
         pathList->append( where ) ;
         where = more ;
         ret = next() ;
      }
      else  // state == RIGHT_TO_DO && more==0 || state == FINISHED
      {
         if( pathList->count() != 0 )
         {
            MAC_BalancedBinaryTreeNode * up = 
            static_cast<MAC_BalancedBinaryTreeNode * >( pathList->at(pathList->count()-1) ) ;
            pathList->remove_at( pathList->count()-1 ) ;
//???????????????????
//            MAC_BalancedBinaryTreeNode * up =
//               static_cast<MAC_BalancedBinaryTreeNode * >( pathList->remove_at( pathList->count()-1) ) ;
            if( up->more == where )
            {
               where = up ;
               state = FINISHED ;
            }
            else
            {
               MAC_ASSERT( up->less == where ) ;
               where = up ;
               state = CURRENT_TO_DO ;
            }
            ret = next() ;
         }
         else
         {
            where = 0 ;
         }
      }
   }
   return ret ;
}



//----------------------------------------------------------------------
void
MAC_BalancedBinaryTreeIterator:: reset( void )
//----------------------------------------------------------------------
{
   init = true ;
   pathList->clear() ;
   where = tree->top ;
   state = LEFT_TO_DO ;
}


//----------------------------------------------------------------------
MAC_Object*
MAC_BalancedBinaryTreeIterator:: item( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BalancedBinaryTreeIterator:: item" ) ;
   MAC_CHECK_PRE( item_PRE() ) ;

   MAC_Object* result = where->val ;

   MAC_CHECK_POST( item_POST( result ) ) ;
   return( result ) ;
}



