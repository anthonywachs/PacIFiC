#ifndef MAC_BALANCED_TREE_NODE_HH
#define MAC_BALANCED_TREE_NODE_HH

#include <MAC_Object.hh>

// Balanced binary tree node class.
// Node realizing a binary tree with balance behaviour.
// This class is used only by MAC_BalancedBinaryTree class.

// A tree node is defined by :
// . a pointer to a MAC_Object
// . the left and right branches such that all elements contained in the
//    left sub-tree are lower than the object of the node and thoses in
//    right sub-tree are larger ( in the sense of three_way_comparison method ).
// . the height (0 for a leap) of the node.
//
// The binary tree is said balanced because it is assumed that left height
// and right one don't differ for more that 2.
// The balance allows then fast inserting and searching.

class MAC_BalancedBinaryTreeNode : public MAC_Object
{
   public: //---------------------------------------------------------------
   protected: //-----------------------------------------------------------

      /* Contructs a leap associated with MAC_Object object. */
      //??????????? mettre le create qui va bien
      MAC_BalancedBinaryTreeNode( MAC_Object * object ) ;
      
      /* Delete self. */ //??????? private.....
      virtual ~MAC_BalancedBinaryTreeNode( void ) ;

      /* Searches in the sub-tree for an item matching value */ 
      MAC_Object * search( const MAC_Object * value ) const ;

      /* Inserts a new item in the sub-tree */
      MAC_Object * addChild( MAC_Object * value,
                             MAC_BalancedBinaryTreeNode *& top ) ;

      /* Delete in the sub-tree the item matching value
         top is the pointer referencing self. */
      MAC_BalancedBinaryTreeNode * remove( const MAC_Object * value,
                                           MAC_BalancedBinaryTreeNode *& top ) ;

      /* Returns cardinal */
      size_t getNbElem( void ) const ;
      
      friend class MAC_BalancedBinaryTree ; //????????? trop de friend
      friend class MAC_BalancedBinaryTreeIterator ;
      
      virtual bool invariant( void ) const ;
     
   private: //--------------------------------------------------------
      
      /* Delete sub-tree */
      void clear( void ) ;
      
      /* Balances sub-tree */ 
      void balance( MAC_BalancedBinaryTreeNode *& top ) ;

      /* Copy items of node in sub-tree */
      void copySubTree( MAC_BalancedBinaryTreeNode * node,
                        MAC_BalancedBinaryTreeNode *& top ) ;

      /* Estimates current height */
      void calculHeight( void ) ;
      
      //--------------------------------------------------------------
      // ATTRIBUTES
      //--------------------------------------------------------------

      size_t height ;
      MAC_Object * val ;
      MAC_BalancedBinaryTreeNode * less ;
      MAC_BalancedBinaryTreeNode * more ;
      
};

#endif // MAC_BALANCED_TREE_NODE_HH
