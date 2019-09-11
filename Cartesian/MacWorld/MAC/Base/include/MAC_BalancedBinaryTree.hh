#ifndef MAC_BALANCED_BINARY_TREE_HH
#define MAC_BALANCED_BINARY_TREE_HH

#include <MAC_Set.hh>

#include <MAC_BalancedBinaryTreeIterator.hh>

//--------------------------------------------------------------------------
// Sorted sets, implemented as balanced binary search trees, to allow
// fast inserting searching operations
//--------------------------------------------------------------------------
// For each item, three pointers and a size_t are stored.
//--------------------------------------------------------------------------

class MAC_BalancedBinaryTreeNode ;

class MAC_BalancedBinaryTree : public MAC_Set
{
      
   public : //---------------------------------------------------------------

   //-- Instance delivery and initialization

      static MAC_BalancedBinaryTree* create( MAC_Object* a_owner ) ;
      
      virtual MAC_BalancedBinaryTree* create_clone( 
                                             MAC_Object* a_owner ) const ;
      
   //-- Measurement   

      virtual size_t count( void ) const ;

   //-- Access

      virtual MAC_Object* item( MAC_Object const* object ) const ;

      virtual MAC_BalancedBinaryTreeIterator* create_iterator(
                                             MAC_Object* a_owner ) const ;
      
   //-- Element change

      virtual void extend( MAC_Object* object ) ;

      
   //-- Removal
      
      virtual void remove( MAC_Object const* object ) ;

      virtual void destroy_item_and_remove( MAC_Object const* object ) ;

      virtual void clear( void ) ;

      virtual void destroy_items_and_clear( void ) ;
     
   protected ://------------------------------------------------------------
      
      virtual ~MAC_BalancedBinaryTree( void ) ;

      MAC_BalancedBinaryTree( MAC_Object* a_owner ) ;

   //-- Preconditions, Postconditions, Invariant
   
      virtual bool invariant( void ) const ;

   private :  //------------------------------------------------------------

      MAC_BalancedBinaryTree( void ) ;
      MAC_BalancedBinaryTree( MAC_BalancedBinaryTree const& other ) ;
      MAC_BalancedBinaryTree const& operator=( 
                              MAC_BalancedBinaryTree const& other ) ;

      friend class MAC_BalancedBinaryTreeIterator ;
      
   //-- Attributes

      MAC_BalancedBinaryTreeNode* top ;
      size_t nbEntries ;
      
};

#endif
