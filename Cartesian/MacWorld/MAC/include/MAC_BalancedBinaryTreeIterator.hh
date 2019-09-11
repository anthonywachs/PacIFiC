#ifndef MAC_BALANCED_BINARY_TREE_ITERATOR_HH
#define MAC_BALANCED_BINARY_TREE_ITERATOR_HH

#include <MAC_Iterator.hh>

class MAC_BalancedBinaryTree ;
class MAC_BalancedBinaryTreeNode ;
class MAC_List ;

//---------------------------------------------------------------------------
// Iterators on items of MAC_BalancedBinaryTree collections
//---------------------------------------------------------------------------
// Items are accessed from the smallest to the biggest, with respect to
// the total order defined by `three_way_comparison'
//---------------------------------------------------------------------------

class MAC_BalancedBinaryTreeIterator : public MAC_Iterator
{
   public: //---------------------------------------------------------------

   //-- Initialization

      // Create and return an iterator over items of `a_tree'.
      static MAC_BalancedBinaryTreeIterator* create( 
                                       MAC_Object* a_owner,
                                       MAC_BalancedBinaryTree const* a_tree ) ;
       
   //-- Status report

      virtual bool is_valid( void ) const ;

   //-- Access

      virtual MAC_Object* item( void ) const ;

   //-- Cursor movement

      virtual void start( void ) ;

      virtual void go_next( void ) ;
             
   protected: //------------------------------------------------------------

      virtual ~MAC_BalancedBinaryTreeIterator( void ) ;

      MAC_BalancedBinaryTreeIterator( MAC_Object* a_owner,
                                      MAC_BalancedBinaryTree const* aTree ) ;
      
   private: //--------------------------------------------------------------

      MAC_BalancedBinaryTreeIterator( void ) ;
      MAC_BalancedBinaryTreeIterator( 
                               MAC_BalancedBinaryTreeIterator const& other ) ;
      MAC_BalancedBinaryTreeIterator const& operator=( 
                               MAC_BalancedBinaryTreeIterator const& other ) ;

      bool next( void ) ;
      void reset( void ) ;

      //--------------------------------------------------------------------
      // ATTRIBUTES
      //--------------------------------------------------------------------
      MAC_BalancedBinaryTree const* tree ;
      MAC_BalancedBinaryTreeNode* where ;
      enum state_type{ LEFT_TO_DO, CURRENT_TO_DO, RIGHT_TO_DO, FINISHED } ;
      state_type state ;
      MAC_List* pathList ;
      bool init ;
      bool valid ;
} ;

#endif
