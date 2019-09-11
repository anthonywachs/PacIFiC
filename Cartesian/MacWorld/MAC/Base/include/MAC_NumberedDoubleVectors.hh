#ifndef MAC_NUMBERED_DOUBLE_VECTOR_S_HH
#define MAC_NUMBERED_DOUBLE_VECTOR_S_HH

#include <MAC_Object.hh>

#include <doubleVector.hh>
#include <doubleArray2D.hh>
#include <size_t_vector.hh>

class MAC_BalancedBinaryTree ;
class MAC_BalancedBinaryTreeIterator ;
class MAC_DoubleComparator ;
class MAC_NumberedDoubleVectorsItem ;

/*
Sets of numbered instances of doubleVector, each of them
having the same size.

PUBLISHED
*/

class MAC_NumberedDoubleVectors : public MAC_Object
{
   public: //------------------------------------------------------------------
   
   //-- Instance delivery and initialization

      // Create, initialize and return an instance devoted to store 
      // doubleVector instances of size `a_size_of_items'.
      static MAC_NumberedDoubleVectors* create(
                            MAC_Object* a_owner,
                            MAC_DoubleComparator const* dbl_comp,
                            size_t a_size_of_items ) ;

   //-- Status

      // number items
      size_t nb_items( void ) const ;

      // number of coordinates of all items
      size_t size_of_items( void ) const ;

   //-- Access with indices

      // Does `self' contains an item comparing equal to `item' ?
      bool has( doubleVector const& item ) const ;

      // index of the point matching `item'
      size_t index( doubleVector const& item ) const ;

   //-- Global access
      
      // array containing all items ordered independently from their 
      // introduction 
      doubleArray2D const& ordered_items( void ) const ;    

      // order of items independent from introduction order
      size_t_vector const& order( void ) const ;    
      
   //-- Element change

      // Ensure that self includes an item equal to `item'.
      void extend( doubleVector const& item ) ;
      
   protected: //---------------------------------------------------------------

   private: //-----------------------------------------------------------------

      MAC_NumberedDoubleVectors( void ) ;
     ~MAC_NumberedDoubleVectors( void ) ;
      MAC_NumberedDoubleVectors( MAC_NumberedDoubleVectors const& other ) ;
      MAC_NumberedDoubleVectors& operator=( 
                                 MAC_NumberedDoubleVectors const& other ) ;
      
      MAC_NumberedDoubleVectors( MAC_Object* a_owner, 
                                 MAC_DoubleComparator const* dbl_comp,
                                 size_t a_size_of_items ) ;

   //-- Attributes

      friend class MAC_NumberedDoubleVectorsItem ;
      
      MAC_DoubleComparator const* const DBL_COMP ;
      
      size_t const DIM ;
      size_t NB_PTS ;
      
      MAC_BalancedBinaryTree* PT_TREE ;
      mutable MAC_BalancedBinaryTreeIterator* PT_TREE_IT ;
      mutable doubleArray2D ALL_ITEMS ;
      mutable size_t_vector ORDER ;

      mutable MAC_NumberedDoubleVectorsItem* TMP ;
} ;

#endif
