#ifndef LA_BLOCK_SPARSE_MATRIX_ITERATOR_HH
#define LA_BLOCK_SPARSE_MATRIX_ITERATOR_HH

#include <LA_MatrixIterator.hh>

class LA_BlockSeqMatrix ;

/*
`LA_MatrixIterator::' objects attached to `LA_BlockSeqMatrix::' objects.

PUBLISHED
*/

class LA_BlockSeqMatrixIterator : public LA_MatrixIterator
{

   public : //--------------------------------------------------------------

   //-- Initialization

      // Create and return an instance attached to `A'.
      static LA_BlockSeqMatrixIterator* create(
                         MAC_Object* a_owner, LA_BlockSeqMatrix const* A ) ;
      
   //-- Status

      virtual bool is_valid( void ) const ;

   //-- Cursor movement

      virtual void start_all_items( void )  ;

      virtual void start_row_items( size_t i_row ) ;
      
      virtual void go_next( void ) ;  

   //-- Access

      virtual LA_Matrix const* matrix( void ) const ;
      
      virtual size_t nb_rows( void ) const ;
      
      virtual size_t row( void ) const ;

      virtual size_t col( void ) const ;

      virtual double item( void ) const ;

   //-- Hidden

      virtual void set_item( double x ) const ;

      virtual void add_to_item( double x ) const ;

   protected: //------------------------------------------------------------

   private : //-------------------------------------------------------------
  
      LA_BlockSeqMatrixIterator( MAC_Object* a_owner,
                                 LA_BlockSeqMatrix const* A ) ;
      
      LA_BlockSeqMatrixIterator( void ) ;
     ~LA_BlockSeqMatrixIterator( void ) ;
      LA_BlockSeqMatrixIterator( LA_BlockSeqMatrixIterator const& other ) ;
      LA_BlockSeqMatrixIterator& operator=( 
                                 LA_BlockSeqMatrixIterator const& other ) ;

      
      bool next( void ) ;

   //-- Attributes
      
      LA_BlockSeqMatrix const* const MAT ;
      LA_MatrixIterator* MAT_IT ;
      bool IS_VALID ;
      size_t IB ;
      size_t JB ;
      bool STAY_IN_ROW ;
      size_t IROW ;
      size_t ROW_SHIFT ;
      size_t COL_SHIFT ;
} ;

#endif
