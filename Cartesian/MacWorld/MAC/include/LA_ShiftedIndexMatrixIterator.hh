#ifndef LA_SHIFT_INDEX_MATRIX_ITERATOR_HH
#define LA_SHIFT_INDEX_MATRIX_ITERATOR_HH

#include <LA_MatrixIterator.hh>

/*
`LA_MatrixIterator::' using inter `LA_MatrixIterator::' with shit on column
and row numbers.

PUBLISHED
*/

class LA_ShiftedIndexMatrixIterator : public LA_MatrixIterator
{

   public : //--------------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return an instance attached to `internal'.
      static LA_ShiftedIndexMatrixIterator* create( MAC_Object* a_owner,
                                                    size_t row_shift,
                                                    size_t col_shift,
                                                    LA_Matrix const* A,
                                                    LA_MatrixIterator* internal ) ;

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
  
      LA_ShiftedIndexMatrixIterator( MAC_Object* a_owner,
                                     size_t row_shift,
                                     size_t col_shift,
                                     LA_Matrix const* A,
                                     LA_MatrixIterator* internal
                                     ) ;
      
      LA_ShiftedIndexMatrixIterator( void ) ;
     ~LA_ShiftedIndexMatrixIterator( void ) ;
      LA_ShiftedIndexMatrixIterator(
                            LA_ShiftedIndexMatrixIterator const& other ) ;
      LA_ShiftedIndexMatrixIterator& operator=(
                            LA_ShiftedIndexMatrixIterator const& other ) ;
    
   //-- Attributes

      LA_Matrix const* const MAT ;
      size_t ROW_SHIFT ;
      size_t COL_SHIFT ;
      LA_MatrixIterator* INTERNAL ;
      bool NULLROW ;
} ;

#endif
