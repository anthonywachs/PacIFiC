#ifndef LA_MATRIX_ITERATOR_HH
#define LA_MATRIX_ITERATOR_HH

#include <MAC_Object.hh>

#include <LA.hh>

class LA_Matrix ;

/*
Iterators to traverse a matrix or a single row of a matrix of type
`LA_Matrix::' and visit all the stored coefficients exactly once.

The provided access to coefficients is efficient, sequential
and one-way.

PUBLISHED
*/

class LA_MatrixIterator : public MAC_Object
{

   public : //--------------------------------------------------------------

   //-- Status

      // Is current position valid ?
      virtual bool is_valid( void ) const = 0 ;

   //-- Cursor movement

      // Initialize for a subsequent traversal of all stored coefficients
      // of the attached matrix. Move to first position.
      virtual void start_all_items( void ) = 0 ;

      // Initialize for a subsequent traversal of the stored coefficients
      // of the `i_row'-th row. Move to first-position.
      virtual void start_row_items( size_t i_row ) = 0 ;

      // Move to next stored coefficient.
      virtual void go_next( void ) = 0 ;

   //-- Access

      // matrix the stored coefficients are visited
      virtual LA_Matrix const* matrix( void ) const = 0 ;

      // number of rows of the attached matrix
      virtual size_t nb_rows( void ) const = 0 ;

      // row index ot current position
      virtual size_t row( void ) const = 0 ;

      // column index at current position
      virtual size_t col( void ) const = 0 ;

      // coefficient at current position
      virtual double item( void ) const = 0 ;

   //-- Hidden

      /** Sets the value of the matrix element at the iterator's
          current position. The result is undefined if the iterator is
          no longer valid.
          @param x : the assigned value. */
      virtual void set_item( double x ) const = 0 ;

      /** Adds a scalar to  the value of the matrix element at the
          iterator's current position. The result is undefined if the
          iterator is no longer valid.
          @param x : the added value. */
      virtual void add_to_item( double x ) const = 0 ;

   protected: //------------------------------------------------------------

      LA_MatrixIterator( MAC_Object* a_owner ) ;

      virtual ~LA_MatrixIterator( void ) ;

      void unsynchronized_matrix_state_for_set( void ) const ;
      void unsynchronized_matrix_state_for_add( void ) const ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

      virtual bool start_all_items_PRE( void ) const ;

      virtual bool start_row_items_PRE( size_t i_row ) const ;

      virtual bool go_next_PRE( void ) const ;

      virtual bool matrix_POST( LA_Matrix const* result ) const ;

      virtual bool row_PRE( void ) const ;
      virtual bool row_POST( size_t result ) const ;

      virtual bool col_PRE( void ) const ;
      virtual bool col_POST( size_t result ) const ;

      virtual double item_PRE( void ) const ;

      virtual bool set_item_PRE( void ) const ;
      virtual bool set_item_POST( LA::SyncState old_state ) const ;

      virtual bool add_to_item_PRE( void ) const ;
      virtual bool add_to_item_POST( LA::SyncState old_state ) const ;

   private : //-------------------------------------------------------------

      LA_MatrixIterator( void ) ;
      LA_MatrixIterator( LA_MatrixIterator const& other ) ;
      LA_MatrixIterator& operator=( LA_MatrixIterator const& other ) ;

} ;

#endif
