#ifndef LA_PAIR_OF_MATRIX_ITERATOR_HH
#define LA_PAIR_OF_MATRIX_ITERATOR_HH

#include <LA_MatrixIterator.hh>

/*
`LA_MatrixIterator::' using two `LA_MatrixIterator::'.

PUBLISHED
*/

class LA_PairOfMatrixIterator : public LA_MatrixIterator
{

   public : //--------------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return an instance attached to `internal'.
      static LA_PairOfMatrixIterator* create( MAC_Object* a_owner,
                                              LA_Matrix const* A,
                                              LA_MatrixIterator* first,
                                              LA_MatrixIterator* second ) ;

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
  
      LA_PairOfMatrixIterator( MAC_Object* a_owner,
                               LA_Matrix const* A,
                               LA_MatrixIterator* first,
                               LA_MatrixIterator* second ) ;
      
      LA_PairOfMatrixIterator( void ) ;
     ~LA_PairOfMatrixIterator( void ) ;
      LA_PairOfMatrixIterator( LA_PairOfMatrixIterator const& other ) ;
      LA_PairOfMatrixIterator& operator=( 
                               LA_PairOfMatrixIterator const& other ) ;
    
   //-- Attributes

      LA_Matrix const* const MAT ;
      size_t NB_ROWS ;
      LA_MatrixIterator* FIRST ;
      LA_MatrixIterator* SECOND ;
      LA_MatrixIterator* CURR ;
} ;

#endif
