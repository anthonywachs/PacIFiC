#ifndef LA_MAC_MATRIX_ITERATOR_HH
#define LA_MAC_MATRIX_ITERATOR_HH

#include <LA_MatrixIterator.hh>

#include <LA_PelMatrix.hh>

/*
`LA_MatrixIterator::' objects attached to `LA_PelMatrix::' objects.

PUBLISHED
*/

class LA_PelMatrixIterator : public LA_MatrixIterator
{

   public : //--------------------------------------------------------------

   //-- Initialization

      // Create and return an instance attached to `A'.
      static LA_PelMatrixIterator* create( MAC_Object* a_owner, 
                                           LA_PelMatrix const* A ) ;

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
  
      LA_PelMatrixIterator( MAC_Object* a_owner, 
                            LA_PelMatrix const* A ) ;
      
      LA_PelMatrixIterator( void ) ;
     ~LA_PelMatrixIterator( void ) ;
      LA_PelMatrixIterator( LA_PelMatrixIterator const& other ) ;
      LA_PelMatrixIterator& operator=( 
                            LA_PelMatrixIterator const& other ) ;

      void go_next_element( void ) ;

   //-- Attributes
      
      LA_PelMatrix const* const MAT ;
      LA_PelMatrix::RowElm*  CUR_ELEM ;
      size_t I_ROW ;
      bool STAY_IN_ROW ;
 
} ;

#endif
