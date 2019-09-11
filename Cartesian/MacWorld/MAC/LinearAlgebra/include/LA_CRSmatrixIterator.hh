#ifndef LA_CRS_MATRIX_ITERATOR_HH
#define LA_CRS_MATRIX_ITERATOR_HH

#include <LA_MatrixIterator.hh>

class LA_CRSmatrix ;

/*
`LA_MatrixIterator::' objects attached to `LA_CRSmatrix::' objects.

PUBLISHED
*/

class LA_CRSmatrixIterator : public LA_MatrixIterator
{

   public : //--------------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return an instance attached to `A'.
      static LA_CRSmatrixIterator* create( MAC_Object* a_owner, 
                                           LA_CRSmatrix const* A ) ;

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
  
      LA_CRSmatrixIterator( MAC_Object* a_owner, 
                            LA_CRSmatrix const* A ) ;
      
      LA_CRSmatrixIterator( void ) ;
     ~LA_CRSmatrixIterator( void ) ;
      LA_CRSmatrixIterator( LA_CRSmatrixIterator const& other ) ;
      LA_CRSmatrixIterator& operator=( LA_CRSmatrixIterator const& other ) ;
    
   //-- Attributes
    
      LA_CRSmatrix const* const MAT ;

      bool VALID ;
      size_t I_ROW ;
      bool STAY_IN_ROW ;
      int I_CUR ;
} ;

#endif
