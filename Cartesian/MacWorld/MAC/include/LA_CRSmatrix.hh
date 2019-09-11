#ifndef LA_CRS_MATRIX_HH
#define LA_CRS_MATRIX_HH

#include <LA_SeqMatrix.hh>
#include <doubleVector.hh>
#include <intVector.hh>
#include <size_t_vector.hh>

class LA_MatrixIterator ;
class LA_CRSmatrixIterator ;
class LA_DiagonalMatrix ;
class MAC_Communicator ;

/*
Sparse matrices with the Compressed Row Storage scheme.

These matrices are designed to be initialized by copying from existing
sparse matrix.

CRS matrix should be optimal for space consumption and product with
vector operator.

Up to now, methods that can extend sparsity pattern can be inhibited by use
of insertion_mode facility (see `::insertion_mode').

Implementation:
   The CRS implementation stores for all element:
      - its value in a table of name `VALUES'
      - its column number in a table of name `COL'
   The values are stored line per line, and the index in `VALUES' and `COL'
   of the first value of each line is stored in a table of name `START'
   (for convenience, the number of elements is stored at the end of `START',
    the index of the first element of the line `i' is then `START'(`i')
    and the last one is `START'(`i+1')-1).

Example:

   matrix:   [ 0 1 2 ]
             [ 0 0 1 ]
             [ 3 0 1 ]

   `VALUES' is [ 1 2 1 3 1 ]
   `COL'    is [ 1 2 2 0 2 ]
   `START'  is [ 0 2 3 5 ]

Hierarchical Data Structure for instantiation :

   MODULE LA_Matrix
      concrete_name = "LA_CRSmatrix"
      [ insertion_mode = <true|false> ]    default is false
   END MODULE LA_Matrix

PUBLISHED
*/

class LA_CRSmatrix : public LA_SeqMatrix
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      // Create a copy of `other' (copying all stored coefficients of `other').
      static LA_CRSmatrix* create( MAC_Object* a_owner,
                                   LA_SeqMatrix const* other ) ;

      virtual LA_CRSmatrix* create_copy( MAC_Object* a_owner,
                                         LA_SeqMatrix const* other ) const ;

      virtual LA_CRSmatrix* create_matrix( MAC_Object* a_owner ) const ;

   //-- Characteristics

      virtual size_t allocated_memory( void ) const ;
      
   //-- Distributed processing

      virtual bool is_desynchronizable( void ) const ;
      
   //-- Access

      virtual LA_MatrixIterator* create_stored_item_iterator(
                                               MAC_Object* a_owner ) const ;

      virtual size_t nb_stored_items( void ) const ;

      virtual size_t nb_rows( void ) const ;

      virtual size_t nb_cols( void ) const ;

      // much slower than accessing through `LA_MatrixIterator::' objects
      virtual double item( size_t i, size_t j ) const ;

      virtual void extract_diag( LA_Vector* diag ) const ;

   //-- Element change

      virtual void nullify_row( size_t i ) ;

      virtual void set_stored_items( double val ) ;

      virtual void set_item( size_t i, size_t j, double x ) ;

      virtual void add_to_item( size_t i, size_t j, double x ) ;

      virtual void scale( double alpha ) ;

      // Can new items be added to existing pattern ?
      bool insertion_mode( void ) const ;

      // Modify insertion mode.
      void set_insertion_mode( bool allowed ) ;

   //-- BLAS level 2 : matrix-vector operators

      virtual void multiply_vec_then_add(
                           LA_Vector const* x, LA_Vector* y,
                           double alpha = 1.0, double beta = 0.0 ) const ;


   //-- BLAS level 3 : matrix-matrix operators

      virtual void set( LA_Matrix const* A, bool same_pattern=false ) ;

   //-- Distributed processing

      void send( MAC_Communicator const* com, size_t dest,
                 bool in_place=false ) const ;

      bool is_sending_in_place( void ) const ;

      void wait_send( MAC_Communicator const* com ) const ;

      static LA_CRSmatrix* receive( MAC_Object* a_owner,
                                    MAC_Communicator const* com,
                                    size_t src ) ;

   //-- System factorization

      virtual void factorize_MILU0( bool modified, double piv_min ) ;

      virtual void solve_LU( LA_SeqVector const* rhs,
                             LA_SeqVector* sol ) const ;

      virtual void relax( double omega,
                          LA_SeqMatrix::relaxation_mode mode,
                          LA_SeqVector const* omega_inv_diag,
                          LA_SeqVector const* rhs,
                          LA_SeqVector* sol ) const ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;

      virtual void readMM( std::string const& file ) ;

      virtual void restore( MAC_ModuleExplorer const* exp ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~LA_CRSmatrix( void ) ;
      LA_CRSmatrix( LA_CRSmatrix const& other ) ;
      LA_CRSmatrix& operator=( LA_CRSmatrix const& other ) ;

      LA_CRSmatrix( MAC_Object* a_owner,
                    size_t a_nb_rows,
                    size_t a_nb_cols) ;

      LA_CRSmatrix( MAC_Object* a_owner,
                    LA_SeqMatrix const* other ) ;

      LA_CRSmatrix( MAC_Object* a_owner,
                    MAC_Communicator const* com,
                    size_t src ) ;

      virtual void re_initialize_with_global_sizes( size_t a_nb_rows,
                                                    size_t a_nb_cols ) ;

   //-- Plug in

      LA_CRSmatrix( void ) ;

      virtual LA_CRSmatrix* create_replica(
                                        MAC_Object* a_owner,
                                        MAC_ModuleExplorer const* exp ) const ;

      virtual void insert( size_t i, size_t j, double x, size_t pos ) ;

   //-- Friends

      friend class LA_CRSmatrixIterator ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

   //-- Class attributes

      static LA_CRSmatrix const* PROTOTYPE ;

   //-- Attributes

      // Matrix dimensions :
      size_t NB_ROWS ;
      size_t NB_COLS ;

      doubleVector VALUES ;
      size_t_vector DIAG ;
      intVector COL ;
      intVector START ;
      size_t NB_ELEMS ;
      bool INSERT ;
      mutable int DIMS[3] ;
      mutable void* REQUEST[4] ;
} ;

#endif
