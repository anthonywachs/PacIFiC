#ifndef LA_BLOCK_SPARSE_MATRIX_HH
#define LA_BLOCK_SPARSE_MATRIX_HH

#include <LA_SeqMatrix.hh>
#include <doubleVector.hh>
#include <size_t_vector.hh>

class MAC_Vector ;

/*
   Block sparse matrices.
   LA_BlockSeqMatrix uses an array of `LA_SeqMatrix::' to form a super-matrix.
   Each sub-matrix can be accessed separatly.
   Such a matrix can be instanciated in MAC datafiles with following module :

   MODULE LA_Matrix
      concrete_name = "LA_BlockSeqMatrix"
      row_partitioning = `vector of row sizes'
      col_partitioning = `vector of column sizes'
      MODULE block_prototype
        `a valid LA_SeqMatrix'
      END MODULE block_prototyp
   END MODULE LA_Matrix

PUBLISHED
*/

class LA_BlockSeqMatrix : public LA_SeqMatrix
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return an instance for whose each block has the
      // dynamic type of `a_block_proto'.
      static LA_BlockSeqMatrix* create(
                                MAC_Object* a_owner,
                                size_t_vector const& a_row_partitioning,
                                size_t_vector const& a_col_partitioning,
                                LA_SeqMatrix const* a_block_proto ) ;

      virtual LA_BlockSeqMatrix* create_matrix( MAC_Object* a_owner ) const ;

   //-- Characteristics

      virtual size_t allocated_memory( void ) const ;

   //-- Distributed processing

      virtual bool is_desynchronizable( void ) const ;

   //-- Access

      virtual size_t nb_stored_items( void ) const ;

      virtual size_t nb_rows( void ) const ;

      virtual size_t nb_cols( void ) const ;

      virtual double item( size_t i, size_t j ) const ;

      virtual LA_MatrixIterator* create_stored_item_iterator(
                                                MAC_Object* a_owner ) const ;

   //-- Block access(45.)

      void set_block_prototype( LA_SeqMatrix const* a_proto ) ;

      LA_SeqMatrix const* block_prototype( void ) const ;

      // row partitioning (the length of the result is the number of row
      // blocks and each coefficient of the result is the number of rows
      // of its associated blocks)
      size_t_vector const& row_partition( void ) const ;

      // column partitioning (the length of the result is the number of column
      // blocks and each coefficient of the result is the number of columns
      // of its associated blocks)
      size_t_vector const& col_partition( void ) const ;

      // Has `::submatrix'( `ib', `jb' ) stored items ?
      bool has_submatrix( size_t ib, size_t jb ) const ;

      // (`ib',`jb') submatrix
      LA_SeqMatrix const* submatrix( size_t ib, size_t jb ) const ;


      // Assign `x' to the (`i',`j') coefficient of the (`ib',`jb') submatrix.
      void set_submatrix_item( size_t ib, size_t jb, size_t i, size_t j,
                               double x ) ;

      // Add `x' to the (`i',`j') coefficient of the (`ib',`jb') submatrix.
      void add_to_submatrix_item( size_t ib, size_t jb, size_t i, size_t j,
                                  double x ) ;

   //-- Element change

      // Reinitialize `self' by modifying its partitioning and making all
      // its coefficients vanish.
      void set_partitioning_then_nullify(
                                  size_t_vector const& row_partitioning,
                                  size_t_vector const& col_partitioning ) ;

      virtual void nullify_row( size_t i ) ;

      void nullify_item( size_t i, size_t j ) ;

      virtual void set_stored_items( double val ) ;

      virtual void set_item( size_t i, size_t j, double x ) ;

      virtual void add_to_item( size_t i, size_t j, double x ) ;

      virtual void scale( double alpha ) ;

   //-- BLAS level 2 : matrix-vector operators

      virtual void multiply_vec_then_add(
                               LA_Vector const* x, LA_Vector* y,
                               double alpha = 1.0, double beta = 0.0 ) const ;

      virtual void tr_multiply_vec_then_add(
                               LA_Vector const* x, LA_Vector* y,
                               double alpha = 1.0, double beta = 0.0 ) const ;


   //-- BLAS level 3 : matrix-matrix operators

      virtual void set( LA_Matrix const* A, bool same_pattern=false ) ;

      virtual void add_Mat( LA_Matrix const* A,
                            double alpha = 1.0,
                            bool same_pattern=false ) ;

      virtual void add_Mat_Mat( LA_Matrix const* A, LA_Matrix const* B,
                                double alpha = 1.0 ) ;

   //-- Input - Output

      virtual void readMM( std::string const& file ) ;

      virtual void restore( MAC_ModuleExplorer const* exp ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

     ~LA_BlockSeqMatrix( void ) ;
      LA_BlockSeqMatrix( LA_BlockSeqMatrix const& other ) ;
      LA_BlockSeqMatrix& operator=( LA_BlockSeqMatrix const& other ) ;

      LA_BlockSeqMatrix( MAC_Object* a_owner,
                         size_t_vector const& a_row_partitioning,
                         size_t_vector const& a_col_partitioning,
                         LA_SeqMatrix const* a_block_proto ) ;

      virtual void re_initialize_with_global_sizes( size_t a_nb_rows,
                                                    size_t a_nb_cols ) ;

      void init( void ) ;

      LA_SeqMatrix* submat( size_t ib, size_t jb ) const ;

      size_t block_idx( size_t ib, size_t jb ) const ;

   //-- Plug in

      LA_BlockSeqMatrix( void ) ;

      virtual LA_BlockSeqMatrix* create_replica(
                                        MAC_Object* a_owner,
                                        MAC_ModuleExplorer const* exp ) const ;

   //-- BLAS level 3 : matrix-matrix operators

      void set_IMP( LA_BlockSeqMatrix const* A, bool same_pattern ) ;

      void add_Mat_IMP( LA_BlockSeqMatrix const* A, double alpha,
                        bool same_pattern ) ;

      void add_Mat_Mat_IMP( LA_BlockSeqMatrix const* A,
                            LA_BlockSeqMatrix const* B,
                            double alpha ) ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

   //-- Friends

      friend class LA_BlockSeqMatrixIterator ;

   //-- Class attributes

      static LA_BlockSeqMatrix const* PROTOTYPE ;

   //-- Attributes

      size_t NB_ROWS ;
      size_t NB_COLS ;

      // number of elements in each block
      size_t_vector ROW_PARTITION ;
      size_t_vector COL_PARTITION ;

      // index of first element of block i
      size_t_vector ROW_BLOCK_2_FIRST_ELEM ;
      size_t_vector COL_BLOCK_2_FIRST_ELEM ;

      // the block to which belongs the element i
      size_t_vector ROW_ELEM_2_BLOCK ;
      size_t_vector COL_ELEM_2_BLOCK ;

      // dummy vectors:
      LA_SeqVector* const VEC_X ;
      LA_SeqVector* const VEC_Y ;

      LA_SeqMatrix const* BLOCK_PROTO ;

      mutable MAC_Vector* BLOCKS ; // LA_SeqMatrix*
} ;

#endif
