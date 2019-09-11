#ifndef EXT_PETSC_MATRIX_HH
#define EXT_PETSC_MATRIX_HH

#include <LA_Matrix.hh>
#include <EXT_PETScAPI.hh>
#include <EXT_PETScVector.hh>

class MAC_ModuleExplorer ;

class LA_SeqMatrix ;
class intVector ;

/*
   PETSc matrices.

   They can be instantiated in MAC datafiles with following modules :

   For sparse sequential matrix :
   MODULE LA_Matrix
      concrete_name = "PETSc_SeqAIJ"
   END MODULE LA_Matrix

   For sparse, block, symmetric and sequential matrix :
   MODULE LA_Matrix
      concrete_name = "PETSc_SeqSBAIJ"
      block_size = `size of blocks'
   END MODULE LA_Matrix

   For sparse distributed matrix :
   MODULE LA_Matrix
      concrete_name = "PETSc_MPIAIJ"
      d_nz = `max nb non-null items on diagonal sub-matrix'
      o_nz = `max nb non-null items on extra diagonal sub-matrices'
   END MODULE LA_Matrix

   For sparse, block, symmetric and distributed matrix :
   MODULE LA_Matrix
      concrete_name = "PETSc_MPISBAIJ"
      d_nz = `max nb non-null items on diagonal sub-matrix'
      o_nz = `max nb non-null items on extra diagonal sub-matrices'
      block_size = `size of blocks'
   END MODULE LA_Matrix

   For sparse, block, and distributed matrix :
   MODULE LA_Matrix
      concrete_name = "PETSc_MPIBAIJ"
      d_nz = `max nb non-null items on diagonal sub-matrix'
      o_nz = `max nb non-null items on extra diagonal sub-matrices'
      block_size = `size of blocks'
    END MODULE LA_Matrix

PUBLISHED
*/

class EXT_PETScMatrix : public LA_Matrix
{
   public: //-----------------------------------------------------------------

   //-- Instance delivery and initialization

      virtual EXT_PETScMatrix* create_matrix( MAC_Object* a_owner ) const ;

      virtual EXT_PETScVector* create_vector( MAC_Object* a_owner ) const ;

      virtual void re_initialize(
                        size_t a_nb_rows, size_t a_nb_cols,
                        size_t a_nb_local_rows = MAC::bad_index(),
                        size_t a_nb_local_cols = MAC::bad_index() ) ;

   //-- Characteristics

      virtual LA_Implementation const* implementation( void ) const ;

      virtual bool is_symmetric( void ) const ;

   //-- Distributed processing

      virtual bool is_desynchronizable( void ) const ;

      virtual void start_local_modifs( void ) ;

      virtual void stop_local_modifs( void ) ;

      virtual void synchronize( void )  ;

      virtual MAC_DistributedPartition const* row_distribution( void ) const ;
      virtual MAC_DistributedPartition const* col_distribution( void ) const ;

   //-- Access

      virtual LA_SeqMatrix* create_local_matrix( MAC_Object* a_owner ) const ;

      virtual size_t nb_stored_items( void ) const ;

      virtual size_t nb_rows( void ) const ;

      virtual size_t nb_cols( void ) const ;

      virtual void extract_diag( LA_Vector* diag ) const ;

   //-- BLAS level 2 : matrix-vector operators

      virtual void multiply_vec_then_add(
                        LA_Vector const* x, LA_Vector* y,
                        double alpha = 1.0, double beta = 0. ) const  ;

      virtual void tr_multiply_vec_then_add(
                        LA_Vector const* x, LA_Vector* y,
                        double alpha = 1.0, double beta = 0. ) const ;

      virtual void scale_as_diag_mat_mat( LA_Vector const* lvec ) ;

      virtual void scale_as_mat_diag_mat( LA_Vector const* rvec ) ;

      virtual void add_to_diag( LA_Vector const* vec ) ;

   //-- BLAS level 3 : matrix-matrix operators

      // Create a PETSc matrix as a copy of `A'.
      void set( LA_SeqMatrix const* A ) ;

      virtual void set( LA_Matrix const* A, bool same_pattern=false ) ;

      virtual void add_Mat( LA_Matrix const* A, double alpha = 1.0,
                            bool same_pattern=false )  ;

      virtual void add_tMat( LA_Matrix const* A, double alpha = 1.0 ) ;

      virtual void add_Mat_Mat( LA_Matrix const* A, LA_Matrix const* B,
                                double alpha = 1.0 ) ;

   //-- Element change

      virtual void set_item( size_t i, size_t j, double x )  ;

      virtual void add_to_item( size_t i, size_t j, double x )  ;

      virtual void nullify( void )  ;

      virtual void scale( double alpha )  ;

   //-- Input - Output

      virtual void readMM( std::string const& file ) ;

      virtual void print( std::ostream& os, size_t indent_width ) const ;

   //-- Hidden

      // pointer to internal data
      // (for time optimization, should be used with care)
      Mat const& matrix( void ) const ;

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

     ~EXT_PETScMatrix( void ) ;
      EXT_PETScMatrix( void ) ;
      EXT_PETScMatrix( EXT_PETScMatrix const& other ) ;
      EXT_PETScMatrix& operator=( EXT_PETScMatrix const& other ) ;

      EXT_PETScMatrix( MAC_Object* a_owner,
                       MAC_ModuleExplorer const* exp,
                       EXT_PETScMatrix const* other ) ;

      EXT_PETScMatrix( MAC_Object* a_owner,
                       EXT_PETScMatrix const* other ) ;

      bool is_assembled( void ) const ;

      void build( intVector const& nnz ) ;

      void destroy_matrix( void ) ;

    //-- Plug in

      EXT_PETScMatrix( std::string const& a_name,
                       bool sequential,
                       bool symmetric,
                       bool block ) ;


      virtual EXT_PETScMatrix* create_replica(
                                   MAC_Object* a_owner,
                                   MAC_ModuleExplorer const* exp ) const ;

   //-- BLAS level 3 : matrix-matrix operators

      void add_Mat_IMP( Mat A, double alpha = 1.0,
                        bool same_pattern=false )  ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

      virtual bool implementation_POST(
                                LA_Implementation const* result ) const ;

   //-- Class attributes

      static EXT_PETScMatrix const* PROTOTYPE_MPIAIJ ;
      static EXT_PETScMatrix const* PROTOTYPE_MPISBAIJ ;
      static EXT_PETScMatrix const* PROTOTYPE_MPIBAIJ ;
      static EXT_PETScMatrix const* PROTOTYPE_SeqSBAIJ ;
      static EXT_PETScMatrix const* PROTOTYPE_SeqAIJ ;
      static EXT_PETScMatrix const* PROTOTYPE_AIJ ;

   //-- Attributes

      MAC_ModuleExplorer* const EXP ;
      bool const DESTROY_ON_EXIT ;
      bool SYMMETRIC ;
      bool SEQ ;
      bool BLOCK ;
      bool HAS_OPT ;
      bool VERB ;
      Mat MATRIX ;
      size_t NB_ROWS ;
      size_t NB_COLS ;
      MAC_DistributedPartition* ROW_DIST ;
      MAC_DistributedPartition* COL_DIST ;
      bool UPPER ;
} ;

#endif
