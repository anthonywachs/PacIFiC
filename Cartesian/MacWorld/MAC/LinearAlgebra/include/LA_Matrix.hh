#ifndef LA_MATRIX_HH
#define LA_MATRIX_HH

#include <MAC.hh>
#include <MAC_Object.hh>
#include <MAC_assertions.hh>

#include <LA.hh>

#include <iosfwd>
#include <string>

class LA_Implementation ;
class LA_MatrixIterator ;
class LA_SeqMatrix ;
class LA_Vector ;

class MAC_ModuleExplorer ;
class MAC_ObjectRegister ;
class MAC_DistributedPartition ;

/*
Matrices of double, for linear algebra computations.

PUBLISHED
*/

class LA_Matrix : public MAC_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static LA_Matrix* make( MAC_Object* a_owner,
                              MAC_ModuleExplorer const* exp ) ;

      // Reinitialize `self' by modifying its dimensions and making all
      // its coefficients vanish.
      virtual void re_initialize(
                      size_t a_nb_rows, size_t a_nb_cols,
                      size_t a_nb_local_rows = MAC::bad_index(),
                      size_t a_nb_local_cols = MAC::bad_index() ) = 0 ;

      // Create a matrix with the same implementation as `self'.
      virtual LA_Matrix* create_matrix( MAC_Object* a_owner ) const = 0 ;

      // Create a vector with the same implementation as `self'.
      virtual LA_Vector* create_vector( MAC_Object* a_owner ) const = 0 ;

   //-- Characteristics

      // built-in name
      std::string const& name( void ) const ;

      // implementation indicator
      virtual LA_Implementation const* implementation( void ) const = 0 ;

      // Might the number of rows or the number of columns be changed ?
      bool is_resizable( void ) const ;

      // Is `self' a symmetric matrix ?
      // IMPLEMENTATION : return false
      virtual bool is_symmetric( void ) const ;

   //-- Characteristics setting

      // Fixe the number of rows and the number of columns.
      void make_non_resizable( void ) ;

   //-- Distributed processing(10.)

      // on-process synchronization state
      LA::SyncState state( void ) const ;

      LA::DistributionStrategy distribution_strategy( void ) const ;

      virtual bool is_desynchronizable( void ) const = 0 ;

      // Define environment between call to "start_local_modifs" and
      // "stop_local_modifs" methods
      // in which "set_item", "add_to_item" methods can be used alternatively
      // to modify on-process items without synchronization.
      virtual void start_local_modifs( void ) ;

      virtual void stop_local_modifs( void ) ;

      // Is in environment defined by "start_local_modifs" and
      // "stop_local_modifs methods" ?
      bool only_local_modifs( void ) const ;

      // In distributed case, are non local items of any process to be copied
      // to other process ?
      bool is_synchronized( void ) const ;

      // Ensure that no dispersed items remain in `self'.
      virtual void synchronize( void ) ;

      // distribution of matrix rows over the processes
      virtual MAC_DistributedPartition const* row_distribution( void ) const = 0 ;

      // distribution of matrix columns over the processes
      virtual MAC_DistributedPartition const* col_distribution( void ) const = 0 ;

   //-- Access(50.)

      // Create and return an iterator on stored coefficients :
      //    - in parallel : current process items only ;
      // Default implementation uses `::create_local_matrix'.
      virtual LA_MatrixIterator* create_stored_item_iterator(
                                              MAC_Object* a_owner ) const ;

      // Create and return a COPY of `self' :
      //    - in parallel : current process rows only ;
      //    - to use with care : potentially greedy in time and memory.
      virtual LA_SeqMatrix* create_local_matrix( MAC_Object* a_owner ) const = 0 ;

      // number of stored coefficients
      //    - in parallel : current process items only
      virtual size_t nb_stored_items( void ) const = 0 ;

      // number of rows
      virtual size_t nb_rows( void ) const = 0 ;
      
      // number of rows handled by current process
      size_t nb_local_rows( void ) const ;

      // number of columns
      virtual size_t nb_cols( void ) const = 0 ;
      
      //number of rows handled by current process
      size_t nb_local_cols( void ) const ;

      // Add diagonal of `self' to `diag'.
      // IMPLEMENTATION : raise fatal error.
      virtual void extract_diag( LA_Vector* diag ) const ;

   //-- Element change(60.)

      // Assign `x' to the coefficient at the `i'-th row and `j'-th column.
      virtual void set_item( size_t i, size_t j, double x ) = 0 ;

      // Add `x' to the coefficient at the `i'-th row and `j'-th column.
      virtual void add_to_item( size_t i, size_t j, double x ) = 0 ;

      // Assign 0.0 to all stored coefficients.
      virtual void nullify( void ) = 0 ;

      // Multiply by `alpha' :
      //    `self' <- `alpha'*`self' .
      virtual void scale( double alpha ) = 0 ;

   //-- BLAS level 2 : matrix-vector operators(110.)

      // Multiply by `x' and `alpha' then add to `y' premultiplied by `beta' :
      //    `y' <- ( `alpha'*`self'*`x' + `beta'*`y' ) .
      virtual void multiply_vec_then_add(
                        LA_Vector const* x, LA_Vector* y,
                        double alpha = 1.0, double beta = 0. ) const = 0 ;

      // Multiply the transpose of `self' by `x' and `alpha'
      // then add to `y' premultiplied by `beta' :
      //    `y' <- `alpha'*transpose(`self')*`x'+`beta'*`y'.
      // IMPLEMENTATION : raise fatal error.
      virtual void tr_multiply_vec_then_add(
                        LA_Vector const* x, LA_Vector* y,
                        double alpha = 1.0, double beta = 0. ) const ;

      // Multiply the diagonal Matrix stored in lvec by 'self'  :
      //     `self'[i][j] <- 'lvec'[i] *`self'[i][j] .
      // IMPLEMENTATION : raise fatal error.
      virtual void scale_as_diag_mat_mat( LA_Vector const* lvec ) ;

      // Multiply 'self' by the diagonal Matrix stored in rvec   :
      //     `self'[i][j] <- `self'[i][j] * 'rvec'[j] .
      // IMPLEMENTATION : raise fatal error.
      virtual void scale_as_mat_diag_mat( LA_Vector const* rvec ) ;

      // Add `vec' to diagonal of 'self'.
      // IMPLEMENTATION : raise fatal error.
      virtual void add_to_diag( LA_Vector const* vec ) ;

   //-- BLAS level 3 : matrix-matrix operators(120.)

      // Set `self' equal to : `A'.
      // If `same_pattern' is true
      // then `A' MUST have the same zero pattern than `self'.
      virtual void set( LA_Matrix const* A, bool same_pattern = false ) ;

      // Multiply `A' by `alpha' then add to `self' :
      //     `self' <- `self' + `alpha' * `A' .
      // If `same_pattern' is true
      // then `A' MUST have the same zero pattern than `self'.
      virtual void add_Mat( LA_Matrix const* A,
                            double alpha = 1.0,
                            bool same_pattern = false ) = 0 ;

      // Multiply the transpose of `A' by `alpha' then add to `self' :
      //     `self' <- `self' + `alpha' * transpose(`A') .
      // IMPLEMENTATION : raise fatal error.
      virtual void add_tMat( LA_Matrix const* A, double alpha = 1.0 ) ;

      // Multiply `A' by `B' and `alpha', then add to `self':
      //      `self' <- `self' + `alpha'*`A'*`B'  .
      // IMPLEMENTATION : raise fatal error.
      virtual void add_Mat_Mat( LA_Matrix const* A, LA_Matrix const* B,
                                double alpha = 1.0 ) ;

   //-- Input - Output

      // Reinitialize with matrix stored in file called `file'
      // under the Matrix Market format.
      virtual void readMM( std::string const& file ) ;

      // Write in file called `file' under the Matrix Market format.
      virtual void writeMM( std::string const& file ) const ;

      virtual void print( std::ostream& os, size_t indent_width ) const ;

      // Print stored items of `self'.
      virtual void print_items( std::ostream& os, size_t indent_width ) const ;

   protected: //--------------------------------------------------------

      virtual ~LA_Matrix( void ) ;

      LA_Matrix( MAC_Object* a_owner, 
                 std::string const& a_name,
                 LA::DistributionStrategy a_dist_strat ) ;

   //-- Plug in

      // for prototype registration only
      LA_Matrix( std::string const& a_name ) ;

      // Is `self' a prototype ?
      bool is_a_prototype( void ) const ;

      virtual LA_Matrix* create_replica(
                                   MAC_Object* a_owner,
                                   MAC_ModuleExplorer const* exp ) const = 0 ;

   //-- Distributed processing

      void set_unsynchronized_state( LA::SyncState new_state ) ;

      void set_only_local_modifs_state( bool only_local ) ;

   //-- Input - Output

      bool readMMHeader( std::ifstream& in,
                         bool& array, bool& sym,
                         size_t& n, size_t& m, size_t& nb ) const ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

      virtual bool create_matrix_PRE( MAC_Object* a_owner ) const ;
      virtual bool create_matrix_POST( LA_Matrix* result,
                                       MAC_Object* a_owner ) const ;

      virtual bool create_vector_PRE( MAC_Object* a_owner ) const ;
      virtual bool create_vector_POST( LA_Vector const* result,
                                       MAC_Object* a_owner ) const ;

      virtual bool re_initialize_PRE(
                     size_t a_nb_rows, size_t a_nb_cols,
                     size_t a_nb_local_rows, size_t a_nb_local_cols ) const ;
      virtual bool re_initialize_POST(
                     size_t a_nb_rows, size_t a_nb_cols,
                     size_t a_nb_local_rows, size_t a_nb_local_cols ) const ;

      virtual bool implementation_POST(
                                LA_Implementation const* result ) const ;

      virtual bool is_symmetric_POST( bool result ) const ;

      virtual bool start_local_modifs_PRE( void ) const ;
      virtual bool start_local_modifs_POST( void ) const ;

      virtual bool stop_local_modifs_PRE( void ) const ;
      virtual bool stop_local_modifs_POST( void ) const ;

      virtual bool synchronize_PRE( void ) const ;
      virtual bool synchronize_POST( void ) const ;

      virtual bool row_distribution_PRE( void ) const ;
      virtual bool row_distribution_POST(
                            MAC_DistributedPartition const* result ) const ;

      virtual bool col_distribution_PRE( void ) const ;
      virtual bool col_distribution_POST(
                            MAC_DistributedPartition const* result ) const ;

      virtual bool create_stored_item_iterator_PRE(
                                        MAC_Object* a_owner ) const ;
      virtual bool create_stored_item_iterator_POST(
                                        MAC_Object* a_owner,
                                        LA_MatrixIterator* result ) const ;

      virtual bool create_local_matrix_PRE( MAC_Object* a_owner ) const ;
      virtual bool create_local_matrix_POST( LA_SeqMatrix const* result,
                                             MAC_Object* a_owner ) const ;

      virtual bool nb_stored_items_PRE( void ) const ;

      virtual bool extract_diag_PRE( LA_Vector const* diag ) const ;
      virtual bool extract_diag_POST( LA_Vector const* diag ) const ;

      virtual bool set_item_PRE( size_t i, size_t j ) const ;
      virtual bool set_item_POST( size_t i, size_t j,
                                  LA::SyncState old_state ) const ;

      virtual bool add_to_item_PRE( size_t i, size_t j ) const ;
      virtual bool add_to_item_POST( size_t i, size_t j,
                                     LA::SyncState old_state ) const ;

      virtual bool nullify_PRE( void ) const ;
      virtual bool nullify_POST( void ) const ;

      virtual bool scale_PRE( double alpha ) const ;
      virtual bool scale_POST( double alpha ) const ;

      virtual bool multiply_vec_then_add_PRE(
                                  LA_Vector const* x, LA_Vector* y,
                                  double alpha, double beta ) const ;
      virtual bool multiply_vec_then_add_POST( LA_Vector* y ) const ;

      virtual bool tr_multiply_vec_then_add_PRE(
                                  LA_Vector const* x, LA_Vector* y,
                                  double alpha, double beta ) const ;
      virtual bool tr_multiply_vec_then_add_POST( LA_Vector* y ) const ;

      virtual bool scale_as_diag_mat_mat_PRE( LA_Vector const* lvec ) const ;
      virtual bool scale_as_diag_mat_mat_POST( void ) const ;

      virtual bool scale_as_mat_diag_mat_PRE( LA_Vector const* rvec ) const ;
      virtual bool scale_as_mat_diag_mat_POST( void ) const ;

      virtual bool add_to_diag_PRE( LA_Vector const* vec ) const ;
      virtual bool add_to_diag_POST( LA::SyncState old_state ) const ;

      virtual bool set_PRE( LA_Matrix const* A, bool same_pattern ) const ;
      virtual bool set_POST( LA_Matrix const* A ) const ;

      virtual bool add_Mat_PRE( LA_Matrix const* A, double alpha,
                                bool same_pattern ) const ;
      virtual bool add_Mat_POST( LA::SyncState old_state ) const ;

      virtual bool add_tMat_PRE( LA_Matrix const* A, double alpha ) const ;
      virtual bool add_tMat_POST( LA::SyncState old_state ) const ;

      virtual bool add_Mat_Mat_PRE( LA_Matrix const* A, LA_Matrix const* B,
                                    double alpha  ) const ;
      virtual bool add_Mat_Mat_POST( LA::SyncState old_state ) const ;

      virtual bool readMM_PRE( std::string const& file ) const ;
      virtual bool readMM_POST( void ) const ;

      virtual bool writeMM_PRE( std::string const& file ) const ;
      virtual bool writeMM_POST( void ) const ;

      virtual bool print_items_PRE(
                             std::ostream& os, size_t indent_width ) const ;

      virtual bool create_replica_PRE( MAC_Object* a_owner,
                                       MAC_ModuleExplorer const* exp ) const ;
      virtual bool create_replica_POST( LA_Matrix const* result,
                                        MAC_Object* a_owner,
                                        MAC_ModuleExplorer const* exp ) const ;

   private: //----------------------------------------------------------

      LA_Matrix( void ) ;
      LA_Matrix( LA_Matrix const& other ) ;
      LA_Matrix& operator=( LA_Matrix const& other ) ;

      static MAC_ObjectRegister* plugins_map( void ) ;

      friend class LA_MatrixIterator ;

   //-- Attributes

      bool const IS_PROTO ;
      std::string const NAME ;
      LA::SyncState SYNC_STATE ;
      bool RESIZABLE ;

      bool ONLY_LOCAL_MODIFS ;

      mutable LA::DistributionStrategy DIST_STRAT ;
} ;

#ifndef OUTLINE
   #include <LA_Matrix.icc>
#endif

#endif
