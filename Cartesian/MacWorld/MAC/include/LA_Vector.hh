#ifndef LA_VECTOR_HH
#define LA_VECTOR_HH

#include <MAC.hh>
#include <MAC_Object.hh>
#include <MAC_ObjectWriter.hh>
#include <MAC_ObjectReader.hh>
#include <MAC_assertions.hh>

#include <LA.hh>

#include <iosfwd>
#include <string>

class MAC_DistributedPartition ;
class size_t_vector ;

class LA_SeqVector ;
class LA_Implementation ;
class LA_Scatter ;

/*
  Vectors of double, for linear algebra computations.

  PUBLISHED
*/

class LA_Vector : public MAC_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      // Reinitialize the internal state, as if it was just created.
      virtual void re_initialize( size_t a_nb_rows,
                        size_t a_nb_local_rows = MAC::bad_index() ) = 0 ;

      // Create a vector with the same implementation as `self'.
      virtual LA_Vector* create_vector( MAC_Object* a_owner ) const = 0 ;

   //-- Characteristics

      // implementation family descriptor
      virtual LA_Implementation const* implementation( void ) const = 0 ;

      // Might the number of rows be changed ?
      bool is_resizable( void ) const ;

   //-- Characteristics setting

      // Fixe the number of rows and the number of columns.
      void make_non_resizable( void ) ;

   //-- Distributed processing(10.)

      // on-process synchronization state
      LA::SyncState state( void ) const ;

      LA::DistributionStrategy distribution_strategy( void ) const ;

      // Is `self' a distributed vector ?
      virtual bool is_desynchronizable( void ) const = 0 ;

      // Define environment between call to "start_local_modifs" and
      // "stop_local_modifs" methods
      // in which "set_item", "add_to_item" or "item" can be used alternatively
      // to modify on-process items without synchronization.
      virtual void start_local_modifs( void ) ;

      virtual void stop_local_modifs( void ) ;

      // Is in environment defined by start_local_modifs and
      // stop_local_modifs methods ?
      bool only_local_modifs( void ) const ;

      // In distributed case, are non local items of any process to be copied
      // to other process ?
      bool is_synchronized( void ) const ;

      // Ensure that no dispersed items are remaining in `self'.
      virtual void synchronize( void ) ;

      // distribution of vector rows over the processes
      virtual MAC_DistributedPartition const* row_distribution( void ) 
      	const = 0 ;

      // Create scatter object compatible with `self'.
      virtual LA_Scatter* create_scatter( MAC_Object* a_owner,
                                          size_t_vector const& from,
                                          size_t_vector const& to ) const = 0 ;

   //-- Access(50.)

      // number of rows
      size_t nb_rows( void ) const ;

      // number of rows handled by current process
      size_t nb_local_rows( void ) const ;
      
      // Create and return a COPY of `self' :
      //    - in parallel : current process rows only ;
      //    - to use with care : potentially greedy in time and memory.
      virtual LA_SeqVector* create_local_vector(
                                           MAC_Object* a_owner ) const = 0 ;

      // `i'-th coefficient
      virtual double item( size_t i ) const = 0 ;

   //-- Element change(60.)

      // Assign `x' to the `i'-th coefficient.
      virtual void set_item( size_t i, double x ) = 0 ;

      // Add `x' to the `i'-th coefficient.
      virtual void add_to_item( size_t i, double x ) = 0 ;

   //-- BLAS level 1 : vector-vector operators(105.)

      // euclidian inner-product with `a'
      virtual double dot( LA_Vector const* a ) const = 0 ;

      // euclidian norm
      virtual double two_norm( void ) const = 0 ;

      // infinity norm
      virtual double max_norm( void ) const = 0 ;

      // Assign 0.0 to all coefficients.
      virtual void nullify( void ) ;

      // Multiply by `alpha' :
      //    `self' <- `alpha'*`self' .
      virtual void scale( double alpha ) = 0 ;

      // Reinitialize by copying all coefficients of `a'.
      virtual void set( LA_Vector const* a ) = 0 ;

      // Reinitialize by setting all coefficients to `value'.
      virtual void set( double value ) = 0 ;

      // Reinitialize by assigning to each coefficient the product
      // of the coefficients of `a' and `b' with the same index :
      //    `self'(i) <- `a'(i) * `b'(i) .
      virtual void set_as_v_product(
                             LA_Vector const* a, LA_Vector const* b ) = 0 ;

      // Reinitialize by assigning to each coefficient the inverse
      // of the coefficients of `a' with the same index.
      // If such coefficient is such that its magnitude is strictly lower than
      // `smallest_inverted_item' then is't value is set to `default_value'.
      // Rem : `a' can be the same than `self'.
      virtual void set_as_reciprocal( LA_Vector const* a,
                                      double smallest_inverted_item = 0.,
                                      double default_value = 1. ) = 0 ;

      // Add `alpha'*`a' to `self' :
      //    `self' <- ( `alpha'*`a' + `self' ) .
      virtual void sum( LA_Vector const* a, double alpha=1.0 ) = 0 ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;

      // Print items of `self'.
      virtual void print_items(
                      std::ostream& os, size_t indent_width ) const = 0 ;

      // Write in file called `filename'.
      virtual void write( std::string const& filename ) const = 0 ;

      // Reinitialize with vector stored in file called `filename'.
      virtual void read( std::string const& filename ) ;

   //-- Persistence

      // Use `writer' to store `self' so that ist can be retrieved
      // with `::restore_state'.
      virtual void save_state( MAC_ObjectWriter* writer ) const ;
      
      // Retrieve `self' from `reader'.
      virtual void restore_state( MAC_ObjectReader* reader ) ;
      
      // Do a fake read state for the reader to proceed to the next object
      static void read_state_nonrestored( MAC_ObjectReader* reader ) ; 
      
      // Set the vector name
      void set_name( std::string const& name );
      
      // Get the vector name
      std::string const& name( void ) const;      

   protected: //--------------------------------------------------------

      virtual ~LA_Vector( void ) ;

      LA_Vector( MAC_Object* a_owner, size_t a_nb_rows ) ;

      // To avoid virtuality of method nb_rows.
      void set_rows_number( size_t a_nb_rows ) ;

   //-- Distributed processing

      void set_distribution_strategy(
                                LA::DistributionStrategy dist_strat ) const ;

      void set_unsynchronized_state( LA::SyncState new_state ) ;

      void set_only_local_modifs_state( bool only_local ) ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

      virtual bool re_initialize_PRE( size_t a_nb_rows,
                                      size_t a_nb_local_rows ) const ;
      virtual bool re_initialize_POST( size_t a_nb_rows,
                                       size_t a_nb_local_rows ) const ;

      virtual bool create_vector_PRE( MAC_Object* a_owner ) const ;
      virtual bool create_vector_POST(
                   LA_Vector* result, MAC_Object* a_owner ) const ;

      virtual bool item_PRE( size_t i ) const ;

      virtual bool implementation_POST(
                                LA_Implementation const* result ) const ;

      virtual bool start_local_modifs_PRE( void ) const ;
      virtual bool start_local_modifs_POST( void ) const ;

      virtual bool stop_local_modifs_PRE( void ) const ;
      virtual bool stop_local_modifs_POST( void ) const ;

      virtual bool synchronize_PRE( void ) const ;
      virtual bool synchronize_POST( void ) const ;

      virtual bool row_distribution_PRE( void ) const ;
      virtual bool row_distribution_POST(
                            MAC_DistributedPartition const* result ) const ;

      virtual bool create_scatter_PRE( MAC_Object* a_owner,
                                       size_t_vector const& from,
                                       size_t_vector const& to ) const ;
      virtual bool create_scatter_POST( LA_Scatter const* result,
                                        MAC_Object const* a_owner,
                                        size_t_vector const& from,
                                        size_t_vector const& to ) const ;

      virtual bool create_local_vector_PRE( MAC_Object* a_owner ) const ;
      virtual bool create_local_vector_POST(
                      LA_SeqVector* result, MAC_Object* a_owner ) const ;

      virtual bool set_item_PRE( size_t i ) const ;
      virtual bool set_item_POST( size_t i, LA::SyncState old_state ) const ;

      virtual bool add_to_item_PRE( size_t i ) const ;
      virtual bool add_to_item_POST( size_t i, LA::SyncState old_state ) const ;

      virtual bool dot_PRE( LA_Vector const* a ) const ;
      virtual bool dot_POST( double result ) const ;

      virtual bool two_norm_PRE( void ) const ;
      virtual bool two_norm_POST( double result ) const ;

      virtual bool max_norm_PRE( void ) const ;
      virtual bool max_norm_POST( double result ) const ;

      virtual bool nullify_PRE( void ) const ;
      virtual bool nullify_POST( void ) const ;

      virtual bool scale_PRE( double value ) const ;
      virtual bool scale_POST( double value ) const ;

      virtual bool set_PRE( LA_Vector const* a ) const ;
      virtual bool set_POST( LA_Vector const* a ) const ;

      virtual bool set_PRE( double value ) const ;
      virtual bool set_POST( double value ) const ;

      virtual bool set_as_v_product_PRE(
                      LA_Vector const* a, LA_Vector const* b ) const ;
      virtual bool set_as_v_product_POST(
                      LA_Vector const* a, LA_Vector const* b ) const ;

      virtual bool set_as_reciprocal_PRE( LA_Vector const* a,
                                          double smallest_inverted_item,
                                          double default_value ) const ;
      virtual bool set_as_reciprocal_POST( LA_Vector const* a ) const ;

      virtual bool sum_PRE( LA_Vector const* a, double alpha ) const ;
      virtual bool sum_POST( LA_Vector const* a, double alpha ) const ;

      virtual bool print_items_PRE(
                          std::ostream& os, size_t indent_width ) const ;

      virtual bool write_PRE( std::string const& filename ) const ;
      virtual bool write_POST( void ) const ;

      virtual bool read_PRE( std::string const& filename ) const ;
      virtual bool read_POST( void ) const ;            

   private: //----------------------------------------------------------

      LA_Vector( void ) ;
      LA_Vector( LA_Vector const& other ) ;
      LA_Vector& operator=( LA_Vector const& other ) ;

   //-- Attributes

      size_t NB_ROWS ;
      LA::SyncState DIST_STATUS ;
      bool RESIZABLE ;

      bool ONLY_LOCAL_MODIFS ;

      mutable LA::DistributionStrategy DIST_STRAT ;
      
      std::string VECNAME;
} ;

#ifndef OUTLINE
   #include <LA_Vector.icc>
#endif

#endif
