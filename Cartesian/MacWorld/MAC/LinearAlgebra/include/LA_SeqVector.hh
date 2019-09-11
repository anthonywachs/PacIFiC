#ifndef LA_SEQ_VECTOR_HH
#define LA_SEQ_VECTOR_HH

#include <LA_Vector.hh>
#include <LA_SeqScatter.hh>

#include <iosfwd>

#include <doubleVector.hh>
#include <MAC_assertions.hh>

class MAC_ModuleExplorer ;

/*
Sequential vectors of double, with a possible block structure in subvectors,
compatibles with `LA_SeqMatrix::' matrices.
*/

class LA_SeqVector : public LA_Vector
{
   public: //-----------------------------------------------------------

  //-- Instance delivery and initialization

      // Create and return an instance.
      static LA_SeqVector* create( MAC_Object* a_owner, size_t a_nb_rows ) ;

      // Create and return an instance initialized with `dvec'.
      static LA_SeqVector* create( MAC_Object* a_owner,
                                   doubleVector const& dvec ) ;

      virtual LA_SeqVector* create_vector( MAC_Object* a_owner ) const ;

      virtual void re_initialize( size_t a_nb_rows,
                                  size_t a_nb_local_rows = MAC::bad_index() ) ;

   //-- Characteristics

      virtual LA_Implementation const* implementation( void ) const ;

   //-- Distributed processing

      virtual bool is_desynchronizable( void ) const ;

      virtual MAC_DistributedPartition const* row_distribution( void ) const ;

      virtual LA_SeqScatter* create_scatter( MAC_Object* a_owner,
                                             size_t_vector const& from,
                                             size_t_vector const& to ) const ;

   //-- Access

      virtual LA_SeqVector* create_local_vector( MAC_Object* a_owner ) const ;

      // `i'-th coefficient
      virtual double item( size_t i ) const ;

   //-- Element change

      virtual void set_item( size_t i, double x ) ;

      virtual void add_to_item( size_t i, double x ) ;

   //-- BLAS level 1 : vector-vector operators

      virtual double dot( LA_Vector const* a ) const ;

      virtual double two_norm( void ) const ;

      virtual double max_norm( void ) const ;

      virtual void nullify( void ) ;

      virtual void scale( double alpha ) ;

      virtual void set( LA_Vector const* a ) ;

      virtual void set( double value ) ;

      virtual void set_as_v_product( LA_Vector const* a, LA_Vector const* b ) ;

      virtual void set_as_reciprocal( LA_Vector const* a,
                                      double smallest_inverted_item=0.0,
                                      double default_value=1.0 ) ;

      virtual void sum( LA_Vector const* a, double alpha=1.0 ) ;

   //-- Input - Output

      virtual void print_items( std::ostream& os, size_t indent_width ) const ;

      virtual void write( std::string const& filename ) const ;

      // Save `self' in `mod'.
      void save( MAC_Module* mod ) const ;

      // Restore `self' from module described by `exp'.
      void restore( MAC_ModuleExplorer const* exp ) ;

   //-- Statistical tools

      // mean value of the elements of `self'
      double mean( void ) const ;

      // standard deviation of the elements of `self'
      double standard_deviation( void ) const ;

   //-- Hidden

      void set( size_t a_nb_rows, double* values ) ;

      // pointer to internal data
      // (for time optimization, should be used with care)
      double const* data( void ) const ;
      double* data( void ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      LA_SeqVector( void ) ;
     ~LA_SeqVector( void ) ;
      LA_SeqVector( LA_SeqVector const& other ) ;
      LA_SeqVector& operator=( LA_SeqVector const& other ) ;

      LA_SeqVector( MAC_Object* a_owner, size_t a_nb_rows ) ;

      friend class LA_SeqMatrix ;
      void make_desynchronizable( void ) ;

      void re_initialize_with_global_sizes( size_t a_nb_rows ) ;
      bool re_initialize_with_global_sizes_PRE( size_t a_nb_rows ) const ;
      bool re_initialize_with_global_sizes_POST(size_t a_nb_rows ) const ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

      virtual bool implementation_POST(
                                LA_Implementation const* result ) const ;

   //-- Attributes

      bool UNSYNCHRO ;
      mutable MAC_DistributedPartition* ROW_DIST ;

      // Vector table :
      bool OWNS_DATA ;
      double* DATA ;
} ;

#ifndef OUTLINE
   #include <LA_SeqVector.icc>
#endif

#endif
