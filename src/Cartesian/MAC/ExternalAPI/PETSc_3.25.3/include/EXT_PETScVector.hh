#ifndef EXT_PETSC_VECTOR_HH
#define EXT_PETSC_VECTOR_HH

#include <MAC_Object.hh>
#include <MAC_ObjectWriter.hh>
#include <MAC_ObjectReader.hh>
#include <iosfwd>

#include <MAC_assertions.hh>
#include <MAC_Communicator.hh>
#include <MAC_DistributedPartition.hh>
#include <doubleVector.hh>

#include <LA_SeqVector.hh>

#include <EXT_PETScScatter.hh>
#include <EXT_PETScAPI.hh>

class size_t_vector ;

/* PETSc vectors.
*/

class EXT_PETScVector : public LA_Vector
{
   public: //-----------------------------------------------------------

  //-- Instance delivery and initialization

      // Create and return an instance.
      static EXT_PETScVector* create( MAC_Object* a_owner,
                                      bool sequential,
                                      size_t a_nb_rows ) ;

      virtual EXT_PETScVector* create_vector( MAC_Object* a_owner ) const ;

      virtual void re_initialize( size_t a_nb_rows,
                                  size_t a_nb_local_rows = MAC::bad_index() ) ;

   //-- Characteristics

      virtual LA_Implementation const* implementation( void ) const ;

   //-- Distributed processing

      virtual bool is_desynchronizable( void ) const ;

      virtual void start_local_modifs( void ) ;

      virtual void stop_local_modifs( void ) ;

      virtual void synchronize( void ) ;

      virtual MAC_DistributedPartition const* row_distribution( void ) const ;

      virtual EXT_PETScScatter* create_scatter( MAC_Object* a_owner,
		size_t_vector const& from,
		size_t_vector const& to ) const ;

   //-- Access

      virtual LA_SeqVector* create_local_vector( MAC_Object* a_owner ) const ;

      // `i'-th coefficient
      virtual double item( size_t i ) const ;

   //-- Element change

      virtual void set_item( size_t i, double x ) ;

      virtual void add_to_item( size_t i, double x ) ;

      void set_verbosity( bool verbosity ) ;

   //-- BLAS level 1 : vector-vector operators

      virtual double dot( LA_Vector const* a ) const ;

      virtual double two_norm( void ) const ;

      virtual double max_norm( void ) const ;

      virtual void scale( double alpha ) ;

      virtual void set( LA_Vector const* a ) ;

      virtual void set( double value ) ;

      virtual void set_as_v_product( LA_Vector const* a,
                                     LA_Vector const* b ) ;

      virtual void set_as_reciprocal( LA_Vector const* a,
                                      double smallest_inverted_item=0.0,
                                      double default_value=1.0 ) ;

      virtual void sum( LA_Vector const* a, double alpha=1.0 ) ;

      LA_SeqVector const* local_vector( void ) const ;

   //-- Input - Output

      virtual void print_items( std::ostream& os, size_t indent_width ) const ;

      virtual void write( std::string const& filename ) const ;

   //-- Hidden

      // pointer to internal data
      // (for time optimization, should be used with care)
      Vec& vector( void ) const ;

   //-- Persistence

      // Use `writer' to store `self' so that ist can be retrieved
      // with `::restore_state'.
      virtual void save_state( MAC_ObjectWriter* writer ) const ;
      
      // Retrieve `self' from `reader'.
      virtual void restore_state( MAC_ObjectReader* reader ) ;          

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      EXT_PETScVector( void ) ;
     ~EXT_PETScVector( void ) ;
      EXT_PETScVector( EXT_PETScVector const& other ) ;
      EXT_PETScVector& operator=(
                            EXT_PETScVector const& other ) ;

      EXT_PETScVector( MAC_Object* a_owner, bool sequential,
                       size_t a_nb_rows ) ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

      virtual bool implementation_POST(
                                LA_Implementation const* result ) const ;

   //-- Attributes

      Vec VECTOR ;
      MAC_DistributedPartition* DIST ;
      bool SEQ ;
      size_t FIRST ;
      size_t LAST ;
} ;


#endif
