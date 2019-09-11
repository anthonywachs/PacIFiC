#ifndef DOUBLE_VECTOR_HH
#define DOUBLE_VECTOR_HH

#include <MAC.hh>
#include <MAC_assertions.hh>

/*
sequences of values, all of type double, ordered according to one index
in a contiguous interval

The values are stored in a block of contiguous storage. Depending on 
the use case, insertion of elements can cause reallocation; that is, the
sequence of values may be stored in a different area after insertion. The
reason is that if there is no room left in the current block for the new
elements, then a larger block is allocated, all of the old elements and the
new elements are copied to the new block, and the old block is deallocated.

Some control can be exercised over when reallocation occurs.

PUBLISHED */

class doubleVector
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      // Contruction of an instance whose index interval admits 0 as an
      // inclusive lower limit and `dim' as an exclusive upper limit.
      doubleVector( size_t dim, double val=0. ) ;

      // Contruction of an instance whose index interval admits 0 as an
      // inclusive lower limit and `dim' as an exclusive upper limit.
      doubleVector( int dim, double val=0. ) ;

      doubleVector( doubleVector const& other ) ;

      // Assign `other' to `self' 
      // (causes reallocation iff `other.size()'>`capacity()').
      doubleVector const& operator=( doubleVector const& other ) ;

      // Reinitialize the internal state, as if the constructor with the
      // same argument list was just completed. 
      void re_initialize( size_t dim, double val=0. ) ;
      
      // Reinitialize the internal state, with the first element being lower,
      // the last element begin upper and the number of elements being ninter
      void re_initialize( size_t ninter, double const& lower, 
      	double const& upper ) ;      

   //-- Termination

     ~doubleVector( void ) ;

   //-- Comparison

      bool operator==( doubleVector const& other ) const ;

      bool operator!=( doubleVector const& other ) const ;
      
   //-- Access

      // size of the currently allocated block
      size_t capacity( void ) const ;

      // exclusive upper limit of the index interval
      size_t size( void ) const ;

      // the smallest index of the itel of `self' whose value
      // is equal to `val' (with respect to `MAC::double_equality')
      // if any, otherwise `MAC::bad_index()'
      size_t index_of( double val, double a_dbl_eps = MAC::epsilon_double(),
                                   double a_dbl_min = MAC::min_double() ) const ;

      // Is there an item of `self' whose value is equal to `val'
      // (with respect to `MAC::double_equality') ?
      bool has( double val, double a_dbl_eps = MAC::epsilon_double(),
                            double a_dbl_min = MAC::min_double() ) const ;

      // item of index i
      double const& operator()( size_t i ) const ;

   //-- Element change

      // Change the exclusive upper limit of the index interval to `dim'
      // (with reallocation if `dim'>`capacity()').
      void resize( size_t dim, double val=0. ) ;

      // Assign `val' to all items.
      void set( double val ) ; 

      // Add `val' at the end of `self' (with possible reallocation).
      void append( double val ) ;

      // Remove the item of index `idx'.
      void remove_at( size_t idx ) ;

      // item of index `i'
      double& operator()( size_t i ) ;

   //-- Input - Output
                  
      friend std::ostream& operator<<( std::ostream& out, 
                                       doubleVector const& vec ) ;
      
      friend std::istream& operator>>( std::istream& in, 
                                       doubleVector& vec ) ;

   //-- Hidden

      // pointer to internal data
      // (for time optimization, should be used with extreme care)
      double const* data( void ) const ;
      
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      doubleVector( void ) ;
      doubleVector( double dim ) ;
      doubleVector( bool dim ) ;
      doubleVector( char dim ) ;

      void allocate( size_t a_size ) ;
      void deallocate( void ) ;

      bool invariant( void ) const ;

   //-- Attributes
      
      double* VECTOR ;
      size_t LENGTH ;
      size_t CAPACITY ;
} ;


#ifndef OUTLINE
#include <doubleVector.icc>
#endif

#endif 
