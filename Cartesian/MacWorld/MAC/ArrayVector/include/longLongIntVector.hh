#ifndef LONGLONGINT_VECTOR_HH
#define LONGLONGINT_VECTOR_HH

#include <MAC_assertions.hh>

/*
sequences of values, all of type long long int, ordered according to one index
in a contiguous interval

The values are stored in a block of contiguous storage. Depending on 
the use case, insertion of elements can cause reallocation; that is, the
sequence of values may be stored in a different area after insertion. The
reason is that if there is no room left in the current block for the new
elements, then a larger block is allocated, all of the old elements and the
new elements are copied to the new block, and the old block is deallocated.

Some control can be exercised over when reallocation occurs.

PUBLISHED */

class longLongIntVector
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      // Contruction of an instance whose index interval admits 0 as an
      // inclusive lower limit and `dim' as an exclusive upper limit.
      longLongIntVector( size_t dim, long long int val=0 ) ;      

      // Contruction of an instance whose index interval admits 0 as an
      // inclusive lower limit and `dim' as an exclusive upper limit.
      longLongIntVector( int dim, long long int val=0 ) ; 

      longLongIntVector( longLongIntVector const& other ) ;

      // Assign `other' to `self' 
      // (causes reallocation iff `other.size()'>`capacity()').
      longLongIntVector const& operator=( longLongIntVector const& other ) ;

      // Reinitialize the internal state, as if the constructor with the
      // same argument list was just completed. 
      void re_initialize( size_t dim, long long int val=0 ) ;
           
   //-- Termination

     ~longLongIntVector( void ) ;

   //-- Comparison

      bool operator==( longLongIntVector const& other ) const ;

      bool operator!=( longLongIntVector const& other ) const ;

   //-- Access

      // size of the currently allocated block
      size_t capacity( void ) const ;
      
      // exclusive upper limit of the index interval
      size_t size( void ) const ;
  
      // the smallest index of the occurences of `val' if any, otherwise
      // a number out of the index interval
      size_t index_of( long long int val ) const ;

      // Does `val' appears in `self' ?
      bool has( long long int val ) const ;

      // item of index i
      long long int const& operator()( size_t i ) const ;

      // sum of all the elements of `self'
      long long int sum( void ) const ;

   //-- Element change

      // Change the exclusive upper limit of the index interval to `dim'
      // (with reallocation if `dim'>`capacity()').
      void resize( size_t dim, long long int val=0 ) ;

      // Assign `val' to all items.
      void set( long long int val ) ; 

      // Add `val' at the end of `self' (with possible reallocation).
      void append( long long int val ) ;

      // Ensure that `self' includes `val' (with possible reallocation).
      void extend( long long int val ) ;
   
      // Remove item at place `idx'.
      void remove_at( size_t idx ) ;

      // item of index i
      long long int& operator()( size_t i ) ;

   //-- Input - Output
     
      friend std::ostream& operator<<( std::ostream& out, 
                                       longLongIntVector const& vec ) ;
      
      friend std::istream& operator>>( std::istream& in, 
                                       longLongIntVector& vec ) ;

   //-- Hidden

      // pointer to internal data
      // (for time optimization, should be used with extreme care)
      long long int const* data( void ) const ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      longLongIntVector( void ) ;
      longLongIntVector( double dim ) ;
      longLongIntVector( bool dim ) ;
      longLongIntVector( char dim ) ;

      void allocate( size_t a_size ) ;
      void deallocate( void ) ;

      bool invariant( void ) const ;

   //-- Attributes
      
      long long int* VECTOR ;
      size_t LENGTH ;
      size_t CAPACITY ;
} ;


#ifndef OUTLINE
#include <longLongIntVector.icc>
#endif

#endif
