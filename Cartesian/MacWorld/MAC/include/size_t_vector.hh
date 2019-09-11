#ifndef SIZE_T_VECTOR_HH
#define SIZE_T_VECTOR_HH

#include <MAC_assertions.hh>

class intVector ;

/*
sequences of values, all of type size_t, ordered according to one index
in a contiguous interval

The values are stored in a block of contiguous storage. Depending on 
the use case, insertion of elements can cause reallocation; that is, the
sequence of values may be stored in a different area after insertion. The
reason is that if there is no room left in the current block for the new
elements, then a larger block is allocated, all of the old elements and the
new elements are copied to the new block, and the old block is deallocated.

Some control can be exercised over when reallocation occurs.

PUBLISHED */

class size_t_vector
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      // Contruction of an instance whose index interval admits 0 as an
      // inclusive lower limit and `dim' as an exclusive upper limit.
      size_t_vector( size_t dim, size_t val=0 ) ;

      // Contruction of an instance whose index interval admits 0 as an
      // inclusive lower limit and `dim' as an exclusive upper limit.
      size_t_vector( int dim, size_t val=0 ) ;

      size_t_vector( intVector const& ivec ) ;

      size_t_vector( size_t_vector const& other ) ;

      // Assign `other' to `self' 
      // (causes reallocation iff `other.size()'>`capacity()').
      size_t_vector const& operator=( size_t_vector const& other ) ;

      // Reinitialize the internal state, as if the constructor with the
      // same argument list was just completed 
      // (causes reallocation iff `size()' and `dim' differ ).
      void re_initialize( size_t dim, size_t val=0 ) ;

   //-- Termination

     ~size_t_vector( void ) ;

   //-- Comparison

      bool operator==( size_t_vector const& other ) const ;

      bool operator!=( size_t_vector const& other ) const ;
     
   //-- Access
      
      // size of the currently allocated block
      size_t capacity( void ) const ;

      // exclusive upper limit of the index interval
      size_t size( void ) const ;

      // the smallest index of the occurences of `val' if any, otherwise
      // a number out of the index interval
      size_t index_of( size_t val ) const ;

      // Is there an item of `self' whose value is equal to `val'
      bool has( size_t val ) const ;

      // item of index i
      size_t const& operator()( size_t i ) const ;

      // sum of all the elements of `self'
      size_t sum( void ) const ;
      
   //-- Element change
      
      // Change the exclusive upper limit of the index interval to `dim'
      // (with reallocation if `dim'>`capacity()').
      void resize( size_t dim, size_t val=0 ) ;
  
      // Assign `val' to all items.
      void set( size_t val ) ;

      // Add `val' at the end of `self' (with possible reallocation).
      void append( size_t val ) ;

      // Ensure that `self' includes `val' (with possible reallocation).
      void extend( size_t val ) ;

      // Remove the item of index `idx'.
      void remove_at( size_t idx ) ;

      // item of index `i'
      size_t& operator()( size_t i ) ;

      // Rearrange all the items into increasing order.
      void sort_increasingly( void ) ;

   //-- Input - Output

      friend std::ostream& operator<<( std::ostream& out, 
                                       size_t_vector const& vec ) ;
      
      friend std::istream& operator>>( std::istream& in, 
                                       size_t_vector& vec ) ;
      
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      size_t_vector( void ) ;
      size_t_vector( double dim ) ;
      size_t_vector( bool dim ) ;
      size_t_vector( char dim ) ;

      void allocate( size_t a_size ) ;
      void deallocate( void ) ;

      bool invariant( void ) const ;

   //-- Attributes
      
      size_t* VECTOR ;
      size_t LENGTH ;   
      size_t CAPACITY ;           
} ;


#ifndef OUTLINE
#include <size_t_vector.icc>
#endif

#endif
