#ifndef BOOL_VECTOR_HH
#define BOOL_VECTOR_HH

#include <MAC_assertions.hh>

/*
sequences of values, all of type bool, ordered according to one index
in a contiguous interval

The values are stored in a block of contiguous storage. Depending on 
the use case, insertion of elements can cause reallocation; that is, the
sequence of values may be stored in a different area after insertion. The
reason is that if there is no room left in the current block for the new
elements, then a larger block is allocated, all of the old elements and the
new elements are copied to the new block, and the old block is deallocated.

Some control can be exercised over when reallocation occurs.

PUBLISHED */

class boolVector
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      // Contruction of an instance whose index interval admits 0 as an
      // inclusive lower limit and `dim' as an exclusive upper limit.
      boolVector( size_t dim, bool val=false ) ;

      // Contruction of an instance whose index interval admits 0 as an
      // inclusive lower limit and `dim' as an exclusive upper limit.
      boolVector( int dim, bool val=false ) ;

      boolVector( boolVector const& other ) ;

      // Assign `other' to `self' 
      // (causes reallocation iff `other.size()'>`capacity()').
      boolVector const& operator=( boolVector const& other ) ;

      // Reinitialize the internal state, as if the constructor with the
      // same argument list was just completed 
      // (causes reallocation iff `size()' and `dim' differ ). 
      void re_initialize( size_t dim, bool val=false ) ;

   //-- Termination

     ~boolVector( void ) ;

   //-- Comparison

     bool operator==( boolVector const& other ) const ;
     
     bool operator!=( boolVector const& other ) const ;
       
   //-- Access

      // size of the currently allocated block
      size_t capacity( void ) const ;
      
      // exclusive upper limit of the index interval
      size_t size( void ) const ;

      // item of index i
      bool const& operator()( size_t i ) const ;
      
      // Does `val' appears in `self' ?
      bool has( bool val ) const ;
      
   //-- Element change

      // Change the exclusive upper limit of the index interval to `dim'
      // (with reallocation if `dim'>`capacity()').
      void resize( size_t dim, bool val=false ) ;
 
      // Assign `val' to all items.
      void set( bool val ) ; 

      // Include `val' at the end of `self'  (with possible reallocation).
      void append( bool val ) ;

      // Remove item at place `idx'.
      void remove_at( size_t idx ) ;

      // item of index i
      bool& operator()( size_t i ) ;

   //-- Input - Output
     
      friend std::ostream& operator<<( std::ostream& out, 
                                       boolVector const& vec ) ;
      
      friend std::istream& operator>>( std::istream& in, 
                                       boolVector& vec ) ;
                  
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      boolVector( void ) ;
      boolVector( double dim ) ;
      boolVector( bool dim ) ;
      boolVector( char dim ) ;

      void allocate( size_t a_size ) ;
      void deallocate( void ) ;

      bool invariant( void ) const ;

   //-- Attributes

      bool* VECTOR ;
      size_t LENGTH ;
      size_t CAPACITY ;              
} ;


#ifndef OUTLINE
#include <boolVector.icc>
#endif

#endif
