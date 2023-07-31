#ifndef LONGLONGINT_ARRAY_3D_HH
#define LONGLONGINT_ARRAY_3D_HH

#include <longLongIntVector.hh>

/* 
sequences of values, all of type int, ordered according to three indices
in contiguous intervals
*/

class longLongIntArray3D
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      // Contruction of an instance whose all index intervals admit 0 as an
      // inclusive lower limit and respectively `dim0', `dim1' and `dim2'  
      // as an exclusive upper limit.
      longLongIntArray3D( size_t dim0, size_t dim1, size_t dim2, 
      	long long int val=0 ) ;

      longLongIntArray3D( longLongIntArray3D const& other ) ;

      longLongIntArray3D& operator=( longLongIntArray3D const& other ) ;

      // Reinitialize the internal state, as if the constructor with the
      // same argument list was just completed. 
      void re_initialize( size_t dim0, size_t dim1, size_t dim2, 
      	long long int val=0 ) ;

   //-- Termination

     ~longLongIntArray3D( void ) ;

   //-- Comparison

      bool operator==( longLongIntArray3D const& other ) const ;
      bool operator!=( longLongIntArray3D const& other ) const ;
        
   //-- Access

      // exclusive upper limit of the `an_index'-th index interval
      size_t index_bound( size_t an_index ) const ;

      // item of indices `i0', 'i1' and `i2'
      long long int const& operator()( size_t i0, size_t i1, size_t i2 ) const ;

   //-- Element change

      // Increase the exclusive upper limit of the first index interval 
      // to `dim0'.
      void raise_first_index_bound( size_t dim0, long long int val=0 ) ;
  
      // Assign `val' to all items.
      void set( long long int val ) ; 

      // item of indices `i0', 'i1' and `i2'
      long long int& operator()( size_t i0, size_t i1, size_t i2 ) ;

   //-- Input - Output
     
      friend std::ostream& operator<<( std::ostream& out, 
                                       longLongIntArray3D const& a ) ;
      
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      longLongIntArray3D( void ) ;

   //-- Attributes
         
      longLongIntVector vector ;
      size_t d0 ; 
      size_t d1 ;
      size_t d2 ;
             
} ;


#ifndef OUTLINE
#include <longLongIntArray3D.icc>
#endif

#endif 
