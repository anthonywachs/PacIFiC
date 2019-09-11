#ifndef INT_ARRAY_3D_HH
#define INT_ARRAY_3D_HH

#include <intVector.hh>

/* 
sequences of values, all of type int, ordered according to three indices
in contiguous intervals
*/

class intArray3D
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      // Contruction of an instance whose all index intervals admit 0 as an
      // inclusive lower limit and respectively `dim0', `dim1' and `dim2'  
      // as an exclusive upper limit.
      intArray3D( size_t dim0, size_t dim1, size_t dim2, int val=0 ) ;

      intArray3D( intArray3D const& other ) ;

      intArray3D& operator=( intArray3D const& other ) ;

      // Reinitialize the internal state, as if the constructor with the
      // same argument list was just completed. 
      void re_initialize( size_t dim0, size_t dim1, size_t dim2, int val=0 ) ;

   //-- Termination

     ~intArray3D( void ) ;

   //-- Comparison

      bool operator==( intArray3D const& other ) const ;
      bool operator!=( intArray3D const& other ) const ;
        
   //-- Access

      // exclusive upper limit of the `an_index'-th index interval
      size_t index_bound( size_t an_index ) const ;

      // item of indices `i0', 'i1' and `i2'
      int const& operator()( size_t i0, size_t i1, size_t i2 ) const ;

   //-- Element change

      // Increase the exclusive upper limit of the first index interval 
      // to `dim0'.
      void raise_first_index_bound( size_t dim0, int val=0 ) ;
  
      // Assign `val' to all items.
      void set( int val ) ; 

      // item of indices `i0', 'i1' and `i2'
      int& operator()( size_t i0, size_t i1, size_t i2 ) ;

   //-- Input - Output
     
      friend std::ostream& operator<<( std::ostream& out, 
                                       intArray3D const& a ) ;
      
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      intArray3D( void ) ;

   //-- Attributes
         
      intVector vector ;
      size_t d0 ; 
      size_t d1 ;
      size_t d2 ;
             
} ;


#ifndef OUTLINE
#include <intArray3D.icc>
#endif

#endif 
