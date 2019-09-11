#ifndef DOUBLE_ARRAY_6D_HH
#define DOUBLE_ARRAY_6D_HH

#include <doubleVector.hh>

/*
sequences of values, all of type double, ordered according to six indices
in contiguous intervals
*/

class doubleArray6D
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      // Contruction of an instance whose all index intervals admit 0 as an
      // inclusive lower limit and respectively `dim0', `dim1', `dim2', 
      // `dim3', `dim4' and `dim5' as an exclusive upper limit.
      doubleArray6D( size_t dim0, size_t dim1, size_t dim2, size_t dim3,
                     size_t dim4, size_t dim5, double val=0. ) ;

      // Reinitialize the internal state, as if the constructor with the
      // same argument list was just completed. 
      void re_initialize( size_t dim0, size_t dim1, size_t dim2, size_t dim3,
                          size_t dim4, size_t dim5, double val=0. ) ;

   //-- Termination

     ~doubleArray6D( void ) ;

   //-- Access

      // exclusive upper limit of the `an_index'-th index interval
      size_t index_bound( size_t an_index ) const ;
  
      // item of indices `i0', 'i1', 'i2', 'i3', 'i4' and `i5'
      double const& operator()( size_t i0, size_t i1, size_t i2, size_t i3, 
                                size_t i4, size_t i5 ) const ;

   //-- Element change

      // Increase the exclusive upper limit of the first index interval 
      // to `dim0'.
      void raise_first_index_bound( size_t dim0, double val=0. ) ;

      // Assign `val' to all items.
      void set( double val ) ; 

      // item of indices `i0', 'i1', 'i2', 'i3', 'i4' and `i5'
      double& operator()( size_t i0, size_t i1, size_t i2, size_t i3, 
                          size_t i4, size_t i5 ) ;

   //-- Input - Output
     
      friend std::ostream& operator<<( std::ostream& out, 
                                       doubleArray6D const& a ) ;
      
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      doubleArray6D( void ) ;

      doubleArray6D( doubleArray6D const& other ) ;
      doubleArray6D& operator=( doubleArray6D const& other ) ;
      
      bool operator==( doubleArray6D const& other ) const ;
      bool operator!=( doubleArray6D const& other ) const ;      

   //-- Attributes
     
      doubleVector vector ;
      size_t d0 ; 
      size_t d1 ;
      size_t d2 ;
      size_t d3 ;
      size_t d4 ;
      size_t d5 ;
             
} ;


#ifndef OUTLINE
#include <doubleArray6D.icc>
#endif

#endif 
