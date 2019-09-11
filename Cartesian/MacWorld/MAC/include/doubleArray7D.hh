#ifndef DOUBLE_ARRAY_7D_HH
#define DOUBLE_ARRAY_7D_HH

#include <doubleVector.hh>

/*
sequences of values, all of type double, ordered according to seven indices
in contiguous intervals
*/

class doubleArray7D
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      // Contruction of an instance whose all index intervals admit 0 as an
      // inclusive lower limit and respectively `dim0', `dim1', `dim2', 
      // `dim3', `dim4', `dim5' and `dim6' as an exclusive upper limit.
      doubleArray7D( size_t dim0, size_t dim1, size_t dim2, size_t dim3,
                     size_t dim4, size_t dim5, size_t dim6,
                     double val=0. ) ;

      // Reinitialize the internal state, as if the constructor with the
      // same argument list was just completed. 
      void re_initialize( size_t dim0, size_t dim1, size_t dim2, size_t dim3,
                          size_t dim4, size_t dim5, size_t dim6,
                          double val=0. ) ;

   //-- Termination

     ~doubleArray7D( void ) ;

   //-- Access

      // exclusive upper limit of the `an_index'-th index interval
      size_t index_bound( size_t an_index ) const ;
  
      // item of indices `i0', 'i1', 'i2', 'i3', 'i4', 'i5' and `i6'
      double const& operator()( size_t i0, size_t i1, size_t i2, size_t i3, 
                                size_t i4, size_t i5, size_t i6 ) const ;

   //-- Element change

      // Increase the exclusive upper limit of the first index interval 
      // to `dim0'.
      void raise_first_index_bound( size_t dim0, double val=0. ) ;

      // Assign `val' to all items.
      void set( double val ) ; 

      // item of indices `i0', 'i1', 'i2', 'i3', 'i4', 'i5' and `i6'
      double& operator()( size_t i0, size_t i1, size_t i2, size_t i3, 
                          size_t i4, size_t i5, size_t i6 ) ;


   //-- Input - Output
     
      friend std::ostream& operator<<( std::ostream& out, 
                                       doubleArray7D const& a ) ;
      
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      doubleArray7D( void ) ;

      doubleArray7D( doubleArray7D const& other ) ;
      doubleArray7D& operator=( doubleArray7D const& other ) ;
      
      bool operator==( doubleArray7D const& other ) const ;
      bool operator!=( doubleArray7D const& other ) const ;

   //-- Attributes      
  
      doubleVector vector ;
      size_t d0 ; 
      size_t d1 ;
      size_t d2 ;
      size_t d3 ;
      size_t d4 ;
      size_t d5 ;
      size_t d6 ;
             
} ;


#ifndef OUTLINE
#include <doubleArray7D.icc>
#endif

#endif 
