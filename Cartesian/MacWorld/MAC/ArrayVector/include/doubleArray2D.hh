#ifndef DOUBLE_ARRAY_2D_HH
#define DOUBLE_ARRAY_2D_HH

#include <doubleVector.hh>

/*
sequences of values, all of type double, ordered according to two indices
in contiguous intervals
*/

class doubleArray2D
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      // Contruction of an instance whose all index intervals admit 0 as an
      // inclusive lower limit and respectively `dim0' and `dim1' as an 
      // exclusive upper limit.
      doubleArray2D( size_t dim0, size_t dim1, double val=0. ) ;

      doubleArray2D( doubleArray2D const& other ) ;

      doubleArray2D const& operator=( doubleArray2D const& other ) ;

      // Reinitialize the internal state, as if the constructor with the
      // same argument list was just completed. 
      void re_initialize( size_t dim0, size_t dim1, double val=0. ) ;

   //-- Termination

     ~doubleArray2D( void ) ;

   //-- Comparison

      bool operator==( doubleArray2D const& other ) const ;

      bool operator!=( doubleArray2D const& other ) const ;
      
   //-- Access

      // exclusive upper limit of the `an_index'-th index interval
      size_t index_bound( size_t an_index ) const ;

      // Reinitialize `x' by copying all items of self whose `an_index'-th 
      // index is equal to `index-value'.
      void extract_section( size_t an_index,
                            size_t index_value,
                            doubleVector& x ) const ;

      // item of indices `i0' and `i1'
      double const& operator()( size_t i0, size_t i1 ) const ;

   //-- Element change

      // Increase the exclusive upper limit of the first index interval 
      // to `dim0'.
      void raise_first_index_bound( size_t dim0, double val=0. ) ;

      // Assign `val' to all items.
      void set( double val ) ; 

      // Assign `x' items to the items of self whose `an_index'-th index is 
      // equal to `index-value' and whose other index is given by `x'.
      void set_section( size_t an_index, 
                        size_t index_value, 
                        doubleVector const& x ) ;

      // item of indices `i0' and `i1'
      double& operator()( size_t i0, size_t i1 ) ;

   //-- Input - Output
     
      friend std::ostream& operator<<( std::ostream& out, 
                                       doubleArray2D const& a ) ;

   //-- Hidden

      // pointer to internal data
      // (for time optimization, should be used with extreme care)
      double const* data( void ) const ;
      
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      doubleArray2D( void ) ;

   //-- Attributes
      
      doubleVector vector ;
      size_t d0 ; 
      size_t d1 ;
             
} ;


#ifndef OUTLINE
#include <doubleArray2D.icc>
#endif

#endif 
