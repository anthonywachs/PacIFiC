#ifndef STRING_ARRAY_2D_HH
#define STRING_ARRAY_2D_HH

#include <stringVector.hh>
#include <string>

/*
sequences of values, all of type string, ordered according to two indices
in contiguous intervals
*/

class stringArray2D
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      // Contruction of an instance whose all index intervals admit 0 as an
      // inclusive lower limit and respectively `dim0' and `dim1' as an 
      // exclusive upper limit.
      stringArray2D( size_t dim0, size_t dim1, std::string val="" ) ;

      stringArray2D( stringArray2D const& other ) ;

      stringArray2D const& operator=( stringArray2D const& other ) ;

      // Reinitialize the internal state, as if the constructor with the
      // same argument list was just completed. 
      void re_initialize( size_t dim0, size_t dim1, std::string val="" ) ;

   //-- Termination

     ~stringArray2D( void ) ;

   //-- Comparison

      bool operator==( stringArray2D const& other ) const ;

      bool operator!=( stringArray2D const& other ) const ;
      
   //-- Access

      // exclusive upper limit of the `an_index'-th index interval
      size_t index_bound( size_t an_index ) const ;

      // Reinitialize `x' by copying all items of self whose `an_index'-th 
      // index is equal to `index-value'.
      void extract_section( size_t an_index,
                            size_t index_value,
                            stringVector& x ) const ;

      // item of indices `i0' and `i1'
      std::string const& operator()( size_t i0, size_t i1 ) const ;

   //-- Element change

      // Increase the exclusive upper limit of the first index interval 
      // to `dim0'.
      void raise_first_index_bound( size_t dim0, std::string val="" ) ;

      // Assign `val' to all items.
      void set( std::string val ) ; 

      // Assign `x' items to the items of self whose `an_index'-th index is 
      // equal to `index-value' and whose other index is given by `x'.
      void set_section( size_t an_index, 
                        size_t index_value, 
                        stringVector const& x ) ;

      // item of indices `i0' and `i1'
      std::string& operator()( size_t i0, size_t i1 ) ;

   //-- Input - Output
     
      friend std::ostream& operator<<( std::ostream& out, 
                                       stringArray2D const& a ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      stringArray2D( void ) ;

   //-- Attributes
      
      stringVector vector ;
      size_t d0 ; 
      size_t d1 ;
             
} ;

#ifndef OUTLINE
#include <stringArray2D.icc>
#endif

#endif 
