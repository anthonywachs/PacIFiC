#ifndef SIZE_T_ARRAY_2D_HH
#define SIZE_T_ARRAY_2D_HH

#include <size_t_vector.hh>

class intArray2D ;

/*
sequences of values, all of type size_t, ordered according to two indices
in contiguous intervals
*/

class size_t_array2D
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      // Contruction of an instance whose all index intervals admit 0 as an
      // inclusive lower limit and respectively `dim0' and `dim1' as an 
      // exclusive upper limit.
      size_t_array2D( size_t dim0, size_t dim1, size_t val=0 ) ;

      size_t_array2D( size_t_array2D const& other ) ;

      size_t_array2D const& operator=( size_t_array2D const& other ) ;

      // Reinitialize the internal state, as if the constructor with the
      // same argument list was just completed. 
      void re_initialize( size_t dim0, size_t dim1, size_t val=0 ) ;

   //-- Termination

     ~size_t_array2D( void ) ;

   //-- Comparison

      bool operator==( size_t_array2D const& other ) const ;
      bool operator!=( size_t_array2D const& other ) const ;
      
   //-- Access

      // exclusive upper limit of the `an_index'-th index interval
      size_t index_bound( size_t an_index ) const ;
  
      // Reinitialize `x' by copying all items of self whose `an_index'-th 
      // index is equal to `index-value'.
      void extract_section( size_t an_index,
                            size_t index_value,
                            size_t_vector& x ) const ;

      // item of indices `i0' and `i1'
      size_t const& operator()( size_t i0, size_t i1 ) const ;

   //-- Element change

      // Increase the exclusive upper limit of the first index interval 
      // to `dim0'.
      void raise_first_index_bound( size_t dim0, size_t val=0 ) ;

      // Assign `val' to all items.
      void set( size_t val ) ; 

      void set( intArray2D const& other ) ;

      // Assign `x' items to the items of self whose `an_index'-th index is 
      // equal to `index-value' and whose other index is given by `x'.
      void set_section( size_t an_index, 
                        size_t index_value, 
                        size_t_vector const& x ) ;

      // item of indices `i0' and `i1'
      size_t& operator()( size_t i0, size_t i1 ) ;

            
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      size_t_array2D( void ) ;

   //-- Attributes      

      size_t_vector vector ;
      size_t d0 ; 
      size_t d1 ;
             
} ;



#ifndef OUTLINE
#include <size_t_array2D.icc>
#endif

#endif
