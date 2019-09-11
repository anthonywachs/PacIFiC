#ifndef SIZE_T_ARRAY_3D_HH
#define SIZE_T_ARRAY_3D_HH

#include <size_t_vector.hh>

/* 
sequences of values, all of type int, ordered according to three indices
in contiguous intervals
*/

class size_t_array3D
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      // Contruction of an instance whose all index intervals admit 0 as an
      // inclusive lower limit and respectively `dim0', `dim1' and `dim2'  
      // as an exclusive upper limit.
      size_t_array3D( size_t dim0, size_t dim1, size_t dim2,
                      size_t val=0 ) ;

      size_t_array3D( size_t_array3D const& other ) ;

      size_t_array3D const& operator=( size_t_array3D const& other ) ;

      // Reinitialize the internal state, as if the constructor with the
      // same argument list was just completed. 
      void re_initialize( size_t dim0, size_t dim1, size_t dim2,
                          size_t val=0 ) ;

   //-- Termination

     ~size_t_array3D( void ) ;

   //-- Comparison

      bool operator==( size_t_array3D const& other ) const ;
      bool operator!=( size_t_array3D const& other ) const ;
      
   //-- Access

      // exclusive upper limit of the `an_index'-th index interval
      size_t index_bound( size_t an_index ) const ;

      // item of indices `i0', 'i1' and `i2'
      size_t const& operator()( size_t i0, size_t i1, size_t i2 ) const ;

   //-- Element change

      // Increase the exclusive upper limit of the first index interval 
      // to `dim0'.
      void raise_first_index_bound( size_t dim0, size_t val=0 ) ;
  
      // Assign `val' to all items.
      void set( size_t val ) ; 

      // item of indices `i0', 'i1' and `i2'
      size_t& operator()( size_t i0, size_t i1, size_t i2 ) ;

   //-- Input - Output
     
      friend std::ostream& operator<<( std::ostream& out, 
                                       size_t_array3D const& a ) ;
      
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      size_t_array3D( void ) ;

   //-- Attributes
      
      size_t_vector vector ;
      size_t d0 ; 
      size_t d1 ;
      size_t d2 ;
             
} ;


#ifndef OUTLINE
#include <size_t_array3D.icc>
#endif

#endif 
