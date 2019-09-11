#ifndef BOOL_ARRAY_3D_HH
#define BOOL_ARRAY_3D_HH

#include <boolVector.hh>

/*
sequences of values, all of type bool, ordered according to three indices
in contiguous intervals
*/

class boolArray3D
{

   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      // Contruction of an instance whose all index intervals admit 0 as an
      // inclusive lower limit and respectively `dim0', `dim1' and `dim2'  
      // as an exclusive upper limit.
      boolArray3D( size_t dim0, size_t dim1, size_t dim2,
                   bool val=false ) ;

      boolArray3D( boolArray3D const& other ) ;

      boolArray3D const& operator=( boolArray3D const& other ) ;

      // Reinitialize the internal state, as if the constructor with the
      // same argument list was just completed. 
      void re_initialize( size_t dim0, size_t dim1, size_t dim2,
                          bool val=false ) ;

   //-- Termination

     ~boolArray3D( void ) ;

   //-- Comparison

     bool operator==( boolArray3D const& other ) const ;
     
     bool operator!=( boolArray3D const& other ) const ;
      
   //-- Access

      // exclusive upper limit of the `an_index'-th index interval
      size_t index_bound( size_t an_index ) const ;

      // item of indices `i0', 'i1' and `i2'
      bool const& operator()( size_t i0, size_t i1, size_t i2 ) const ;

   //-- Element change

      // Increase the exclusive upper limit of the first index interval 
      // to `dim0'.
      void raise_first_index_bound( size_t dim0, bool val=false ) ;
  
      // Assign `val' to all items.
      void set( bool val ) ; 

      // item of indices `i0', 'i1' and `i2'
      bool& operator()( size_t i0, size_t i1, size_t i2 ) ;

   //-- Input - Output
     
      friend std::ostream& operator<<( std::ostream& out, 
                                       boolArray3D const& a ) ;
      
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      boolArray3D( void ) ;
      
   //-- Attributes
    
      boolVector vector ;
      size_t d0 ; 
      size_t d1 ;
      size_t d2 ;
             
} ;


#ifndef OUTLINE
#include <boolArray3D.icc>
#endif

#endif 
