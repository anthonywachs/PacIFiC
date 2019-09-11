#ifndef LA_DIST_IMPLEMENTATION_HH
#define LA_DIST_IMPLEMENTATION_HH

#include <LA_Implementation.hh>

/*
Implementation of MAC distributed matrices.

PUBLISHED
*/

class LA_DistImplementation : public LA_Implementation
{
   public: //-----------------------------------------------------------
      
   //-- Instance delivery and initialization
      
      static LA_DistImplementation const* object( void ) ;
      
   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //--------------------------------------------------------
      
   private: //----------------------------------------------------------

     ~LA_DistImplementation( void ) ;
      LA_DistImplementation( void ) ;
      
      LA_DistImplementation( LA_DistImplementation const& other ) ; 
      LA_DistImplementation& operator=( LA_DistImplementation const& other ) ;
      
} ;

#endif

