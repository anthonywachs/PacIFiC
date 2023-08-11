#ifndef EXT_PETSc_IMPLEMENTATION_HH
#define EXT_PETSc_IMPLEMENTATION_HH

#include <LA_Implementation.hh>

/*
Implementation of PETSc matrices.

PUBLISHED
*/

class EXT_PETScImplementation : public LA_Implementation
{
   public: //-----------------------------------------------------------
      
   //-- Instance delivery and initialization
      
      static EXT_PETScImplementation const* object( void ) ;
      
   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;
      
   protected: //--------------------------------------------------------
      
   private: //----------------------------------------------------------

     ~EXT_PETScImplementation( void ) ;
      EXT_PETScImplementation( void ) ;
      
      EXT_PETScImplementation( EXT_PETScImplementation const& other ) ; 
      EXT_PETScImplementation& operator=(
                               EXT_PETScImplementation const& other ) ;
} ;

#endif

