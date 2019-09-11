#ifndef LA_Implementation_HH
#define LA_Implementation_HH

#include <MAC_Object.hh>

#include <iosfwd>

class MAC_ObjectRegister ;

/*
Class of singletons referring to a particular implementation of matrices and vectors.
This class is used to verify that matrix-vector operations are consistents.
*/

class LA_Implementation : public MAC_Object
{
   public: //-----------------------------------------------------------
     
   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const = 0 ;
      
   protected: //--------------------------------------------------------

      virtual ~LA_Implementation( void ) ;

   //-- Plug in

      LA_Implementation( std::string const& name ) ;
      
   private: //----------------------------------------------------------

      LA_Implementation( void ) ;
      LA_Implementation( LA_Implementation const& other ) ; 
      LA_Implementation& operator=( LA_Implementation const& other ) ;

      static MAC_ObjectRegister* plugins_map( void ) ;

   //-- Attributes
      
} ;

#endif

