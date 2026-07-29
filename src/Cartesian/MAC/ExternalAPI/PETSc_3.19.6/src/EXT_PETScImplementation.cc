#include <EXT_PETScImplementation.hh>

#include <MAC_assertions.hh>

#include <string>
#include <iostream>

//----------------------------------------------------------------------
EXT_PETScImplementation const*
EXT_PETScImplementation:: object( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScImplementation:: object" ) ;
   
   static EXT_PETScImplementation const* result =
                                          new EXT_PETScImplementation() ;
   return( result ) ;
}

//----------------------------------------------------------------------
EXT_PETScImplementation:: EXT_PETScImplementation( void )
//----------------------------------------------------------------------
  : LA_Implementation( "EXT_PETScImplementation" )
{
}

//----------------------------------------------------------------------
EXT_PETScImplementation:: ~EXT_PETScImplementation( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
EXT_PETScImplementation:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   std::string s( indent_width, ' ' ) ;
   os << s << "PETSc implementation" << std::endl ;
}

