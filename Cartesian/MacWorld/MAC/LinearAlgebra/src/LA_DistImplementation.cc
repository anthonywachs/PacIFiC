#include <LA_DistImplementation.hh>

#include <MAC_assertions.hh>

#include <string>
#include <iostream>

//----------------------------------------------------------------------
LA_DistImplementation const*
LA_DistImplementation:: object( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_DistImplementation:: object" ) ;
   
   static LA_DistImplementation const* result = new LA_DistImplementation() ;

   MAC_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
LA_DistImplementation:: LA_DistImplementation( void )
//----------------------------------------------------------------------------
   : LA_Implementation( "LA_DistImplementation" )
{
}

//----------------------------------------------------------------------
LA_DistImplementation:: ~LA_DistImplementation( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
LA_DistImplementation:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   std::string s( indent_width, ' ' ) ;
   os << s << "MAC distributed implementation" << std::endl ;
}
