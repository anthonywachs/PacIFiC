#include <LA_SeqImplementation.hh>

#include <MAC_assertions.hh>

#include <string>
#include <iostream>

//----------------------------------------------------------------------
LA_SeqImplementation const*
LA_SeqImplementation:: object( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_SeqImplementation:: object" ) ;
   
   static LA_SeqImplementation const* result = new LA_SeqImplementation() ;

   MAC_CHECK_POST( result != 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
LA_SeqImplementation:: LA_SeqImplementation( void )
//----------------------------------------------------------------------------
   : LA_Implementation( "LA_SeqImplementation" )
{
}

//----------------------------------------------------------------------
LA_SeqImplementation:: ~LA_SeqImplementation( void )
//----------------------------------------------------------------------
{
}

//----------------------------------------------------------------------
void
LA_SeqImplementation:: print( std::ostream& os, size_t indent_width ) const
//----------------------------------------------------------------------
{
   std::string s( indent_width, ' ' ) ;
   os << s << "MAC sequetial implementation" << std::endl ;
}
