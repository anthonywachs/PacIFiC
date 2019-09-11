#include <LA_Implementation.hh>

#include <MAC_assertions.hh>
#include <MAC_ObjectRegister.hh>
#include <MAC_Root.hh>

using std::string ;

//----------------------------------------------------------------------------
LA_Implementation:: LA_Implementation( std::string const& name )
//----------------------------------------------------------------------------
   : MAC_Object( plugins_map() )
{
   MAC_LABEL( "LA_Implementation:: LA_Implementation" ) ;
   
   plugins_map()->register_item( name, this ) ;
}

//----------------------------------------------------------------------
LA_Implementation:: ~LA_Implementation( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_Implementation:: ~LA_Implementation" ) ;
}

//----------------------------------------------------------------------
MAC_ObjectRegister*
LA_Implementation:: plugins_map( void )
//----------------------------------------------------------------------
{
   static MAC_ObjectRegister* result =
         MAC_ObjectRegister::create( MAC_Root::object(),
                                     "LA_Implementation descendant" ) ;
   return( result ) ;
}
