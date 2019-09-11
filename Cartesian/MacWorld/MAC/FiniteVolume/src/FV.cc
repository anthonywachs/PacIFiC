#include <FV.hh>
#include <MAC_Exec.hh>
#include <MAC_Error.hh>
#include <MAC_assertions.hh>


//---------------------------------------------------------------------------
std::ostream&
FV:: out( void )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "FV:: out" ) ;
   
   return MAC::out() ;

}
