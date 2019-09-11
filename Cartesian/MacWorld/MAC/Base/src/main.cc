#include <MAC.hh>
#include <MAC_Exec.hh>
#include <MAC_Module.hh>
#include <MAC_Root.hh>
#include <stringVector.hh>

//--------------------------------------------------------------------------
int main( int argc, char * argv[] )
//--------------------------------------------------------------------------
{
   stringVector args(0) ;
   int ok = MAC_Exec::initialize( argc, argv, args ) ;

   unsigned long long int result = MAC::used_memory() ;
   MAC::display_memory( MAC::out(), result ) ;
   
   if( ok == 0 )
   {
      MAC_Module const* m =
              MAC_Exec::create_module_with_data( MAC_Root::object(), args ) ;

      MAC_Application* app =
                 MAC_Exec::create_application( MAC_Root::object(), m, args ) ;
      MAC_Exec::run_application( app ) ;
   }
   
   result = MAC::used_memory() ;
   MAC::display_memory( MAC::out(), result ) ;   

   MAC_Exec::terminate() ;

   return( MAC_Exec::exit_code() ) ;
}






