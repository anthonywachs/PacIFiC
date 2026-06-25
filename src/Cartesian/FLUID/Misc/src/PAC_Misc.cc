#include <PAC_Misc.hh>
#include <FV_DiscreteField.hh>
#include <FV_DomainBuilder.hh>
#include <MAC.hh>
#include <MAC_Error.hh>
#include <MAC_Utility.hh>
#include <cstdlib>
#include <sstream>

using std::ostringstream ;


//---------------------------------------------------------------------------
double 
PAC_Misc::compute_flow_rate( FV_DiscreteField const* FF,
      	size_t level, string const& boundary_name )
//--------------------------------------------------------------------------- 
{
   MAC_LABEL( "MAC_Utils:: compute_flow_rate" ) ;
   MAC_ASSERT( boundary_name == "left" || boundary_name == "right" ||
   	boundary_name == "top" || boundary_name == "bottom" || 
	boundary_name == "behind" || boundary_name == "front" ) ;
	
   size_t component = 0;
   double coef = 1. ;
   if ( boundary_name == "left" || boundary_name == "right" )
     component = 0 ;
   else if ( boundary_name == "top" || boundary_name == "bottom" )
     component = 1 ; 
   else component = 2 ;  
   
   if ( boundary_name == "left" || boundary_name == "bottom" ||
   	boundary_name == "behind" ) coef = -1. ;    
   
   double flow_rate = coef * FF->compute_boundary_cell_centered_DOF_integral( 
      	component, level, boundary_name ) ;

   return flow_rate ; 
   
}	




//---------------------------------------------------------------------------
bool 
PAC_Misc::is_main_boundary( string const& boundary_name, 
                            const size_t& dim, string const& message,
                            MAC_ModuleExplorer const* exp )
//--------------------------------------------------------------------------- 
{
   MAC_LABEL( "PAC_Misc:: is_main_boundary" ) ;

   bool result = false ;

   if( FV_DomainBuilder:: does_primary_color_exist( boundary_name ) )
     if( FV_DomainBuilder:: is_main_color( FV_DomainBuilder:: get_color_number( boundary_name ) ) )
     {
       if( dim == 3 ) result = true ;
       else if( boundary_name != "front" && boundary_name != "behind" )
         result = true ;
     }

   if ( exp )
     if ( !result )
     {
       string error_message = "   - left\n   - right\n   - bottom\n   - top\n";
       if ( dim == 3 )
         error_message += "   - behind\n   - front\n";
       MAC_Error::object()->raise_bad_data_value( exp, 
		message, error_message );     
     }
   
   return result ;

}	




//---------------------------------------------------------------------------
void 
PAC_Misc:: clearAllFiles( string const& resultsDirectory,
      	string const& savingsDirectory,
	const size_t& process_rank )
//---------------------------------------------------------------------------
{ 
   MAC_LABEL( "REG_UzawaNavierStokes:: clearAllFiles" ) ;
   
   if ( process_rank == 0 )
   {
     MAC::out() << "    CLEAR ALL FILES IN DIRECTORY " << resultsDirectory 
     	<< " AND " << savingsDirectory << endl << endl;

    //  char* path = getenv("FLUID_MAC_HOME");
    //  string exe(path);
    //  exe = exe+"/ExeUtils/Pacific_clear_exec "+resultsDirectory+" "+savingsDirectory;
    //  int sys_res = system( exe.c_str() ) ;


    if( pacific_clear_exec(resultsDirectory, savingsDirectory) != 0 )
     {
        ostringstream mesg ;
        mesg << "Unable to clear files in directories "
             << resultsDirectory << " and " << savingsDirectory << "."
             << endl ;
        MAC_Error::object()->raise_plain( mesg.str() ) ;
     }
   }   

}
