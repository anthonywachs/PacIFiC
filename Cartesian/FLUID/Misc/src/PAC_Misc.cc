#include <PAC_Misc.hh>
#include <FV_DiscreteField.hh>
#include <MAC.hh>
#include <cstdlib>


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

     char* path = getenv("PACIFIC_HOME");
     string exe(path);
     exe = exe+"/ExeUtils/Pacific_clear_exec "+resultsDirectory+" "
     	+savingsDirectory;
     system(exe.c_str()); 
   }   

}
