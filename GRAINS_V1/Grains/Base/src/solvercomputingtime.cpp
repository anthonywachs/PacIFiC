#include <solvercomputingtime.hh>
#include <cstdlib>

  
// ----------------------------------------------------------------------------
// Default constructor
SolverComputingTime::SolverComputingTime()
{}




// ----------------------------------------------------------------------------
// Constructor with a list of the application names as input parameter 
SolverComputingTime::SolverComputingTime( list<string>& all_apps_ )
{
  list<string>::iterator il;
  for (il=all_apps_.begin();il!=all_apps_.end();il++)
  {
    ComputingTime cpt(*il);
    m_all_apps.push_back(cpt);
  }  
}    




// ----------------------------------------------------------------------------
// Copy constructor
SolverComputingTime::SolverComputingTime( ComputingTime const& CT ) 
{}




// ----------------------------------------------------------------------------
// Destructor 
SolverComputingTime::~SolverComputingTime()
{
  m_all_apps.clear();
} 




// ----------------------------------------------------------------------------
// Set start
void SolverComputingTime::SCT_set_start( string const& app_name_ )       
{
  list<ComputingTime>::iterator il=m_all_apps.begin();
  bool found=false;
  
  while( ( il != m_all_apps.end() ) && ( !found ) )
  {
    if ( il->CT_get_app_name() == app_name_ ) 
    {
      il->CT_set_start(); 
      found = true;
    }
    else il++;
  }
  
  if ( !found )
  { 
    cout << "Wrong application name \"" << app_name_ 
    	<< "\" in SolverComputingTime::set_start" << endl;
    cout << "Available applications are :" << endl;
    for (il=m_all_apps.begin();il!=m_all_apps.end();il++)
      cout << "  - " << il->CT_get_app_name() << endl;
    exit(0);
  } 
}


    
 
// ----------------------------------------------------------------------------
// Get the elapsed time 
double SolverComputingTime::SCT_get_elapsed_time( string const& app_name_ )
{
  list<ComputingTime>::iterator il=m_all_apps.begin();
  bool found=false;
  double elapsed_time=0.;
  
  while( ( il != m_all_apps.end() ) && ( !found ) )
  {
    if ( il->CT_get_app_name() == app_name_ ) 
    {
      elapsed_time=il->CT_get_elapsed_time(); 
      found = true;
    }
    else il++;
  }
  
  if ( !found )
  { 
    cout << "Wrong application name in SolverComputingTime::"
    	<< "get_elapsed_time" << endl;
    exit(0);
  }
  
  return ( elapsed_time ); 
}




// ----------------------------------------------------------------------------
// Add the elapsed time without incrementing the counter
void SolverComputingTime::SCT_add_elapsed_time( string const& app_name_ )
{
  list<ComputingTime>::iterator il=m_all_apps.begin();
  bool found=false;
  
  while( ( il != m_all_apps.end() ) && ( !found ) )
  {
    if ( il->CT_get_app_name() == app_name_ ) 
    {
      il->CT_add_elapsed_time(); 
      found = true;
    }
    else il++;
  }
  
  if ( !found )
  { 
    cout << "Wrong application name in SolverComputingTime::"
    	<< "get_elapsed_time" << endl;
    exit(0);
  }
}




// ----------------------------------------------------------------------------
// Get the total elapsed time 
double SolverComputingTime::SCT_get_total_elapsed_time( string const& app_name_ )
	const
{
  list<ComputingTime>::const_iterator il=m_all_apps.begin();
  bool found=false;
  double total_elapsed_time=0.;
  
  while( ( il != m_all_apps.end() ) && ( !found ) )
  {
    if ( il->CT_get_app_name() == app_name_ ) 
    {
      total_elapsed_time=il->CT_get_total_elapsed_time(); 
      found = true;
    }
    else il++;
  }
  
  if ( !found )
  { 
    cout << "Wrong application name in SolverComputingTime::"
    	<< "get_elapsed_time" << endl;
    exit(0);
  }
  
  return ( total_elapsed_time ); 
}



    
// ----------------------------------------------------------------------------
// Get summary 
void SolverComputingTime::SCT_get_summary( ostream& f,
	double const& solver_total_time ) const
{
  list<ComputingTime>::const_iterator il;
  double total_time_apps=0,time_other;
  
  if ( !m_all_apps.empty() )
  {
    f.precision(3);
    f << "----------------------------------------------------------" << 
    	"--------------------" << endl;
    f << "                          Name          CT      Number "
    	<< "    Mean CT        % of" << endl;
    f << "                                              of calls   "
    	<< "             total CT" << endl;  	 
    for (il=m_all_apps.begin();il!=m_all_apps.end();il++)
    {
      f << setw(30) << il->CT_get_app_name() 
      	<< setw(12) << il->CT_get_total_elapsed_time() 
	<< setw(12) << il->CT_get_counter()   
	<< setw(12) << il->CT_get_mean_elapsed_time() 
	<< setw(12) << 100*il->CT_get_total_elapsed_time()/solver_total_time 
	<< endl;
      total_time_apps+=il->CT_get_total_elapsed_time();	
    }
    time_other=solver_total_time-total_time_apps;
    f << setw(30) << "Other (outputs, ...)" 
      	<< setw(12) << time_other
	<< setw(12) << " "  
	<< setw(12) << " "
	<< setw(12) << 100*time_other/solver_total_time 
	<< endl;    
    f << "----------------------------------------------------------" << 
    	"--------------------" << endl;
    f << setw(30) << "Total" << setw(12) << solver_total_time << endl;
    f << "----------------------------------------------------------" << 
    	"--------------------" << endl;    	
  }     
}         




// ----------------------------------------------------------------------------
// Insert a new application 
void SolverComputingTime::SCT_insert_app( string const& app_name_ )
{
  list<ComputingTime>::iterator il=m_all_apps.begin();
  bool found=false;
  
  while( ( il != m_all_apps.end() ) && ( !found ) )
  {
    if ( il->CT_get_app_name() == app_name_ ) found = true;
    else il++;
  }
  
  if ( !found )
  { 
    ComputingTime cpt(app_name_);
    m_all_apps.push_back(cpt);
  }
}
