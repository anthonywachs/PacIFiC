#include <mpi.h>
#include <computingtime.hh>


// ----------------------------------------------------------------------------
// Default constructor
ComputingTime::ComputingTime()
{}



  
// ----------------------------------------------------------------------------
// Constructor the application name as the input parameter 
ComputingTime::ComputingTime( string const& app_name_ )
{
  m_app_name = app_name_;
  m_total_elapsed_time = 0.;
  m_counter = 0;  
}    




// ----------------------------------------------------------------------------
// Copy constructor
ComputingTime::ComputingTime( ComputingTime const& CT )
{
  m_start = CT.m_start;
  m_end = CT.m_end;
  m_app_name = CT.m_app_name;
  m_total_elapsed_time = CT.m_total_elapsed_time;
  m_counter = CT.m_counter;
}




// ----------------------------------------------------------------------------
// Destructor
ComputingTime::~ComputingTime()
{}        




// ----------------------------------------------------------------------------
// Set start 
void ComputingTime::CT_set_start()
{
//  gettimeofday(&start,&tz);
  m_start = MPI_Wtime();    
}        




// ----------------------------------------------------------------------------
// Get the elapsed time 
double ComputingTime::CT_get_elapsed_time()
{
  double elapsed_time;
  
//  gettimeofday(&end,&tz);
  m_end = MPI_Wtime();    
//  elapsed_time=end.tv_sec-start.tv_sec+double(end.tv_usec-start.tv_usec)/1e6;
  elapsed_time = m_end - m_start;
  m_total_elapsed_time += elapsed_time;
  m_counter++;

  return( elapsed_time );  
}




// ----------------------------------------------------------------------------
// Add the elapsed time without incrementing the counter 
void ComputingTime::CT_add_elapsed_time()
{
  double elapsed_time;
  
//  gettimeofday(&end,&tz);  
  m_end = MPI_Wtime();    
//  elapsed_time=end.tv_sec-start.tv_sec+double(end.tv_usec-start.tv_usec)/1e6;
  elapsed_time = m_end - m_start;  
  m_total_elapsed_time += elapsed_time;  
}



    
// ----------------------------------------------------------------------------
// Get the application name 
string ComputingTime::CT_get_app_name() const
{
  return ( m_app_name );           
}




// ----------------------------------------------------------------------------
// Get the total elasped time 
double ComputingTime::CT_get_total_elapsed_time() const
{
  return ( m_total_elapsed_time );
}



    
// ----------------------------------------------------------------------------
// Get the mean elasped time 
double ComputingTime::CT_get_mean_elapsed_time() const
{
  return ( m_total_elapsed_time / max(m_counter,1) );
}



    
// ----------------------------------------------------------------------------
// Get the number of calls to this application  
int ComputingTime::CT_get_counter() const
{
  return ( m_counter );  
}




// ----------------------------------------------------------------------------
// Write elapsed time in seconds, minutes, hours and days
void write_elapsed_time_smhd( ostream& f, double const& elapsed_time,
    	string const& name )
{
  int days=int(elapsed_time/86400.),hours=int((elapsed_time-86400.*days)/3600.),
  	minutes=int((elapsed_time-86400.*days-3600.*hours)/60.);
  double seconds=elapsed_time-86400.*days-3600.*hours-60.*minutes;
  
  f << name << " (seconds) = " << elapsed_time << endl
	<< name << " = " << days << "d " << hours << "h " 
  	<< minutes << "m " << seconds << "s" << endl;
}
