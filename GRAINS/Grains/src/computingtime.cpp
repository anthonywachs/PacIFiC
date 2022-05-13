#include <mpi.h>
#include <computingtime.hh>

  
/* Constructor with argument the application name  
-------------------------------------------------*/
ComputingTime::ComputingTime(const string &app_name_)
{
  app_name=app_name_;
  total_elapsed_time=0;
  counter=0;
  
}    




/* Copy constructor 
-------------------*/
ComputingTime::ComputingTime(const ComputingTime &CT)
{
  start=CT.start;
  end=CT.end;
  app_name=CT.app_name;
  total_elapsed_time=CT.total_elapsed_time;
  counter=CT.counter;

}




/* Destructor 
-------------*/
ComputingTime::~ComputingTime()
{}        




/* Set start 
------------*/
void ComputingTime::CT_set_start()
{
//  gettimeofday(&start,&tz);
  start = MPI_Wtime();  
  
}        




/* Get the elapsed time 
-----------------------*/
double ComputingTime::CT_get_elapsed_time()
{
  double elapsed_time;
  
//  gettimeofday(&end,&tz);
  end = MPI_Wtime();    
//  elapsed_time=end.tv_sec-start.tv_sec+double(end.tv_usec-start.tv_usec)/1e6;
  elapsed_time = end - start;
  total_elapsed_time+=elapsed_time;
  counter++;

  return(elapsed_time);
  
}




/* Add the elapsed time without incrementing the counter 
--------------------------------------------------------*/
void ComputingTime::CT_add_elapsed_time()
{
  double elapsed_time;
  
//  gettimeofday(&end,&tz);  
  end = MPI_Wtime();    
//  elapsed_time=end.tv_sec-start.tv_sec+double(end.tv_usec-start.tv_usec)/1e6;
  elapsed_time = end - start;  
  total_elapsed_time+=elapsed_time;
  
}



    
/* Get the application name 
---------------------------*/
string ComputingTime::CT_get_app_name() const
{
  return app_name;
           
}




/* Get the total elasped time 
-----------------------------*/
double ComputingTime::CT_get_total_elapsed_time() const
{
  return total_elapsed_time;

}



    
/* Get the mean elasped time 
----------------------------*/
double ComputingTime::CT_get_mean_elapsed_time() const
{
  return(total_elapsed_time/max(counter,1));

}



    
/* Get the number of calls to this application 
----------------------------------------------*/    
int ComputingTime::CT_get_counter() const
{
  return counter;
  
}




/* Write elapsed time in seconds, minutes, hours and days
---------------------------------------------------------*/
void write_elapsed_time_smhd(ostream &f,const double &elapsed_time,
	const string &name)
{
  int days=int(elapsed_time/86400.),hours=int((elapsed_time-86400.*days)/3600.),
  	minutes=int((elapsed_time-86400.*days-3600.*hours)/60.);
  double seconds=elapsed_time-86400.*days-3600.*hours-60.*minutes;
  
  f << name << " (seconds) = " << elapsed_time << endl
	<< name << " = " << days << "d " << hours << "h " 
  	<< minutes << "m " << seconds << "s" << endl;

}
