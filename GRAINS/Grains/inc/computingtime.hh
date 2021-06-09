#ifndef _computingtime
#define _computingtime

#include <time.h>
#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
using namespace std;


/** @brief The Class ComputingTime.

Use for the definition of computation time (CT) measurements the precision of
which is the microsecond.

@author A.Wachs - Particulate flow project 2003-2005 */
class ComputingTime
{
  private :
//    struct timeval start; /* start time */
//    struct timeval end; /* end time */
//    struct timezone tz; /* time zone (obsolete, see sys/time.h) */
    double start;
    double end;
    string app_name; /**< application name */
    double total_elapsed_time; /**< total elapsed time */
    int counter; /**< number of calls to this application */ 
	
  private :
    /** @name Constructors & Destructor */
    //@{
    /** @brief Constructor without argument */
    ComputingTime() {}
    //@}  

  public :
    /** @name Constructors & Destructor */
    //@{  
    /** @brief Constructor with argument the application name  
    @param app_name_ application name */
    ComputingTime(const string &app_name_);    

    /** @brief Copy constructor */
    ComputingTime(const ComputingTime &CT);

    /** @brief Destructor */
    ~ComputingTime();        
    //@}    

    /** @brief Write elapsed time in seconds, minutes, hours and days
    @param f output stream
    @param elapsed_time elapsed time 
    @param name defines the elapsed time (ex : "Computation time") */
    friend void write_elapsed_time_smhd(ostream &f,const double &elapsed_time,
    	const string &name);
    
    /** @name SET methods */
    //@{
    /** @brief Set start */
    void CT_set_start();        
    //@}
    
    /** @name GET methods */
    //@{ 
    /** @brief Get the elapsed time */
    double CT_get_elapsed_time();
    
    /** @brief Add the elapsed time without incrementing the counter */
    void CT_add_elapsed_time();      
    
    /** @brief Get the application name */
    string CT_get_app_name() const;
    
    /** @brief Get the total elapsed time */
    double CT_get_total_elapsed_time() const;
    
    /** @brief Get the mean elasped time */
    double CT_get_mean_elapsed_time() const;
    
    /** @brief Get the number of calls to this application */    
    int CT_get_counter() const;            
    //@} 
           
};

/** @brief Write elapsed time in seconds, minutes, hours and days
@param f output stream
@param elapsed_time elapsed time 
@param name defines the elapsed time (ex : "Computation time") */
void write_elapsed_time_smhd(ostream &f,const double &elapsed_time,
    	const string &name);
	
#endif
