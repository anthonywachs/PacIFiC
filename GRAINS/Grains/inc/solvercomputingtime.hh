#ifndef _solvercomputingtime
#define _solvercomputingtime

#include <string>
#include <list>
#include <time.h>
#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <iomanip>
#include <computingtime.hh>
using namespace std;


/** @brief The Class SolverComputingTime.

Use for the definition of solver computation time (CT) measurements.

@author A.Wachs - Particulate flow project 2003-2005 */
class SolverComputingTime
{
  private :
    list<ComputingTime> all_apps; /**< list of all applications timed in this
    	solver */
	
  private :
    /** @name Constructors & Destructor */
    //@{
    /** @brief Copy constructor */
    SolverComputingTime(const ComputingTime &CT) {}
    //@}  

  public :
    /** @name Constructors & Destructor */
    //@{  
    /** @brief Constructor without argument */
    SolverComputingTime();

    /** @brief Constructor with argument the application name  
    @param all_apps_ list of application names */
    SolverComputingTime(list<string> &all_apps_);    

    /** @brief Destructor */
    ~SolverComputingTime();        
    //@}    

    /** @brief Insert a new application 
    @param app_name_ application name */
    void SCT_insert_app(const string &app_name_);

    /** @name SET methods */
    //@{
    /** @brief Set start 
    @param app_name_ application name */
    void SCT_set_start(const string &app_name_);       
    //@}
    
    /** @name GET methods */
    //@{ 
    /** @brief Get the elapsed time 
    @param app_name_ application name */
    double SCT_get_elapsed_time(const string &app_name_);
    
    /** @brief Add the elapsed time without incrementing the counter  
    @param app_name_ application name */
    void SCT_add_elapsed_time(const string &app_name_);    
    
    /** @brief Get the total elapsed time 
    @param app_name_ application name */
    double SCT_get_total_elapsed_time(const string &app_name_) const;    
    
    /** @brief Get summary 
    @param f output stream 
    @param solver_total_time solver total computation time */
    void SCT_get_summary(ostream &f,const double &solver_total_time) const;
    //@} 
           
};

#endif
