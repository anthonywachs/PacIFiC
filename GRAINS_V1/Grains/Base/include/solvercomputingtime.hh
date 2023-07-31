#ifndef _SOLVERCOMPUTINGTIME_HH_
#define _SOLVERCOMPUTINGTIME_HH_

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


/** @brief The class SolverComputingTime.

    Use for the definition of solver computing time (CT) measurements.

    @author A.Wachs - Particulate flow project 2003-2005
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================    
class SolverComputingTime
{
  public :
    /** @name Constructors */
    //@{  
    /** @brief Default constructor */
    SolverComputingTime();

    /** @brief Constructor with a list of the application names as input 
    parameter  
    @param all_apps_ list of application names */
    SolverComputingTime( list<string>& all_apps_ );    

    /** @brief Destructor */
    ~SolverComputingTime();        
    //@}    


    /** @name Methods */
    //@{
    /** @brief Insert a new application 
    @param app_name_ application name */
    void SCT_insert_app( string const& app_name_ );
    //@}
        

    /** @name Set methods */
    //@{
    /** @brief Set start 
    @param app_name_ application name */
    void SCT_set_start( string const& app_name_ );       
    //@}

    
    /** @name Get methods */
    //@{ 
    /** @brief Get the elapsed time 
    @param app_name_ application name */
    double SCT_get_elapsed_time( string const& app_name_ );
    
    /** @brief Add the elapsed time without incrementing the counter  
    @param app_name_ application name */
    void SCT_add_elapsed_time( string const& app_name_ );    
    
    /** @brief Get the total elapsed time 
    @param app_name_ application name */
    double SCT_get_total_elapsed_time( string const& app_name_ ) const;    
    
    /** @brief Get summary 
    @param f output stream 
    @param solver_total_time solver total computation time */
    void SCT_get_summary( ostream& f, double const& solver_total_time ) const;
    //@}  

    
  private :
    /** @name Parameters */
    //@{
    list<ComputingTime> m_all_apps; /**< list of all applications timed in this
    	solver */
    //@} 
    	

    /** @name Constructors */
    //@{
    /** @brief Copy constructor */
    SolverComputingTime( SolverComputingTime const& SCT );
    //@}              
};

#endif
