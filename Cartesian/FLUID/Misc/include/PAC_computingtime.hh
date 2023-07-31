#ifndef PAC_COMPUTINGTIME_HH_
#define PAC_COMPUTINGTIME_HH_

#include <time.h>
#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
using namespace std;


/** @brief The class PAC_ComputingTime.

    Use for the definition of computing time (CT) measurements the precision of
    which is the microsecond.

    @author A.Wachs - Particulate flow project 2003-2005
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class PAC_ComputingTime
{
  public :
    /** @name Constructors */
    //@{  
    /** @brief Constructor the application name as the input parameter 
    @param app_name_ application name */
    PAC_ComputingTime( string const& app_name_ );    

    /** @brief Copy constructor */
    PAC_ComputingTime( PAC_ComputingTime  const& CT );

    /** @brief Destructor */
    ~PAC_ComputingTime();        
    //@}    


    /** @name Friend methods */
    //@
    /** @brief Write elapsed time in seconds, minutes, hours and days
    @param f output stream
    @param elapsed_time elapsed time 
    @param name defines the elapsed time (ex : "Computation time") */
    friend void write_elapsed_time_smhd( ostream& f, 
    	double const& elapsed_time,
    	string const& name );
    //@} 

    
    /** @name Set methods */
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
    
    
  private :
    /** @name Parameters */
    //@{     
    double m_start; /**< start of timing */
    double m_end; /**< end of timing */
    string m_app_name; /**< application name */
    double m_total_elapsed_time; /**< total elapsed time */
    int m_counter; /**< number of calls to this application */
    //@}     


    /** @name Constructors */
    //@{
    /** @brief Default constructor */
    PAC_ComputingTime();
    //@}     
           
};

/** @brief Write elapsed time in seconds, minutes, hours and days
@param f output stream
@param elapsed_time elapsed time 
@param name defines the elapsed time (ex : "Computation time") */
void write_elapsed_time_smhd( ostream& f, double const& elapsed_time,
    	string const& name );
	
#endif
