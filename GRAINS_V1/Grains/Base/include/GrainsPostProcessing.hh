#ifndef _GRAINSPOSTPROCESSING_HH_
#define _GRAINSPOSTPROCESSING_HH_

#include "Grains.hh"
#include "ReaderXML.hh"
#include <list>
#include <string>
using namespace std;


/** @brief Global porosity */
struct GlobalPorosity
{
  struct Window domain; /**< porosity domain */
  size_t nintervals[3]; /**< Number of intervals in each direction */   
};


/** @brief The class GrainsPostProcessing.

    Grains3D application to postprocess the granular medium at a given time
    through reloading the system as in a restarted simulation.

    @author A.WACHS - 2022 - Major cleaning & refactoring */
// ============================================================================
class GrainsPostProcessing : virtual public Grains
{
  public:
    /** @name Contructors & Destructor */
    //@{
    /** @brief Default constructor 
    @param fluid_density_ fluid density */
    GrainsPostProcessing();

    /** @brief Destructor */
    virtual ~GrainsPostProcessing();
    //@}


    /** @name Core methods */
    //@{
    /** @brief Tasks to perform before time-stepping 
    @param rootElement XML root */
    virtual void do_before_time_stepping( DOMElement* rootElement );

    /** @brief Tasks to perform after time-stepping */
    virtual void do_after_time_stepping();

    /** @brief Runs the simulation over the prescribed time interval
    @param time_interval runs the simulation on a time interval different than 
    	m_tend - m_tstart  */
    virtual void Simulation( double time_interval = 0. ) ; 	
    //@}


    /** @name Methods */
    //@{  
    //@}
    
    
    /** @name I/O methods */
    //@{  
    /** @brief Writes an initial message in the standard output only
    on the process ranked 0 */
    virtual void initialOutputMessage(); 	          
    //@}
    

  protected:
    /** @name Parameters */
    //@{ 
    struct GlobalPorosity* m_global_porosity;    	       
    //@}


    /**@name Methods */
    //@{
    /** @brief Construction of the simulation: linked cell, particles &
    obstacles, domain decomposition 
    @param rootElement XML root */
    virtual void Construction( DOMElement* rootElement );

    /** @brief Additional features of the simulation: time features, insertion,
    post-processing
    @param rootElement XML root */
    virtual void AdditionalFeatures( DOMElement* rootElement );

    /** @brief External force definition
    @param rootElement XML root */
    virtual void Forces( DOMElement* rootElement );
    //@}    
    
    
  private:
};

#endif
