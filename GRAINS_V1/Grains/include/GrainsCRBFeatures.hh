#ifndef _GRAINSCRBFEATURES_HH_
#define _GRAINSCRBFEATURES_HH_

#include "Grains.hh"
#include "AllComponents.hh"
#include "App.hh"
#include "LinkedCell.hh"

#include <list>
#include <vector>
#include <string>
#include <iostream>
using namespace std;

#include "ReaderXML.hh"


/** @brief The class GrainsCRBFeatures.

    Application that computes the volume, center of mass coordinates and moment
    of inertia tensor of a composite rigid body.

    @author A.WACHS - Institut Francais du Petrole - 2010 - Creation
    @author A.WACHS - 2021 - Major cleaning & refactoring */
// ============================================================================
class GrainsCRBFeatures : public Grains
{
  public:
    /** @name Contructors & Destructor */
    //@{
    /** @brief Default constructor */
    GrainsCRBFeatures();

    /** @brief Destructor */
    virtual ~GrainsCRBFeatures();
    //@}


    /** @name Core methods */
    //@{
    /** @brief Tasks to perform before time-stepping 
    @param rootElement XML root */
    virtual void do_before_time_stepping( DOMElement* rootElement );

    /** @brief Tasks to perform after time-stepping */
    virtual void do_after_time_stepping();

    /** @brief Runs the simulation over the prescribed time interval
    @param predict time algorithm mode 
    @param isPredictorCorrector true if the time algorithm is
    	predictor-corrector
    @param explicit_added_mass true if added mass is treated explicitly 
    @param time_interval runs the simulation on a time interval different than 
    	m_tend - m_tstart */
    virtual void Simulation( double time_interval = 0., 
    	bool predict = true,
  	bool isPredictorCorrector = false,
	bool explicit_added_mass = false ) ;	
    //@}
    
    
    /** @name I/O methods */
    //@{  
    /** @brief Writes an initial message in the standard output only
    on the process ranked 0 */
    virtual void initialOutputMessage();        
    //@}       

  
  private:
    /** @name Parameters */
    //@{  
    string m_outputfilename; /**< output file name */
    int m_geomtypeID; /**< geometric type number of the composite particle */
    double m_eqsphrad; /**< equivalent sphere radius of the composite 
    	particle */
    double m_dimensionless_delta; /**< dimensionless grid cell size */
    size_t m_Nmax; /**< maximum number of cells of the grid in all 
    	directions */         
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
};

#endif
  
