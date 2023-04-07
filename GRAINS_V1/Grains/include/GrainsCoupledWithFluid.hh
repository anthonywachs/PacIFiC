#ifndef _GRAINSCOUPLEDWITHFLUID_HH_
#define _GRAINSCOUPLEDWITHFLUID_HH_

#include "Grains.hh"
#include "ReaderXML.hh"
#include "AppPRSHydroFT.hh"
#include <list>
#include <string>
using namespace std;


/** @brief The class GrainsCoupledWithFluid.

    Grains3D application coupled to a fluid flow to describe the dynamics
    of immersed rigid bodies.

    @author G.FERRER - Institut Francais du Petrole - 2003 - Creation 
    @author A.WACHS - 2021 - Major cleaning & refactoring */
// ============================================================================
class GrainsCoupledWithFluid : virtual public Grains
{
  public:
    /** @name Contructors & Destructor */
    //@{
    /** @brief Default constructor 
    @param fluid_density_ fluid density */
    GrainsCoupledWithFluid( double fluid_density_ );

    /** @brief Destructor */
    virtual ~GrainsCoupledWithFluid();
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


    /** @name Set methods */
    //@{
    /** @brief Sets the initial physical time
    @param time0 initial physical time */
    void setInitialTime( double const& time0 ) ;

    /** @brief Sets the initial postprocessing cycle number
    @param cycle0 initial postprocessing cycle number */    
    void setInitialCycleNumber( int const& cycle0 );

    /** @brief Sets the boolean m_forceReloadSame to true. This forces the code 
    to restart a simulation as a continuation of a previous simulation */
    void setReloadSame();
    
    /** @brief Sets the boolean Particle::setFluidCorrectedAcceleration. Default
    value is True, i.e., the particle acceleration is corrected by
    the factor ( 1 - fluid_density / particle_density )
    @param correct particle acceleration correction factor */
    void setFluidCorrectedAcceleration( bool correct );
    //@}  


    /** @name Methods */
    //@{
    /** @brief Writes features of moving rigid bodies in a stream to be used
    by the fluid flow solver
    @param is output stream */
    virtual void GrainsToFluid( istringstream &is ) const;
    
    /** @brief Updates particles velocity with data from the fluid solver
    @param velocity_data_array velocity data array
    @param b_set_velocity_nm1_and_diff updates the velocity at the previous time
    and the explicit velocity difference */
    virtual void updateParticlesVelocity( 
  	vector< vector<double> > const& velocity_data_array,
  	bool const& b_set_velocity_nm1_and_diff );
	
    /** @brief Updates particles hydro force and torque with data from the 
    fluid solver
    @param hydroft_data_array hydro force and torque data array */
    virtual void updateParticlesHydroFT( 
  	vector< vector<double> > const* hydroft_data_array );	   
    //@}
    
    
    /** @name I/O methods */
    //@{  
    /** @brief Writes an initial message in the standard output only
    on the process ranked 0 */
    virtual void initialOutputMessage(); 
    
    /** @brief Checks that the Paraview post-processing writer exists, otherwise
    creates it
    @param name_ files name
    @param root_ file root name
    @param isBinary whether to write in binary */
    void checkParaviewPostProcessing( string const& name_, string const& root_,
	bool const& isBinary );
    void checkParaviewPostProcessing( char const* name_, char const* root_,
  	bool const& isBinary ); 
	
    /** @brief Initially writes postprocessing files    
    @param indent_width output message indentation width */
    void InitialPostProcessing( size_t indent_width = 0 );

    /** @brief Writes postprocessing files    
    @param indent_width output message indentation width */
    void doPostProcessing( size_t indent_width = 0 );
    
    /** @brief Sets the Paraview post-processing translation vector in case of
    projection-translation
    @param tvx x coordinate
    @param tvy y coordinate
    @param tvz z coordinate */
    void setParaviewPostProcessingTranslationVector( 
      	double const& tvx, double const& tvy, double const& tvz );
    //@}
    

  protected:
    /** @name Parameters */
    //@{
    double m_fluid_density; /**< fluid density */
    double m_fluidflow_dt; /**< fluid flow simulation time step */ 
    bool m_forceReloadSame; /**< forces the code to restart a simulation as 
    	a continuation of a previous simulation */
    double m_min_dt; /**< minimum granular simulation time step */
    double m_max_dt; /**< minimum granular simulation time step */
    size_t m_ndt; /**< number of granular simulation time steps over 
    	a fluid flow simulation time step */ 
    AppPRSHydroFT* m_PRSHydroFT; /**< explicit pointer to the PRS hydro force
    	and torque application */	       	       
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
