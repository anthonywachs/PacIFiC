#ifndef _GRAINSMPI_HH_
#define _GRAINSMPI_HH_

#include "Grains.hh"


/** @brief The class GrainsMPI.

    Standard Grains3D application running in parallel mode with domain
    decomposition.

    @author G.FERRER - Institut Francais du Petrole - 2003 - Creation 
    @author A.WACHS - 2021 - Major cleaning & refactoring 
    @author A.WACHS - 2024 - Reimplementation */
// ============================================================================
class GrainsMPI : virtual public Grains
{
  public:
    /** @name Contructors & Destructor */
    //@{
    /** @brief Default constructor */
    GrainsMPI();

    /** @brief Destructor */
    virtual ~GrainsMPI();
    //@}



    /** @name Core methods */
    //@{
    /** @brief Runs the simulation over the prescribed time interval
    @param time_interval runs the simulation on a time interval different than 
    	m_tend - m_tstart  */
    virtual void Simulation( double time_interval = 0. ) ;  	
    //@}
    
    
    /** @name I/O methods */
    //@{  
    /** @brief Writes an initial message in the standard output only
    on the process ranked 0 */
    virtual void initialOutputMessage();        
    //@}    
  
   
  protected:
    /**@name Methods */
    //@{
    /** @brief Reads data for MPI simulations and creates and sets the MPI
    wrapper
    @param lx global domain size in the X direction
    @param ly global domain size in the Y direction
    @param lz global domain size in the Z direction  
    @param root XML node */
    virtual void readDomainDecomposition( DOMNode* root,
  	double const& lx, double const& ly, double const& lz );	
	
    /** @brief Returns the full result file name
    @param rootname root file name */
    virtual string fullResultFileName( string const& rootname ) const;
    
    /** @brief Sets the linked cell grid
    @param radius maximum circumscribed radius of particles 
    @param oshift empty string to shift the output */
    virtual void defineLinkedCell( double const& radius, string const& oshift );

    /** @brief Emergency termination in case of an issue */
    virtual void grainsAbort() const;
    
    /** @brief Returns the maximum particle ID number */
    virtual int getMaxParticleIDnumber() const;
    
    /** @brief Displays the memory used by the simulation */
    virtual void display_used_memory() const; 
    
    /** @brief Synchronizes the PPWindow boolean relative to each sub-domain */
    virtual void synchronize_PPWindow();
    
    /** @brief Outputs timer summary */
    virtual void display_timer_summary();                  
    //@}
  
};

#endif
  
