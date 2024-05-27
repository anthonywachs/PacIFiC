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
    /** @brief Attempts to insert a particle in the simulation
    @param mode insertion order */
    virtual bool insertParticle( PullMode const& mode );

    /** @brief Sets particle initial positions from a file */
    virtual size_t setPositionParticlesFromFile();
  
    /** @brief Sets particle initial position with a structured array */
    virtual size_t setPositionParticlesArray();

    /** @brief Reads data for MPI simulations and creates and sets the MPI
    wrapper
    @param lx global domain size in the X direction
    @param ly global domain size in the Y direction
    @param lz global domain size in the Z direction  
    @param root XML node */
    virtual void readDomainDecomposition( DOMNode* root,
  	double const& lx, double const& ly, double const& lz );	
    
    /** @brief Sets the linked cell grid
    @param radius maximum circumscribed radius of particles 
    @param oshift empty string to shift the output */
    virtual void defineLinkedCell( double const& radius, string const& oshift );

    /** @brief Emergency termination in case of an issue */
    virtual void grainsAbort() const;
    
    /** @brief Returns the maximum particle ID number */
    virtual int getMaxParticleIDnumber() const;

    /** @brief Creates, inserts and links new particles in the simulation */
    virtual void InsertCreateNewParticles(); 
    
    /** @brief Displays the memory used by the simulation */
    virtual void display_used_memory() const; 
    
    /** @brief Synchronizes the PPWindow boolean relative to each sub-domain */
    virtual void synchronize_PPWindow();
    
    /** @brief Outputs timer summary */
    virtual void display_timer_summary();
    
    /** @brief Returns a particle class among the classes of new particles to
    insert
    @param mode particle insertion order
    @param ParticleClassesForCreation classes of new particles to insert
    @param random_local true if random is on the process, otherwise random over
    all processes */
    Particle* getParticleClassForCreation( PullMode const& mode,
  	list< pair<Particle*,size_t> >& ParticleClassesForCreation,
	bool const& random_local );
	
    /** @brief Returns the number of insertion positions */
    virtual size_t getNbInsertionPositions() const;
    
    /** @brief Checks the clones when a simulation is reloaded */
    virtual void checkClonesReload();      	                      
    //@}
  
};

#endif
  
