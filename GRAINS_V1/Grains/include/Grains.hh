#ifndef _GRAINS_HH_
#define _GRAINS_HH_

#include <mpi.h>
#include "GrainsMPIWrapper.hh"
#include "GrainsExec.hh"
#include "AllComponents.hh"
#include "AllInsertionWindows.hh"
#include "App.hh"
#include "solvercomputingtime.hh"
#include "computingtime.hh"
#include "LinkedCell.hh"
#include <list>
#include <vector>
#include <string>
#include <iostream>
using namespace std;

#include "ReaderXML.hh"


/** @brief Insertion mode */
enum InsertMode 
{
  IM_INITIALTIME, /**< Insertion at initial time */    
  IM_OVERTIME, /**< Insertion over time */    
  IM_NOINSERT /**< No insertion */
}; 


/** @brief Particle initial velocity */
enum InitialVelocity
{
  IV_ZERO, /**< initialisation to zero */    
  IV_CONSTANT, /**< initialisation to constant velocity */
  IV_RANDOM /**< initialisation to random velocity */
};


/** @brief Initial angular position of inserted particles */
enum InitialAngularPosition 
{
  IAP_FIXED, /**< fixed as defined in the particle class in the input file */
  IAP_RANDOM, /**< randomly assigned */
  IAP_FILE /**< fixed defined by an external file */
};


/** @brief Random generator seed */
enum RandomGeneratorSeed 
{
  RGS_DEFAULT, /**< initialized to default value (i.e., 1) */
  RGS_UDEF, /**< initialized to a user provided value */
  RGS_RANDOM /**< randomly initialized */
};


/** @brief The class Grains.

    Standard Grains3D application.

    @author Institut Francais du Petrole - 2003 - Creation 
    @author A.WACHS - 2021 - Major cleaning & refactoring */
// ============================================================================
class Grains : public ComputingTime, public SolverComputingTime
{
  public:
    /** @name Contructors & Destructor */
    //@{
    /** @brief Default constructor */
    Grains();

    /** @brief Destructor */
    virtual ~Grains();
    //@}


    /** @name Static methods */
    //@{  
    /** @brief Sets the time algorithm to predictor or corrector mode 
    @param predictor true if predictor */
    static void setMode( bool const& predictor ); 
  
    /** @brief Returns whether the time algorithm is in predictor mode */
    static bool isModePredictor();        
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
    
    
    /** @name I/O methods */
    //@{  
    /** @brief Writes an initial message in the standard output only
    on the process ranked 0 */
    virtual void initialOutputMessage();        
    //@}    

 
    /** @name Parameters */
    //@{
    static Vector3 m_vgravity; /**< gravity vector */	 
    //@} 
  
   
  protected:
    /** @name Parameters */
    //@{
    AllComponents m_allcomponents; /**< all components */
    list<App*> m_allApp; /**< list of applications, the first application is
    	always the linked cell collision detection application */  
    LinkedCell* m_collision; /**< explicit pointer to the linked cell collision 
    	detection application */ 
    double m_tstart; /**< initial simulation time */  
    double m_tend; /**< end simulation time */  
    double m_dt; /**< simulation time step */
    double m_time; /**< Physical time */   
    list<double> m_save; /**< list of restart and post-processing times */  
    list<double>::iterator m_timeSave; /**< iterator in the list m_save */
    ofstream fVitMax; /**< output file to store max and average component
    	velocity over the simulation */
    bool m_lastTime_save; /**< true if the last time was a time to write output
    	files */
    bool m_error_occured; /**< true if an error occured over the simulation */	
    string m_fileSave; /**< Root name of restart files */
    CloneInReload m_clonesInReloadFile; /**< clone (periodic, parallel or both)
    	writing mode */
    size_t m_dimension; /**< space dimension */
    bool m_periodic; /**< true if the domain is periodic in at least one
    	direction */
    vector<bool> m_periodicity; /**< vector of periodicity (3 booleans) */	
    bool m_allProcTiming; /**< whether all processes return data on the time
    	consumption in different parts of the code */
    bool m_restart; /**< is this simulation a continuation of a previous
    	simulation ? */	
    //@}


    /** @name Particle insertion parameters */
    //@{
    PullMode m_insertion_order; /**< Particle insertion order */  
    InsertMode m_insertion_mode; /**< Particle insertion mode */
    InitialVelocity m_initvit_mode; /**< Particle velocity initialization 
    	mode */ 
    InitialAngularPosition m_init_angpos; /**< Initial angular position of 
    	inserted particles */
    RandomGeneratorSeed m_randomseed; /**< Random generator seed */
    AllInsertionWindows m_insertion_windows; /**< Insertion windows */  
    string m_position; /**< External position file name or lattice type */
    string m_angular_position; /**< External angular position file name */    
    struct InsertionLattice* m_InsertionLattice; /**< Insertion lattice 
	features */   
    list< pair<Particle*,size_t> > m_newParticles; /**< types of new particles 
    	to be inserted */
    vector<Point3>* m_insertion_position; /**< list of insertion positions */
    list<Matrix>* m_insertion_angular_position; /**< list of insertion angular 
    	positions */    
    size_t m_i_sp; /**< index of the selected position in 
    	m_insertion_position */
    list<Matrix>::iterator m_il_sap; /**< iterator on the selected angular 
    	position in m_insertion_angular_position */	
    size_t m_insertion_frequency; /**< Insertion attempted every 
    	m_insertion_frequency time steps */
    bool m_force_insertion; /**< Force insertion even in case of contact with
    	other components */
    double m_RandomMotionCoefTrans; /**< maximum magnitude of translational 
    	random velocity */
    double m_RandomMotionCoefRot; /**< maximum magnitude of angular 
    	random velocity */		
    Vector3 m_InitVtrans; /**< Initial translational velocity */
    Vector3 m_InitVrot; /**< Initial angular velocity */
    size_t m_npwait_nm1; /**< number of particles not yet inserted at the
    	previous discrete time */
    //@}


    /** @name Fluid coupling parameters */
    //@{  
    static bool m_predictor_mode; /**< time algorithm mode */
    //@}


    /** @name Generic MPI parameters */
    //@{  
    int m_rank; /**< Rank of process in the MPI_COMM_activProc communicator
  	 (=0 in serial) */
    int m_nprocs; /**< Total number of processes in the MPI_COMM_activProc 
    	communicator (=1 in serial) */	 
    bool m_processorIsActive; /**< true if the process is active
  	(=true in serial) */
    GrainsMPIWrapper* m_wrapper; /**< manages MPI communications */		
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

    /** @brief Returns a point randomly selected in one of the insertion 
    windows or in the list of positions. Note: list of positions has priority
    over windows until it is empty */
    Point3 getInsertionPoint();
  
    /** @brief Attempts to insert a particle in the simulation
    @param mode insertion order */
    virtual bool insertParticle( PullMode const& mode );
  
    /** @brief Sets particle initial positions from a file */
    virtual size_t setPositionParticlesFromFile();

    /** @brief Sets angular particle initial positions from a file */
    virtual size_t setAngularPositionParticlesFromFile();    
  
    /** @brief Sets particle initial position with a structured array */
    virtual size_t setPositionParticlesStructuredArray();
    
    /** @brief Sets particle initial position through filling a cylinder with a
    regular array */
    virtual size_t setPositionParticlesCyl();    
  
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
  
    /** @brief Writes reload files 
    @param time physical time */
    virtual void saveReload( double const& time );	
  
    /** @brief Returns the maximum particle ID number */
    virtual int getMaxParticleIDnumber() const; 
  
    /** @brief Creates, inserts and links new particles in the simulation */
    virtual void InsertCreateNewParticles(); 
  
    /** @brief Deletes .result and .xml result files */
    void clearResultXmlFiles() const;
  
    /** @brief Displays the memory used by the simulation */
    virtual void display_used_memory() const;          

    /** @brief Synchronizes the PPWindow boolean relative to each sub-domain */
    virtual void synchronize_PPWindow();

    /** @brief Computes the initial velocity of a particle
    @param vtrans particle translational velocity
    @param vrot particle rotational velocity */
    void computeInitVit( Vector3& vtrans, Vector3& vrot ) const;
    
    /** @brief Computes particle forces and torques */
    void computeParticlesForceAndTorque();
	
    /** @brief Moves particles and obstacles
    @param dt_particle_vel time step to advance particle velocity 
    @param dt_particle_disp time step to advance particle position     
    @param dt_obstacle time step to advance obstacle velocity and position */
    void moveParticlesAndObstacles( double const& dt_particle_vel, 
    	double const& dt_particle_disp,
	double const& dt_obstacle );
	
    /** @brief Outputs timer summary */
    virtual void display_timer_summary();
    
    /** @brief Returns the number of insertion positions */
    virtual size_t getNbInsertionPositions() const;
    
    /** @brief Checks the periodic clones when a simulation is reloaded */
    virtual void checkClonesReload(); 
    
    /** @brief Reads time output frequency in MPI 
    @param root XML node */
    virtual void read_timeouputfrequency( DOMNode* root );            	
    //@}
  
};

#endif
  
