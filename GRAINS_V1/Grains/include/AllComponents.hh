#ifndef _ALLCOMPONENTS_HH_
#define _ALLCOMPONENTS_HH_

#include "GrainsMPIWrapper.hh"
#include "Basic.hh"
#include "Error.hh"
#include <fstream>
#include <sstream>
#include <list>
#include <map>
using namespace std;
#include "Particle.hh"
#include "CompositeParticle.hh"
#include "Obstacle.hh"
#include "SimpleObstacle.hh"
#include "WriterXML.hh"


class App;
class AppCollision;
class RigidBody;
class ObstacleImposedVelocity;
class ObstacleImposedForce;
class GrainsMPIWrapper;
class PostProcessingWriter;
class AllInsertionWindows;


/** @brief Insertion order */
enum PullMode
{
  PM_ORDERED, /**< in the order of the input file */
  PM_RANDOM /**< in a random order */
};


/** @brief Insertion order */
enum CloneInReload
{
  CIR_NONE, /**< no clone in the reload file */
  CIR_NOPERIODIC, /**< no periodic clone in the reload file */
  CIR_ALL, /**< all clones in the reload file */  
};


/** @brief The class AllComponents.

    Manages all rigid components in the simulation.

    @author G.FERRER - Institut Francais du Petrole - 1999 - Creation
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class AllComponents
{
  public:
    /** @name Constructors & Destructor */
    //@{
    /** @brief Default constructor */
    AllComponents();

    /** @brief Destructor */
    ~AllComponents();
    //@}


    /** @name Methods */
    //@{
    /** @brief Updates particle activity */
    void UpdateParticleActivity();

    /** @brief Adds a reference particle
    @param particle the reference particle to be added 
    @param n number of such particles to be inserted in the simulation */
    void AddReferenceParticle( Particle *particle, size_t const &n );

    /** @brief Adds an obstacle
    @param obstacle_ the obstacle to be added */
    void AddObstacle( Obstacle* obstacle_ );

    /** @brief Associates the imposed velocity to the obstacle
    @param impvel imposed velocity */
    void LinkImposedMotion( ObstacleImposedVelocity* impvel );

    /** @brief Associates the imposed force to the obstacle
    @param load imposed force */
    void LinkImposedMotion( ObstacleImposedForce* load );

    /** @brief Initializes forces exerted on all components and set coordination
    number to 0
    @param time physical time
    @param dt time step magnitude
    @param withWeight with or without particle weight */
    void InitializeForces( double time, double dt, bool const& withWeight );

    /** @brief Initializes the transformation with crust of all components to not
    computed
    @param time physical time
    @param dt time step magnitude */
    void InitializeRBTransformWithCrustState( double time, double dt );

    /** @brief Computes all particles weight
    @param time physical time
    @param dt time step magnitude */
    void computeWeight( double time, double dt );

    /** @brief Moves all components
    @param time physical time
    @param dt_particle_vel velocity time step magnitude 
    @param dt_particle_disp motion time step magnitude
    @param dt_obstacle obstacle velocity and motion time step magnitude 
    @param LC linked cell grid */
    list<SimpleObstacle*> Move( double time,
	double const& dt_particle_vel, 
    	double const& dt_particle_disp,
	double const& dt_obstacle,
	LinkedCell const* LC );
    
    /** @brief Advances particles velocity over dt_particle_vel
    @param time physical time 
    @param dt_particle_vel velocity time step magnitude */
    void advanceParticlesVelocity( double time, 
    	double const& dt_particle_vel );        

    /** @brief Links all active particles to the application
    @param app application */
    void Link( AppCollision& app );

    /** @brief Writes components for Post-Processing at the start of the
    simulation
    @param time physical time
    @param dt time step magnitude
    @param LC linked cell grid
    @param insert_windows insertion windows
    @param rank process rank
    @param nprocs number of processes
    @param wrapper MPI wrapper */
    void PostProcessing_start( double time, double dt,
	LinkedCell const* LC, AllInsertionWindows const& insert_windows,
	int rank = 0, int nprocs = 1,
  	GrainsMPIWrapper const* wrapper = NULL );

    /** @brief Writes components for Post-Processing over the simulation
    @param time physical time
    @param dt time step magnitude
    @param LC linked cell grid
    @param insert_windows insertion windows    
    @param rank process rank
    @param nprocs number of processes
    @param wrapper MPI wrapper */
    void PostProcessing( double time, double dt,
	LinkedCell const* LC, AllInsertionWindows const& insert_windows,
	int rank = 0, int nprocs = 1,
	GrainsMPIWrapper const* wrapper = NULL );

    /** @brief Finalizes Post-Processing at the end of the simulation */
    void PostProcessing_end();

    /** @brief Writes components for Post-Processing in case of an error in
    contact or motion
    @param filename file root name
    @param errcomposants list of components involved in the error */
    void PostProcessingErreurComponents( string const& filename,
  	list<Component*> const& errcomposants );

    /** @brief Sets all active particles velocity to 0 if reset == "Reset"
    @param reset keyword to reset the particle velocity */
    void resetKinematics( string const& reset );

    /** @brief Adds the waiting particle to the list of active particles 
    @param parallel true if Grains runs in parallel mode */
    void WaitToActive( bool const& parallel = false );

    /** @brief Destroys the waiting particle */
    void DeleteAndDestroyWait();
  
    /** @brief Computes and returns the maximum and mean translational velocity
    of all components i.e. particles and obstacles
    @param vmax maximum translational velocity
    @param vmean mean translational velocity
    @param wrapper MPI wrapper */
    void ComputeMaxMeanVelocity( double& vmax, double& vmean,
  	GrainsMPIWrapper const* wrapper = NULL ) const;

    /** @brief Computes and writes in a file the minimum, maximum and mean
    translational velocity of all particles
    @param time physical time
    @param fileOut output file name
    @param rank process rank
    @param wrapper MPI wrapper */
    void monitorParticlesVelocity( double time, ofstream& fileOut,
  	int rank = 0, GrainsMPIWrapper const* wrapper = NULL ) const;

    /** @brief Updates obstacles' velocity without actually moving them
    @param time physical time
    @param dt time step magnitude */
    void setKinematicsObstacleWithoutMoving( double time, double dt );
    
    /** @brief Updates list of particles in parallel
    @param time physical time 
    @param newBufPart list of new buffer particles */
    void updateParticleLists( double time, list<Particle*>* newBufPart );    
    //@}


    /**@name Accessors */
    //@{
    /** @brief Returns an obstacle using its name as an input
    @param name obstacle name */
    Obstacle const* getObstacle( string const& name ) const;

    /** @brief Returns the root obstacle */
    Obstacle* getObstacles();

    /** @brief Returns particle with ID number id
    @param id ID number */
    Particle* getParticle( int id );

    /** @brief Returns component with ID number id */
    Component* getComponent( int id );

    /** @brief Returns a pointer to the particle to be inserted
    @param mode insertion mode */
    Particle* getParticleToInsert( PullMode mode );
    
    /** @brief Returns a pointer to the particle to be inserted
    @param geomtype particle geometric type */
    Particle* getParticleToInsert( int const& geomtype );    

    /** @brief Returns a pointer to the list of active particles */
    list<Particle*>* getActiveParticles();

    /** @brief Returns a const pointer to the list of active particles */
    list<Particle*> const* getActiveParticles() const;

    /** @brief Returns a pointer to the list of removed particles */
    list<Particle*>* getRemovedParticles();

    /** @brief Returns a pointer to the list of particles in the buffer zone */
    list<Particle*>* getParticlesInBufferzone();
    
    /** @brief Returns a pointer to the list of particles in the buffer zone */
    list<Particle*> const* getParticlesInBufferzone() const;    

    /** @brief Returns a pointer to the list of clone particles */
    list<Particle*>* getCloneParticles();

    /** @brief Returns a pointer to the map of serial clone particles */
    multimap<int,Particle*>* getPeriodicCloneParticles();
    
    /** @brief Returns a pointer to the map of serial clone particles */
    multimap<int,Particle*> const* getPeriodicCloneParticles() const;    

    /** @brief Returns a pointer to the vector of reference particles */
    vector<Particle*>* getReferenceParticles();

    /** @brief Returns a const pointer to the vector of reference particles */
    vector<Particle*> const* getReferenceParticles() const;

    /** @brief Returns the maximum particle circumscribed radius */
    double getCircumscribedRadiusMax();

    /** @brief Returns the minimum particle circumscribed radius */
    double getCircumscribedRadiusMin();

    /** @brief Returns the maximum particle crust thickness */
    double getCrustThicknessMax();

    /** @brief Returns the minimum particle crust thickness */
    double getCrustThicknessMin();

    /** @brief Returns the cumulative volume of all particles, both active and
    to be inserted */
    double getVolume() const;

    /** @brief Returns the cumulative volume of all actives particles */
    double getVolumeIn() const;

    /** @brief Returns the cumulative volume of all particles to be inserted */
    double getVolumeOut() const;

    /** @brief Returns the number of particles in this process */
    size_t getNumberParticles() const;
    
    /** @brief Returns the total number of particles in the
    system (i.e. on all subdomains/processes) */
    size_t getTotalNumberParticles() const;    

    /** @brief Returns the number of active particles in this process */
    size_t getNumberActiveParticles() const;
    
    /** @brief Returns the total number of active particles in the system (i.e. 
    on all subdomains/processes) */
    size_t getTotalNumberActiveParticles() const;    

    /** @brief Returns the number of active particles in this process with tag 
    0 ou 1 */
    size_t getNumberActiveParticlesOnProc() const;

    /** @brief Returns the total number of active particles with tag 
    0 ou 1 in the system (i.e. on all subdomains/processes) */
    size_t getNumberActiveParticlesOnAllProc() const;   
    
    /** @brief Returns the total number of particles in the physical system 
    (i.e. on all subdomains/processes), i.e. sum of total number of active 
    particles with tag 0 or 1 and particles to be inserted */
    size_t getTotalNumberPhysicalParticles() const; 
    
    /** @brief Returns the number of particles to insert in the physical 
    system */
    size_t getNumberPhysicalParticlesToInsert() const;        

    /** @brief Returns the highest particle ID number */
    int getMaxParticleIDnumber() const;

    /** @brief Returns the list of obstacles to send to the fluid */
    list<Obstacle*> getObstaclesToFluid() const;
    //@}


    /**@name Set methods */
    //@{
    /** @brief Computes and sets the numbers of particles in the system 
    @param wrapper MPI wrapper */
    void computeNumberParticles( GrainsMPIWrapper const* wrapper );    

    /** @brief Sets the frequency at which the relationship between obstacles
    and linked cell grid is updated
    @param updateFreq update frequency */
    void setObstaclesLinkedCellUpdateFrequency( int const& updateFreq );

    /** @brief Sets the parameters to output load exerted on obstacles
    @param root_ output directory name
    @param freq_ output frequency
    @param ObsNames liste de noms des obstacles */
    void setOutputObstaclesLoadParameters( string const& root_,
  	int const& freq_,
	list<string> const& ObsNames );

    /** @brief Sets a random translational and angular velocity to all particles
    @param coefTrans translational velocity amplitude
    @param coefRot angular velocity amplitude */
    void setRandomMotion( double const& coefTrans, double const& coefRot );

    /** @brief Set all contact map entries to false in all particles
    and all elementary obstacles */
    void setAllContactMapToFalse();

    /** @brief Update all contact map entries in all particles
    and all elementary obstacles */
    void updateAllContactMaps();
    
    /** @brief Set all contact map cumulative features to zero in all particles
    and all elementary obstacles */
    void setAllContactMapCumulativeFeaturesToZero(); 
    
    /** @brief Sets time integration scheme in all particles using the macro 
    variable GrainsExec::m_TIScheme */
    void setTimeIntegrationScheme();     
    //@}


    /**@name Methods I/O */
    //@{
    /** @brief Reloads components from an input stream
    @param fileSave input stream
    @param filename file name corresponding to the input stream 
    @param wrapper MPI wrapper */
    void read_pre2024( istream& fileSave, string const& filename,
    	GrainsMPIWrapper const* wrapper );
	
    /** @brief Reloads reference particles and obstacles from an input stream
    @param fileSave input stream 
    @param known_positions the list of remaining insertion positions    
    @param rank process rank    
    @param nprocs number of processes */
    size_t read( istream& fileSave, list<Point3>* known_positions,
    	int const& rank, int const& nprocs );
    
    /** @brief Reloads particles from an input stream 
    @param filename file name 
    @param npart number of particles to be read
    @param LC linked cell grid 
    @param rank process rank
    @param nprocs number of processes     
    @param wrapper MPI wrapper */
    void read_particles( string const& filename, size_t const& npart,
    	LinkedCell const* LC, int const& rank,
  	int const& nprocs, GrainsMPIWrapper const* wrapper );     	   

    /** @brief Writes components to an output stream. Only active particles are
    written to the stream
    @param fileSave output stream
    @param filename file name corresponding to the output stream 
    @param known_positions the list of remaining insertion positions
    @param cir clone (periodic, parallel or both) writing mode
    @param LC linked cell grid 
    @param rank process rank
    @param nprocs number of processes     
    @param wrapper MPI wrapper */
    void write( ostream &fileSave, string const& filename, 
    	list<Point3> const* known_positions, CloneInReload cir, 
	LinkedCell const* LC, int const& rank,
  	int const& nprocs, GrainsMPIWrapper const* wrapper ) const;
    
    /** @brief Writes components to a single MPI File in parallel. Only active 
    particles are written to the stream
    @param fileSave output stream
    @param filename file name corresponding to the output stream 
    @param known_positions the list of remaining insertion positions
    @param cir clone (periodic, parallel or both) writing mode
    @param LC linked cell grid    
    @param wrapper MPI wrapper 
    @param periodic true if the domain is periodic */
    void write_singleMPIFile( ostream &fileSave, string const& filename,
    	list<Point3> const* known_positions, CloneInReload cir, 
	LinkedCell const* LC, GrainsMPIWrapper const* wrapper, 
	bool periodic ) const;

    /** @brief Output operator
    @param f output stream
    @param EC the all components object */
    friend ostream& operator << ( ostream& f, AllComponents const& EC );

    /** @brief Debugging method
    @param s message de debug */
    void debug( char* s );

    /** @brief Adds a post-processing writer
    @param ppw post-processing writer */
    void addPostProcessingWriter( PostProcessingWriter* ppw );

    /** @brief Sets the initial post-processing cycle number
    @param cycle0 initial cycle number */
    void setInitialCycleNumber( int const& cycle0 );

    /** @brief Checks that the Paraview post-processing writer exists, and if
    not creates it
    @param rank process rank
    @param nprocs number of processes
    @param name_ files name
    @param root_ root file name
    @param isBinary whether to write in binary mode */
    void checkParaviewPostProcessing( int const& rank,
  	int const& nprocs,
  	string const& name_,
	string const& root_,
  	const bool& isBinary );

    /** @brief Writes load on obstacles in a file
    @param time physical time
    @param dt time step magnitude
    @param enforceOutput force writing
    @param increaseCounterOnly increases the writing counter only
    @param rank process rank */
    void outputObstaclesLoad( double time, double dt,
    	bool enforceOutput, bool increaseCounterOnly,
	int rank );

    /** @brief Computes load on obstacles
    @param time physical time
    @param dt time step magnitude
    @param wrapper MPI wrapper */
    void computeObstaclesLoad( double time, double dt,
      	GrainsMPIWrapper const* wrapper = NULL );

    /** @brief Initialises output files to write loads on obstacles
    @param rank process rank
    @param coupledFluid whether the simulation is coupled to a fluid solver
    @param time physical time */
    void initialiseOutputObstaclesLoadFiles( int rank,
      	bool coupledFluid, double time );
    //@}


  private:
    /** @name Parameters */
    //@{
    vector<Particle*> m_ReferenceParticles; /**< reference particle for each
    	class of particles */
    vector<size_t> m_NbRemainingParticlesToInsert; /**< number of remaining 
    	particles in each class to be inserted */
    size_t m_nb_physical_particles_to_insert; /**< number of 
    	remaining particles to insert in the physical system */		
    Particle* m_wait; /**< Next particle to be inserted */
    list<Particle*> m_RemovedParticles; /**< Particles removed from the
  	simulation */	
    list<Particle*> m_ActiveParticles; /**< All active particles in the
  	simulation */
    list<Particle*> m_ParticlesInBufferzone; /**< Active particles in a
	buffer zone */
    list<Particle*> m_CloneParticles; /**< Active particles that are parallel
  	clones of other active particles located in another subdomain/process */
    multimap<int,Particle*> m_PeriodicCloneParticles; /**< Periodic clone
    	particles and their relation to their master particle through its ID
	number */
    size_t m_nb_particles; /**< number of particles in this process */
    size_t m_total_nb_particles; /**< total number of particles in the system
  	(i.e. on all subdomains/processes) */
    size_t m_nb_active_particles; /**< number of active particles in this 
    	process */
    size_t m_total_nb_active_particles; /**< total number of active particles 
    	in the system (i.e. on all subdomains/processes) */
    size_t m_nb_active_particles_on_proc; /**< number of active particles in 
    	this process with tag 0 or 1 */
    size_t m_total_nb_active_particles_on_all_procs; /**< total number of 
    	active particles with tag 0 or 1 in the system 
	(i.e. on all subdomains/processes) */	
    size_t m_total_nb_physical_particles; /**< total number of particles 
    	in the physical system (i.e. on all subdomains/processes), i.e. sum of 
	total number of active particles with tag 0 or 1 and physical
	particles to be inserted */			
    Obstacle *m_obstacle; /**< Root obstacle */
    list<PostProcessingWriter*> m_postProcessors; /**< list of
  	post-processors */
    list<ObstacleImposedVelocity*> m_AllImposedVelocitiesOnObstacles; /**< list
  	of all imposed velocities on obstacles */
    list<ObstacleImposedForce*> m_AllImposedForcesOnObstacles; /**< list of all
  	imposed forces on obstacles */
    list<Obstacle*> m_outputTorsorObstacles; /**< Obstacles for which exerted
  	force and torque are written in a file */
    string m_outputTorsorObstacles_dir; /**< directory name where to write force
  	and torque on obstacles files */
    int m_outputTorsorObstacles_counter; /**< counter for force and torque on
  	obstacles output  */
    int m_outputTorsorObstacles_frequency; /**< frequency of force and torque on
  	obstacles output */
    //@}
    
    
    /**@name Methods */
    //@{
    /** @brief Computes and sets the maximum particle ID number and minimum 
    obstacle ID number
    @param wrapper MPI wrapper */
    void setParticleMaxIDObstacleMinID( GrainsMPIWrapper const* wrapper );    
    //@}    
};


/** @name Useful tools */
//@{
/** @brief Delete the fist instance of a pointer to a particle in a list of
pointers to particle
@param pointerslist particle pointer list
@param value pointer to delete */
bool removeParticleFromList( list<Particle*>& pointerslist,
	Particle* value );

/** @brief Delete the fist instance of a pointer to a particle in a set of
pointers to particle
@param pointersSet particle pointer set
@param value pointer to delete */
bool removeParticleFromSet( set<Particle*>& pointersSet, Particle* value );

/** @brief Delete the fist instance of a pointer to a simple obstacle in a list
of pointers to simple obstacle
@param pointerslist simple obstaclr pointer list
@param value pointer to delete */
bool removeObstacleFromList( list<SimpleObstacle*>& pointerslist,
	SimpleObstacle* value );
//@}

#endif
