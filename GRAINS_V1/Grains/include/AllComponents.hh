#ifndef _ALLCOMPONENTS_HH_
#define _ALLCOMPONENTS_HH_

#include "Basic.hh"
#include "Error.hh"
#include <fstream>
#include <sstream>
#include <list>
#include <map>
using namespace std;
#include "Particle.hh"
#include "CompositeParticle.hh"
#include "SpheroCylinder.hh"
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
struct Window;


/** @brief Insertion order */
enum PullMode
{
  PM_ORDERED, /**< in the order of the input file */
  PM_RANDOM /**< in a random order */
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

    /** @brief Adds a particle
    @param particle the particle to be added */
    void AddParticle( Particle* particle );

    /** @brief Adds a reference particle
    @param particle the reference particle to be added */
    void AddReferenceParticle( Particle *particle );

    /** @brief Adds an obstacle
    @param obstacle_ L'obstacle a ajouter. */
    void AddObstacle( Obstacle* obstacle_ );

    /** @brief Associates the imposed velocity to the obstacle
    @param impvel imposed velocity */
    void LinkImposedMotion( ObstacleImposedVelocity* impvel );

    /** @brief Associates the imposed velocity to the obstacle
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
    @param dt_particle_disp displacement time step magnitude
    @param dt_obstacle obstacle velocity and displacement time step magnitude */
    list<SimpleObstacle*> Move( double time,
	double const& dt_particle_vel, 
    	double const& dt_particle_disp,
	double const& dt_obstacle );
	
    /** @brief Computes particles acceleration
    @param time physical time */
    void computeParticlesAcceleration( double time );
    
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
	LinkedCell const* LC, vector<Window> const& insert_windows,
	int rank = 0,
	int nprocs = 1,
  	GrainsMPIWrapper const* wrapper = NULL );

    /** @brief Writes components for Post-Processing over the simulation
    @param time physical time
    @param dt time step magnitude
    @param LC linked cell grid
    @param rank process rank
    @param nprocs number of processes
    @param wrapper MPI wrapper */
    void PostProcessing( double time, double dt,
	LinkedCell const* LC,
	int rank = 0,
  	int nprocs = 1,
	GrainsMPIWrapper const* wrapper = NULL );

    /** @brief Finalizes Post-Processing at the end of the simulation */
    void PostProcessing_end();

    /** @brief Writes components for Post-Processing in case of an error in
    contact or displacement
    @param filename file root name
    @param errcomposants list of components involved in the error */
    void PostProcessingErreurComponents( string const& filename,
  	list<Component*> const& errcomposants );

    /** @brief Set all active particles velocity to 0 if reset == "Reset"
    @param reset keyword to reset the particle velocity */
    void resetKinematics( string const& reset );

    /** @brief Transfer the inactive particle waiting to be inserted to the list
    of active particles */
    void ShiftParticleOutIn();

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

    /** @brief Updates geolocalization of particles in the halozone */
    void updateGeoPositionParticlesHalozone();

    /** @brief Updates obstacles' velocity without actually moving them
    @param time physical time
    @param dt time step magnitude */
    void setKinematicsObstacleWithoutMoving( double time, double dt );
    //@}


    /**@name Methods Get */
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

    /** @brief Returns a particle from the list of inactive particles
    @param mode insertion mode
    @param wrapper MPI wrapper */
    Particle* getParticle( PullMode mode,
  	GrainsMPIWrapper const* wrapper = NULL );

    /** @brief Returns a pointer to the list of active particles */
    list<Particle*>* getActiveParticles();

    /** @brief Returns a const pointer to the list of active particles */
    list<Particle*> const* getActiveParticles() const;

    /** @brief Returns a pointer to the list of inactive particles */
    list<Particle*>* getInactiveParticles();

    /** @brief Returns a const pointer to the list of inactive particles */
    list<Particle*> const* getInactiveParticles() const;

    /** @brief Returns a pointer to the list of particles in the halozone */
    list<Particle*>* getParticlesInHalozone();

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
    inactive */
    double getVolume() const;

    /** @brief Returns the cumulative volume of all actives particles */
    double getVolumeIn() const;

    /** @brief Returns the cumulative volume of all inactives particles */
    double getVolumeOut() const;

    /** @brief Returns the total number of particles, both active and
    inactive */
    size_t getNumberParticles() const;

    /** @brief Returns the number of active particles */
    size_t getNumberActiveParticles() const;

    /** @brief Returns the number of active particles with tag 0 ou 1 */
    size_t getNumberActiveParticlesOnProc() const;

    /** @brief Returns the total number of particles on all processes */
    size_t getNumberParticlesOnAllProc() const;

    /** @brief Returns the number of inactive particles */
    size_t getNumberInactiveParticles() const;

    /** @brief Returns the highest particle ID number */
    int getMaxParticleIDnumber() const;

    /** @brief Returns the list of obstacles to send to the fluid */
    list<Obstacle*> getObstaclesToFluid() const;
    //@}


    /**@name Methods Set */
    //@{
    /** @brief Sets the total number of particles on all processes
    @param nb_ total number of particles on all processes */
    void setNumberParticlesOnAllProc( size_t const& nb_ );

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
    
    /** @brief Set all contact map entry features to zero in all particles
    and all elementary obstacles */
    void setAllContactMapFeaturesToZero();    
    //@}


    /**@name Methods I/O */
    //@{
    /** @brief Reloads components from an input stream
    @param fileSave input stream
    @param filename file name corresponding to the input stream */
    void read( istream& fileSave, string const& filename );

    /** @brief Writes components to an output stream
    @param fileSave output stream
    @param filename file name corresponding to the output stream */
    void write( ostream &fileSave, string const& filename ) const;

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
    @param rank process rank
    @param nprocs number of processes
    @param wrapper MPI wrapper */
    void outputObstaclesLoad( double time, double dt,
    	bool enforceOutput = false,
      	bool increaseCounterOnly = false,
      	int rank = 0, int nprocs = 1, GrainsMPIWrapper const* wrapper = NULL );

    /** @brief Initialises output files to write loads on obstacles
    @param rank process rank
    @param coupledFluid whether the simulation is coupled to a fluid solver
    @param time physical time */
    void initialiseOutputObstaclesLoadFiles( int rank = 0,
      	bool coupledFluid = false, double time = 0. );
    //@}


  private:
    /** @name Parameters */
    //@{
    vector<Particle*> m_ReferenceParticles; /**< reference particle for each
    	class of particles */
    list<Particle*> m_InactiveParticles; /**< Inactive particles, either deleted
  	from the simulation or waiting to be inserted */
    Particle* m_wait; /**< Particle from the inactive particle list waiting to
    	be inserted */
    list<Particle*> m_ActiveParticles; /**< All active particles in the
  	simulation */
    list<Particle*> m_ParticlesInHalozone; /**< Active particles in a
    	halozone (i.e.a buffer zone) */
    list<Particle*> m_CloneParticles; /**< Active particles that are parallel
  	clones of other active particles located in another subdomain/process */
    multimap<int,Particle*> m_PeriodicCloneParticles; /**< Periodic clone
    	particles and their relation to their master particle through its ID
	number */
    size_t m_total_nb_particles; /**< total number of particles in the system
  	(i.e. on all subdomains/processes) */
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
