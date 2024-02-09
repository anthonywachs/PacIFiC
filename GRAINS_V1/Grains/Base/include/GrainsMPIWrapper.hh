#ifndef _GRAINSMPIWRAPPER_HH_
#define _GRAINSMPIWRAPPER_HH_

#include <mpi.h>
#include "MPINeighbors.hh"
#include "solvercomputingtime.hh"
#include "computingtime.hh"
#include "Particle.hh"
#include "LinkedCell.hh"
#include "Point3.hh"
#include "Vector3.hh"
#include "Matrix.hh"
#include "AppCollision.hh"  
#include <fstream>
#include <list>
#include <map>
using namespace std;


/** @brief The class GrainsMPIWrapper.

    Manages MPI communications in the parallel version of Grains3D.

    @author A.WACHS - 2021 - Major cleaning & refactoring */
// ============================================================================
class GrainsMPIWrapper : public SolverComputingTime
{
  public:
    /** @name Constructors & Destructor */
    //@{
    /** @brief Constructor with domain decomposition and periodicity as input 
    parameters 
    @param NX number of processes (=subdomains) in the X direction
    @param NY number of processes (=subdomains) in the Y direction
    @param NZ number of processes (=subdomains) in the Z direction
    @param PERX periodicity in the X direction (1 if periodic, 0 otherwise) 
    @param PERY periodicity in the Y direction (1 if periodic, 0 otherwise) 
    @param PERZ periodicity in the Z direction (1 if periodic, 0 otherwise) 
    @param oshift empty string to shift the output */
    GrainsMPIWrapper( int NX, int NY, int NZ,
  	int PERX, int PERY, int PERZ, string const& oshift );
  
    /** @brief Destructor */
    ~GrainsMPIWrapper();
    //@}

  
    /** @name Methods Get */
    //@{
    /** @brief Returns the number of processes in one direction
    @param i direction */
    int get_nb_procs_direction( int i ) const;
  
    /** @brief Returns the number of processes in all directions */
    int const* get_nb_procs_direction() const;   
  
    /** @brief Returns the MPI cartesian coordinates of the process */
    int const* get_MPI_coordinates() const;
  
    /** @brief Returns the MPINeighbors of the process */
    MPINeighbors const* get_MPI_neighbors() const;
  
    /** @brief Returns the periodicity of the domain that is the same as the
    MPI periodicity of the MPI cartesian topology */
    int const* get_MPI_periodicity() const;  
  
    /** @brief Returns the total number of processes in the MPI_COMM_WORLD 
    communicator */
    int get_total_number_of_processes() const; 
  
    /** @brief Returns the total number of active processes in the
    MPI_COMM_activProc communicator */
    int get_total_number_of_active_processes() const;   
  
    /** @brief Returns the process rank in the MPI_COMM_WORLD communicator */
    int get_rank_world() const; 
  
    /** @brief Returns the process rank in the MPI_COMM_activProc 
    communicator */
    int get_rank_active() const;   
  
    /** @brief Returns whether the process is active */
    bool isActive() const;                  
    //@}  


    /** @name Methods */
    //@{
    /** @brief Creates and updates clone particles using a Send-Recv strategy
    with neighboring processes in the MPI cartesian topology
    @param time physical time
    @param particles list of active particles
    @param particlesHalozone list of active particles in the halo (buffer) zone 
    @param particlesClones list of active clone particles
    @param referenceParticles reference particles
    @param LC linked cell grid */
    void UpdateOrCreateClones_SendRecvLocal_GeoLoc( double time,  	
	list<Particle*>* particles,
  	list<Particle*> const* particlesHalozone,
  	list<Particle*>* particlesClones,
	vector<Particle*> const* referenceParticles,
	LinkedCell* LC );
	
    /** @brief Gathers all particles on the master process for post-processing
    purposes
    @param particles list of active particles
    @param pwait list of non-active particles 
    @param referenceParticles reference particles
    @param nb_total_particles total number of particles on all processes */	
    vector<Particle*>* GatherParticles_PostProcessing(
  	list<Particle*> const& particles,
	list<Particle*> const& pwait,
	vector<Particle*> const& referenceParticles,
	size_t const& nb_total_particles ) const;	

    /** @brief Gathers all particle velocity-position data on the master 
    process for post-processing purposes
    @param particles list of active particles
    @param nb_total_particles total number of particles on all processes */
    vector< vector<double> >* GatherPositionVelocity_PostProcessing(
  	list<Particle*> const& particles,
	size_t const& nb_total_particles ) const;	

    /** @brief Gathers the class of all particles on the master process 
    @param particles list of active particles
    @param nb_total_particles total number of particles on all processes */
    vector< vector<double> >* GatherParticlesClass_PostProcessing(
  	list<Particle*> const& particles,
	size_t const& nb_total_particles ) const;	

    /** @brief Gathers all periodic clone particles on the master 
    process for post-processing purposes
    @param periodicCloneParticles list of periodic clone particles
    @param referenceParticles reference particles */	
    list<Particle*>* GatherPeriodicClones_PostProcessing(
  	list<Particle*> const& periodicCloneParticles,
	vector<Particle*> const& referenceParticles ) const;	
		
    /** @brief Broadcasts an integer from the master to all processes within the
    MPI_COMM_activProc communicator
    @param i integer */
    int Broadcast_INT( int const& i ) const;

    /** @brief Broadcasts a double from the master to all processes within the 
    MPI_COMM_activProc communicator
    @param d double */
    double Broadcast_DOUBLE( double const& d ) const;
  
    /** @brief Broadcasts an unsigned integer from the master to all processes 
    within the MPI_COMM_activProc communicator
    @param i unsigned integer */
    size_t Broadcast_UNSIGNED_INT( size_t const& i ) const;  
  
    /** @brief Sums an integer from all processes on all processes 
    within the MPI_COMM_activProc communicator
    @param i integer */
    int sum_INT( int const& i ) const;
  
    /** @brief Sums an integer from all processes on all processes 
    within the MPI_COMM_activProc communicator
    @param i unsigned integer */
    size_t sum_UNSIGNED_INT( size_t const& i ) const;  
  
    /** @brief Sums an integer from all processes on the master process 
    within the MPI_COMM_activProc communicator
    @param i integer */
    int sum_INT_master( int const& i ) const; 
  
    /** @brief Sums an unsigned integer from all processes on the master process 
    within the MPI_COMM_activProc communicator
    @param i unsigned integer */
    size_t sum_UNSIGNED_INT_master( size_t const& i ) const;     

    /** @brief Performs a "logical and" operation on input boolean value from 
    all processes on all processes
    @param input boolean value on which we want to perform the operation */
    bool logical_and( bool const& input ) const;     

    /** @brief Sums a double from all processes on all processes 
    within the MPI_COMM_activProc communicator
    @param x double */
    double sum_DOUBLE( double const& x ) const;
  
    /** @brief Sums a double from all processes on the master process 
    within the MPI_COMM_activProc communicator
    @param x double */
    double sum_DOUBLE_master( double const& x ) const;     
  
    /** @brief Minimum of an integer from all processes on all processes 
    within the MPI_COMM_activProc communicator
    @param i integer */
    int min_INT( int const& i ) const;
  
    /** @brief Minimum of an unsigned integer from all processes on all 
    processes within the MPI_COMM_activProc communicator
    @param i unsigned integer */
    size_t min_UNSIGNED_INT( size_t const& i ) const;  
  
    /** @brief Maximum of an integer from all processes on all processes 
    within the MPI_COMM_activProc communicator
    @param i integer */
    int max_INT( int const& i ) const;
  
    /** @brief Maximum of an unsigned integer from all processes on all 
    processes within the MPI_COMM_activProc communicator
    @param i unsigned integer */
    size_t max_UNSIGNED_INT( size_t const& i ) const;   
  
    /** @brief Maximum of a double from all processes on all processes 
    within the MPI_COMM_activProc communicator
    @param x double */
    double max_DOUBLE( double const& x ) const;
  
    /** @brief Maximum of a double from all processes on the master process 
    within the MPI_COMM_activProc communicator
    @param x double */
    double max_DOUBLE_master( double const& x ) const;  
  
    /** @brief Minimum of a double from all processes on all processes 
    within the MPI_COMM_activProc communicator
    @param x double */
    double min_DOUBLE( double const& x ) const;
  
    /** @brief Minimum of a double from all processes on the master process
    within the MPI_COMM_activProc communicator
    @param x double */
    double min_DOUBLE_master( double const& x ) const;
  
    /** @brief AllGather of an unsigned integer from all processes on all 
    processes within the MPI_COMM_activProc communicator
    @param i unsigned integer */
    size_t* AllGather_UNSIGNED_INT( size_t const& i ) const;          
  
    /** @brief Broadcasts a 3D point from the master to all processes within 
    the MPI_COMM_activProc communicator
    @param pt 3D point */
    Point3 Broadcast_Point3( Point3 const& pt ) const;

    /** @brief Broadcasts a 3D vector from the master to all processes within 
    the MPI_COMM_activProc communicator
    @param v 3D vector */
    Vector3 Broadcast_Vector3( Vector3 const& v ) const;    

    /** @brief Broadcasts a 3D vector from the master to all processes within 
    the MPI_COMM_activProc communicator
    @param mat matrix */
    Matrix Broadcast_Matrix( Matrix const& mat ) const;   
  
    /** @brief Sets the local MPI communicators involving the MPI cartesian 
    neighbors */
    void setCommLocal(); 
  
    /** @brief Outputs timer summary */
    void timerSummary() const;
  
    /** @brief Outputs the MPI log string per process and reinitialize it to
    empty
    @param f output stream */
    void writeAndFlushMPIString( ostream &f );
  
    /** @brief MPI_Barrier for all active processes */
    void MPI_Barrier_ActivProc() const; 

    /** @brief Sets periodicity vectors
    @param lx length of the domain in the x direction
    @param ly length of the domain in the y direction
    @param lz length of the domain in the z direction */
    void setMPIperiodicVectors( double const& lx, double const& ly, 
  	double const& lz ); 
  
    /** @brief Shares contact features among all active processes
    @param overlap_max maximum overlap
    @param overlap_mean average overlap    
    @param time_overlapMax time of maximum overlap
    @param nbIterGJK_mean average number of iterations of GJK for convergence */
    void ContactsFeatures( double& overlap_max,
	double& overlap_mean,
	double& time_overlapMax,
	double& nbIterGJK_mean ) const;
	
    /** @brief Sums force & torque exerted on obstacles on the master process 
    @param allMyObs list of simple obstacles */
    void sumObstaclesLoad( list<SimpleObstacle*> const& allMyObs ) const;
  
    /** @brief Distributes the number of particles in each class and on each
    process in the case of the block structured insertion   
    @param newPart total number of particles per class
    @param newPartProc number of particles per class on this process
    @param npartproc total number of particles to insert on this process
    @param ntotalinsert total number of particles to insert on all processes */
    void distributeParticlesClassProc( 
  	list< pair<Particle*,int> > const& newPart,
	list< pair<Particle*,int> >& newPartProc,
	size_t const& npartproc,
	size_t const& ntotalinsert ) const; 	
    //@}  


    /** @name Test methods */
    //@{
    /** @brief AllGather of a 1D array of n integers within the 
    MPI_COMM_activProc communicator
    @param n number of elements of the 1D array */
    void test_AllGatherv_INT( int const& n ) const;
  
    /** @brief AllGather of a 1D array of n integers within the 
    m_commMPINeighbors local communicator
    @param n number of elements of the 1D array */
    void testCommLocal_AllGatherv_INT( int const& n ) const;  

    /** @brief AllGather of a 1D array of n doubles within the 
    MPI_COMM_activProc communicator
    @param n number of elements of the 1D array */
    void test_AllGatherv_DOUBLE( int const& n ) const;
  
    /** @brief AllGather of a 1D array of n doubles within the 
    m_commMPINeighbors local communicator
    @param n number of elements of the 1D array */
    void testCommLocal_AllGatherv_DOUBLE( int const& n ) const;  

    /** @brief Send-Recv of a 1D array of n doubles within the 
    m_commMPINeighbors local communicator
    @param n number of elements of the 1D array */
    void testCommLocal_SendRecv_DOUBLE( int const& n ) const;  

    /** @brief Gather of an integer within the m_commMPINeighbors local 
    communicator */
    void testCommLocal_Gather_INT() const;
    //@}  


    /** @name Methodes static */
    //@{
    /** @brief Adds a string to the MPI log string
    @param add string to be added to the MPI log string */
    static void addToMPIString( string const& add ); 
  
    /** @brief Returns the GeoPosition as a function of the relative
    position in the MPI cartesian topology
    @param i relative position in X direction (-1,0 ou 1) 
    @param j relative position in Y direction (-1,0 ou 1)   
    @param k relative position in Z direction (-1,0 ou 1) */
    static GeoPosition getGeoPosition( int i, int j, int k );
    //@}  


    /** @name I/O methods */
    //@{  
    /** @brief Writes the MPI wrapper features in a stream
    @param f output stream 
    @param oshift empty string to shift the output */
    void display( ostream& f, string const& oshift ) const;

    /** @brief Writes the memory consumption per process in a stream
    @param f output flux */
    void display_used_memory( ostream& f ) const;  
    //@}
    

  private:
    /** @name Parameters */
    //@{
    MPI_Group m_MPI_GROUP_activProc; /**< active process group */    
    MPI_Comm m_MPI_COMM_activProc; /**< active process communicator */
    int *m_coords; /**< coordinates in the MPI cartesian topology */
    int *m_dim; /**< number of processes in each direction of the MPI cartesian
	topology */
    int *m_period; /**< Periodicity in each direction */
    bool m_isMPIperiodic; /**< true if at least one direction is periodic */
    int m_rank; /**< rank in the MPI_COMM_activProc communicator */
    int m_rank_world; /**< rank in the MPI_COMM_WORLD communicator */
    int m_rank_masterWorld; /**< rank in the MPI_COMM_activProc communicator of
    the process that has rank 0 in the MPI_COMM_WORLD communicator */
    int m_nprocs; /**< number of active processes */
    int m_nprocs_world; /**< total number of processes */
    bool m_is_activ; /**< is this process active ? */  
    MPINeighbors *m_neighbors; /**< neighbors of the process in the MPI cartesian 
    	topology */
    MPI_Comm *m_commgrainsMPI_3D; /**< MPI cartesian communicator */
    vector<MPI_Group*> m_groupMPINeighbors; /**< local groups involving 
    	neighbors in the MPI cartesian topology */	
    vector<MPI_Comm*> m_commMPINeighbors; /**< local MPI communicators involving
     	neighbors in the MPI cartesian topology */
    vector<bool> m_isInCommMPINeighbors; /**< local MPI communicators the
    	process is part of */
    vector<int> m_masterGeoPos; /**< GeoPosition of the master process
    	in local MPI communicators the process is part of */
    int m_nprocs_localComm; /**< number of processes in the local MPI 
    	communicator */
    int m_rank_localComm; /**< rank in the local MPI communicator */	
    int *m_master_localComm; /**< rank of master process in local MPI 
    	communicators */
    static string *m_MPILogString; /**< MPI log string */
    static vector< vector<int> > m_particleHalozoneToNeighboringProcs; /**< 
  	relationship between the GeoPosition in a halo zone from which 
	data are sent and the GeoPosition of the neighboring processes 
	that receive the data */
    static vector<int> m_GeoLocReciprocity; /**< reciprocal correspondence of
  	GeoPosition (ex: GEOPOS_BEHIND -> GEOPOS_FRONT ) */
    vector<Vector3> m_MPIperiodes; /**< periodic vectors */
    multimap<int,Particle*> m_AccessToClones; /**< facilitates access to clone
    	particles via their ID number */
    //@}


    /** @name Constructors & Destructor */
    //@{
    /** @brief Constructor */
    GrainsMPIWrapper();
    //@}


    /** @name Methods  */
    //@{
    /** @brief Sets the relationship between the GeoPosition in a halo
    zone from which data are sent and the GeoPosition of the neighboring
    processes that receive the data */
    void setParticleHalozoneToNeighboringProcs(); 
  
    /** @brief Sets the GeoPosition of the master process in the local
    communicators involving neighbors only to which this process belongs to */
    void setMasterGeoLocInLocalComm(); 
  
    /** @brief Creates and updates clones with the data sent by the neighboring
    processes 
    @param time physical time
    @param recvsize number of particles received
    @param recvbuf_DOUBLE array of double containing the data received
    @param NB_DOUBLE_PART number of doubles per particle       
    @param particlesClones list of active clone particles
    @param particles list of active particles
    @param particlesHalozone list of active particles in the halo (buffer) zone
    @param referenceParticles vector of reference particles
    @param LC linked cell grid */
    void UpdateOrCreateClones( double time,
 	int const& recvsize, double const* recvbuf_DOUBLE,
	int const& NB_DOUBLE_PART, 
  	list<Particle*>* particlesClones,
	list<Particle*>* particles,
  	list<Particle*> const* particlesHalozone,
	vector<Particle*> const* referenceParticles,
	LinkedCell* LC );               
    //@}  
};

#endif
