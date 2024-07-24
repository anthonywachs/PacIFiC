#ifndef _LINKEDCELL_HH_
#define _LINKEDCELL_HH_

#include "AppCollision.hh"
#include "BBox.hh"
#include "RigidBodyWithCrust.hh"

class Cell;
class MPINeighbors;
class GrainsMPIWrapper;

#include <vector>
using namespace std;


/** @brief The class LinkedCell.

    Primary application to detect collisions and compute collision force
    & torques on rigid bodies. Broad phase detection is based on a linked-cell
    grid.

    @author G.FERRER - Institut Francais du Petrole - 2003 - Creation
    @author A.WACHS - 2021 - Major cleaning & refactoring */
// ============================================================================
class LinkedCell : public AppCollision
{
  public:
    /** @name Contructors & Destructor */
    //@{
    /** @brief Default constructor */
    LinkedCell();

    /** @brief Destructor */
    ~LinkedCell();
    //@}


    /** @name Methods */
    //@{
    /** @brief Computes forces and torques exerted on rigid bodies
    @param time physical time
    @param dt time step magnitude
    @param particles active particles */
    void ComputeForces( double time, double dt,
  	list<Particle*> const* particles );

    /** @brief Returns whether a particle is in contact with another component
    using the method Component::isContact
    @param particle particle */
    bool isContact( Particle const* particle ) const;

    /** @brief Returns whether a particle is in contact with another component
    using the method Component::isContactWithCrust
    @param particle particle */
    bool isContactWithCrust( Particle const* particle ) const;

    /** @brief Returns whether a particle is close to another component
    using the method Component::isClose
    @param particle particle */
    bool isClose( Particle const* particle ) const ;

    /** @brief Returns whether a particle is close to another component
    using the method Component::isCloseWithCrust
    @param particle particle */
    bool isCloseWithCrust( Particle const* particle ) const;

    /** @brief Links a particle with the linked cell grid without checking if
    the particle overlaps with another rigid body
    @param particle particle */
    void Link( Particle* particle );

    /** @brief Links the root obstacle with the linked cell grid at the start
    of the simulation
    @param root_obstacle root obstacle */
    void Link( Obstacle* root_obstacle );

    /** @brief Updates links between particles & obstacles and the linked cell
    grid
    @param time physical time
    @param dt time step magnitude
    @param particles active particles */
    void LinkUpdate( double time, double dt,
  	list<Particle*>* particles );

    /** @brief Updates the link of an active particle and the linked cell grid
    @param particle particle */
    void LinkUpdateActiveParticle( Particle* particle );

    /** @brief Output operator
    @param f output stream
    @param LC LinkedCell object */
    friend ostream& operator << ( ostream& f, LinkedCell const& LC );

    /** @brief Returns whether a point lies inside the linked cell grid
    @param position point coordinates */
    bool isInLinkedCell( Point3 const& position ) const;

    /** @brief Returns whether a point lies inside the linked cell grid
    @param gx x-ccordinate of the point
    @param gy y-ccordinate of the point
    @param gz z-ccordinate of the point */
    bool isInLinkedCell( double const& gx, double const& gy,
	double const& gz ) const;

    /** @brief Removes a particle from the linked cell grid
    @param particle particle to be deleted */
    void remove( Particle* particle );

    /** @brief Removes an obstacle from the linked cell grid
    @param obs obstacle to be deleted */
    void remove( SimpleObstacle* obs );

    /** @brief Removes clone particles that exited the local linked cell grid
    @param time physical time
    @param particles list of active particles
    @param particlesClones list of active clone particles
    @param wrapper MPI wrapper */
    void DestroyOutOfDomainClones( double time,
	list<Particle*>* particles,
	list<Particle*>* particlesClones,
	GrainsMPIWrapper const* wrapper = NULL );

    /** @brief Sets the list of all neighboring cells */
    void setCellCompleteNeighborhood();

    /** @brief Attempts to insert a particle in serial mode
    @param particle particle
    @param particles list of active particles
    @param particlesPeriodicClones map of active periodic clone particles
    @param ReferenceParticles vector of reference particles
    @param periodic true if the domain is at least mono-periodic
    @param force_insertion force insertion regardless of potential contacts with
    other particles or obstacles */
    bool insertParticleSerial( Particle* particle, list<Particle*>* particles,
    	multimap<int,Particle*>* particlesPeriodicClones,
	vector<Particle*> const* ReferenceParticles,
	bool const& periodic,
    	bool const& force_insertion );

    /** @brief Updates periodic clones and destroys those out of the linked cell
    grid in serial mode
    @param particles list of active particles
    @param particlesPeriodicClones map of active periodic clone particles 
    @param destroyOnly performs periodic clones destruction only when true (i.e.
    does not update periodic clones in the linked cell grid) */
    void updateDestroyPeriodicClones( list<Particle*>* particles,
    	multimap<int,Particle*>* particlesPeriodicClones,
	bool destroyOnly = false );

    /** @brief Creates/destroys periodic clones after LinkUpdate in serial mode
    @param particles list of active particles
    @param particlesPeriodicClones map of active periodic clone particles
    @param ReferenceParticles vector of reference particles */
    void createDestroyPeriodicClones( list<Particle*>* particles,
    	multimap<int,Particle*>* particlesPeriodicClones,
	vector<Particle*> const* ReferenceParticles );
	
    /** @brief Checks periodic clones in serial mode when a simulation is
    reloaded
    @param particles list of active particles
    @param particlesPeriodicClones map of active periodic clone particles
    @param ReferenceParticles vector of reference particles 
    @param time physical time */
    void checkPeriodicClonesReload( list<Particle*>* particles,
    	multimap<int,Particle*>* particlesPeriodicClones,
	vector<Particle*> const* ReferenceParticles, 
	double const& time );	
	
    /** @brief Attempts to insert a particle in parallel mode
    @param time physical time    
    @param particle particle
    @param particles list of active particles
    @param particlesClones list of active clone particles 
    @param ReferenceParticles vector of reference particles       
    @param periodic true if the domain is at least mono-periodic
    @param force_insertion force insertion regardless of potential contacts with
    other particles or obstacles 
    @param wrapper MPI wrapper */
    pair<bool,bool> insertParticleParallel( double time,
    	Particle* particle, 
    	list<Particle*>* particles,
	list<Particle*>* particlesClones,
	vector<Particle*> const* ReferenceParticles,
    	bool const& periodic, bool const& force_insertion,
	GrainsMPIWrapper const* wrapper = NULL );
	
    /** @brief Returns an array of point coordinates of the local grid in a 
    direction
    @param dir direction */
    vector<double> local_coordinates( size_t const& dir ) const;
    
    /** @brief Returns an array of point coordinates of the global grid in a 
    direction
    @param dir direction */
    vector<double> global_coordinates( size_t const& dir ) const; 
    
    /** @brief Checks that none of the structured array positions is exactly 
    at a limit of the linked cell grid, otherwise shift by 1e-12 
    @param InsertionArray structured array positions 
    @param wrapper MPI wrapper */
    void checkStructuredArrayPositionsMPI( struct StructArrayInsertion* 
    	InsertionArray, GrainsMPIWrapper const* wrapper ) const;   		
    //@}


    /**@name Methods Get */
    //@{
    /** @brief Returns a pointer to the cell that contains a point
    @param position the point coordinates */
    Cell* getCell( Point3 const& position ) const;

    /** @brief Returns the cell edge length in a direction
    @param dir direction */
    double getCellSize( int const& dir ) const;

    /** @brief Returns a pointer to the vector of all cells of the linked cell
    grid */
    vector<Cell*> const* getAllCells() const;
    
    /** @brief Returns a list of pointers to the cell that contains a point
    and the neighboring cells to that cell
    @param position the point coordinates */
    list<Cell*> getCellAndCellNeighborhood( Point3 const& position ) const;    
    //@}


    /**@name Methods Set */
    //@{
    /** @brief Sets the linked cell grid in serial mode
    @param cellsize_ minimum cell edge length in each direction
    @param oshift empty string to shift the output */
    size_t set( double cellsize_, string const& oshift );

    /** @brief Sets the linked cell grid in parallel mode
    @param cellsize_ minimum cell edge length in each direction
    @param nprocsdir number of subdomains in each direction
    @param MPIcoords subdomain coordinates in the MPI Cartesian topology
    @param voisins neighboring subdomain in the MPI Cartesian topology
    @param oshift empty string to shift the output */
    size_t set( double cellsize_, int const* nprocsdir, int const* MPIcoords,
  	MPINeighbors const* voisins, string const& oshift );
    //@}


  private:
    /** @name Methods */
    //@{
    /** @brief Returns a pointer to the cell given its ijk indexing
    @param i index in the X direction
    @param j index in the Y direction
    @param k index in the Z direction */
    Cell* getCell( int i, int j, int k ) const;

    /** @brief Returns the cell number given its ijk indexing
    @param i index in the X direction
    @param j index in the Y direction
    @param k index in the Z direction */
    int getCellNumber( int i, int j, int k ) const;

    /** @brief Sets the list of neighboring cells over which broad phase
    contact detection is performed (3 cells above, 1 cell to the right, 9 cells
    behind) */
    void setCellContactNeighborhood();

    /** @brief Updates the link between the cells and a simple obstacle
    @param time physical time
    @param dt time step magnitude
    @param myObs simple obstacle */
    void LinkUpdate( double time, double dt, SimpleObstacle *myObs );
    //@}


    /** @name Parameters */
    //@{
    vector<Cell*> m_allcells; /**< vector of cells */
    list<Cell*> m_buffer_cells; /**< list of buffer cells (i.e. tag = 1) */
    int m_nb; /**< total number of cells */
    int m_nbi; /**< number of cells in X direction */
    int m_nbj; /**< number of cells in Y direction */
    int m_nbk; /**< number of cells in Z direction */
    double m_cellsize_X; /**< cell size in X direction */
    double m_cellsize_Y; /**< cell size in Y direction */
    double m_cellsize_Z; /**< cell size in Z direction */
    Point3 m_LC_global_origin; /**< Linked cell global origin */
    Point3 m_LC_global_max; /**< Linked cell global max point */    
    Point3 m_LC_local_origin; /**< Linked cell local origin */
    Point3 m_LC_local_max; /**< Linked cell local max point */      
    BBox* m_extendedBBox; /**< Bounding Box of the local Linked cell extended by
    	half a cell in each direction */
    //@}
};

#endif
