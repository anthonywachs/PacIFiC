#ifndef _CELL_HH_
#define _CELL_HH_


#include "Basic.hh"
#include <list>
using namespace std;


/** @brief Geographic position in the linked-cell grid
    @author GRAINS Project - IFP - 2009 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
enum GeoPosition
  {
    GEOPOS_NORTH		= 0,
    GEOPOS_NORTH_EAST		= 1,
    GEOPOS_NORTH_WEST		= 2,
    GEOPOS_NORTH_FRONT		= 3,
    GEOPOS_NORTH_BEHIND		= 4,
    GEOPOS_NORTH_EAST_FRONT	= 5,
    GEOPOS_NORTH_EAST_BEHIND	= 6, 
    GEOPOS_NORTH_WEST_FRONT	= 7,
    GEOPOS_NORTH_WEST_BEHIND	= 8, 
    GEOPOS_SOUTH		= 9,
    GEOPOS_SOUTH_EAST		= 10,
    GEOPOS_SOUTH_WEST		= 11,
    GEOPOS_SOUTH_FRONT		= 12,
    GEOPOS_SOUTH_BEHIND		= 13,
    GEOPOS_SOUTH_EAST_FRONT	= 14,
    GEOPOS_SOUTH_EAST_BEHIND	= 15, 
    GEOPOS_SOUTH_WEST_FRONT	= 16,
    GEOPOS_SOUTH_WEST_BEHIND	= 17,
    GEOPOS_EAST			= 18,
    GEOPOS_WEST			= 19,
    GEOPOS_EAST_FRONT		= 20,
    GEOPOS_EAST_BEHIND		= 21,
    GEOPOS_WEST_FRONT		= 22,
    GEOPOS_WEST_BEHIND		= 23,
    GEOPOS_FRONT		= 24,
    GEOPOS_BEHIND		= 25,
    GEOPOS_NONE			= 26,
    // Special cases for periodicity with a single cell in the periodic
    // direction
    // Note: all cases are not considered so far, only those with practical
    // applications
    GEOPOS_EASTWEST		= 27,
    GEOPOS_NORTHSOUTH		= 28,
    GEOPOS_EASTWESTNORTHSOUTH   = 29                              
  };
  
    
#include "Particle.hh"
#include "SimpleObstacle.hh"  


/** @brief The class Cell.

    Defines a cell of the linked-cell grid. The cell has a parallelepipedic
    shape by construction.
    
    @author G.FERRER - Institut Francais du Petrole - 1999 - Creation 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class Cell
{
  public:
    /**@name Constructors */
    //@{
    /** @brief Default constructor */
    Cell();

    /** @brief Constructor with input parameters
    @param id1 cell number
    @param x index in the x direction
    @param y index in the y direction
    @param z index in the z direction
    @param arete_X Edge length in the x direction 
    @param arete_Y Edge length in the y direction      
    @param arete_Z Edge length in the z direction 
    @param tag_ cell tag 
    @param OL local origin of the linked-cell grid 
    @param xmax_ maximum coordinate of the linked-cell grid in the x direction
    @param ymax_ maximum coordinate of the linked-cell grid in the y direction
    @param zmax_ maximum coordinate of the linked-cell grid in the z 
    	direction
    @param geoloc_ geographic position */
    Cell( int id1, int x, int y, int z, Point3 const& OL,
    	double arete_X, double arete_Y, double arete_Z, 
	double xmax_, double ymax_, double zmax_,
	int tag_ = 0, GeoPosition geoloc_ = GEOPOS_NONE );  

    /** @brief Destructor */
    ~Cell();
    //@}


    /**@name Methods Get */
    //@{  
    /** @brief Returns the complete neighborhood of the cell */
    list<Cell*> const* getCompleteNeighborhood() const;  

    /** @brief Returns the number of obstacles in the vicinity of the cell */
    int numberOfObstacles() const;

    /** @brief Returns the volume of all particles whose center of mass is in
    the cell */
    double getVolumeParticles();

    /** @brief Returns the cell center */
    Point3 const* getCentre() const;
  
    /** @brief Returns a pointer to the list of particles in the cell */
    list<Particle*>* getParticles();
  
    /** @brief Returns the cell geographic position */
    GeoPosition getGeoPosition() const;  
    //@}


    /**@name Methods */
    //@{
    /** @brief Adds a particle to the list of particles in the cell. No
    geometric test is performed, the particle pointer is simply added to the
    list
    @param particle_ particle to be added */
    void add( Particle* particle_ );

    /** @brief Adds a cell to the list of neighboring cells for contact 
    detection
    @param neighbor neighboring cell */
    void addNeighboringCellContact( Cell* neighbor );

    /** @brief Adds a cell to the list of all neighboring cells
    @param neighbor neighboring cell */
    void addNeighboringCell( Cell* neighbor );
  
    /** @brief Adds an obstacle to the list of obstacles in the vicinity of the 
    cell. No geometric test is performed, the obstacle pointer is simply added 
    to the list
    @param obstacle_ obstacle to be added */
    void addObstacle( SimpleObstacle* obstacle_ );  

    /** @brief Clears the list of particles that belong to the cell  */
    void clearParticles();

    /** @brief Returns whether a particle belongs to the list of particles in
    the cell
    @param particle_ particle */
    bool contains( Particle* particle_ );

    /** @brief Returns whether a particle is in contact with another component
    in the vicinity of the cell
    @param particle_ particle */
    bool isContact( Particle const* particle_ ) const;
  
    /** @brief Returns whether a particle is in contact with another component
    in the vicinity of the cell. The contact detection is performed with the
    crust width
    @param particle_ particle */
    bool isContactWithCrust( Particle const* particle_ ) const;   
  
    /** @brief Returns whether a particle is close to another component
    in the vicinity of the cell 
    @param particle_ particle */
    bool isClose( Particle const* particle_ ) const; 
  
    /** @brief Returns whether a particle is close to another component
    in the vicinity of the cell. The closeness detection is performed with the
    crust width
    @param particle_ particle */
    bool isCloseWithCrust( Particle const* particle_ ) const;      

    /** @brief Returns whether the cell does not contain any particle */
    bool isEmpty() const;

    /** @brief Returns a list of particles that do not belong to the cell
    anymore
    @param particlesExit list of particles that do not belong to the cell */
    void LinkUpdate( list<Particle*>& particlesExit );

    /** @brief Removes a particle from the list of particles in the cell
    @param particle_ particle to be removed */
    void remove( Particle* particle_ );
  
    /** @brief Removes an obstacle from the list of obstacles in the vicinity of
    the cell
    @param obs obstacle to be removed */
    void remove( SimpleObstacle* obs ); 
    
    /** @brief Returns whether a point belongs to any particle that belongs to 
    the cell
    @param position the point coordinates */
    bool isInParticle( Point3 const& position ) const;
    
    /** @brief Returns whether a point belongs to a specific particle that 
    belongs to the cell
    @param position the point coordinates 
    @param particle_ the specific particle */
    bool isInParticle( Point3 const& position, Particle const* particle_ ) 
    	const;                 
  
    /** @brief Output operator
    @param f output stream
    @param C the cell */
    friend ostream& operator << ( ostream& f, Cell const& C );  
    //@}


    /** @name Methods Static */
    //@{
    /** @brief Returns the cell ijk indices that contains a point
    @param position the point coordinates
    @param id cell ijk indices */
    static void GetCell( Point3 const& position, int* id );

    /** @brief Sets the number of cells of the linked-cell in each direction
    @param nbX number of cells in the x direction
    @param nbY number of cells in the y direction
    @param nbZ number of cells in the z direction */
    static void setNbCellsPerDirection( int nbX, int nbY, int nbZ );
  
    /** @brief Returns the geographic position name
    @param geoloc_ geographic position */
    static string getGeoPositionName( GeoPosition geoloc_ );
  
    /** @brief Returns the geographic position name
    @param geoloc_ geographic position */
    static string getGeoPositionName( int geoloc_ );   
    //@}


    /**@name Operators */
    //@{
    /** @brief Returns the cell index in direction dir in the linked-cell
    @param dir space direction 0, 1 or 2 */
    int operator [] ( int dir ) const;

    /** @brief Comparison operator based on addresses
    @param cell_ the other cell */
    bool operator == ( Cell const& cell_ ) const;
    //@}


    /**@name Class Friend */
    //@{
    friend class LinkedCell; /**< linked-cell grid the cell is part of */
    //@}
  

  private:
    /**@name Parameters */
    //@{  
    int m_number; /**< Cell number */
    int m_tag;   /**< Cell tag:
    <ul> 
      <li> 0=interior, 
      <li> 1=buffer zone, 
      <li> 2=clone zone
    </ul> */
    Point3 m_centre; /**< cell center coordinates */
    GeoPosition m_GeoPosCell; /**< geographic position in the linked-cell
    	grid, useful for cells tagged 1 */
    int m_cel[3]; /**< ijk indexing of the cell in the linked-cell */
    list<Cell*> m_neighborsContact; /**< Neighboring cells for contact 
    	detection: 3 above, 1 on the right, 9 behind */
    list<Cell*> m_allNeighbors; /**< All neighboring cells (in the general case
    	one cell has 26 neighboring cells  */
    list<SimpleObstacle*> m_obstacles; /**< obstacles in the vicinity of the
    	cell */		
    list<Particle*> m_particles; /**< Particles whose center of mass belongs to
    	the cell */
    static int m_nbi; /**< number of cells in the linked-cell in the 
    	x direction */
    static int m_nbj; /**< number of cells in the linked-cell in the 
    	y direction */
    static int m_nbk; /**< number of cells in the linked-cell in the 
    	z direction */  
    static double m_edge_X; /**< Edge length in the x direction */
    static double m_edge_Y; /**< Edge length in the y direction */
    static double m_edge_Z; /**< Edge length in the z direction */
    static Point3 m_LC_local_origin; /**< Local origin of the linked-cell 
    	grid */
    static double m_LC_local_xmax; /**< maximum coordinate of the linked-cell 
    	grid in the x direction */
    static double m_LC_local_ymax; /**< maximum coordinate of the linked-cell 
    	grid in the y direction */	
    static double m_LC_local_zmax; /**< maximum coordinate of the linked-cell 
    	grid in the z direction */
    //@}

  
    /**@name Methods Set */
    //@{  
    /** @brief Computes the cell center coordinates */
    void setCentre(); 
    
    /** @brief Returns the geographic position name
    @param geoloc_ geolocalisation */
    static string getGeoPositionName_generic( int geoloc_ );
    //@}  	
};

#endif
