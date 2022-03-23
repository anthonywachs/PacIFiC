#ifndef _SIMPLEOBSTACLE_HH_
#define _SIMPLEOBSTACLE_HH_

#include "Obstacle.hh"
#include <list>
#include <set>
using namespace std;

#include "ReaderXML.hh"

class Particle;
class App;
class Cell;
class LinkedCell;
class PeriodicObstacle;


/** @brief The class SimpleObstacle.

    A simple obstacle is an obstacle composed of a single convex rigid body.

    @author G.FERRER - Institut Francais du Petrole - 2002 - Creation 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class SimpleObstacle : public Obstacle
{
  public:
    /** @name Constructeurs */
    //@{
    /** @brief Constructor with name as input parameter
    @param s obstacle name */ 
    SimpleObstacle( string const& s = "" );

    /** @brief Constructor with an XML node as an input parameter
    @param root XML node */
    SimpleObstacle( DOMNode* root );

    /** @brief Copy constructor from a Component
    @param copy copied Component 
    @param s obstacle name */
    SimpleObstacle( Component& copy, char const* s = "obstacle" );
  
    /** @brief Constructor with a rigid body, a name and a material as input
    parameters
    @param geoRBWC rigid body
    @param name obstacle name
    @param materialName material name
    @param transferToFluid_ whether the obstacle is transferred to the fluid */
    SimpleObstacle( RigidBodyWithCrust* geoRBWC, string const& name = "",
      string const& materialName = "", bool const& transferToFluid_ = false );
  
    /** @brief Destructor */
    ~SimpleObstacle();
    //@}


    /** @name Accessors */
    //@{
    /** @brief Returns a pointer to the obstacle if the name matches
    @param nom_ obstacle name */
    const Obstacle* getObstacleFromName( string const& nom_ ) const;

    /** @brief Returns a list of simple obstacles that belong to the obstacle.
    Here returns this (i.e. itself) */
    list<SimpleObstacle*> getObstacles();

    /** @brief Returns a list of obstacles to be sent to the fluid solver
    in case of coupling with a fluid. Here returns this (i.e. itself) */
    list<Obstacle*> getObstaclesToFluid() ;
  
    /** @brief Returns a pointer to the list of cells the obstacle is linked 
    to */
    list<Cell*> const* getInCells() const;
  
    /** @brief Returns the bounding box of the obstacle */
    BBox const* getObstacleBox() const;  
  
    /** @brief Returns the frequency at which the obstacle link to the cells of
    the linked-cell grid is upadted */
    int getObstacleLinkedCellUpdateFrequency() const 
  	{ return m_LinkUpdate_frequency; }

    /** @brief Returns obstacle type */
    string getObstacleType() ; 
    //@}


    /** @name Methods */
    //@{
    /** @brief Adds an obstacle (single or composite) to the composite obstacle
    tree. Should not be called by SimpleObstacle.
    @param obstacle obstacle to be added */
    void append( Obstacle* obstacle ) ;

    /** @brief Moves the simple obstacle and returns a list of moved 
    obstacles (here itself)
    @param time physical time
    @param dt time step magnitude
    @param b_deplaceCine_Comp whether to move the composite that the composite 
    obstacle belongs to (imposed velocity)
    @param b_deplaceF_Comp whether to move the composite that the composite 
    obstacle belongs to (imposed force) */
    list<SimpleObstacle*> Move( double time,
	double dt, bool const& b_deplaceCine_Comp, 
        bool const& b_deplaceF_Comp ) ;
	
    /** @brief Contact between a simple obstacle and a component. If contact 
    exists, computes the contact force and torque and adds to each component
    @exception ContactError if overlapping distance is larger than the sum of
    the crust thicknesses of the components
    @param voisin the other component
    @param dt time step magnitude
    @param time physical time
    @param LC linked-cell grid */
    virtual void InterAction( Component* voisin, 
	double dt, double const& time, LinkedCell* LC ) throw (ContactError);	

    /** @brief Rotates the obstacle with a quaternion
    @param rotation the quaternion defining the rotation */
    void Rotate( Quaternion const& rotation );

    /** @brief Deletes the obstacle if it belongs to a prescribed box
    @param box the box */
    void Suppression( BBox const& box );

    /** @brief Translates the obstacle
    @param translation translation vector */
    void Translate( Vector3 const& translation );
  
    /** @brief Deletes an obstacle in the obstacle tree
    @param name_ name of obstacle to be deleted */
    void DestroyObstacle( string const& name_ ); 
  
    /** @brief Deletes an obstacle in the obstacle tree and removes it from the
    LinkedCell
    @param name_ name of obstacle to be deleted
    @param LC linked-cell grid */
    void ClearObstacle( string const& name_, LinkedCell* LC );
		  
    /** @brief Adds a cell to the list of cells the obstacle is linked to
    @param cel_ the cell */
    void add( Cell* cel_ ); 
  
    /** @brief Empties the list of cells the obstacle is linked to and deletes
    the pointer to the obstacle is these cells */
    void resetInCells();
	
    /** @brief Returns whether an update of the link between the obstacle and the
    linked-cell grid is required. If yes, returns true and sets the counter to
    0, if no, returns false and increments the counter */
    bool performLinkUpdate(); 
  
    /** @brief Returns the maximum of the absolute value of the obstacle 
    velocity in each direction */
    Vector3 vitesseMaxPerDirection() const;
  
    /** @brief Updates contact map */
    void updateContactMap(); 
  
    /** @brief Does the contact exist in the map, if yes return the pointer to 
    the cumulative tangential displacement 
    @param tangentialDepl pointer to the cumulative tangential displacement 
    @param id id number of the other component */
    bool ContactInMapIsActive( double*& tangentialDepl, int const& id );
  
    /** @brief Adds new contact in the map
    @param tangentialDepl initial tangential displacement 
    @param id id number of the other component */
    void addNewContactInMap( double const& tangentialDepl, 
  	int const& id ); 

    /** @brief Increases cumulative tangential displacement with component id
    @param tangentialDepl additional tangential displacement 
    @param id id number of the other component */
    void addDeplContactInMap( double const& tangentialDepl, 
  	int const& id );  	     		    
    //@}


    /** @name Set methods */
    //@{
    /** @brief Initializes all contact map entries to false */
    void setContactMapToFalse();
  
    /** @brief Sets the frequency at which the obstacle link to 
    the cells of the linked-cell grid is updated 
    @param updateFreq updating frequency */
    void setObstacleLinkedCellUpdateFrequency( int const& updateFreq );      
    //@}


    /** @name State storing/restoring methods */
    //@{
    /** @brief Creates obstacle state and adds state to the list of states of
    all obstacles 
    @param obsStates list of states of all obstacles  */
    void createState( list<struct ObstacleState*>& obsStates ) ; 
  
    /** @brief Restores obstacle state
    @param obsStates list of states of all obstacles */
    void restoreState( list<struct ObstacleState*>& obsStates ) ;  
    //@}


    /** @name I/O methods */
    //@{
    /** @brief Reloads the simple obstacle and links it to the higher level 
    obstacle in the obstacle tree
    @param mother higher level obstacle
    @param file input stream */
    virtual void reload( Obstacle& mother, istream& file ) ;

    /** @brief Outputs the simple obstacle for reload
    @param fileSave output stream */
    virtual void write( ostream& fileSave ) const;

    /** @brief Returns the number of points to write the simple obstacle in a
    Paraview format */
    virtual int numberOfPoints_PARAVIEW() const ;

    /** @brief Returns the number of elementary polytopes to write the 
    simple obstacle shape in a Paraview format */
    virtual int numberOfCells_PARAVIEW() const;

    /** @brief Returns a list of points describing the simple obstacle in a
    Paraview format 
    @param translation additional center of mass translation */
    virtual list<Point3> get_polygonsPts_PARAVIEW(
  	Vector3 const* translation = NULL ) const ;

    /** @brief Writes the points describing the simple obstacle in a
    Paraview format 
    @param f output stream
    @param translation additional center of mass translation */
    virtual void write_polygonsPts_PARAVIEW( ostream& f, 
  	Vector3 const* translation = NULL ) const ; 
  
    /** @brief Writes the simple obstacle in a Paraview format
    @param connectivity connectivity of Paraview polytopes
    @param offsets connectivity offsets
    @param cellstype Paraview polytopes type
    @param firstpoint_globalnumber global number of the 1st point
    @param last_offset last offset used for the previous convex shape */
    virtual void write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const ;

    /**  @brief Outputs information to be transferred to the fluid
    @param fluid output stream */
    virtual void writePositionInFluid( ostream& fluid );
    //@}

  
  protected:
    /**@name Parameters */
    //@{  
    BBox m_obstacleBox; /**< bounding box of the obstacle */
    list<Cell*> m_inCells; /**< cells the obstacle is linked to */
    int m_LinkUpdate_frequency; /**< frequency at which the obstacle link to 
    	the cells of the linked-cell grid is upadted */
    int m_LinkUpdate_counter; /**< counter of updates between the obstacle and 
    	the linked-cell grid */
    bool m_transferToFluid; /**< whether to transfer the obstacle to the fluid
  	solver or not in case of coupling to a fluid solver */
    //@}
  	  
};

#endif

