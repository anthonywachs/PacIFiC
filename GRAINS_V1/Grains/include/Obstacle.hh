#ifndef _OBSTACLE_HH_
#define _OBSTACLE_HH_

#include "Component.hh"
#include "ObstacleKinematicsVelocity.hh"
#include "ObstacleKinematicsForce.hh"
#include "Basic.hh"
#include "Vector3.hh"
using namespace solid;

class SimpleObstacle;
class ObstacleImposedVelocity;
class Quaternion;
class LinkedCell;
class LinkedBox;
class App;


/** @brief Data to store the state of an obstacle */
// ============================================================================
struct ObstacleState
{
  string nom; /**< obstacle name */
  ConfigurationMemento* memento_config; /**< obstacle configuration */
  ObstacleKinematicsMemento* memento_cine; /**< obstacle kinematics */
};    


/** @brief The class Obstacle.

    Rigid body with a prescribed motion or a motion controlled by a prescribed
    force load.

    @author G.FERRER - Institut Francais du Petrole - 1999 - Creation 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class Obstacle : public Component
{
  public:
    /** @name Constructors */
    //@{
    /** @brief Constructor with name and autonumbering as input parameters
    @param s obstacle name
    @param autonumbering obstacle autonumbering */
    Obstacle( string const& s = "obstacle", bool const& autonumbering = true );

    /** @brief Copy constructor from a Component
    @param copy copied Component 
    @param s obstacle name */
    Obstacle( Component& copy, char const* s = "obstacle" );

    /** @brief Destructor */
    virtual ~Obstacle();
    //@}


    /** @name Purely virtual methods */
    //@{
    /** @brief Adds an obstacle (single or composite) to the obstacle tree
    @param obstacle obstacle to be added */
    virtual void append( Obstacle* obstacle ) = 0;

    /** @brief Moves the obstacle and returns a list of moved obstacles
    @param time physical time
    @param dt time step magnitude
    @param b_deplaceCine_Comp whether to move the composite that the obstacle
    belongs to (imposed velocity)
    @param b_deplaceF_Comp whether to move the composite that the obstacle
    belongs to (imposed force) */
    virtual list<SimpleObstacle*> Move( double time,
	double dt, bool const& b_deplaceCine_Comp, 
        bool const& b_deplaceF_Comp ) = 0;

    /** @brief Returns a pointer to the obstacle if the name matches
    @param nom_ obstacle name */
    virtual Obstacle const* getObstacleFromName( string const& nom_ ) const = 0;

    /** @brief Returns a list of simple obstacles that belong to the obstacle
    (itself in the case of a single obstacle) */
    virtual list<SimpleObstacle*> getObstacles() = 0;
  
    /** @brief Returns a list of obstacles to be sent to the fluid solver
    in case of coupling with a fluid (itself in the case of a single 
    obstacle) */
    virtual list<Obstacle*> getObstaclesToFluid() = 0;

    /** @brief Rotates the obstacle with a quaternion
    @param rotation the quaternion defining the rotation */
    virtual void Rotate( Quaternion const& rotation ) = 0;

    /** @brief Deletes the obstacle if it belongs to a prescribed box
    @param box the box */
    virtual void Suppression( BBox const& box ) = 0;
  
    /** @brief Deletes an obstacle in the obstacle tree
    @param name_ name of obstacle to be deleted */
    virtual void DestroyObstacle( string const& name_ ) = 0; 
  
    /** @brief Deletes an obstacle in the obstacle tree and removes it from the
    LinkedCell
    @param name_ name of obstacle to be deleted
    @param LC linked-cell grid */
    virtual void ClearObstacle( string const& name_, LinkedCell* LC ) = 0;
    
    /** @brief Returns obstacle type */
    virtual string getObstacleType() = 0;    
    //@}


    /** @name Set methods */
    //@{
    /** @brief Sets kinematics
    @param kine_ the new kinematics */
    void setKinematics( ObstacleKinematicsVelocity& kine_ );
  
    /** @brief Sets kinematics using translational and angular velocities
    @param vtrans translational velocity
    @param vrot angular velocity */
    void setVelocity( Vector3 const* vtrans, Vector3 const* vrot );  
  
    /** @brief Sets indicator for Paraview post-processing */
    void setIndicator( double const& value );
  
    /** @brief Initializes all contact map entries to false */
    virtual void setContactMapToFalse();       
    //@}


    /** @name Accessors */
    //@{
    /** @brief Returns obstacle name */
    string getName() const;
  
    /** @brief Returns the velocity at a point in space based on the
    translational and angular velocity of the obstacle. This method assumes
    that the point belongs to the obstacle but this assumption is not verified.
    @param pt point where the velocity is computed */
    Vector3 getVelocityAtPoint( Point3 const& pt ) const;

    /** @brief Returns the angular velocity */
    Vector3 const* getAngularVelocity() const ;

    /** @brief Returns the translational velocity */
    Vector3 const* getTranslationalVelocity() const ;
  
    /** @brief Returns indicator for Paraview post-processing */
    double getIndicator() const;  
  
    /** @brief Returns a pointer to the torsor exerted on the obstacle */
    virtual Torsor const* getTorsor();   
    //@}

  
    /** @name State storing/restoring methods */
    //@{
    /** @brief Saves obstacle state */
    virtual void saveState();

    /** @brief Creates obstacle state and adds state to the list of states of
    all obstacles 
    @param obsStates list of states of all obstacles  */
    virtual void createState( list<struct ObstacleState*>& obsStates ) = 0; 
  
    /** @brief Restores obstacle state */
    virtual void restoreState();
  
    /** @brief Restores obstacle state
    @param obsStates list of states of all obstacles */
    virtual void restoreState( list<struct ObstacleState*>& obsStates ) = 0;  
    //@}


    /** @name Virtual methods */
    //@{
    /** @brief Links imposed kinematics to the obstacle and returns true if the
    linking process is successful
    @param imposed the imposed kinematics */
    virtual bool LinkImposedMotion( ObstacleImposedVelocity& imposed );
 
    /** @brief Links imposed force kinematics to the obstacle and returns true 
    if the linking process is successful
    @param imposed the imposed force */
    virtual bool LinkImposedMotion( ObstacleImposedForce& imposed );

    /** @brief Resets kinematics to 0 */
    virtual void resetKinematics();
  
    /** @brief Updates contact map */
    virtual void updateContactMap(); 
  
    /** @brief Does the contact exist in the map, if yes return the pointer to 
    the cumulative tangential displacement 
    @param tangentialDepl pointer to the cumulative tangential displacement 
    @param id id number of the other component */
    virtual bool ContactInMapIsActive( double*& tangentialDepl, int const& id );
  
    /** @brief Adds new contact in the map
    @param tangentialDepl initial tangential displacement 
    @param id id number of the other component */
    virtual void addNewContactInMap( double const& tangentialDepl, 
  	int const& id ); 

    /** @brief Increases cumulative tangential displacement with component id
    @param tangentialDepl additional tangential displacement 
    @param id id number of the other component */
    virtual void addDeplContactInMap( double const& tangentialDepl, 
  	int const& id ); 	
    //@}


    /** @name Methods */
    //@{
    /** @brief Composes the obstacle kinematics with another "higher level"
    velocity kinematics
    @param other the higher level kinematics
    @param lever lever arm of the higher level kinematics applied to the
    obstacle */
    void Compose( ObstacleKinematicsVelocity const& other, 
    	Vector3 const& lever );
  
    /** @brief Composes the obstacle kinematics with another "higher level"
    force kinematics
    @param other the higher level kinematics
    @param centre the obstacle center of mass */
    void Compose( ObstacleKinematicsForce const& other, Point3 const& centre );
  
    /** @brief Returns whether the obstacle has moved over the last time step */
    bool hasMoved() const; 
    
    /** @brief Contact between an obstacle and a component. If contact exists,
    computes the contact force and torque and adds to each component
    @exception ContactError if overlapping distance is larger than the sum of
    the crust thicknesses of the components
    @param voisin the other component
    @param dt time step magnitude
    @param time physical time
    @param LC linked-cell grid */
    virtual void InterAction( Component* voisin, 
	double dt, double const& time, LinkedCell* LC ) throw (ContactError);
    //@}


    /**@name Static methods */
    //@{
    /** @brief Returns total number of single obstacles */
    static int getTotalNbSingleObstacles();
  
    /** @brief Sets the boolean to actually move obstacles 
    param depObs boolean to actually move obstacles  */
    static void setMoveObstacle( bool const& depObs );
  
    /** @brief Gets the boolean to actually move obstacles */
    static bool getMoveObstacle();    
    //@}

  
    /**@name I/O methods */
    //@{
    /** @brief Reloads the obstacle and links it to the higher level obstacle in
    the obstacle tree
    @param mother higher level obstacle
    @param file input stream */
    virtual void reload( Obstacle& mother, istream& file ) = 0;

    /** @brief Writes the identity of the obstacle i.e its id number and 
    name */
    void writeIdentity( ostream& file ) const;

    /** @brief Writes the obstacle position
    @param position output stream */
    virtual void writePosition( ostream& position );

    /** @brief Writes the obstacle's "static" data
    @param fileOut output stream */
    virtual void writeStatic( ostream& fileOut ) const;
  
    /** @brief Updates indicator for Paraview post-processing
    @param time physical time
    @param dt time step magnitude */
    virtual void updateIndicator( double time, double dt ) {}

    /** @brief Returns the number of points to write the obstacle in a
    Paraview format */
    virtual int numberOfPoints_PARAVIEW() const ;

    /** @brief Returns the number of elementary polytopes to write the obstacle 
    shape in a Paraview format */
    virtual int numberOfCells_PARAVIEW() const;

    /** @brief Returns a list of points describing the obstacle in a
    Paraview format 
    @param translation additional center of mass translation */
    virtual list<Point3> get_polygonsPts_PARAVIEW(
  	Vector3 const* translation = NULL ) const ;

    /** @brief Writes the points describing the obstacle in a
    Paraview format 
    @param f output stream
    @param translation additional center of mass translation */
    virtual void write_polygonsPts_PARAVIEW( ostream& f, 
  	Vector3 const* translation = NULL ) const ; 
  
    /** @brief Writes the obstacle in a Paraview format
    @param connectivity connectivity of Paraview polytopes
    @param offsets connectivity offsets
    @param cellstype Paraview polytopes type
    @param firstpoint_globalnumber global number of the 1st point
    @param last_offset last offset used for the previous convex shape */
    virtual void write_polygonsStr_PARAVIEW(list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset) const ;

    /**  @brief Outputs information to be transferred to the fluid
    @param fluid output stream */
    virtual void writePositionInFluid( ostream& fluid );
    //@}


  protected:
    /** @name Parameters */
    //@{
    ObstacleKinematicsVelocity m_kinematics; /**< obstacle kinematics with 
    	imposed velocity */
    ObstacleKinematicsForce m_confinement; /**< obstacle kinematics with imposed
  	force */
    string m_name; /**< obstacle name */
    bool m_ismoving; /**< whether the obstacle moves or not */
    double m_indicator; /**< post-processing indicator for the rotation of
  	composite obstacle in Paraview */
    string m_ObstacleType; /**< obstacle type */
    //@}

    /** @name Parameters Static */
    //@{
    static int m_totalNbSingleObstacles; /**< Total number of single 
    	obstacles */
    static bool m_MoveObstacle; /**< enables to imposed a kinematics while 
    	not actually moving obstacles if set to false, useful in periodic 
	cases */
    static bool m_isConfinement; /**< true if imposed force */
  //@}
      
};

#endif
