#ifndef _COMPOSITEOBSTACLE_HH_
#define _COMPOSITEOBSTACLE_HH_

#include "Obstacle.hh"
#include "SimpleObstacle.hh"

#include <list>
using namespace std;

#include "ReaderXML.hh"


/** @brief The class CompositeObstacle

    Cluster of obstacles that can be either single obstacles or composite
    obstacles. Useful to impose the same motion to the whole cluster or to
    compose motions.

    @author G.FERRER - Institut Francais du Petrole - 2002 - Creation
    @author D. RAKOTONIRINA - IFP Energies nouvelles - 2014 - Modification 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class CompositeObstacle : public Obstacle
{
  public:
    /** @name Constructors */
    //@{
    /** @brief Constructor with name as input parameter
    @param s obstacle name */
    CompositeObstacle( string const& s = "" );

    /** @brief Constructor with an XML node as an input parameter
    @param root XML node */
    CompositeObstacle( DOMNode* root );

    /** @brief Destructor */
    ~CompositeObstacle();
    //@}


    /** @name Accessors */
    //@{
    /**
    @brief Returns a pointer to the obstacle if the name matches. Searches all
    leaves of the composite obstacle tree
    @param nom_ obstacle name */
    Obstacle const* getObstacleFromName( string const& nom_ ) const;

    /** @brief Returns a list of simple obstacles that belong to the composite 
    obstacle */
    list<SimpleObstacle*> getObstacles() ;

    /** @brief Returns a list of obstacles to be sent to the fluid solver
    in case of coupling with a fluid */
    list<Obstacle*> getObstaclesToFluid() ;

    /** @brief Returns obstacle type */
    string getObstacleType() ;
    
    /** @brief Returns a pointer to the torsor exerted on the composite 
    obstacle */
    Torsor const* getTorsor();       
    
    /** @brief Returns the volume of the composite obstacle */
    double getVolume() const; 
    
    /** @brief Returns the radius of the sphere of volume equivalent to that of 
    the composite obstacle */ 
    double getEquivalentSphereRadius() const;       
    //@}


    /** @name Set methods */
    //@{  
    /** @brief Initializes all contact map entries to false */
    virtual void setContactMapToFalse();
    
    /** @brief Set contact map cumulative features to zero */
    virtual void setContactMapCumulativeFeaturesToZero();            
    //@}


    /** @name Methods */
    //@{
    /** @brief Links imposed kinematics to the obstacle and returns true if the
    linking process is successful
    @param imposed the imposed kinematics */
    virtual bool LinkImposedMotion( ObstacleImposedVelocity* imposed );
 
    /** @brief Links imposed force kinematics to the obstacle and returns true 
    if the linking process is successful
    @param imposed the imposed force */
    virtual bool LinkImposedMotion( ObstacleImposedForce* imposed );

    /** @brief Adds an obstacle (single or composite) to the composite obstacle
    tree
    @param obstacle obstacle to be added */
    void append( Obstacle* obstacle ) ;

    /** @brief Moves the composite obstacle and returns a list of moved 
    obstacles
    @param time physical time
    @param dt time step magnitude
    @param motherCompositeHasImposedVelocity whether the composite that the 
    obstacle belongs to has a non-zero imposed velocity
    @param motherCompositeHasImposedForce whether the composite that the 
    obstacle belongs to has a non-zero imposed force */
    list<SimpleObstacle*> Move( double time,
	double dt, bool const& motherCompositeHasImposedVelocity, 
        bool const& motherCompositeHasImposedForce ) ;

    /** @brief Returns whether the component is a composite obstacle ? */
    virtual bool isCompositeObstacle() const;

    /** @brief Returns whether there is geometric contact with another
    component 
    @param voisin the other component */
    virtual bool isContact( Component const* voisin ) const;

    /** @brief Returns whether there is geometric contact with another
    component accounting for crust thickness 
    @param voisin the other component */
    virtual bool isContactWithCrust( Component const* voisin ) const;

    /** @brief Returns whether there is geometric proximity with another
    component in the sense of whether their respective bounding boxes overlap 
    @param voisin the other component */
    virtual bool isClose( Component const* voisin ) const;

    /** @brief Returns whether there is geometric proximity with another
    component in the sense of whether their respective bounding boxes minus 
    their crust thickness overlap 
    @param voisin the other component */
    virtual bool isCloseWithCrust( Component const* voisin ) const;  

    /** @brief Resets kinematics to 0 */
    void resetKinematics();

    /** @brief Rotates the composite obstacle with a quaternion about its 
    center of mass. Warning: This method is not defined and must not be called
    @param rotation the quaternion defining the rotation */
    void Rotate( Quaternion const& rotation ) ;

    /** @brief Deletes the composite obstacle if it belongs to a prescribed box
    @param box the box */
    void Suppression( BBox const& box ) ;

    /** @brief Translates the composite obstacle
    @param translation translation vector */
    void Translate( Vector3 const& translation ) ;

    /** @brief Initializes the composite obstacle's torsor
    @param withWeight intializes the force to the composite obstacle weight if 
    true or to 0 is false */
    void InitializeForce( bool const& withWeight ); 

    /** @brief Deletes an obstacle in the obstacle tree
    @param name_ name of obstacle to be deleted */
    void DestroyObstacle( string const& name_ ) ; 
  
    /** @brief Deletes an obstacle in the obstacle tree and removes 
    it from the LinkedCell
    @param name_ name of obstacle to be deleted
    @param LC linked-cell grid */
    void ClearObstacle( string const& name_, LinkedCell* LC ) ;
      
    /** @brief Update contact map */
    virtual void updateContactMap();

    /** @brief Does the contact exist in the map? If so, return true and make
    kdelta, prev_normal and cumulSpringTorque point to the memorized info. 
    Otherwise, return false and set those pointers to NULL.
    @param id key in the map
    @param kdelta pointer to the memory of the vector kt * delta_t
    @param prev_normal pointer to the previous vector normal to the contact 
    plane
    @param cumulSpringTorque pointer to the memory of the spring-like component 
    of the friction torque 
    @param createContact when true, create contact if it does not exist */
    virtual bool getContactMemory( std::tuple<int,int,int> const& id,
  	Vector3* &kdelta, Vector3* &prev_normal, Vector3* &cumulSpringTorque,
  	bool createContact );

    /** @brief Adds new contact in the map
    @param id key in the map
    @param kdelta kt * delta_t vector
    @param prev_normal pointer to the previous vector normal to the contact 
    plane
    @param cumulSpringTorque pointer to the memory of the spring-like component 
    of the friction torque */
    virtual void addNewContactInMap( std::tuple<int,int,int> const& id,
  	Vector3 const& kdelta, Vector3 const& prev_normal,
  	Vector3 const& cumulSpringTorque );

    /** @brief Stores memory of the contact with component id: increase 
    cumulative tangential motion and cumulative spring torque, remember 
    contact normal.
    @param id key in the map
    @param kdelta kt * delta_t vector
    @param prev_normal pointer to the previous vector normal to the contact 
    plane
    @param cumulSpringTorque pointer to the memory of the spring-like component 
    of the friction torque */
    virtual void addDeplContactInMap( std::tuple<int,int,int> const& id,
  	Vector3 const& kdelta, Vector3 const& prev_normal,
  	Vector3 const& cumulSpringTorque );

    /** @brief Writes the contact map information in an array of doubles
    @param destination the array of double where the contact map should be 
    stored
    @param start_index the index of destination where the copy should start */
    virtual void copyContactMap( double* destination, int start_index );

    /** @brief Adds a single contact info to the contact map
    @param id key in the map
    @param isActive boolean: true if the contact is active, false otherwise
    @param kdelta kt * delta_t vector
    @param prev_normal pointer to the previous vector normal to the contact 
    plane
    @param cumulSpringTorque pointer to the memory of the spring-like component 
    of the friction torque */
    virtual void copyContactInMap( std::tuple<int,int,int> const& id,
  	bool const& isActive, Vector3 const& kdelta, Vector3 const& prev_normal,
  	Vector3 const& cumulSpringTorque );

    /** @brief Returns the number of contacts in the contact map */
    virtual int getContactMapSize();

    /** @brief Displays the active neighbours in the 
    format "my_elementary_id/neighbour_id/neightbout_elementary_id ; ...". 
    Useful for debugging only.
    @param id id of this component */
    virtual void printActiveNeighbors( int const& id );
	
    /** @ brief Returns whether a point lies inside the composite obstacle
    @param pt point */
    bool isIn( Point3 const& pt ) const;
            
    /** @brief Resets the minimum ID number of an obstacle for autonumbering */
    virtual void setMinIDnumber(); 
    
    /** @brief Checks if there is anything special to do about periodicity and
    if there is applies periodicity 
    @param LC linked-cell grid */
    virtual void periodicity( LinkedCell* LC );
    
    /** @brief Empties the list of cells the obstacle is linked to and deletes
    the pointer to the obstacle is these cells */
    virtual void resetInCells();              
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
    /** @brief Reloads the composite obstacle and links it to the higher level 
    obstacle in the obstacle tree
    @param mother higher level obstacle
    @param file input stream */
    virtual void reload( Obstacle& mother, istream& file ) ;    
    
    /** @brief Outputs the composite obstacle for reload
    @param fileSave output stream */
    virtual void write( ostream& fileSave ) const;

    /** @brief Writes the composite obstacle position
    @param position output stream */
    virtual void writePosition( ostream& position );

    /** @brief Writes the composite obstacle's "static" data
    @param fileOut output stream */
    virtual void writeStatic( ostream& fileOut ) const;
  
    /** @brief Updates indicator for Paraview post-processing
    @param time physical time
    @param dt time step magnitude */
    virtual void updateIndicator( double time, double dt );

    /** @brief Returns the number of points to write the composite obstacle in a
    Paraview format */
    virtual int numberOfPoints_PARAVIEW() const ;

    /** @brief Returns the number of elementary polytopes to write the 
    composite obstacle shape in a Paraview format */
    virtual int numberOfCells_PARAVIEW() const;

    /** @brief Returns a list of points describing the composite obstacle in a
    Paraview format 
    @param translation additional center of mass translation */
    virtual list<Point3> get_polygonsPts_PARAVIEW(
  	Vector3 const* translation = NULL ) const ;

    /** @brief Writes the points describing the composite obstacle in a
    Paraview format 
    @param f output stream
    @param translation additional center of mass translation */
    virtual void write_polygonsPts_PARAVIEW( ostream& f, 
  	Vector3 const* translation = NULL ) const ; 
  
    /** @brief Writes the composite obstacle in a Paraview format
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
    /** @name Methods */
    //@{  
    /** @brief Computes center of mass position */
    pair<Point3,double> computeVolumeCenterOfMass();
    //@}
    

    /** @name Parameters */
    //@{  
    list<Obstacle*> m_obstacles; /**< list of obstacles in the composite
    	obstacle */
    int m_CompositeObstacle_id; /**< standard ID number is 0, the composite
    	obstacle ID number is used to distinguish the highest level composite
    	obstacle that is artificially created as the root of the obstacle tree 
	from all other physically relevant and user created composite 
	obstacles */    	 
    static int m_minCompositeObstacleID; /**< Minimum composite
    	obstacle ID number, composite obstacle ID numbers range
    	from 0 (the root) down to m_min_CompositeObstacle_id and are 
	therefore always negative except the root */
    string m_type; /**< type of composite obstacle: standard or one of the
    	specific derived types such as CylindricalShell or RoughWall */
    double m_volume; /**< approximate volume of the composite obstacles that 
    	ignores geometric overlaps between obstacles */	 
    //@}


  private:
    /** @name Constructors */
    //@{
    /** @brief Copy constructor
    @param copy copied CompositeObstacle */
    CompositeObstacle( CompositeObstacle const& copy );
    //@}   
};

#endif
