#ifndef _COMPONENT_HH_
#define _COMPONENT_HH_

#include "Basic.hh"
#include "RigidBodyWithCrust.hh"
#include "Point3.hh"
#include "PointContact.hh"
#include "Torsor.hh"
#include "CollisionHistory.hh"
#include <list>
#include <string>
#include <map>
#include <set>
using namespace std;

class BBox;
class ContactError;
class ConfigurationMemento;
class LinkedCell;


struct ContactInfos
{
  PointContact ContactPoint; /** contact point */
  Component* p0; /** 1st component involved in contact */
  Component* p1; /** 2nd component involved in contact */
};


/** @brief The class Component.

    Physical component of a granular simulation (can be particle or obstacle).

    @author G.FERRER - Institut Francais du Petrole - 2000 - Creation
    @author A.WACHS - 2019 - Major cleaning & refactoring 
    @author D.HUET - 2022 - Contact force model with memory */
// ============================================================================
class Component
{
  public:
    /** @name Constructors */
    //@{
    /** @brief Destructor */
    virtual ~Component();
    //@}


    /** @name Virtual methods */
    //@{
    /** @brief Returns the velocity at a point in space based on the
    translational and angular velocity of the component. This method assumes
    that the point belongs to the component but this assumption is not verified.
    @param pt point where the velocity is computed */
    virtual Vector3 getVelocityAtPoint( Point3 const& pt ) const = 0;

    /** @brief Returns the angular velocity */
    virtual Vector3 const* getAngularVelocity() const = 0;

    /** @brief Returns the translational velocity */
    virtual Vector3 const* getTranslationalVelocity() const = 0;

    /** @brief Returns whether there is geometric contact with another
    component. Note: the other component must not be of the derived type
    CompositeObstacle
    @param voisin the other component */
    virtual bool isContact( Component const* voisin ) const;

    /** @brief Returns whether there is geometric contact with another
    component accounting for crust thickness. Note: the other component must
    not be of the derived type CompositeObstacle
    @param voisin the other component */
    virtual bool isContactWithCrust( Component const* voisin ) const;

    /** @brief Returns whether there is geometric proximity with another
    component in the sense of whether their respective bounding boxes overlap.
    Note: the other component must not be of the derived type CompositeObstacle
    @param voisin the other component */
    virtual bool isClose( Component const* voisin ) const;

    /** @brief Returns whether there is geometric proximity with another
    component in the sense of whether their respective bounding boxes minus
    their crust thickness overlap. Note: the other component must not be of the
    derived type CompositeObstacle
    @param voisin the other component */
    virtual bool isCloseWithCrust( Component const* voisin ) const;

    /** @brief Rotates the component using a quaternion about its center of mass
    @param rotation quaternion representing the rotation */
    virtual void Rotate( Quaternion const& rotation );

    /** @brief Translates the component
    @param translation translation vector */
    virtual void Translate( Vector3 const& translation );

    /** @brief Initializes the component's torsor
    @param withWeight intializes the force to the component weight if true
    or to 0 is false */
    virtual void InitializeForce( bool const& withWeight );

    /** @brief Updates contact map */
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

    /** @brief Returns whether a point lies inside the component
    @param pt point */
    virtual bool isIn( Point3 const& pt ) const;
    
    /** @brief Returns whether to store the contact force for post-processing 
    @param othercomp the other component involved in the contact */
    virtual bool storePPForce( Component const* othercomp ) const = 0;
    
    /** @brief Returns whether the bounding volumes of two components overlap
    @param othercomp the other component involved in the overlap test */
    virtual bool doBVolumeOverlap( Component const* othercomp ) const;    
    //@}


    /** @name Set methods */
    //@{
    /** @brief Sets the component ID number
    @param id_ component ID number */
    void setID( int const& id_ );

    /** @brief Sets the component's transformation with an 1D array of 12
    values (see class Transform for details)
    @param pos 1D array of values containing the tranformation coefficients */
    virtual void setPosition( double const* pos );

    /** @brief Sets the origin of the component's transformation
    @param centre origin coordinates as a Point3 */
    virtual void setPosition( Point3 const& centre );

    /** @brief Sets the material type
    @param mat material type */
    void setMaterial( string const& mat );

    /** @brief Initializes all contact map entries to false */
    virtual void setContactMapToFalse();
    
    /** @brief Sets contact map cumulative features to zero */
    virtual void setContactMapCumulativeFeaturesToZero();
    
    /** @brief Sets the contact map 
    @param othermap the contact map to be copied */
    void setContactMap( map< std::tuple<int,int,int>,
     	std::tuple<bool, Vector3, Vector3, Vector3> > const& othermap );
	
    /** @brief Sets the total force exerted on the component
    @param force new total force */
    void setForce( Vector3 const& force );
    
    /** @brief Sets the total torque exerted on the component
    @param torque new total torque */
    void setTorque( Vector3 const& torque ); 	         
    //@}


    /** @name Accessors */
    //@{
    /** @brief Returns the bounding box of the component */
    virtual BBox BoundingBox() const;

    /** @brief Returns the component ID number */
    int getID() const;

    /** @brief Returns a point to the rigid body with crust */
    RigidBodyWithCrust const* getRigidBody() const;

    /** @brief Returns a point to the rigid body with crust */
    RigidBodyWithCrust* getRigidBody();

    /** @brief Returns the mass of the component */
    double getMass() const;

    /** @brief Returns a pointer to the center of mass of the component */
    Point3 const* getPosition() const;

    /** @brief Copy the center of mass of the component in the array pos */
    void getPosition( double* pos ) const;

    /** @brief Returns the circumscribed radius of the rigid body */
    double getCircumscribedRadius() const;

    /** @brief Returns the radius of the sphere of volume equivalent to that of
    the rigid body */
    virtual double getEquivalentSphereRadius() const;

    /** @brief Returns the crust thickness of the rigid body */
    double getCrustThickness() const;

    /** @brief Returns the volume of the rigid body */
    virtual double getVolume() const;

    /** @brief Returns the component material name  */
    string getMaterial() const;

    /** @brief Returns a pointer to the total force exerted on the component */
    Vector3 const* getForce() const;

    /** @brief Returns a pointer to the total torque exerted on the component */
    Vector3 const* getTorque() const;

    /** @brief Returns a pointer to the torsor exerted on the component */
    virtual Torsor const* getTorsor();

    /** @brief Returns the fluid velocity interpolated at the center of mass of
    the component */
    virtual Vector3 const* getTranslationalVelocity_fluide() const;

    /** @brief Returns a pointer to the master component of the component:
    this in general and the CompositeParticle for an elementary particle */
    virtual Component* getMasterComponent();

    /** @brief Returns the particle class */
    virtual int getGeometricType() const;
    
    /** @brief Returns a pointer to the contact map */
    map< std::tuple<int,int,int>,
     	std::tuple<bool, Vector3, Vector3, Vector3> > const* getContactMap()
	const;
	
    /** @brief Returns component tag at the current discrete time (default 
    is 0) */
    virtual int getTag() const;	    
    //@}


    /** @name Methods */
    //@{
    /** @brief Adds a force exerted at a point to the torsor (torsor adds torque
    automatically)
    @param force force
    @param point point where the force is exerted 
    @param tagSecondComp tag of the other compoenent in case of a contact 
    force */
    virtual void addForce( Point3 const& point, Vector3 const& force,
    	int tagSecondComp );

    /** @brief Adds a body force exerted at the center of mass of the component
    to the torsor (associated torque is 0)
    @param force force */
    virtual void addBodyForce( Vector3 const& force );

    /** @brief Adds a torque to the torsor
    @param torque torque 
    @param tagSecondComp tag of the other compoenent in case of a contact 
    torque */
    virtual void addTorque( Vector3 const& torque, int tagSecondComp );

    /** @brief Copies the center of mass of the component in a 1D array
    @param pos 1D array where center of mass is copied
    @param i start index to copy in the 1D array */
    void copyPosition( double* pos, int i ) const;

    /** @brief Copies the component transformation in a 1D array
    @param vit 1D array where transformation is copied
    @param i start index to copy in the 1D array */
    void copyTransform( double* vit, int i ) const;

    /** @brief Copies the component transformation in a 1D array with the center
    of mass translated
    @param vit 1D array where transformation is copied
    @param i start index to copy in the 1D array
    @param vec translation vector */
    void copyTransform( double* vit, int i, Vector3 const& vec ) const;

    /** @brief Returns whether the center of mass of the component is located in
    a geometric box */
    bool isIn( BBox const& boite ) const;

    /** @brief Checks whether 2 components are in contact. If contact, computes
    contact force & torque and adds to each component torsor
    @exception ContactError in case of an excessice overlap
    @param voisin the other component
    @param dt time step magnitude
    @param time physical time
    @param LC linked-cell grid */
    virtual void InterAction( Component* voisin,
	double dt, double const& time, LinkedCell* LC );

    /** @brief Searches and stores all contact points between two components
    @exception ContactError if overlapping distance is larger than the sum of
    the crust thicknesses of the components
    @param voisin the other component
    @param dt time step magnitude
    @param time physical time
    @param LC linked-cell grid
    @param listContact list of information about contacts */
    virtual void SearchContact( Component* voisin, double dt,
      double const& time, LinkedCell *LC,
      list<ContactInfos*>& listContact );
    
    /** @brief Searches for the previous direction of contact with component
    with id = _id
    @param _id id of the component to search for */
    Vector3 lookupCollision( int _id ) const;

    /** @brief Adds a contacting component ID to the set of contacting 
    component IDs
    @param id contacting component ID */
    virtual void addContactingComponentID( int const& id );

    /** @brief Returns whether the component is a composite particle */
    virtual bool isCompositeParticle() const;

    /** @brief Returns whether the component is a composite obstacle */
    virtual bool isCompositeObstacle() const;

    /** @brief Returns whether the component is an obstacle ? (use the fact
    that all obstacles have the same very large mass by convention) */
    bool isObstacle() const;

    /** @brief Returns whether the component is an elementary particle */
    virtual bool isElementaryParticle() const;
    
    /** @brief Returns whether the component is an STL obstacle */
    virtual bool isSTLObstacle() const;    
    //@}


    /**@name I/O methods */
    //@{
    /** @brief Outputs the component for reload
    @param fileSave output stream */
    virtual void write( ostream& fileSave ) const = 0;

    /** @brief Outputs information to be transferred to the fluid
    @param fluid output stream */
    virtual void writePositionInFluid( ostream& fluid ) = 0;

    /** @brief Writes the identity of a component */
    virtual void writeIdentity( ostream& file ) const = 0;

    /** @brief Returns the number of points to write the component in a
    Paraview format */
    virtual int numberOfPoints_PARAVIEW() const = 0;

    /** @brief Returns the number of elementary polytopes to write the
    component shape in a Paraview format */
    virtual int numberOfCells_PARAVIEW() const = 0;

    /** @brief Returns a list of points describing the component in a
    Paraview format
    @param translation additional center of mass translation */
    virtual list<Point3> get_polygonsPts_PARAVIEW(
  	Vector3 const* translation = NULL ) const = 0;

    /** @brief Writes the points describing the component in a
    Paraview format
    @param f output stream
    @param translation additional center of mass translation */
    virtual void write_polygonsPts_PARAVIEW( ostream& f,
  	Vector3 const* translation = NULL ) const = 0;

    /** @brief Writes the component in a Paraview format
    @param connectivity connectivity of Paraview polytopes
    @param offsets connectivity offsets
    @param cellstype Paraview polytopes type
    @param firstpoint_globalnumber global number of the 1st point
    @param last_offset last offset used for the previous convex shape */
    virtual void write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const = 0;
    //@}

    /** @brief Writes the contact map to file in plain 2014 format
    @param fileSave output file stream */
    void writeContactMemory2014( ostream &fileSave ) const;

    /** @brief Writes the contact map to file in binary format
    @param fileOut output file stream */
    void writeContactMemory2014_binary( ostream &fileOut );

    /** @brief Reads the contact map to file in plain 2014 format
    @param fileSave input file stream */
    void readContactMap2014( istream &fileSave );

    /** @brief Reads the contact map to file in binary format
    @param fileSave input file stream */
    void readContactMap2014_binary( istream &fileSave );
    //@}


    /** @name Static methods */
    //@{    
    /** @brief Returns the number of created components  */
    static size_t getNbCreatedComponents();
    //@}


    /** @name State restoring methods */
    //@{
    /** @brief Saves the component state */
    void saveConfigState();

    /** @brief Creates and returns the component state */
    ConfigurationMemento* createConfigState();
    //@}


    /**@name Parameters */
    //@{
    static size_t m_sizeofContactMemory; /** binary size of one 
    	contact memory tuple */
    //@}
    
    
  protected:
    /** @name Constructors */
    //@{
    /** @brief Default constructor */
    Component();

    /** @brief Copy constructor 
    @param copy copied Component object */
    Component( Component const& copy );
    //@}


    /** @name I/O methods */
    //@{
    /** @brief Outputs the data related to the position of the component
    @param position output stream */
    void writePosition( ostream& position ) const;

    /** @brief Writes the component's "static" data
    @param fileOut output stream */
    virtual void writeStatic( ostream& fileOut ) const;    
    //@}


    /** @name Parameters */
    //@{
    int m_id; /**< ID number, particle ID numbers are always positive and
    	obstacles ID numbers are always negative, reference particle ID number 
	and composite obstacle ID number are always 0 */
    string m_materialName; /**< Material name */
    double m_mass; /**< Mass */
    RigidBodyWithCrust *m_geoRBWC; /**< geometric shape with crust */
    Torsor m_torsor; /**< Torsor of forces exerted on the component at its 
    	center of mass */
    ConfigurationMemento *m_memento; /**< To store the component features */
    map< std::tuple<int,int,int>,
     	std::tuple<bool, Vector3, Vector3, Vector3> > m_contactMap; /** List of 
     	active contacts with other components. It reads as follows:
    	map<tuple<own elementary particle id, neighbour id, neighbour elementary
    	particle id>, tuple<isContactActive, kt * cumulative tangential 
	dispacement, previous normal vector, kr * cumulative rotational 
	motion> > */
    CollisionHistory m_collisionHistory; /**< keep track of collisions */
    static size_t m_nb; /**< Number of created components */
    static double m_massObstacle; /**< Mass of obstacle assumed to be infinite
    */
    //@}
};

#endif
