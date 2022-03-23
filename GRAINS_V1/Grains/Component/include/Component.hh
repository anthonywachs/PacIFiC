#ifndef _COMPONENT_HH_
#define _COMPONENT_HH_

#include "Basic.hh"
#include "RigidBodyWithCrust.hh"
#include "Point3.hh"
#include "PointContact.hh"
#include "Torsor.hh"
#include <list>
#include <string>
#include <map>
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
    @author A.WACHS - 2019 - Major cleaning & refactoring */
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

    /** @brief Rotates the component using a quaternion
    @param rotation quaternion representing the rotation */
    virtual void Rotate( Quaternion const& rotation );

    /** @brief Translates the component
    @param translation translation vector */
    virtual void Translate( Vector3 const& translation );
  
    /** @brief Initializes the component's torsor
    @param withWeight intializes the force to the component weight if true 
    or to 0 is false */
    virtual void InitializeForce( bool const& withWeight ); 
    
    /** @brief Update contact map */
    virtual void updateContactMap();
  
    /** @brief Does the contact exist in the map, if yes return the pointer to 
    the cumulative tangential displacement 
    @param tangentialDepl pointer to the cumulative tangential displacement 
    @param id id number of the other component */
    virtual bool ContactInMapIsActive( double* &tangentialDepl, int const& id );
  
    /** @brief Add new contact in the map
    @param tangentialDepl initial tangential displacement 
    @param id id number of the other component */
    virtual void addNewContactInMap( double const& tangentialDepl, 
  	int const& id ); 
	
    /** @brief Increase cumulative tangential displacement with component id
    @param tangentialDepl additional tangential displacement 
    @param id id number of the other component */
    virtual void addDeplContactInMap( double const& tangentialDepl, 
  	int const& id );
		
    /** @ brief Returns whether a point lies inside the component
    @param pt point */
    virtual bool isIn( Point3 const& pt ) const;	 	   
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
  
    /** @brief Initialize all contact map entries to false */
    virtual void setContactMapToFalse();     
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

    /** @brief Returns the radius ofthe sphere of equivalent volume as the rigid
    body */
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
    //@}


    /** @name Methods */
    //@{
    /** @brief Adds a force exerted at a point to the torsor (torsor adds torque
    automatically)
    @param force force 
    @param point point where the force is exerted */
    virtual void addForce( Point3 const& point, Vector3 const& force );
  
    /** @brief Adds a body force exerted at the center of mass of the component
    to the torsor (associated torque is 0)
    @param force force */
    virtual void addBodyForce( Vector3 const& force );  
  
    /** @brief Adds a torque to the torsor 
    @param torque torque */
    virtual void addTorque( Vector3 const& torque );  
  
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
	double dt, double const& time, LinkedCell* LC ) 
	throw (ContactError) = 0;

    /** @brief Searches and stores all contact points between a composite 
    particle and a component.
    @exception ContactError if overlapping distance is larger than the sum of
    the crust thicknesses of the components    
    @param voisin the other component
    @param dt time step magnitude
    @param time physical time
    @param LC linked-cell grid
    @param listContact list of information about contacts */
    virtual void SearchContact( Component* voisin, double dt,
      double const& time, LinkedCell *LC,
      list<ContactInfos*>& listContact ) throw (ContactError);

    /** @brief Increments the coordination number by nc
    @param nc increment of the coordination number */
    virtual void addToCoordinationNumber( int const& nc );

    /** @brief Returns whether the component is a composite particle */
    virtual bool isCompositeParticle() const;
  
    /** @brief Returns whether the component is a composite obstacle */
    virtual bool isCompositeObstacle() const;

    /** @brief Returns whether the component is an obstacle ? (use the fact
    that obstacle have a zero mass by convention) */
    bool isObstacle() const;

    /** @brief Returns whether the component is an elementary particle */
    virtual bool isElementaryParticle() const;
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


    /** @name Static methods */
    //@{
    /** @brief Resets the number of created components to nb_
    @param nb_ number of created components */
    static void setNbCreatedComponents( const int &nb_ );
  
    /** @brief Returns the number of created components  */
    static int getNbCreatedComponents();
    //@}


    /** @name State restoring methods */
    //@{
    /** @brief Saves the component state */
    void saveConfigState();

    /** @brief Creates and returns the component state */
    ConfigurationMemento* createConfigState(); 
    //@}   
  
  

  protected:
    /** @name Constructors */
    //@{
    /** @brief Default constructor
    @param autonumbering whether to increase the number of created components or
    not */
    Component( bool const& autonumbering = true );

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
    int m_id; /**< ID number */
    string m_materialName; /**< Material name */
    double m_mass; /**< Mass */
    RigidBodyWithCrust *m_geoRBWC; /**< geometric shape with crust */
    Torsor m_torsor; /**< Torsor of forces exerted on the component at its center
    	of mass */
    ConfigurationMemento *m_memento; /**< To store the component features */
    map< int, pair<double,bool> > m_contactMap; /**< List of active contacts 
    	with other components */
    static int m_nb; /**< Number of created components */
    //@}
};

#endif

