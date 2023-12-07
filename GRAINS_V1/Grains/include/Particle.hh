#ifndef _PARTICLE_HH_
#define _PARTICLE_HH_

#include "Component.hh"
#include "ParticleKinematics.hh"
#include "PointContact.hh"
#include "Vector3.hh"
#include "Cell.hh"
#include "ReaderXML.hh"
#include "WriterXML.hh"
using namespace solid;


class Cell;
class Obstacle;
class SimpleObstacle;
class PeriodicObstacle;
class Quaternion;
class LinkedCell;
struct PointForcePostProcessing;
ostream& operator << ( ostream& f, Particle const& P );


/** @brief Particle activity */
enum ParticleActivity
  {
    /** @brief Not inserted yet */
    WAIT,
    /** @brief Active in the simulation */
    COMPUTE,
    /** @brief Extracted from simulation and wait */
    CLEARandWAIT
  };


/** @brief Data to compute part of the particle acceleration explicitly.
Note: generally used for particles lighther than the fluid */
struct VelocityInfosNm1
{
  Vector3 TranslationalVelocity_nm1; /**< Translation velocity at the previous
  	fluid discrete time */
  Vector3 TranslationalVelocity_difference; /**< Translational velocity
  	difference at the previous fluid discrete time */
  Vector3 RotationalVelocity_nm1; /**< Angular velocity at the previous
  	fluid discrete time */
  Vector3 RotationalVelocity_difference; /**< Angular velocity
  	difference at the previous fluid discrete time */
};


/** @brief The class Particle.

    A freely moving particle.

    @author F.PRADEL - Institut Francais du Petrole - 1999 - Creation
    @author G.FERRER - Institut Francais du Petrole - 2000 - Creation
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class Particle : public Component
{
  public:
    /**@name Constructors & Destructor */
    //@{
    /** @brief Constructor with autonumbering as input parameter
    @param autonumbering whether to increment the component indexing */
    Particle( bool const& autonumbering = true );

    /** @brief Constructor with an XML node as an input parameter. This
    constructor is expected to be used for reference particles
    @param root XML node
    @param autonumbering whether to increment the component indexing
    @param pc particle class */
    Particle( DOMNode* root, bool const& autonumbering = true,
  	int const& pc = 0 );
	
    /** @brief Constructor with input parameters. This
    constructor is expected to be used for reference particles
    @param georbwc pointer to a rigid body with crust object
    @param density particle density
    @param mat particle material
    @param autonumbering whether to increment the component indexing
    @param pc particle class */
    Particle( RigidBodyWithCrust* georbwc, double const& density,
    	string const& mat,
    	bool const& autonumbering = true,
  	int const& pc = 0 );	

    /** @brief Constructor with input parameters
    @param id_ ID number
    @param ParticleRef reference particle
    @param vx x translational velocity component
    @param vy y translational velocity component
    @param vz z translational velocity component
    @param rx x angular velocity component
    @param ry y angular velocity component
    @param rz z angular velocity component
    @param qrotationx x rotation quaternion component
    @param qrotationy y rotation quaternion component
    @param qrotationz z rotation quaternion component
    @param qrotations scalar rotation quaternion component
    @param m particle position & configuration as a 1D array
    @param activ particle activity
    @param tag_ tag of the cell the particle belongs to
    @param coordination_number_ particle coordination number */
    Particle( int const& id_, Particle const* ParticleRef,
	double const& vx, double const& vy, double const& vz,
	double const& qrotationx, double const& qrotationy,
	double const& qrotationz, double const& qrotations,
	double const& rx, double const& ry, double const& rz,
	const double m[12],
	ParticleActivity const& activ,
	int const& tag_,
	int const& coordination_number_ = 0 );

    /** @brief Constructor with input parameters. This constructor is expected
    to be used for periodic clone particle
    @param id_ ID number
    @param ParticleRef reference particle
    @param vtrans translational velocity
    @param vrot angular velocity
    @param qrot rotation quaternion
    @param config particle transformation
    @param activ particle activity */
    Particle( int const& id_, Particle const* ParticleRef,
	Vector3 const& vtrans,
	Quaternion const& qrot,
	Vector3 const& vrot,
	Transform const& config,
	ParticleActivity const& activ );

    /** @brief Copy constructor (the torsor is initialized to 0)
    @param other copied Particle object */
    Particle( Particle const& other );

    /** @brief Destructor */
    virtual ~Particle();
    //@}


    /** @name Static methods */
    //@{
    /** @brief Returns the fluid density */
    static double getFluidDensity();

    /** @brief Sets the fluid density
    @param rho fluid density */
    static void setFluidDensity( double rho );

    /** @brief Defines whether the particle acceleration (i.e. change of
    momentum) is corrected by the fluid density in case of immersed rigid
    bodies
    @param correct true if fluid corrected */
    static void setFluidCorrectedAcceleration( bool correct );

    /** @brief Returns whether the particle acceleration (i.e. change of
    momentum) is corrected by the fluid density in case of immersed rigid
    bodies */
    static bool getFluidCorrectedAcceleration();

    /** @brief Returns the viscosity of the surrounding fluid */
    static double getFluidViscosity();

    /** @brief Sets the viscosity of the surrounding fluid
    @param mu fluid viscosity */
    static void setFluidViscosity( double mu );
    //@}


    /** @name Set methods */
    //@{
    /** @brief Resets kinematics and transformation to 0 */
    void reset();

    /** @brief Resets kinematics to 0 */
    void resetKinematics();

    /** @brief Sets the angular velocity
    @param vrot angular velocity */
    virtual void setAngularVelocity( Vector3 const& vrot );

    /** @brief Sets the translation velocity
    @param vtrans translation velocity */
    virtual void setTranslationalVelocity( Vector3 const& vtrans );

    /** @brief Sets particle activity
    @param activity particle activity */
    void setActivity( ParticleActivity activity );

    /** @brief Updates velocity difference and velocity at previous discrete
    time */
    void setVelocityAndVelocityDifferencePreviousTime();

    /** @brief Sets the velocity at the previous discrete time when a simulation
    is restarted (as Grains3D does not save this information)
    @param vx x translational velocity component
    @param vy y translational velocity component
    @param vz z translational velocity component
    @param omx x angular velocity component
    @param omy y angular velocity component
    @param omz z angular velocity component */
    void setVelocityPreviousTimeRestart(
      double const& vx, double const& vy, double const& vz,
      double const& omx, double const& omy, double const& omz );

    /** @brief Sets and returns the particle tag
    @param tag_ tag */
    int setTag( int tag_ );

    /** @brief Sets the geographic position of the particle in the LinkedCell
    @param geoloc_ geographic position of the particle in the LinkedCell */
    void setGeoPosition( GeoPosition const& geoloc_ );

    /** @brief Sets the cell the particle belonged to at the previous discrete
    time
    @param cel the cell the particle belonged to at the previous discrete
    time */
    void setCellNm1( Cell* cel );

    /** @brief Sets the cell the particle belongs to, the particle tag and
    the geographic location of the particle at the current time
    @param cell_ the cell containing the particle
    @param tag_ tag
    @param geoloc_ geographic position of the particle in the LinkedCell */
    void setCellTagGeoPosition( Cell* cell_, int tag_,
    	GeoPosition const& geoloc_ );

    /** @brief Sets the cell the particle belonged to, the particle tag and
    the geographic location of the particle at the previous time
    @param cell_ the cell containing the particle
    @param tag_ tag
    @param geoloc_ geographic position of the particle in the LinkedCell */
    void setCellTagGeoPosition_nm1( Cell* cell_, int tag_,
    	GeoPosition const& geoloc_ );

    /** @brief Sets the rotation quaternion
    @param vecteur0 x component of the quaternion
    @param vecteur1 y component of the quaternion
    @param vecteur2 z component of the quaternion
    @param scalaire scalar component of the quaternion */
    void setQuaternionRotation( double const& vecteur0,
	double const& vecteur1,
	double const& vecteur2,
	double const& scalaire );

    /** @brief Sets the rotation quaternion
    @param qrot rotation quaternion */
    void setQuaternionRotation( Quaternion const& qrot );

    /** @brief Sets the particle geometric type
    @param pc particle class */
    void setGeometricType( int const& pc );

    /** @brief Sets kinematics at time t-2dt from a 1D array of 12 scalars
    (translational velocity, angular velocity, variation of translational
    velocity, variation of angular velocity)
    @param tab 1D array of 4 vectors containing translational velocity, angular
    velocity, variation of translational velocity and variation of angular
    velocity */
    void setKinematicsNm2( double const* tab );

    /** @brief Sets the velocity with a random motion
    @param coefTrans translational random velocity amplitude
    @param coefRot angular random velocity amplitude */
    void setRandomMotion( double const& coefTrans,
	double const& coefRot );

    /** @brief Sets the pointer to the master particle of the particle (
    in general the CompositeParticle for an elementary particle)
    @param master_ pointer to the master particle */
    void setMasterParticle( Particle* master_ );

    /** @brief Sets the particle density */
    void setDensity( double const& density_ );

    /** @brief Sets the particle's transformation with a transformation
    @param transform_ transformation */
    virtual void setTransform( Transform const& transform_ );

    /** @brief Sets the pointer to the particle's kinematics
    @param pkine transformation */
    void setKinematics( ParticleKinematics* pkine );
    //@}


    /** @name Methods */
    //@{
    /** @brief Solves the Newton's law and move particle to their new position
    @exception DisplacementError displacement is larger than crust thickness
    @param time physical time
    @param dt time step magnitude */
    virtual void Move( double time, double dt );

    /** @brief Contact between a particle and a component. If contact exists,
    computes the contact force and torque and adds to each component
    @exception ContactError if overlapping distance is larger than the sum of
    the crust thicknesses of the components
    @param voisin the other component
    @param dt time step magnitude
    @param time physical time
    @param LC linked-cell grid */
    virtual void InterAction( Component* voisin,
	double dt, double const& time, LinkedCell* LC );

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
      list<ContactInfos*>& listContact );

    /** @brief Copies angular velocity in a 1D array
    @param vit 1D array where angular velocity is copied
    @param i start index to copy in the 1D array */
    void copyAngularVelocity( double* vit, int i ) const;

    /** @brief Copies translational velocity in a 1D array
    @param vit 1D array where translational velocity is copied
    @param i start index to copy in the 1D array */
    void copyTranslationalVelocity( double* vit, int i ) const;

    /** @brief Copies rotation quaternion in a 1D array
    @param vit 1D array where rotation quaternion is copied
    @param i start index to copy in the 1D array */
    void copyQuaternionRotation( double* vit, int i ) const;

    /** @brief Copies force and torque in a 1D array
    @param fm 1D array where force and torque are copied
    @param i start index to copy in the 1D array */
    void copyForceTorque( double* fm, int i ) const;

    /** @brief Copies kinematics at time t-2dt (translational velocity, angular
    velocity, variation of translational velocity, variation of angular
    velocity) in a 1D array
    @param vit 1D array where kinematics at time t-2dt is copied
    @param i start index to copy in the 1D array */
    void copyKinematicsNm2( double* vit, int i ) const;

    /** @brief Adds a force whose point of application is different from the
    reference point of the torsor (additional torque contribution)
    @param f_ the added force
    @param point point of application of the force */
    void addForce( Point3 const& point, Vector3 const& f_ );

    /** @brief Updates geographic localisation in the LinkedCell. Note that
    this method uses the cell from the previous time m_cellule_nm1 */
    void updateGeoPosition();

    /** @brief Creates the VelocityInfosNm1 structure */
    void createVelocityInfosNm1();

    /** @brief Returns an orientation vector to describe the angular position of
    the particle */
    Vector3 computeOrientationVector() const;

    /** @brief Initializes the particle torsor
    @param withWeight intializes the force to the particle weight if true
    or to 0 is false */
    void InitializeForce( bool const& withWeight );

    /** @brief Computes particle net weight */
    void computeWeight();

    /** @brief Creates a clone of the particle. This method calls the standard
    copy constructor and is used for new particles to be inserted in the
    simulation. Numbering is automatic, total number of components is
    incremented by 1 and activity is set to WAIT. The calling object is
    expected to be a reference particle */
    virtual Particle* createCloneCopy() const ;

    /** @brief Creates a clone of the particle. This method calls the
    constructor Particle( int const& id_, Particle const* ParticleRef, Vector3
    const& vtrans, Quaternion const& qrot, Vector3 const& vrot,	Transform
    const& config, ParticleActivity const& activ ) and is used for periodic
    clone particles to be inserted in the simulation. Numbering is set with the
    parameter id_ and total number of components left unchanged.
    @param id_ ID number
    @param ParticleRef reference particle
    @param vtrans translational velocity
    @param vrot angular velocity
    @param qrot rotation quaternion
    @param config particle transformation
    @param activ particle activity */
    virtual Particle* createCloneCopy( int const& id_,
    	Particle const* ParticleRef, Vector3 const& vtrans,
	Quaternion const& qrot,	Vector3 const& vrot,
	Transform const& config, ParticleActivity const& activ ) const ;

    /** @brief Sets the boolean that tells that the rigid body's transformation
    with the scaling by the crust thickness to shrink the rigid bodies has
    already been computed to false */
    virtual void initialize_transformWithCrust_to_notComputed();

    /** @brief Copies the cell the particle belonged to, the particle tag and
    the geographic location of the particle from current time to previous
    time */
    void copyCellTagGeoPosition_n_to_nm1();

    /** @brief Increments the coordination number by nc
    @param nc increment of the coordination number */
    virtual void addToCoordinationNumber( int const& nc );

    /** @brief Returns whether the particle is an elementary particle */
    virtual bool isElementaryParticle() const;

    /** @brief Compose the component transformation on the left by another
    transformation: this = t o this (this first followed by t)
    @param t the other affine transformation */
    virtual void composePositionLeftByTransform( Transform const& t );

    /** @brief Compose the component transformation on the right by another
    transformation: this = this o t (t first followed by this)
    @param t the other affine transformation */
    virtual void composePositionRightByTransform( Transform const& t );

    /** @brief Computes particle inertia tensor in the space fixed coordinate
    frame
    @param inertia inertia tensor in the space fixed coordinate frame */
    void computeInertiaTensorSpaceFixed( vector<double>& inertia ) const;
    //@}


    /** @name Accessors */
    //@{
    /** @brief Returns a pointer to the particle kinematics */
    ParticleKinematics const* getKinematics() const;

    /** @brief Returns particle activity */
    ParticleActivity getActivity() const;

    /** @brief Returns particle tag at the current discrete time */
    int getTag() const;

    /** @brief Returns particle tag at the previous discrete time */
    int getTagNm1() const;

    /** @brief Returns translational velocity difference at the previous
    discrete time */
    Vector3 getTranslationalVelocityDifferencePreviousTime() const;

    /** @brief Returns angular velocity difference at the previous
    discrete time */
    Vector3 getRotationalVelocityDifferencePreviousTime() const;

    /** @brief Returns particle inertia tensor in the body fixed coordinate
    frame */
    double const* getInertiaTensorBodyFixed() const;

    /** @brief Returns inverse of particle inertia tensor in the body fixed
    coordinate frame */
    double const* getInverseInertiaTensorBodyFixed() const;

    /** @brief Returns particle density */
    double getDensity() const;

    /** @brief Returns the radius of the sphere of same volume */
    double getEquivalentSphereRadius() const;

    /** @brief Returns the velocity at a point in space based on the
    translational and angular velocity of the particle. This method assumes
    that the point belongs to the particle but this assumption is not verified.
    @param pt the point at which velocity is computed */
    Vector3 getVelocityAtPoint( Point3 const& pt ) const;

    /** @brief Returns the rotation quaternion */
    Quaternion const* getQuaternionRotation() const;

    /** @brief Returns angular velocity */
    virtual Vector3 const* getAngularVelocity() const;

    /** @brief Returns translational velocity */
    virtual Vector3 const* getTranslationalVelocity() const;

    /** @brief Returns total force exerted on the particle */
    Vector3 const* getForce() const;

    /** @brief Returns the cell the particle belonged to at the previous
    discrete time */
    Cell* getCellNm1() const;

    /** @brief Returns the cell the particle belongs to at the current
    discrete time */
    Cell* getCell() const;

    /** @brief Returns particle class */
    int getGeometricType() const;

    /** @brief Returns particle geographic localisation at the current
    discrete time */
    GeoPosition getGeoPosition() const;

    /** @brief Returns particle geographic localisation at the previous
    discrete time */
    GeoPosition getGeoPositionNm1() const;

    /** @brief Returns particle coordination number */
    virtual int getCoordinationNumber() const;

    /** @brief Returns the number of corners of the rigib body shape and a code
    describing the rigid body shape */
    virtual int getNbCorners() const;

    /** @brief Returns a pointer to the reference component of the component:
    this in general and the CompositeParticle for an elementary particle */
    virtual Component* getMasterComponent();
    //@}


    /** @name State storing/restoring methods */
    //@{
    /** @brief Saves particle state */
    void saveState();

    /** @brief Creates and returns particle state */
    pair<ConfigurationMemento*,ParticleKinematicsMemento*> createState();

    /** @brief Restores particle state */
    void restoreState();
    //@}


    /**@name I/O methods */
    //@{
    /** @brief Reads a (in practice reference) particle data from a stream
    @param fileIn input stream
    @param elemPart true if the particle is an elementary particle of a
    composite particle, false otherwise */
    void read( istream& fileIn, bool elemPart = false );

    /** @brief Reads particle data from a stream. Usage: for standard particles
    in the 2014 reload format
    @param fileIn input stream
    @param referenceParticles reference particles for each class of
    particles */
    virtual void read2014( istream& fileIn, vector<Particle*> const*
  	referenceParticles );

    /** @brief Reads particle data from a stream in a binary form. Usage: for
    standard particles in the 2014 reload format
    @param fileIn input stream
    @param referenceParticles reference particles for each class of
    particles */
    virtual void read2014_binary( istream& fileIn, vector<Particle*> const*
  	referenceParticles );

    /** @brief Saves a (in practice reference) particle for reload
    @param fileSave output stream */
    void write( ostream& fileSave ) const;

    /** @brief Saves particle for reload in 2014 format
    @param fileSave output stream */
    void write2014( ostream& fileSave ) const;

    /** @brief Saves particle for reload in 2014 and binary format
    @param fileSave output stream */
    void write2014_binary( ostream& fileSave );

    /** @brief Saves the identity of the particle i.e its id number
    @param file output stream */
    void writeIdentity( ostream& file ) const;

    /** @brief Writes the particle's "static" data
    @param fileOut output stream */
    void writeStatic( ostream& fileOut ) const;

    /**  @brief Outputs information to be transferred to the fluid
    @param fluid output stream */
    virtual void writePositionInFluid( ostream& fluid );

    /** @brief Output operator
    @param f output stream
    @param P Particle object */
    friend ostream& operator << ( ostream& f, Particle const& P );

    /** @brief Writes the particle in a Paraview format
    @param connectivity connectivity of Paraview polytopes
    @param offsets connectivity offsets
    @param cellstype Paraview polytopes type
    @param firstpoint_globalnumber global number of the 1st point
    @param last_offset last offset used for the previous convex shape */
    virtual void write_polygonsStr_PARAVIEW(list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset) const ;

    /** @brief Writes the points describing the particle in a
    Paraview format with a transformation that may be different than the current
    transformation of the particle
    @param f output stream
    @param transform transformation
    @param translation additional center of mass translation */
    virtual void write_polygonsPts_PARAVIEW( ostream& f,
  	Transform const& transform,
  	Vector3 const* translation = NULL ) const;

    /** @brief Returns the number of points to write the particle in a
    Paraview format */
    virtual int numberOfPoints_PARAVIEW() const ;

    /** @brief Returns the number of elementary polytopes to write the particle
    shape in a Paraview format */
    virtual int numberOfCells_PARAVIEW() const;

    /** @brief Returns a list of points describing the component in a
    Paraview format
    @param translation additional center of mass translation */
    virtual list<Point3> get_polygonsPts_PARAVIEW(
  	Vector3 const* translation = NULL ) const ;

    /** @brief Writes the points describing the particle in a
    Paraview format
    @param f output stream
    @param translation additional center of mass translation */
    virtual void write_polygonsPts_PARAVIEW( ostream& f,
  	Vector3 const* translation = NULL ) const ;
    //@}


  protected:
    /**@name Parameters */
    //@{
    Particle* m_masterParticle; /**< master particle (the particle itself in
    	general and the actual master for elementary particles of a
	composite particle */
    ParticleKinematics* m_kinematics; /**< particle kinematics */
    double m_density; /**< Density */
    static double m_fluidDensity; /**< Surrounding fluid density */
    static double m_fluidViscosity; /**< Surrounding fluid viscosity */
    static bool m_fluidCorrectedAcceleration; /**< Whether the particle
    	acceleration is corrected by the fluid density */
    static bool m_splitExplicitAcceleration; /**< Whether part of the particle
    	acceleration is computed explicitly */
    double m_inertia[6]; /**< Inertia tensor I={I(1,1), I(1,2), I(1,3),
  	I(2,2), I(2,3), I(3,3)} */
    double m_inertia_1[6]; /**< Inverse inertia tensor */
    ParticleActivity m_activity; /**< Particle activity */
    struct VelocityInfosNm1* m_VelocityInfosNm1; /**< data to compute the
    	contribution of the particle acceleration treated explicitly (used for
	neutrally buoyant or lighter particles than the fluid) */
    int m_tag; /**< tag of the cell the particle belongs to at the
    	current time: 0=interior, 1=buffer zone, 2=halo zone */
    GeoPosition m_GeoLoc; /**< geographic position of the particle in the
    	Linked cell, i.e. geographic position of the cell the particle belongs
	to at the current time */
    Cell* m_cellule; /**< Cell that the particle belongs to at the
    	current time */
    int m_tag_nm1; /**< tag of the cell the particle belonged to at the
    	previous time: 0=interior, 1=buffer zone, 2=halo zone */
    GeoPosition m_GeoLoc_nm1; /**< geographic position of the particle in the
    	Linked cell, i.e. geographic position of the cell the particle belonged
	to at the previous time */
    Cell* m_cellule_nm1; /**< Cell that the particle belonged to at the
    	previous time */
    int m_GeomType; /**< particle geometric type */
    int m_coordination_number; /**< coordination number */
    Vector3 m_weight; /**< particle weight */
    string m_specific_composite_shape; /**< specific composite particle 
    	shape, e.g., SpheroCylinder */  
    //@}


    /**@name I/O methods */
    //@{
    /** @brief Saves additional features of a (in practice reference) particle
    for reload
    @param fileSave output stream */
    virtual void writeAdditionalFeatures( ostream& fileSave ) const;

    /** @brief Reads additional features of a (in practice reference) particle
    data from a stream
    @param fileIn input stream */
    virtual void readAdditionalFeatures( istream& fileIn );
    //@}


  private:
    /**@name Constructors  */
    //@{
    /** @brief Default constructor */
    Particle();
    //@}
};

#endif
