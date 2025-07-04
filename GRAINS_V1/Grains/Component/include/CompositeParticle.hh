#ifndef _COMPOSITEPARTICLE_HH_
#define _COMPOSITEPARTICLE_HH_

#include "Particle.hh"
#include "RigidBody.hh"
#include "PointContact.hh"


/** @brief The class CompositeParticle.

    A freely moving particle made of multiple glued particles.

    @author D.RAKOTONIRINA - IFP Energies nouvelles - 2014 - Creation
    @author A.WACHS - 2021 - Major cleaning & refactoring */
// ============================================================================
class CompositeParticle : public Particle
{
  public:
    /**@name Constructors & Destructor */
    //@{
    /** @brief Constructor with autonumbering as input parameter
    @param autonumbering whether to increment the component indexing */
    CompositeParticle( bool const& autonumbering );

    /** @brief Constructor with an XML node as an input parameter. This
    constructor is expected to be used for reference composite particles.
    Autonumbering is set to false
    @param root XML node
    @param pc particle class */
    CompositeParticle( DOMNode* root, int const& pc );

    /** @brief Constructor with input parameters. Autonumbering
    is set to false and numbering is set with the parameter id_
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
    CompositeParticle( int const& id_, Particle const* ParticleRef,
	double const& vx, double const& vy, double const& vz,
	double const& qrotationx, double const& qrotationy,
	double const& qrotationz, double const& qrotations,
	double const& rx, double const& ry, double const& rz,
	const double m[12],
	ParticleActivity const& activ,
	int const& tag_,
	int const& coordination_number_ = 0 );

    /** @brief Constructor with input parameters. This constructor is expected
    to be used for periodic clone particle. Autonumbering
    is set to false and numbering is set with the parameter id_
    @param id_ ID number
    @param ParticleRef reference particle
    @param vtrans translational velocity
    @param vrot angular velocity
    @param qrot rotation quaternion
    @param config particle transformation
    @param activ particle activity 
    @param contactMap contact map */
    CompositeParticle( int const& id_, Particle const* ParticleRef,
	Vector3 const& vtrans,
	Quaternion const& qrot,
	Vector3 const& vrot,
	Transform const& config,
	ParticleActivity const& activ,
     	map< std::tuple<int,int,int>,
     	std::tuple<bool, Vector3, Vector3, Vector3> > const* contactMap );

    /** @brief Copy constructor (the torsor is initialized to 0)
    @param other copied CompositeParticle object
    @param autonumbering whether to increment the component indexing */
    CompositeParticle( CompositeParticle const& other, 
    	bool const& autonumbering );

    /** @brief Destructor */
    virtual ~CompositeParticle();
    //@}


    /**@name Methods */
    //@{
    /** @brief Returns whether the component is a composite particle
    (here true) */
    bool isCompositeParticle() const;

    /** @brief Creates a clone of the composite particle. This method calls
    the standard copy constructor and is used for new composite particles to be
    inserted in the simulation. Activity is set to WAIT. The calling
    object is expected to be a reference composite particle 
    @param autonumbering whether to increment the component indexing */
    Particle* createCloneCopy( bool const& autonumbering ) const ;

    /** @brief Creates a clone of the composite particle. This method calls the
    constructor CompositeParticle( int const& id_, Particle const* ParticleRef,
    Vector3 const& vtrans, Quaternion const& qrot, Vector3 const& vrot,
    Transform const& config, ParticleActivity const& activ ) and is used for
    periodic clone composite particles to be inserted in the simulation.
    Autonumbering is set to false and numbering is set with the parameter id_
    @param id_ ID number
    @param ParticleRef reference particle
    @param vtrans translational velocity
    @param vrot angular velocity
    @param qrot rotation quaternion
    @param config particle transformation
    @param activ particle activity 
    @param contactMap contact map */
    Particle* createCloneCopy( int const& id_,
    	Particle const* ParticleRef, Vector3 const& vtrans,
	Quaternion const& qrot,	Vector3 const& vrot,
	Transform const& config, ParticleActivity const& activ,
	map< std::tuple<int,int,int>,
     	std::tuple<bool, Vector3, Vector3, Vector3> > const* contactMap ) 
	const ;

    /** @brief Sets the boolean that tells that the rigid body's transformation
    with the scaling by the crust thickness to shrink the rigid bodies has
    already been computed to false */
    void initialize_transformWithCrust_to_notComputed();

    /** @brief Returns whether there is geometric contact with another
    component. Note: the other component must not be of the derived type
    CompositeObstacle
    @param voisin the other component */
    bool isContact( Component const* voisin ) const;

    /** @brief Returns whether there is geometric contact with another
    component accounting for crust thickness. Note: the other component must
    not be of the derived type CompositeObstacle
    @param voisin the other component */
    bool isContactWithCrust( Component const* voisin ) const;

    /** @brief Returns whether there is geometric proximity with another
    component in the sense of whether their respective bounding boxes overlap.
    Note: the other component must not be of the derived type CompositeObstacle
    @param voisin the other component */
    bool isClose( Component const* voisin ) const;

    /** @brief Returns whether there is geometric proximity with another
    component in the sense of whether their respective bounding boxes minus
    their crust thickness overlap. Note: the other component must not be of the
    derived type CompositeObstacle
    @param voisin the other component */
    bool isCloseWithCrust( Component const* voisin ) const;

    /** @brief Contact between a composite particle and a component. If contact
    exists, computes the contact force and torque and adds to each component
    @exception ContactError if overlapping distance is larger than the sum of
    the crust thicknesses of the components
    @param voisin the other component
    @param dt time step magnitude
    @param time physical time
    @param LC linked-cell grid */
    void InterAction( Component* voisin,
	double dt, double const& time, LinkedCell* LC );

    /** @brief Searches and stores all contact points between two components
    @exception ContactError if overlapping distance is larger than the sum of
    the crust thicknesses of the components
    @param voisin the other component
    @param dt time step magnitude
    @param time physical time
    @param LC linked-cell grid
    @param listContact list of information about contacts */
    void SearchContact( Component* voisin, double dt,
      double const& time, LinkedCell *LC,
      list<ContactInfos*>& listContact );

    /** @brief Compose the component transformation on the left by another
    transformation: this = t o this (this first followed by t)
    @param t the other affine transformation */
    void composePositionLeftByTransform( Transform const& t );

    /** @brief Compose the component transformation on the right by another
    transformation: this = this o t (t first followed by this)
    @param t the other affine transformation */
    void composePositionRightByTransform( Transform const& t );

    /** @brief Solves the Newton's law and move particle to their new position
    @exception MotionError motion is larger than crust thickness
    @param time physical time
    @param dt_particle_vel velocity time step magnitude 
    @param dt_particle_disp motion time step magnitude */
    void Move( double time, double const& dt_particle_vel, 
    	double const& dt_particle_disp );

    /** @brief Translates the composite particle
    @param translation translation vector */
    void Translate( Vector3 const& translation );

    /** @ brief Returns whether a point lies inside the composite particle
    @param pt point */
    bool isIn( Point3 const& pt ) const;

    /** @brief Returns whether two particles are of the same type
    @param other the other particle */
    virtual bool equalType( Particle const* other ) const;    
    //@}


    /**@name Methods I/O */
    //@{
    /** @brief Returns the number of points to write the composite particle in a
    Paraview format */
    virtual int numberOfPoints_PARAVIEW() const;

    /** @brief Returns the number of elementary polytopes to write the
    composite particle shape in a Paraview format */
    virtual int numberOfCells_PARAVIEW() const;

    /** @brief Writes the points describing the composite particle in a
    Paraview format
    @param f output stream
    @param translation additional center of mass translation */
    virtual void write_polygonsPts_PARAVIEW( ostream& f,
  	Vector3 const* translation = NULL ) const;

    /** @brief Returns a list of points describing the component in a
    Paraview format
    @param translation additional center of mass translation */
    virtual list<Point3> get_polygonsPts_PARAVIEW(
  	Vector3 const* translation = NULL ) const;

    /** @brief Writes the composite particle in a Paraview format
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
    
    /** @brief Writes the particle in an OBJ format
    @param f output stream 
    @param group_number elementary particle number in case of composite 
    particle 
    @param firstpoint_number number of the 1st point */
    virtual void write_OBJ( ostream& f, size_t const& group_number,
    	size_t& firstpoint_number ) const;  
    //@}


    /**@name Accessors */
    //@{
    /** @brief Returns the bounding box of the composite particle */
    BBox BoundingBox() const;

    /** @brief Returns the number of elementary particles */
    size_t getNbElemPart() const;

    /** @brief Returns the number of corners of the rigib body shape and a code
    describing the rigid body shape */
    virtual int getNbCorners() const;

    /** @brief Returns the volume of the composite particle */
    double getVolume() const;

    /** @brief Returns a pointer to the vector of initial relative positions */
    vector<Vector3> const* getInitialRelativePositions() const;
    
    /** @brief Returns the radius of the sphere of volume equivalent to that of 
    the composite particle */ 
    double getEquivalentSphereRadius() const;    
    //@}


    /**@name Methods set */
    //@{
    /** @brief Sets the angular velocity
    @param vrot angular velocity */
    void setAngularVelocity( Vector3 const& vrot );

    /** @brief Sets the angular velocity
    @param omx x-angular velocity component 
    @param omy y-angular velocity component     
    @param omz z-angular velocity component */
    void setAngularVelocity( double const& omx, double const& omy,
	double const& omz ); 

    /** @brief Sets the translation velocity
    @param vtrans translation velocity */
    void setTranslationalVelocity( Vector3 const& vtrans );

    /** @brief Sets the translation velocity
    @param vx x-translation velocity component 
    @param vy y-translation velocity component     
    @param vz z-translation velocity component */
    void setTranslationalVelocity( double const& vx, double const& vy,
	double const& vz ); 

    /** @brief Sets the origin of the composite particle's transformation
    @param centre origin coordinates as a Point3 */
    void setPosition( Point3 const& centre );

    /** @brief Sets the composite particle's transformation with an 1D array
    of 12 values (see class Transform for details)
    @param pos 1D array of values containing the tranformation coefficients */
    void setPosition( double const* pos );

    /** @brief Sets the composite particle's transformation with a
    transformation
    @param transform_ transformation */
    void setTransform( Transform const& transform_ );
    
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
    //@}


    /**@name I/O methods */
    //@{
    /** @brief Reads composite particle data from a stream. Usage: for standard
    composite particles in the 2014 reload format
    @param fileIn input stream
    @param referenceParticles reference particles for each class of
    particles */
    virtual void read2014( istream& fileIn, vector<Particle*> const*
  	referenceParticles );

    /** @brief Reads composite particle data from a stream in a binary form.
    Usage: for standard composite particles in the 2014 reload format
    @param fileIn input stream
    @param referenceParticles reference particles for each class of
    particles */
    virtual void read2014_binary( istream& fileIn, vector<Particle*> const*
  	referenceParticles );
    //@}


  protected:
    /** @name Parameters */
    //@{
    size_t m_nbElemPart; /**< Number of elementary (glued) particles */
    vector<Particle*> m_elementaryParticles; /**< vector of elementary particles
    	of the composite */
    vector<Vector3> m_InitialRelativePositions; /**< vector of relative
    	positions of elementary particles with respect to (0,0,0) required to
	be the center of mass position of the composite in its initial
	reference position */
    vector<Matrix> m_InitialRotationMatrices; /**< vector of initial rotation
    	matrices of elementary particles in their initial reference position */
    //@}


    /**@name I/O methods */
    //@{
    /** @brief Saves additional features of a (in practice reference) composite
    particle for reload
    @param fileSave output stream */
    virtual void writeAdditionalFeatures( ostream& fileSave ) const;

    /** @brief Reads additional features of a (in practice reference) particle
    data from a stream
    @param fileIn input stream */
    virtual void readAdditionalFeatures( istream& fileIn );
    //@}


    /** @name Methods */
    //@{
    /** @brief Computes and sets the moment of inertia tensor and its inverse */
    void BuildInertia();

    /** @brief Computes and sets the circumscribed radius */
    virtual void setCircumscribedRadius();

    /** @brief Sets the elementary particle position given the position of the
    composite particle */
    void setElementaryParticlesPosition();

    /** @brief Creates and sets the elementary particles
    @param CompParticleRef reference composite particle */
    void createSetElementaryParticles(
    	CompositeParticle const* CompParticleRef );
	
    /** @brief Creates the non-spherical bounding volume */
    void createBoundingVolume();	
    //@}


  private:
    /**@name Constructors */
    //@{
    /** @brief Default constructor */
    CompositeParticle();
    //@}
};

#endif
