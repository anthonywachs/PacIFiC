#ifndef _SPHEROCYLINDER_HH_
#define _SPHEROCYLINDER_HH_

#include "CompositeParticle.hh"


/** @brief The class SpheroCylinder.

    A freely moving sphero-cylinder.

    @author A.WACHS - 2023 - Creation */
// ============================================================================
class SpheroCylinder : public CompositeParticle
{
  public:
    /**@name Constructors & Destructor */
    //@{
    /** @brief Constructor with autonumbering as input parameter
    @param autonumbering whether to increment the component indexing */
    SpheroCylinder( bool const& autonumbering = true );

    /** @brief Constructor with an XML node as an input parameter. This
    constructor is expected to be used for reference composite particles
    @param root XML node
    @param autonumbering whether to increment the component indexing
    @param pc particle class */
    SpheroCylinder( DOMNode* root, bool const& autonumbering = true,
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
    @param coordination_number_ particle coordination number
    @param updatePosition whether we update position or not */
    SpheroCylinder( int const& id_, Particle const* ParticleRef,
	double const& vx, double const& vy, double const& vz,
	double const& qrotationx, double const& qrotationy,
	double const& qrotationz, double const& qrotations,
	double const& rx, double const& ry, double const& rz,
	const double m[12],
	ParticleActivity const& activ,
	int const& tag_,
	int const& coordination_number_ = 0,
 	bool const& updatePosition = false );

    /** @brief Constructor with input parameters. This constructor is expected
    to be used for periodic clone particle
    @param id_ ID number
    @param ParticleRef reference particle
    @param vtrans translational velocity
    @param vrot angular velocity
    @param qrot rotation quaternion
    @param config particle transformation
    @param activ particle activity */
    SpheroCylinder( int const& id_, Particle const* ParticleRef,
	Vector3 const& vtrans,
	Quaternion const& qrot,
	Vector3 const& vrot,
	Transform const& config,
	ParticleActivity const& activ );

    /** @brief Copy constructor (the torsor is initialized to 0)
    @param other copied SpheroCylinder object */
    SpheroCylinder( SpheroCylinder const& other );

    /** @brief Destructor */
    virtual ~SpheroCylinder();
    //@}


    /**@name Methods */
    //@{
    /** @brief Creates a clone of the composite particle. This method calls
    the standard copy constructor and is used for new composite particles to be
    inserted in the simulation. Numbering is automatic, total number of
    components is incremented by 1 and activity is set to WAIT. The calling
    object is expected to be a reference composite particle */
    Particle* createCloneCopy() const ;

    /** @brief Creates a clone of the composite particle. This method calls the
    constructor SpheroCylinder( int const& id_, Particle const* ParticleRef,
    Vector3 const& vtrans, Quaternion const& qrot, Vector3 const& vrot,
    Transform const& config, ParticleActivity const& activ ) and is used for
    periodic clone composite particles to be inserted in the simulation.
    Numbering is set with the parameter id_ and total number of components left
    unchanged.
    @param id_ ID number
    @param ParticleRef reference particle
    @param vtrans translational velocity
    @param vrot angular velocity
    @param qrot rotation quaternion
    @param config particle transformation
    @param activ particle activity */
    Particle* createCloneCopy( int const& id_,
    	Particle const* ParticleRef, Vector3 const& vtrans,
	Quaternion const& qrot,	Vector3 const& vrot,
	Transform const& config, ParticleActivity const& activ ) const ;
    //@}


    /**@name Methods I/O */
    //@{
    /** @brief Returns the number of points to write the composite particle in a
    Paraview format */
    int numberOfPoints_PARAVIEW() const;

    /** @brief Returns the number of elementary polytopes to write the
    composite particle shape in a Paraview format */
    int numberOfCells_PARAVIEW() const;

    /** @brief Writes the points describing the composite particle in a
    Paraview format
    @param f output stream
    @param translation additional center of mass translation */
    void write_polygonsPts_PARAVIEW( ostream& f,
  	Vector3 const* translation = NULL ) const;

    /** @brief Returns a list of points describing the component in a
    Paraview format
    @param translation additional center of mass translation */
    list<Point3> get_polygonsPts_PARAVIEW(
  	Vector3 const* translation = NULL ) const;

    /** @brief Writes the composite particle in a Paraview format
    @param connectivity connectivity of Paraview polytopes
    @param offsets connectivity offsets
    @param cellstype Paraview polytopes type
    @param firstpoint_globalnumber global number of the 1st point
    @param last_offset last offset used for the previous convex shape */
    void write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const ;

    /**  @brief Outputs information to be transferred to the fluid
    @param fluid output stream */
    void writePositionInFluid( ostream& fluid );
    //@}


    /**@name Accessors */
    //@{
    /** @brief Returns the number of corners of the rigib body shape and a code
    describing the rigid body shape */
    int getNbCorners() const;   
    //@}
    
    
    /**@name I/O methods */
    //@{
    /** @brief Reads composite particle data from a stream. Usage: for standard
    composite particles in the 2014 reload format
    @param fileIn input stream
    @param referenceParticles reference particles for each class of
    particles */
    void read2014( istream& fileIn, vector<Particle*> const*
  	referenceParticles );

    /** @brief Reads composite particle data from a stream in a binary form.
    Usage: for standard composite particles in the 2014 reload format
    @param fileIn input stream
    @param referenceParticles reference particles for each class of
    particles */
    void read2014_binary( istream& fileIn, vector<Particle*> const*
  	referenceParticles );
    //@}    


  protected:


  private:
    /** @name Parameters */
    //@{
    double m_height; /**< Height of the elementary cylinder */
    double m_radius; /**< Radius of the elementary cylinder and the two
    	spherical caps */
    static int m_visuNodeNbPerQar; /**< number of points over a quarter of 
    	the circular perimeter of the sphero-cylinder for Paraview 
	post-processing */	
    //@}

    /**@name Constructors */
    //@{
    /** @brief Default constructor */
    SpheroCylinder();
    //@}
};

#endif
