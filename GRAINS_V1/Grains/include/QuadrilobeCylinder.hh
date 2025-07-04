#ifndef _QUADRILOBECYLINDER_HH_
#define _QUADRILOBECYLINDER_HH_

#include "CompositeParticle.hh"


/** @brief The class QuadrilobeCylinder.

    A freely moving quadrilobe cylinder.

    @author A.WACHS - 2023 - Creation */
// ============================================================================
class QuadrilobeCylinder : public CompositeParticle
{
  public:
    /**@name Constructors & Destructor */
    //@{
    /** @brief Constructor with autonumbering as input parameter
    @param autonumbering whether to increment the component indexing */
    QuadrilobeCylinder( bool const& autonumbering );

    /** @brief Constructor with an XML node as an input parameter. This
    constructor is expected to be used for reference composite particles.
    Autonumbering is set to false
    @param root XML node
    @param pc particle class */
    QuadrilobeCylinder( DOMNode* root, int const& pc );

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
    QuadrilobeCylinder( int const& id_, Particle const* ParticleRef,
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
    QuadrilobeCylinder( int const& id_, Particle const* ParticleRef,
	Vector3 const& vtrans,
	Quaternion const& qrot,
	Vector3 const& vrot,
	Transform const& config,
	ParticleActivity const& activ,
     	map< std::tuple<int,int,int>,
     	std::tuple<bool, Vector3, Vector3, Vector3> > const* contactMap );

    /** @brief Copy constructor (the torsor is initialized to 0)
    @param other copied QuadrilobeCylinder object 
    @param autonumbering whether to increment the component indexing */
    QuadrilobeCylinder( QuadrilobeCylinder const& other, 
    	bool const& autonumbering );

    /** @brief Destructor */
    virtual ~QuadrilobeCylinder();
    //@}


    /**@name Methods */
    //@{
    /** @brief Creates a clone of the composite particle. This method calls
    the standard copy constructor and is used for new composite particles to be
    inserted in the simulation. Activity is set to WAIT. The calling
    object is expected to be a reference composite particle */
    Particle* createCloneCopy( bool const& autonumbering ) const ;

    /** @brief Creates a clone of the composite particle. This method calls the
    constructor QuadrilobeCylinder( int const& id_, Particle const* ParticleRef,
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
     	std::tuple<bool, Vector3, Vector3, Vector3> > const* contactMap ) const ;
    //@}


    /**@name Methods I/O */
    //@{
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
    /**@name I/O methods */
    //@{
    /** @brief Saves additional features of a (in practice reference) composite
    particle for reload
    @param fileSave output stream */
    void writeAdditionalFeatures( ostream& fileSave ) const;

    /** @brief Reads additional features of a (in practice reference) particle
    data from a stream
    @param fileIn input stream */
    void readAdditionalFeatures( istream& fileIn );
    //@}
    
    
    /** @name Methods */
    //@{
    /** @brief Computes and sets the circumscribed radius */
    void setCircumscribedRadius();
    //@}    
    

  private:
    /** @name Parameters */
    //@{
    double m_radius; /**< Radius of the two half cylinders */
    double m_armLength; /**< Length of each of the 4 arms of the quadrilobe,
    	from the center of mass to end of each arm minus the radius of the 
	cylindrical cap */ 
    double m_height; /**< Height */       
    static int m_visuNodeNbPerHalf; /**< number of points over 
    	the circular perimeter of each half cylinder (half of a cricle ) 
	for Paraview post-processing */	
    //@}

    /**@name Constructors */
    //@{
    /** @brief Default constructor */
    QuadrilobeCylinder();
    //@}
};

#endif
