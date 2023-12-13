#ifndef _PARTICLEKINEMATICS2D_HH_
#define _PARTICLEKINEMATICS2D_HH_

#include "ParticleKinematics.hh"

class ParticleKinematics2D;
ostream& operator << ( ostream& fileOut, ParticleKinematics2D const& kine_ );
istream& operator >> ( istream& fileIn, ParticleKinematics2D& kine_ );


/** @brief The class ParticleKinematics2D.
    
    Specific methods for the kinematics of 2D particles.

    @author F.PRADEL - Institut Francais du Petrole - 2000 - Creation 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class ParticleKinematics2D : public ParticleKinematics
{
  public:
    /**@name Constructors */
    //@{
    /** @brief Default constructor */
    ParticleKinematics2D();

    /** @brief Copy constructor
    @param copy the ParticleKinematics2D copied object */
    ParticleKinematics2D( ParticleKinematics2D const& copy );

    /** @brief Destructor */
    ~ParticleKinematics2D();
    //@}


    /**@name Methods */
    //@{
    /** @brief Returns the class name */
    string className() const;

    /** @brief Creates and returns a clone of the object */
    ParticleKinematics* clone() const;

    /** @brief Computes the momentum change over dt 
    @param torseur particle torsor 
    @param particle particle related to the kinematics */
    void computeAcceleration( Torsor const& torseur,
	Particle const* particle ) ;
			   
    /** @brief Computes explicitly Idw/dt
    @param dw explicit change of angular velocity
    @param inertie particle inertia tensor */
    Vector3 computeExplicitDJomDt( Vector3 const& dw,
	double const* inertie ) const ;
    //@}


    /**@name Friend methods */
    //@{
    /** @brief Output operator
    @param fileOut output stream
    @param kine_ ParticleKinematics2D object */
    friend ostream& operator << ( ostream& fileOut, 
	ParticleKinematics2D const& kine_ );
	
    /** @brief Input operator
    @param fileIn input stream
    @param kine_ ParticleKinematics2D object */
    friend istream& operator >> ( istream& fileIn, 
    	ParticleKinematics2D& kine_ );
    //@}
};

#endif
