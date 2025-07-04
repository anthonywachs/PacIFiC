#ifndef _PARTICLEKINEMATICS2D_HH_
#define _PARTICLEKINEMATICS2D_HH_

#include "ParticleKinematics.hh"

class ParticleKinematics2D;
ostream& operator << ( ostream& fileOut, ParticleKinematics2D const& kine_ );
istream& operator >> ( istream& fileIn, ParticleKinematics2D& kine_ );


/** @brief The class ParticleKinematics2D.
    
    Specific methods for the kinematics of 2D particles.

    @author Institut Francais du Petrole - 2000 - Creation 
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

    /** @brief Computes the angular acceleration in body fixed space
    @param particle particle related to the kinematics 
    @param torque_bf torque in body fixed space
    @param om_bf angular velocity in body-fixed coordinates system 
    @param dOmdt_bf angular acceleration in body fixed space */
    virtual void computeAngularAccelerationBodyFixed( Particle const* particle,
    	Vector3 const& torque_bf, Vector3 const& om_bf, Vector3& dOmdt_bf );
			   
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
