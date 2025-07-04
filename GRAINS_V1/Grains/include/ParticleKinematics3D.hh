#ifndef _PARTICLEKINEMATICS3D_HH_
#define _PARTICLEKINEMATICS3D_HH_

#include "ParticleKinematics.hh"

class ParticleKinematics3D;
ostream& operator << ( ostream& fileOut, ParticleKinematics3D const& kine_ );
istream& operator >> ( istream& fileIn, ParticleKinematics3D& kine_ );


/** @brief The class ParticleKinematics3D.
    
    Specific methods for the kinematics of 3D particles.

    @author Institut Francais du Petrole - 2000 - Creation 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class ParticleKinematics3D : public ParticleKinematics
{
  public:
    /**@name Constructors */
    //@{
    /** @brief Default constructor */
    ParticleKinematics3D();

    /** @brief Copy constructor
    @param copy the ParticleKinematics3D copied object */
    ParticleKinematics3D( ParticleKinematics3D const& copy );

    /** @brief Destructor */
    ~ParticleKinematics3D();
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
    @param kine_ ParticleKinematics3D object */
    friend ostream& operator << ( ostream& fileOut, 
	ParticleKinematics3D const& kine_ );
	
    /** @brief Input operator
    @param fileIn input stream
    @param kine_ ParticleKinematics3D object */
    friend istream& operator >> ( istream& fileIn, 
    	ParticleKinematics3D& kine_ );
    //@}
};

#endif
