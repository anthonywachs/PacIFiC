#ifndef _PARTICLEKINEMATICSSPHERE_HH_
#define _PARTICLEKINEMATICSSPHERE_HH_

#include "ParticleKinematics.hh"
#include "ParticleKinematics3D.hh"
#include "Vector3.hh"
#include "Quaternion.hh"
using namespace solid;


/** @brief The class ParticleKinematicsSphere
    
    Specific methods for the kinematics of 3D spheres.

    @author A.WACHS - Institut Francais du Petrole - 2009 - Creation 
    @author A.WACHS - 2019 - Major cleaning & refactoring */
// ============================================================================
class ParticleKinematicsSphere : public ParticleKinematics3D
{
  public:
    /**@name Constructors */
    //@{
    /** @brief Default constructor */
    ParticleKinematicsSphere();

    /** @brief Copy constructor
    @param copy the ParticleKinematics3D copied object */
    ParticleKinematicsSphere( ParticleKinematicsSphere const& copy );

    /** @brief Destructor */
    ~ParticleKinematicsSphere();
    //@}


    /** @name Methods */
    //@{
    /** @brief Creates and returns a clone of the object */
    ParticleKinematics* clone() const;
    
    /** @brief Computes the angular acceleration in body fixed space
    @param particle particle related to the kinematics 
    @param torque_bf torque in body fixed space
    @param om_bf angular velocity in body-fixed coordinates system 
    @param dOmdt_bf angular acceleration in body fixed space */
    virtual void computeAngularAccelerationBodyFixed( Particle const* particle,
    	Vector3 const& torque_bf, Vector3 const& om_bf, Vector3& dOmdt_bf );  
    //@}

};

#endif
