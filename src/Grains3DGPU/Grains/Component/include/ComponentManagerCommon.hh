#ifndef _COMPONENTMANAGERCOMMON_HH_
#define _COMPONENTMANAGERCOMMON_HH_

#include "Basic.hh"
#include "GrainsParameters.hh"
#include "Kinematics.hh"
#include "Quaternion.hh"
#include "QuaternionMath.hh"
#include "RigidBody.hh"
#include "TimeIntegrator.hh"
#include "Torce.hh"
#include "Vector3.hh"
#include "VectorMath.hh"

// =================================================================================================
/** @brief ComponentManager common functions between host and device.

    This is a header-only file that contains common functions between CPU/GPU
    for the ComponentManager class. The functions are templated to allow
    for flexibility in usage. The functions are marked as inline to allow for
    better optimization by the compiler.

    Force computation helpers (computeContactForces_common, reduceTorces_common,
    addExternalForces_common) have been moved to ForceModuleCommon.hh.

    @author A.Yazdani - 2025 - Construction */
// =================================================================================================
/** @name ComponentManager common functions between host and device */
//@{

// -------------------------------------------------------------------------------------------------
/** @brief Moves a component using the given time integration method
    @param TI the time integrator
    @param rigidBody the rigid body of the components
    @param position the position of the component
    @param quaternion the quaternion of the component
    @param kinematics the kinematics of the component
    @param torce the torce acting on the component
    @param cID the ID of the component */
template <typename T>
__HOSTDEVICE__ static INLINE void moveParticles_common(const TimeIntegrator<T>* const* TI,
                                                       const RigidBody<T>* const*      rigidBody,
                                                       Vector3<T>*                     position,
                                                       Quaternion<T>*                  quaternion,
                                                       Kinematics<T>*                  kinematics,
                                                       Torce<T>*                       torce,
                                                       const uint                      cID)
{
    // Rigid body
    const RigidBody<T>* rb = rigidBody[cID];
    // Computing momentums in the space-fixed coordinate
    const Kinematics<T>& momentum
        = rb->computeMomentum(kinematics[cID].getAngularComponent(), torce[cID], quaternion[cID]);
    // Reset torces
    torce[cID].reset();
    // Finally, we move particles using the given time integration
    Vector3<T>    transMotion;
    Quaternion<T> rotMotion;
    TI[0]->Move(momentum, kinematics[cID], transMotion, rotMotion);

    position[cID] += transMotion;
    quaternion[cID] = rotMotion * quaternion[cID];
    T qn            = norm(quaternion[cID]);
    if(qn > EPS<T>)
        quaternion[cID] *= (T(1) / qn);
}

// -------------------------------------------------------------------------------------------------
/** @brief Performs the second velocity half-kick for split-step schemes (KDK leapfrog Step 3).
    Computes the acceleration from the current torce (forces at x_{n+1}) and delegates to
    TI[0]->AdvanceVelocity. The torce is intentionally NOT reset here so that it remains
    available as a_n for the next call to moveParticles_common.
    @param TI the time integrator
    @param rigidBody the rigid body of the components
    @param quaternion the quaternion of the component
    @param kinematics the kinematics of the component (velocity updated in-place)
    @param torce the accumulated torce at x_{n+1} (read but not reset)
    @param cID the ID of the component */
template <typename T>
__HOSTDEVICE__ static INLINE void advanceVelocity_common(const TimeIntegrator<T>* const* TI,
                                                         const RigidBody<T>* const*      rigidBody,
                                                         const Quaternion<T>*            quaternion,
                                                         Kinematics<T>*                  kinematics,
                                                         const Torce<T>*                 torce,
                                                         const uint                      cID)
{
    const RigidBody<T>* rb = rigidBody[cID];
    // Compute acceleration from forces at x_{n+1} -- torce is NOT reset
    const Kinematics<T> acceleration
        = rb->computeMomentum(kinematics[cID].getAngularComponent(), torce[cID], quaternion[cID]);
    TI[0]->AdvanceVelocity(acceleration, kinematics[cID]);
}
//@}

#endif
