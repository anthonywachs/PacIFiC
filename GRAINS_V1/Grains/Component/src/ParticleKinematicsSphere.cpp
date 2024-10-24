#include "Grains.hh"
#include "ParticleKinematicsSphere.hh"
#include "Particle.hh"
#include "Torsor.hh"
#include <math.h>


// ----------------------------------------------------------------------------
// Default constructor
ParticleKinematicsSphere::ParticleKinematicsSphere() 
  : ParticleKinematics3D()
{}




// ----------------------------------------------------------------------------
// Copy constructor
ParticleKinematicsSphere::ParticleKinematicsSphere( 
	ParticleKinematicsSphere const& copy ) 
  : ParticleKinematics3D( copy )
{}




// ----------------------------------------------------------------------------
// Destructor
ParticleKinematicsSphere::~ParticleKinematicsSphere()
{}




// ----------------------------------------------------------------------------
// Creates and returns a clone of the object
ParticleKinematics* ParticleKinematicsSphere::clone() const
{
  return ( new ParticleKinematicsSphere( *this ) );
}




// ----------------------------------------------------------------------------
// Computes the angular acceleration in body fixed space
void ParticleKinematicsSphere::computeAngularAccelerationBodyFixed( 
	Particle const* particle, Vector3 const& torque_bf,
	Vector3 const& om_bf, Vector3& dOmdt_bf )
{
  const double *inertia = particle->getInertiaTensorBodyFixed();
  dOmdt_bf[X] = torque_bf[X] / ( m_coupling_factor * inertia[0] );
  dOmdt_bf[Y] = torque_bf[Y] / ( m_coupling_factor * inertia[3] );
  dOmdt_bf[Z] = torque_bf[Z] / ( m_coupling_factor * inertia[5] );
}
