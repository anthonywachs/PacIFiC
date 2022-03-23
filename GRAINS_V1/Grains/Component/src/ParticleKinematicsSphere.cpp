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
// Computes the momentum change over dt 
void ParticleKinematicsSphere::computeMomentumChangeOverDt( 
	Torsor const& torseur,
	double dt, Particle const* particle )
{
  double couplingFactor = 1.;

  // Translational momentum
  m_dUdt = *(torseur.getForce()) / ( particle->getMass() * couplingFactor );
  m_dUdt.round();

  // Angular momentum
  Vector3 const* Moment = torseur.getTorque();

  double const* inertieInverse = particle->getInverseInertiaTensor();
  m_dOmegadt[0] = inertieInverse[0] * (*Moment)[0] / couplingFactor;
  m_dOmegadt[1] = inertieInverse[3] * (*Moment)[1] / couplingFactor;
  m_dOmegadt[2] = inertieInverse[5] * (*Moment)[2] / couplingFactor;
}
