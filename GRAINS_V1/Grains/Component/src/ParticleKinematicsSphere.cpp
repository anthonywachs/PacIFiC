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
void ParticleKinematicsSphere::computeAcceleration( 
	Torsor const& torseur, Particle const* particle )
{
  // Values of the coupling factor:
  // 1) purely granular
  //    FluidCorrectedAcceleration is true and fluid density is 0
  //    so coupling factor = 1
  // 2) coupled to fluid, fluid density is not 0
  //    a. FluidCorrectedAcceleration is true so coupling factor = 1 - 
  //    fluid density / particle density
  //    b. FluidCorrectedAcceleration is false so coupling factor = 1  
  double couplingFactor = 1.;
  if ( Particle::getFluidCorrectedAcceleration() )
    couplingFactor -=
    	Particle::getFluidDensity() / particle->getDensity(); 

  // Translational momentum
  m_dUdt = *(torseur.getForce()) / ( particle->getMass() * couplingFactor );
  m_dUdt.round();

  // Angular momentum
  Vector3 const* Moment = torseur.getTorque();

  double const* inverseInertia = particle->getInverseInertiaTensorBodyFixed();
  m_dOmegadt[0] = inverseInertia[0] * (*Moment)[0] / couplingFactor;
  m_dOmegadt[1] = inverseInertia[3] * (*Moment)[1] / couplingFactor;
  m_dOmegadt[2] = inverseInertia[5] * (*Moment)[2] / couplingFactor;
}
