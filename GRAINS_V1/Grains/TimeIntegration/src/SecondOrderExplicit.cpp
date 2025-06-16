#include "SecondOrderExplicit.hh"
#include "ParticleKinematics.hh"
#include "Particle.hh"


// ----------------------------------------------------------------------------
// Default constructeur
SecondOrderExplicit::SecondOrderExplicit() 
  : TimeIntegrator()
{}




// ----------------------------------------------------------------------------
// Copy constructor
SecondOrderExplicit::SecondOrderExplicit( SecondOrderExplicit const& copy ) 
  : TimeIntegrator( copy )
{}




// ----------------------------------------------------------------------------
// Destructor
SecondOrderExplicit::~SecondOrderExplicit()
{}




// ----------------------------------------------------------------------------
// Creates and returns a clone of the time integrator
TimeIntegrator* SecondOrderExplicit::clone() const
{
  return ( new SecondOrderExplicit(*this) );
}




// ----------------------------------------------------------------------------
// Computes the new velocity and position at time t+dt
void SecondOrderExplicit::Move( Particle* particle, ParticleKinematics* kine,
	double const& coupling_factor, Vector3 const& torque_bf,
	Vector3& vtrans, Vector3& transMotion, 
	Vector3& vrot, Vector3& meanVRot, double const& dt_particle_vel, 
    	double const& dt_particle_disp )
{
  // Translational velocity and motion
  Vector3 dUdt = *(particle->getForce()) / 
  	( particle->getMass() * coupling_factor );
  transMotion = vtrans * dt_particle_disp 
  	+ 0.5 * dUdt * dt_particle_disp * dt_particle_disp;
  vtrans += dUdt * dt_particle_vel;

  // Angular velocity
  Vector3 dOmegadt_bf;
  kine->computeAngularAccelerationBodyFixed( particle, torque_bf, vrot, 
  	dOmegadt_bf );
  meanVRot = vrot + 0.5 * dOmegadt_bf * dt_particle_disp;
  vrot += dOmegadt_bf * dt_particle_vel;   
}
