#include "FirstOrderExplicit.hh"
#include "ParticleKinematics.hh"
#include "Particle.hh"


// ----------------------------------------------------------------------------
// Default constructor
FirstOrderExplicit::FirstOrderExplicit() 
  : TimeIntegrator()
{}




// ----------------------------------------------------------------------------
// Copy constructor
FirstOrderExplicit::FirstOrderExplicit( FirstOrderExplicit const& copy ) 
  : TimeIntegrator( copy )
{}




// ----------------------------------------------------------------------------
// Destructor
FirstOrderExplicit::~FirstOrderExplicit()
{}




// ----------------------------------------------------------------------------
// Creates and returns a clone of the time integrator
TimeIntegrator* FirstOrderExplicit::clone() const
{
  return ( new FirstOrderExplicit(*this) );
}




// ----------------------------------------------------------------------------
// Computes the new velocity and position at time t+dt
void FirstOrderExplicit::Move( Particle* particle, ParticleKinematics* kine,
	double const& coupling_factor, Vector3 const& torque_bf,
	Vector3& vtrans, Vector3& transMotion, 
	Vector3& vrot, Vector3& meanVRot, double const& dt_particle_vel, 
    	double const& dt_particle_disp )
{
  // Translational velocity and motion
  transMotion = vtrans * dt_particle_disp;   
  vtrans += ( *(particle->getForce()) / 
  	( particle->getMass() * coupling_factor ) ) * dt_particle_vel;

  // Angular velocity
  Vector3 dOmegadt_bf_nm1;
  kine->computeAngularAccelerationBodyFixed( particle, torque_bf, vrot, 
  	dOmegadt_bf_nm1 );
  meanVRot = vrot;
  vrot += dOmegadt_bf_nm1 * dt_particle_vel;     
}
