#include "SecondOrderLeapFrog.hh"


// ----------------------------------------------------------------------------
// Default constructeur
SecondOrderLeapFrog::SecondOrderLeapFrog() 
  : TimeIntegrator()
{}




// ----------------------------------------------------------------------------
// Copy constructor
SecondOrderLeapFrog::SecondOrderLeapFrog( SecondOrderLeapFrog const& copy ) 
  : TimeIntegrator( copy )
{}




// ----------------------------------------------------------------------------
// Destructor
SecondOrderLeapFrog::~SecondOrderLeapFrog()
{}




// ----------------------------------------------------------------------------
// Creates and returns a clone of the time integrator
TimeIntegrator* SecondOrderLeapFrog::clone() const
{
  return ( new SecondOrderLeapFrog(*this) );
}




// ----------------------------------------------------------------------------
// Computes the new velocity and position at time t+dt
void SecondOrderLeapFrog::Move( Vector3 const& dUdt, Vector3& vtrans, 
	Vector3& transDisplacement, Vector3 const& dOmegadt,
	Vector3& vrot, Vector3& meanVRot, double const& dt_particle_vel, 
    	double const& dt_particle_disp )
{
  // Translational velocity and displacement
  vtrans += dUdt * dt_particle_vel;
  transDisplacement = vtrans * dt_particle_disp;

  // Angular velocity and displacement
  vrot += dOmegadt * dt_particle_vel;
  meanVRot = vrot;    
}




// ----------------------------------------------------------------------------
// Advances velocity over dt_particle_vel
void SecondOrderLeapFrog::advanceVelocity( Vector3 const& dUdt, Vector3& vtrans,
	Vector3 const& dOmegadt, Vector3& vrot, double const& dt_particle_vel )
{
  // Translational velocity and displacement
  vtrans += dUdt * dt_particle_vel;

  // Angular velocity and displacement
  vrot += dOmegadt * dt_particle_vel;
}

