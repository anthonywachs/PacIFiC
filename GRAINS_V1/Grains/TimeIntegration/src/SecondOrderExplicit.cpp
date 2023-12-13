#include "SecondOrderExplicit.hh"


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
void SecondOrderExplicit::Move( Vector3 const& dUdt, Vector3& vtrans, 
	Vector3& transDisplacement, Vector3 const& dOmegadt,
	Vector3& vrot, Vector3& meanVRot, double const& dt_particle_vel, 
    	double const& dt_particle_disp )
{
  // Translational velocity and displacement
  transDisplacement = vtrans * dt_particle_disp 
  	+ 0.5 * dUdt * dt_particle_disp * dt_particle_disp;
  vtrans += dUdt * dt_particle_vel;

  // Angular velocity and displacement
  meanVRot = vrot + 0.5 * dOmegadt * dt_particle_disp;  
  vrot += dOmegadt * dt_particle_vel;  
}
