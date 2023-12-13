#include "FirstOrderExplicit.hh"


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
void FirstOrderExplicit::Move( Vector3 const& dUdt, Vector3& vtrans, 
	Vector3& transDisplacement, Vector3 const& dOmegadt,
	Vector3& vrot, Vector3& meanVRot, double const& dt_particle_vel, 
    	double const& dt_particle_disp )
{
  // Translational velocity and displacement
  transDisplacement = vtrans * dt_particle_disp;
  vtrans += dUdt * dt_particle_vel;

  // Angular velocity and displacement
  meanVRot = vrot;  
  vrot += dOmegadt * dt_particle_vel;  
}
