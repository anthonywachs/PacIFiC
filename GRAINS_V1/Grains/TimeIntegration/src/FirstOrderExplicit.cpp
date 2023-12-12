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
void FirstOrderExplicit::Move( Vector3& vtrans, Vector3 const& dUdt,
	Vector3& transDisplacement, Vector3 const& dOmegadt,
	Vector3& vrot, Vector3& meanVRot, double dt )
{
  // Translational velocity and displacement
  transDisplacement = vtrans * dt;
  vtrans += dUdt * dt;

  // Angular velocity and displacement
  meanVRot = vrot;  
  vrot += dOmegadt * dt;  
}
