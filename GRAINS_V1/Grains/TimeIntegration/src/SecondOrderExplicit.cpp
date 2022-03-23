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
void SecondOrderExplicit::Move( Vector3& vtrans, Vector3 const& dUdt,
	Vector3& transDisplacement, Vector3 const& dOmegadt,
	Vector3& vrot, Vector3& meanVRot, double dt )
{
  // Velocity et deplacement translationnels
  transDisplacement = vtrans * dt + 0.5 * dUdt * dt * dt;
  vtrans += dUdt * dt;

  // Velocity et deplacement rotationnels
  meanVRot = vrot + 0.5 * dOmegadt * dt;  
  vrot += dOmegadt * dt;  
}
