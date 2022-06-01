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
void SecondOrderLeapFrog::Move( Vector3& vtrans, Vector3 const& dUdt,
	Vector3& transDisplacement, Vector3 const& dOmegadt,
	Vector3& vrot, Vector3& meanVRot, double dt )
{
  // Velocity et deplacement translationnels
  vtrans += dUdt * dt;
  transDisplacement = vtrans * dt;

  // Velocity et deplacement rotationnels
  vrot += dOmegadt * dt;
  meanVRot = vrot;    
}
