#include "TimeIntegrator.hh"


// ----------------------------------------------------------------------------
// Default constructor
TimeIntegrator::TimeIntegrator()
{}




// ----------------------------------------------------------------------------
// Copy constructor
TimeIntegrator::TimeIntegrator( TimeIntegrator const& copy )
{}




// ----------------------------------------------------------------------------
// Destructor
TimeIntegrator::~TimeIntegrator()
{}




// ----------------------------------------------------------------------------
// Advances velocity over dt_particle_vel
void TimeIntegrator::advanceVelocity( Vector3 const& dUdt, Vector3& vtrans, 
	Vector3 const& dOmegadt, Vector3& vrot, double const& dt_particle_vel )
{}
