#include "TimeIntegrator.hh"
#include "ParticleKinematics.hh"
#include "Particle.hh"


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
void TimeIntegrator::advanceVelocity( Particle* particle, 
	ParticleKinematics* kine, double const& coupling_factor, 
	Vector3 const& torque_bf, Vector3& vtrans, Vector3& vrot, 
	double const& dt_particle_vel )
{}




// ----------------------------------------------------------------------------
// Copies kinematics at time t-2dt (translational velocity, angular 
// velocity, variation of translational velocity, variation of angular 
// velocity) in a 1D array 
void TimeIntegrator::copyKinematicsNm2( double* vit, int i ) const 
{}





// ----------------------------------------------------------------------------
// Sets kinematics at time t-2dt from a 1D array of 12 scalars
// (translational velocity, angular velocity, variation of translational
// velocity, variation of angular velocity)
void TimeIntegrator::setKinematicsNm2( double const* tab )
{}



// ----------------------------------------------------------------------------
// Writes time integrator data in an output stream with a high
// precision and 2014 format
void TimeIntegrator::writeParticleKinematics2014( ostream& fileOut, 
    	Particle const* particle ) const
{} 



  
// ----------------------------------------------------------------------------
// Writes time integrator data in an output stream with a binary and 2014 format
void TimeIntegrator::writeParticleKinematics2014_binary( ostream& fileOut, 
    	Particle const* particle )
{}




// ----------------------------------------------------------------------------
// Reads time integrator data from a stream in the 2014 format 
void TimeIntegrator::readParticleKinematics2014( istream& StreamIN,
    	Particle* particle )
{} 



  
// ----------------------------------------------------------------------------
// Reads time integrator data from a stream in a binary form in the 2014 format 
void TimeIntegrator::readParticleKinematics2014_binary( istream& StreamIN,
    	Particle* particle )
{}




// ----------------------------------------------------------------------------
// Returns the number of bytes of the time integrator data when written in a 
// binary format to an output stream
size_t TimeIntegrator::get_numberOfBytes() const
{
  return ( 0 );
}
