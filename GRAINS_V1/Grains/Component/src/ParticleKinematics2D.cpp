#include "ParticleKinematics2D.hh"
#include "Particle.hh"
#include "Torsor.hh"
#include <math.h>


// ----------------------------------------------------------------------------
//Default constructor
ParticleKinematics2D::ParticleKinematics2D() 
  : ParticleKinematics()
{}




// ----------------------------------------------------------------------------
// Copy constructor
ParticleKinematics2D::ParticleKinematics2D( ParticleKinematics2D const& copy ) 
  : ParticleKinematics( copy )
{}




// ----------------------------------------------------------------------------
// Destructor
ParticleKinematics2D::~ParticleKinematics2D()
{}




// ----------------------------------------------------------------------------
// Returns the class name
string ParticleKinematics2D::className() const
{
  return ( "*ParticleKinematics2D" );
}




// ----------------------------------------------------------------------------
// Creates and returns a clone of the object
ParticleKinematics* ParticleKinematics2D::clone() const
{
  return ( new ParticleKinematics2D( *this ) );
}




// ----------------------------------------------------------------------------
// Computes the angular acceleration in body fixed space
void ParticleKinematics2D::computeAngularAccelerationBodyFixed( 
	Particle const* particle, Vector3 const& torque_bf,
	Vector3 const& om_bf, Vector3& dOmdt_bf )
{
  const double *inertia = particle->getInertiaTensorBodyFixed();
  dOmdt_bf[Z] = torque_bf[Z] / ( m_coupling_factor * inertia[5] );
}




// ----------------------------------------------------------------------------
// Output operator
ostream& operator << ( ostream& fileOut, ParticleKinematics2D const& kine_ )
{
  fileOut << "*ParticleKinematics2D\n";
  fileOut << kine_.m_translationalVelocity
	<< kine_.m_QuaternionRotation
	<< kine_.m_angularVelocity;
	
  return ( fileOut );
}




// ----------------------------------------------------------------------------
// Input operator
istream& operator >> ( istream& fileIn, ParticleKinematics2D& kine_ )
{
  fileIn >> kine_.m_translationalVelocity
	>> kine_.m_QuaternionRotation
	>> kine_.m_angularVelocity;
  kine_.m_angularVelocity[X] = kine_.m_angularVelocity[Y] = 0.;	
		
  return ( fileIn );
}




// ----------------------------------------------------------------------------
// Computes explicitly Idw/dt
Vector3 ParticleKinematics2D::computeExplicitDJomDt( Vector3 const & dw,
	double const* inertie ) const
{
  Vector3 Idw;

  // Calcul de I.dw
  Idw[Z] = inertie[5] * dw[Z];  
  
  return ( Idw );
}
