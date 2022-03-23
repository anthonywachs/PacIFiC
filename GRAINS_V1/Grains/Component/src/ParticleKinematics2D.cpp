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
// Computes the momentum change over dt
void ParticleKinematics2D::computeMomentumChangeOverDt( Torsor const& torseur,
	double dt, Particle const* particle )
{
  double couplingFactor = 1.;

  // Translational momentum
  m_dUdt = *torseur.getForce() / ( particle->getMass() * couplingFactor );
  m_dUdt.round();
  m_dUdt[Z] = 0.;

  // Angular momentum
  double const* inertieInverse = particle->getInverseInertiaTensor();
  m_dOmegadt[Z] = (*torseur.getTorque())[Z] * inertieInverse[5] 
  	/ couplingFactor;
}




// ----------------------------------------------------------------------------
// Output operator
ostream& operator << ( ostream& fileOut, ParticleKinematics2D const& kine_ )
{
  fileOut << "*ParticleKinematics2D\n";
  fileOut << kine_.m_translationalVelocity
	<< kine_.m_QuaternionRotation
	<< kine_.m_dQuaternionRotationdt;
	
  return ( fileOut );
}




// ----------------------------------------------------------------------------
// Input operator
istream& operator >> ( istream& fileIn, ParticleKinematics2D& kine_ )
{
  fileIn >> kine_.m_translationalVelocity
	>> kine_.m_QuaternionRotation
	>> kine_.m_dQuaternionRotationdt;
  kine_.m_angularVelocity = 
  	2.0 * kine_.m_dQuaternionRotationdt.multConjugateToVector3( 
  	kine_.m_QuaternionRotation );
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
