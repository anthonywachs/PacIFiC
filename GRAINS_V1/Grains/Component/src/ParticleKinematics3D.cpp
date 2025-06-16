#include "Grains.hh"
#include "ParticleKinematics3D.hh"
#include "Particle.hh"
#include "Torsor.hh"
#include <math.h>


// ----------------------------------------------------------------------------
// Default constructor
ParticleKinematics3D::ParticleKinematics3D() 
  : ParticleKinematics()
{}




// ----------------------------------------------------------------------------
// Copy constructor
ParticleKinematics3D::ParticleKinematics3D( ParticleKinematics3D const& copy ) 
  : ParticleKinematics( copy )
{}




// ----------------------------------------------------------------------------
// Destructeur
ParticleKinematics3D::~ParticleKinematics3D()
{}




// ----------------------------------------------------------------------------
// Returns the class name
string ParticleKinematics3D::className() const
{
  return ( "*ParticleKinematics3D" );
}




// ----------------------------------------------------------------------------
// Creates and returns a clone of the object
ParticleKinematics* ParticleKinematics3D::clone() const
{
  return ( new ParticleKinematics3D( *this ) );
}




// ----------------------------------------------------------------------------
// Computes the angular acceleration in body fixed space
void ParticleKinematics3D::computeAngularAccelerationBodyFixed( 
	Particle const* particle, Vector3 const& torque_bf,
	Vector3 const& om_bf, Vector3& dOmdt_bf )
{
  // Note: in the body-fixed space, the moment of inertia tensor is diagonal
  double const* inertia = particle->getInertiaTensorBodyFixed();    
  if ( Grains::isModePredictor() )
  {   
    dOmdt_bf[X] = torque_bf[X] / ( m_coupling_factor * inertia[0] )
    	+ om_bf[Y] * om_bf[Z] * ( inertia[3] - inertia[5] ) / inertia[0];
    dOmdt_bf[Y] = torque_bf[Y] / ( m_coupling_factor * inertia[3] )
    	+ om_bf[Z] * om_bf[X] * ( inertia[5] - inertia[0] ) / inertia[3];
    dOmdt_bf[Z] = torque_bf[Z] / ( m_coupling_factor * inertia[5] )
    	+ om_bf[X] * om_bf[Y] * ( inertia[0] - inertia[3] ) / inertia[5];	
  }
  else
  {
    dOmdt_bf[X] = torque_bf[X] / ( m_coupling_factor * inertia[0] );
    dOmdt_bf[Y] = torque_bf[Y] / ( m_coupling_factor * inertia[3] );
    dOmdt_bf[Z] = torque_bf[Z] / ( m_coupling_factor * inertia[5] );
  }  
}  




// ----------------------------------------------------------------------------
// Output operator
ostream& operator << ( ostream& fileOut, ParticleKinematics3D const& kine_ )
{
  fileOut << "*ParticleKinematics3D\n";
  fileOut << kine_.m_translationalVelocity
	<< kine_.m_QuaternionRotation
	<< kine_.m_angularVelocity;

  return ( fileOut );
}




// ----------------------------------------------------------------------------
// Input operator
istream& operator >> ( istream& fileIn, ParticleKinematics3D& kine_ )
{
  fileIn >> kine_.m_translationalVelocity
	>> kine_.m_QuaternionRotation
	>> kine_.m_angularVelocity;
	
  return ( fileIn );
}




// ----------------------------------------------------------------------------
// Computes explicitly Idw/dt
Vector3 ParticleKinematics3D::computeExplicitDJomDt( Vector3 const& dw,
	double const* inertie ) const
{
  Vector3 Idw;
  Quaternion pConjugue = m_QuaternionRotation.Conjugate();
  
  // Write dw in body-fixed coordinates system
  Vector3 vTmp = pConjugue.multToVector3( ( dw , m_QuaternionRotation ) );
  
  // Compute I.dw in body-fixed coordinates system
  Vector3 work; 
  work[0] = inertie[0]*vTmp[0] + inertie[1]*vTmp[1] + inertie[2]*vTmp[2];
  work[1] = inertie[1]*vTmp[0] + inertie[3]*vTmp[1] + inertie[4]*vTmp[2];
  work[2] = inertie[2]*vTmp[0] + inertie[4]*vTmp[1] + inertie[5]*vTmp[2];
  
  // Write I.dw back in the space-fixed coordinates system
  Idw = m_QuaternionRotation.multToVector3( ( work , pConjugue ) );  
     
  return ( Idw );
}
