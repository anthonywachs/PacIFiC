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
// Computes the momentum change over dt 
void ParticleKinematics3D::computeMomentumChangeOverDt( const Torsor &torseur,
	double dt, Particle const* particle )
{
  double couplingFactor = 1.;

  // Translational momentum
  m_dUdt = *(torseur.getForce()) / ( particle->getMass() * couplingFactor );
  m_dUdt.round();

  // Angular momentum
  Vector3 work, vTmp;  

  // Rotation quaternion conjugate
  Quaternion pConjugue = m_QuaternionRotation.Conjugate(); 
	
  // Inverse inertia tensor
  double const* inertieInverse = particle->getInverseInertiaTensor(); 
  
  if ( Grains::isModePredictor() )
  {      
    // Write omega in body-fixed coordinates system
    Vector3 omb = pConjugue.multToVector3( 
    	( m_angularVelocity , m_QuaternionRotation ) );
	
    // Compute I.w in body-fixed coordinates system
    const double *inertie = particle->getInertiaTensor(); 
    vTmp[0] = inertie[0]*omb[0] + inertie[1]*omb[1] + inertie[2]*omb[2];
    vTmp[1] = inertie[1]*omb[0] + inertie[3]*omb[1] + inertie[4]*omb[2];
    vTmp[2] = inertie[2]*omb[0] + inertie[4]*omb[1] + inertie[5]*omb[2];
    
    // Compute I.w ^ w in body-fixed coordinates system 
    work = vTmp ^ omb;
    
    // Write torque in body-fixed coordinates system
    vTmp = pConjugue.multToVector3( 
    	( *(torseur.getTorque()) , m_QuaternionRotation ) );
	
    // Compute T + I.w ^ w in body-fixed coordinates system     
    work += vTmp; 
  }
  else
    // Write torque in body-fixed coordinates system
    work = pConjugue.multToVector3( 
    	( *(torseur.getTorque()) , m_QuaternionRotation ) );
	 
  // Compute I^-1.(T + I.w ^ w) in body-fixed coordinates system 
  vTmp[0] = inertieInverse[0]*work[0] + inertieInverse[1]*work[1] 
    + inertieInverse[2]*work[2];
  vTmp[1] = inertieInverse[1]*work[0] + inertieInverse[3]*work[1] 
    + inertieInverse[4]*work[2];
  vTmp[2] = inertieInverse[2]*work[0] + inertieInverse[4]*work[1] 
    + inertieInverse[5]*work[2];
  
  // Write I^-1.(T + I.w ^ w) in space-fixed coordinates system
  work = m_QuaternionRotation.multToVector3( ( vTmp , pConjugue ) );

  // Compute m_dOmegadt
  m_dOmegadt = work / couplingFactor; 
}  




// ----------------------------------------------------------------------------
// Output operator
ostream& operator << ( ostream& fileOut, ParticleKinematics3D const& kine_ )
{
  fileOut << "*ParticleKinematics3D\n";
  fileOut << kine_.m_translationalVelocity
	<< kine_.m_QuaternionRotation
	<< kine_.m_dQuaternionRotationdt;

  return ( fileOut );
}




// ----------------------------------------------------------------------------
// Input operator
istream& operator >> ( istream& fileIn, ParticleKinematics3D& kine_ )
{
  fileIn >> kine_.m_translationalVelocity
	>> kine_.m_QuaternionRotation
	>> kine_.m_dQuaternionRotationdt;
  kine_.m_angularVelocity = 
  	2.0 * kine_.m_dQuaternionRotationdt.multConjugateToVector3( 
  	kine_.m_QuaternionRotation );
	
  return ( fileIn );
}




// ----------------------------------------------------------------------------
// Computes explicitly Idw/dt
Vector3 ParticleKinematics3D::computeExplicitDJomDt( Vector3 const & dw,
	double const* inertie ) const
{
  Vector3 Idw;
  Quaternion pConjugue = m_QuaternionRotation.Conjugate();
  
  // Changement de repère de R->R'
  Vector3 vTmp = pConjugue.multToVector3( ( dw , m_QuaternionRotation ) );
  
  // Calcul de I.pConjugue.dw.p dans R'
  Vector3 work; 
  work[0] = inertie[0]*vTmp[0] + inertie[1]*vTmp[1] + inertie[2]*vTmp[2];
  work[1] = inertie[1]*vTmp[0] + inertie[3]*vTmp[1] + inertie[4]*vTmp[2];
  work[2] = inertie[2]*vTmp[0] + inertie[4]*vTmp[1] + inertie[5]*vTmp[2];
  
  // Changement de repère R'->R
  Idw = m_QuaternionRotation.multToVector3( ( work , pConjugue ) );  
     
  return ( Idw );
}
