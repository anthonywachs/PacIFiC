#include "ParticleKinematics.hh"
#include "Particle.hh"
#include "Torsor.hh"
#include "Memento.hh"
#include "TimeIntegratorBuilderFactory.hh"
#include "TimeIntegrator.hh"


// ----------------------------------------------------------------------------
// Constructeur par defaut
ParticleKinematics::ParticleKinematics() : 
  m_timeIntegrationScheme( NULL ),
  m_memento( NULL ) 
{
  m_timeIntegrationScheme = TimeIntegratorBuilderFactory::create();
  m_coupling_factor = 1.;    
}




// ----------------------------------------------------------------------------
// Copy constructor
ParticleKinematics::ParticleKinematics( ParticleKinematics const& copy ) : 
  Kinematics( copy ), 
  m_memento( NULL ) 
{
  m_translationalVelocity = copy.m_translationalVelocity;
  m_angularVelocity = copy.m_angularVelocity;
  m_angularVelocity_bf = copy.m_angularVelocity_bf;    
  m_QuaternionRotation = copy.m_QuaternionRotation; 
  m_translationalMotionOverDt = copy.m_translationalMotionOverDt;  
  m_averageAngularVelocityOverDt = copy.m_averageAngularVelocityOverDt;  
  m_QuaternionRotationOverDt = copy.m_QuaternionRotationOverDt;   
  m_timeIntegrationScheme = copy.m_timeIntegrationScheme->clone();
  m_coupling_factor = copy.m_coupling_factor;
  if ( !copy.m_memento )
  {
    m_memento = new ParticleKinematicsMemento();
    m_memento->m_QuaternionRotation = m_QuaternionRotation;
    m_memento->m_translationalVelocity = m_translationalVelocity;
    m_memento->m_angularVelocity = m_angularVelocity;  
  }
}




// ----------------------------------------------------------------------------
// Destructor
ParticleKinematics::~ParticleKinematics()
{
  if ( m_timeIntegrationScheme ) delete m_timeIntegrationScheme;
  if ( m_memento ) delete m_memento;
}




// ----------------------------------------------------------------------------
// Integrates Newton's law and moves the particle
double ParticleKinematics::Move( Particle* particle, 
	double const& dt_particle_vel, 
    	double const& dt_particle_disp )
{
  // Quaternion and torque in body-fixed coordinates system
  Quaternion pConjugue = m_QuaternionRotation.Conjugate();
  Vector3 torque_bf = pConjugue.multToVector3( 
    	( *(particle->getTorque()) , m_QuaternionRotation ) ); 
//  Vector3 torque_bf( 0., 0.5, 0.);	 

  // Time integration
  m_timeIntegrationScheme->Move( particle, this, m_coupling_factor, torque_bf, 
  	m_translationalVelocity, m_translationalMotionOverDt,
	m_angularVelocity_bf, m_averageAngularVelocityOverDt, 
	dt_particle_vel, dt_particle_disp );

  // Rotate back angular velocity to space-fixed coordinates system
  m_angularVelocity = m_QuaternionRotation.multToVector3( 
  	( m_angularVelocity_bf , pConjugue ) );
  m_averageAngularVelocityOverDt = m_QuaternionRotation.multToVector3( 
  	( m_averageAngularVelocityOverDt , pConjugue ) );	  	
  
  // Translational motion
  particle->Translate( m_translationalMotionOverDt );
  
  // Angular motion 
  double nOmega = Norm( m_averageAngularVelocityOverDt );
  if ( nOmega > EPSILON ) 
  {
    double c = cos( nOmega * dt_particle_disp / 2. );
    double s = sin( nOmega * dt_particle_disp / 2. );
    Vector3 t;
    t = ( s * 1. / nOmega ) * m_averageAngularVelocityOverDt;
    m_QuaternionRotationOverDt.setQuaternion( t, c );
  } 
  else 
     m_QuaternionRotationOverDt.setQuaternion( 0., 0., 0., 1. );
   
  m_QuaternionRotation = m_QuaternionRotationOverDt * m_QuaternionRotation;
  particle->Rotate( m_QuaternionRotationOverDt );

  return Norm( m_translationalMotionOverDt );
}




// ----------------------------------------------------------------------------
// Advances velocity over dt_particle_vel
void ParticleKinematics::advanceVelocity( Particle* particle, 
	double const& dt_particle_vel )
{
  // Quaternion and torque in body-fixed coordinates system
  Quaternion pConjugue = m_QuaternionRotation.Conjugate();
  Vector3 torque_bf = pConjugue.multToVector3( 
    	( *(particle->getTorque()) , m_QuaternionRotation ) ); 
//  Vector3 torque_bf( 0., 0.5, 0.);

  // Time integration
  m_timeIntegrationScheme->advanceVelocity( particle, this, 
  	m_coupling_factor, torque_bf, m_translationalVelocity, 
	m_angularVelocity_bf, dt_particle_vel );

  // Rotate back angular velocity to space-fixed coordinates system
  m_angularVelocity = m_QuaternionRotation.multToVector3( 
  	( m_angularVelocity_bf , pConjugue ) );  	
}




// ----------------------------------------------------------------------------
// Returns the rotation quaternion
Quaternion const* ParticleKinematics::getQuaternionRotation() const
{
  return ( &m_QuaternionRotation );
}




// ----------------------------------------------------------------------------
// Returns angular velocity
Vector3 const* ParticleKinematics::getAngularVelocity() const
{
  return ( &m_angularVelocity ); 
}




// ----------------------------------------------------------------------------
// Returns translational velocity
Vector3 const* ParticleKinematics::getTranslationalVelocity() const
{
  return ( &m_translationalVelocity );
}




// ----------------------------------------------------------------------------
// Resets kinematics to 0
void ParticleKinematics::reset()
{
  m_translationalVelocity = 0.;
  m_angularVelocity = 0.; 
  m_QuaternionRotation.setQuaternion( 0., 0., 0., 1. );
}




// ----------------------------------------------------------------------------
// Sets the translation velocity
void ParticleKinematics::setTranslationalVelocity( Vector3 const& vtrans )
{
  m_translationalVelocity = vtrans;
}




// ----------------------------------------------------------------------------
// Sets the translation velocity
void ParticleKinematics::setTranslationalVelocity( double const& vx, 
	double const& vy, double const& vz )
{
  m_translationalVelocity[X] = vx;
  m_translationalVelocity[Y] = vy;  
  m_translationalVelocity[Z] = vz;  
}




// ----------------------------------------------------------------------------
// Sets the angular velocity
void ParticleKinematics::setAngularVelocity( Vector3 const& omega )
{
  m_angularVelocity = omega;
  m_angularVelocity_bf = m_angularVelocity;   
}




// ----------------------------------------------------------------------------
// Sets the angular velocity
void ParticleKinematics::setAngularVelocity( double const& omx, 
	double const& omy, double const& omz )
{
  m_angularVelocity[X] = omx;
  m_angularVelocity[Y] = omy;  
  m_angularVelocity[Z] = omz;
  m_angularVelocity_bf = m_angularVelocity; 
}




// ----------------------------------------------------------------------------
// Returns the total velocity U + om x R given R 
Vector3 ParticleKinematics::Velocity( Vector3 const& lev ) const
{
  return ( m_translationalVelocity + ( m_angularVelocity ^ lev ) );
}




// ----------------------------------------------------------------------------
// Saves particle kinematics state
void ParticleKinematics::saveState()
{
  if ( !m_memento ) m_memento = new ParticleKinematicsMemento();

  m_memento->m_QuaternionRotation = m_QuaternionRotation;
  m_memento->m_translationalVelocity = m_translationalVelocity;
  m_memento->m_angularVelocity = m_angularVelocity;
}




// ----------------------------------------------------------------------------
// Creates and returns particle kinematics state
ParticleKinematicsMemento* ParticleKinematics::createState()
{
  ParticleKinematicsMemento* memento_ = new ParticleKinematicsMemento();
  
  memento_->m_QuaternionRotation = m_QuaternionRotation;
  memento_->m_translationalVelocity = m_translationalVelocity;
  memento_->m_angularVelocity = m_angularVelocity;
  
  return ( memento_ );
}




// ----------------------------------------------------------------------------
// Restores particle kinematics state
void ParticleKinematics::restoreState()
{
  m_QuaternionRotation = m_memento->m_QuaternionRotation;
  m_translationalVelocity = m_memento->m_translationalVelocity;
  m_angularVelocity = m_memento->m_angularVelocity;
  Quaternion pConjugue = m_QuaternionRotation.Conjugate();
  m_angularVelocity_bf = pConjugue.multToVector3( 
	( m_angularVelocity , m_QuaternionRotation ) );     
}




// ----------------------------------------------------------------------------
// Restores particle kinematics state
void ParticleKinematics::restoreState( 
	ParticleKinematicsMemento const* memento_ )
{
  m_QuaternionRotation = memento_->m_QuaternionRotation;
  m_translationalVelocity = memento_->m_translationalVelocity;
  m_angularVelocity = memento_->m_angularVelocity;
  Quaternion pConjugue = m_QuaternionRotation.Conjugate();
  m_angularVelocity_bf = pConjugue.multToVector3( 
	( m_angularVelocity , m_QuaternionRotation ) );   
}




// ----------------------------------------------------------------------------
// Output operator
ostream& operator << ( ostream& fileOut, ParticleKinematics const& cinematique )
{
  fileOut << cinematique.className() << '\n';
  fileOut << cinematique.m_translationalVelocity
	  << cinematique.m_QuaternionRotation
	  << cinematique.m_angularVelocity;
	  
  return ( fileOut );
}




// ---------------------------------------------------------------------------
// Writes particle kinematics in an output stream with a high precision format
void ParticleKinematics::writeParticleKinematics( ostream& fileOut ) const
{
  fileOut << className() << endl;
  m_translationalVelocity.writeGroup3( fileOut ); 
  fileOut << endl;  
  m_QuaternionRotation.writeQuaternion( fileOut ); 
  fileOut << endl; 
  m_angularVelocity.writeGroup3( fileOut ); 
}




// ---------------------------------------------------------------------------
// Writes particle kinematics in an output stream with a high
// precision and 2014 format
void ParticleKinematics::writeParticleKinematics2014( ostream& fileOut, 
    	Particle const* particle ) const
{
  m_translationalVelocity.writeGroup3( fileOut ); 
  fileOut << " ";  
  m_QuaternionRotation.writeQuaternion( fileOut ); 
  fileOut << " "; 
  m_angularVelocity.writeGroup3( fileOut );
  m_timeIntegrationScheme->writeParticleKinematics2014( fileOut, particle ); 
}




// ---------------------------------------------------------------------------
// Writes particle kinematics in an output stream with a binary and 2014 format
void ParticleKinematics::writeParticleKinematics2014_binary( ostream &fileOut, 
    	Particle const* particle )
{
  m_translationalVelocity.writeGroup3_binary( fileOut ); 
  m_QuaternionRotation.writeQuaternion_binary( fileOut ); 
  m_angularVelocity.writeGroup3_binary( fileOut );
  m_timeIntegrationScheme->writeParticleKinematics2014_binary( fileOut,
  	particle );      
}




// ----------------------------------------------------------------------------
// Input operator
istream& operator >> ( istream& fileIn, ParticleKinematics& cinematique )
{
  fileIn >> cinematique.m_translationalVelocity
	 >> cinematique.m_QuaternionRotation
	 >> cinematique.m_angularVelocity;
  Quaternion pConjugue = cinematique.m_QuaternionRotation.Conjugate();
  cinematique.m_angularVelocity_bf = pConjugue.multToVector3( 
	( cinematique.m_angularVelocity , cinematique.m_QuaternionRotation ) );
  return ( fileIn );
}




// ----------------------------------------------------------------------------
// Reads particle kinematics from a stream in a binary form in the 2014 format
void ParticleKinematics::readParticleKinematics2014( istream &StreamIN,
    	Particle* particle )
{
  StreamIN >> m_translationalVelocity
	 >> m_QuaternionRotation
	 >> m_angularVelocity;
  Quaternion pConjugue = m_QuaternionRotation.Conjugate();
  m_angularVelocity_bf = pConjugue.multToVector3( 
	( m_angularVelocity , m_QuaternionRotation ) ); 	
  m_timeIntegrationScheme->readParticleKinematics2014( StreamIN, particle );
} 




// ----------------------------------------------------------------------------
// Reads particle kinematics from a stream in a binary form in the 2014 format
void ParticleKinematics::readParticleKinematics2014_binary( istream &StreamIN,
    	Particle* particle )
{
  m_translationalVelocity.readGroup3_binary( StreamIN );
  m_QuaternionRotation.readQuaternion_binary( StreamIN );
  m_angularVelocity.readGroup3_binary( StreamIN );
  Quaternion pConjugue = m_QuaternionRotation.Conjugate();
  m_angularVelocity_bf = pConjugue.multToVector3( 
	( m_angularVelocity , m_QuaternionRotation ) ); 	
  m_timeIntegrationScheme->readParticleKinematics2014_binary( StreamIN,
  	particle );   
} 




// ----------------------------------------------------------------------------
// Sets the rotation quaternion
void ParticleKinematics::setQuaternionRotation( double const& vecteur0, 
	double const& vecteur1, 
	double const& vecteur2, 
	double const& scalaire )
{
  m_QuaternionRotation.setQuaternion( vecteur0, vecteur1, vecteur2, scalaire );
}




// ----------------------------------------------------------------------------
// Sets the rotation quaternion
void ParticleKinematics::setQuaternionRotation( Quaternion const& qrot )
{
  m_QuaternionRotation = qrot;
}




// ----------------------------------------------------------------------------
// Copies kinematics at time t-2dt as a 1D array of 12 scalars
// (translational velocity, angular velocity, variation of translational
// momentum, variation of angular momentum)
void ParticleKinematics::copyKinematicsNm2( double *vit, int i ) const
{
  m_timeIntegrationScheme->copyKinematicsNm2( vit, i );
}




// ----------------------------------------------------------------------------
// Sets kinematics at time t-2dt as a 1D array of 12 scalars
// (translational velocity, angular velocity, variation of translational
// momentum, variation of angular momentum)
void ParticleKinematics::setKinematicsNm2( double const* tab )
{
  m_timeIntegrationScheme->setKinematicsNm2( tab );
} 




// ----------------------------------------------------------------------------
// Returns the number of bytes of the ParticleKinematics when written in a 
// binary format to an output stream
size_t ParticleKinematics::get_numberOfBytes() const
{
  return ( 2 * solid::Group3::m_sizeofGroup3 + Quaternion::m_sizeofQuaternion
  	+  m_timeIntegrationScheme->get_numberOfBytes() );
}




// ----------------------------------------------------------------------------
// Sets time integration scheme using the macro variable GrainsExec::m_TIScheme
void ParticleKinematics::setTimeIntegrationScheme()
{
  if ( m_timeIntegrationScheme ) delete m_timeIntegrationScheme;  
  m_timeIntegrationScheme = TimeIntegratorBuilderFactory::create();
}




// ----------------------------------------------------------------------------
// Sets momentum equation coupling factor
void ParticleKinematics::setCouplingFactor( double const& particle_density )
{
  m_coupling_factor = 1.;
  if ( Particle::getFluidCorrectedAcceleration() )
    m_coupling_factor -= Particle::getFluidDensity() / particle_density;   
}
