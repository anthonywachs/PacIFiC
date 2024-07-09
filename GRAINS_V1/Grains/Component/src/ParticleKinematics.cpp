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
}




// ----------------------------------------------------------------------------
// Copy constructor
ParticleKinematics::ParticleKinematics( ParticleKinematics const& copy ) : 
  Kinematics( copy ), 
  m_memento( NULL ) 
{
  m_translationalVelocity = copy.m_translationalVelocity;
  m_angularVelocity = copy.m_angularVelocity;  
  m_dQuaternionRotationdt = copy.m_dQuaternionRotationdt;
  m_QuaternionRotation = copy.m_QuaternionRotation;
  m_timeIntegrationScheme = copy.m_timeIntegrationScheme->clone();
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
  // Time integration
  m_timeIntegrationScheme->Move( m_dUdt, m_translationalVelocity,
  	m_translationalMotionOverDt,
	m_dOmegadt, m_angularVelocity, m_averageAngularVelocityOverDt, 
	dt_particle_vel, dt_particle_disp );
  
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
  m_dQuaternionRotationdt = 0.5 * ( m_angularVelocity , m_QuaternionRotation );
  particle->Rotate( m_QuaternionRotationOverDt );

  return Norm( m_translationalMotionOverDt );
}




// ----------------------------------------------------------------------------
// Advances velocity over dt_particle_vel
void ParticleKinematics::advanceVelocity( double const& dt_particle_vel )
{
  // Time integration
  m_timeIntegrationScheme->advanceVelocity( m_dUdt, m_translationalVelocity,
  	m_dOmegadt, m_angularVelocity, dt_particle_vel );
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
  m_dQuaternionRotationdt = 0.;
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
  m_dQuaternionRotationdt = 0.5 * ( m_angularVelocity , m_QuaternionRotation );
}




// ----------------------------------------------------------------------------
// Sets the angular velocity
void ParticleKinematics::setAngularVelocity( double const& omx, 
	double const& omy, double const& omz )
{
  m_angularVelocity[X] = omx;
  m_angularVelocity[Y] = omy;  
  m_angularVelocity[Z] = omz; 
  m_dQuaternionRotationdt = 0.5 * ( m_angularVelocity , m_QuaternionRotation );
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
  m_memento->m_dQuaternionRotationdt = m_dQuaternionRotationdt;
}




// ----------------------------------------------------------------------------
// Creates and returns particle kinematics state
ParticleKinematicsMemento* ParticleKinematics::createState()
{
  ParticleKinematicsMemento* memento_ = new ParticleKinematicsMemento();
  
  memento_->m_QuaternionRotation = m_QuaternionRotation;
  memento_->m_translationalVelocity = m_translationalVelocity;
  memento_->m_dQuaternionRotationdt = m_dQuaternionRotationdt;
  
  return ( memento_ );
}




// ----------------------------------------------------------------------------
// Restores particle kinematics state
void ParticleKinematics::restoreState()
{
  m_QuaternionRotation = m_memento->m_QuaternionRotation;
  m_translationalVelocity = m_memento->m_translationalVelocity;
  m_dQuaternionRotationdt = m_memento->m_dQuaternionRotationdt;
  m_angularVelocity = 2.0 * m_dQuaternionRotationdt.multConjugateToVector3( 
  	m_QuaternionRotation );
}




// ----------------------------------------------------------------------------
// Restores particle kinematics state
void ParticleKinematics::restoreState( 
	ParticleKinematicsMemento const* memento_ )
{
  m_QuaternionRotation = memento_->m_QuaternionRotation;
  m_translationalVelocity = memento_->m_translationalVelocity;
  m_dQuaternionRotationdt = memento_->m_dQuaternionRotationdt;
  m_angularVelocity = 2.0 * m_dQuaternionRotationdt.multConjugateToVector3( 
  	m_QuaternionRotation );
}




// ----------------------------------------------------------------------------
// Output operator
ostream& operator << ( ostream& fileOut, ParticleKinematics const& cinematique )
{
  fileOut << cinematique.className() << '\n';
  fileOut << cinematique.m_translationalVelocity
	  << cinematique.m_QuaternionRotation
	  << cinematique.m_dQuaternionRotationdt;
	  
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
  m_dQuaternionRotationdt.writeQuaternion( fileOut ); 
}




// ---------------------------------------------------------------------------
// Writes particle kinematics in an output stream with a high
// precision and 2014 format
void ParticleKinematics::writeParticleKinematics2014( ostream& fileOut ) const
{
  m_translationalVelocity.writeGroup3( fileOut ); 
  fileOut << " ";  
  m_QuaternionRotation.writeQuaternion( fileOut ); 
  fileOut << " "; 
  m_dQuaternionRotationdt.writeQuaternion( fileOut );
  m_timeIntegrationScheme->writeParticleKinematics2014( fileOut,
    	m_dUdt, m_dOmegadt ); 
}




// ---------------------------------------------------------------------------
// Writes particle kinematics in an output stream with a binary and 2014 format
void ParticleKinematics::writeParticleKinematics2014_binary( ostream &fileOut )
{
  m_translationalVelocity.writeGroup3_binary( fileOut ); 
  m_QuaternionRotation.writeQuaternion_binary( fileOut ); 
  m_dQuaternionRotationdt.writeQuaternion_binary( fileOut );
  m_timeIntegrationScheme->writeParticleKinematics2014_binary( fileOut,
    	m_dUdt, m_dOmegadt );      
}




// ----------------------------------------------------------------------------
// Input operator
istream& operator >> ( istream& fileIn, ParticleKinematics& cinematique )
{
  fileIn >> cinematique.m_translationalVelocity
	 >> cinematique.m_QuaternionRotation
	 >> cinematique.m_dQuaternionRotationdt;
  cinematique.m_angularVelocity = 
  	2.0 * cinematique.m_dQuaternionRotationdt.multConjugateToVector3( 
  	cinematique.m_QuaternionRotation );
	  
  return ( fileIn );
}




// ----------------------------------------------------------------------------
// Reads particle kinematics from a stream in a binary form in the 2014 format
void ParticleKinematics::readParticleKinematics2014( istream &StreamIN )
{
  StreamIN >> m_translationalVelocity
	 >> m_QuaternionRotation
	 >> m_dQuaternionRotationdt;

  m_angularVelocity = 2.0 * m_dQuaternionRotationdt.multConjugateToVector3( 
  	m_QuaternionRotation );
	
  m_timeIntegrationScheme->readParticleKinematics2014( StreamIN,
    	m_dUdt, m_dOmegadt );     
} 




// ----------------------------------------------------------------------------
// Reads particle kinematics from a stream in a binary form in the 2014 format
void ParticleKinematics::readParticleKinematics2014_binary( istream &StreamIN )
{
  m_translationalVelocity.readGroup3_binary( StreamIN );
  m_QuaternionRotation.readQuaternion_binary( StreamIN );
  m_dQuaternionRotationdt.readQuaternion_binary( StreamIN );

  m_angularVelocity = 2.0 * m_dQuaternionRotationdt.multConjugateToVector3( 
  	m_QuaternionRotation );
	
  m_timeIntegrationScheme->readParticleKinematics2014_binary( StreamIN,
    	m_dUdt, m_dOmegadt );   
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
  return ( solid::Group3::m_sizeofGroup3 + 2 * Quaternion::m_sizeofQuaternion
  	+  m_timeIntegrationScheme->get_numberOfBytes() ); 
}
