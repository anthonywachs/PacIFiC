#include "ObstacleKinematicsVelocity.hh"
#include "SimpleObstacle.hh"
#include "Torsor.hh"
#include "Memento.hh"


// ----------------------------------------------------------------------------
// Default constructor
ObstacleKinematicsVelocity::ObstacleKinematicsVelocity():
  m_memento( NULL )
{
}




// ----------------------------------------------------------------------------
// Destructor
ObstacleKinematicsVelocity::~ObstacleKinematicsVelocity()
{
  if ( m_memento ) delete m_memento;
  clearAndDestroy();
}




// ----------------------------------------------------------------------------
// Adds an imposed velocity motion to the obstacle kinematics
void ObstacleKinematicsVelocity::append( ObstacleImposedVelocity* chargement )
{
  if ( m_imposedVelocities.empty( ))
    m_imposedVelocities.push_back( chargement );
  else 
  {
    list<ObstacleImposedVelocity*>::iterator c;
    for (c=m_imposedVelocities.begin(); 
    	c!=m_imposedVelocities.end() && **c<*chargement; c++)
      {}
    m_imposedVelocities.insert( c, chargement );
  }
}




// ----------------------------------------------------------------------------
// Deletes all imposed motions
void ObstacleKinematicsVelocity::clearAndDestroy()
{
  list<ObstacleImposedVelocity*>::iterator chargement;

  for (chargement=m_imposedVelocities.begin(); 
       chargement!=m_imposedVelocities.end(); chargement++) 
    delete *chargement;

  m_imposedVelocities.clear();
}




// ----------------------------------------------------------------------------
// Composes the obstacle kinematics with another "higher level" velocity 
// kinematics
void ObstacleKinematicsVelocity::Compose( 
	ObstacleKinematicsVelocity const& other, 
    	Vector3 const& lever )
{
  Matrix mat( other.m_QuaternionRotationOverDt );
  Vector3 rota = other.m_rotationOverTimeStep;
  double  d = Norm(rota);

  if ( d != 0. ) 
  {
    Vector3 vect = ( sin( d / 2.) / d ) * rota;
    m_QuaternionRotationOverDt = Quaternion( vect, cos( d / 2. ) );
  }
  else 
    m_QuaternionRotationOverDt = Quaternion( 0., 0., 0., 1. );
  
  // Pour des composites de composites de ... etc, vérifier
  // que l'écriture correcte ne serait pas:
  // m_rotationOverTimeStep += voisine.m_rotationOverTimeStep;
  // pour conserver l'aspect récursif
  // A discuter avec Gilles
  m_rotationOverTimeStep = other.m_rotationOverTimeStep;
  m_translationOverTimeStep = other.m_translationOverTimeStep 
  	+ ( ( mat * lever ) - lever );
	
  m_angularVelocity += other.m_angularVelocity;
  m_translationalVelocity += other.m_translationalVelocity 
  	+ ( other.m_angularVelocity ^ lever );
}




// ----------------------------------------------------------------------------
// Returns whether the obstacle moved from t to t+dt
bool ObstacleKinematicsVelocity::Deplacement( double time, double dt )
{
  Vector3 depl, rota, vt, vr;

  // Chargements de l'obstacle
  list<ObstacleImposedVelocity*>::iterator chargement;
  for (chargement=m_imposedVelocities.begin(); 
  	chargement!=m_imposedVelocities.end();chargement++) 
    if ( (*chargement)->isActif(time, dt) ) 
    {
      depl += (*chargement)->translationalDisplacement( time, dt ); 
      rota += (*chargement)->angularDisplacement( time, dt );
      vt += *(*chargement)->translationalVelocity( time, dt );
      vr += *(*chargement)->angularVelocity( time, dt );      
    } 

  double d = Norm( rota );
  if ( d != 0. ) 
  {
    Vector3 vect = ( sin( d /2. ) / d ) * rota;
    m_QuaternionRotationOverDt = Quaternion( vect, cos( d / 2. ) );
  }

  // Sur le pas de time, le mouvement de l'obstacle correspond à celui du
  // composite dont il fait partie plus son mouvement propre
  // Si le composite dont il fait partie n'a pas de mouvement imposé, 
  // m_translationOverTimeStep et m_rotationOverTimeStep sont mis à 0 par
  // la méthode ObstacleKinematicsVelocity::Compose et seuls les chargements 
  // propres de
  // l'obstacle sont pris en compte, ce qui justifie l'utilisation du +=  
  m_translationOverTimeStep += depl;
  m_rotationOverTimeStep += rota;
  
  m_angularVelocity += vr;
  m_translationalVelocity += vt; 

  return ( Norm(depl) != 0. || Norm(rota) != 0. );
}




// ----------------------------------------------------------------------------
// Returns a pointer to the rotation quaternion over dt
Quaternion const* ObstacleKinematicsVelocity::getQuaternionRotationOverDt() 
	const
{
  return ( &m_QuaternionRotationOverDt );
}




// ----------------------------------------------------------------------------
// Returns a pointer to the translation vector over dt
Vector3 const* ObstacleKinematicsVelocity::getTranslation() const
{
  return &m_translationOverTimeStep;
}




// ----------------------------------------------------------------------------
// Returns a pointer to the current angular velocity vector
Vector3 const* ObstacleKinematicsVelocity::getAngularVelocity() const
{
  return ( &m_angularVelocity );
}




// ----------------------------------------------------------------------------
// Returns a pointer to the current translational velocity vector
Vector3 const* ObstacleKinematicsVelocity::getTranslationalVelocity() const
{
  return ( &m_translationalVelocity );
}




// ----------------------------------------------------------------------------
// Resets kinematics to 0
void ObstacleKinematicsVelocity::reset()
{
  m_translationOverTimeStep = 0.;
  m_rotationOverTimeStep = 0.;
  m_translationalVelocity = 0.;
  m_angularVelocity = 0.;
}




// ----------------------------------------------------------------------------
// Sets the velocity and displacement using another velocity kinematics
void ObstacleKinematicsVelocity::set( ObstacleKinematicsVelocity& kine_ )
{
  m_translationOverTimeStep = kine_.m_translationOverTimeStep;
  m_rotationOverTimeStep = kine_.m_rotationOverTimeStep;
  m_translationalVelocity = kine_.m_translationalVelocity;
  m_angularVelocity = kine_.m_angularVelocity;
}




// ----------------------------------------------------------------------------
// Sets velocity
void ObstacleKinematicsVelocity::setVelocity( Vector3 const* vtrans,
  	Vector3 const* vrot )
{
  m_translationalVelocity = *vtrans ;
  m_angularVelocity = *vrot ;
}   




// ----------------------------------------------------------------------------
// Computes the total velocity of the obstacle using the arm lever
Vector3 ObstacleKinematicsVelocity::Velocity( Vector3 const& om ) const
{
  return ( m_translationalVelocity + (m_angularVelocity ^ om) );
}




// ----------------------------------------------------------------------------
// Output operator
ostream& operator << ( ostream& fileOut, 
    	ObstacleKinematicsVelocity const& kine_ )
{
  fileOut << "*ObstacleKinematicsVelocity\n";
  fileOut << kine_.m_translationOverTimeStep;
   fileOut << "*Chargement\n";
  fileOut << kine_.m_imposedVelocities.size() << '\n';
  list<ObstacleImposedVelocity*>::const_iterator chargement;
  for (chargement=kine_.m_imposedVelocities.begin(); 
       chargement!=kine_.m_imposedVelocities.end(); chargement++)
    fileOut << **chargement;

 return ( fileOut );
}




// ----------------------------------------------------------------------------
// Input operator
istream& operator >> ( istream& fileIn, ObstacleKinematicsVelocity& kine_ )
{
  string cle;
  fileIn >> cle;
  fileIn >> kine_.m_translationOverTimeStep;

  return ( fileIn );
}




// ----------------------------------------------------------------------------
// Returns whether there is an active angular motion imposed from t to t+dt
bool ObstacleKinematicsVelocity::activAngularMotion( double time, double dt ) 
	const
{
  bool rotation = false ;
  list<ObstacleImposedVelocity*>::const_iterator chargement;
  for (chargement=m_imposedVelocities.begin(); 
  	chargement!=m_imposedVelocities.end() && !rotation; chargement++)
    if ( (*chargement)->isActif( time, dt ) &&
    	( (*chargement)->getType() == "Rotation" ||
	(*chargement)->getType() == "RotationSinusoidale" ) ) rotation = true;
  
  return ( rotation ); 
}




// ----------------------------------------------------------------------------
// Saves obstacle kinematics state
void ObstacleKinematicsVelocity::saveState()
{
  if (!m_memento) m_memento = new ObstacleKinematicsMemento();

  m_memento->m_translationalVelocity = m_translationalVelocity;
  m_memento->m_angularVelocity = m_angularVelocity;
}




// ----------------------------------------------------------------------------
// Creates and returns obstacle kinematics state
ObstacleKinematicsMemento* ObstacleKinematicsVelocity::createState()
{
  ObstacleKinematicsMemento* memento_ = new ObstacleKinematicsMemento();
  
  memento_->m_translationalVelocity = m_translationalVelocity;
  memento_->m_angularVelocity = m_angularVelocity;
  
  return ( memento_ );
}




// ----------------------------------------------------------------------------
// Restores obstacle kinematics state
void ObstacleKinematicsVelocity::restoreState()
{
  m_translationalVelocity = m_memento->m_translationalVelocity;
  m_angularVelocity = m_memento->m_angularVelocity;
}




// ----------------------------------------------------------------------------
// Restores obstacle kinematics state
void ObstacleKinematicsVelocity::restoreState( 
	ObstacleKinematicsMemento const* memento_ )
{
  m_translationalVelocity = memento_->m_translationalVelocity;
  m_angularVelocity = memento_->m_angularVelocity;
} 




// ----------------------------------------------------------------------------
// Returns the list of velocity motions imposed on the obstacle
list<ObstacleImposedVelocity*> ObstacleKinematicsVelocity::getChargements() 
	const
{
  return ( m_imposedVelocities );
}
