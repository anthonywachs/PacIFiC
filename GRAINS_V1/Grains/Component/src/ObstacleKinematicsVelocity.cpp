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
void ObstacleKinematicsVelocity::append( ObstacleImposedVelocity* oiv )
{
  if ( m_imposedVelocities.empty( ))
    m_imposedVelocities.push_back( oiv );
  else 
  {
    list<ObstacleImposedVelocity*>::iterator c;
    for (c=m_imposedVelocities.begin(); 
    	c!=m_imposedVelocities.end() && **c<*oiv; c++) {}
    m_imposedVelocities.insert( c, oiv );
  }
}




// ----------------------------------------------------------------------------
// Deletes all imposed motions
void ObstacleKinematicsVelocity::clearAndDestroy()
{
  list<ObstacleImposedVelocity*>::iterator il;

  for (il=m_imposedVelocities.begin(); il!=m_imposedVelocities.end(); il++) 
    delete *il;

  m_imposedVelocities.clear();
}




// ----------------------------------------------------------------------------
// Composes the obstacle kinematics with another "higher level" velocity 
// kinematics
void ObstacleKinematicsVelocity::Compose( 
	ObstacleKinematicsVelocity const& other, 
    	Vector3 const& lever )
{  
  m_rotationOverTimeStep += other.m_rotationOverTimeStep;
  m_QuaternionRotationOverDt = 
  	other.m_QuaternionRotationOverDt * m_QuaternionRotationOverDt;
  m_translationOverTimeStep = other.m_translationOverTimeStep 
  	+ ( other.m_QuaternionRotationOverDt.rotateVector( lever ) - lever );
	
  m_angularVelocity += other.m_angularVelocity;
  m_translationalVelocity += other.m_translationalVelocity 
  	+ ( other.m_angularVelocity ^ lever );
}




// ----------------------------------------------------------------------------
// Updates the obstacle translational and angular velocity at time time
// and translational and angular motion from time - dt to time and returns 
// whether the obstacle moved from time - dt to time
bool ObstacleKinematicsVelocity::ImposedMotion( double time, double dt, 
	Point3 const& cg )
{
  Vector3 depl, rota, vt, vr;
  Quaternion qrot;
  double dtt = 0.;

  // Imposed kinematics
  list<ObstacleImposedVelocity*>::iterator il;
  for (il=m_imposedVelocities.begin();il!=m_imposedVelocities.end();il++) 
    if ( (*il)->isActif( time - dt, time, dt, dtt ) ) 
    {
      depl += (*il)->translationalMotion( time, dt, dtt, cg ); 
      rota += (*il)->angularMotion( time, dt, dtt );
      vt += *(*il)->translationalVelocity( time, dt, cg );
      vr += *(*il)->angularVelocity( time, dt );      
    } 

  // Rotation quaternion over dt
  double d = Norm( rota );
  if ( d != 0. ) 
  {
    Vector3 vect = ( sin( d / 2. ) / d ) * rota;
    qrot.setQuaternion( vect, cos( d / 2. ) );
  }

  // The obstacle motion is the sum of its own and that of the composite
  // obstacle it belongs. The motion of the composite obstacle it belongs
  // to is added by the Compose method 
  m_translationOverTimeStep += depl;
  m_rotationOverTimeStep += rota; 
  m_QuaternionRotationOverDt = qrot * m_QuaternionRotationOverDt;   
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
  m_QuaternionRotationOverDt.setQuaternion( 0., 0., 0., 1. );
  m_translationalVelocity = 0.;
  m_angularVelocity = 0.;
}




// ----------------------------------------------------------------------------
// Sets the velocity and motion using another velocity kinematics
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
  return ( m_translationalVelocity + ( m_angularVelocity ^ om ) );
}




// ----------------------------------------------------------------------------
// Returns whether there is an active angular motion imposed from time - dt 
// to time
bool ObstacleKinematicsVelocity::activeAngularMotion( double time, double dt ) 
	const
{
  bool rotation = false ;
  double dtt = 0.;
  list<ObstacleImposedVelocity*>::const_iterator il;
  for (il=m_imposedVelocities.cbegin(); 
  	il!=m_imposedVelocities.cend() && !rotation; il++)
    if ( (*il)->isActif( time - dt, time, dt, dtt ) 
    	&& ( (*il)->getType() == "Rotation" ) ) 
      rotation = true;
  
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
