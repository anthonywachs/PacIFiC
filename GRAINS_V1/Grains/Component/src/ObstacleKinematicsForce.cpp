#include "ObstacleKinematicsForce.hh"
#include "GrainsMPIWrapper.hh"
#include "GrainsExec.hh"
#include <limits>
#include "Obstacle.hh"
#include "Torsor.hh"


// ----------------------------------------------------------------------------
// Default constructor
ObstacleKinematicsForce::ObstacleKinematicsForce() 
  : m_currentImposedForce( NULL )
{}




// ----------------------------------------------------------------------------
// Destructor
ObstacleKinematicsForce::~ObstacleKinematicsForce()
{}




// ----------------------------------------------------------------------------
// Adds an imposed force load to the obstacle kinematics
void ObstacleKinematicsForce::append( ObstacleImposedForce* oif )
{
  m_imposedForces.push_back( oif );
}




// ----------------------------------------------------------------------------
// Deletes all imposed force loads
void ObstacleKinematicsForce::clearAndDestroy()
{
  list<ObstacleImposedForce*>::iterator iter;

  for (iter=m_imposedForces.begin(); iter!=m_imposedForces.end(); iter++)
    delete *iter;
  m_imposedForces.clear();

  if ( m_currentImposedForce ) delete m_currentImposedForce;
}




// ----------------------------------------------------------------------------
// Composes the obstacle kinematics with another "higher level"
// force kinematics
void ObstacleKinematicsForce::Compose( ObstacleKinematicsForce const& other, 
    	Point3 const& centre )
{
  m_translationalVelocity += other.m_translationalVelocity;
}




// ----------------------------------------------------------------------------
// Computes the obstacle velocity and returns whether the obstacle 
// moved from t to t+dt
bool ObstacleKinematicsForce::ImposedMotion( double time, double dt, 
	Obstacle* obstacle )
{
  // Force load over [t,t+dt]
  double fin = time + dt;
  double dtt = 0.;
  if ( m_currentImposedForce ) 
  {
    if ( m_currentImposedForce->isActif( time, fin ) ) 
    {
      dtt = m_currentImposedForce->getTime( time, fin );
      m_translationalVelocity = *m_currentImposedForce->translationalVelocity( 
      	time, dtt, obstacle );
    }
    else 
    {
      delete m_currentImposedForce;
      m_currentImposedForce = NULL;
    }
  }
  else 
  {
    list<ObstacleImposedForce*>::iterator iter;
    for (iter=m_imposedForces.begin(); iter!=m_imposedForces.end(); iter++)
      if ( (*iter)->isActif( time, fin ) ) 
      {
	m_currentImposedForce = *iter;
	iter = m_imposedForces.erase( iter );
	iter = m_imposedForces.end();
      }
  }
  
  m_vitesseD = Norm( m_translationalVelocity );
  
  return ( m_vitesseD != 0. );
}




// ----------------------------------------------------------------------------
// Returns translational motion over dt 
Vector3 ObstacleKinematicsForce::getTranslation( double dt ) const
{
  return ( m_translationalVelocity * dt );
}




// ----------------------------------------------------------------------------
// Resets kinematics to 0
void ObstacleKinematicsForce::reset()
{
  m_translationalVelocity = 0.;
  m_vitesseD = 0.;
}




// ----------------------------------------------------------------------------
// Computes the total velocity of the obstacle using the arm lever
Vector3 ObstacleKinematicsForce::Velocity( const Vector3 &om ) const
{
  return ( m_translationalVelocity );
}
