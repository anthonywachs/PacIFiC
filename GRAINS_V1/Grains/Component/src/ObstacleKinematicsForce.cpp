#include "ObstacleKinematicsForce.hh"
#include "GrainsMPIWrapper.hh"
#include "GrainsExec.hh"
#include <limits>
#include "Obstacle.hh"
#include "Torsor.hh"


// ----------------------------------------------------------------------------
// Default constructor
ObstacleKinematicsForce::ObstacleKinematicsForce() 
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
}




// ----------------------------------------------------------------------------
// Composes the obstacle kinematics with another "higher level"
// force kinematics
void ObstacleKinematicsForce::Compose( ObstacleKinematicsForce const& other )
{
  m_translationalVelocity += other.m_translationalVelocity;
  m_translationOverTimeStep += other.m_translationOverTimeStep;
}




// ----------------------------------------------------------------------------
// Computes the obstacle velocity and returns whether the obstacle 
// moved from time - dt to time
bool ObstacleKinematicsForce::ImposedMotion( double time, double dt, 
	Obstacle* obstacle )
{
  // Force load over [time-dt,time]
  double dtt = 0.;
  bool found = false;
  ObstacleImposedForce* currentImposedForce = NULL;
  list<ObstacleImposedForce*>::iterator iter;
  
  // Note: only a single imposed force can be active over [time-dt,time]
  // Therefore we impose the first force that is active in the list, it is up 
  // to the user to impose forces properly
  for (iter=m_imposedForces.begin(); iter!=m_imposedForces.end() && !found; 
  	iter++)
    if ( (*iter)->isActif( time - dt, time, dt, dtt ) ) 
    {
      currentImposedForce = *iter;
      currentImposedForce->translationalVelocity( 
      	time, dt, dtt, m_translationalVelocity, m_translationOverTimeStep,
	obstacle );
    }  
  
  m_vitesseD = Norm( m_translationalVelocity );
  
  return ( m_vitesseD != 0. );
}




// ----------------------------------------------------------------------------
// Returns translational motion over dt 
Vector3 ObstacleKinematicsForce::getTranslation( double dt ) const
{
  return ( m_translationOverTimeStep );
}




// ----------------------------------------------------------------------------
// Resets kinematics to 0
void ObstacleKinematicsForce::reset()
{
  m_translationalVelocity = 0.;
  m_translationOverTimeStep = 0.;
  m_vitesseD = 0.;
}




// ----------------------------------------------------------------------------
// Computes the total velocity of the obstacle using the arm lever
Vector3 ObstacleKinematicsForce::Velocity( const Vector3 &om ) const
{
  return ( m_translationalVelocity );
}




// ----------------------------------------------------------------------------
// Returns a pointer to the current translational velocity vector
Vector3 const* ObstacleKinematicsForce::getTranslationalVelocity() const
{
  return ( &m_translationalVelocity );
}
