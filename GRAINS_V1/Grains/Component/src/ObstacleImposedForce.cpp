#include "ObstacleImposedForce.hh"
#include "GrainsMPIWrapper.hh"
#include "GrainsExec.hh"
#include "Obstacle.hh"
#include <stdlib.h>


// ----------------------------------------------------------------------------
// Default constructor
ObstacleImposedForce::ObstacleImposedForce()
{
  m_ObstacleName = "Undefined";
  m_type = "Undefined";
  m_tstart = 0.; 
  m_tend = 0.;
  m_force = Vector3Nul;
  m_prev = Vector3Nul;
  m_mass = 0.;
  m_direction = Vector3Nul;
  m_translationalVelocity = Vector3Nul;
  m_freqX = 0.;
  m_freqY = 0.;
  m_freqZ = 0.;
  m_phase = 0.;
  m_prev = Vector3Nul;
}




// ----------------------------------------------------------------------------
// Constructor with an XML node as input parameter
ObstacleImposedForce::ObstacleImposedForce( DOMNode* root, double dt, 
	int rank, size_t& error )
{
  m_ObstacleName = ReaderXML::getNodeAttr_String( root, "ObstacleName" );
  
  DOMNode* nTimeInterval = ReaderXML::getNode( root, "TimeInterval" );
  m_tstart = ReaderXML::getNodeAttr_Double( nTimeInterval, "Start" );
  m_tend = ReaderXML::getNodeAttr_Double( nTimeInterval, "End" );
  
  DOMNode* force    = ReaderXML::getNode( root, "Amplitude" );
  DOMNode* nVector3 = ReaderXML::getNode( root, "Vector3" );
  DOMNode* property = ReaderXML::getNode( root, "Property" );
  m_mass	    = ReaderXML::getNodeAttr_Double( property, "Masse" );

  m_direction[X] = ReaderXML::getNodeAttr_Double( nVector3, "X" );
  m_direction[Y] = ReaderXML::getNodeAttr_Double( nVector3, "Y" );
  m_direction[Z] = ReaderXML::getNodeAttr_Double( nVector3, "Z" );

  m_type            = ReaderXML::getNodeAttr_String( root, "Type" );

  m_force[X]   	    = ReaderXML::getNodeAttr_Double( force, "AX" );
  m_force[Y]   	    = ReaderXML::getNodeAttr_Double( force, "AY" );
  m_force[Z]   	    = ReaderXML::getNodeAttr_Double( force, "AZ" );

  if ( m_type == "Cyclic")
  {
    DOMNode* frequence  = ReaderXML::getNode( root, "Frequence" );
    m_phase             = ReaderXML::getNodeAttr_Double( frequence, "Phi" );
    m_freqX        = ReaderXML::getNodeAttr_Double( frequence, "FX" );
    m_freqY        = ReaderXML::getNodeAttr_Double( frequence, "FY" );
    m_freqZ        = ReaderXML::getNodeAttr_Double( frequence, "FZ" );
    m_phase            *= PI / 180.;
  }

  if ( rank == 0 )
  {
    cout << "Chargement en force sur " << m_ObstacleName << endl;
    cout << "Type de chargement : " << m_type << endl;
    cout << "Amplitude de la force = " << m_force[X] << "\t" << m_force[Y] 
	 << "\t" << m_force[Z] << endl;
    cout << "Direction de la force = " << m_direction[X] << "\t"
	 << m_direction[Y] << "\t" << m_direction[Z] << endl;
    if ( m_type == "Cyclic" )
      cout << "Frequence : FX = " << m_freqX << "\tFY = " << m_freqY
           << "\tFZ = " << m_freqZ << endl;
    cout << "   Temps de depart = " << m_tstart << endl;
    cout << "   Temps de fin = " << m_tend << endl;
  }
}




// ----------------------------------------------------------------------------
// Destructor
ObstacleImposedForce::~ObstacleImposedForce()
{}




// ----------------------------------------------------------------------------
// Returns obstacle name
string ObstacleImposedForce::getNom() const
{
  return ( m_ObstacleName );
}




// ----------------------------------------------------------------------------
// Returns the remaining active time interval of the imposed motion
double ObstacleImposedForce::getTime( double debut, double fin ) const
{
  double activtimeint = fin - debut;

  if ( debut < m_tstart ) activtimeint -= ( m_tstart - debut );
  if ( m_tend < fin ) activtimeint -= ( fin - m_tend );

  return ( activtimeint );
}




// ----------------------------------------------------------------------------
// Returns whether the imposed motion is activ at time t
bool ObstacleImposedForce::isActif( double t, double dt ) const 
{
  return ( t > m_tstart - dt * 1.e-5  && t < m_tend + dt * 1.e-5 );
}




// ----------------------------------------------------------------------------
// Returns whether the imposed motion is completed at time t
bool ObstacleImposedForce::isCompleted( double t, double dt ) const 
{
  return ( t > m_tend + dt * 1.e-5 );
}




// ----------------------------------------------------------------------------
// Creates and reads the imposed force features from an input stream
ObstacleImposedForce* ObstacleImposedForce::read( istream& fileIn )
{
  ObstacleImposedForce* chargement;
  chargement = new ObstacleImposedForce();

  fileIn >> chargement->m_ObstacleName
	>> chargement->m_tstart >> chargement->m_tend
	>> chargement->m_force 
	>> chargement->m_mass
	>> chargement->m_direction;
  
  return ( chargement );
}




// ----------------------------------------------------------------------------
// Returns the imposed force
Vector3 ObstacleImposedForce::getForce() const
{
  return ( m_force );
}




// ----------------------------------------------------------------------------
// Returns the obstacle virtual mass
double ObstacleImposedForce::getMass() const
{
  return ( m_mass );
}




// ----------------------------------------------------------------------------
// Returns the direction of displacement
Vector3 const* ObstacleImposedForce::getDirection() const
{
  return ( &m_direction );
}




// ----------------------------------------------------------------------------
// Returns the imposed force type
string ObstacleImposedForce::getType() const
{
  return ( m_type );
}




// ----------------------------------------------------------------------------
// Returns the translational velocity at time t 
Vector3 const* ObstacleImposedForce::translationalVelocity( double time, 
	double dt, Obstacle* obstacle )
{
  Vector3 center = *obstacle->getPosition();
  GrainsMPIWrapper* wrapper = GrainsExec::getComm();
  // Somme des forces sur l'obstacle
  Torsor const* somme  = obstacle->getTorsor();
  Vector3 const* forces = somme->getForce();
  Vector3 force = *forces;
  force[X] = wrapper->sum_DOUBLE_master( force[X] ); 
  force[Y] = wrapper->sum_DOUBLE_master( force[Y] ); 
  force[Z] = wrapper->sum_DOUBLE_master( force[Z] ); 

  force[X] = wrapper->Broadcast_DOUBLE( force[X] );
  force[Y] = wrapper->Broadcast_DOUBLE( force[Y] );
  force[Z] = wrapper->Broadcast_DOUBLE( force[Z] );
  
  Vector3 dforce;
  Vector3 depl;
  Vector3 trans;
  if ( m_type == "Translation" )
  {
    dforce = m_force - force;
    dforce[X] *= m_direction[X]; 
    dforce[Y] *= m_direction[Y]; 
    dforce[Z] *= m_direction[Z]; 
    depl = 0.5 *( dt * dt / m_mass ) * dforce; 
    m_translationalVelocity = depl / dt;
  }
  else if ( m_type == "Cyclic" )
  {
    dforce = cyclicForce( time ) - force;
    dforce[X] *= m_direction[X]; 
    dforce[Y] *= m_direction[Y]; 
    dforce[Z] *= m_direction[Z]; 
    trans = 0.5 *( dt * dt / m_mass ) * dforce; 
    depl = trans - m_prev;
    m_translationalVelocity = depl / dt;
    // t-dt 
    m_prev = trans;
  }
  
  return ( &m_translationalVelocity );
}




// ----------------------------------------------------------------------------
// Returns the imposed force in cyclic mode 
Vector3 ObstacleImposedForce::cyclicForce( double time ) const
{
  Vector3 cycForce;
  cycForce[X] = m_force[X] * sin( 2. * PI * m_freqX * 
      	( time - m_tstart ) );
  cycForce[Y] = m_force[Y] * sin( 2. * PI * m_freqY *
      	( time - m_tstart ) + m_phase );
  cycForce[Z] = m_force[Z] * sin( 2. * PI * m_freqZ *
      	( time - m_tstart ) + m_phase );

  return ( cycForce );
}
