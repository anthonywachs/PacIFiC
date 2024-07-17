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
  m_force_amplitude = Vector3Null;
  m_force = Vector3Null;
  m_prev = Vector3Null;
  m_mass = 0.;
  m_automass = false;
  m_direction = Vector3Null;
  m_translationalVelocity = Vector3Null;
  m_SinCyclic_period = Vector3Null;
  m_SinCyclic_phase_shift = Vector3Null;
  m_prev = Vector3Null;
  m_vmaxzeroforce = 0.;
}




// ----------------------------------------------------------------------------
// Constructor with an XML node as input parameter
ObstacleImposedForce::ObstacleImposedForce( DOMNode* root, double dt, 
	int rank, size_t& error )
{
  m_ObstacleName = "Undefined";
  m_type = "Undefined";
  m_tstart = 0.; 
  m_tend = 0.;
  m_force_amplitude = Vector3Null;
  m_force = Vector3Null;
  m_prev = Vector3Null;
  m_mass = 0.;
  m_automass = false;  
  m_direction = Vector3Null;
  m_translationalVelocity = Vector3Null;
  m_SinCyclic_period = Vector3Null;
  m_SinCyclic_phase_shift = Vector3Null;
  m_prev = Vector3Null;
  m_vmaxzeroforce = 0.;  

  m_ObstacleName = ReaderXML::getNodeAttr_String( root, "ObstacleName" );
  
  DOMNode* nTimeInterval = ReaderXML::getNode( root, "TimeInterval" );
  m_tstart = ReaderXML::getNodeAttr_Double( nTimeInterval, "Start" );
  m_tend = ReaderXML::getNodeAttr_Double( nTimeInterval, "End" );
  
  // Constant translation
  if ( ReaderXML::getNode( root, "ConstantTranslation" ) )
  {  
    m_type = "ConstantTranslation";
    DOMNode* nCT = ReaderXML::getNode( root, "ConstantTranslation" );      
    DOMNode* force = ReaderXML::getNode( nCT, "Amplitude" );
    m_force_amplitude[X] = ReaderXML::getNodeAttr_Double( force, "AX" );
    m_force_amplitude[Y] = ReaderXML::getNodeAttr_Double( force, "AY" );
    m_force_amplitude[Z] = ReaderXML::getNodeAttr_Double( force, "AZ" ); 
    m_direction = m_force_amplitude / Norm( m_force_amplitude ); 
    m_force = m_force_amplitude; 
    DOMNode* property = ReaderXML::getNode( nCT, "Property" );
    if ( ReaderXML::hasNodeAttr( property, "Mass" ) )
      m_mass = ReaderXML::getNodeAttr_Double( property, "Mass" );
    else m_automass = true;
    m_vmaxzeroforce = ReaderXML::getNodeAttr_Double( property, "Vmax" );    
    if ( rank == 0 )
    {
      cout << GrainsExec::m_shift12 << "Obstacle name = " << m_ObstacleName 
      	<< endl;      
      cout << GrainsExec::m_shift12 << "Time interval = [" 
      	<< m_tstart << "," << m_tend << "]" << endl;
      cout << GrainsExec::m_shift12 << "Type = " << m_type << endl;
      cout << GrainsExec::m_shift12 << "Force = " << m_force[X] << " " << 
      	m_force[Y] << " " << m_force[Z] << endl;
      cout << GrainsExec::m_shift12 << "Mass = ";
      if ( m_automass ) cout << "auto" << endl;
      else cout << m_mass << endl;
      cout << GrainsExec::m_shift12 << "Vmax = " << m_vmaxzeroforce << endl;
    }
  } 
  else if ( ReaderXML::getNode( root, "SinCyclicTranslation" ) )
  {
    m_type = "SinCyclicTranslation";
    DOMNode* nCyclic = ReaderXML::getNode( root, "SinCyclicTranslation" );        
    DOMNode* force = ReaderXML::getNode( nCyclic, "Amplitude" );
    m_force_amplitude[X] = ReaderXML::getNodeAttr_Double( force, "AX" );
    m_force_amplitude[Y] = ReaderXML::getNodeAttr_Double( force, "AY" );
    m_force_amplitude[Z] = ReaderXML::getNodeAttr_Double( force, "AZ" ); 
    m_direction = m_force_amplitude / Norm( m_force_amplitude );   
    DOMNode* property = ReaderXML::getNode( nCyclic, "Property" );
    m_mass = ReaderXML::getNodeAttr_Double( property, "Mass" );
    DOMNode* nPer = ReaderXML::getNode( nCyclic, "Period" );
    m_SinCyclic_period[X] = ReaderXML::getNodeAttr_Double( nPer, "PX" );
    m_SinCyclic_period[Y] = ReaderXML::getNodeAttr_Double( nPer, "PY" );
    m_SinCyclic_period[Z] = ReaderXML::getNodeAttr_Double( nPer, "PZ" ); 
    DOMNode* nPhaseShift = ReaderXML::getNode( nCyclic, "PhaseShift" );    
    m_SinCyclic_phase_shift[X] = ReaderXML::getNodeAttr_Double( nPhaseShift, 
    	"PhiX" ) * PI / 180.;
    m_SinCyclic_phase_shift[Y] = ReaderXML::getNodeAttr_Double( nPhaseShift, 
    	"PhiY" ) * PI / 180.;	
    m_SinCyclic_phase_shift[Z] = ReaderXML::getNodeAttr_Double( nPhaseShift, 
    	"PhiZ" ) * PI / 180.;    
    if ( rank == 0 )
    {
      cout << GrainsExec::m_shift12 << "Obstacle name = " << m_ObstacleName 
      	<< endl;      
      cout << GrainsExec::m_shift12 << "Time interval = [" 
      	<< m_tstart << "," << m_tend << "]" << endl;
      cout << GrainsExec::m_shift12 << "Type = " << m_type << endl;
      cout << "Force = " << m_force[X] << " " << m_force[Y] << " " << 
    	m_force[Z] << endl; 
      cout << GrainsExec::m_shift12 << "Period = " << 
      	m_SinCyclic_period[X] << " " << m_SinCyclic_period[Y] << " " <<
	m_SinCyclic_period[Z] << endl;
      cout << GrainsExec::m_shift12 << "Phase shift in rad = " << 
      	m_SinCyclic_phase_shift[X] << " " << m_SinCyclic_phase_shift[Y] << " " 
	<< m_SinCyclic_phase_shift[Z] << endl;	
    }    
  }
}




// ----------------------------------------------------------------------------
// Destructor
ObstacleImposedForce::~ObstacleImposedForce()
{}




// ----------------------------------------------------------------------------
// Returns obstacle name
string ObstacleImposedForce::getObstacleName() const
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
// Returns the imposed force
Vector3 ObstacleImposedForce::getForce( double time )
{
  if ( m_type == "SinCyclicTranslation" ) SinCyclicForce( time );
  return ( m_force );
}




// ----------------------------------------------------------------------------
// Returns the obstacle virtual mass
double ObstacleImposedForce::getMass() const
{
  return ( m_mass );
}




// ----------------------------------------------------------------------------
// Returns the direction of motion
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
  Vector3 center = *(obstacle->getPosition());
  Vector3 force = *(obstacle->getForce());

  // Mass coeficient
  // If automatically determined, it is set such that the obstacle motion
  // over [t,t+dt] is equal to its crust thickness divided by 10^7 when
  // the controller is simply a proportional controller and df is equal 
  // to the imposed force  
  if ( m_automass )
    m_mass = dt * dt * Norm( m_force_amplitude ) * 1.e7 /
    	( 2. * obstacle->getCrustThickness() );

  // PID controller coefficients (Ziegler-Nichols method) 
  // The proportionality coefficient Kp is the solution of the simple
  // ODE m_mass * d^22 x/dt^2 = df with df constant over [t,t+dt] where
  // x is the obstacle motion, i.e. x = dt^2 * df / ( 2 * m_mass )  
  double Kp = 0.5 * dt * dt / m_mass ;
  double Ki = 2. * Kp / dt ;
  double Kd = 3. * Kp * dt / 24.; 
     
  // Translational velocity of the obstacle over [t,t+dt]
  Vector3 dforce, depl, trans;
  static Vector3 dforce_nm1;
  static Vector3 dforce_nm2;
  static Vector3 depl_nm1;
  
  static bool nonzeroforce = false;
  if ( nonzeroforce == false ) nonzeroforce = Norm( force ) > EPSILON;

  if ( m_type == "ConstantTranslation" )
  {
    m_force = m_force_amplitude;
    dforce = force - m_force;
    for (size_t i=0;i<3;++i) dforce[i] *= m_direction[i]; 
    if ( nonzeroforce )
    { 
      depl = depl_nm1 + ( Kp + Ki * dt + Kd / dt ) * dforce
      	- ( Kp + 2. * Kd / dt ) * dforce_nm1 + ( Kd / dt ) * dforce_nm2; 
      m_translationalVelocity = depl / dt;
      dforce_nm2 = dforce_nm1;
      dforce_nm1 = dforce;
      depl_nm1 = depl;
    }
    else
      m_translationalVelocity = - m_vmaxzeroforce * m_direction;           
  }
  else if ( m_type == "SinCyclicTranslation" )
  {
    // TO DO
    
//     SinCyclicForce( time );
//     dforce = force - m_force;
//     for (size_t i=0;i<3;++i) dforce[i] *= m_direction[i]; 
//     trans = 0.5 *( dt * dt / m_mass ) * dforce; 
//     depl = trans - m_prev;
//     m_translationalVelocity = depl / dt;
//     // t-dt 
//     m_prev = trans;
  }
  
  return ( &m_translationalVelocity );
}




// ----------------------------------------------------------------------------
// Sets the sinusoidal cyclic force at a given time 
void ObstacleImposedForce::SinCyclicForce( double time )
{
  for (size_t i=0;i<3;++i)
    m_force[i] = m_force_amplitude[i] * sin( 2. * PI * ( time - m_tstart ) 
    	/ m_SinCyclic_period[i] );
}
