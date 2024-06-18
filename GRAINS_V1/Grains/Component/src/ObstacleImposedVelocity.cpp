#include "ObstacleImposedVelocity.hh"
#include "Obstacle.hh"
#include "GrainsExec.hh"
#include <stdlib.h>


// ----------------------------------------------------------------------------
// Default constructor
ObstacleImposedVelocity::ObstacleImposedVelocity()
{
  m_ObstacleName = "Undefined";
  m_type = "Undefined";
  m_tstart = 0.; 
  m_tend = 0.;
  m_translationalVelocity = Vector3Null;
  m_angularVelocity = Vector3Null;
  m_rotationCenterIsCenterOfMass = true;
  m_rotationCenter = Point3Null;
  m_Sin_amplitude = 0.;
  m_Sin_period = 0.;
  m_Sin_phase_shift = 0.;  
  m_unit_vitRef = Vector3Null;
  m_SinCyclic_period = Vector3Null;
  m_SinCyclic_amplitude = Vector3Null;
  m_SinCyclic_phase_shift = Vector3Null;  
}




// ----------------------------------------------------------------------------
// Constructor with an XML node as input parameter
ObstacleImposedVelocity::ObstacleImposedVelocity( DOMNode* root, 
	double dt, int rank, size_t& error )
{
  m_type = "Undefined";
  m_translationalVelocity = Vector3Null;
  m_angularVelocity = Vector3Null;
  m_rotationCenterIsCenterOfMass = true;
  m_rotationCenter = Point3Null;
  m_Sin_amplitude = 0.;
  m_Sin_period = 0.;
  m_Sin_phase_shift = 0.;  
  m_unit_vitRef = Vector3Null;
  m_SinCyclic_period = Vector3Null;
  m_SinCyclic_amplitude = Vector3Null;
  m_SinCyclic_phase_shift = Vector3Null;
  
  error = 0;   
    
  m_ObstacleName = ReaderXML::getNodeAttr_String( root, "ObstacleName" );
  
  DOMNode* nTimeInterval = ReaderXML::getNode( root, "TimeInterval" );
  m_tstart = ReaderXML::getNodeAttr_Double( nTimeInterval, "Start" );
  m_tend = ReaderXML::getNodeAttr_Double( nTimeInterval, "End" );

  // Constant translation
  if ( ReaderXML::getNode( root, "ConstantTranslation" ) )
  {
    m_type = "ConstantTranslation";
    DOMNode* nVTranslation = ReaderXML::getNode( root, "ConstantTranslation" );
    DOMNode* nTV = ReaderXML::getNode( nVTranslation, "TranslationalVelocity" );
    m_translationalVelocity[X] = ReaderXML::getNodeAttr_Double( nTV, "VX" );
    m_translationalVelocity[Y] = ReaderXML::getNodeAttr_Double( nTV, "VY" );    
    m_translationalVelocity[Z] = ReaderXML::getNodeAttr_Double( nTV, "VZ" ); 
    if ( rank == 0 )
    {
      cout << GrainsExec::m_shift12 << "Obstacle name = " << m_ObstacleName 
      	<< endl;      
      cout << GrainsExec::m_shift12 << "Time interval = [" 
      	<< m_tstart << "," << m_tend << "]" << endl;
      cout << GrainsExec::m_shift12 << "Type = " << m_type << endl;
      cout << GrainsExec::m_shift12 << "Translational velocity = " << 
      	m_translationalVelocity << endl;
    }   
  }
  // Sinusoidal translation  
  else if ( ReaderXML::getNode( root, "SinTranslation" ) )
  {
    m_type = "SinTranslation";
    DOMNode* nVTranslation = ReaderXML::getNode( root, "SinTranslation" ); 
    DOMNode* nTV = ReaderXML::getNode( nVTranslation, "Direction" );
    m_unit_vitRef[X] = ReaderXML::getNodeAttr_Double( nTV, "VX" );
    m_unit_vitRef[Y] = ReaderXML::getNodeAttr_Double( nTV, "VY" );    
    m_unit_vitRef[Z] = ReaderXML::getNodeAttr_Double( nTV, "VZ" ); 
    m_unit_vitRef.normalize();
    DOMNode* nPar = ReaderXML::getNode( nVTranslation, "Parameters" );        
    m_Sin_amplitude = ReaderXML::getNodeAttr_Double( nPar, "Amplitude" ); 
    m_Sin_period = ReaderXML::getNodeAttr_Double( nPar, "Period" ); 
    m_Sin_phase_shift = ReaderXML::getNodeAttr_Double( nPar, "PhaseShift" ) 
    	* PI / 180.;
    if ( rank == 0 )
    {
      cout << GrainsExec::m_shift12 << "Obstacle name = " << m_ObstacleName 
      	<< endl;      
      cout << GrainsExec::m_shift12 << "Time interval = [" 
      	<< m_tstart << "," << m_tend << "]" << endl;
      cout << GrainsExec::m_shift12 << "Type = " << m_type << endl;
      cout << GrainsExec::m_shift12 << "Unit direction vector = " << 
      	m_unit_vitRef << endl;
      cout << GrainsExec::m_shift12 << "Amplitude = " << 
      	m_Sin_amplitude << endl;
      cout << GrainsExec::m_shift12 << "Period = " << 
      	m_Sin_period << endl;
      cout << GrainsExec::m_shift12 << "Phase shift in rad = " << 
      	m_Sin_phase_shift << endl;	
      cout << GrainsExec::m_shift12 << "Maximum motion = " << 
      	m_Sin_amplitude * m_Sin_period / ( 2. * PI ) << endl;
      cout << GrainsExec::m_shift12 << "Maximum acceleration = " << 
      	m_Sin_amplitude * 2. * PI / m_Sin_period << endl; 
    }             
  }
  // Constant rotation
  else if ( ReaderXML::getNode( root, "ConstantRotation" ) )
  {
    m_type = "ConstantRotation";
    DOMNode* nVRotation = ReaderXML::getNode( root, "ConstantRotation" );
    DOMNode* nRV = ReaderXML::getNode( nVRotation, "AngularVelocity" );
    m_angularVelocity[X] = ReaderXML::getNodeAttr_Double( nRV, "WX" );
    m_angularVelocity[Y] = ReaderXML::getNodeAttr_Double( nRV, "WY" );    
    m_angularVelocity[Z] = ReaderXML::getNodeAttr_Double( nRV, "WZ" );
    DOMNode* nRC = ReaderXML::getNode( nVRotation, "RotationCenter" );
    if ( nRC )
    {
      m_rotationCenterIsCenterOfMass = false;
      m_rotationCenter[X] = ReaderXML::getNodeAttr_Double( nRC, "CX" );
      m_rotationCenter[Y] = ReaderXML::getNodeAttr_Double( nRC, "CY" );    
      m_rotationCenter[Z] = ReaderXML::getNodeAttr_Double( nRC, "CZ" );      
    }      
    if ( rank == 0 )
    {
      cout << GrainsExec::m_shift12 << "Obstacle name = " << m_ObstacleName 
      	<< endl;      
      cout << GrainsExec::m_shift12 << "Time interval = [" 
      	<< m_tstart << "," << m_tend << "]" << endl;
      cout << GrainsExec::m_shift12 << "Type = " << m_type << endl;
      cout << GrainsExec::m_shift12 << "Angular velocity = " << 
      	m_angularVelocity << endl;
      cout << GrainsExec::m_shift12 << "Center of rotation = ";
      if ( m_rotationCenterIsCenterOfMass ) cout << "center of mass";
      else cout << m_rotationCenter[X] << " " << m_rotationCenter[Y] <<
	" " <<  m_rotationCenter[Z]; 
      cout << endl;	
    }   
  }
  // Cyclic translation
  else if ( ReaderXML::getNode( root, "SinCyclicTranslation" ) )
  {
    m_type = "SinCyclicTranslation";
    DOMNode* nCyclic = ReaderXML::getNode( root, "SinCyclicTranslation" );    
    DOMNode* nTV = ReaderXML::getNode( nCyclic, "Direction" );
    m_unit_vitRef[X] = ReaderXML::getNodeAttr_Double( nTV, "VX" );
    m_unit_vitRef[Y] = ReaderXML::getNodeAttr_Double( nTV, "VY" );    
    m_unit_vitRef[Z] = ReaderXML::getNodeAttr_Double( nTV, "VZ" ); 
    m_unit_vitRef.normalize();    
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
    DOMNode* nAmp = ReaderXML::getNode( nCyclic, "Amplitude" );
    m_SinCyclic_amplitude[X] = ReaderXML::getNodeAttr_Double( nAmp, "AX" ) 
    	* m_unit_vitRef[X];
    m_SinCyclic_amplitude[Y] = ReaderXML::getNodeAttr_Double( nAmp, "AY" ) 
    	* m_unit_vitRef[Y];
    m_SinCyclic_amplitude[Z] = ReaderXML::getNodeAttr_Double( nAmp, "AZ" ) 
    	* m_unit_vitRef[Z];
    if ( rank == 0 )
    {
      cout << GrainsExec::m_shift12 << "Obstacle name = " << m_ObstacleName 
      	<< endl;      
      cout << GrainsExec::m_shift12 << "Time interval = [" 
      	<< m_tstart << "," << m_tend << "]" << endl;
      cout << GrainsExec::m_shift12 << "Type = " << m_type << endl;
      cout << GrainsExec::m_shift12 << "Amplitude = " << 
      	 m_SinCyclic_amplitude[X] << " " << m_SinCyclic_amplitude[Y] << " " << 
	 m_SinCyclic_amplitude[Z] << endl;
      cout << GrainsExec::m_shift12 << "Period = " << 
      	m_SinCyclic_period[X] << " " << m_SinCyclic_period[Y] << " " <<
	m_SinCyclic_period[Z] << endl;
      cout << GrainsExec::m_shift12 << "Phase shift in rad = " << 
      	m_SinCyclic_phase_shift[X] << " " << m_SinCyclic_phase_shift[Y] << " " 
	<< m_SinCyclic_phase_shift[Z] << endl;	 
    }   
  }    
  else 
  {
    error = 1;
    if ( rank == 0 ) cout << GrainsExec::m_shift6 << 
	"Unknown or missing obstacle imposed velocity node !!" << endl;
  }
}




// ----------------------------------------------------------------------------
// Destructor
ObstacleImposedVelocity::~ObstacleImposedVelocity()
{}




// ----------------------------------------------------------------------------
// Returns obstacle name
string ObstacleImposedVelocity::getObstacleName() const
{
  return ( m_ObstacleName );
}




// ----------------------------------------------------------------------------
// Returns the remaining active time interval of the imposed motion
double ObstacleImposedVelocity::getTime( double debut, double fin ) const
{
  double activtimeint = fin - debut;

  if ( debut < m_tstart ) activtimeint -= ( m_tstart - debut );
  if ( m_tend < fin ) activtimeint -= ( fin - m_tend );

  return ( activtimeint );
}




// ----------------------------------------------------------------------------
// Returns whether the imposed motion is activ at time t
bool ObstacleImposedVelocity::isActif( double t, double dt ) const 
{
  return ( t > m_tstart - dt * 1.e-5  && t < m_tend + dt * 1.e-5 );
}




// ----------------------------------------------------------------------------
// Returns whether the imposed motion is completed at time t
bool ObstacleImposedVelocity::isCompleted( double t, double dt ) const 
{
  return ( t > m_tend + dt * 1.e-5 );
}




// ----------------------------------------------------------------------------
// Returns the translational velocity at time t 
Vector3 const* ObstacleImposedVelocity::translationalVelocity( double time, 
	double dt, Point3 const& cg )
{
  if ( m_type == "SinTranslation" )
    m_translationalVelocity = m_Sin_amplitude * 
    	sin( 2. * PI * ( time - m_tstart ) / m_Sin_period + m_Sin_phase_shift )
	* m_unit_vitRef ;
  else if ( m_type == "SinCyclicTranslation" )
  {
    for (size_t i=0;i<3;++i)
      m_translationalVelocity[i] = m_SinCyclic_amplitude[i] * 
    	sin( 2. * PI * ( time - m_tstart ) / m_SinCyclic_period[i] 
	+ m_SinCyclic_phase_shift[i] ) * m_unit_vitRef[i] ;	    	
  }
  else if ( m_type == "ConstantRotation" && !m_rotationCenterIsCenterOfMass )
    m_translationalVelocity = m_angularVelocity ^ ( cg - m_rotationCenter );

  return ( &m_translationalVelocity );
}




// ----------------------------------------------------------------------------
// Returns the angular velocity at time t 
Vector3 const* ObstacleImposedVelocity::angularVelocity( double time, 
	double dt )
{    
  return ( &m_angularVelocity );
}




// ----------------------------------------------------------------------------
// Returns the translational motion over dt at time t
Vector3 ObstacleImposedVelocity::translationalMotion( double time, 
	double dt, Point3 const& cg ) 
{
  // We integrate the velocity over [t,t+dt] analytically to be more accurate
  Vector3 translation;
  
  if ( m_type == "ConstantTranslation" ) 
    translation = m_translationalVelocity * dt;
  else if ( m_type == "SinTranslation" )
    translation = ( m_Sin_amplitude * m_Sin_period / ( 2. * PI ) )
    	* ( cos( 2. * PI * ( time - m_tstart ) / m_Sin_period 
		+ m_Sin_phase_shift ) 
	- cos( 2. * PI * ( time + dt - m_tstart ) / m_Sin_period 
		+ m_Sin_phase_shift ) ) * m_unit_vitRef ;
  else if ( m_type == "SinCyclicTranslation" )
  {
    for (size_t i=0;i<3;++i)
      translation[i] = 
      	( m_SinCyclic_amplitude[i] * m_SinCyclic_period[i] / ( 2. * PI ) )
    	* ( cos( 2. * PI * ( time - m_tstart ) / m_SinCyclic_period[i] 
		+ m_SinCyclic_phase_shift[i] ) 
	- cos( 2. * PI * ( time + dt - m_tstart ) / m_SinCyclic_period[i] 
		+ m_SinCyclic_phase_shift[i] ) ) * m_unit_vitRef[i] ;	    	
  }
  else if ( m_type == "ConstantRotation" && !m_rotationCenterIsCenterOfMass )
  {
    Vector3 rota = angularMotion( time, dt );
    double d = Norm(rota);
    Quaternion q;

    if ( d != 0. ) 
    {
      Vector3 vect = ( sin( d / 2. ) / d ) * rota;
      q = Quaternion( vect, cos( d / 2. ) );
    }
    else 
      q = Quaternion( 0., 0., 0., 1. );
    translation = q.rotateVector( cg - m_rotationCenter ) 
    	- ( cg - m_rotationCenter );    
  }
    
  return ( translation );
}  




// ----------------------------------------------------------------------------
// Returns the angular motion over dt at time t 
Vector3 ObstacleImposedVelocity::angularMotion( double time, double dt )
{
  // We only consider contant angular velocity so far
  return ( m_angularVelocity * dt );
}  




// ----------------------------------------------------------------------------
// Operator == based on the object address
bool ObstacleImposedVelocity::operator == ( 
	ObstacleImposedVelocity const& other ) const
{
  return ( this == &other );
}




// ----------------------------------------------------------------------------
// Operator < based on the start time of the imposed motion.
// Returns true if c0.tdebut < c1.tdebut
bool operator < ( ObstacleImposedVelocity const& c0,
	ObstacleImposedVelocity const& c1 )
{
  return ( c0.m_tstart < c1.m_tstart );
}




// ----------------------------------------------------------------------------
// Debug
void ObstacleImposedVelocity::debug( char* c )
{
  cout << m_tstart << '-' << m_tend << '\t';
}




// ----------------------------------------------------------------------------
// Returns the imposed motion type
string ObstacleImposedVelocity::getType() const
{
  return ( m_type );
}
 
