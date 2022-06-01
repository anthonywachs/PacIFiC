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
  m_translationalVelocity = Vector3Nul;
  m_angularVelocity = Vector3Nul;
  m_Sin_amplitude = 0.;
  m_Sin_period = 0.;
  m_Sin_vitRef = Vector3Nul;
  m_freqX = 0.;
  m_freqY = 0.;
  m_freqZ = 0.;
  m_phase = 0.;
  m_ampX = 0.;
  m_ampY = 0.;  
  m_ampZ = 0.;  
  m_prev = Vector3Nul;
}




// ----------------------------------------------------------------------------
// Constructor with an XML node as input parameter
ObstacleImposedVelocity::ObstacleImposedVelocity( DOMNode* root, 
	double dt, int rank, size_t& error )
{
  error = 0;
  m_ObstacleName = ReaderXML::getNodeAttr_String( root, "ObstacleName" );
  
  DOMNode* nTimeInterval = ReaderXML::getNode( root, "TimeInterval" );
  m_tstart = ReaderXML::getNodeAttr_Double( nTimeInterval, "Start" );
  m_tend = ReaderXML::getNodeAttr_Double( nTimeInterval, "End" );

  // Constant translation
  if ( ReaderXML::getNode( root, "VTranslation" ) )
  {
    DOMNode* nVTranslation = ReaderXML::getNode( root, "VTranslation" );
    DOMNode* nTV = ReaderXML::getNode( nVTranslation, "TranslationalVelocity" );
    m_translationalVelocity[X] = ReaderXML::getNodeAttr_Double( nTV, "VX" );
    m_translationalVelocity[Y] = ReaderXML::getNodeAttr_Double( nTV, "VY" );    
    m_translationalVelocity[Z] = ReaderXML::getNodeAttr_Double( nTV, "VZ" ); 
    m_type = "ConstantTranslation";
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
  else
  {
    error = 1;
    if ( rank == 0 ) cout << GrainsExec::m_shift6 << 
	"Unknown or missing obstacle imposed velocity node !!" << endl;
  }

// 
// 
//   DOMNode* nVector3 = ReaderXML::getNode( root, "Vector3" );
//   Vector3 vecteur;
//   vecteur[X] = ReaderXML::getNodeAttr_Double( nVector3, "X" );
//   vecteur[Y] = ReaderXML::getNodeAttr_Double( nVector3, "Y" );
//   vecteur[Z] = ReaderXML::getNodeAttr_Double( nVector3, "Z" );
// 
//   double deltaT = m_tend - m_tstart;
// 
//   string mode = ReaderXML::getNodeAttr_String( root, "Mode" );
// 
//   bool istype = ReaderXML::hasNodeAttr( root, "Type" ); 
// 
//   if ( mode == "Translation" ) 
//   {
//     if ( !istype )
//     {
//       m_translationalVelocity = vecteur / deltaT;
//       m_type = mode;
//   
//       cout << endl << "Chargement translationnel sur " << m_ObstacleName << endl;
//       cout << "   Vitese constante de translation = " 
//       	<< m_translationalVelocity[X] << " " << m_translationalVelocity[Y] << " " 
//   	<< m_translationalVelocity[Z] << endl;
//     }
//     else
//     {
//       m_type = "Cyclic";
//       DOMNode* freq = ReaderXML::getNode( root, "Frequence" );
//       DOMNode* amp = ReaderXML::getNode( root, "Amplitude" );
//       m_freqX = ReaderXML::getNodeAttr_Double( freq, "FX" );
//       m_freqY = ReaderXML::getNodeAttr_Double( freq, "FY" );
//       m_freqZ = ReaderXML::getNodeAttr_Double( freq, "FZ" );
//       m_phase = ReaderXML::getNodeAttr_Double( freq, "Phi" ) * PI / 180.;
//       m_ampX = ReaderXML::getNodeAttr_Double( amp, "AX" ) * vecteur[X];
//       m_ampY = ReaderXML::getNodeAttr_Double( amp, "AY" ) * vecteur[Y];
//       m_ampZ = ReaderXML::getNodeAttr_Double( amp, "AZ" ) * vecteur[Z];
//       cout << "Chargement translationnel cyclic sur " << m_ObstacleName << endl;
//       cout << "   Amplitude de translation = " 
//       	<< m_ampX << "\t" << m_ampY << "\t" << m_ampZ << endl;
//       cout << "   Frequence : FX = " << m_freqX << "\tFY = " << m_freqY 
// 	   << "\tFZ = " << m_freqZ << endl;
//     }
//   }
//   else if ( mode == "Rotation" ) 
//   {
//     m_angularVelocity = vecteur / deltaT;
//     m_type = mode;
//     
//     cout << endl << "Chargement rotationnel sur " << m_ObstacleName << endl;
//     cout << "   Vitese constante de rotation = " 
//     	<< m_angularVelocity[X] << " " << m_angularVelocity[Y] << " " 
// 	<< m_angularVelocity[Z] << endl;    
//   }
//   else if ( mode == "RotationSinusoidale" ) 
//   {
//     m_Sin_vitRef = vecteur / Norm(vecteur);
//     m_angularVelocity.reset();    
//     m_type = mode;
//     m_Sin_amplitude = ReaderXML::getNodeAttr_Double( root, "A" );
//     m_Sin_period = ReaderXML::getNodeAttr_Double( root, "P" ); 
//     
//     cout << endl << "Chargement rotationnel sinusoidal sur " 
//     	<< m_ObstacleName << endl;
//     cout << "   Periode du mouvement = " << m_Sin_period << endl;
//     cout << "   Amplitude angulaire max du mouvement = " << m_Sin_amplitude * 
//     	m_Sin_period / ( 2. * PI ) << endl;
//     cout << "   Acceleration angulaire max du mouvement = " << m_Sin_amplitude *
//     	2. * PI / m_Sin_period << endl;       
//   } 
//   else if ( mode == "TranslationSinusoidale" ) 
//   {
//     m_Sin_vitRef = vecteur / Norm(vecteur);
//     m_translationalVelocity.reset();
//     m_type = mode;
//     m_Sin_amplitude = ReaderXML::getNodeAttr_Double( root, "A" );
//     m_Sin_period = ReaderXML::getNodeAttr_Double( root, "P" );
//     cout << endl << "Chargement translationnel sinusoidal sur " 
//     	<< m_ObstacleName << endl;
//     cout << "   Periode du mouvement = " << m_Sin_period << endl;
//     cout << "   Amplitude max du mouvement = " << m_Sin_amplitude * 
//     	m_Sin_period / ( 2. * PI ) << endl;
//     cout << "   Acceleration max du mouvement = " << m_Sin_amplitude *
//     	2. * PI / m_Sin_period << endl; 
//   }
// 
//   if ( rank == 0 )
//   {
//     cout << "   Temps de depart = " << m_tstart << endl;
//     cout << "   Temps de fin = " << m_tend << endl << endl;
//   }     
}




// ----------------------------------------------------------------------------
// Destructor
ObstacleImposedVelocity::~ObstacleImposedVelocity()
{}




// ----------------------------------------------------------------------------
// Returns obstacle name
string ObstacleImposedVelocity::getNom() const
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
	double dt )
{
  if ( m_type == "TranslationSinusoidale" )
    m_translationalVelocity = m_Sin_amplitude * 
    	sin( 2. * PI * ( time - m_tstart ) / m_Sin_period )
	* m_Sin_vitRef ;
  else if ( m_type == "Cyclic" )
  {
    Vector3 trans, dx;
    trans[X] = m_ampX * sin( 2. * PI * m_freqX 
	* ( time - m_tstart ) ); 
    trans[Y] = m_ampY * sin( 2. * PI * m_freqY 
	* ( time - m_tstart ) + m_phase ); 
    trans[Z] = m_ampZ * sin( 2. * PI * m_freqZ 
	* ( time - m_tstart ) + m_phase ); 
    dx = trans - m_prev; 
    m_prev = trans;
    m_translationalVelocity = dx / dt;
  }

  return ( &m_translationalVelocity );
}




// ----------------------------------------------------------------------------
// Returns the angular velocity at time t 
Vector3 const* ObstacleImposedVelocity::angularVelocity( double time, 
	double dt )
{
  if ( m_type == "RotationSinusoidale" )
    m_angularVelocity = m_Sin_amplitude * 
    	sin( 2. * PI * ( time - m_tstart ) / m_Sin_period )
	* m_Sin_vitRef ;
     
  return ( &m_angularVelocity );
}




// ----------------------------------------------------------------------------
// Returns the translational displacement over dt at time t
Vector3 ObstacleImposedVelocity::translationalDisplacement( double time, 
	double dt ) 
{
  return ( m_translationalVelocity * dt );
}  




// ----------------------------------------------------------------------------
// Returns the angular displacement over dt at time t 
Vector3 ObstacleImposedVelocity::angularDisplacement( double time, double dt )
{
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
// Output operator
ostream &operator << ( ostream& fileOut, 
	ObstacleImposedVelocity const& motion )
{
  fileOut << motion.m_ObstacleName << '\n';
  fileOut << motion.m_tstart << '\t' << motion.m_tend << '\n';
  fileOut << motion.m_translationalVelocity;
  fileOut << motion.m_angularVelocity;

  return ( fileOut );
}




// ----------------------------------------------------------------------------
// Input operator
istream &operator >> ( istream& fileIn, 
	ObstacleImposedVelocity& motion )
{
  fileIn >> motion.m_ObstacleName;
  fileIn >> motion.m_tstart >> motion.m_tend;
  fileIn >> motion.m_translationalVelocity;
  if ( motion.m_translationalVelocity != Vector3Nul ) 
    motion.m_type = "Translation";
  fileIn >> motion.m_angularVelocity;
  if ( motion.m_angularVelocity != Vector3Nul ) 
    motion.m_type = "Rotation";  

  return ( fileIn );
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
 
