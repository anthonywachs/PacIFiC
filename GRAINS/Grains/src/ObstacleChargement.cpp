#include "ObstacleChargement.H"
#include "Obstacle.H"
#include <stdlib.h>



Vecteur ObstacleChargement::m_prev = VecteurNul;
// ----------------------------------------------------------------------------
// G.FERRER - Dece.1999 - Creation
// Constructeur par defaut
ObstacleChargement::ObstacleChargement()
{
  m_nomObstacle = "Undefined";
  m_type = "Undefined";
  m_tdebut = 0.; 
  m_tfin = 0.;
  m_Sin_amplitude = 0.;
  m_Sin_periode = 0.;
}




// ----------------------------------------------------------------------------
// Constructeur par copie
// G.FERRER - Dece.2001 - Creation
ObstacleChargement::ObstacleChargement( const ObstacleChargement &copie ) :
  m_nomObstacle( copie.m_nomObstacle ), 
  m_type( copie.m_type ),
  m_tdebut( copie.m_tdebut ), 
  m_tfin( copie.m_tfin ),
  m_vitesse_translation( copie.m_vitesse_translation ), 
  m_vitesse_rotation( copie.m_vitesse_rotation ),
  m_Sin_amplitude( copie.m_Sin_amplitude ),
  m_Sin_periode( copie.m_Sin_periode ),
  m_Sin_vitRef( copie.m_Sin_vitRef )  
{}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Creation d'un chargement d'Obstacle
ObstacleChargement::ObstacleChargement( DOMNode* root, double dt, int rank )
{
  // Remarque: pour l'instant, valable pour des translations et des
  // rotations à vitesse constante uniquement 

  m_nomObstacle = ReaderXML::getNodeAttr_String( root, "NomObstacle" );
  
  DOMNode* temps = ReaderXML::getNode( root, "Temps" );
  m_tdebut      = ReaderXML::getNodeAttr_Double( temps, "Debut" );
  m_tfin        = ReaderXML::getNodeAttr_Double( temps, "Fin" );

  DOMNode* nVecteur = ReaderXML::getNode( root, "Vecteur" );
  Vecteur vecteur;
  vecteur[X] = ReaderXML::getNodeAttr_Double( nVecteur, "X" );
  vecteur[Y] = ReaderXML::getNodeAttr_Double( nVecteur, "Y" );
  vecteur[Z] = ReaderXML::getNodeAttr_Double( nVecteur, "Z" );

  double deltaT = m_tfin - m_tdebut;

  string mode = ReaderXML::getNodeAttr_String( root, "Mode" );

  bool istype = ReaderXML::hasNodeAttr_String( root, "Type" ); 

  if ( mode == "Translation" ) 
  {
    if ( !istype )
    {
      m_vitesse_translation = vecteur / deltaT;
      m_type = mode;
  
      cout << endl << "Chargement translationnel sur " << m_nomObstacle << endl;
      cout << "   Vitese constante de translation = " 
      	<< m_vitesse_translation[X] << " " << m_vitesse_translation[Y] << " " 
  	<< m_vitesse_translation[Z] << endl;
    }
    else
    {
      m_type = "Cyclic";
      DOMNode* freq = ReaderXML::getNode( root, "Frequence" );
      DOMNode* amp = ReaderXML::getNode( root, "Amplitude" );
      m_freqX = ReaderXML::getNodeAttr_Double( freq, "FX" );
      m_freqY = ReaderXML::getNodeAttr_Double( freq, "FY" );
      m_freqZ = ReaderXML::getNodeAttr_Double( freq, "FZ" );
      m_phase = ReaderXML::getNodeAttr_Double( freq, "Phi" ) * PI / 180.;
      m_ampX = ReaderXML::getNodeAttr_Double( amp, "AX" ) * vecteur[X];
      m_ampY = ReaderXML::getNodeAttr_Double( amp, "AY" ) * vecteur[Y];
      m_ampZ = ReaderXML::getNodeAttr_Double( amp, "AZ" ) * vecteur[Z];
      cout << "Chargement translationnel cyclic sur " << m_nomObstacle << endl;
      cout << "   Amplitude de translation = " 
      	<< m_ampX << "\t" << m_ampY << "\t" << m_ampZ << endl;
      cout << "   Frequence : FX = " << m_freqX << "\tFY = " << m_freqY 
	   << "\tFZ = " << m_freqZ << endl;
    }
  }
  else if ( mode == "Rotation" ) 
  {
    m_vitesse_rotation = vecteur / deltaT;
    m_type = mode;
    
    cout << endl << "Chargement rotationnel sur " << m_nomObstacle << endl;
    cout << "   Vitese constante de rotation = " 
    	<< m_vitesse_rotation[X] << " " << m_vitesse_rotation[Y] << " " 
	<< m_vitesse_rotation[Z] << endl;    
  }
  else if ( mode == "RotationSinusoidale" ) 
  {
    m_Sin_vitRef = vecteur / Norm(vecteur);
    m_vitesse_rotation.reset();    
    m_type = mode;
    m_Sin_amplitude = ReaderXML::getNodeAttr_Double( root, "A" );
    m_Sin_periode = ReaderXML::getNodeAttr_Double( root, "P" ); 
    
    cout << endl << "Chargement rotationnel sinusoidal sur " 
    	<< m_nomObstacle << endl;
    cout << "   Periode du mouvement = " << m_Sin_periode << endl;
    cout << "   Amplitude angulaire max du mouvement = " << m_Sin_amplitude * 
    	m_Sin_periode / ( 2. * PI ) << endl;
    cout << "   Acceleration angulaire max du mouvement = " << m_Sin_amplitude *
    	2. * PI / m_Sin_periode << endl;       
  } 
  else if ( mode == "TranslationSinusoidale" ) 
  {
    m_Sin_vitRef = vecteur / Norm(vecteur);
    m_vitesse_translation.reset();
    m_type = mode;
    m_Sin_amplitude = ReaderXML::getNodeAttr_Double( root, "A" );
    m_Sin_periode = ReaderXML::getNodeAttr_Double( root, "P" );
    cout << endl << "Chargement translationnel sinusoidal sur " 
    	<< m_nomObstacle << endl;
    cout << "   Periode du mouvement = " << m_Sin_periode << endl;
    cout << "   Amplitude max du mouvement = " << m_Sin_amplitude * 
    	m_Sin_periode / ( 2. * PI ) << endl;
    cout << "   Acceleration max du mouvement = " << m_Sin_amplitude *
    	2. * PI / m_Sin_periode << endl; 
  }

  if ( rank == 0 )
  {
    cout << "   Temps de depart = " << m_tdebut << endl;
    cout << "   Temps de fin = " << m_tfin << endl << endl;
  }     
}




// ----------------------------------------------------------------------------
// G.FERRER - Dece.1999 - Creation
// Destructeur
ObstacleChargement::~ObstacleChargement()
{
}




// ----------------------------------------------------------------------------
// G.FERRER - Dece.1999 - Creation
// Identificateur de la obstacle
string ObstacleChargement::getNom() const
{
  return m_nomObstacle;
}




// ----------------------------------------------------------------------------
// Temps de chargement dans l'espace indique
// G.FERRER - Juil.2003 - Creation
Scalar ObstacleChargement::getTime( Scalar debut, Scalar fin ) const
{
  Scalar dt = fin - debut;

  if ( debut < m_tdebut )
    dt -= ( m_tdebut - debut );
  if ( m_tfin < fin )
    dt -= ( fin - m_tfin );

  return dt;
}




// ----------------------------------------------------------------------------
// Le chargement est il actif ?
// G.FERRER - Mars.2000 - Creation
// G.FERRER - Juil.2003 - Doit renvoyer un bool
bool ObstacleChargement::isActif( Scalar t, Scalar dt ) const 
{
  return ( t > m_tdebut- dt * 1.e-5  && t < m_tfin + dt * 1.e-5 );
}




// ----------------------------------------------------------------------------
// Le chargement est il termine?
// A.WACHS - Mars.2012 - Creation
bool ObstacleChargement::isCompleted( Scalar t, Scalar dt ) const 
{
  return ( t > m_tfin + dt * 1.e-5 );
}




// ----------------------------------------------------------------------------
// Vitesse de deplacement
// G.FERRER - Juin.2000 - Creation
// D. RAKOTONIRINA - Mai. 2017 - Modification
Vecteur const* ObstacleChargement::VitesseTranslation( Scalar temps, 
	Scalar dt )
{
  if ( m_type == "TranslationSinusoidale" )
    m_vitesse_translation = m_Sin_amplitude * 
    	sin( 2. * PI * ( temps - m_tdebut ) / m_Sin_periode )
	* m_Sin_vitRef ;
  else if ( m_type == "Cyclic" )
  {
    Vecteur trans, dx;
    trans[X] = m_ampX * sin( 2. * PI * m_freqX 
	* (temps - m_tdebut) ); 
    trans[Y] = m_ampY * sin( 2. * PI * m_freqY 
	* (temps - m_tdebut) + m_phase ); 
    trans[Z] = m_ampZ * sin( 2. * PI * m_freqZ 
	* (temps - m_tdebut) + m_phase ); 
    dx = trans - m_prev; 
    // t-dt
    m_prev = trans;
    m_vitesse_translation = dx / dt;
  }

  return &m_vitesse_translation;
}




// ----------------------------------------------------------------------------
// Renvoie le signe d'un double
// D. RAKOTONIRINA - Mai. 2017 - Modification
Scalar ObstacleChargement::sign( Scalar v ) const
{
  Scalar ttt = (v > 0) - (v < 0);
  return ttt;
}




// ----------------------------------------------------------------------------
// Vitesse de rotation instantanee
// F.PRADEL - Nov.2001 - Creation
Vecteur const* ObstacleChargement::VitesseRotation( Scalar temps, 
	Scalar dt )
{
  if ( m_type == "RotationSinusoidale" )
    m_vitesse_rotation = m_Sin_amplitude * 
    	sin( 2. * PI * ( temps - m_tdebut ) / m_Sin_periode )
	* m_Sin_vitRef ;
     
  return &m_vitesse_rotation;
}




// ----------------------------------------------------------------------------
// Deplacement en translation sur un pas de temps à un temps donné
Vecteur ObstacleChargement::DeplacementTranslation( Scalar temps, Scalar dt ) 
{
  return m_vitesse_translation * dt;
}  




// ----------------------------------------------------------------------------
// Deplacement en rotation sur un pas de temps à un temps donné
Vecteur ObstacleChargement::DeplacementRotation( Scalar temps, Scalar dt )
{
  return m_vitesse_rotation * dt;
}  




// ----------------------------------------------------------------------------
// Decalage dans le temps du chargement
// F.PRADEL - Nov.2001 - Creation
void ObstacleChargement::TempsSlip( Scalar decalage )
{
  m_tdebut += decalage;
  m_tfin   += decalage;
}




// ----------------------------------------------------------------------------
// Identite entre deux chargements
// G.FERRER - Fevr.2002 - Creation
bool ObstacleChargement::operator == ( const ObstacleChargement &chargement ) 
  const
{
  return ( this == &chargement );
}




// ----------------------------------------------------------------------------
// Ordre entre deux chargements
// G.FERRER - Juil.2003 - Creation
bool operator < (const ObstacleChargement &c0,
		 const ObstacleChargement &c1)
{
  return c0.m_tdebut < c1.m_tdebut;
}




// ----------------------------------------------------------------------------
// Ecriture du chargement d'un obstacle
// G.FERRER - Dece.2000 - Creation
ostream &operator << ( ostream &fileOut, 
	ObstacleChargement &chargement )
{
  fileOut << chargement.m_nomObstacle << '\n';
  fileOut << chargement.m_tdebut << '\t' << chargement.m_tfin << '\n';
  fileOut << chargement.m_vitesse_translation;
  fileOut << chargement.m_vitesse_rotation;

  return( fileOut );
}




// ----------------------------------------------------------------------------
// Lecture du chargement d'un obstacle
// G.FERRER - Mars.2000 - Creation
istream &operator >> ( istream &fileIn, 
	ObstacleChargement &chargement )
{
  fileIn >> chargement.m_nomObstacle;
  fileIn >> chargement.m_tdebut >> chargement.m_tfin;
  fileIn >> chargement.m_vitesse_translation;
  if ( chargement.m_vitesse_translation != VecteurNul ) 
    chargement.m_type = "Translation";
  fileIn >> chargement.m_vitesse_rotation;
  if ( chargement.m_vitesse_rotation != VecteurNul ) 
    chargement.m_type = "Rotation";  

  return ( fileIn );
}




// ----------------------------------------------------------------------------
// Debug
void ObstacleChargement::debug( char* c )
{
  cout << m_tdebut << '-' << m_tfin << '\t';
}




// ----------------------------------------------------------------------------
// Type de chargement
// A.WACHS - Fev.2011 - Creation
string ObstacleChargement::getType() const
{
  return m_type;
}
 
