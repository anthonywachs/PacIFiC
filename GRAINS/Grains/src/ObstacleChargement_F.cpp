#include "ObstacleChargement_F.H"
#include "MPIWrapperGrains.hh"
#include "Grains_Exec.hh"
#include "Obstacle.H"
#include <stdlib.h>

Vecteur ObstacleChargement_F::m_prev = VecteurNul;
// ----------------------------------------------------------------------------
// Constructeur par defaut
// D. RAKOTONIRINA - Fev.2017 - Creation
ObstacleChargement_F::ObstacleChargement_F( DOMNode* root, double dt, int rank )
{
  m_nomObstacle = ReaderXML::getNodeAttr_String( root, "NomObstacle" );
  
  DOMNode* temps    = ReaderXML::getNode( root, "Temps" );
  DOMNode* force    = ReaderXML::getNode( root, "Amplitude" );
  DOMNode* nVecteur = ReaderXML::getNode( root, "Vecteur" );
  DOMNode* property = ReaderXML::getNode( root, "Property" );
  m_masse	    = ReaderXML::getNodeAttr_Double( property, "Masse" );

  m_direction[X] = ReaderXML::getNodeAttr_Double( nVecteur, "X" );
  m_direction[Y] = ReaderXML::getNodeAttr_Double( nVecteur, "Y" );
  m_direction[Z] = ReaderXML::getNodeAttr_Double( nVecteur, "Z" );

  m_tdebut          = ReaderXML::getNodeAttr_Double( temps, "Debut" );
  m_tfin            = ReaderXML::getNodeAttr_Double( temps, "Fin" );

  m_type            = ReaderXML::getNodeAttr_String( root, "Type" );

  m_force[X]   	    = ReaderXML::getNodeAttr_Double( force, "AX" );
  m_force[Y]   	    = ReaderXML::getNodeAttr_Double( force, "AY" );
  m_force[Z]   	    = ReaderXML::getNodeAttr_Double( force, "AZ" );

  if ( m_type == "Cyclic")
  {
    DOMNode* frequence  = ReaderXML::getNode( root, "Frequence" );
    m_phase             = ReaderXML::getNodeAttr_Double( frequence, "Phi" );
    m_frequenceX        = ReaderXML::getNodeAttr_Double( frequence, "FX" );
    m_frequenceY        = ReaderXML::getNodeAttr_Double( frequence, "FY" );
    m_frequenceZ        = ReaderXML::getNodeAttr_Double( frequence, "FZ" );
    m_phase            *= PI / 180.;
  }

  if ( rank == 0 )
  {
    cout << "Chargement en force sur " << m_nomObstacle << endl;
    cout << "Type de chargement : " << m_type << endl;
    cout << "Amplitude de la force = " << m_force[X] << "\t" << m_force[Y] 
	 << "\t" << m_force[Z] << endl;
    cout << "Direction de la force = " << m_direction[X] << "\t"
	 << m_direction[Y] << "\t" << m_direction[Z] << endl;
    if ( m_type == "Cyclic" )
      cout << "Frequence : FX = " << m_frequenceX << "\tFY = " << m_frequenceY
           << "\tFZ = " << m_frequenceZ << endl;
    cout << "   Temps de depart = " << m_tdebut << endl;
    cout << "   Temps de fin = " << m_tfin << endl;
  }
}




// ----------------------------------------------------------------------------
// Constructeur par defaut
// G.FERRER - Aout.2003 - Creation
ObstacleChargement_F::ObstacleChargement_F()
{
}




// ----------------------------------------------------------------------------
// Destructeur
// G.FERRER - Aout.2003 - Creation
ObstacleChargement_F::~ObstacleChargement_F()
{}




// ----------------------------------------------------------------------------
// Identification de l'obstacle
// G.FERRER - Aout.2003 - Creation
string ObstacleChargement_F::getNom() const
{
  return m_nomObstacle;
}




// ----------------------------------------------------------------------------
// Temps de chargement dans l'espace indique
// G.FERRER - aout.2003 - Creation
Scalar ObstacleChargement_F::getTime( Scalar debut, Scalar fin ) const
{
  Scalar dt = fin - debut;

  if ( debut < m_tdebut )
    dt -= ( m_tdebut - debut );
  if ( m_tfin < fin )
    dt -= ( fin - m_tfin );

  return (dt);
}




// ----------------------------------------------------------------------------
// Le chargement est il actif ?
// G.FERRER - Aout.2003 - Creation
bool ObstacleChargement_F::isActif( Scalar t, Scalar dt ) const
{
  return ( t > m_tdebut - dt * 1.e-5  && t < m_tfin + dt * 1.e-5 );
}




// ----------------------------------------------------------------------------
// Le chargement est il termine?
// A.WACHS - Mars.2012 - Creation
bool ObstacleChargement_F::isCompleted( Scalar t, Scalar dt ) const 
{
  return ( t > m_tfin + dt * 1.e-5 );
}




// ----------------------------------------------------------------------------
// Construction du chargement par lecture des donnees
// G.FERRER - Aout.2003 - Creation
ObstacleChargement_F* ObstacleChargement_F::read( istream &fileIn )
{
  ObstacleChargement_F* chargement;
  chargement = new ObstacleChargement_F();

  fileIn >> chargement->m_nomObstacle
	>> chargement->m_tdebut >> chargement->m_tfin
	>> chargement->m_force 
	>> chargement->m_masse
	>> chargement->m_direction;
  
  return ( chargement );
}




// ----------------------------------------------------------------------------
// Renvoie de la force
// D. RAKOTONIRINA - Fev.2017 - Creation
Vecteur ObstacleChargement_F::getForce() const
{
  return m_force;
}




// ----------------------------------------------------------------------------
// Renvoie de la masse fictive
// D. RAKOTONIRINA - Fev.2017 - Creation
Scalar ObstacleChargement_F::getMasse() const
{
  return m_masse;
}




// ----------------------------------------------------------------------------
// Renvoie de la direction
// D. RAKOTONIRINA - Fev.2017 - Creation
Vecteur const* ObstacleChargement_F::getDirection() const
{
  return &m_direction;
}




// ----------------------------------------------------------------------------
// Type de chargement
// D. RAKOTONIRINA - Mai.2017 - Creation
string ObstacleChargement_F::getType() const
{
  return m_type;
}




// ----------------------------------------------------------------------------
// Vitesse de deplacement
// D.RAKOTONIRINA - Mai.2017 - Creation
Vecteur const* ObstacleChargement_F::VitesseTranslation( Scalar temps, 
	Scalar dt, Obstacle* obstacle )
{
  Vecteur center = *obstacle->getPosition();
  MPIWrapperGrains* wrapper = Grains_Exec::getComm();
  // Somme des forces sur l'obstacle
  Torseur const* somme  = obstacle->getTorseur();
  Vecteur const* forces = somme->getForce();
  Vecteur force = *forces;
  force[X] = wrapper->sum_DOUBLE_master( force[X] ); 
  force[Y] = wrapper->sum_DOUBLE_master( force[Y] ); 
  force[Z] = wrapper->sum_DOUBLE_master( force[Z] ); 

  force[X] = wrapper->Broadcast_DOUBLE( force[X] );
  force[Y] = wrapper->Broadcast_DOUBLE( force[Y] );
  force[Z] = wrapper->Broadcast_DOUBLE( force[Z] );
  
  Vecteur dforce;
  Vecteur depl;
  Vecteur trans;
  if ( m_type == "Translation" )
  {
    dforce = m_force - force;
    dforce[X] *= m_direction[X]; 
    dforce[Y] *= m_direction[Y]; 
    dforce[Z] *= m_direction[Z]; 
    depl = 0.5*( dt * dt / m_masse ) * dforce; 
    m_vitesse_translation = depl / dt;
  }
  else if ( m_type == "Cyclic" )
  {
    dforce = cyclicForce( temps ) - force;
    dforce[X] *= m_direction[X]; 
    dforce[Y] *= m_direction[Y]; 
    dforce[Z] *= m_direction[Z]; 
    trans = 0.5*( dt * dt / m_masse ) * dforce; 
    depl = trans - m_prev;
    m_vitesse_translation = depl / dt;
    // t-dt 
    m_prev = trans;
  }
  return &m_vitesse_translation;
}




// ----------------------------------------------------------------------------
// Calcul de la force appliquee a l'obstacle
Vecteur ObstacleChargement_F::cyclicForce( Scalar temps ) const
{
  Vecteur cycForce;
  cycForce[X] = m_force[X] * sin( 2. * PI * m_frequenceX * 
      	( temps - m_tdebut ) );
  cycForce[Y] = m_force[Y] * sin( 2. * PI * m_frequenceY *
      	( temps - m_tdebut ) + m_phase );
  cycForce[Z] = m_force[Z] * sin( 2. * PI * m_frequenceZ *
      	( temps - m_tdebut ) + m_phase );

  return cycForce;
}
