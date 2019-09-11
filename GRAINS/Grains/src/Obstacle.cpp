#include "Obstacle.H"
#include "CompObstacle.H"
#include "MonObstacle.H"
#include "LinkedCell.H"
#include "ObstacleChargement.H"
#include "Quaternion.H"
#include "App.H"
#include "Memento.hh"
#include <math.h>


int Obstacle::m_nbreObstacles = 0;
bool Obstacle::m_DeplaceObstacle = true ;
bool Obstacle::m_isConfinement = false;


//-----------------------------------------------------------------------------
// Constructeur par defaut
// G.FERRER - Nove.1999 - Creation
// D.RAKOTONIRINA - juin 2014 - Modification
Obstacle::Obstacle( const string &s, const bool &autonumbering ) :
  Composant( autonumbering ),
  m_nom( s ),
  m_deplace( false ),
  m_indicator( 0. ),
  m_ObstacleType ( "0" )
{}




// ----------------------------------------------------------------------------
// Constructeur par copie d'un Composant
// G.FERRER - Juil.2003 - Creation
Obstacle::Obstacle( Composant &copie, char const* s ) :
  Composant( copie ),
  m_nom( s ),
  m_deplace( false ),
  m_indicator( 0. ),
  m_ObstacleType ( "0" )
{}




//-----------------------------------------------------------------------------
// G.FERRER - Nove.1999 - Creation
// Destructeur
Obstacle::~Obstacle()
{}




// ----------------------------------------------------------------------------
// Associaton du chargement a l'obstacle
// G.FERRER - Nove.2000 - Creation
bool Obstacle::Associer( ObstacleChargement &chargement )
{
  bool status = false;
  if ( m_nom == chargement.getNom() )
  {
    m_cinematique.append( chargement );
    status = true;
  }

  return status;
}




// ----------------------------------------------------------------------------
// Renvoie le torseur applique
Torseur const* Obstacle::getTorseur()
{
  return &m_somme;
}




// ----------------------------------------------------------------------------
// Associaton du chargement a l'obstacle
// G.FERRER - Aout.2003 - Creation
bool Obstacle::Associer( ObstacleChargement_F &chargement )
{
  bool status = false;
  if ( m_nom == chargement.getNom() )
  {
    m_confinement.append( chargement );
    Obstacle::m_isConfinement = status = true;
  }
  return status;
}




// ----------------------------------------------------------------------------
// Description de la cinematique a partir de la cinematique en argument,
// appliquee au bras de levier.
// G.FERRER - Octo.2002 - Creation
void Obstacle::Decompose( const CineObstacle &voisine,
	const Vecteur &levier )
{
  m_cinematique.Decompose( voisine, levier );
}




// ----------------------------------------------------------------------------
// Description de la cinematique a partir de la cinematique en argument,
// G.FERRER - Aout.2003 - Creation
void Obstacle::Decompose( const CineObstacle_F &voisine,
	const Point &centre )
{
  m_confinement.Decompose( voisine, centre );
}




// ----------------------------------------------------------------------------
// Nombre d'obstacles unitaires
// G.FERRER - Octo.2002 - Creation
int Obstacle::getNombre()
{
  return Obstacle::m_nbreObstacles;
}




// ----------------------------------------------------------------------------
// Vitesse en un point de l'obstacle
// G.FERRER - Janv.2002 - Creation
// A.WACHS - Octo.2010 - Modif
Vecteur Obstacle::getVitesse( const Point &pt ) const
{
  Vecteur levier = pt - *m_geoFormeVdw->getCentre();
  // if m_confinement
  if ( Obstacle::m_isConfinement )
    return ( m_confinement.VitesseRelative( levier ) );
  else
    return m_cinematique.Vitesse( levier );
}




// ----------------------------------------------------------------------------
// Vitesse de rotation
// G.FERRER - Janv.2002 - Creation
Vecteur const* Obstacle::getVitesseRotation() const
{
  return m_cinematique.getVitesseRotation();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @brief Vitesse de translation. */
Vecteur const* Obstacle::getVitesseTranslation() const
{
  return m_cinematique.getVitesseTranslation();
}




// ----------------------------------------------------------------------------
// Remise a 0 de la cinematique
void Obstacle::resetCinematique()
{
  m_cinematique.reset();
  m_confinement.reset();
}




// ----------------------------------------------------------------------------
// Affectation d'une vitesse de Rotation.
// G.FERRER - Octo.2002 - Creation
void Obstacle::setCinematique( CineObstacle &cinematique_ )
{
  m_cinematique.set( cinematique_ );
}




// ----------------------------------------------------------------------------
// Affectation de la vitesse de l'obstacle
// A.WACHS - Fev.2012 - Creation
void Obstacle::setVitesse( Vecteur const* vitesseTranslation,
  	Vecteur const* vitesseRotation )
{
  m_cinematique.setVitesse( vitesseTranslation, vitesseRotation );
}




// ----------------------------------------------------------------------------
// Ecriture de l'information de position
// G.FERRER - Janv.2001 - Creation
void Obstacle::writePosition( ostream &position )
{
  position << "*Obstacle\n" << m_nom << '\n';
  Composant::writePosition( position );
  position << "*EndObstacle\n\n";
}




// ----------------------------------------------------------------------------
// Ecriture de l'information statique
// G.FERRER - Janv.2001 - Creation
void Obstacle::writeStatique( ostream &statique )
{
  statique << "*Obstacle\n" << m_nom << '\n';
  Composant::writeStatique( statique );
  statique << "*EndObstacle\n\n";
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecriture de l'identite de l'obstacle
// G.FERRER - Janv.2004 - Creation
void Obstacle::writeIdentity( ostream &file ) const
{
  file << m_nom << '\t' << this << '\n';
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie l'indicateur pour post-traitement Paraview
double Obstacle::getIndicator() const
{
  return m_indicator;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// L'obstacle s'est il deplace ?
bool Obstacle::hasMoved() const
{
  return m_deplace;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Sauvegarde de l'etat
void Obstacle::saveState()
{
  if (!m_memento)
    m_memento = new ConfigurationMemento();
  m_memento->m_position = *m_geoFormeVdw->getTransform();
  m_cinematique.saveState();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Restauration de l'etat
void Obstacle::restaureState()
{
  m_geoFormeVdw->setTransform( m_memento->m_position );
  m_cinematique.restaureState();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Deplacement geometrique des obstacles
void Obstacle::setDeplaceObstacle( bool const& depObs )
{
  Obstacle::m_DeplaceObstacle = depObs ;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Deplacement geometrique des obstacles
bool Obstacle::getDeplaceObstacle()
{
  return Obstacle::m_DeplaceObstacle ;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Contact entre un composant et un composant.
// D. RAKOTONIRINA - Fev. 2014 - Creation
void Obstacle::InterAction( Composant* voisin,
	double dt, double const& temps, LinkedCell *LC ) throw (ErreurContact)
{
  cout << "Warning when calling Obstacle::InterAction() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Contact entre un composant et un composant.
// D. RAKOTONIRINA - Fev. 2014 - Creation
void Obstacle::InterActionPostProcessing( Composant* voisin, Scalar dt,
	list<struct PointForcePostProcessing>* listOfContacts )
  throw (ErreurContact)
{
  cout << "Warning when calling Obstacle::InterActionPostProcessing() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Nombre de points pour post-processing avec Paraview
// D. RAKOTONIRINA - Dec. 2014 - Creation
int Obstacle::numberOfPoints_PARAVIEW() const
{
  cout << "Warning when calling Obstacle::numberOfPoints_PARAVIEW() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Nombre de points pour post-processing avec Paraview
// D. RAKOTONIRINA - Dec. 2014 - Creation
int Obstacle::numberOfCells_PARAVIEW() const
{
  cout << "Warning when calling Obstacle::numberOfCells_PARAVIEW() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Nombre de polygones elementaires pour post-processing avec Paraview
// D. RAKOTONIRINA - Dec. 2014 - Creation
list<Point> Obstacle::get_polygonsPts_PARAVIEW( Vecteur const* translation )
const
{
  cout << "Warning when calling Obstacle::get_polygonsPts_PARAVIEW() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Ecrit les points du convexe pour post-processing avec Paraview
// D. RAKOTONIRINA - Dec. 2014 - Creation
void Obstacle::write_polygonsPts_PARAVIEW( ostream &f,
  	Vecteur const* translation ) const
{
  cout << "Warning when calling Obstacle::write_polygonsPts_PARAVIEW() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Ecrit les points du convexe pour post-processing avec Paraview
// D. RAKOTONIRINA - Dec. 2014 - Creation
void Obstacle::write_polygonsStr_PARAVIEW(list<int> &connectivity,
    	list<int> &offsets, list<int> &cellstype, int& firstpoint_globalnumber,
	int& last_offset) const
{
  cout << "Warning when calling Obstacle::write_polygonsStr_PARAVIEW() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecriture de la position de la particule pour le Fluide
// D. RAKOTONIRINA - Avril. 2015 - Creation
void Obstacle::writePositionInFluid( ostream &fileOut )
{
  cout << "Warning when calling Obstacle::writePositionInFluid() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Initialize all contact map entries to false
void Obstacle::setContactMapToFalse()
{
  cout << "Warning when calling Obstacle::setContactMapToFalse() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Update contact map
void Obstacle::updateContactMap()
{
  cout << "Warning when calling Obstacle::updateContactMap() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Does the contact exist in the map, if yes return the pointer to the
// cumulative tangential displacement
bool Obstacle::ContactInMapIsActive( double* &tangentialDepl, int const& id )
{
  cout << "Warning when calling Obstacle::ContactInMapIsActive() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Add new contact in the map
void Obstacle::addNewContactInMap( double const& tangentialDepl, int const& id )
{
  cout << "Warning when calling Obstacle::addNewContactInMap() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Increase cumulative tangential displacement with component id
void Obstacle::addDeplContactInMap( double const& tangentialDepl,
	int const& id )
{
  cout << "Warning when calling Obstacle::addDeplContactInMap() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}
