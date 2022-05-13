#include "MPIWrapperGrains.hh"
#include "MonObstacle.H"
#include "LinkedCell.H"
#include "Grains.H"
#include "Contact_BuilderFactory.hh"
#include "GrainsCoupledWithFluid.hh"
#include "AppFluide_LubricationCorrection.H"
#include "Particule.H"
#include "PointContact.hh"
#include "SaveTable.H"
#include "Cellule.H"
#include "Grains_Exec.hh"
#include "Memento.hh"
#include "ContactLaw.hh"
#include <sstream>
#include <limits>
using namespace std;

#include <assert.h>


// ----------------------------------------------------------------------------
// Constructeur par defaut
// G.FERRER - Octo.2002 - Creation
MonObstacle::MonObstacle( const string &s ) :
  Obstacle( s ),
  m_transferToFluid( false )
{
  m_ObstacleType = "MonObstacle";
  Obstacle::m_nbreObstacles++;
  m_LinkUpdate_frequency = 1;
  m_LinkUpdate_counter = 0;
}




// ----------------------------------------------------------------------------
// Constructeur par decodage d'un noeud XML
MonObstacle::MonObstacle( DOMNode *root ) :
  Obstacle(),
  m_transferToFluid( false )
{
  m_ObstacleType = "MonObstacle";

  Obstacle::m_nbreObstacles++;

  m_nom = ReaderXML::getNodeAttr_String( root, "name" );

  // Convex - Position & Orientation
  m_geoFormeVdw = new FormeVdW( root );

  // Materiau
  DOMNode* materiau_ = ReaderXML::getNode( root, "Materiau" );
  m_nomMateriau = ReaderXML::getNodeValue_String( materiau_ );
  Contact_BuilderFactory::defineMaterial( m_nomMateriau, true );

  // Obstacle � transf�rer au fluide
  DOMNode* statut = ReaderXML::getNode( root, "Statut" );
  if ( statut )
    m_transferToFluid = ReaderXML::getNodeAttr_Int( statut, "ToFluid" );


  m_obstacleBox = Composant::Boite();
  m_LinkUpdate_frequency = 1;
  m_LinkUpdate_counter = 0;
}




// ----------------------------------------------------------------------------
// Constructeur par copie d'un Composant
// G.FERRER - Juil.2003 - Creation
MonObstacle::MonObstacle( Composant &copie, char const* s ) :
  Obstacle( copie, s )
{
  m_ObstacleType = "MonObstacle";

  Obstacle::m_nbreObstacles++;
  m_obstacleBox = Composant::Boite();
  m_LinkUpdate_frequency = 1;
  m_LinkUpdate_counter = 0;
  m_transferToFluid = false;
}




// ----------------------------------------------------------------------------
// Constructeur a partir d'une forme
// M.BERNARD - April 2015 - Creation
MonObstacle::MonObstacle( FormeVdW *geoFormeVdw, const string &name,
      const string &materialName, const bool &transferToFluid_ ) :
  Obstacle()
{
  m_ObstacleType = "MonObstacle";

  Obstacle::m_nbreObstacles++;

  m_nom = name;
  m_geoFormeVdw = geoFormeVdw;
  m_nomMateriau = materialName;
  m_transferToFluid = transferToFluid_;
  Contact_BuilderFactory::defineMaterial( materialName, true );

  m_obstacleBox = Composant::Boite();
  m_LinkUpdate_frequency = 1;
  m_LinkUpdate_counter = 0;
}




// ----------------------------------------------------------------------------
// Destructeur
// G.FERRER - Octo.2002 - Creation
MonObstacle::~MonObstacle()
{
  m_inCells.clear();
}




// ----------------------------------------------------------------------------
// Ajout d'un obstacle (unique ou composite) a l'arborescence Obstacle
// En fait, sans action dans le cas d'un obstacle unique.
void MonObstacle::append( Obstacle* obstacle )
{
  cout << "Warning : Obstacle append on unique obstacle !\n";
}




// ----------------------------------------------------------------------------
// Deplacement de l'obstacle du chargement indiquee
// G.FERRER - Dece.1999 - Creation
// G.FERRER - Fevr.2002 - Nouvelle version de cinematique
// D. RAKOTONIRINA - Apr. 2017 - Modification
list<MonObstacle*> MonObstacle::Deplacer( Scalar temps, Scalar dt,
	const bool &b_deplaceCine_Comp,
	const bool &b_deplaceF_Comp )
{

  m_deplace = m_cinematique.Deplacement( temps, dt );
  m_deplace = m_deplace || b_deplaceCine_Comp;

  if ( m_deplace && Obstacle::m_DeplaceObstacle )
  {
    Vecteur const* translation = m_cinematique.getTranslation();
    *m_geoFormeVdw += *translation;
    Quaternion const* w = m_cinematique.getQuaternionRotationOverDt();
    Rotate(*w);
  }

  bool deplaceF = m_confinement.Deplacement( temps, dt, this );
  deplaceF = deplaceF || b_deplaceF_Comp;

  if ( deplaceF && Obstacle::m_DeplaceObstacle )
  {
    Vecteur translation = m_confinement.getTranslation( dt );
    *m_geoFormeVdw += translation;
    //    Quaternion w = cinematique.getRotation(dt);
    //    Rotate(w);
  }
  m_deplace = m_deplace || deplaceF;

  list<MonObstacle*> obstacleEnDeplacement;
  if ( m_deplace )
  {
    m_obstacleBox = Composant::Boite();
    obstacleEnDeplacement.push_back( this );
  }

  return obstacleEnDeplacement;
}




// ----------------------------------------------------------------------------
// L'obstacle a t'il le nom indique.
// G.FERRER - Juil.2003 - Creation
const Obstacle* MonObstacle::getNom( const string &nom_ ) const
{
  if ( m_nom == nom_ )
    return (Obstacle*)this;
  else
    return NULL;
}




// ============================================================================
// Liste des obstacles.
list<MonObstacle*> MonObstacle::getObstacles()
{
  list<MonObstacle*> liste;
  liste.push_back(this);
  return liste;
}




// ============================================================================
// Liste des obstacles.
//list<MonObstacle*> MonObstacle::getObstaclesToFluid() // Modif Manu 04/2015
list<Obstacle*> MonObstacle::getObstaclesToFluid()
{
//  list<MonObstacle*> liste;
  list<Obstacle*> liste;
  if ( m_transferToFluid == true ) liste.push_back(this);
  return liste;
}




// ----------------------------------------------------------------------------
// Determination du contact avec une particule
// DGR : Le contact est dirige Particule->Obstacle cf. Obstacle::getContact
void MonObstacle::InterAction( Composant* voisin,
	double dt, double const& temps, LinkedCell *LC )
  throw (ErreurContact)
{
  bool contactProbable = m_obstacleBox.InZone( voisin->getPosition(),
    	voisin->getRayon() );

  PointContact closestPoint;
  if ( contactProbable )
  {
    try {
      closestPoint = voisin->getForme()->ClosestPoint( *m_geoFormeVdw );
    }
    catch (ErreurContact &erreur)
    {
      try {
	closestPoint = voisin->getForme()->ClosestPoint_ErreurHandling(
		*m_geoFormeVdw, 10., m_id, voisin->getID() );
      }
      catch (ErreurContact &erreur_level2)
      {
        cout << endl << "Processor = " <<
    		(Grains_Exec::m_MPI ? Grains_Exec::getComm()->rank_ACTIV() : 0 )
		<< " has thrown an ErreurContact exception" <<  endl;
        erreur_level2.setMessage(
		"MonObstacle::InterAction : choc de croute ! a t="
		+Grains_Exec::doubleToString(temps,TIMEFORMAT));
        erreur_level2.setComposants( this, voisin, temps );
        Grains_Exec::m_exception_Contact = true;
        throw(erreur_level2);
      }
    }
  }
  else
    closestPoint = PointNoContact;

  LC->addToContactsFeatures( temps, closestPoint );

  // cout << "\tDistance GJK = " << closestPoint.getDistance() << endl ;

  if ( closestPoint.getDistance() < 0. )
  {
    if ( Contact_BuilderFactory::contactForceModel(
    	m_nomMateriau, voisin->materiau() )
    	->computeForces( voisin, this, closestPoint, LC, dt ) )
      voisin->addToCoordinationNumber( 1 );
  }
  else if ( Grains_Exec::m_withlubrication )
    (GrainsCoupledWithFluid::LubricationCorrection())->computeforce( voisin,
		     this, LC, dt );
 }



// ----------------------------------------------------------------------------
// Determination du contact avec une particule
void MonObstacle::InterActionPostProcessing( Composant* voisin, Scalar dt,
	list<struct PointForcePostProcessing>* listOfContacts )
  throw (ErreurContact)
{
  bool contactProbable = m_obstacleBox.InZone( voisin->getPosition(),
    	voisin->getRayon() );

  PointContact closestPoint;
  if ( contactProbable )
  {
    try {
      closestPoint = voisin->getForme()->ClosestPoint( *m_geoFormeVdw );
    }
    catch (ErreurContact &erreur)
    {
      closestPoint = voisin->getForme()->ClosestPoint_ErreurHandling(
		*m_geoFormeVdw, 10., m_id, voisin->getID() );
    }
  }
  else
    closestPoint = PointNoContact;

  if ( closestPoint.getDistance() < 0. )
    Contact_BuilderFactory::contactForceModel(
    	m_nomMateriau, voisin->materiau() )
    	->computeForcesPostProcessing( voisin, this, dt, closestPoint,
	listOfContacts );
}




// ----------------------------------------------------------------------------
// Lecture de l'obstacle
// G.FERRER - Octo.2002 - Creation
void MonObstacle::reload( istream &fileSave )
{
  string buffer, adresse;
  Scalar buf = 0.;

  fileSave >> m_nom >> adresse >> buffer;
  SaveTable::create( adresse, this );
  fileSave >> buffer >> m_nomMateriau;
  m_geoFormeVdw = new FormeVdW( fileSave );
  fileSave >> buffer;
  if ( buffer == "*Couleur" )
    fileSave >> buf >> buf >> buf;
  else if ( buffer == "*ToFluid" )
    fileSave >> m_transferToFluid;
  m_geoFormeVdw->readPosition( fileSave );
  m_obstacleBox = Composant::Boite();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Reload de l'obstacle & publication dans son referent
void MonObstacle::reload( Obstacle &obstacle, istream &file )
{
  string buffer, adresse;
  Scalar buf = 0.;
  file >> adresse >> m_id;
  SaveTable::create( adresse, this );
  file >> buffer >> m_nomMateriau;
  m_geoFormeVdw = new FormeVdW(file);
  file >> buffer;
  if ( buffer == "*Couleur" )
    file >> buf >> buf >> buf;
  else if ( buffer == "*ToFluid" )
    file >> m_transferToFluid;
  m_geoFormeVdw->readPosition( file );
  obstacle.append( this );
  file >> buffer;
  assert(buffer == "</Simple>");
  m_obstacleBox = Composant::Boite();
}




// ----------------------------------------------------------------------------
// Rotation de l'obstacle
// D.PETIT - Aout.2000 - Creation
void MonObstacle::Rotate( const Quaternion &rotation )
{
  m_geoFormeVdw->Rotate( rotation );
  m_obstacleBox = Composant::Boite();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Suppression de l'obstacle dans la zone specifie.
// Pour un monolithe -> pas d'action.
// G.FERRER - Fevr.2004 - Creation
void MonObstacle::Suppression( const BBox &box )
{
}




// ----------------------------------------------------------------------------
// Translation de l'obstacle
// D.PETIT - Aout.2000 - Creation
void MonObstacle::Translate( const Vecteur &translation )
{
  *m_geoFormeVdw += translation;
  m_obstacleBox = Composant::Boite();
}




// ----------------------------------------------------------------------------
// Sauvegarde de l'obstacle pour Reload
// G.FERRER - Octo.2002 - Creation
// D. RAKOTONIRINA - Oct. 2014 - Modification
void MonObstacle::write( ostream &fileSave, Composant const* composant ) const
{
  fileSave << "<Simple>\t" << m_nom << '\n'
	   << this << '\t' << m_id << '\n';
  fileSave << "*Materiau \n"
	   << m_nomMateriau << '\n';
  m_geoFormeVdw->writeStatique( fileSave );
  fileSave << "*ToFluid\n" << m_transferToFluid << '\n';
  m_geoFormeVdw->writePosition( fileSave );
  fileSave << endl;
  fileSave << "</Simple>\n";
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ajout d'une cellule � la liste de cellules auxquelles l'obstacle est lie
// A.WACHS - Aout.2009 - Creation
void MonObstacle::add( Cellule *cel_ )
{
  m_inCells.push_back(cel_);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Vide la liste des cellules auxquelles l'obstacle est lie et
// supprime le pointeur sur l'obstacle dans les cellules de la liste
// A.WACHS - Aout.2009 - Creation
void MonObstacle::resetInCells()
{
  list<Cellule*>::iterator il;
  for (il=m_inCells.begin();il!=m_inCells.end();il++) (*il)->remove( this );
  m_inCells.clear();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie la liste des cellules auxquelles l'obstacle est lie
// A.WACHS - Aout.2009 - Creation
const list<Cellule*>* MonObstacle::getInCells() const
{
  return &m_inCells;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie la BBox de l'obstacle
// A.WACHS - Aout.2009 - Creation
const BBox* MonObstacle::getObstacleBox() const
{
  return &m_obstacleBox;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie le type de l'obstacle
//
string MonObstacle::getObstacleType()
{
  return m_ObstacleType;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Destruction de l'obstacle dans l'arborescence des obstacles
// A.WACHS - Aout.2010 - Creation
void MonObstacle::DestroyObstacle( const string &name_ )
{
  if ( m_nom == name_ || name_ == "ToBeDestroyed" )
    Obstacle::m_nbreObstacles--;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Supprime l'obstacle de nom name_
// A.WACHS - Aout.2010 - Creation
void MonObstacle::SupprimeObstacle( const string &name_, LinkedCell *LC )
{
  if ( m_nom == name_ || name_ == "ToBeErased" )
    // Suppression de l'obstacle du LinkedCell
    LC->remove( this );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Indique si il y a lieu de mettre a jour le lien entre l'obstacle et
// le LinkedCell: si oui, la methode renvoie true et met le compteur a zero, si
// non la methode renvoie false et met le compteur est incremente
bool MonObstacle::performLinkUpdate()
{
  bool b_doit = false;

  if ( m_LinkUpdate_counter == m_LinkUpdate_frequency )
  {
    m_LinkUpdate_counter = 1;
    b_doit = true;
  }
  else ++m_LinkUpdate_counter;

  return b_doit;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Maximum de la valeur absolue de vitesse de l'obstacle
// (en parcourant les points de sa surface) dans chaque direction
Vecteur MonObstacle::vitesseMaxPerDirection() const
{
  list<Point> surface = m_geoFormeVdw->get_polygonsPts_PARAVIEW();
  list<Point>::iterator iv;
  Vecteur vmax, v;

  for (iv=surface.begin();iv!=surface.end();iv++)
  {
    v = getVitesse(*iv);
    vmax[X] = fabs(v[X]) > vmax[X] ? fabs(v[X]) : vmax[X];
    vmax[Y] = fabs(v[Y]) > vmax[Y] ? fabs(v[Y]) : vmax[Y];
    vmax[Z] = fabs(v[Z]) > vmax[Z] ? fabs(v[Z]) : vmax[Z];
  }

  return vmax;

}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Fr�quence de mise � jour du lien entre l'obstacle et le LinkedCell
void MonObstacle::setObstacleLinkedCellUpdateFrequency( int const &updateFreq )
{
  m_LinkUpdate_frequency = updateFreq;
  m_LinkUpdate_counter = m_LinkUpdate_frequency;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Cree et ajoute l'etat
void MonObstacle::createState( list<struct ObstacleState*> &obsStates )
{
  struct ObstacleState* obss = new ObstacleState;
  obss->nom = m_nom;
  obss->memento_config = new ConfigurationMemento();
  obss->memento_config->m_position = *m_geoFormeVdw->getTransform();
  obss->memento_cine = m_cinematique.createState();
  obsStates.push_back( obss );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Restauration de l'etat
void MonObstacle::restaureState( list<struct ObstacleState*>& obsStates )
{
  list<struct ObstacleState*>::iterator il;
  bool found = false ;
  for (il=obsStates.begin();il!=obsStates.end() && !found; )
    if ( (*il)->nom == m_nom )
    {
      m_geoFormeVdw->setTransform((*il)->memento_config->m_position);
      m_cinematique.restaureState((*il)->memento_cine);
      delete (*il)->memento_config;
      delete (*il)->memento_cine;
      delete *il;
      il = obsStates.erase(il) ;
      found = true ;
    }
    else il++;
}




// ----------------------------------------------------------------------------
// Nombre de points pour post-processing avec Paraview
// D. RAKOTONIRINA - Dec. 2014 - Creation
int MonObstacle::numberOfPoints_PARAVIEW() const
{
  return ( m_geoFormeVdw->getConvex()->numberOfPoints_PARAVIEW() );
}




// ----------------------------------------------------------------------------
// Nombre de points pour post-processing avec Paraview
// D. RAKOTONIRINA - Dec. 2014 - Creation
int MonObstacle::numberOfCells_PARAVIEW() const
{
  return ( m_geoFormeVdw->getConvex()->numberOfCells_PARAVIEW() );
}




// ----------------------------------------------------------------------------
// Nombre de polygones elementaires pour post-processing avec Paraview
// D. RAKOTONIRINA - Dec. 2014 - Creation
list<Point> MonObstacle::get_polygonsPts_PARAVIEW( Vecteur const* translation )
const
{
  return ( getForme()->get_polygonsPts_PARAVIEW( translation ) );
}




// ----------------------------------------------------------------------------
// Ecrit les points du convexe pour post-processing avec Paraview
// D. RAKOTONIRINA - Dec. 2014 - Creation
void MonObstacle::write_polygonsPts_PARAVIEW( ostream &f,
  	Vecteur const* translation ) const
{
  getForme()->write_polygonsPts_PARAVIEW( f, translation );
}




// ----------------------------------------------------------------------------
// Ecrit les points du convexe pour post-processing avec Paraview
// D. RAKOTONIRINA - Dec. 2014 - Creation
void MonObstacle::write_polygonsStr_PARAVIEW(list<int> &connectivity,
    	list<int> &offsets, list<int> &cellstype, int& firstpoint_globalnumber,
	int& last_offset) const
{
  m_geoFormeVdw->getConvex()->write_polygonsStr_PARAVIEW(connectivity,
	offsets, cellstype, firstpoint_globalnumber, last_offset);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecriture de la position de la particule pour le Fluide
// D. RAKOTONIRINA - Avril. 2015 - Creation
void MonObstacle::writePositionInFluid( ostream &fileOut )
{
  m_geoFormeVdw->writePositionInFluid( fileOut );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Initialize all contact map entries to false
void MonObstacle::setContactMapToFalse()
{
  Composant::setContactMapToFalse();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Update contact map
void MonObstacle::updateContactMap()
{
  Composant::updateContactMap();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Does the contact exist in the map, if yes return the pointer to the
// cumulative tangential displacement
bool MonObstacle::getContactMemory( std::tuple<int,int,int> const& id,
  Vecteur* &tangent, Vecteur* &prev_normal, Vecteur* &cumulSpringTorque )
{
  return ( Composant::getContactMemory( id, tangent, prev_normal,
    cumulSpringTorque) );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Add new contact in the map
void MonObstacle::addNewContactInMap( std::tuple<int,int,int> const& id,
  Vecteur const& tangent, Vecteur const& prev_normal,
  Vecteur const& cumulSpringTorque )
{
  Composant::addNewContactInMap( id, tangent, prev_normal, cumulSpringTorque );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Increase cumulative tangential displacement with component id
void MonObstacle::addDeplContactInMap( std::tuple<int,int,int> const& id,
  Vecteur const& tangent, Vecteur const& prev_normal,
  Vecteur const& cumulSpringTorque )
{
  Composant::addDeplContactInMap( id, tangent, prev_normal, cumulSpringTorque );
}


void MonObstacle::copyHistoryContacts( double* &destination, int start_index )
{
  Composant::copyHistoryContacts( destination, start_index ) ;
}

// ----------------------------------------------------------------------------
// Copy existing contact in the map
void MonObstacle::copyContactInMap( std::tuple<int,int,int> const& id,
  bool const& isActive, Vecteur const& tangent, Vecteur const& prev_normal,
  Vecteur const& cumulSpringTorque )
{
  Composant::copyContactInMap( id, isActive, tangent, prev_normal,
    cumulSpringTorque ) ;
}

int MonObstacle::getContactMapSize()
{
  return ( Composant::getContactMapSize() );
}

void MonObstacle::updateContactMapId( int prev_id, int new_id)
{
  Composant::updateContactMapId( prev_id, new_id);
}
