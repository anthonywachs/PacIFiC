// G.FERRER - Octo.2002 - Creation
// ============================================================================
#include "CompObstacle.H"
#include "Obstacle_BuilderFactory.H"
#include "Memento.hh"
#include "LinkedCell.H"
#include "Torseur.H"
#include "PointC.H"


// ----------------------------------------------------------------------------
// Constructeur par defaut
// G.FERRER - Octo.2002 - Creation
// D.RAKOTONIRINA - juin 2014 - Modification
CompObstacle::CompObstacle( const string &s ) :
  Obstacle( s, false )
{
  m_id = -4;
  m_geoFormeVdw = new FormeVdW( new PointC(), Transform() );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur du groupe d'obstacles
CompObstacle::CompObstacle( DOMNode* root ) :
  Obstacle( "obstacle", false )
{
  m_id = -4;
  assert(root != NULL);

  m_geoFormeVdw = new FormeVdW( new PointC(), Transform() );

  m_nom = ReaderXML::getNodeAttr_String( root, "name" );

  Obstacle *obstacle = NULL;
  DOMNodeList* allObstacles = ReaderXML::getNodes(root);
  for (XMLSize_t i=0; i<allObstacles->getLength(); i++)
  {
    obstacle = Obstacle_BuilderFactory::create( allObstacles->item( i ) );
    m_obstacles.push_back(obstacle);
  }
  EvalPosition();
}




// ----------------------------------------------------------------------------
// Destructeur
// G.FERRER - Octo.2002 - Creation
CompObstacle::~CompObstacle()
{
  Obstacle *obstacle;
  list<Obstacle*>::iterator iter;
  for (iter=m_obstacles.begin(); iter!=m_obstacles.end(); iter++)
  {
    obstacle = *iter;
    delete obstacle;
  }
}




// ----------------------------------------------------------------------------
// Ajout d'un obstacle (unique ou composite) a l'arborescence Obstacle.
// G.FERRER - Juil.2003 - Creation
void CompObstacle::append( Obstacle* obstacle )
{
  m_obstacles.push_back(obstacle);
}




// ----------------------------------------------------------------------------
// Associaton du chargement a l'obstacle
// Le chargement est lie a l'obstacle composite ou a un des composants.
// Cette association est definie par le nom de l'obstacle dans le chargement.
// G.FERRER - Octo.2002 - Creation
bool CompObstacle::Associer( ObstacleChargement &chargement )
{
  bool status = false;
  if ( m_nom == chargement.getNom() )
  {
    m_cinematique.append( chargement );
    status = true;
  }
  else
  {
    list<Obstacle*>::iterator obstacle;
    for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); obstacle++)
      status = (*obstacle)->Associer( chargement );
  }

  return status;
}




// ----------------------------------------------------------------------------
// Associaton du chargement a l'obstacle
// Le chargement est lie a l'obstacle composite ou a un des composants.
// Cette association est definie par le nom de l'obstacle dans le chargement.
// G.FERRER - Octo.2002 - Creation
bool CompObstacle::Associer( ObstacleChargement_F &chargement )
{
  bool status = false;
  if ( m_nom == chargement.getNom() )
  {
    m_confinement.append( chargement );
    status = true;
  }
  else
  {
    list<Obstacle*>::iterator obstacle;
    for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); obstacle++)
      status = (*obstacle)->Associer( chargement );
  }

  return status;
}




// ----------------------------------------------------------------------------
// Deplacement de l'obstacle ou de ses sous-composants dans l'espace de temps
// On evalue la vitesse de l'obstacle composite
// On affecte cette vitesse aux composants
// On evalue la vitesse de chaque composant dans l'espace temps
// G.FERRER - Octo.2002 - Creation
list<MonObstacle*> CompObstacle::Deplacer( Scalar temps, Scalar dt,
	const bool &b_deplaceCine_Comp,
	const bool &b_deplaceF_Comp )
{
  m_deplace = m_cinematique.Deplacement( temps, dt );
  m_deplace = m_deplace || b_deplaceCine_Comp;

  // Deplacement du centre du composite
  if ( m_deplace && Obstacle::m_DeplaceObstacle )
  {
    Vecteur const* translation = m_cinematique.getTranslation();
    *m_geoFormeVdw += *translation;
    /* Rotation non utilisee pour le centre du composite.
      Quaternion w = cinematique.getQuaternionRotationOverDt();
      Rotate(w);
    */
  }
  Point centre = *getPosition();

  // Deplacement des obstacles
  list<Obstacle*>::iterator obstacle;
  if ( m_deplace )
  {
    Vecteur levier;
    for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); obstacle++)
    {
      levier = *(*obstacle)->getPosition() - centre;
      (*obstacle)->Decompose( m_cinematique, levier );
    }
  }


  bool deplaceF = m_confinement.Deplacement( temps, dt, this );
  deplaceF = deplaceF || b_deplaceF_Comp;

  if ( deplaceF )
    for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); obstacle++)
      (*obstacle)->Decompose( m_confinement, *(*obstacle)->getPosition() );

  m_deplace = m_deplace || deplaceF; // ??? demander � Gillos !!

  list<MonObstacle*> obstacleEnDeplacement;
  list<MonObstacle*>::iterator ilo;
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); obstacle++) {
    list<MonObstacle*> lod = (*obstacle)->Deplacer( temps, dt, m_deplace,
    	deplaceF );
    for (ilo=lod.begin();ilo!=lod.end();ilo++)
      obstacleEnDeplacement.push_back(*ilo);
  }

  return obstacleEnDeplacement;
}




// ----------------------------------------------------------------------------
// Definition du centre de gravite
// G.FERRER - Octo.2002 - Creation
void CompObstacle::EvalPosition()
{
  Point centre;
  int nbre=0;
  list<Obstacle*>::iterator obstacle;
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end();
      nbre++, obstacle++)
    centre += *(*obstacle)->getPosition();
  centre /= nbre;
  setPosition(centre);

}




// ----------------------------------------------------------------------------
// L'obstacle a t'il le nom indique.
// G.FERRER - Juil.2003 - Creation
const Obstacle* CompObstacle::getNom( const string &nom_ ) const
{
  const Obstacle *obst = NULL;
  if ( m_nom == nom_ ) obst = (Obstacle*)this;
  else
  {
    list<Obstacle*>::const_iterator obstacle;
    for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end() &&
    	obst == NULL; obstacle++)
      obst = (*obstacle)->getNom( nom_ );
  }

  return obst;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Liste des obstacles.
list<MonObstacle*> CompObstacle::getObstacles()
{
  list<MonObstacle*> liste;

  list<Obstacle*>::iterator obstacle;
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); obstacle++)
  {
    list<MonObstacle*> tmp = (*obstacle)->getObstacles();
    liste.insert(liste.end(), tmp.begin(), tmp.end());
  }

  return liste;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Liste des obstacles a transmettre au fluide
// Modified by Manuel, May 2015, from list<MonObstable> to list<Obstacle>
//list<MonObstacle*> CompObstacle::getObstaclesToFluid()
list<Obstacle*> CompObstacle::getObstaclesToFluid()
{
//  list<MonObstacle*> liste;
  list<Obstacle*> liste;

  list<Obstacle*>::iterator obstacle;
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); obstacle++)
  {
//    list<MonObstacle*> tmp = (*obstacle)->getObstaclesToFluid();
    list<Obstacle*> tmp = (*obstacle)->getObstaclesToFluid();
    liste.insert(liste.end(), tmp.begin(), tmp.end());
  }

  return liste;
}




// ----------------------------------------------------------------------------
// Les deux composants sont ils en contact ?
// A.WACHS - aout 2009 - Creation
// D.RAKOTONIRINA - janv 2014 - Modification
bool CompObstacle::isContact( const Composant* voisin ) const
{
  bool contact = false;

  list<Obstacle*>::const_iterator obstacle;
  for (obstacle=m_obstacles.begin();
       obstacle!=m_obstacles.end() && !contact; obstacle++)
  {
    if ( voisin->isCompObstacle() )
      contact = voisin->isContact( *obstacle );
    else
      contact = (*obstacle)->isContact( voisin );
  }

  return contact;
}




// ----------------------------------------------------------------------------
// Les deux composants sont ils en contact ? Variante avec les rayons de VdW.
// A.WACHS - aout 2009 - Creation
// D.RAKOTONIRINA - janv 2014 - Modification
bool CompObstacle::isContactVdW( const Composant* voisin ) const
{
  bool contact = false;

  list<Obstacle*>::const_iterator obstacle;
  for (obstacle=m_obstacles.begin();
       obstacle!=m_obstacles.end() && !contact; obstacle++)
  {
    if ( voisin->isCompObstacle() )
      contact = voisin->isContactVdW( *obstacle );
    else
      contact = (*obstacle)->isContactVdW( voisin );
  }

  return contact;
}




// ----------------------------------------------------------------------------
// Les deux composants sont ils proches ?
// A.WACHS - aout 2009 - Creation
// D.RAKOTONIRINA - janv 2014 - Modification
bool CompObstacle::isProche( const Composant* voisin ) const
{
  bool contact = false;

  list<Obstacle*>::const_iterator obstacle;
  for (obstacle=m_obstacles.begin();
       obstacle!=m_obstacles.end() && !contact; obstacle++)
    if ( voisin->isCompObstacle() )
      contact = voisin->isProche( *obstacle );
    else
      contact = (*obstacle)->isProche( voisin );

  return contact;
}




// ----------------------------------------------------------------------------
// Les deux composants sont ils proches ? Variante avec les rayons de VdW.
// A.WACHS - aout 2009 - Creation
// D.RAKOTONIRINA - janv 2014 - Modification
bool CompObstacle::isProcheVdW( const Composant* voisin ) const
{
  bool contact = false;

  list<Obstacle*>::const_iterator obstacle;
  for (obstacle=m_obstacles.begin();
       obstacle!=m_obstacles.end() && !contact; obstacle++)
    if ( voisin->isCompObstacle() )
      contact = voisin->isProcheVdW( *obstacle );
    else
      contact = (*obstacle)->isProcheVdW( voisin );

  return contact;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Reload du groupe d'Obstacle & publication dans son referent
void CompObstacle::reload( Obstacle &obstacle, istream &file )
{
  string ttag;
  file >> ttag;
  while ( ttag != "</Composite>" )
  {
    Obstacle_BuilderFactory::reload( ttag, *this, file );
    file >> ttag;
  }
  EvalPosition();
  obstacle.append(this);
}




// ----------------------------------------------------------------------------
// Remise a 0 de la cinematique
// G.FERRER - Octo.2002 - Creation
void CompObstacle::resetCinematique()
{
  m_cinematique.reset();
  m_confinement.reset();

  list<Obstacle*>::iterator obstacle;
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); obstacle++)
    (*obstacle)->resetCinematique();
}




// ----------------------------------------------------------------------------
// Rotation de l'obstacle
// G.FERRER - Octo.2002 - Creation
void CompObstacle::Rotate( const Quaternion &rotation )
{
  list<Obstacle*>::iterator obstacle;
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); obstacle++)
    (*obstacle)->Rotate(rotation);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Suppression de l'obstacle dans la zone specifie.
// Supprime les obstacles contenus se trouvant dans la zone.
// G.FERRER - Fevr.2004 - Creation
void CompObstacle::Suppression( const BBox &box )
{
  list<Obstacle*>::iterator obstacle;
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); )
  {
    (*obstacle)->Suppression( box );
    if ( (*obstacle)->isIn( box ) ) obstacle = m_obstacles.erase(obstacle);
    else obstacle++;
  }
}




// ----------------------------------------------------------------------------
// Translation de l'obstacle
// G.FERRER - Octo.2002 - Creation
void CompObstacle::Translate( const Vecteur &translation )
{
  list<Obstacle*>::iterator obstacle;
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); obstacle++)
    (*obstacle)->Translate( translation );
}




// ----------------------------------------------------------------------------
// Ecriture de la position
// G.FERRER - Octo.2002 - Creation
void CompObstacle::writePosition( ostream &position )
{
  list<Obstacle*>::iterator obstacle;
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); obstacle++)
    (*obstacle)->writePosition( position );
}




// ----------------------------------------------------------------------------
// Ecriture du composant
// G.FERRER - Octo.2002 - Creation
// D. RAKOTONIRINA - Oct. 2014 - Modification
void CompObstacle::writeStatique( ostream &statique, Composant* composant )
{
  list<Obstacle*>::iterator obstacle;
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); obstacle++)
    (*obstacle)->writeStatique( statique );
}




// ----------------------------------------------------------------------------
// Sauvegarde de l'obstacle compose pour Reload
// G.FERRER - Octo.2002 - Creation
// D. RAKOTONIRINA - Oct. 2014 - Modification
void CompObstacle::write( ostream &fileSave, Composant const* composant ) const
{
  fileSave << "<Composite>\t" << m_nom << '\n';
  list<Obstacle*>::const_iterator obstacle;
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); obstacle++)
    (*obstacle)->write( fileSave );
  fileSave << "</Composite>\n";
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecriture de l'objet dans sa configuration courante pour
// post-processing avec GMV
void CompObstacle::GMVoutput( ostream &fileOut ) const
{
  list<Obstacle*>::const_iterator obstacle;
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); obstacle++)
    (*obstacle)->GMVoutput( fileOut );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Destruction de l'obstacle dans l'arborescence des obstacles
// A.WACHS - Mars.2010 - Creation
void CompObstacle::DestroyObstacle( const string &name_ )
{
  list<Obstacle*>::iterator obstacle;
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); )
  {
    if ( (*obstacle)->getName() == name_ ||
    	 (*obstacle)->getName() == "ToBeDestroyed" )
    {
      (*obstacle)->DestroyObstacle( "ToBeDestroyed" );
      delete *obstacle;
      obstacle = m_obstacles.erase(obstacle);
    }
    else
    {
      (*obstacle)->DestroyObstacle( name_ );
      obstacle++;
    }
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Supprime l'obstacle de nom name_
// A.WACHS - Mars.2010 - Creation
void CompObstacle::SupprimeObstacle( const string &name_, LinkedCell *LC )
{
  list<Obstacle*>::iterator obstacle;
  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end();obstacle++)
  {
    if ( (*obstacle)->getName() == name_ )
    {
      list<MonObstacle*> allObs = (*obstacle)->getObstacles();
      list<MonObstacle*>::iterator il;
      for(il=allObs.begin(); il!=allObs.end(); il++)
        (*il)->SupprimeObstacle( "ToBeErased", LC );
    }
    else
      (*obstacle)->SupprimeObstacle( name_, LC );
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Mise a jour de l'indicateur pour post-traitement Paraview
// A.WACHS - Fev.2011 - Creation
void CompObstacle::updateIndicator( Scalar temps, Scalar dt )
{
  list<Obstacle*>::iterator obstacle;

  if ( m_cinematique.rotationEnCours( temps, dt ) )
      getObstacles().front()->setIndicator( 1. );

  for (obstacle=m_obstacles.begin(); obstacle!=m_obstacles.end(); obstacle++)
    (*obstacle)->updateIndicator( temps, dt );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Cree et ajoute l'etat
void CompObstacle::createState( list<struct ObstacleState*> &obsStates )
{
  struct ObstacleState* obss = new ObstacleState;
  obss->nom = m_nom;
  obss->memento_config = new ConfigurationMemento();
  obss->memento_config->m_position = *m_geoFormeVdw->getTransform();
  obss->memento_cine = m_cinematique.createState();
  obsStates.push_back( obss );

  for (list<Obstacle*>::const_iterator obstacle=m_obstacles.begin();
  	obstacle!=m_obstacles.end();obstacle++)
    (*obstacle)->createState( obsStates );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Restauration de l'etat
void CompObstacle::restaureState( list<struct ObstacleState*>& obsStates )
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

  for (list<Obstacle*>::iterator obstacle=m_obstacles.begin();
  	obstacle!=m_obstacles.end();obstacle++)
    (*obstacle)->restaureState( obsStates );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Initialise le torseur des efforts sur le composant
// A.WACHS - Mai.2012 - Creation
void CompObstacle::InitializeForce( bool const& withWeight )
{
  for (list<Obstacle*>::iterator obstacle=m_obstacles.begin();
  	obstacle!=m_obstacles.end(); obstacle++)
    (*obstacle)->InitializeForce( false );
}




// ----------------------------------------------------------------------------
// Renvoie le torseur applique
// A.WACHS - Mai.2012 - Creation
Torseur const* CompObstacle::getTorseur()
{
  m_somme.setToBodyForce( *getPosition(), VecteurNul );

  for (list<Obstacle*>::iterator obstacle=m_obstacles.begin();
  	obstacle!=m_obstacles.end(); obstacle++)
    m_somme += *(*obstacle)->getTorseur();

  return &m_somme;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie le type de l'obstacle
string CompObstacle::getObstacleType()
{
  return m_ObstacleType;
}




// ----------------------------------------------------------------------------
// Nombre de points pour post-processing avec Paraview
// D. RAKOTONIRINA - Dec. 2014 - Creation
int CompObstacle::numberOfPoints_PARAVIEW() const
{
  cout << "Warning when calling CompObstacle::numberOfPoints_PARAVIEW() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Nombre de polygones elementaires pour post-processing avec Paraview
// D. RAKOTONIRINA - Dec. 2014 - Creation
list<Point> CompObstacle::get_polygonsPts_PARAVIEW( Vecteur const* translation )
const
{
  cout << "Warning when calling CompObstacle::get_polygonsPts_PARAVIEW() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Nombre de points pour post-processing avec Paraview
// D. RAKOTONIRINA - Dec. 2014 - Creation
int CompObstacle::numberOfCells_PARAVIEW() const
{
  cout << "Warning when calling CompObstacle::numberOfCells_PARAVIEW() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Ecrit les points du convexe pour post-processing avec Paraview
// D. RAKOTONIRINA - Dec. 2014 - Creation
void CompObstacle::write_polygonsPts_PARAVIEW( ostream &f,
  	Vecteur const* translation ) const
{
  cout << "Warning when calling CompObstacle::write_polygonsPts_PARAVIEW() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Ecrit les points du convexe pour post-processing avec Paraview
void CompObstacle::write_polygonsStr_PARAVIEW(list<int> &connectivity,
    	list<int> &offsets, list<int> &cellstype, int& firstpoint_globalnumber,
	int& last_offset) const
{
  cout << "Warning when calling CompObstacle::write_polygonsStr_PARAVIEW() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecriture de la position de la particule pour le Fluide
// D. RAKOTONIRINA - Avril. 2015 - Creation
void CompObstacle::writePositionInFluid( ostream &fileOut )
{
  cout << "Warning when calling CompObstacle::writePositionInFluid() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Initialize all contact map entries to false
void CompObstacle::setContactMapToFalse()
{
  cout << "Warning when calling CompObstacle::setContactMapToFalse() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Update contact map
void CompObstacle::updateContactMap()
{
  cout << "Warning when calling CompObstacle::updateContactMap() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Does the contact exist in the map, if yes return the pointer to the
// cumulative tangential displacement
bool CompObstacle::ContactInMapIsActive( double* &tangentialDepl, int const& id )
{
  cout << "Warning when calling CompObstacle::ContactInMapIsActive() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Add new contact in the map
void CompObstacle::addNewContactInMap( double const& tangentialDepl,
	int const& id )
{
  cout << "Warning when calling CompObstacle::addNewContactInMap() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Increase cumulative tangential displacement with component id
void CompObstacle::addDeplContactInMap( double const& tangentialDepl,
	int const& id )
{
  cout << "Warning when calling CompObstacle::addDeplContactInMap() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}

void CompObstacle::copyHistoryContacts( double* &destination, int start_index )
{
  cout << "Warning when calling CompObstacle::copyHistoryContacts() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}

// ----------------------------------------------------------------------------
// Copy existing contact in the map
void CompObstacle::copyContactInMap( std::tuple<int,int,int> const& id,
  bool const& isActive, int const& nbDtRemember, int const& nbCumulTangent,
  Vecteur const& tangent, Vecteur const& prev_normal,
  Vecteur const& cumulSpringTorque )
{
  cout << "Warning when calling CompObstacle::copyContactInMap() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}

int CompObstacle::getContactMapSize()
{
  cout << "Warning when calling CompObstacle::getContactMapSize() "
       << "\nShould not go into this class !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}
