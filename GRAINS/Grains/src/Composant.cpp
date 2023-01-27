#include "Composant.H"
#include "SaveTable.H"
#include "Particule.H"
#include "Obstacle.H"
#include "Torseur.H"
#include "Erreur.H"
#include "Memento.hh"
#include <algorithm>
using namespace std;


// ----------------------------------------------------------------------------
// Initialisation des attributs Static
int Composant::m_nb = 0;
// double Composant::m_heatCapacity = 0.;



// ----------------------------------------------------------------------------
// Constructeur par defaut
Composant::Composant( const bool &autonumbering ) :
  m_nomMateriau( "" ),
  m_masse( 0. ),
  m_geoFormeVdw( NULL ),
  m_vectIdParticle(10, -10),
  m_vectFmaxDist(10, 0.),
  m_vectKnElast(10, 0.),
  m_vectKtElast(10, 0.),
  m_vectInitialOverlap(10 ,0.),
  m_memento( NULL )
{
  if ( autonumbering )
  {
    m_id = Composant::m_nb;
    Composant::m_nb++;
  }
  else m_id = -1;
}




// ----------------------------------------------------------------------------
// Constructeur par copie
Composant::Composant( const Composant &copie ) :
  m_id( Composant::m_nb ),
  m_nomMateriau( copie.m_nomMateriau ),
  m_masse( copie.m_masse ),
  m_vectIdParticle(10, -10),
  m_vectFmaxDist(10, 0.),
  m_vectKnElast(10, 0.),
  m_vectKtElast(10, 0.),
  m_vectInitialOverlap(10, 0.),
  m_memento( NULL )
{
  Composant::m_nb++;

  m_geoFormeVdw = new FormeVdW( *copie.m_geoFormeVdw );
}




// ----------------------------------------------------------------------------
// Destructeur
Composant::~Composant()
{
  Composant::m_nb--;
  delete m_geoFormeVdw;
  if ( m_memento ) delete m_memento;
}




// ----------------------------------------------------------------------------
// Ajout d'une force au torseur des efforts
void Composant::addForce( const Point &point, const Vecteur &force )
{
  m_somme.addForce( point, force );
}




// ----------------------------------------------------------------------------
// Force at the contact point for the stress tensor
// D. RAKOTONIRINA - Fev. 2017 - Creation
void Composant::computeInternalMoments( const Point &point, const Vecteur &force )
{
  cout << "WARNING!!! Composant::computeInternalMoments() "
       << "is no yet implemented\n"
       << "Need for an assistance! Stop running!" << endl;
  exit(10);
}




// ----------------------------------------------------------------------------
// Internal moment of particles
// D. RAKOTONIRINA - Fev. 2017 - Creation
vector<Scalar> const* Composant::getStressTensor( )
{
  cout << "WARNING!!! Composant::getInternalMoment() "
       << "is no yet implemented\n"
       << "Need for an assistance! Stop running!" << endl;
  exit(10);
  return( NULL );
}




// ----------------------------------------------------------------------------
// Internal moment divided by the volume of the system for the stress tensor
// D. RAKOTONIRINA - Fev. 2017 - Creation
vector<Scalar> const* Composant::getInternalMoment( )
{
  cout << "WARNING!!! Composant::getInternalMoment() "
       << "is no yet implemented\n"
       << "Need for an assistance! Stop running!" << endl;
  exit(10);
  return( NULL );
}




// ----------------------------------------------------------------------------
// Add contact force on each particle for postprocessing purposes
void Composant::addContactForcePP( const Vecteur &force )
{
  Vecteur absforce;
  absforce[X] = fabs(force[X]);
  absforce[Y] = fabs(force[Y]);
  absforce[Z] = fabs(force[Z]);
  m_ForceContactPP += absforce;
}




// ----------------------------------------------------------------------------
// Add contact force on each particle for postprocessing purposes
void Composant::addContactForcePP_instantaneous( const Vecteur &force )
{
  m_ForceContactPP_instantaneous += force;
}





// ----------------------------------------------------------------------------
// Add lubrication force on each particle for postprocessing purposes
void Composant::addLubriForcePP( const Vecteur &force )
{
  Vecteur absforce;
  absforce[X] = fabs(force[X]);
  absforce[Y] = fabs(force[Y]);
  absforce[Z] = fabs(force[Z]);
  m_ForceLubriPP += absforce;
}




// ----------------------------------------------------------------------------
// Ajout d'une force s'exercant au centre de gravite au torseur des efforts
void Composant::addBodyForce( const Vecteur &force )
{
  m_somme.addForce( force );
}




// ----------------------------------------------------------------------------
// Ajouter un moment
void Composant::addMoment( const Vecteur &moment )
{
  m_somme.addMoment( moment );
}




// ----------------------------------------------------------------------------
// Determination de la "Bouding Box" du composant
BBox Composant::Boite() const
{
  return m_geoFormeVdw->BoxForme();
}




// ----------------------------------------------------------------------------
// Acces a la forme du composant.
// G.FERRER - Octo.2000 - Creation
const FormeVdW* Composant::getForme() const
{
  return ( m_geoFormeVdw );
}




// ----------------------------------------------------------------------------
// Acces a la forme du composant.
// AWACHS - Janv.2011 - Creation
FormeVdW* Composant::getForme()
{
  return ( m_geoFormeVdw );
}




// ----------------------------------------------------------------------------
// Acces a l'ID
// G.FERRER - Nove.2000 - Creation
int Composant::getID() const
{
  return ( m_id );
}




// ----------------------------------------------------------------------------
// Masse du Composant.
Scalar Composant::getMasse() const
{
  return m_masse;
}




// ----------------------------------------------------------------------------
// Renvoi la valeur du point position de la particule.
Point const* Composant::getPosition() const
{
  return m_geoFormeVdw->getCentre();
}




// ----------------------------------------------------------------------------
// Renvoi le vecteur position de la particule dans le pointeur pos.
void Composant::getPosition( Scalar *pos ) const
{
  Point const* pot = m_geoFormeVdw->getCentre();
  for (int i=0; i<3; i++) pos[i] = (*pot)[i];
}




// ----------------------------------------------------------------------------
// Rayon du composant
// F.PRADEL - Fevr.2000 - Creation
Scalar Composant::getRayon() const
{
  return ( m_geoFormeVdw->getRayon() );
}




// ----------------------------------------------------------------------------
// retourne l'information si la particule est retrecissante
// M. SULAIMAN- Nov.2015 - Creation
int Composant::getShrinkingMode()
{
  return m_geoFormeVdw->getShrinkingMode();
}




// ----------------------------------------------------------------------------
// Rayon courant de la particule
// M. SULAIMAN- Nov.2015 - Creation
void Composant::set_shrinking_radius(Scalar R)
{
  m_geoFormeVdw->set_shrinking_radius(R);
}




// ----------------------------------------------------------------------------
// Rayon de la sphere de meme volume que le composant
// A. WACHS- Sept.2014 - Creation
Scalar Composant::getRayonSphereEquivalente() const
{
  return ( 0. );
}




// ----------------------------------------------------------------------------
// Rayon d'interaction du composant
// A.WACHS - Fevr.2009 - Creation
Scalar Composant::getRayonInteraction() const
{
  return ( m_geoFormeVdw->getRayonInterAction() );
}




// ----------------------------------------------------------------------------
// Volume du composant
Scalar Composant::getVolume() const
{
  return m_geoFormeVdw->getVolume();
}




// ----------------------------------------------------------------------------
// Copy le vecteur position de la particule dans le pointeur pos
// en d�butant � la position i
void Composant::copyPosition( double *pos, int i ) const
{
  Point const* pot = m_geoFormeVdw->getCentre();
  for (int j=0 ;j<3; j++) pos[i+j] = (*pot)[j];
}





// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Copie de la transformation dans le vecteur vit
// en d�butant � la position i
void Composant::copyTransform( double *vit, int i ) const
{
  m_geoFormeVdw->copyTransform( vit, i );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Copie de la transformation dans le vecteur vit
// en d�butant � la position i, avec une translation du centre de gravite
void Composant::copyTransform( double *vit, int i, Vecteur const& vec ) const
{
  m_geoFormeVdw->copyTransform( vit, i, vec );
}




// ----------------------------------------------------------------------------
// Les deux composants sont ils en contact ?
bool Composant::isContact( const Composant* voisin ) const
{
  if ( this != voisin )
    return m_geoFormeVdw->Forme::isContact( *voisin->m_geoFormeVdw );
  else return false;
}




// ----------------------------------------------------------------------------
// Les deux composants sont ils en contact ? Variante avec les rayons de VdW
bool Composant::isContactVdW( const Composant* voisin ) const
{
  if ( this != voisin )
    return m_geoFormeVdw->isContact( *voisin->m_geoFormeVdw );
  else return false;
}




// ----------------------------------------------------------------------------
// Les deux composants sont ils proches ?
bool Composant::isProche( const Composant* voisin ) const
{
  if ( this != voisin )
    return m_geoFormeVdw->Forme::isProche( *voisin->m_geoFormeVdw );
  else return false;
}




// ----------------------------------------------------------------------------
// Les deux composants sont ils proches ? Variante avec les rayons de VdW
bool Composant::isProcheVdW( const Composant* voisin ) const
{
  if ( this != voisin )
    return m_geoFormeVdw->isProche( *voisin->m_geoFormeVdw );
  else return false;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Le Composant est il dans la boite ?
bool Composant::isIn( const BBox &boite ) const
{
  const Point   &origin = boite.getCenter();
  const Vecteur &extent = boite.getExtent();

  Point const* centre = m_geoFormeVdw->getCentre();
  Vecteur dist = *centre - origin;

  bool status =
    fabs( dist[X] ) <= extent[X] &&
    fabs( dist[Y] ) <= extent[Y] &&
    fabs( dist[Z] ) <= extent[Z];

  return status;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Nom du materiau
const string& Composant::materiau() const
{
  return m_nomMateriau;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Definit le type de materiau
void Composant::setMateriau( const string& mat )
{
  m_nomMateriau = mat ;
}




// ----------------------------------------------------------------------------
// Rotation de la particule
void Composant::Rotate( const Quaternion &rotation )
{
  m_geoFormeVdw->Rotate( rotation );
}




// ----------------------------------------------------------------------------
// Translation de la particule
// D.PETIT - Aout.2000 - Creation
void Composant::Translate( const Vecteur &translation )
{
  *m_geoFormeVdw += translation;
}




// ----------------------------------------------------------------------------
// Positionne un composant
void Composant::setPosition( const Scalar *pos )
{
  m_geoFormeVdw->setPosition( pos );
}




// ----------------------------------------------------------------------------
// Positionne un composant
void Composant::setPosition( const Point &centre )
{
  m_geoFormeVdw->setOrigin( (Scalar *) &centre );
}




/*
// ----------------------------------------------------------------------------
// Set solid-body heat capacity
void Composant::set_heatCapacity( const Scalar &heatCapacity_ )
{
  m_heatCapacity = heatCapacity_;
}
*/



// ----------------------------------------------------------------------------
// Applique la transformation trot au composant
void Composant::composePosition( const Transform &trot )
{
  m_geoFormeVdw->composeTransform( trot );
}





// ----------------------------------------------------------------------------
// Comparaison entre deux composants
// G.FERRER - Octo.2000 - Creation
bool Composant::operator == ( const Composant &composant ) const
{
  return ( this == &composant );
}




// ----------------------------------------------------------------------------
// Comparaison entre deux composants
// G.FERRER - Octo.2000 - Creation
bool Composant::operator != ( const Composant &composant ) const
{
  return ( this != &composant );
}





// ----------------------------------------------------------------------------
// Ecriture de l'information de position
// G. FERRER - Janv.2001 - Creation
// D. RAKOTONIRINA - Sept. 2014 - Modification
void Composant::writePosition( ostream &position, Composant* composant ) const
{
  Scalar buf = 0. ;
  position << "*Couleur\n"
	   << buf << '\t'
	   << buf << '\t'
	   << buf << '\n';
  if ( composant ) composant->getForme()->writePosition( position );
  else m_geoFormeVdw->writePosition( position );
}




// ----------------------------------------------------------------------------
// Ecriture de l'information statique
// G.FERRER - Janv.2001 - Creation
// D. RAKOTONIRINA - Sept. 2014 - Modification
void Composant::writeStatique( ostream &statique,
	  Composant const* composant ) const
{
  statique << this << '\t'
	   << m_id   << '\n';
  statique << "*Materiau\n"
	   << m_nomMateriau << '\n';
  m_geoFormeVdw->writeStatique( statique, composant );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecriture de l'objet dans sa configuration courante pour
// post-processing avec GMV
void Composant::GMVoutput( ostream &fileOut ) const
{
  m_geoFormeVdw->GMVoutput( fileOut );
}




// ----------------------------------------------------------------------------
// Renvoie la force appliquee
Vecteur const* Composant::getForce() const
{
  return m_somme.getForce();
}




// ----------------------------------------------------------------------------
// Renvoie le moment applique
Vecteur const* Composant::getMoment() const
{
  return m_somme.getMoment();
}




// ----------------------------------------------------------------------------
// Renvoie le torseur applique
Torseur const* Composant::getTorseur()
{
  return &m_somme;
}




// ----------------------------------------------------------------------------
// Initialise the PostProcessing Contact Force
void Composant::InitializePostProcessingForce()
{
  m_ForceContactPP = VecteurNul;
  m_ForceContactPP_instantaneous = VecteurNul;
  m_ForceLubriPP= VecteurNul;
  m_ForceContactPP_instantaneous = VecteurNul;
}




// ----------------------------------------------------------------------------
// Integrate the PostProcessing contact force in time
void Composant::IntegrateContactPostProcessingForce(Scalar nb)
{
  m_ForceContactPP /= nb;
  m_ForceLubriPP /= nb;
}




// ----------------------------------------------------------------------------
// Initialise force at the contact point for the post-processing
// of the stress tensor
// D. RAKOTONIRINA - Fev. 2017 - Creation
void Composant::InitializeForceAtContactPoint()
{
  cout << "WARNING!!! Composant::InitializeForceAtContactPoint()"
       << "is no yet implemented\n"
       << "Need for an assistance! Stop running!" << endl;
  exit(10);
}




// ----------------------------------------------------------------------------
// Initialise le torseur des efforts sur le composant
void Composant::InitializeForce( bool const& withWeight )
{
  m_somme.setToBodyForce( *m_geoFormeVdw->getCentre(), VecteurNul );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Sauvegarde de l'etat
void Composant::saveConfigState()
{
  if ( !m_memento ) m_memento = new ConfigurationMemento();
  m_memento->m_position = *m_geoFormeVdw->getTransform();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Cree et renvoie l'etat
ConfigurationMemento* Composant::createConfigState()
{
  ConfigurationMemento* Pmemento_ = new ConfigurationMemento();
  Pmemento_->m_position = *m_geoFormeVdw->getTransform();

  return Pmemento_;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Determine s'il y a intersection entre une Forme et un Composant
// Vrai s'il y a intersection
// D. RAKOTONIRINA - Mars 2014 - Creation
bool Composant::intersectWithForme( const Forme &b )
{
  Vecteur v = *(b.getTransform()->getOrigin()) - *m_geoFormeVdw->getCentre();
  if ( Norm(v) < EPSILON ) return true;
  return intersect( *m_geoFormeVdw->getConvex(), *b.getConvex(),
  	*m_geoFormeVdw->getTransform(), *b.getTransform(), v );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie this pour une particule convex et m_masterComposite pour une
// particule elemenatire de CompParticule
// D. RAKOTONIRINA - Sept 2014 - Creation
Composant* Composant::ReferenceComposant()
{
  return ( this );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// GET vecteur vitesse de translation du fluide interpole
// au centre de gravite de la particule
// A. Esteghamatian - Janvier 2015 - Creation
Vecteur const* Composant::getVitesseTranslation_fluide() const
{
  cout << "WARNING!!!!Composant::getVitesseTranslation_fluide()"
  "This should not be called here, refer to Particule.cpp" << endl;
  return ( NULL );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie le nom associe a la particule composite
// D. RAKOTONIRINA - Sept 2014 - Creation
string Composant::getPartName() const
{
  cout << "WARNING!!! Composant::getPartName()"
       << "is no yet implemented\n"
       << "Need for an assistance! Stop running!" << endl;
  exit(10);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ajoute un nombre de contacts au nombre de contacts de la particule;
//  Utilisation: ajoute a la particule de reference periodique les contacts
//  de son clone periodique
//  D. RAKOTONIRINA - Nov 2014 - Creation
void Composant::addToCoordinationNumber( int const& nc )
{
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Only in case of demcfd temperature coupling for the moment, add a heatFlux
// to component's heatFlux sum
void Composant::add_heatFlux( const Scalar &heatFlux_ )
{
  m_sum_HeatFlux += heatFlux_ ;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Returns maxdist
Scalar Composant::get_FmaxDist(int id_particle) const
{
  size_t size = m_vectIdParticle.size();
  Scalar res = 0.;
  for( size_t i=0; i<size; i++ )
    if( m_vectIdParticle[i] == id_particle )
      res =  m_vectFmaxDist[i];
  return res;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Returns normal cohesive elastic coef
Scalar Composant::get_KnElast(int id_particle) const
{
  size_t size = m_vectIdParticle.size();
  Scalar res = 0.;
  for( size_t i=0; i<size; i++ )
    if( m_vectIdParticle[i] == id_particle )
      res =  m_vectKnElast[i];
  return res;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Returns tangential cohesive elastic coef
Scalar Composant::get_KtElast(int id_particle) const
{
  size_t size = m_vectIdParticle.size();
  Scalar res = 0.;
  for( size_t i=0; i<size; i++ )
    if( m_vectIdParticle[i] == id_particle )
      res =  m_vectKtElast[i];
  return res;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Initialization of Fmaxdist for particles in contact
void Composant::initializeCohesiveProperties( int id_particle, Scalar distance,
    Scalar kn_elast, Scalar kt_elast, Scalar overlap )
{
  size_t size = m_vectIdParticle.size();
  for( size_t i=0; i<size; i++ )
    if( m_vectIdParticle[i]==-10 && id_particle!=m_vectIdParticle[i] )
    {
      m_vectIdParticle[i] = id_particle;
      m_vectFmaxDist[i] = distance;
      m_vectKnElast[i] = kn_elast;
      m_vectKtElast[i] = kt_elast;
      m_vectInitialOverlap[i] = overlap;
      break;
    }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Modify the cohesive properties of a chosen particle ( the Fmax distance and
// all K coefficient )
void Composant::set_CohesiveProperties( int id_particle, Scalar distance,
    Scalar kn_elast, Scalar kt_elast )
{
  size_t size = m_vectIdParticle.size();
  for( size_t i=0; i<size; i++ )
    if( m_vectIdParticle[i] == id_particle )
    {
      m_vectFmaxDist[i] = distance;
      m_vectKnElast[i] = kn_elast;
      m_vectKtElast[i] = kt_elast;
   }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Copy the list of particles in cohesion with
void Composant::copy_VectIdParticle( double* id_part, int i) const
{
  size_t size = m_vectIdParticle.size();
  for( size_t j=0; j<size; j++ )
    id_part[i+j] = (double)m_vectIdParticle[j];
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Copy the vector containing maxdist corresponding to IDvector in cohesion with
void Composant::copy_VectFmaxDist( double* fmax_dist, int i) const
{
  size_t size = m_vectFmaxDist.size();
  for( size_t j=0; j<size; j++ )
    fmax_dist[i+j] = m_vectFmaxDist[j];
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Copy the vector containing K coef corresponding to IDvector in cohesion with
void Composant::copy_VectKnElast( double* kn_elast, int i) const
{
  size_t size = m_vectKnElast.size();
  for( size_t j=0; j<size; j++ )
    kn_elast[i+j] = m_vectKnElast[j];
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Copy the vector containing K coef corresponding to IDvector in cohesion with
void Composant::copy_VectKtElast( double* kt_elast, int i) const
{
  size_t size = m_vectKtElast.size();
  for( size_t j=0; j<size; j++ )
    kt_elast[i+j] = m_vectKtElast[j];
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Set the list of particles in cohesion with
void Composant::set_VectIdParticle( double id_part, int i)
{
  size_t size = m_vectIdParticle.size();
  for( size_t j=0; j<size; j++ )
    if( j == size_t(i) )
      m_vectIdParticle[j] = int(id_part);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Set the vector containing maxdist corresponding to
// ID vector in cohesion with
void Composant::set_VectFmaxDist( double fmax_dist, int i)
{
  size_t size = m_vectFmaxDist.size();
  for( size_t j=0 ;j<size; j++ )
    if( j == size_t(i) )
      m_vectFmaxDist[j] = fmax_dist;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Set the vector containing Ke coef  corresponding to
// ID vector in cohesion with
void Composant::set_VectKnElast( double kn_elast, int i)
{
  size_t size = m_vectKnElast.size();
  for (size_t j=0 ;j<size; j++)
    if( j == size_t(i) )
      m_vectKnElast[j] = kn_elast;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Set the vector containing Ke coef  corresponding to
// ID vector in cohesion with
void Composant::set_VectKtElast( double kt_elast, int i)
{
  size_t size = m_vectKtElast.size();
  for( size_t j=0; j<size; j++ )
    if( j == size_t(i) )
      m_vectKtElast[j] = kt_elast;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Copy the vector containing initial overlap corresponding to
// ID vector in cohesion with
void Composant::copy_VectInitialOverlap( double* initialOverlap, int i) const
{
  size_t size = m_vectInitialOverlap.size();
  for( size_t j=0; j<size; j++ )
    initialOverlap[i+j] = m_vectInitialOverlap[j];
}





// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Set the vector containing initial overlap coef  corresponding to
// ID vector in cohesion with
void Composant::set_VectInitialOverlap( double Kelast, int i)
{
  size_t size = m_vectInitialOverlap.size();
  for (size_t j=0; j<size; j++)
    if( j == size_t(i) )
      m_vectInitialOverlap[j] = Kelast;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Returns cohesive elastic coef
Scalar Composant::get_InitialOverlap(int id_particle) const
{
  size_t size = m_vectIdParticle.size();
  Scalar res = 0.;
  for( size_t i=0; i<size; i++)
    if( m_vectIdParticle[i] == id_particle )
      res = m_vectInitialOverlap[i];
  return res;
}




// ----------------------------------------------------------------------------
// Initialize all contact map entries to false
void Composant::setContactMapToFalse()
{
  map<int,std::tuple<bool, int, int, Vecteur, double, Vecteur>,bool >::iterator it;

  for (it=m_contactMap.begin();it!=m_contactMap.end();++it)
    get<0>(it->second) = false;
}




// ----------------------------------------------------------------------------
// Update contact map
void Composant::updateContactMap()
{
  map<int,std::tuple<bool, int, int, Vecteur, double, Vecteur>>::iterator it;
  list<int> keywithfalse;
  list<int>::const_iterator il;

  for (it=m_contactMap.begin();it!=m_contactMap.end();++it)
    if ( !(get<0>(it->second)) ) keywithfalse.push_back(it->first);

  for (il=keywithfalse.begin();il!=keywithfalse.end();il++)
  {
    if ( get<1>(m_contactMap[*il])==0 )
    {
      m_contactMap.erase(*il);
    }
    else
    {
      get<1>(m_contactMap[*il]) --;
    }
  }
}



// ----------------------------------------------------------------------------
// Does the contact exist in the map, if yes return the pointer to the
// cumulative tangential displacement
bool Composant::ContactInMapIsActive( int const& id, int* &nbCumulTangent,
  Vecteur* &tangent, double* &tangentialDepl, Vecteur* &cumulSpringTorque )
{
  bool active = false;
  map<int,std::tuple<bool, int, int, Vecteur, double, Vecteur> >
    ::iterator it = m_contactMap.find(id);

  if ( it != m_contactMap.end() )
  {
    active = true;
    get<0>(it->second) = true;
    get<1>(it->second) = NB_STEPS_REMEMBER_ENDED_CONTACT;
    nbCumulTangent = &(get<2>(it->second));
    tangent = &(get<3>(it->second));
    tangentialDepl = &(get<4>(it->second));
    cumulSpringTorque = &(get<5>(it->second));
  }
  else{
    nbCumulTangent = NULL;
    tangent = NULL;
    tangentialDepl = NULL;
    cumulSpringTorque = NULL;
  }
  return active;
}




// ----------------------------------------------------------------------------
// Add new contact in the map
void Composant::addNewContactInMap( int const& id, int const& nbCumulTangent,
  Vecteur const& tangent, double const& tangentialDepl,
  Vecteur const& cumulSpringTorque )
{
  // pair<double,bool> pp (tangentialDepl,true);
  // pair<int,pair<double,bool> > ppp (id,pp);
 // cout << "Create (" << m_id << "," << id << ")" << endl;
  m_contactMap.insert(std::make_pair( id, std::make_tuple( true,
    NB_STEPS_REMEMBER_ENDED_CONTACT,
    nbCumulTangent, tangent, tangentialDepl, cumulSpringTorque) ));
}




// ----------------------------------------------------------------------------
// Increase cumulative tangential displacement with component id
void Composant::addDeplContactInMap( int const& id, int const& nbCumulTangent,
  Vecteur const& tangent, double const& tangentialDepl,
  Vecteur const& cumulSpringTorque )
{
  get<0>(m_contactMap[id]) = true;
  get<1>(m_contactMap[id]) = NB_STEPS_REMEMBER_ENDED_CONTACT;
  get<2>(m_contactMap[id]) = nbCumulTangent;
  get<3>(m_contactMap[id]) = tangent;
  get<4>(m_contactMap[id]) += tangentialDepl;
  get<5>(m_contactMap[id]) = cumulSpringTorque;
}

// ---------------------------------------------------------------------------
// Print active neighbors of the particle
void Composant::printActiveNeighbors(int const& id )
{
    map<int,std::tuple<bool, int, int, Vecteur, double, Vecteur> >
      ::iterator it;
    if (m_contactMap.begin() != m_contactMap.end())
    {
        cout << "Neighbors of #" << id << ": ";
        for (it=m_contactMap.begin();it!=m_contactMap.end();++it){
            if (get<0>(it->second)){
                cout << it->first  << "  ";
            }
        }
        cout << endl;
    }
}
