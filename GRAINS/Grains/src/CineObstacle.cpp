// Gestion de la cinematique pour des obstacles
// G.FERRER - Juin.2000 - Creation
// ============================================================================
#include "CineObstacle.H"
#include "MonObstacle.H"
#include "Torseur.H"
#include "Memento.hh"


// ----------------------------------------------------------------------------
// Constructeur par defaut
// G.FERRER - Juin.2000 - Creation
CineObstacle::CineObstacle():
  m_memento( NULL )
{
}




// ----------------------------------------------------------------------------
// Destructeur
// G.FERRER - Juin.2000 - Creation
CineObstacle::~CineObstacle()
{
  if ( m_memento ) delete m_memento;
  clearAndDestroy();
}




// ----------------------------------------------------------------------------
// Ajout du chargement impose a la cinematique
// G.FERRER - Fevr.2002 - Creation
void CineObstacle::append( ObstacleChargement &chargement )
{
  if ( m_chargements.empty( ))
    m_chargements.push_back(&chargement);
  else 
  {
    list<ObstacleChargement*>::iterator c;
    for (c=m_chargements.begin(); c!=m_chargements.end() && **c<chargement; c++)
	{}
    m_chargements.insert(c, &chargement);
  }
}




// ----------------------------------------------------------------------------
// Destruction de l'ensemble des m_chargements d'obstacle
// G.FERRER - Juin.2003 - Creation
void CineObstacle::clearAndDestroy()
{
  list<ObstacleChargement*>::iterator chargement;

  for (chargement=m_chargements.begin(); 
       chargement!=m_chargements.end(); chargement++) 
    delete *chargement;

  m_chargements.clear();
}




// ----------------------------------------------------------------------------
// Decomposition de la cinematique de l'Obstacle avec une autre cinematique
// Si l'autre cinématique est nulle, soit si le composite dont l'obstacle 
// fait partie n'a pas de mouvement imposé, cette méthode initialise la
// cinématique de l'objet à nulle
// G.FERRER - Octo.2002 - Creation
void CineObstacle::Decompose( const CineObstacle &voisine,
	const Vecteur &levier )
{
  Matrix mat( voisine.m_QuaternionRotationOverDt );
  Vecteur rota = voisine.m_rotationOverTimeStep;
  Scalar  d = Norm(rota);

  if ( d != 0. ) 
  {
    Vecteur vect = ( sin( d / 2.) / d ) * rota;
    m_QuaternionRotationOverDt = Quaternion( vect, cos( d / 2. ) );
  }
  else 
    m_QuaternionRotationOverDt = Quaternion( 0., 0., 0., 1. );
  
  // Pour des composites de composites de ... etc, vérifier
  // que l'écriture correcte ne serait pas:
  // m_rotationOverTimeStep += voisine.m_rotationOverTimeStep;
  // pour conserver l'aspect récursif
  // A discuter avec Gilles
  m_rotationOverTimeStep = voisine.m_rotationOverTimeStep;
  m_translationOverTimeStep = voisine.m_translationOverTimeStep 
  	+ ( ( mat * levier ) - levier );
	
  m_vitesseRotation += voisine.m_vitesseRotation;
  m_vitesseTranslation += voisine.m_vitesseTranslation 
  	+ ( voisine.m_vitesseRotation ^ levier );
}




// ----------------------------------------------------------------------------
// Evaluation du deplacement en fonction des m_chargements 
// dans l'intervalle de temps specifie et du chargement du composite dont
// l'obstacle fait partie
// G.FERRER - Fevr.2002 - Creation
bool CineObstacle::Deplacement( Scalar temps, Scalar dt )
{
  Vecteur depl, rota, vt, vr;

  // Chargements de l'obstacle
  list<ObstacleChargement*>::iterator chargement;
  for (chargement=m_chargements.begin(); chargement!=m_chargements.end();
  	chargement++) 
    if ( (*chargement)->isActif(temps, dt) ) 
    {
      depl += (*chargement)->DeplacementTranslation( temps, dt ); 
      rota += (*chargement)->DeplacementRotation( temps, dt );
      vt += *(*chargement)->VitesseTranslation( temps, dt );
      vr += *(*chargement)->VitesseRotation( temps, dt );      
    } 

  Scalar d = Norm( rota );
  if ( d != 0. ) 
  {
    Vecteur vect = ( sin( d /2. ) / d ) * rota;
    m_QuaternionRotationOverDt = Quaternion( vect, cos( d / 2. ) );
  }

  // Sur le pas de temps, le mouvement de l'obstacle correspond à celui du
  // composite dont il fait partie plus son mouvement propre
  // Si le composite dont il fait partie n'a pas de mouvement imposé, 
  // m_translationOverTimeStep et m_rotationOverTimeStep sont mis à 0 par
  // la méthode CineObstacle::Decompose et seuls les chargements propres de
  // l'obstacle sont pris en compte, ce qui justifie l'utilisation du +=  
  m_translationOverTimeStep += depl;
  m_rotationOverTimeStep += rota;
  
  m_vitesseRotation += vr;
  m_vitesseTranslation += vt; 

  return ( Norm(depl) != 0. || Norm(rota) != 0. );
}




// ----------------------------------------------------------------------------
// Valeur du quaternion de rotation sur le pas de temps
// G.FERRER - Octo.2002 - Creation
// A.WACHS - Octo.2010 - Modif
Quaternion const* CineObstacle::getQuaternionRotationOverDt() const
{
  return &m_QuaternionRotationOverDt;
}




// ----------------------------------------------------------------------------
// Valeur de la translation pendant un pas de temps
// G.FERRER - Octo.2002 - Creation
Vecteur const* CineObstacle::getTranslation() const
{
  return &m_translationOverTimeStep;
}




// ----------------------------------------------------------------------------
// Vitesse de rotation de l'obstacle
// G.FERRER - Janv.2002 - Creation
Vecteur const* CineObstacle::getVitesseRotation() const
{
  return &m_vitesseRotation;
}




// ----------------------------------------------------------------------------
// Vitesse de translation de l'obstacle
// G.FERRER - Juin.2000 - Creation
Vecteur const* CineObstacle::getVitesseTranslation() const
{
  return &m_vitesseTranslation;
}




// ----------------------------------------------------------------------------
// Mise a zero de la cinematique
// G.FERRER - Juin.2000 - Creation
void CineObstacle::reset()
{
  m_translationOverTimeStep = 0.;
  m_rotationOverTimeStep = 0.;
  m_vitesseTranslation = 0.;
  m_vitesseRotation = 0.;
}




// ----------------------------------------------------------------------------
// Affectation d'une valeur des vitesses initiale a la cinematique
// G.FERRER - Octo.2002 - Creation
void CineObstacle::set( CineObstacle &cinematique )
{
  m_translationOverTimeStep = cinematique.m_translationOverTimeStep;
  m_rotationOverTimeStep = cinematique.m_rotationOverTimeStep;
  m_vitesseTranslation = cinematique.m_vitesseTranslation;
  m_vitesseRotation = cinematique.m_vitesseRotation;
}




// ----------------------------------------------------------------------------
// Affectation de la vitesse de l'obstacle
// A.WACHS - Fev.2012 - Creation
void CineObstacle::setVitesse( Vecteur const* vitesseTranslation,
  	Vecteur const* vitesseRotation )
{
  m_vitesseTranslation = *vitesseTranslation ;
  m_vitesseRotation = *vitesseRotation ;
}   




// ----------------------------------------------------------------------------
// Vitesse relative de l'obstacle
// D.PETIT  - Aout.2000 - Creation
// G.FERRER - Octo.2002 - Prise en compte de la rotation
Vecteur CineObstacle::Vitesse( const Vecteur &om ) const
{
  return ( m_vitesseTranslation + (m_vitesseRotation ^ om) );
}




// ----------------------------------------------------------------------------
// Operateur d'ecriture
// G.FERRER - Juin.2000 - Creation
ostream &operator << ( ostream &fileOut, 
	CineObstacle &cinematique )
{
  fileOut << "*CineObstacle\n";
  fileOut << cinematique.m_translationOverTimeStep;
   fileOut << "*Chargement\n";
  fileOut << cinematique.m_chargements.size() << '\n';
  list<ObstacleChargement*>::iterator chargement;
  for (chargement=cinematique.m_chargements.begin(); 
       chargement!=cinematique.m_chargements.end(); chargement++)
    fileOut << **chargement;

 return ( fileOut );
}




// ----------------------------------------------------------------------------
// Operateur de lecture
// G.FERRER - Jun.2000 - Creation
istream &operator >> ( istream &fileIn, 
	CineObstacle &cinematique )
{
  string cle;
  fileIn >> cle;
  fileIn >> cinematique.m_translationOverTimeStep;

  return ( fileIn );
}




// ----------------------------------------------------------------------------
// L'obstacle a t il un chargement en rotation actif
// A.WACHS - Fev.2011 - Creation
bool CineObstacle::rotationEnCours( Scalar temps, Scalar dt ) const
{
  bool rotation = false ;
  list<ObstacleChargement*>::const_iterator chargement;
  for (chargement=m_chargements.begin(); chargement!=m_chargements.end() 
  	&& !rotation; chargement++)
    if ( (*chargement)->isActif( temps, dt ) &&
    	( (*chargement)->getType() == "Rotation" ||
	(*chargement)->getType() == "RotationSinusoidale" ) ) rotation = true;
  
  return rotation; 
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Sauvegarde de l'etat
void CineObstacle::saveState()
{
  if (!m_memento) m_memento = new CineObstacleMemento();

  m_memento->m_vitesseTranslation = m_vitesseTranslation;
  m_memento->m_vitesseRotation = m_vitesseRotation;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Cree et renvoie l'etat
CineObstacleMemento* CineObstacle::createState()
{
  CineObstacleMemento* memento_ = new CineObstacleMemento();
  
  memento_->m_vitesseTranslation = m_vitesseTranslation;
  memento_->m_vitesseRotation = m_vitesseRotation;
  
  return memento_;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Restauration de l'etat
void CineObstacle::restaureState()
{
  m_vitesseTranslation = m_memento->m_vitesseTranslation;
  m_vitesseRotation = m_memento->m_vitesseRotation;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Restauration de l'etat
void CineObstacle::restaureState( CineObstacleMemento const* memento_ )
{
  m_vitesseTranslation = memento_->m_vitesseTranslation;
  m_vitesseRotation = memento_->m_vitesseRotation;
} 




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Acces au chargement
list<ObstacleChargement*> CineObstacle::getChargements() const
{
  return m_chargements;
}
