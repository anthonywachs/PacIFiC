#include "ParticulePeriodique.hh"
#include "Memento.hh"


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur avec arguments
// !!! utilise en SEQUENTIEL !!!
ParticulePeriodique::ParticulePeriodique( Particule* reference,
	Particule const* particuleClasseReference,
  	ObstaclePeriodique const* contactReferenceObstacle,
	Vecteur const& trans,
	const int &nbper ) :
  Particule( -2, particuleClasseReference, 
	*(reference->getVitesseTranslation()),
	*(reference->getCinematique()->getRotation()),		
	*(reference->getVitesseRotation()),		 
	*(reference->getForme()->getTransform()),
	COMPUTE, 
	0 ),
  m_reference( reference ),
  m_referenceID( reference->getID() ),
  m_contactReferenceObstacle( contactReferenceObstacle ),
  m_nperiodes( nbper ),
  m_restaureState_done( true )  
{
  string ptype;
  switch( nbper )
  {
    case 1: ptype = "PC";
      break;
    case 2: ptype = "BIPC";
      break;      
    case 3: ptype = "TRIPC";
      break;
  }                  
  setType( ptype );
  
  // Si la particule est uni-périodique, on affecte le pointeur sur le vecteur
  // de periodicite de l'obstacle périodique
  // Si la particule est bi- ou tri-periodique, on crée un nouveau vecteur,
  // supposé être la somme des vecteurs de périodicité des obstacles concernés
  if ( m_nperiodes > 1 ) m_periodicTranslation = new Vecteur(trans);
  else m_periodicTranslation = &trans;
  *m_geoFormeVdw += *m_periodicTranslation;
  
  // Initialise les efforts
  InitializeForce( false ) ;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur avec arguments
// !!! utilise en PARALLELE !!!
// Rem: la position passee en arguments dans m[16] est deja celle du clone
// periodique
ParticulePeriodique::ParticulePeriodique( const int &IDref, 
	Particule const* particuleClasseReference,
  	ObstaclePeriodique const* contactReferenceObstacle,
	Vecteur const& trans,
	const int &nbper,	
	const double &vx, const double &vy, const double &vz,
	const double &qrotationx, const double &qrotationy, 
	const double &qrotationz, const double &qrotations,	 
	const double &rx, const double &ry, const double &rz,	 
	const Scalar m[16],
	const ParticuleActivity &activ, 
	const int &tag_ ) :
  Particule( -2, particuleClasseReference, vx, vy, vz,
	qrotationx, qrotationy, qrotationz, qrotations,	 
	rx, ry, rz, m, COMPUTE, 0 ),
  m_reference( NULL ),
  m_referenceID( IDref ),
  m_contactReferenceObstacle( contactReferenceObstacle ),
  m_nperiodes( nbper ),
  m_restaureState_done( true ) 
{
  m_activity = COMPUTE;
  string ptype;
  switch( nbper )
  {
    case 1: ptype = "PC";
      break;
    case 2: ptype = "BIPC";
      break;      
    case 3: ptype = "TRIPC";
      break;
  }                  
  setType( ptype );
  
  // Si la particule est uni-périodique, on affecte le pointeur sur le vecteur
  // de periodicite de l'obstacle périodique
  // Si la particule est bi- ou tri-periodique, on crée un nouveau vecteur,
  // supposé être la somme des vecteurs de périodicité des obstacles concernés
  if ( m_nperiodes > 1 ) m_periodicTranslation = new Vecteur(trans);
  else m_periodicTranslation = &trans;

  // Contrairement au 1er constructeur utilise en sequentiel, l'operation
  // *m_geoFormeVdw += *m_periodicTranslation n'est pas realise dans celui ci 
  // utilise en parallele car la position est deja celle du clone periodique
  // soit celle de la particule de ref translate du vecteur de periodicite
  
  // Initialise les efforts
  InitializeForce( false ) ;
}	




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Destructeur
ParticulePeriodique::~ParticulePeriodique()
{
  if ( m_nperiodes > 1 ) delete m_periodicTranslation;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Mise à jour de la vitesse & position de la particule periodique
void ParticulePeriodique::updateVitessePositionPeriodique()
{
  m_geoFormeVdw->setTransform( *(m_reference->getForme()->getTransform()) );
  *m_geoFormeVdw += *m_periodicTranslation;
  setVitesseTranslation( *(m_reference->getVitesseTranslation()) ); 
  setQuaternionRotation( *(m_reference->getCinematique()->getRotation()) );
  setVitesseRotation( *(m_reference->getVitesseRotation()) );    
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie l'obstacle en contact avec la reference 
ObstaclePeriodique const* ParticulePeriodique::getObstacle() const
{
  return m_contactReferenceObstacle;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie le nombre de periodes de la particule periodique
int ParticulePeriodique::getNbPeriodes() const
{
  return m_nperiodes;
}
  



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie le vecteur de periodicite */
Vecteur const* ParticulePeriodique::getVecteurPeriodique() const
{
  return m_periodicTranslation;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie la particule de reference
Particule* ParticulePeriodique::getPeriodicReference() const
{
  return m_reference;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  Renvoie le numero de la particule de reference
int ParticulePeriodique::getPeriodicReferenceID() const
{
  return m_referenceID;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ajout des forces & moments de contact à la particule maitre du clone 
void ParticulePeriodique::AddForcesFromPeriodicCloneToParticule(
	Scalar temps, Scalar dt, bool ContactforceOutput  ) const
{
  if ( ContactforceOutput )
  {
    m_reference->addBodyForce( *m_somme.getForce() );
    m_reference->addContactForcePP( *m_somme.getForce() );
  }  
  else  
    m_reference->addBodyForce( *m_somme.getForce() );
  m_reference->addMoment( *m_somme.getMoment() ); 
  m_reference->addToCoordinationNumber( m_coordination_number ); 
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Restauration de l'etat
void ParticulePeriodique::restaureState( ConfigurationMemento const* Pmemento_,
  	CineParticuleMemento const* Cmemento_,
	ObstaclePeriodique const* contactReferenceObstacle_,
	Particule* reference_ )
{
  m_geoFormeVdw->setTransform( Pmemento_->m_position );
  m_cinematique->restaureState( Cmemento_ );
  m_contactReferenceObstacle = contactReferenceObstacle_ ;
  m_periodicTranslation = m_contactReferenceObstacle->getPeriode() ;
  m_reference = reference_,
  m_referenceID = m_reference->getID();  
}
