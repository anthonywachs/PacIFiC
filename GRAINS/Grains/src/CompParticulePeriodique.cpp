#include "CompParticulePeriodique.hh"
#include "Memento.hh"


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur avec arguments
// !!! utilise en SEQUENTIEL !!!
CompParticulePeriodique::CompParticulePeriodique( Particule* reference,
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
  ParticulePeriodique( reference, particuleClasseReference, 
	contactReferenceObstacle, trans, nbper ),	
  CompParticule( -2, particuleClasseReference, 
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
  m_id = -2;
  updateVitessePositionPeriodique();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Destructeur
CompParticulePeriodique::~CompParticulePeriodique()
{
//  if ( m_nperiodes > 1 ) delete m_periodicTranslation;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Mise à jour de la vitesse & position de la particule periodique
void CompParticulePeriodique::updateVitessePositionPeriodique()
{
  m_geoFormeVdw->setTransform( *(m_reference->getForme()->getTransform()) );
  *m_geoFormeVdw += *m_periodicTranslation;
  setVitesseTranslation( *(m_reference->getVitesseTranslation()) ); 
  setQuaternionRotation( *(m_reference->getCinematique()->getRotation()));
  setVitesseRotation( *(m_reference->getVitesseRotation()) ); 
  setElementPosition();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie l'obstacle en contact avec la reference 
ObstaclePeriodique const* CompParticulePeriodique::getObstacle() const
{
  return m_contactReferenceObstacle;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie le nombre de periodes de la particule periodique
int CompParticulePeriodique::getNbPeriodes() const
{
  return m_nperiodes;
}
  



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie le vecteur de periodicite */
Vecteur const* CompParticulePeriodique::getVecteurPeriodique() const
{
  return m_periodicTranslation;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie la particule de reference
Particule* CompParticulePeriodique::getPeriodicReference() const
{
  return m_reference;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  Renvoie le numero de la particule de reference
int CompParticulePeriodique::getPeriodicReferenceID() const
{
  return m_referenceID;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ajout des forces & moments de contact à la particule maitre du clone 
void CompParticulePeriodique::AddForcesFromPeriodicCloneToParticule(
	Scalar temps, Scalar dt, bool ContactforceOutput  ) const
{
  ParticulePeriodique::AddForcesFromPeriodicCloneToParticule( temps, dt, 
  	ContactforceOutput );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Restauration de l'etat
void CompParticulePeriodique::restaureState( 
	ConfigurationMemento const* Pmemento_,
  	CineParticuleMemento const* Cmemento_,
	ObstaclePeriodique const* contactReferenceObstacle_,
	Particule* reference_ )
{
}
