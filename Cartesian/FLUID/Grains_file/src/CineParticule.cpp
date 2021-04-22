#include "CineParticule.H"
#include "Particule.H"
#include "Torseur.H"
#include "Memento.hh"
#include "TimeIntegrator_BuilderFactory.hh"
#include "TimeIntegrator.hh"


// ----------------------------------------------------------------------------
// Constructeur par defaut
CineParticule::CineParticule() :
  m_timeIntegrationScheme( NULL ),
  m_memento( NULL )
{
  m_timeIntegrationScheme = TimeIntegrator_BuilderFactory::create();
}




// ----------------------------------------------------------------------------
// Constructeur par copie
CineParticule::CineParticule( const CineParticule &copie ) :
  Cinematique( copie ),
  m_memento( NULL )
{
  m_vitesseT = copie.m_vitesseT;
  m_vitesseR = copie.m_vitesseR;
  m_QuaternionVitesseR = copie.m_QuaternionVitesseR;
  m_QuaternionRotation = copie.m_QuaternionRotation;
  m_timeIntegrationScheme = copie.m_timeIntegrationScheme->clone();
}




// ----------------------------------------------------------------------------
// Destructeur
CineParticule::~CineParticule()
{
  if ( m_timeIntegrationScheme ) delete m_timeIntegrationScheme;
  if ( m_memento ) delete m_memento;
}




// ----------------------------------------------------------------------------
// Deplacement de la particule
Scalar CineParticule::Deplacer( Particule* particule, double dt )
{
  // Integration en temps
  m_timeIntegrationScheme->Deplacer( m_vitesseT,
  	m_dUdt, deplacementTranslationnelOverDt,
	m_dOmegadt, m_vitesseR, vitesseRMoyenne, dt );

  particule->Translate( deplacementTranslationnelOverDt );

   // Deplacement rotationnel
   double nOmega = Norm( vitesseRMoyenne );
   if ( nOmega > EPS )
   {
     double c = cos( nOmega * dt / 2. );
     double s = sin( nOmega * dt / 2. );
     Vecteur t;
     t = ( s * 1. / nOmega ) * vitesseRMoyenne;
     QuaternionRotationOverDt.setQuaternion( t, c );
   }
   else
     QuaternionRotationOverDt.setQuaternion( 0., 0., 0., 1. );

   m_QuaternionRotation = QuaternionRotationOverDt * m_QuaternionRotation;
   m_QuaternionVitesseR = 0.5 * ( m_vitesseR , m_QuaternionRotation );
   particule->Rotate( QuaternionRotationOverDt );

  return Norm( deplacementTranslationnelOverDt );
}




// ----------------------------------------------------------------------------
// Renvoi de la rotation
Quaternion const* CineParticule::getRotation() const
{
  return &m_QuaternionRotation;
}




// ----------------------------------------------------------------------------
// Vitesse de rotation
Vecteur const* CineParticule::getVitesseRotation() const
{
  return &m_vitesseR;
}




// ----------------------------------------------------------------------------
// Vitesse de translation
Vecteur const* CineParticule::getVitesseTranslation() const
{
  return &m_vitesseT;
}




// ----------------------------------------------------------------------------
// Mise a zero de la cinematique
void CineParticule::reset()
{
  m_vitesseT = 0.;
  m_vitesseR = 0.;
  m_QuaternionVitesseR = 0.;
  m_QuaternionRotation.setQuaternion( 0. , 0., 0., 1. );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Modification de la vitesse de translation
void CineParticule::setVitesseTranslation( const Vecteur &vitesse )
{
  m_vitesseT = vitesse;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Modification de la vitesse de rotation
void CineParticule::setVitesseRotation( const Vecteur &omega )
{
  m_vitesseR = omega;
  m_QuaternionVitesseR = 0.5 * ( omega , m_QuaternionRotation );
}




// ----------------------------------------------------------------------------
// Vitesse en un point de la particule
// A.WACHS - Octo.2010 - Modif
Vecteur CineParticule::Vitesse( const Vecteur &om ) const
{
  return ( m_vitesseT + ( m_vitesseR ^ om ) );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Sauvegarde de l'etat
void CineParticule::saveState()
{
  if ( !m_memento ) m_memento = new CineParticuleMemento();

  m_memento->m_QuaternionRotation = m_QuaternionRotation;
  m_memento->m_vitesseT = m_vitesseT;
  m_memento->m_QuaternionVitesseR = m_QuaternionVitesseR;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Cree et renvoie l'etat
CineParticuleMemento* CineParticule::createState()
{
  CineParticuleMemento* memento_ = new CineParticuleMemento();

  memento_->m_QuaternionRotation = m_QuaternionRotation;
  memento_->m_vitesseT = m_vitesseT;
  memento_->m_QuaternionVitesseR = m_QuaternionVitesseR;

  return memento_;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Restauration de l'etat
void CineParticule::restaureState()
{
  m_QuaternionRotation = m_memento->m_QuaternionRotation;
  m_vitesseT = m_memento->m_vitesseT;
  m_QuaternionVitesseR = m_memento->m_QuaternionVitesseR;
  m_vitesseR = 2.0 * m_QuaternionVitesseR.multConjugateToVecteur(
  	m_QuaternionRotation );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Restauration de l'etat
void CineParticule::restaureState( CineParticuleMemento const* memento_ )
{
  m_QuaternionRotation = memento_->m_QuaternionRotation;
  m_vitesseT = memento_->m_vitesseT;
  m_QuaternionVitesseR = memento_->m_QuaternionVitesseR;
  m_vitesseR = 2.0 * m_QuaternionVitesseR.multConjugateToVecteur(
  	m_QuaternionRotation );
}




// ----------------------------------------------------------------------------
// Operateur d'ecriture
ostream &operator << ( ostream &fileOut, const CineParticule &cinematique )
{
  fileOut << cinematique.className() << '\n';
  fileOut << cinematique.m_vitesseT
	  << cinematique.m_QuaternionRotation
	  << cinematique.m_QuaternionVitesseR;

  return fileOut;
}




// ----------------------------------------------------------------------
// Operateur d'ecriture avec un format de precision elevee
// A.WACHS - Janv.2011 - Creation.
void CineParticule::writeCineParticule( ostream &fileOut ) const
{
  fileOut << className() << endl;
  m_vitesseT.writeGroup3( fileOut );
  fileOut << endl;
  m_QuaternionRotation.writeQuaternion( fileOut );
  fileOut << endl;
  m_QuaternionVitesseR.writeQuaternion( fileOut );
}




// ----------------------------------------------------------------------
// Operateur d'ecriture avec un format de precision elevee
// A.WACHS - Aout 2014 - Creation.
void CineParticule::writeCineParticule2014( ostream &fileOut ) const
{
  m_vitesseT.writeGroup3( fileOut );
  fileOut << " ";
  m_QuaternionRotation.writeQuaternion( fileOut );
  fileOut << " ";
  m_QuaternionVitesseR.writeQuaternion( fileOut );
}




// ----------------------------------------------------------------------
// Ecriture de l'objet en binaire sur le flux de sortie
// A.WACHS - Aout 2014 - Creation.
void CineParticule::writeCineParticule2014_binary( ostream &fileOut )
{
  m_vitesseT.writeGroup3_binary( fileOut );
  m_QuaternionRotation.writeQuaternion_binary( fileOut );
  m_QuaternionVitesseR.writeQuaternion_binary( fileOut );
}




// ----------------------------------------------------------------------------
// Operateur de lecture
istream &operator >> ( istream &fileIn, CineParticule &cinematique )
{
  fileIn >> cinematique.m_vitesseT
	 >> cinematique.m_QuaternionRotation
	 >> cinematique.m_QuaternionVitesseR;
  cinematique.m_vitesseR =
  	2.0 * cinematique.m_QuaternionVitesseR.multConjugateToVecteur(
  	cinematique.m_QuaternionRotation );

  return fileIn;
}




// ----------------------------------------------------------------------------
// Lecture de l'objet en binaire sur le flux d'entree
// avec le format de reload 2014
void CineParticule::readCineParticule2014_binary( istream &StreamIN )
{
  m_vitesseT.readGroup3_binary( StreamIN );
  m_QuaternionRotation.readQuaternion_binary( StreamIN );
  m_QuaternionVitesseR.readQuaternion_binary( StreamIN );

  m_vitesseR = 2.0 * m_QuaternionVitesseR.multConjugateToVecteur(
  	m_QuaternionRotation );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Modification du quaternion de rotation.
void CineParticule::setQuaternionRotation( const Scalar &vecteur0,
	const Scalar &vecteur1,
	const Scalar &vecteur2,
	const Scalar &scalaire )
{
  m_QuaternionRotation.setQuaternion( vecteur0, vecteur1, vecteur2, scalaire );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Modification du quaternion de rotation.
void CineParticule::setQuaternionRotation( const Quaternion &qrot )
{
  m_QuaternionRotation = qrot;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Copie de la cinematique au temps t-2dt: vitesse de translation,
// vitese de rotation, variation de QDM translationnelle, variation de QDM
// rotationnalle
void CineParticule::copyCinematiqueNm2( double *vit, int i ) const
{
  m_timeIntegrationScheme->copyCinematiqueNm2( vit, i );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Cinematique au temps t-2dt: vitesse de translation,
// vitese de rotation, variation de QDM translationnelle, variation de QDM
// rotationnalle
void CineParticule::setCinematiqueNm2( double const* tab )
{
  m_timeIntegrationScheme->setCinematiqueNm2( tab );
}
