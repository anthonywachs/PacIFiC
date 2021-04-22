#include "LeapFrog_2D.H"
#include "Particule.H"
#include "Torseur.H"
#include <math.h>


// ----------------------------------------------------------------------------
// Constructeur par defaut
LeapFrog_2D::LeapFrog_2D() : 
  CineParticule()
{}




// ----------------------------------------------------------------------------
// Constructeur par copie
LeapFrog_2D::LeapFrog_2D( const LeapFrog_2D &copie ) : 
  CineParticule( copie )
{}




// ----------------------------------------------------------------------------
// Destructeur
LeapFrog_2D::~LeapFrog_2D()
{}




// ----------------------------------------------------------------------------
// Classname for BuilderFactory
string LeapFrog_2D::className() const
{
  return ( "*LeapFrog_2D" );
}




// ----------------------------------------------------------------------------
// Construction d'un clone de la cinematique
CineParticule* LeapFrog_2D::clone() const
{
  return ( new LeapFrog_2D( *this ) );
}




// ----------------------------------------------------------------------------
// Evaluation du deplacement  & de la rotation de la particule
void LeapFrog_2D::CalculerVariationQDM( const Torseur &torseur,
	Scalar dt,const Particule* particule )
{
  Scalar facteurCouplage = 1.;
  if ( !(Particule::getExplicitMassCorrection()) )
    facteurCouplage -= 
    	Particule::getFluideMasseVolumique() / particule->getMasseVolumique();
  Scalar masseReduite = particule->getMasse() * facteurCouplage;

  // Quantite de mouvement translationnelle
  m_dUdt = *torseur.getForce() / masseReduite;
  m_dUdt.round();
  m_dUdt[Z] = 0.;

  // Quantite de mouvement rotationnelle
  const double *inertieInverse = particule->getInertieInverse();
  m_dOmegadt[Z] = (*torseur.getMoment())[Z] * inertieInverse[5] 
  	/ facteurCouplage;
}




// ----------------------------------------------------------------------------
// Operateur d'ecriture
// G.FERRER - Juil.2000 - Creation
ostream &operator << ( ostream &fileOut,
	const LeapFrog_2D &cinematique )
{
  fileOut << "*LeapFrog_2D\n";
  fileOut << cinematique.m_vitesseT
	<< cinematique.m_QuaternionRotation
	<< cinematique.m_QuaternionVitesseR;
	
  return fileOut;
}




// ----------------------------------------------------------------------------
// Operateur de lecture
// G.FERRER - Juil.2000 - Creation
istream &operator >> ( istream &fileIn,
	LeapFrog_2D &cinematique )
{
  fileIn >> cinematique.m_vitesseT
	>> cinematique.m_QuaternionRotation
	>> cinematique.m_QuaternionVitesseR;
  cinematique.m_vitesseR = 
  	2.0 * cinematique.m_QuaternionVitesseR.multConjugateToVecteur( 
  	cinematique.m_QuaternionRotation );
  cinematique.m_vitesseR[X] = cinematique.m_vitesseR[Y] = 0.;	
		
  return fileIn;
}




// ----------------------------------------------------------------------------
// Calcul du terme de masse ajouté explicite sur les moments
Vecteur LeapFrog_2D::calculIdwExplicite( const Vecteur &dw,
  	const double *inertie ) const
{
  Vecteur Idw;

  // Calcul de I.dw
  Idw[Z] = inertie[5] * dw[Z];  
  
  return Idw;
}
