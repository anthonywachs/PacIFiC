#include "Grains.H"
#include "LeapFrog_Sphere.hh"
#include "Particule.H"
#include "Torseur.H"
#include <math.h>


// ----------------------------------------------------------------------------
// Constructeur par defaut
LeapFrog_Sphere::LeapFrog_Sphere() : 
  LeapFrog_3D()
{}




// ----------------------------------------------------------------------------
// Constructeur par copie
LeapFrog_Sphere::LeapFrog_Sphere( const LeapFrog_Sphere &copie ) : 
	LeapFrog_3D( copie )
{}




// ----------------------------------------------------------------------------
// Destructeur
LeapFrog_Sphere::~LeapFrog_Sphere()
{}




// ----------------------------------------------------------------------------
// Construction d'un clone de la cinematique
CineParticule* LeapFrog_Sphere::clone() const
{
  return ( new LeapFrog_Sphere( *this ) );
}




// ----------------------------------------------------------------------------
// Evaluation du deplacement  & de la rotation de la particule
void LeapFrog_Sphere::CalculerVariationQDM( const Torseur &torseur,
	Scalar dt, const Particule* particule )
{
  Scalar facteurCouplage = 1.;
  if ( Particule::getMassCorrection() && 
  				!Particule::getExplicitMassCorrection() )
    facteurCouplage -= 
    	Particule::getFluideMasseVolumique() / particule->getMasseVolumique();

  // Quantite de mouvement translationnelle
  Scalar masseReduite = particule->getMasse() * facteurCouplage;
  m_dUdt = *torseur.getForce() / masseReduite;
  m_dUdt.round();

  // Quantite de mouvement rotationnelle
  Vecteur const* Moment = torseur.getMoment();

  const double *inertieInverse = particule->getInertieInverse();
  m_dOmegadt[0] = inertieInverse[0] * (*Moment)[0] / facteurCouplage;
  m_dOmegadt[1] = inertieInverse[3] * (*Moment)[1] / facteurCouplage;
  m_dOmegadt[2] = inertieInverse[5] * (*Moment)[2] / facteurCouplage;
}
