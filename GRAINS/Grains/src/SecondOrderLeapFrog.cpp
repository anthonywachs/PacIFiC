#include "SecondOrderLeapFrog.hh"


// ----------------------------------------------------------------------------
// Constructeur par defaut
SecondOrderLeapFrog::SecondOrderLeapFrog() :
  TimeIntegrator()
{}




// ----------------------------------------------------------------------------
// Constructeur par copie.
SecondOrderLeapFrog::SecondOrderLeapFrog(
	const SecondOrderLeapFrog &copie) :
  TimeIntegrator(copie)
{}




// ----------------------------------------------------------------------------
// Destructeur
SecondOrderLeapFrog::~SecondOrderLeapFrog()
{}




// ----------------------------------------------------------------------------
// Contruction d'un clone de l'integrateur en temps
TimeIntegrator* SecondOrderLeapFrog::clone() const
{
  return new SecondOrderLeapFrog(*this);
}




// ----------------------------------------------------------------------------
// Integration du mouvement: calcul des nouvelles vitesses & position
void SecondOrderLeapFrog::Deplacer(Vecteur& vitesseT,
  	Vecteur const &dUdt,
	Vecteur &deplacementTranslationnel,
	Vecteur const &dOmegadt,
	Vecteur &vitesseR,
	Vecteur &vitesseRMoyenne, 
	double dt)
{
  // Vitesse et deplacement translationnels
  vitesseT += dUdt * dt;
  deplacementTranslationnel = vitesseT * dt;

  // Vitesse et deplacement rotationnels
  vitesseR += dOmegadt * dt;
  vitesseRMoyenne = vitesseR;    
}
