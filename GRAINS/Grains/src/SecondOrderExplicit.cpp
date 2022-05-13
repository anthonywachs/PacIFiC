#include "SecondOrderExplicit.hh"


// ----------------------------------------------------------------------------
// Constructeur par defaut
SecondOrderExplicit::SecondOrderExplicit() :
  TimeIntegrator()
{}




// ----------------------------------------------------------------------------
// Constructeur par copie.
SecondOrderExplicit::SecondOrderExplicit(
	const SecondOrderExplicit &copie) :
  TimeIntegrator(copie)
{}




// ----------------------------------------------------------------------------
// Destructeur
SecondOrderExplicit::~SecondOrderExplicit()
{}




// ----------------------------------------------------------------------------
// Contruction d'un clone de l'integrateur en temps
TimeIntegrator* SecondOrderExplicit::clone() const
{
  return new SecondOrderExplicit(*this);
}




// ----------------------------------------------------------------------------
// Integration du mouvement: calcul des nouvelles vitesses & position
void SecondOrderExplicit::Deplacer(Vecteur& vitesseT,
  	Vecteur const &dUdt,
	Vecteur &deplacementTranslationnel,
	Vecteur const &dOmegadt,
	Vecteur &vitesseR,
	Vecteur &vitesseRMoyenne, 
	double dt)
{
  // Vitesse et deplacement translationnels
  deplacementTranslationnel = vitesseT * dt + 0.5 * dUdt * dt * dt;
  vitesseT += dUdt * dt;

  // Vitesse et deplacement rotationnels
  vitesseRMoyenne = vitesseR + 0.5 * dOmegadt * dt;  
  vitesseR += dOmegadt * dt;  
}
