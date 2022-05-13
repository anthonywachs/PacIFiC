#include "FirstOrderExplicit.hh"


// ----------------------------------------------------------------------------
// Constructeur par defaut
FirstOrderExplicit::FirstOrderExplicit() :
  TimeIntegrator()
{}




// ----------------------------------------------------------------------------
// Constructeur par copie.
FirstOrderExplicit::FirstOrderExplicit(
	const FirstOrderExplicit &copie) :
  TimeIntegrator(copie)
{}




// ----------------------------------------------------------------------------
// Destructeur
FirstOrderExplicit::~FirstOrderExplicit()
{}




// ----------------------------------------------------------------------------
// Contruction d'un clone de l'integrateur en temps
TimeIntegrator* FirstOrderExplicit::clone() const
{
  return new FirstOrderExplicit(*this);
}




// ----------------------------------------------------------------------------
// Integration du mouvement: calcul des nouvelles vitesses & position
void FirstOrderExplicit::Deplacer(Vecteur& vitesseT,
  	Vecteur const &dUdt,
	Vecteur &deplacementTranslationnel,
	Vecteur const &dOmegadt,
	Vecteur &vitesseR,
	Vecteur &vitesseRMoyenne, 
	double dt)
{
  // Vitesse et deplacement translationnels
  deplacementTranslationnel = vitesseT * dt;
  vitesseT += dUdt * dt;

  // Vitesse et deplacement rotationnels
  vitesseRMoyenne = vitesseR;  
  vitesseR += dOmegadt * dt;  
}
