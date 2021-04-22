#include "SecondOrderRungeKutta.hh"

// ----------------------------------------------------------------------------
// Constructeur par defaut
SecondOrderRungeKutta::SecondOrderRungeKutta() :
  TimeIntegrator()
{}




// ----------------------------------------------------------------------------
// Constructeur par copie.
SecondOrderRungeKutta::SecondOrderRungeKutta(
	const SecondOrderRungeKutta &copie) :
  TimeIntegrator(copie)
{
  intermediate_step = copie.intermediate_step;
  prev_translational_displacement_contribution = copie.prev_translational_displacement_contribution;
  prev_translational_velocity_contribution = copie.prev_translational_velocity_contribution;
  prev_acc = copie.prev_acc;
  prev_vel = copie.prev_vel;
}




// ----------------------------------------------------------------------------
// Destructeur
SecondOrderRungeKutta::~SecondOrderRungeKutta()
{}




// ----------------------------------------------------------------------------
// Contruction d'un clone de l'integrateur en temps
TimeIntegrator* SecondOrderRungeKutta::clone() const
{
  return new SecondOrderRungeKutta(*this);
}




// ----------------------------------------------------------------------------
// Integration du mouvement: calcul des nouvelles vitesses & position
void SecondOrderRungeKutta::Deplacer(Vecteur& vitesseT,
  	Vecteur const &dUdt,
	Vecteur &deplacementTranslationnel,
	Vecteur const &dOmegadt,
	Vecteur &vitesseR,
	Vecteur &vitesseRMoyenne,
	double dt)
{
  // Vitesse et deplacement translationnels
  // // Here we code the Raslton method
  // if (intermediate_step==false)
  // {
  //     deplacementTranslationnel= 2.*dt/3 * vitesseT;
  //     vitesseT += 2.*dt/3 * dUdt;
  //     prev_translational_displacement_contribution = 2.*dt/3 * vitesseT;
  //     prev_translational_velocity_contribution = 2.*dt/3 * dUdt;
  //     prev_vel=vitesseT;
  //     prev_acc=dUdt;
  // }
  // else
  // {
  //     deplacementTranslationnel = -prev_translational_displacement_contribution
  //                                 + dt*(0.25*prev_vel+0.75*vitesseT);
  //    vitesseT += -prev_translational_velocity_contribution + dt*(0.25*prev_acc + 0.75*dUdt);
  // }
  if (intermediate_step==false)
  {
      deplacementTranslationnel= dt/2. * vitesseT;
      prev_translational_displacement_contribution = dt/2. * vitesseT;
      prev_translational_velocity_contribution = dt/2. * dUdt;
      vitesseT += dt/2. * dUdt;
      prev_vel=vitesseT;
      prev_acc=dUdt;
  }
  else
  {
      deplacementTranslationnel = -prev_translational_displacement_contribution
                                  + dt*vitesseT;
     vitesseT += -prev_translational_velocity_contribution + dt*dUdt;
  }

  // Vitesse et deplacement rotationnels
  // No rotation for now
  vitesseR += Vecteur(0.,0.,0.);
  vitesseRMoyenne = Vecteur(0.,0.,0.);

  intermediate_step = !intermediate_step;
}
