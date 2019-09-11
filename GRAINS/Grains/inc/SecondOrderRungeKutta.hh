#ifndef _SecondOrderRungeKutta__
#define _SecondOrderRungeKutta__

#include "TimeIntegrator.hh"

// WARNING: NEED TO INVESTIGATE IF THERE IS A MEMORY LEAK

class SecondOrderRungeKutta : public TimeIntegrator
{
public:
  ~SecondOrderRungeKutta();

  SecondOrderRungeKutta();

  virtual TimeIntegrator* clone() const;

  virtual void Deplacer(Vecteur& vitesseT,
  	Vecteur const &dUdt,
	  Vecteur &deplacementTranslationnel,
    Vecteur const &dOmegadt,
    Vecteur &vitesseR,
    Vecteur &vitesseRMoyenne,
    double dt);

protected:
  SecondOrderRungeKutta(const SecondOrderRungeKutta &copie);

private:
  bool intermediate_step=false; // True if the current time is between two time steps
  Vecteur prev_translational_displacement_contribution;
  Vecteur prev_translational_velocity_contribution;
  Vecteur prev_acc;
  Vecteur prev_vel;
};

#endif
