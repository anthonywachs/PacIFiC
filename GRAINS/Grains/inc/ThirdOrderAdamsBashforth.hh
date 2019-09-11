#ifndef _ThirdOrderAdamsBashforth__
#define _ThirdOrderAdamsBashforth__

#include "TimeIntegrator.hh"



class ThirdOrderAdamsBashforth : public TimeIntegrator
{
public:
  ~ThirdOrderAdamsBashforth();

  ThirdOrderAdamsBashforth();

  virtual TimeIntegrator* clone() const;

  //To be added: previous acceleration
  virtual void Deplacer(Vecteur& vitesseT,
  	Vecteur const &dUdt,
	  Vecteur &deplacementTranslationnel,
    Vecteur const &dOmegadt,
    Vecteur &vitesseR,
    Vecteur &vitesseRMoyenne,
    double dt);

protected:
  ThirdOrderAdamsBashforth(const ThirdOrderAdamsBashforth &copie);

private:
  Vecteur m_prev_dUdt;
  Vecteur m_2prev_dUdt;
  Vecteur m_prev_v;
  Vecteur m_2prev_v;
  Vecteur prev_translational_displacement_contribution;
  Vecteur prev_translational_velocity_contribution;
  int remaining_use_rk2=4;
};

#endif
