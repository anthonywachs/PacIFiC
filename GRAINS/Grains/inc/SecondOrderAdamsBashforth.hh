#ifndef _SecondOrderAdamsBashforth__
#define _SecondOrderAdamsBashforth__

#include "TimeIntegrator.hh"


/** @brief Sch�ma d'int�gration en temps 2� ordre Adams-Bashforth
   x(t+dt)=x(t)+0.5*(3*v(t)-v(t-dt)) puis v(t+dt)=v(t)+0.5*(3*a(t)-a(t-dt))

    @author A.WACHS - Institut Francais du Petrole - 2011 - Creation */
// ============================================================================
class SecondOrderAdamsBashforth : public TimeIntegrator
{
public:
  /**@name Constructors */
  //@{
  /** @brief Destructeur. */
  ~SecondOrderAdamsBashforth();

  /** @brief Constructeur par defaut */
  SecondOrderAdamsBashforth();
  //@}


  /** @name Methods */
  //@{
  /** @brief Contruction d'un clone de l'integrateur en temps
  @return Le clone construit. */
  virtual TimeIntegrator* clone() const;

  /** @brief Integration du mouvement: calcul des nouvelles vitesses & position
  */
  virtual void Deplacer(Vecteur& vitesseT,
  	Vecteur const &dUdt,
	Vecteur &deplacementTranslationnel,
	Vecteur const &dOmegadt,
	Vecteur &vitesseR,
	Vecteur &vitesseRMoyenne,
	double dt);

  /** @brief Copie de la cinematique au temps t-2dt: vitesse de translation,
  vitese de rotation, variation de QDM translationnelle, variation de QDM
  rotationnalle
  @param vit vecteur de copie
  @param i position dans le vecteur vit */
  virtual void copyCinematiqueNm2(double *vit,int i) const;
  //@}


  /** @name Methods Set */
  //@{
  /** @brief Cinematique au temps t-2dt: vitesse de translation,
  vitese de rotation, variation de QDM translationnelle, variation de QDM
  rotationnalle
  @param tab tableau de 12 doubles contenant les valeurs des 4 vecteurs vitesse
  de translation, vitese de rotation, variation de QDM translationnelle,
  variation de QDM rotationnalle */
  virtual void setCinematiqueNm2(double const* tab);


protected:
  /**@name Constructors */
  //@{
  /** @brief Constructeur par copie.
  @param copie La cinematique a copier. */
  SecondOrderAdamsBashforth(const SecondOrderAdamsBashforth &copie);
  //@}


private:
  /** @name Parameters */
  //@{
  Vecteur m_vitesseT_nm2; /**< Vitesse de translation a t-2*dt */
  Vecteur m_vitesseR_nm2; /**< Vitesse de rotation a t-2*dt */
  Vecteur m_dUdt_nm2; /**< variation de quantite de mouvement translationnelle
	a t-2*dt */
  Vecteur m_dOmegadt_nm2; /**< variation de quantite de mouvement rotationnelle
	a t-2*dt */
  bool first_step=true;
  //@}
};

#endif
