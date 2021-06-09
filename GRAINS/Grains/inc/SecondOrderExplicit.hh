#ifndef _SecondOrderExplicit__
#define _SecondOrderExplicit__

#include "TimeIntegrator.hh"


/** @brief Schéma d'intégration en temps 2è ordre explicite classique:
   x(t+dt)=x(t)+dt*v(t)+0.5*a(t)*dt^2 puis v(t+dt)=v(t)+dt*a(t). 

    @author A.WACHS - Institut Francais du Petrole - 2011 - Creation */
// ============================================================================
class SecondOrderExplicit : public TimeIntegrator
{
public:
  /**@name Constructors */
  //@{
  /** @brief Destructeur. */
  ~SecondOrderExplicit();
  
  /** @brief Constructeur par defaut */
  SecondOrderExplicit();  
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
  //@}


protected:
  /**@name Constructors */
  //@{
  /** @brief Constructeur par copie.
  @param copie La cinematique a copier. */
  SecondOrderExplicit(const SecondOrderExplicit &copie);
  //@}
};

#endif
