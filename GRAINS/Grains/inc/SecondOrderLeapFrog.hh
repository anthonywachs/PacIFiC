#ifndef _SecondOrderLeapFrog__
#define _SecondOrderLeapFrog__

#include "TimeIntegrator.hh"


/** @brief Schéma d'intégration en temps 2e ordre Leap frog
    v(t+dt/2)=v(t-dt/2)+dt*a(t) puis x(t+dt)=x(t)+dt*v(t+dt/2).

    @author A.WACHS - Institut Francais du Petrole - 2011 - Creation */
// ============================================================================
class SecondOrderLeapFrog : public TimeIntegrator
{
public:
  /**@name Constructors */
  //@{
  /** @brief Destructeur. */
  ~SecondOrderLeapFrog();
  
  /** @brief Constructeur par defaut */
  SecondOrderLeapFrog();  
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
  SecondOrderLeapFrog(const SecondOrderLeapFrog &copie);
  //@}
};

#endif
