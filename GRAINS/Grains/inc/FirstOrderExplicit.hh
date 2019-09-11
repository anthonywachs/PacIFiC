#ifndef _FirstOrderExplicit__
#define _FirstOrderExplicit__

#include "TimeIntegrator.hh"


/** @brief Schéma d'intégration en temps 1er ordre explicite classique:
   x(t+dt)=x(t)+dt*v(t) puis v(t+dt)=v(t)+dt*a(t). 

    @author A.WACHS - Institut Francais du Petrole - 2011 - Creation */
// ============================================================================
class FirstOrderExplicit : public TimeIntegrator
{
public:
  /**@name Constructors */
  //@{
  /** @brief Destructeur. */
  ~FirstOrderExplicit();
  
  /** @brief Constructeur par defaut */
  FirstOrderExplicit();  
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
  FirstOrderExplicit(const FirstOrderExplicit &copie);
  //@}
};

#endif
