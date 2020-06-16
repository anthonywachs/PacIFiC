#ifndef _TimeIntegrator__
#define _TimeIntegrator__

#include "Vecteur.H"
#include "Quaternion.H"


/** @brief Schéma d'intégration en temps.

    @author A.WACHS - Institut Francais du Petrole - 2011 - Creation */
// ============================================================================
class TimeIntegrator
{
public:
  /**@name Constructors */
  //@{
  /** @brief Destructeur. */
  virtual ~TimeIntegrator();
  //@}


  /** @name Methods */
  //@{
  /** @brief Contruction d'un clone de l'integrateur en temps
  @return Le clone construit. */
  virtual TimeIntegrator* clone() const = 0;

  /** @brief Integration du mouvement: calcul des nouvelles vitesses & position
  */
  virtual void Deplacer(Vecteur& vitesseT,
  	Vecteur const &dUdt,
	Vecteur &deplacementTranslationnel,
	Vecteur const &dOmegadt,
	Vecteur &vitesseR,
	Vecteur &vitesseRMoyenne, 
	double dt) = 0;

  /** @brief Copie de la cinematique au temps t-2dt: vitesse de translation, 
  vitese de rotation, variation de QDM translationnelle, variation de QDM
  rotationnalle
  @param vit vecteur de copie 
  @param i position dans le vecteur vit */
  virtual void copyCinematiqueNm2(double *vit,int i) const {}
  //@}


  /** @name Methods Set */
  //@{  
  /** @brief Cinematique au temps t-2dt: vitesse de translation, 
  vitese de rotation, variation de QDM translationnelle, variation de QDM
  rotationnalle
  @param tab tableau de 12 doubles contenant les valeurs des 4 vecteurs vitesse 
  de translation, vitese de rotation, variation de QDM translationnelle, 
  variation de QDM rotationnalle */
  virtual void setCinematiqueNm2(double const* tab) {}   
  //@}


protected:
  /**@name Constructors */
  //@{
  /** @brief Constructeur par defaut */
  TimeIntegrator();

  /** @brief Constructeur par copie.
  @param copie La cinematique a copier. */
  TimeIntegrator(const TimeIntegrator &copie);
  //@}
};

#endif
