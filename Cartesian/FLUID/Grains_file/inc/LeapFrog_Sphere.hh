#ifndef _LeapFrog_Sphere
#define _LeapFrog_Sphere

#include "CineParticule.H"
#include "LeapFrog_3D.hh"
#include "Vecteur.H"
#include "Quaternion.H"
using namespace solid;


/** @brief Cinematique de translation et rotation par schema de Leap_Frog 
    sur des particules spheriques.

    @author A.WACHS - Institut Francais du Petrole - 2009 - Creation */
// ============================================================================
class LeapFrog_Sphere : public LeapFrog_3D
{
public:
  /**@name Constructors */
  //@{
  /** @brief Constructeur par defaut */
  LeapFrog_Sphere();

  /** @brief Constructeur par copie.
  @param copie La cinematique a copier. */
  LeapFrog_Sphere(const LeapFrog_Sphere &copie);

  /** @brief Destructeur. */
  ~LeapFrog_Sphere();
  //@}


  /** @name Methods */
  //@{
  /** @brief Contruction d'un clone de la cinematique.
  @return Le clone construit. */
  CineParticule* clone() const;

  /** @brief Calcul de la variation de quantite de mouvement au temps precedent
  @param torseur Torseur somme des forces appliquees sur la particule 
  @param dt increment de temps. 
  @param particule Particule associee a la cinematique. */
  void CalculerVariationQDM( const Torseur &torseur,
	Scalar dt,const Particule* particule );
  //@}

};

#endif
