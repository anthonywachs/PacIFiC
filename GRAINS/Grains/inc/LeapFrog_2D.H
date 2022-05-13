#ifndef _LeapFrog_2D
#define _LeapFrog_2D

#include "CineParticule.H"


class LeapFrog_2D;
ostream &operator << ( ostream &fileOut, const LeapFrog_2D &cinematique );
istream &operator >> ( istream &fileIn, LeapFrog_2D &cinematique );


/** @brief Cinematique de tanslation et rotation par schema de Leap_Frog 
    sur les particules.

    @author F.PRADEL - Institut Francais du Petrole - 2000 - Creation */
// ============================================================================
class LeapFrog_2D : public CineParticule
{
public:
  /**@name Constructors */
  //@{
  /** @brief Constructeur par defaut */
  LeapFrog_2D();

  /** @brief Constructeur par copie.
  @param copie La cinematique a copier. */
  LeapFrog_2D( const LeapFrog_2D &copie );

  /** @brief Destructeur. */
  ~LeapFrog_2D();
  //@}


  /**@name Methods */
  //@{
  /** @brief Classname for BuilderFactory */
  string className() const;

  /** @brief Contruction d'un clone de la cinematique.
  @return Le clone construit. */
  CineParticule* clone() const;

  /** @brief Calcul de la variation de quantite de mouvement au temps precedent
  @param torseur Torseur somme des forces appliquees sur la particule 
  @param dt increment de temps. 
  @param particule Particule associee a la cinematique. */
  void CalculerVariationQDM( const Torseur &torseur,
	Scalar dt,const Particule* particule );

  /** @brief Calcul du terme de masse ajouté explicite sur les moments
  @param dw difference explicite de vitesse de rotation 
  @param inertie tenseur d'inertie de la particule */
  virtual Vecteur calculIdwExplicite( const Vecteur &dw,
  	const double *inertie ) const;
  //@}


  /**@name Methods Friends */
  //@{
  /** @brief Operateur d'ecriture.
  @return Flux recepteur
  @param fileOut Flux recepteur
  @param cinematique Objet courant. */
  friend ostream &operator << ( ostream &fileOut, 
	const LeapFrog_2D &cinematique );
	
  /** @brief Operateur de lecture.
  @return Flux emetteur
  @param fileIn Flux emetteur
  @param cinematique Objet courant */
  friend istream &operator >> ( istream &fileIn, 
	LeapFrog_2D &cinematique );
  //@}
};

#endif
