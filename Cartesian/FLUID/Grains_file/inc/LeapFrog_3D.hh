#ifndef _LeapFrog_3D
#define _LeapFrog_3D

#include "CineParticule.H"
#include "Vecteur.H"
#include "Quaternion.H"
using namespace solid;


class LeapFrog_3D;
ostream &operator << ( ostream &fileOut, const LeapFrog_3D &cinematique );
istream &operator >> ( istream &fileIn, LeapFrog_3D &cinematique );


/** @brief Cinematique de tanslation et rotation par schema de Leap_Frog 
    sur les particules.

    @author F.PRADEL - Institut Francais du Petrole - 2000 - Creation */
// ============================================================================
class LeapFrog_3D : public CineParticule
{
public:
  /**@name Constructors */
  //@{
  /** @brief Constructeur par defaut */
  LeapFrog_3D();

  /** @brief Constructeur par copie.
  @param copie La cinematique a copier. */
  LeapFrog_3D( const LeapFrog_3D &copie );

  /** @brief Destructeur. */
  ~LeapFrog_3D();
  //@}


  /** @name Methods */
  //@{
  /** @brief Classname for BuilderFactory */
  string className() const;

  /** @brief Contruction d'un clone de la cinematique.
  @return Le clone construit. */
  virtual CineParticule* clone() const;

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
	const LeapFrog_3D &cinematique );
	
  /** @brief Operateur de lecture.
  @return Flux emetteur
  @param fileIn Flux emetteur
  @param cinematique Objet courant */
  friend istream &operator >> ( istream &fileIn, 
	LeapFrog_3D &cinematique );
  //@}

};

#endif
