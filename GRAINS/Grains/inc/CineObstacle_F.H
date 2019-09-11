#ifndef _CineObstacle_F
#define _CineObstacle_F

#include "Cinematique.H"
#include "ObstacleChargement_F.H"
#include "Point.H"
#include "Vecteur.H"
using namespace solid;


class Obstacle;


/** @brief Gestion d'une cinematique d'obstacle par un chargement en force.

    @author G.FERRER - Aout.2003 - Institut Francais du Petrole - Creation */
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class CineObstacle_F : public Cinematique
{
public:
  /**@name Constructors */
  //@{
  /** @brief Constructeur par defaut */
  CineObstacle_F();

  /** @brief Destructeur */
  ~CineObstacle_F();
  //@}

  /**@name Methods */
  //@{
  /** @brief Ajout du chargement a la cinematique de l'obstacle
  @param chargement Le chargement lie a la cinematique */
  void append( ObstacleChargement_F &chargement );
  
  /** @brief Decomposition de la cinematique vers la cinematique courante.
  En fonction de la cinematique indiquee, les valeurs de vitesse de
  rotation et translation ont une valeur initiale.
  @param voisine La cinematique a appliquer. 
  @param centre Centre de gravite de l'obstacle sur lequel est applique
      la cinematique courante. */
  void Decompose( const CineObstacle_F &voisine, const Point& centre );

  /** @brief Evaluation des vitesses Translation & Rotation dans l'intervalle 
  de temps specifie.
  @return TRUE su un deplacement est actif
  @param temps Temps de debut
  @param dt pas de temps granulaire
  @param obstacle Obstacle sur lequel s'applique la cinematique courante*/
  bool Deplacement( Scalar temps, Scalar dt, Obstacle* obstacle );

  /** @brief Valeur de la translation dans l'intervalle de temps dt 
  @param dt temps de temps granulaire */
  Vecteur getTranslation( const Scalar dt ) const;

  /** @brief Mise a zero de la cinematique. */
  void reset();

  /** @brief Calcul de la vitesse relative de l'obstacle
  @return La vitesse relative
  @param om Le vecteur permettant de definir la partie rotation. */
  Vecteur VitesseRelative( const Vecteur &om ) const;
  //@}

private:
  /**@name Parameters */
  //@{
  Vecteur m_vitesseT; /**< Composante de vitesse de translation. */
  Vecteur m_vitesseR; /**< Composante de vitesse de rotation. */
  Scalar m_vitesseD; /**< Vitesse de deplacement impose a l'obstacle (norme) */  
  list<ObstacleChargement_F*> m_chargements; /**< Ensemble des chargements 
  	lies a l'obstacle */
  ObstacleChargement_F* m_chargement; /**< Chargement en cours sur l'obstacle */
  Vecteur m_translationOverTimeStep; /**< Translation imposee a l'obstacle 
  	pendant un pas de temps */
  //@}


  /**@name Methods */
  //@{
  /** @brief Destruction de l'ensemble des chargements */
  void clearAndDestroy();
  //@}
};

#endif
