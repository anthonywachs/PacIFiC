#ifndef _ObstacleChargement_F
#define _ObstacleChargement_F

#include "Vecteur.H"
#include "ObstacleChargement.H"
using namespace solid;
#include <list>
#include <string>
#include <vector>
#include <iostream>
using namespace std;

class Obstacle;

/** @brief Description du chargement en force sur un obstacle.

    Il constitue la definition sur un espace de temps du chargement en force 
    (~confinement) a appliquer sur un obstacle composite ou unitaire.

    @author G.FERRER - Institut Francais du Petrole - 2003 - Creation */
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class ObstacleChargement_F
{
public:
  /**@name Constructeur */
  //@{
  /** @brief Constructeur par defaut */
  ObstacleChargement_F();

  ObstacleChargement_F( DOMNode* root, double dt, int rank );

  /** @brief Destructeur */
  ~ObstacleChargement_F();
  //@}


  /**@name Methods */
  //@{
  /** @brief Nom de l'obstacle associe au chargement */
  string getNom() const;

  /** @brief Temps de chargement actif dans l'espace indique.
  @return Le temps de chargement (<= fin - debut). */
  Scalar getTime( Scalar debut, Scalar fin ) const;

  /** @brief Le chargement d'obstacle est il actif ?
  @return Vrai si termine
  @param t temps physique
  @param dt Increment de temps */
  bool isActif( Scalar t, Scalar dt ) const;
  
  /** @brief Le chargement d'obstacle est il termine ?
  @return Vrai si actif
  @param t temps physique
  @param dt Increment de temps */
  bool isCompleted( Scalar t, Scalar dt ) const;  

  /** @brief Vitesse de deplacement en force 
  @param temps Temps de debut
  @param dt Increment de temps
  @param obstacle Obstacle en deplacement*/
  Vecteur const* VitesseTranslation( Scalar temps, Scalar dt, 
      Obstacle* obstacle );

  /** @brief Renvoie la force en mode cyclique */
  Vecteur cyclicForce( Scalar temps ) const ;

  /** @brief Force imposee */
  Vecteur getForce() const ;

  /** @brief Masse affectee a l'obstacle. */
  Scalar getMasse() const ;

  /** @brief Point de direction du deplacement. */
  Vecteur const* getDirection() const ;

  /** @brief Type de chargement */
  string getType() const; 
  //@}


  /**@name Methods Static */
  //@{
  /** @brief Construction du chargement par lecture des donnees 
  @param fileIn flux d'entrée */
  static ObstacleChargement_F* read( istream &fileIn );
  //@}


private:
  /**@name Parameters */
  //@{  
  static Vecteur m_prev; /**< Force sur l'osbtacle a t et t-dt */
  string m_nomObstacle; /**< Nom de l'obstacle associe au chargement */  
  double m_tdebut; /**< Temps de debut */
  double m_tfin; /**< Temps de fin */
  Vecteur m_force; /**< Force de reaction a atteindre */
  string m_type; /**< Type de chargement */
  Scalar m_masse; /**< Masse de l'obstacle */
  Vecteur m_direction; /**< Direction du deplacement (~de la force) */
  Vecteur m_vitesse_translation; /**< Vitesse de translation */
  Scalar m_frequenceX; /**< Frequence du cisaillement dans la direction X */ 
  Scalar m_frequenceY; /**< Frequence du cisaillement dans la direction Y */ 
  Scalar m_frequenceZ; /**< Frequence du cisaillement dans la direction Z */ 
  Scalar m_phase; /**< Dephasage*/ 
  Vecteur m_origin; /**< Origin du repere de l'obstacle par rapport au repere
      global */
  //@}
};

#endif
