// Gestion de la cinematique pour des obstacles
// G.FERRER - Juin.2000 - Creation
// ============================================================================
#ifndef _CineObstacle
#define _CineObstacle

#include "Cinematique.H"
#include "ObstacleChargement.H"
#include "Quaternion.H"
#include <stdlib.h>

class MonObstacle;
class Torseur;
class CineObstacleMemento;


class CineObstacle;
ostream &operator << ( ostream &fileOut, CineObstacle &cinematique );
istream &operator >> ( istream &fileIn,  CineObstacle &cinematique );


/** @brief Gestion de la cinematique des obstacles.

    Les obstacles ont deux possibilites de mouvement. Il peut etre impose par
    l'intermediaire d'un programme de chargement : vitesse de deplacement
    constante entre temps-debut et temps-fin. Il peut etre le resultat d'un
    effort sur celles-ci du au contact des particules sur celle-ci.

    @author G.FERRER - Institut Francais du Petrole - 2000 - Creation */
// ============================================================================
class CineObstacle : public Cinematique
{
public:
  /**@name Constructors */
  //@{
  /** @brief Constructeur par defaut */
  CineObstacle();

  /** @brief Destructeur */
  virtual ~CineObstacle();
  //@}

  /**@name Methods */
  //@{
  /** @brief Ajout du chargement a la cinematique de l'obstacle
  @param chargement Le chargement lie a la cinematique */
  void append( ObstacleChargement &chargement );

  /** @brief Decomposition de la cinematique sur la Cinematique.
  @param voisine Cinematique a appliquer
  @param levier Bras de levier de la cinematique a appliquer */
  void Decompose( const CineObstacle &voisine, const Vecteur &levier );

  /** @brief Evaluation du deplacement en fonction des chargements 
  dans l'intervalle de temps specifie.
  @return TRUE su un deplacement est actif
  @param temps Temps de debut
  @param dt Increment de temps. */
  bool Deplacement( Scalar temps, Scalar dt );
  
  /** @brief Calcul de la vitesse relative de l'obstacle
  @param om Le bras de levier permettant de definir la partie rotation. */
  Vecteur Vitesse( const Vecteur &om ) const; 
  
  /** @brief L'obstacle a t il un chargement en rotation actif 
  @param temps temps physique
  @param dt pas de temps */
  bool rotationEnCours( Scalar temps,Scalar dt ) const; 
  //@}
  

  /** @name Accessors */
  //@{
  /** @brief Valeur du quaternion de rotation sur le pas de temps */
  Quaternion const* getQuaternionRotationOverDt() const;

  /** @brief Valeur de la translation imposee a l'obstacle 
  pendant un pas de temps */
  Vecteur const* getTranslation() const;

  /** @brief Vitesse de rotation de l'obstacle
  pendant un pas de temps */  
  Vecteur const* getVitesseRotation() const;

  /** @brief Vitesse de translation de l'obstacle */
  Vecteur const* getVitesseTranslation() const;

  /** @brief Acces au chargement */
  list<ObstacleChargement*> getChargements() const;
  //@}
  

  /** @name Methods Set */
  //@{
  /** @brief Mise a zero de la cinematique. */
  void reset();

  /** @brief Affectation d'une valeur des vitesses initiale a la cinematique */
  void set( CineObstacle &cinematique );
  
  /** @brief Affectation de la vitesse
  @param vitesseTranslation vitesse de translation 
  @param vitesseRotation vitesse de rotation */
  void setVitesse( Vecteur const* vitesseTranslation,
  	Vecteur const* vitesseRotation );   
  //@}


  /** @name Processus de restauration */
  //@{
  /** @brief Sauvegarde de l'etat */
  void saveState();
  
  /** @brief Cree et renvoie l'etat */
  CineObstacleMemento* createState();  
  
  /** @brief Restauration de l'etat */
  void restaureState();
  
  /** @brief Restauration de l'etat 
  @param memento_ cinematique */
  void restaureState( CineObstacleMemento const* memento_ );  
  //@}
  

  /**@name Methods Friends */
  //@{
  /** @brief Operateur d'ecriture.
  @return Flux recepteur
  @param fileOut Flux recepteur
  @param cinematique Objet courant. */
  friend ostream &operator << ( ostream &fileOut, 
	CineObstacle &cinematique );
	
  /** @brief Operateur de lecture.
  @return Flux emetteur
  @param fileIn Flux emetteur
  @param cinematique Objet courant */
  friend istream &operator >> ( istream &fileIn, 
	CineObstacle &cinematique );
  //@}


private:
  /**@name Methods */
  //@{
  /** @brief Destruction de l'ensemble des chargements d'obstacle */
  void clearAndDestroy();
  //@}


  /**@name Parameters */
  //@{
  Vecteur m_translationOverTimeStep; /**< Translation imposee a l'obstacle 
  	pendant un pas de temps */
  Vecteur m_rotationOverTimeStep; /**< Rotation imposee a l'obstacle 
  	pendant un pas de temps  */
  Quaternion m_QuaternionRotationOverDt; /**< Quaternion associe a la rotation 
  	pendant un pas de temps (rotation instantanée) */
  list<ObstacleChargement*> m_chargements; /**< Ensemble des chargements 
  	lies a l'obstacle : translation ou rotation */  
  Vecteur m_vitesseTranslation; /**< vitesse de translation de l'obstacle */
  Vecteur m_vitesseRotation; /**< vitesse de rotation de l'obstacle */
  CineObstacleMemento *m_memento; /**< pour sauver l'état */    
  //@}
};

#endif
  
