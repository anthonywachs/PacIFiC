#ifndef _AppSec
#define _AppSec

#include "App.H"
#include "EnsComposant.H"
#include <iostream>
#include <list>
using namespace std;

class ContactMemento;
class ObstaclePeriodique;
class MonObstacle;


/** @brief Definition des methodes Sec decrivant la resolution
    du probleme sec pour les particules.

    Package APP.

    @author G.FERRER - Institut Francais du Petrole - 2003 - Creation */
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class AppSec : public App
{
public:
  /** @name Contructors & Destructor*/
  //@{
  /** @brief Destructeur */
  virtual ~AppSec();
  //@}


  /** @name Virtual Methods */
  //@{
  /** @brief Calcul des forces d'interaction
  @param time temps physique
  @param dt pas de temps granulaire
  @param particules liste de particules */
  virtual void CalculerForces( Scalar time, Scalar dt,
  	list<Particule*> const* particules ) = 0 ;

  /** @brief Calcul des forces d'interaction pour l'initialisation des
  conditions de cohesion.
  @param dt pas de temps granulaire
  @param time temps physique
  @param particules liste de particules */
  virtual void InitializeCohesiveForces( Scalar time, Scalar dt,
  	list<Particule*> const* particules ) {}

  /** @brief Indication du diametre max des particules a l'algorithme.
  @param dist diametre max des particules */
  virtual void define( double dist ) = 0;

  /** @brief La particule est elle en contact avec un composant de l'App ?
  @return TRUE si contact
  @param particule particule */
  virtual bool isContact( const Particule *particule ) const {return false;};

  /** @brief La particule est elle en contact avec un composant de l'App ?
  variante VdW
  @return TRUE si contact
  @param particule particule */
  virtual bool isContactVdW( const Particule *particule ) const {return false;};

  /** @brief La particule est elle proche d'un composant de l'App ?
  @return TRUE si contact
  @param particule particule */
  virtual bool isProche( const Particule *particule ) const {return false;};

  /** @brief La particule est elle proche d'un composant de l'App ?
  variante VdW
  @return TRUE si contact
  @param particule particule */
  virtual bool isProcheVdW( const Particule *particule ) const {return false;};

  /** @brief Association de la particule avec l'algorithme .
  sans controle de contact.
  @param particule particule */
  virtual void Link( Particule* particule )=0;

  /** @brief Association de l'obstacle avec l'algorithme.
  @param obstacle obstacle */
  virtual void Link( Obstacle* obstacle );

  /** @brief Actualisation des associations des particules avec l'algorithme.
  Une particule dont l'etat devient inactif est supprimee de la gestion.
  @param time temps physique
  @param dt pas de temps granulaire
  @param particules liste de particules */
  virtual void LinkUpdate( Scalar time, Scalar dt,
  	list<Particule*>* particules )=0;

  /** @brief Suppression de la particule de l'APPlication
  @param particule particule a supprimer */
  virtual void remove( Particule* particule )=0;

  /** @brief Suppression de l'obstacle de l'APPlication
  @param obs obstacle a supprimer */
  virtual void remove( MonObstacle* obs );
  //@}


  /** @name Methods */
  //@{
  /** @brief Caracteristiques des contacts
  @param nbPart nombre de particules */
  void computeMeanNbParticules( size_t const& nbPart );

  /** @brief Caracteristiques des contacts
  @param time temps
  @param contactPoint point de contact */
  void addToContactsFeatures( Scalar time, PointContact & contactPoint );

  /** @brief Renvoie un obstacle en fonction de son numero
  @param num numero de l'obstacle */
  MonObstacle* getMonObstacle( int const& num ) const;

  /** @brief Renvoie un obstacle periodique en fonction de son numero (si la
  simulation n'est pas periodique, renvoie NULL)
  @param num numero de l'obstacle */
  ObstaclePeriodique const* getObstaclePeriodique( int const& num ) const;
  //@}


  /** @name Methods static */
  //@{
  /** @brief Renvoie la penetration maximale */
  static Scalar getOverlapMax();

  /** @brief Renvoie la penetration moyenne */
  static Scalar getOverlapMean();

  /** @brief Renvoie l'instant de l'overlap max */
  static Scalar getTimeOverlapMax();

  /** @brief Renvoie le nombre d'iterations moyen de GJK */
  static Scalar getNbIterGJKMean();

  /** @brief Renvoie le nombre moyen de particules par proc */
  static Scalar getNbParticulesPerProcMean();

  /** @brief Renvoie la proprotion de régimes de Coulomb*/
  static Scalar getNbCoulombRegimes();

  /** @brief Affecte les caracteristiques globales des contacts
  @param overlap_max_ penetration maximale
  @param overlap_mean_ penetration moyenne
  @param time_overlapMax_ instant de l'overlap max
  @param nbIterGJK_ le nombre d'iterations moyen de GJK */
  static void setContactsFeatures( Scalar const& overlap_max_,
	  Scalar const& overlap_mean_,
	  Scalar const& time_overlapMax_,
	  Scalar const& nbIterGJK_ );
  //@}

protected:
  /** @name Constructors */
  //@{
  /** @brief Constructeur par defaut */
  AppSec();
  //@}

  /** @name Parameters */
  //@{
  Obstacle* m_obstacles; /**< Obstacles presents dans la resolution */
  list<MonObstacle*> m_allObstacles; /**< Liste des obstacles primaires */
  static Scalar m_overlap_max; /**< Penetration maximale */
  static Scalar m_overlap_mean; /**< Penetration moyenne */
  static Scalar m_time_overlapMax; /**< Instant de l'overlap max */
  static Scalar m_nbIterGJK_mean; /**< Nombre d'iterations moyen de GJK */
  static Scalar m_nbParticules_mean; /**< Nombre moyen de particules par proc */
  static Scalar m_propCoulombFriction; /**< Proportion de frottements de coulomb */
  //@}
};

#endif
