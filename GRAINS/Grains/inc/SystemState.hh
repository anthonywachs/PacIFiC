#ifndef _SystemState
#define _SystemState

#include <list>
#include <string>
#include <iostream>
#include "Memento.hh"
using std::list;
using std::string;
using std::cout;
using std::endl;

class EnsComposant;
class LinkedCell;
class App;


/** @brief Etat d'un clone periodique */
struct EtatCloneMonoPeriodique
{
  ConfigurationMemento* memento_config; /**< configuration de la particule */
  Particule* pPerClone; /**< pointeur sur le clone periodique */
  ObstaclePeriodique const* contactReferenceObstacle; /**< Obstacle en 
  	contact avec la particule de reference */
};


/** @brief Etat d'une particule */
struct EtatParticule
{
  int PartID; /**< ID de la particule */
  ConfigurationMemento* memento_config; /**< configuration de la particule */
  CineParticuleMemento* memento_cine; /**< cinematique de la particule */
  int PartClasse; /**< classe de la particule */
  list<struct EtatCloneMonoPeriodique*>* EtatPeriodicClones; /**< liste des 
  	etats des clones mono-periodiques */
};



/** @brief Classe pour la sauvegarde de l'etat du systeme a un temps donne

@author A.WACHS - Institut Francais du Petrole - 2010 - Creation */
//=============================================================================
class SystemState
{
public:
  /**@name Constructors */
  //@{
  /** @brief Constructeur par defaut
  @param temps_sauvegarde temps physique de sauvegarde
  @param ENScomposants ensemble des composants */
  SystemState( Scalar const& temps_sauvegarde,
  	EnsComposant const& ENScomposants );

  /** @brief Destructeur */
  ~SystemState();
  //@}


  /**@name Methods */
  //@{
  /** @brief Retablit l'etat du systeme tel qu'il a ete sauvegarde
  @param temps_sauvegarde temps physique de sauvegarde
  @param ENScomposants ensemble des composants 
  @param LC grille de cellules 
  @param rank rang du processus */    
  void RestaureState( Scalar& temps_sauvegarde, 
  	EnsComposant& ENScomposants,
	LinkedCell* LC,
	const int &rank );  
	
  /** @brief Detruit un etat
  @param etat etat de la particule */
  static void emptyEtatParticule( struct EtatParticule* etat ) ;		 
  //@}


private:
  /** @name Parameters */
  //@{
  Scalar m_temps_sauvegarde; /**< temps physique de sauvegarde */
  list<struct EtatParticule*> allEtats; /**< etat de l'ensemble des particules
  	sur ce processeur */
  list<struct ObstacleState*> obsStates; /**< etat de l'ensemble des 
  	obstacles */
  //@}


  /**@name Methods */
  //@{

  //@}


  /**@name Constructors */
  //@{
  /** @brief Constructeur par defaut */
  SystemState() {}
  
  /** @brief Constructeur de recopie */
  SystemState(const SystemState& Sys) {}  
  //@}

};

#endif
  
