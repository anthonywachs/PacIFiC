#ifndef _AddedMass
#define _AddedMass

#include "App.H"

#include "ReaderXML.hh"


/** @brief Terme de masse ajoutée

  Calcul du terme de masse ajoutée pour un traitement explicite

  @author A.WACHS - Institut Francais du Petrole - 2009 - Creation */
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class AddedMass : public App
{
public:
  /** @name Constructors */
  //@{
  /** @brief Constructeur 
  @param rhoFluide_ masse volumique du fluide
  @param simulTime_ pas de temps de simulation du fluide */
  AddedMass( Scalar rhoFluide_, Scalar simulTime_ );

  /** @brief Destructeur */
  ~AddedMass();
  //@}


  /**@name Methods */
  //@{
  /** @brief Calcul des forces d'interaction 
  @param time temps physique 
  @param dt pas de temps granulaire 
  @param particules liste de particules */
  virtual void CalculerForces( Scalar time, Scalar dt,
  	list<Particule*> const* particules );
	
  void setsimulTime( double simulTime_ );	
  //@}


private:
  /**@name Parameters */
  //@{
  Scalar rhoFluideInExplicitMass; /**< Masse volumique du fluide dans le calcul 
  	de la masse ajoutée explicite: si rhoFluideInExplicitMass = 0, alors la 
	masse ajoutée est traitée de manière implicite */
  Scalar simulTime; /**< Pas de temps de simulation fluide */ 
  //@}
  
};

#endif
