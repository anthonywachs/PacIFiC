#ifndef _AppFluide_Puits
#define _AppFluide_Puits

#include "App.H"

#include "Point.H"

#include <list>
#include <iostream>
using namespace std;


/** @brief Application avec une force fluide de type puits c.a.d. 
    dirigee vers un point et qui decroit lineairement avec la distance à ce 
    point pour finalement s'annuler au delà d'une distance critique. 

    Package APP.

    @author G.FERRER - Institut Francais du Petrole - 2003 - Creation */
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class AppFluide_Puits : public App
{
public:
  /**@name Constructors */
  //@{
  /** @brief Constructeur */
  AppFluide_Puits();

  /** @brief Constructeur 
  @param fileIn flux d'entree */
  AppFluide_Puits( istream &fileIn );

  /** @brief Constructeur
  @param root Le noeud "<FluidePuits>" */
  AppFluide_Puits( DOMNode* root );

  /** @brief Destructeur */
  ~AppFluide_Puits();
  //@}


  /**@name Methods */
  //@{
  /** @brief Calcul des forces d'interaction 
  @param time temps physique 
  @param dt pas de temps granulaire 
  @param particules liste de particules */
  virtual void CalculerForces( Scalar time, Scalar dt,
  	list<Particule*> const* particules );
  //@}


protected:
  /**@name Parameters */
  //@{  
  Point puits; /**< Puits: point de convergence du fluide */
  Scalar f0; /**< Force fluide au point de convergence */
  Scalar distance; /**< Distance critique du gradient de force fluide */
  Scalar deltaP; /**< Gradient de force fluide (pression) */
  //@}


private:
  /**@name Methods */
  //@{
  /** @brief Lecture des parametres d'application 
  @param fileIn flux d'entrée */
  void read( istream &fileIn );
  //@}
};

#endif
