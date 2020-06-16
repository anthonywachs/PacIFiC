#ifndef _AppFluide_Drag
#define _AppFluide_Drag

#include "App.H"
#include "ReaderXML.hh"

#include "EnsComposant.H"

/** @brief Force de trainee hydrodynamique

    Force hydrodynamique opposee au mouvement et proportionnelle a la vitesse
    de la particule via le coefficient de trainee.

    @author G.FERRER - Institut Francais du Petrole - 2004 - Creation */
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class AppFluide_Drag : public App
{
public:
  /** @name Constructors */
  //@{
  /** @brief Constructeur */
  AppFluide_Drag();

  /** @brief Constructeur 
  @param fileIn flux d'entree */
  AppFluide_Drag( istream &fileIn );

  /** @brief Constructeur
  @param root Le noeud "<Drag>" */
  AppFluide_Drag( DOMNode* root );

  /** @brief Destructeur */
  ~AppFluide_Drag();
  //@}


  /**@name Methods */
  //@{
  /** @brief Calcul des forces d'interaction 
  @param time temps physique 
  @param dt pas de temps granulaire 
  @param particules liste de particules */
  virtual void CalculerForces( Scalar time, Scalar dt,
  	list<Particule*> const* particules );

  /**brief Set if slip velocity is on output
  @param slipveloutput */
  void set_slipveloutput( bool slipveloutput );
  //@}

  /**brief Set fluid time step
  @param fluid time step */
  void set_simultime( double dt );
  //@}

  /** @name Parameters */
  //@{
  bool b_slipveloutput; /**< with slip velocity in output ? */	
  //@}

protected:
  /**@name Parameters */
  //@{  
  string type_Fhydro; /**< modele de force de trainee: Ergun, DiFelice, etc */
  bool b_withLiftForces;        /**< are lift forces On or Off ? */
  bool b_withPressureGradient;  /**< do we add the pressure gradient around
    the particle to the total hydro force or not (for model A or B) ? */
  bool b_withAddedMass; /** do we consider added mass effect? */
  double meanD; /**< Sauter Mean Diameter */
  double sum_XY; /** sum_i(Xi*Yi) for polydispersity*/
  //@}


private:
  /**@name Parameters */
  //@{
  bool  b_fluctdrag; /** is with stochastic drag law? */
  double tau; /** Stochastic time scale */
  double b_fluidTimeStep; /** fluid time step */

  //@}

  /**@name Methods */
  //@{
  /** @brief Lecture des parametres d'application 
  @param fileIn flux d'entree */
  void read( istream &fileIn );

  /** @brief Generate random number using Box Muller dist
  @param mean
  @param stddev standard deviation */
  double  rand_normal(double mean, double stddev );

  //@}  
};

#endif
