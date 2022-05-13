#ifndef _AppFluide_Temperature
#define _AppFluide_Temperature

#include "App.H"
#include "ReaderXML.hh"

#include "EnsComposant.H"

/** @brief Force de trainee hydrodynamique

    Force hydrodynamique opposee au mouvement et proportionnelle a la vitesse
    de la particule via le coefficient de trainee.

    @author G.FERRER - Institut Francais du Petrole - 2004 - Creation */
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class AppFluide_Temperature : public App
{
public:
  /** @name Constructors */
  //@{
  /** @brief Constructeur */
  AppFluide_Temperature();

  /** @brief Constructeur
  @param root Le noeud "<Drag>" */
  AppFluide_Temperature( DOMNode* nTemperature );

  /** @brief Destructeur */
  ~AppFluide_Temperature();
  //@}


  /**@name Methods */
  //@{
  /** @brief Calcul des forces d'interaction 
  @param time temps physique 
  @param dt pas de temps granulaire 
  @param particules liste de particules */
  virtual void CalculerForces( Scalar time, Scalar dt,
  	list<Particule*> const* particules );

  /**
    @brief Initialization step
    @param time temps physique 
    @param dt pas de temps granulaire 
    @param particules liste de particules
  */
  void InitializeTemperature( Scalar time, Scalar dt,
      list<Particule*> const* particules,
      const Scalar fluidThermalConductivity_ );
  /**brief Set fluid time step
  @param fluid time step */
  //@{
  void set_simultime( double dt );
  //@}  
 
  //@}
  /**@name Static Parameters */
  //@{
  static double m_heatCapacityS; /** Solid-bodies heat capacity */
  //@}

protected:
  /**@name Parameters */
  //@{  
  Scalar m_Prandtl;
  //@}


private:
  /**@name Parameters */
  //@{
    string type_Nusselt; /**< modele de Nusselt: RanzMarshall, Gunn, etc */
    Scalar m_cellDiamRatio; /** Size ratio between cell and particle diameter */
    Scalar m_diffCoefS; /** Solid-bodies diffusion coefficient */
    Scalar m_initialTempS; /** Solid-bodies thermal conductivity */
    
    Scalar m_diffCoefF; /** Fluid diffusion coefficient (in case of
        non-resolved fluid dynamics) */
    Scalar m_heatCapacityF; /** Fluid heat capacity (in case of
        non-resolved fluid dynamics) */
    Scalar m_initialTempF; /** Fluid initial temperature (in case of
        non-resolved fluid dynamics) */
    bool  b_fluctNu; /** is with stochastic drag law? */
    double tau; /** Stochastic time scale */
    double b_fluidTimeStep; /** fluid time step */
  //@}

  /**@name Methods */
  //@{
  /** @brief Generate random number using Box Muller dist
  @param mean
  @param stddev standard deviation */
  double  rand_normal(double mean, double stddev );


  //@}  
};

#endif
