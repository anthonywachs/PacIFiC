#ifndef _Text_PostProcessingWriter__
#define _Text_PostProcessingWriter__

#include "PostProcessingWriter.hh"
#include <iostream>
#include <fstream>
using std::ofstream;

/**
  @brief Classe Text_PostProcessingWriter
  @details Ecriture de la position (P), des vitesses (V) et des nombres
      de contact des particules dans des fichiers.
      Possibilite de sortir les forces de contact, les vitesses relatives F/P,
      la temperature, la force Hydro
  @remarks Inspired from previous PositionVitesse class
  @remarks Here particles are not reconstructed on the master proc anymore,
      only their informations are gathered, using an MPI buffer
  @author M.BERNARD - IFPEN - 2012 - Creation
*/
//=============================================================================
class Text_PostProcessingWriter : public PostProcessingWriter
{
public:
  /** @name Constructors */
  //@{
  /** Constructeur par defaut 
  @param dn noeud XML 
  @param rank_ rang du processus 
  @param nbranks_ nombre total de processus */
  Text_PostProcessingWriter(DOMNode* dn,int const& rank_,int const& nbranks_);

  /** Destructeur */
  virtual ~Text_PostProcessingWriter();
  //@}


  /** @name Methods */
  //@{
  /** @brief Initialisation du post processeur
  @param temps Temps de sauvegarde
  @param dt pas de temps 
  @param particules particules actives
  @param pwait particules en attente
  @param pperiodiques particules periodiques
  @param ParticuleClassesReference classes de reference de particule
  @param obstacle obstacles 
  @param LC grille de cellules
  @param insert_windows fenetres d'insertion */
  void PostProcessing_start(Scalar const& temps,
  	Scalar const& dt,
  	list<Particule*> const* particules,
	list<Particule*> const* pwait,
	list<Particule*> const* pperiodiques,	
	vector<Particule*> const* ParticuleClassesReference,
	Obstacle* obstacle,
	LinkedCell const* LC,
	vector<Fenetre> const& insert_windows );

  /** @brief Ecriture d'evolution 
  @param temps Temps de sauvegarde  
  @param dt pas de temps  
  @param particules particules actives 
  @param pwait particules en attente
  @param pperiodiques particules periodiques  
  @param ParticuleClassesReference classes de reference de particule
  @param obstacle obstacles 
  @param LC grille de cellules */
  void PostProcessing(Scalar const& temps,
  	Scalar const& dt,
  	list<Particule*> const* particules,
	list<Particule*> const* pwait,
	list<Particule*> const* pperiodiques,		
	vector<Particule*> const* ParticuleClassesReference,
	Obstacle *obstacle,
	LinkedCell const* LC );

  /** @brief Clot les ecritures */
  void PostProcessing_end();

  /** @brief Type de post-processing writer */
  string getPostProcessingWriterType() const {return "PV";};

  /** @brief Le post-processing writer est il parallele ? */
  bool isParallelWriter() const { return true; }

  //@}
  /** @name Static variables */
  //@{
  static bool b_totalForce; /** Write down total force experienced by
      each particle */
  static bool b_cumulatedContactForce; /** Write down sum of contact forces
      on particles cumulated between each outputs */
  static bool b_instantaneousContactForce; /** Write down sum of contact forces
      on particles instantaneously at each output */
  static bool b_cumulatedLubriForce; /** Write down sum of lubrication forces
      on particles cumulated between each outputs */
  static bool b_hydroForce; /** Write down hydro force (drag+potential lift)
      experienced by each particle (for DEMCFD) */
  static bool b_slipVelocity; /** Write down fluid/solid slip velocity
      (for DEMCFD) */
  static bool b_temperature; /** Write down particles temperature */
  static bool b_stressTensor; /** Write the stress tensor */
  static bool b_particleStressTensor; /** Write the individual stress tensor */
  //@}
private:
  /** @name Methods */
  //@{
  /** @brief Write-down method for MPI jobs (even on 1 core)
  @param temps Temps de sauvegarde
  @param nb_total_part Total number of particles
  @param cinematique_Global Vector containing the global kinematic */
  void one_output_MPI(Scalar const& temps, size_t &nb_total_part,
  	vector< vector<double> > const* cinematique_Global);
	
  /** @brief Write-down method for sequential jobs (without MPI pattern)
  @param temps Temps de sauvegarde  
  @param particules particules actives
  @param pwait particules en attente */
  void one_output_Standard(Scalar const& temps,
  	list<Particule*> const* particules,
  	list<Particule*> const* pwait);
	
  /** @brief Efface les fichiers resultats */
  void clearResultFiles() const ; 	  
  
  /** @brief Creates output files and open streams
  @param mode File opening mode (here : ios::app) */
  void prepareResultFiles( ios_base::openmode mode ) ;
  
  /**
    @brief Assemble the cinematique_Global object with particles contained
        only on this proc (used in sequential)
    @param particules
    @param pwait 
  */
  vector< vector<double> >* build_myCinematiqueGlobal( 
	const list<Particule*> &particules,
	const list<Particule*> &pwait ) ;
  //@}
  

private:
  ofstream gc_coordinates_x; /** Flux de coordonnees X du centre de gravite */
  ofstream gc_coordinates_y; /** Flux de coordonnees Y du centre de gravite */ 
  ofstream gc_coordinates_z; /** Flux de coordonnees Z du centre de gravite */
  ofstream gc_velocity_x; /** Flux de coordonnees X de vitesse du centre de 
      gravite */ 
  ofstream gc_velocity_y; /** Flux de coordonnees Y de vitesse du centre de 
      gravite */    
  ofstream gc_velocity_z; /** Flux de coordonnees Z de vitesse du centre de 
      gravite */
  ofstream gc_rotation_x; /** Flux de coordonnees X de vitesse de rotation 
      du centre de gravite */ 
  ofstream gc_rotation_y; /** Flux de coordonnees Y de vitesse de rotation 
      du centre de gravite */    
  ofstream gc_rotation_z; /** Flux de coordonnees Z de vitesse de rotation 
      du centre de gravite */
  ofstream coordination_number; /** Flux de nombre de contacts sur chaque
      particule */
  ofstream demcfd_HydroForce_x; /** Flux de coordonnees X de la drag force */
  ofstream demcfd_HydroForce_y; /** Flux de coordonnees Y de la drag force */ 
  ofstream demcfd_HydroForce_z; /** Flux de coordonnees Z de la drag force */
  ofstream demcfd_SlipVel_x; /** Flux de coordonnees X de la slip velocity */
  ofstream demcfd_SlipVel_y; /** Flux de coordonnees Y de la slip velocity */ 
  ofstream demcfd_SlipVel_z; /** Flux de coordonnees Z de la slip velocity */  
  ofstream particle_class; /** Flux de classe des particules */
  ofstream demcfd_ParticleTemperature; /** Particle temperature flux */
  ofstream demcfd_FluidTemperature; /** Fluid temperature interpolated on the particle */
  ofstream demcfd_ParticleHeatflux; /** Particle heatflux */
  ofstream demcfd_ParticleNusselt; /** Particle Nusselt number */
  ofstream contact_force_x; /**< Flux de coordonnees X de force de contact */ 
  ofstream contact_force_y; /**< Flux de coordonnees Y de force de contact */    
  ofstream contact_force_z; /**< Flux de coordonnees Z de force de contact */
  ofstream contact_force_inst_x; /**< Flux de coordonnees X de force de contact */ 
  ofstream contact_force_inst_y; /**< Flux de coordonnees Y de force de contact */    
  ofstream contact_force_inst_z; /**< Flux de coordonnees Z de force de contact */
  ofstream lubri_force_x; /**< Flux de coordonnees X de force de lubrication */ 
  ofstream lubri_force_y; /**< Flux de coordonnees Y de force de lubrication */    
  ofstream lubri_force_z; /**< Flux de coordonnees Z de force de lubrication */
  ofstream gc_force_x; /** Flux de coordonnees X de force */
  ofstream gc_force_y; /** Flux de coordonnees Y de force */
  ofstream gc_force_z; /** Flux de coordonnees Z de force */
  ofstream fIntMmt; /** Flux somme des moments internes des particules */
  ofstream fIntMmt0; /** Flux element de la matrice de contrainte individuelle */
  ofstream fIntMmt1; /** Flux element de la matrice de contrainte individuelle */
  ofstream fIntMmt2; /** Flux element de la matrice de contrainte individuelle */
  ofstream fIntMmt3; /** Flux element de la matrice de contrainte individuelle */
  ofstream fIntMmt4; /** Flux element de la matrice de contrainte individuelle */
  ofstream fIntMmt5; /** Flux element de la matrice de contrainte individuelle */
  ofstream fIntMmt6; /** Flux element de la matrice de contrainte individuelle */
  ofstream fIntMmt7; /** Flux element de la matrice de contrainte individuelle */
  ofstream fIntMmt8; /** Flux element de la matrice de contrainte individuelle */
  string simul; /** racine des fichiers de sortie */
protected:
  /** @name Constructors */
  //@{
  /** @brief Constructeur par defaut */
  Text_PostProcessingWriter() {};  
  //@}  
};

#endif
  
