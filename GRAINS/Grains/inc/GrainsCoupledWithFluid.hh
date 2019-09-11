#ifndef _GrainsCoupledWithFluid
#define _GrainsCoupledWithFluid

#include "Grains.H"
#include "InterfaceFluide.hh"
#include "SystemState.hh"
#include "ReaderXML.hh"
#include "AppFluide_LubricationCorrection.H"
#include <list>
#include <string>
using namespace std;



/** @brief Interface entre la modelisation fluide et la modelisation sec.

    @author G.FERRER - Institut Francais du Petrole - 2004 - Creation */
// ============================================================================
class GrainsCoupledWithFluid : virtual public Grains
{
public:
  /** @name Constuctors & Destructor */
  //@{
  /** @brief Constructeur
  @param rhoFluide Masse volumique du fluide dans lequel baignent les
  particules
  @param grid_size size of the smallest grid cell */
  GrainsCoupledWithFluid( Scalar rhoFluide, double grid_size );

  /** @brief Destructeur */
  virtual ~GrainsCoupledWithFluid();
  //@}


  /** @name Methods */
  //@{
  /** @brief Construction de la simulation.
  @param rootElement Le noeud racine */
  virtual void Chargement( DOMElement* rootElement );

  /** @brief Defini le pas temps (couplage avec Basilisk).
  @param dtfluid Le pas de temps fluide */
  virtual void set_timeStep( double const dtfluid );
  virtual void set_ExplicitAddedMasstimeStep( double dtfluid ); 
  
  /** @brief Construction du probleme.
  @param rootElement Le noeud racine */
  virtual void Construction( DOMElement* rootElement );

  /** @brief Construction des forces actives
  @param rootElement Le noeud racine */
  virtual void Forces( DOMElement* rootElement );

  /** @brief Sauvegarde de l'etat de simulation */
  virtual void Save( const string &ext ) const;

  /** @brief Appel a la simulation granulaire.
  @param predict Simulation en etat de prediction
  @param isPredictorCorrector le sch�ma de couplage avec le fluide est il de
  type pr�dicteur-correcteur
  @param explicit_added_mass mass reduite explicite */
  virtual void Simulation( bool predict = true,
  	bool isPredictorCorrector = false,
	bool explicit_added_mass = false );

  /** @brief Mise � jour de la vitesse des particules.
  @param b_set_velocity_nm1_and_diff mise ou non a jour des vitesses et des
  differences de vitesse du pas de temps flduide precedent */
  virtual void UpdateParticulesVelocities(
  	const bool &b_set_velocity_nm1_and_diff );

  /** @brief Ecriture des particules dans un fichier
  @param filename nom du fichier */
  virtual void WriteParticulesInFluid( const string &filename ) const;

  /** @brief Mise � jour de la vitesse des particules.
  @param b_set_velocity_nm1_and_diff mise ou non a jour des vitesses et des
  differences de vitesse du pas de temps flduide precedent
  @param velocities tableau contenant les vitesse de translation & rotation
  @param b_set_velocity_nm1_and_diff mise a jour de la vitesse au pas de temps
  precedent ainsi que la difference de vitesse explicite */
  virtual void UpdateParticulesVelocities(
  	const vector<vector<double> > &velocities,
  	const bool &b_set_velocity_nm1_and_diff );

  /** @brief Recuperation de l'etat du pas de temps precedent et mise a jour
  des vitesses venant du fluide
  @param velocities tableau contenant les vitesse de translation & rotation */
  virtual void InitializeCorrectorStep(
  	const vector<vector<double> > &velocities );

  /** @brief Ecriture des particules dans un flux
  @param is flux d'entr�e */
  virtual void WriteParticulesInFluid( istringstream &is ) const;


  /** @brief Ecriture des particules dans une structure Basilisk
  @param is structure Basilisk d'entr�e */
  virtual void WriteParticulesInFluid (BasiliskDataStructure * is) const;

  /** @brief Modification de la vitesse des particules de Grains par Basilisk
      @param is structure de donnee basilisk contenant les vitesses de translation et rotation
*/
  virtual void UpdateParticulesVelocities (
	BasiliskDataStructure * is,
	const bool &b_set_velocity_nm1_and_diff);


  /** @brief Ecriture de la vitesse et du centre de gravite des particules dans
  un fichier
  @param filename nom du fichier */
  virtual void WritePVGCInFluid( const string &filename ) const;

  /** @brief Ecriture de la vitesse et du centre de gravite des particules dans
  un flux
  @param is flux d'entr�e */
  virtual void WritePVGCInFluid( istringstream &is ) const;

  /** @brief Ajout de l'application masse ajout�e explicite
  @param restart reprise d'un calcul couple PeliGRIFF
  @param PelDirRes repertoire de resultats PeliGRIFF */
  void AddExplicitAddedMass( bool const& restart, string const& PelDirRes );

  /** @brief Initialisation de l'application masse ajout�e explicite en cas de
  reprise de calcul PeliGRIFF
  @param restart reprise d'un calcul couple avec un solveur fluide
  @param dirRes repertoire de resultats
  @param fluidsolver nom du solveur fluide (PeliGRIFF ou Basilisk) */
  void InitializeExplicitAddedMassRestart( bool const& restart,
  	string const& dirRes, string const& fluidsolver = "PeliGRIFF" );

  /** @brief Sauvegarde par defaut de l'etat initial pour post-processing */
  virtual void InitialPostProcessing();

  /** @brief Sauvegarde pour post-processing et restart */
  virtual void doPostProcessing();

  /** @brief Sauvegarde pour post-processing evolution et restart */
  virtual void doPostProcessingEvo( const double &time );

  /** @brief Operations a effectuer avant appel du destructeur */
  virtual void BeforeDestructor();

  /** @brief Temps de simulation (pas de temps fluide) */
  double getSimulTime() const;

  /**
    @brief Is the lift force active in DEMCFD ?
    @return The boolean value Grains_Exec::m_withLiftForce
  */
  bool getIsLiftActive() const;

  /**
    @brief Is the DEMCFD temperature module active ?
    @return The boolean value Grains_Exec::m_withFluidTemperature
  */
  bool getIsTemperature() const;

  /** @brief Renvoi le vecteur gravite */
  vector<double> getVecteurGravite() const;

  /** @brief Send number of particles on all proc*/
  size_t getNumberOfParticleOnAllProc() const;

  /** @brief Renvoi la dimension de l'espace physique */
  int getDimension() const;

  /** @brief Renvoi la liste des particules */
  list<Particule*>* getListOfParticles();

  /** @brief Definit le temps initial
  @param time0 temps initial */
  void setInitialTime( const Scalar &time0 ) { m_temps = time0; } ;

  /** @brief Definit le numero initial de cycle d'ecriture des fichiers de post
  processing
  @param cycle0 numero de cycle initial */
  void setInitialCycleNumber( const int& cycle0 );

  /** @brief Verifie que le post processing Paraview est actif, sinon le cree
  @param name_ nom des fichiers
  @param root_ racine du nom des fichiers
  @param isBinary ecriture en mode binaire */
  void checkParaviewPostProcessing( const string &name_, const string &root_,
  	const bool &isBinary );

  /** @brief Verifie que le post processing Paraview est actif pour Basilisk, sinon le cree
  @param name_ nom des fichiers
  @param root_ racine du nom des fichiers
  @param isBinary ecriture en mode binaire */
  void checkParaviewPostProcessing( const char * name_, const char * root_,
  	const bool &isBinary );


  /** @brief Verifie que le post processing Matlab est actif, sinon le cree
  @param name_ nom des fichiers
  @param root_ racine du nom des fichiers
  @param isBinary ecriture en mode binaire */
  void checkMatlabPostProcessing( const string &name_, const string &root_,
  	const bool &isBinary );


  /** @brief Force le code a activer le mode reload "same" */
  void setReloadSame() { m_forceReloadSame = true ; } ;

  /** @brief Definit la valeur de la translation de post processing dans
  Paraview pour le cas projection-translation
  @param tvx coordonnee x du vecteur de translation
  @param tvy coordonnee y du vecteur de translation
  @param tvz coordonnee z du vecteur de translation */
  void setParaviewPostProcessingTranslationVector( double const& tvx,
  	double const& tvy, double const& tvz );

  /** @brief Allocate DEM-CFD fluid informations */
  void allocateDEMCFD_FluidInfos();

  /** @brief Definit la viscosite du fluide environnant */
  void setFluidViscosity( double mu );

  /** @brief Set boolean "b_fixed_particles" to true in grains3D.
    Then, we still go through GRAINS3D, we still compute forces, but
    we dont move  particles */
  void set_fixed_particles( );

  /** @brief Set boolean "b_fixed_particles_temp" to true in grains3D.
    Then, we still go through GRAINS3D, we still compute forces, but
    we dont move  particles */
  void set_fixed_particles_temp( );

  /** @brief Get processor's origin coordinates */
  void getProcessorOrigin( double& x_solid, double& y_solid,
  	double &z_solid ) ;

  /** @brief Initialise les clones periodiques */
  virtual void initializeClonesPeriodiques();

  /** @brief Calcule I.w et I.dw pour chaque particule pour evaluer le moment
  hydrodynamique exerce par le fluide
  @param w_dw tableau contenant w et dw pour chaque particule
  @param Iw_Idw tableau contenant I.w et I.dw pour chaque particule */
  void compute_Iw_Idw( vector< vector<double> > const* w_dw,
  	vector< vector<double> >* Iw_Idw );

  /** @brief Ecrit les efforts sur les obstacles dans des fichiers
  @param temps temps physique */
  virtual void outputObstaclesLoad( double temps );

  /** @brief Initialise les fichiers de sortie des efforts sur les obstacles
  @param temps temps physique */
  virtual void initialiseOutputObstaclesLoadFiles( double temps );

  /**
    @brief Forces initialization step
  */
  void InitializeForcesT0( const double &fluidThermalConductivity_ );

  /** @brief Return pointer to the lubrication class  */
  static AppFluide_LubricationCorrection* LubricationCorrection() ;
  //@}


protected:
  /** @name Methods */
  //@{
  /** @brief Sauvegarde de l'etat de simulation */
  void saveState();

  /** @brief Restauration de l'etat de simulation */
  void restaureState();
  //@}


  /** @name Parameters */
  //@{
  double m_gridsize; /** Smallest grid size used for lubrication model */

  //~~~ Variables ajoutees pour couplage avec Basilisk
  double m_N; /** Nombre d'iterations granulaire */
  double m_dtgmin; /** Pas de temps granulaire minimum */
  double m_dtgmax; /** Pas de temps granulaire maximum */
  //~~~ Variables ajoutees pour couplage avec Basilisk

  Scalar m_rayon;  /** Maximum radius */
  InterfaceFluide *m_InterfaceFluide; /** Couplage avec le probleme fluide */
  App* m_explicitAddedMass; /** Masse ajout�e explicite */
  Scalar m_simulTime; /** Temps de simulation du probleme de contacts */
  bool m_bd_completed; /** si les operations a effectuer avant appel du
  	destructeur ont effectivement ete realisees */
  SystemState* m_etatSysteme; /** sauvegarde de l'etat du systeme pour schema
  	predictor-corrector */
  bool m_forceReloadSame; /** Force le code a activer le reload de type "same",
  	soit la suite d'une precedente simu */
//  App* m_hydroForce; /** Force hydro dans le couplage DEM-CFD */
//  App* m_demcfdTemperature; /** Temperature dans le couplage DEM-CFD */
  bool b_fixed_particles; /** we still compute forces, but
    we dont move  particles */
  bool b_fixed_particles_temp; /** Taking into account particle temperature
    moving */
  static AppFluide_LubricationCorrection* m_lubricationForce;
      /** Force de correction de lubrication */
  //@}


  /**@name Methods */
  //@{
  /** @brief Insertion d'une particule dans les algorithmes
  @param mode mode d'insertion */
  virtual bool insertParticule( const PullMode& mode );

  /** @brief Positionne les particules en attente pour la simulation a partir
  d'un fichier de positions
  @param mode mode d'insertion */
  virtual void setPositionParticulesFichier( const PullMode& mode = ORDER );

  /** @brief Positionne les particules en attente pour la simulation sous forme
  d'un bloc structure de positions
  @param mode mode d'insertion */
  virtual void setPositionParticulesBloc( const PullMode& mode = ORDER );

  /** @brief Add added mass info structure to the particles when demcfd model
  used */
  virtual void AddAddedMassInfos2Particles_demcfd();
  //@}
};

#endif
