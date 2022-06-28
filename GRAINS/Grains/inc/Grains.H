#ifndef _Grains
#define _Grains

#include "Grains_Exec.hh"
#include "EnsComposant.H"
#include "App.H"
#include "solvercomputingtime.hh"
#include "computingtime.hh"
#include "LinkedCell.H"
#include "AppFluide_Temperature.H"
#include "AppFluide_Drag.H"
#include <list>
#include <vector>
#include <string>
#include <iostream>
using namespace std;

#include "ReaderXML.hh"



/** @brief Classe principale de gestion du probleme GRAINS.

    Construction du probleme en fonction du fichier de description.
    Gestion de la simulation.

    @author G.FERRER - Institut Francais du Petrole - 1999 - Creation */
//=============================================================================
class Grains : public ComputingTime, public SolverComputingTime
{
public:
  /** @name Constructors */
  //@{
  /** Constructeur par defaut */
  Grains();

  /** Destructeur */
  virtual ~Grains();
  //@}


  /** @name Methods Static */
  //@{
  /** @brief Schema de resolution: predicteur ou correcteur
  @param predictor vrai si le schema est predicteur, sinon correcteur */
  static void setMode( const bool &predictor );

  /** @brief Schema de resolution: predicteur ou correcteur */
  static bool isModePredictor();
  //@}


  /** @name Methods Virtual */
  //@{
  /** @brief Construction de la simulation
  @param rootElement Le noeud racine */
  virtual void Chargement( DOMElement* rootElement );

  /**
    @brief Construction du probleme.
    @param rootElement Le noeud racine
    @remark Here we read "Construction" module information from the input file,
      namely :
        - Domain dimension
        - Eventual domain origine (otherwise [0., 0., 0.])
        - Domain decomposition
        - Eventual TimeIntegration scheme (otherwise SecondOrderLeapFrog)
        - Eventual reload type (otherwise start from scratch)
        - Particles geometrical properties
        - Contacts properties (calls Contact_BuilderFactory::create)
        - Obstacles informations (calls Obstacle_BuilderFactory::create)
  */
  virtual void Construction( DOMElement* rootElement );

  /** @brief Construction des forces actives
  @param rootElement Le noeud racine */
  virtual void Forces( DOMElement* rootElement );

  /** @brief Appel a la simulation granulaire
  @param predict Simulation en etat de prediction
  @param isPredictorCorrector le sch�ma de couplage avec le fluide est il de
  type pr�dicteur-correcteur
  @param explicit_added_mass la masse ajout�e est elle trait�e de mani�re
  explicite */
  virtual void Simulation( bool predict = true,
  	bool isPredictorCorrector = false,
	bool explicit_added_mass = false );

  /** @brief Operations particulieres pour la sauvegarde de l'etat initial
  pour post-processing */
  virtual void InitialPostProcessing() {};

  /** @brief Operations particulieres pour la sauvegarde pour post-processing
  et restart */
  virtual void doPostProcessing() {};
  //@}

  /** @name Methods */
  //@{
  /** @brief Allocate fluid informations */
  void allocateDEMCFD_FluidInfos();
  //@}


  /** @name Parameters */
  //@{
  static Vecteur m_vgravite; /**< vecteur gravite */
  //@}


protected:
  /** @name Mode d'insertion (position) des particules */
  //@{
  enum InsertMode {
    BEFORE, /** Insertion avant simulation */
    INLINE, /** Insertion au cours de la simulation */
    NOINSERT /** Pas d'insertion */
  };  //@}

  /** @name Mode d'initialisation de la vitesse des particules */
  //@{
  enum InitVitMode {
    ZERO, /** initialisation a zero */
    CONSTANT, /** Vitesse constante (translation et rotation) */
    RANDOM_MOTION /** Vitesse aleatoire (translation et rotation) */
  };  //@}


  /** @name Parameters */
  //@{
  EnsComposant m_composants; /** Liste des composants :
  	particules & obstacles */
  list<App*> m_allApp; /** Liste des applications utilisees dans la simulation;
      	la premiere application correspond au probleme sec */
  LinkedCell* m_sec; /** App Sec : Particule-Particule */
  AppFluide_Drag* app_HydroForce; /** Drag and potential lift in
     DEMCFD coupling */
  AppFluide_Temperature* app_FluidTemperature; /** Temperature dans le couplage
      DEM-CFD */
  AppFluide_Temperature* app_SolidTemperature; /** Heat transfert between
      particles */
  Scalar m_tdeb; /** Temps de debut de simulation */
  Scalar m_tfin; /** Temps de fin de simulation */
  Scalar m_dt; /** Pas de temps de simulation */
  Scalar m_temps; /** Temps physique */
  list<Scalar> m_save; /** Temps de Sauvegarde */
  string m_fileSave; /** Nom du fichier de sauvegarde en fin de simulation */
  int m_dimension; /** dimension du domaine = 2 ou 3 */
  bool m_allProcTiming; /** chaque processus mesure les temps pass�s dans
      les diff�rentes parties du code */
  bool m_new_reload_format; /** format des fichiers de reload (nouveau ou
      ancien) */
//  bool m_dragForce_with_fluidAtRest; /** activation de la force hydro
//      lorsque le fluide n'est pas resolu (fluide suppose au repos) */
  //@}

  /** @name Parameters pour la gestion des insertions */
  //@{
  PullMode m_mode_insertion_particules; /** Mode d'insertion des particules */
  InsertMode m_methode_insertion_particules; /**< Methode d'insertion
      des particules */
  InitVitMode m_methode_initvit_particules; /** Mode d'initialisation de la
      vitesse des particules */
  vector<Fenetre> m_fenetres; /** Fenetres d'insertion des particules */
  vector<Fenetre> m_constrainteFenetres; /** Fenetres de calcul des contraintes */
  string m_position; /** Fichier des positions */
  struct BlocInsertion* m_blocInsert; /** Bloc d'insertion structure */
  bool m_configAleatoire; /** Si vrai, choisit aleatoirement la configuration
  	des particules */
  list< pair<Particule*,int> > m_newParticules; /** classes de nouvelles
  	particules a inserer */
  int m_insertion_frequency; /** Insertion tous les m_insertion_frequency
  	pas de temps */
  bool m_force_insertion; /** Force l'insertion meme en cas de contact, option
  	utilisable avec la methode Preinstall uniquement */
  vector<Vecteur const*> m_periodicVectors; /** liste des vecteurs periodiques
  	des obstacles periodiques pour test d'insertion */
  double m_RandomMotionCoefTrans; /** coefficient pour l'amplitude du mouvement
  	aleatoire initial en translation */
  double m_RandomMotionCoefRot; /** coefficient pour l'amplitude du mouvement
  	aleatoire initial en rotation */
  Vecteur m_InitVtrans; /** Initialisation de la vitesse de translation de la
	particule a l'insertion */
  Vecteur m_InitVrot; /** Initialisation de la vitesse de rotation de la
	particule a l'insertion */
  Scalar m_VolumeIn; /** Current volume in simulation (porosity) */
  Scalar m_VolumeInLimit; /** Limit of void ratio due to compression */
  //@}

  /** @name Parameters pour le couplage avec le fluide */
  //@{
  static bool m_predictor_mode; /**< Mode de calcul: predicteur ou correcteur */
  //@}

  /** @name Parameters MPI generique */
  //@{
  int m_rank; /** Rang du processus dans le communicateur MPI_COMM_activProc
  	 (=0 en sequentiel) */
  int m_nprocs; /** nombre total de processus dans le communicateur
  	MPI_COMM_activProc (=1 en sequentiel) */
  bool m_processorIsActiv; /** le processus est il actif ?
  	(=true en sequentiel) */
  //@}

  /**@name Methods */
  //@{
  /** @brief Selection d'un point dans les fenetres d'insertion */
  Point getPoint() const;

  /** @brief Selection aleatoire d'une transformation de type rotation */
  Matrix getRandomRotation() const;

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

  /** @brief Lecture des donnees pour les simulations MPI
  @param lx taille du domaine global en X
  @param ly taille du domaine global en Y
  @param lz taille du domaine global en Z
  @param root le noeud XML */
  virtual void readDomainDecomposition( DOMNode* root,
  	const Scalar& lx, const Scalar& ly, const Scalar& lz );

  /** @brief Renvoie le nom complet du fichier de resultat .result
  @param rootname nom de la racine */
  virtual string fullResultFileName( const string &rootname ) const;

  /** @brief Lecture des donnees pour les simulations periodiques
  @param rootElement le noeud XML */
  virtual void readPeriodic( DOMElement* rootElement );

  /** @brief Definition du LinkedCell
  @param rayon rayon max des particules */
  virtual void defineLinkedCell( Scalar const& rayon );

  /** @brief Arret du code en cas de probleme */
  virtual void grainsAbort() const;

  /** @brief Sauvegarde des fichiers de reload
  @param temps temps physique */
  virtual void saveReload( Scalar const& temps );

  /** @brief Numero maximum de particules */
  virtual int numeroMaxParticules() const;

  /** @brief Cree, insert et link les nouvelles particules */
  virtual void InsertCreateNewParticules();

  /** @brief Efface les fichiers .result et .xml */
  void clearResultXmlFiles() const;

  /** @brief Ecriture de la memoire utilisee par la simulation */
  virtual void display_used_memory() const;

  /** Synchronize the PPWindow boolean relative to each sub-domain.
  vector of bool is then correct on every processor*/
  virtual void synchronize_PPWindow() {};

  /**
    @brief Renvoie la vitesse initiale d'une particule a partir du mode
        d'initialisation
    @param vtrans Particle translational velocity
    @param vrot Particle rotational velocity
  */
  void computeInitVit( Vecteur& vtrans, Vecteur& vrot ) const;
  //@}

};

#endif
