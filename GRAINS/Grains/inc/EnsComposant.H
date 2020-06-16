// G.FERRER - Dece.1999 - Creation
// G.FERRER - Aout.2002 - Recherche des voisins
// ============================================================================
#ifndef _EnsComposant
#define _EnsComposant

#include "Basic.H"
#include "Erreur.H"

#include <fstream>
#include <sstream>
#include <list>
#include <set>
using namespace std;
#include "Particule.H"
#include "Obstacle.H"
#include "MonObstacle.H"
#include "WriterXML.hh"


class App;
class AppSec;
class Forme;
class ObstacleChargement;
class ObstacleChargement_F;
class MPIWrapperGrains;
class PostProcessingWriter;
struct Fenetre;


/** @name Enumeration */
//@{
enum PullMode {
  // Tirage des particules dans l'ordre d'attente.
  ORDER,
  // Tirage aleatoire des particules.
  RANDOM,
  // Tirage non effectif.
  NONE
};
//@}


/** @brief Gestion de l'ensemble des composants : particules & obstacles..

    Au cours de la simulation, les particules peuvent changer de physionomie
    pour passer d'une particule elementaire a un amas de particules.
    Ce passage est lie a la forme de la fonction de contact utilisee pour
    les particules. Toutefois, les particules elementaires restent toujours
    accessibles.

    @author G.FERRER - Institut Francais du Petrole - 1999 - Creation */
// ============================================================================
class EnsComposant
{
public:
  /** @name Constructors */
  //@{
  /// @brief Constructeur par defaut.
  EnsComposant();

  /// @brief Destructeur.
  ~EnsComposant();
  //@}


  /** @name Methods */
  //@{
  /** @brief Actualisation des relations entre composants.
  Verification des composants actifs-inactifs. */
  void Actualiser();

  /** @brief Ajout d'une particule a l'ensemble.
  @param particule La particule a ajouter */
  void Ajouter( Particule* particule );

  /** @brief Ajout d'une classe de particules.
  @param particule La particule a ajouter */
  void AjouterClasseParticules( Particule *particule );

  /** @brief Adjust the size of the vector nbParticlesPerClass to
  the number of particle class (used for dry granular simulation, when
  we know initialy the number of class)
  @param NbClass Number of class  */
  void prepare_Polydisperse_vectors( int NbClass );

  /** @brief Set the number of particles per class to the vector
  m_nbParticlesPerClass
  @param i class of particles
  @param nb Number of particles from the class i  */
  void set_NbParticlesPerClass( int i, int nb );

  /** @brief Set the number of particles per class to the vector
  m_nbParticlesPerClass. Used for restarts or coupled simulations */
  void set_NbParticlesPerClass( void );

  /** @brief compute Sauter mean diameter */
  void compute_SauterMeanDiameter( void );

  /** @brief compute Sauter mean diameter */
  void compute_ParticleClassesConcentration( void );

  /** @brief Ajout d'un obstacle a l'ensemble.
  @param obstacle_ L'obstacle a ajouter. */
  void Ajouter( Obstacle* obstacle_ );

  /** @brief Association du chargement a son obstacle.
  @param chargement Le chargement candidat. */
  void Associer( ObstacleChargement &chargement );

  /** @brief Association du chargement a son obstacle.
  @param chargement Le chargement candidat. */
  void Associer( ObstacleChargement_F &chargement );

  /** @brief Somme les forces sur les particules
  @param temps Temps de simulation.
  @param dt Pas de temps */
  void AddForcesFromPeriodicClonesToParticules( Scalar temps, Scalar dt );

  /** @brief Initialise les forces sur les composants et l'indicateur de calcul
  de la transformation avec epaiseur de croute a faux
  @param temps Temps de simulation.
  @param dt Pas de temps
  @param withWeight avec ou sans le poids des particules */
  void InitializeForces( Scalar temps, Scalar dt,
	bool const& withWeight );

  /** @brief Initialize PP forces */
  void InitializePostProcessingForces();

  /** @brief Initialise force at the contact point for the post-processing
  of the stress tensor */
  void InitializeForcesAtContactPoint();

  /** @brief Integrate contact PP forces in time
  @param nb number of analyzed solid time steps between two saving output  */
  void IntegrateContactPostProcessingForces(Scalar nb);

  /** @brief Initialise l'indicateur de calcul de la transformation
  avec epaiseur de croute a faux pour tous les composants
  @param temps Temps de simulation.
  @param dt Pas de temps */
  void InitializeVdWState( Scalar temps, Scalar dt );

  /** @brief Calcul le poids de toutes les particules */
  void computeWeight( Scalar temps, Scalar dt );

  /**
    @brief Deplacement des particules avec prise en compte de
        la force de gravite.
    @exception Si le deplacement admissible d'une particule est depasse.
    @param temps Temps de simulation.
    @param dt Pas de temps.
  */
  list<MonObstacle*> Deplacer( Scalar temps,
		Scalar dt ) throw(ErreurDeplacement);

  /**
    @brief Solve temperature problem
    @param temps Temps de simulation.
    @param dt Pas de temps.
  */
  void ComputeTemperature( Scalar temps, Scalar dt ) ;

  /**
    @brief Post-processing contraintes
    @param fenetre Fenetre de calcul.
  */
  void setStressTensorDomain( vector<Fenetre> fenetres ) ;

  /**
    @brief Compute stress tensor
    @param stress  Stress Tensor.
  */
  void computeStressTensor( vector<Scalar> &stressTensor,
    MPIWrapperGrains const* wrapper ) ;

  /** @brief Association des particules et obstacles a l'algorithme
  @param app algorithme */
  void Link( AppSec &app );

  /** @brief Ecriture des composants pour Post-Processing (Ecriture de
  reference)
  @param temps Temps de sauvegarde
  @param dt pas de temps
  @param LC grille de cellules
  @param insert_windows fenetres d'insertion
  @param rang rang du processeur
  @param nprocs nombre de processeurs
  @param wrapper wrapper MPI */
  void PostProcessing_start( Scalar temps, Scalar dt,
	LinkedCell const* LC, vector<Fenetre> const& insert_windows,
	int rang = 0,
	int nprocs = 1,
  	MPIWrapperGrains const* wrapper = NULL );

  /** @brief Ecriture des composants pour Post-Processing (Ecriture
  d'evolution)
  @param temps Temps de sauvegarde
  @param dt pas de temps
  @param LC grille de cellules
  @param rang rang du processeur
  @param nprocs nombre de processeurs
  @param wrapper wrapper MPI */
  void PostProcessing( Scalar temps, Scalar dt,
	LinkedCell const* LC,
	int rang = 0,
  	int nprocs = 1,
	MPIWrapperGrains const* wrapper = NULL );

  /** @brief Ecriture des composants pour Post-Processing GrainsCompFeatures */
  void PostProcessing();

  /** Fermeture du Post-Processing */
  void PostProcessing_end();

  /** @brief Ecriture des composants pour Post-Processing d'une erreur de
  contact ou de deplacement
  @param filename racine du nom du fichier
  @param errcomposants liste contenant les 2 composants */
  void PostProcessingErreurComposants( string const& filename,
  	list<Composant*> const& errcomposants );

  /** @brief Mise a zero de la cinematique si necessaire.
  @param reset Cinematique = 0. si reset= "Reset". */
  void ResetCinematique( string &reset );

  /** @brief Transfert de la particule inactive en attente vers les actives.
  La particule inactive est definie par la fonction getParticule(PullMode) */
  void ShiftParticuleOutIn();

  /** @brief Cas o� la particule � ins�rer n'est pas dans le sous-domaine:
  on efface le pointeur wait des listes reference et pwait et on d�truit l'objet
  point� par wait */
  void DeleteAndDestroyWait();

  /** @brief Mise � jour des particules clones periodiques
  @param LC grille de cellules */
  void updateClonesPeriodiques( LinkedCell* LC );

  void updateContactMapId( int prev_id, int new_id);

  /** @brief Mise � jour du statut des particules clones paralleles qui sont
  egalement des references periodiques
  @param temps temps physique
  @param PartRefPerHalozone liste (ID,obstacle ID) de particules de
  	reference periodique de la zone de recouvrement
  @param PartRefPerOutDomainHalozone list ID de particules de
  	reference periodique hors du domaine de la zone de recouvrement
  @param InNotRefPerHalozone list (ID,obstacle ID)  de particules ayant perdu
  	leur statut de reference periodique mais toujours dans le domaine
	et la zone de recouvrement
  @param LC grille de cellules */
  void statutClonesReferencesPeriodiques_MPI_Step2( Scalar temps,
  	list<int>& PartRefPerHalozone,
  	list<int>& PartRefPerOutDomainHalozone,
  	list<int>& InNotRefPerHalozone,
	LinkedCell* LC );

  /** @brief Mise � jour du statut des particules clones periodiques en
  parallele
  @param temps temps physique
  @param ClonestoDestroy liste des IDs des clones a detruire
  @param ClonestoParticules liste des IDs des clones a transformer en particules
  	standards
  @param LC grille de cellules */
  void statutClonesPeriodiques_MPI_Step3( Scalar temps,
  	list<int>& ClonestoDestroy,
  	list<int>& ClonestoParticules,
	LinkedCell* LC );

  /** @brief Renvoie les vitesses max et moyenne de l'ensemble des composants
  @param vmax vitesse max
  @param vmean vitesse moyenne
  @param wrapper wrapper MPI */
  void ComputeMaxMeanVelocity( Scalar &vmax, Scalar &vmean,
  	MPIWrapperGrains const* wrapper = NULL ) const;

  /** @brief Compute the void ratio : Volume of void / Volume of particles
  @param vmax vitesse max
  @param vmean vitesse moyenne
  @param wrapper wrapper MPI */
  void ComputeVoidRatio( Scalar &volume, vector<Point> &obsPos,
  	MPIWrapperGrains const* wrapper = NULL ) const;

  /** @brief Ecrit dans un fichier les vitesses min, max et moy des particles en
  fonction du temps
  @param temps Temps de simulation
  @param fileOut fichier de sortie
  @param rang rang du processus
  @param wrapper wrapper MPI */
  void monitorParticlesVelocity( Scalar temps, ofstream& fileOut,
  	int rang = 0, MPIWrapperGrains const* wrapper = NULL ) const;

  /** @brief Mise � jour de la localisation geographique des particules de
  la zone de recouvrement */
  void updateGeoLocalisationParticulesHalozone();

  /** @brief Cree une sauvegarde de l'etat des obstacles
  @param obsStates liste des etats de tous les obstacles */
  void createStateObstacles( list<struct ObstacleState*> &obsStates ) const;

  /** @brief Restaure l'etat des obstacles
  @param obsStates liste des etats de tous les obstacles */
  void restaureStateObstacles( list<struct ObstacleState*> &obsStates );

  /** @brief Affecte la cinematique des obstacles sans les d�placer
  @param temps Temps de simulation.
  @param dt Pas de temps. */
  void setCinematiqueObstacleSansDeplacement( Scalar temps, Scalar dt );
  //@}


  /**@name Methods Get */
  //@{
  /** @brief Recherche d'un obstacle par son nom.
  @return L'obstacle : unitaire ou composite,
  NULL si pas d'obstacle de ce nom.
  @param nom Le nom de l'obstacle. */
  const Obstacle* getObstacle( const string &nom ) const;

  /** @brief Acces aux obstacles */
  Obstacle* getObstacles();

  /** @brief Particule d'indice id. Attention, methode peu efficace.
  @return La particule a l'indice id, NULL si l'indice n'est pas correct.*/
  Particule* getParticule( int id );

  /** @brief Composant d'indice id. Attention, methode peu efficace.
  @return Le composant a l'indice id, NULL si l'indice n'est pas correct.*/
  Composant* getComposant( int id );

  /** @brief Particule cliente en fonction du mode de tirage.
  @return La particule cliente,
  NULL si le mode est NONE ou si pas de particules en attente
  @param mode mode d'insertion
  @param wrapper wrapper MPI*/
  Particule* getParticule( PullMode mode,
  	MPIWrapperGrains const* wrapper = NULL );

  /** @brief Liste des particules actives.
  @return La liste des particules. */
  list<Particule*>* getParticulesActives();

  /** @brief Liste des particules actives.
  @return La liste des particules. */
  list<Particule*> const* getParticulesActives() const;

  /** @brief Liste des particules inactives.
  @return La liste des particules. */
  list<Particule*>* getParticulesWait();

  /** @brief Liste des particules inactives.
  @return La liste des particules. */
  list<Particule*> const* getParticulesWait() const;

  /** @brief Liste des particules dans la zone de recouvrement d'un autre
  processeur.
  @return La liste des particules. */
  list<Particule*>* getParticulesHalozone();

  /** @brief Liste des particules clones.
  @return La liste des particules. */
  list<Particule*>* getParticulesClones();

  /** @brief Liste des particules clones periodiques
  @return La liste des particules. */
  list<Particule*>* getParticulesClonesPeriodiques();

  /** @brief Liste des particules possedant des clones periodiques
  @return La liste des particules. */
  set<Particule*>* getParticulesReferencesPeriodiques();

  /** @brief Vecteur des particules de r�f�rence des classes.
  @return Le vecteur des particules de r�f�rence. */
  vector<Particule*>* getParticuleClassesReference();

  /** @brief Vecteur des particules de r�f�rence des classes.
  @return Le vecteur des particules de r�f�rence. */
  vector<Particule*> const* getParticuleClassesReference() const;

  /** @brief Rayon maximal des particules.
  @return Le plus grand rayon */
  Scalar getRayonMax();

  /** @brief updates the shrinking radius.
  @return the updated radius
  MOSTAFA SULAIMAN ADDED HERE,02-11-2015 */
 Scalar update_shrinking_radius(Scalar dt);

  /** @brief gets the initial radius.
  @return the intial radius
  MOSTAFA SULAIMAN ADDED HERE,02-11-2015*/
 Scalar get_initial_radius();

  /** @brief gets the intial mass.
  @return the intitial mass
  MOSTAFA SULAIMAN ADDED HERE,02-11-2015*/
 Scalar get_initial_mass();

  /** @brief tells if the particle is shrinking
  or not.
  @return Shrinking
  MOSTAFA SULAIMAN ADDED HERE,02-11-2015*/
  bool IsShrinking();


  /** @brief computes the shrinking mass
    using the shrinking radius .
  @return the shrinking mass
  MOSTAFA SULAIMAN ADDED HERE,02-11-2015*/
  Scalar compute_shrinking_weight();
//---------------------------------------------------------------------------------------------------------------------------
  /** @brief Rayon minimal des particules.
  @return Le plus petit rayon. */
  Scalar getRayonMin();

  /** @brief Rayon d'interaction maximal des particules.
  @return Le plus grand rayon */
  Scalar getRayonInteractionMax();

  /** @brief Rayon d'interaction minimal des particules.
  @return Le plus petit rayon. */
  Scalar getRayonInteractionMin();

  /** @brief Volume des particules.
  @return Le volume global des particules : inserees et en attente. */
  Scalar getVolume() const;

  /** @brief Volume des particules inserees.*/
  Scalar getVolumeIn() const;

  /** @brief Volume des particules en attente. */
  Scalar getVolumeOut() const;

  /** @brief Nombre de particules totales */
  size_t nbreParticules() const { return ( m_particulesActives.size()
  	+ m_pwait.size() ); }

  /** @brief Nombre de particules actives */
  size_t nbreParticulesActives() const { return m_particulesActives.size(); }

  /** @brief Nombre de particules actives de tag 0 ou 1 */
  size_t nbreParticulesActivesOnProc() const;

  /** @brief Nombre total de particules sur tous les procs */
  size_t nbreParticulesOnAllProc() const { return m_nb_total_particules; }

  /** @brief Nombre de particules en attente */
  size_t nbreParticulesWait() const { return m_pwait.size(); }

  /** @brief Numero maximum de particules */
  int numeroMaxParticules() const;

  /** @brief Liste des obstacles � transmettre au fluide */
//  list<MonObstacle*> getObstaclesToFluid() const;
  list<Obstacle*> getObstaclesToFluid() const;

  /** @brief Return Sauter Mean diameter */
  static double get_SauterMeanDiameter( );

  /** @brief Return sum_i(Xi*Yi) for polydispersity */
  static double get_sumXY( );
  //@}


  /**@name Methods Set */
  //@{
  /** @brief Nombre total de particules sur tous les procs
  @param nb_ nombre total de particules sur tous les procs */
  void setNbreParticulesOnAllProc( const size_t &nb_ );

  /** @brief Fr�quence de mise � jour du lien entre les obstacles
  et le LinkedCell
  @param updateFreq la frequence */
  void setObstaclesLinkedCellUpdateFrequency( int const &updateFreq );

  /** @brief Parametres de post-processing des efforts sur les obstacles
  @param root_ nom du repertoire de sortie
  @param freq_ frequence de sortie
  @param ObsNames liste de noms des obstacles */
  void setOutputObstaclesLoadParameters( string const& root_,
  	int const& freq_,
	list<string> const& ObsNames );

  /** @brief Initialisation de la vitesse translationnelle et rotationnelle des
  particules au pas de temps precedent en cas de restart (lecture dans les
  fichiers resultats generes par PeliGRIFF, solution temporaire)
  @param PelDirRes repertoire de resultats PeliGRIFF */
  void setVelocityAndVelocityDifferencePreviousTimeRestart(
  	string const& PelDirRes );
  double setVelocityAndVelocityDifferencePreviousTimeRestart_Basilisk(
  	string const& dirRes );

  /** @brief Mouvement aleatoire sur les particules actives
  @param coefTrans amplitude max du mouvement translationnel
  @param coefRot amplitude max du mouvement rotationnel	*/
  void setRandomMotion( double const& coefTrans, double const& coefRot );

  /** @brief Set all contact map entries to false in all particles
  and all elementary obstacles */
  void setAllContactMapToFalse();

  /** @brief Update all contact map entries in all particles
  and all elementary obstacles */
  void updateAllContactMaps();
  //@}


  /** @name Processus de restauration */
  //@{
  /** @brief Sauvegarde de l'etat */
  void saveState();

  /** @brief Restauration de l'etat */
  void restaureState();
  //@}

   /** @brief Calls other shrinking
      functions if shrinking choice is 1
      in insert.xml */
  void ShrinkingRate(Scalar CurrentTime);
  //@}


  /**@name Methods I/O */
  //@{
  /** @brief Lecture des composants pour Reload.
  @param fileSave Flux de lecture
  @param filename nom du fichier correspondant au flux de lecture
  @param new_reload_format reload avec le format de reload 2014 */
  void read( istream &fileSave, string const& filename,
  	bool const& new_reload_format = false,
    bool const& contact_history_storage = false );

  /** @brief Sauvegarde des composants pour Reload.
  @param fileSave Flux de sauvegarde
  @param filename nom du fichier correspondant au flux de sauvegarde */
  void write( ostream &fileSave, string const& filename ) const;

  /** @brief Operateur << */
  friend ostream& operator << ( ostream &f, const EnsComposant &EC );

  /** @brief Debug
  @param s message de debug */
  void debug( char *s );

  /** @brief Ajoute un post-processing writer � la liste
  @param ppw post-processing writer */
  void addPostProcessingWriter( PostProcessingWriter* ppw );

  /** @brief Numero de cycle initial
  @param cycle0 numero de cycle initial */
  void setInitialCycleNumber( const int& cycle0 );

  /** @brief Verifie que le post processing Paraview est actif, sinon le cree
  @param rang rang du processeur
  @param nprocs nombre de processeurs
  @param name_ nom des fichiers
  @param root_ racine du nom des fichiers
  @param isBinary ecriture en mode binaire */
  void checkParaviewPostProcessing( int const& rang,
  	int const& nprocs,
  	const string &name_,
	const string &root_,
  	const bool &isBinary );

  /** @brief Verifie que le post processing Matlab est actif, sinon le cree
  @param rang rang du processeur
  @param nprocs nombre de processeurs
  @param name_ nom des fichiers
  @param root_ racine du nom des fichiers
  @param isBinary ecriture en mode binaire */
  void checkMatlabPostProcessing( int const& rang,
  	int const& nprocs,
  	const string &name_,
	const string &root_,
  	const bool &isBinary );

  /** @brief Post-processing des efforts sur les obstacles
  @param temps Temps de sauvegarde
  @param dt pas de temps
  @param enforceOutput oblige l'ecriture
  @param increaseCounterOnly augmente simplement le compteur pour ecriture
  @param rang rang du processeur
  @param nprocs nombre de processeurs
  @param wrapper wrapper MPI */
  void outputObstaclesLoad( Scalar temps, Scalar dt, bool enforceOutput = false,
      bool increaseCounterOnly = false,
      int rang = 0, int nprocs = 1, MPIWrapperGrains const* wrapper = NULL );

  /** @brief Initialise les fichiers de sortie des efforts sur les obstacles
  @param rang Rang du processeur
  @param coupledFluid Avec ou sans fluide
  @param temps Temps physique */
  void initialiseOutputObstaclesLoadFiles( int rang = 0,
      bool coupledFluid = false, Scalar temps = 0. );
  //@}


private:
  /** @name Parameters */
  //@{
  vector<Particule*> m_ParticuleClassesReference; /**< particules de r�f�rence
  	pour chaque classe de particules */
  vector<int> m_nbParticlesPerClass; /**< Number of particles per class (used
  for demcfd polydisperse) */
  static double SauterMeanDiameter; /** Sauter Mean Diameter (used for demcfd
  polydisperse) \<d\> = sum(ni di^3)/sum(ni di^2) */
  vector<double> ParticleClassesConcentration; /** Volumetric
  concentration of particles per class (used for demcfd polydisperse)
  xi = (ni di^3)/sum(ni di^3) */
  vector<double> diameterRatio; /** Contains Di / SauterMeanDiameter */
  static double sumXY; /**  sum_i(Xi*Yi) for polydispersity */
  list<Particule*> m_pwait; /**< Particules inactives dans la simulation.
      Elle sont supprimees du probleme ou/et en attente d'insertion. */
  Particule* m_wait; /**< Particule de la liste inactive en attente
  	d'insertion.*/
  list<Particule*> m_particulesActives; /**< Particules actives dans la
  	simulation. */
  list<Particule*> m_particulesHalozone; /**< Particules actives dans la zone de
  	recouvrement */
  list<Particule*> m_particulesClones; /**< Particules actives representant des
  	clones de particules sur un autre processeur (sous domaine) */
  list<Particule*> m_particulesClonesPeriodiques; /**< Particules actives
  	representant des clones periodiques de particules */
  set<Particule*> m_particulesReferencesPeriodiques; /**< Particules actives
  	possedant des clones periodiques */
  size_t m_nb_total_particules; /**< nombre total de particules dans le syst�me
  	(actives ou en attente) i.e. sur tous les processeurs */
  Obstacle *m_obstacle; /**< Obstacles decrivant le probleme */
  list<PostProcessingWriter*> m_postProcessors; /**< liste des
  	post-processors */
  bool m_hasSerialPostProcessors; /**< un des post-processors est il
  	sequentiel ? */
  list<ObstacleChargement*> m_ChargementsCinematiques; /**< liste de l'ensemble
  	des chargements cinematiques */
  list<ObstacleChargement_F*> m_ChargementsForces; /**< liste de
  	l'ensemble des chargements en force */
  list<Obstacle*> m_outputTorseurObstacles; /**< liste des obstacles pour
  	lesquels les efforts sont ecrits dans un fichier */
  string m_outputTorseurObstacles_dir; /**< nom du repertoire dans lequel sont
  	stockes les fichiers de sortie des efforts sur les obstacles */
  int m_outputTorseurObstacles_counter; /**< compteur pour la sortie des
  	efforts sur les obstacles */
  int m_outputTorseurObstacles_frequency; /**< frequence de sortie des
  	efforts sur les obstacles */
  Scalar m_initial_time; /** Initial time used in the integration of the contact
        force */
  //@}
};


/** @name Useful tools */
//@{
/** @brief Supprime la 1ere instance d'une valeur, si elle existe, dans une
liste de pointeurs de particule
@param pointerslist liste de pointeurs de particule
@param value la valeur � supprimer */
bool removeParticuleFromList( list<Particule*> &pointerslist,
	Particule* value );

/** @brief Supprime l'instance d'une valeur, si elle existe, dans un
set de pointeurs de particule
@param pointersSet liste de pointeurs de particule
@param value la valeur � supprimer */
bool removeParticuleFromSet( set<Particule*> &pointersSet, Particule* value );

/** @brief Supprime la 1ere instance d'une valeur, si elle existe, dans une
liste de pointeurs d'obstacle
@param pointerslist liste de pointeurs d'obstacle
@param value la valeur � supprimer */
bool removeObstacleFromList( list<MonObstacle*> &pointerslist,
	MonObstacle* value );
//@}

#endif
