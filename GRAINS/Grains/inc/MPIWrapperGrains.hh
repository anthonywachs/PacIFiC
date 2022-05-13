#ifndef _MPIWrapperGrains
#define _MPIWrapperGrains

#include <mpi.h>
#include "Voisins.hh"
#include "solvercomputingtime.hh"
#include "computingtime.hh"
#include "Particule.H"
#include "CompParticule.hh"
#include "LinkedCell.H"
#include "Point.H"
#include "Vecteur.H"
#include "Matrix.H"
#include "AppSec.H"  
#include <fstream>
#include <list>
#include <map>
using namespace std;


/** @brief Gestion des communications MPI
       
    @author L. GIROLAMI - Institut Francais du Petrole - 2009 - Creation */
// ============================================================================
class MPIWrapperGrains : public SolverComputingTime
{
private:
  /** @name Constructors & Destructor */
  //@{
  /** @brief Constructeur */
  MPIWrapperGrains();
  //@}

public:
  /** @name Constructors & Destructor */
  //@{
  /** @brief Constructeur 
  @param NX nombre de domaines dans la direction X 
  @param NY nombre de domaines dans la direction Y   
  @param NZ nombre de domaines dans la direction Z 
  @param PERX periodicite en x traite par MPI 
  @param PERY periodicite en y traite par MPI  
  @param PERZ periodicite en z traite par MPI */
  MPIWrapperGrains( int NX, int NY, int NZ,
  	int PERX, int PERY, int PERZ );
  
  /** @brief Destructeur. */
  ~MPIWrapperGrains();
  //@}

  
  /** @name Methodes d'accès */
  //@{
  /** @brief Nombre de processus dans la direction i
  @param i direction */
  int nb_procs_direction(int i) const;
  
  /** @brief Nombre de processus dans chaque direction */
  int const* nb_procs_direction() const;   
  
  /** @brief Coordonnees du processus */
  int const* MPI_coordonnees() const;
  
  /** @brief Voisins du processus */
  Voisins const* MPI_voisins() const;
  
  /** @brief Periodicite de la topologie MPI */
  int const* MPI_periodicite() const;  
  
  /** @brief Nombre total de processus dans le communicateur MPI_COMM_WORLD */
  int nombre_total_procs() const;
  
  /** @brief Nombre total de processus dans le communicateur MPI_COMM_WORLD */
  static int nombreTotalProcs(); 
  
  /** @brief Nombre de processus dans le communicateur MPI_COMM_activProc */
  int nombre_total_procs_ACTIV() const;   
  
  /** @brief Rang de processus dans le communicateur MPI_COMM_WORLD */
  int rank_WORLD() const; 
  
  /** @brief Rang de processus dans le communicateur MPI_COMM_activProc */
  int rank_ACTIV() const;   
  
  /** @brief Le processus est il actif ? */
  bool isActiv() const;   
  
  /** @brief Rang de processus dans le communicateur MPI_COMM_WORLD */
  static int rankOf_WORLD();                
  //@}  


  /** @name Methodes MPI pour les particules */
  //@{
  /** @brief Création & mise à jour des clones, en utilisant une stratégie de
  comm MPI de type AllGatherv sur tous les procs 
  @param time temps physique 
  @param particulesClones clones
  @param particules particules actives 
  @param particulesHalozone particules dans la zone de recouvrement
  @param ParticuleClassesReference particules de référence
  @param LC grille de cellules */
  void UpdateOrCreateClones_AllGatherGlobal( Scalar time,
  	list<Particule*>* particulesClones,
	list<Particule*>* particules,
  	list<Particule*> const* particulesHalozone,
	vector<Particule*> const* ParticuleClassesReference,
	LinkedCell* LC );

  /** @brief Création & mise à jour des clones, en utilisant une stratégie de
  comm MPI de type AllGatherv sur les voisins 
  @param time temps physique 
  @param particulesClones clones
  @param particules particules actives 
  @param particulesHalozone particules dans la zone de recouvrement
  @param ParticuleClassesReference particules de référence
  @param LC grille de cellules */
  void UpdateOrCreateClones_AllGatherLocal( Scalar time,
  	list<Particule*>* particulesClones,
	list<Particule*>* particules,
  	list<Particule*> const* particulesHalozone,
	vector<Particule*> const* ParticuleClassesReference,
	LinkedCell* LC );

  /** @brief Création & mise à jour des clones, en utilisant une stratégie de
  comm MPI de type Send-Recv entre voisins et en traitant les données recues à
  chaque réception en provenance d'un processeur
  @param time temps physique 
  @param particulesClones clones
  @param particules particules actives 
  @param particulesHalozone particules dans la zone de recouvrement
  @param ParticuleClassesReference particules de référence
  @param LC grille de cellules */
  void UpdateOrCreateClones_SendRecvLocal( Scalar time,
  	list<Particule*>* particulesClones,
	list<Particule*>* particules,
  	list<Particule*> const* particulesHalozone,
	vector<Particule*> const* ParticuleClassesReference,
	LinkedCell* LC );

  /** @brief Création & mise à jour des clones, en utilisant une stratégie de
  comm MPI de type Send-Recv entre voisins, en envoyant que les donnees utiles
  aux voisins en terme de particules partagees et en traitant les données 
  recues à chaque réception en provenance d'un processeur
  @param time temps physique 
  @param particulesClones clones
  @param particules particules actives 
  @param particulesHalozone particules dans la zone de recouvrement
  @param ParticuleClassesReference particules de référence
  @param LC grille de cellules */
  void UpdateOrCreateClones_SendRecvLocal_GeoLoc( Scalar time,  	
  	list<Particule*>* particulesClones,
	list<Particule*>* particules,
  	list<Particule*> const* particulesHalozone,
	vector<Particule*> const* ParticuleClassesReference,
	LinkedCell* LC );
	
  /** @brief Collecte sur le processeur master de l'ensemble des particules 
  sur les différents processeurs pour post-processing
  @param particules liste des particules actives dans le système
  @param pwait liste des particules en attente dans le système  
  @param ParticuleClassesReference particules de référence
  @param nb_total_particules nombre total de particules sur l'ensemble des
  processeurs */	
  vector<Particule*>* GatherParticules_PostProcessing(
  	const list<Particule*>& particules,
	const list<Particule*>& pwait,
	vector<Particule*> const& ParticuleClassesReference,
	const size_t& nb_total_particules ) const;	

  /** @brief Collecte sur le processeur master la cinématique de l'ensemble 
  des particules sur les différents processeurs pour post-processing
  @param particules liste des particules actives dans le système
  @param nb_total_particules nombre total de particules sur l'ensemble des
  processeurs */	
  vector< vector<double> >* GatherPositionVitesse_PostProcessing(
  	const list<Particule*>& particules,
	const size_t& nb_total_particules ) const;	

  /** @brief Collecte sur le masterProc de la classe des particules
  de tous les coeurs.
  @param particules liste des particules actives dans le système
  @param nb_total_particules nombre total de particules sur l'ensemble des
  processeurs */	
  vector< vector<double> >* GatherParticlesClass_PostProcessing(
  	const list<Particule*>& particules,
	const size_t& nb_total_particules ) const;	

  /** @brief Collecte sur le processeur master de l'ensemble des clones
  periodiques sur les différents processeurs pour post-processing
  @param particulesClonesPeriodiques particules clones periodiques
  @param ParticuleClassesReference particules de référence */	
  list<Particule*>* GatherClonesPeriodiques_PostProcessing(
  	const list<Particule*> &particulesClonesPeriodiques,
	vector<Particule*> const& ParticuleClassesReference ) const;

  /** @brief Création & mise à jour des nouveaux clones periodiques par comm 
  des particules periodiques de reference, en 
  utilisant une stratégie de comm MPI de type AllGatherv sur tous les procs du 
  communicateur periodique
  @param time temps physique 
  @param particulesReferencesPeriodiques particules periodiques de reference
  @param particulesClonesPeriodiques particules clones periodiques
  @param ParticuleClassesReference particules de référence
  @param LC grille de cellules */
  void commParticulesRefPer_AllGatherGlobal_UpdateOrCreate( Scalar time,
  	set<Particule*>* particulesReferencesPeriodiques,
  	list<Particule*>* particulesClonesPeriodiques,
	vector<Particule*> const* ParticuleClassesReference,
	LinkedCell* LC );
	
  /** @brief AllGather de listes d'entiers en utilisant une stratégie de comm 
  MPI de type AllGatherv sur tous les procs du communicateur periodique; 
  Utilisation principale: comm des listes des ID des particules periodiques de 
  reference pour lesquelles le statut des clones change (destruction ou 
  transformation en particules classiques)
  @param time temps physique 
  @param IDs liste d'entiers */
  void commPer_AllGatherGlobal_listINT( Scalar time,
  	list<int> &IDs );
	
  /** @brief AllGather de listes d'entiers en utilisant une stratégie de comm 
  MPI de type AllGatherv sur tous les procs du communicateur periodique; 
  Utilisation principale: comm des listes des ID des particules periodiques de 
  reference pour lesquelles le statut des clones change (destruction ou 
  transformation en particules classiques). En sortie la liste est filtree pour
  qu'une seule occurence de chaque valeur ne subsiste.
  @param time temps physique 
  @param IDs liste d'entiers */
  void commPer_AllGatherGlobal_listINT_unique( Scalar time,
  	list<int> &IDs );			
	
  /** @brief Envoi des forces de contact des clones periodiques a leur 
  reference, en utilisant une stratégie de comm MPI de type AllGatherv sur tous 
  les procs du communicateur periodique
  @param time temps physique 
  @param particulesReferencesPeriodiques particules periodiques de reference
  @param particulesClonesPeriodiques particules clones periodiques */
  void commParticulesClonesPer_AllGatherGlobal_Forces( Scalar time,
  	set<Particule*>* particulesReferencesPeriodiques,
  	list<Particule*>* particulesClonesPeriodiques ) const;		
		
  /** @brief Broadcast un entier à partir du master
  @param i entier */
  int Broadcast_INT( const int &i ) const;

  /** @brief Broadcast un double à partir du master
  @param i double */
  double Broadcast_DOUBLE( const double &i ) const;
  
  /** @brief Broadcast un entier non signe à partir du master
  @param i entier */
  size_t Broadcast_UNSIGNED_INT( const size_t &i ) const;  
  
  /** @brief Somme un entier sur tous les proc
  @param i entier */
  int sum_INT( int i ) const;
  
  /** @brief Somme un entier non signe sur tous les proc
  @param i entier */
  size_t sum_UNSIGNED_INT( size_t i ) const;  
  
  /** @brief Somme un entier de tous les proc sur le master
  @param i entier */
  int sum_INT_master( int i ) const; 
  
  /** @brief Somme un entier non signe de tous les proc sur le master
  @param i entier */
  size_t sum_UNSIGNED_INT_master( size_t i ) const;     

  /** @brief perform a "logical and" operation on input boolean value
  @param input boolean value on which we want to perform the operation */
  bool logical_and( bool input ) const;     

  /** @brief Somme d'un double sur tous les proc
  @param x double */
  double sum_DOUBLE( double x ) const;
  
  /** @brief Somme d'un double de tous les proc sur le master
  @param x double */
  double sum_DOUBLE_master( double x ) const;     
  
  /** @brief Min d'un entier sur tous les proc
  @param i entier */
  int min_INT( int i ) const;
  
  /** @brief Min d'un entier non signe sur tous les proc
  @param i entier */
  size_t min_UNSIGNED_INT( size_t i ) const;  
  
  /** @brief Max d'un entier sur tous les proc
  @param i entier */
  int max_INT( int i ) const;
  
  /** @brief Max d'un entier non signe sur tous les proc
  @param i entier */
  size_t max_UNSIGNED_INT( size_t i ) const;   
  
  /** @brief Max d'un double sur tous les proc
  @param x double */
  double max_DOUBLE( double x ) const;
  
  /** @brief Max d'un double de tous les proc sur le master
  @param x double */
  double max_DOUBLE_master( double x ) const;  
  
  /** @brief Min d'un double sur tous les proc
  @param x double */
  double min_DOUBLE( double x ) const;
  
  /** @brief Min d'un double de tous les proc sur le master
  @param x double */
  double min_DOUBLE_master( double x ) const;
  
  /** @brief AllGather d'un entier non signe
  @param i entier */
  size_t* AllGather_UNSIGNED_INT( size_t i ) const;          
  
  /**
    @brief Broadcast un point a partir du master
    @param pt point
  */
  Point Broadcast_Point( const Point &pt ) const;

  /**
    @brief Broadcast un vecteur a partir du master
    @param pt Vecteur
  */
  Vecteur Broadcast_Vecteur( const Vecteur &v ) const;    

  /** @brief Broadcast une matrice a partir du master
  @param mat matrice */
  Matrix Broadcast_Matrix( const Matrix &mat ) const;  
  
  /** @brief Test: AllGather d'un vecteur de n entiers
  @param n nb d'elements du vecteur */
  void test_AllGatherv_INT( const int &n ) const;
  
  /** @brief Test avec communicateur local: AllGather d'un vecteur de n entiers
  @param n nb d'elements du vecteur */
  void testCommLocal_AllGatherv_INT( const int &n ) const;  

  /** @brief Test: AllGather d'un vecteur de n doubles
  @param n nb d'elements du vecteur */
  void test_AllGatherv_DOUBLE( const int &n ) const;
  
  /** @brief Test avec communicateur local: AllGather d'un vecteur de n doubles
  @param n nb d'elements du vecteur */
  void testCommLocal_AllGatherv_DOUBLE( const int &n ) const;  

  /** @brief Test avec communicateur local: Send-Recv d'un vecteur de n doubles
  @param n nb d'elements du vecteur */
  void testCommLocal_SendRecv_DOUBLE( const int &n ) const;  

  /** @brief Test avec communicateur local: validation */
  void testCommLocal() const; 
  
  /** @brief Definition des communicateurs locaux lies aux voisins */
  void setCommLocal(); 
  
  /** @brief Bilan timer */
  void bilanTimer() const;
  
  /** @brief Affiche la chaine de log des comm MPI par proc et la reinitialise 
  @param f flux de sortie */
  void writeAndFlushMPIString( ostream &f );
  
  /** @brief MPI_Barrier pour les processus actifs uniquement */
  void MPI_Barrier_ActivProc() const; 
  
  /** @brief Definition du communicateur pour les simulations periodiques 
  @param hasPeriodicObstacle indique si le sous-domaine du processeur contient 
  	des obstacles periodiques */
  void setCommPeriodic( bool const& hasPeriodicObstacle ); 

  /** @brief Definition des vecteurs de periodicite MPI dans le communicateur
  standard commgrainsMPI_3D
  @param lx longueur du domaine global en x
  @param ly longueur du domaine global en y  
  @param lz longueur du domaine global en z */
  void setMPIperiodicVectors( const Scalar& lx, const Scalar& ly, 
  	const Scalar& lz ); 
  
  /** @brief Caracteristiques globales des contacts 
  @param overlap_max penetration maximale
  @param overlap_mean penetration moyenne    
  @param time_overlapMax instant de l'overlap max 
  @param nbIterGJK_mean nombre d'iterations moyen de GJK */
  void ContactsFeatures( Scalar& overlap_max,
	Scalar& overlap_mean,
	Scalar& time_overlapMax,
	Scalar& nbIterGJK_mean ) const;
	
  /** @brief Somme des efforst sur les obstacles sur le master 
  @param allMyObs liste des obstacles elementaires */
  void sumObstaclesLoad( list<MonObstacle*> const& allMyObs ) const;
  
  /** @brief Repartition des nombres de particules par proc et par classe 
  dans l'insertion par bloc structure
  @param newPart nombre total de particules par classe
  @param newPartProc nombre de particules par classe sur le proc
  @param npartproc nombre total de particules a insere sur le proc */
  void distributeParticulesClassProc( 
  	list< pair<Particule*,int> > const& newPart,
	list< pair<Particule*,int> >&newPartProc,
	size_t const& npartproc,
	size_t const& ntotalinsert ) const; 	
  //@}  


  /** @name Methodes static */
  //@{
  /** @brief Complete la chaine contenant les comm MPI
  @param add chaine a ajouter a la chaine de comm MPI */
  static void addToMPIString( const string &add ); 
  
  /** @brief Renvoie la MPIGeoLocalisation en fonction de la position relative
  @param i index dans la direction X (-1,0 ou 1) 
  @param j index dans la direction Y (-1,0 ou 1)   
  @param k index dans la direction Z (-1,0 ou 1) */
  static MPIGeoLocalisation getMPIGeoLocalisation( int i, int j, int k );
  //@}  


  /** @name Methodes i/o */
  //@{  
  /** @brief Affichage des caracteristiques */
  void display( ostream &f ) const;

  /** @brief Ecriture de la memoire utilisee par la simulation sur chaque 
  processeur  
  @param f output flux */
  void display_used_memory( ostream &f ) const;  
  //@}
    

private:
  /** @name Parameters */
  //@{
  MPI_Group m_MPI_GROUP_activProc; /**< groupe des processus actifs */    
  MPI_Comm m_MPI_COMM_activProc; /**< communicateur des processus actifs */
  int *m_coords; /**< coordonnees dans la topologie cartesienne MPI */
  int *m_dim; /**< nombre de domaines dans chaque direction de topologie 
  	cartesienne MPI */
  int *m_period; /**< periodicite dans chaque direction */
  bool m_isMPIperiodic; /**< au moins une direction du pattern MPI est 
  	periodique */
  int m_rank; /**< rang du processus dans le communicateur MPI_COMM_activProc */
  int m_rank_world; /**< rang du processus dans le communicateur
  	MPI_COMM_WORLD */
  int m_rank_masterWorld; /**< rang dans le communicateur MPI_COMM_activProc du
  	processus ayant le rang 0 (master) dans le communicateur 
	MPI_COMM_WORLD */  
  int m_nprocs; /**< nombre de processus actifs */
  int m_nprocs_world; /**< nombre total de processus */
  bool m_is_activ; /**< le processus est il actif ? */  
  Voisins *m_voisins; /**< voisins du processus */
  MPI_Comm *m_commgrainsMPI_3D; /**< communicateur lie a la topologie 
  	cartesienne */
  vector<MPI_Group*> m_groupVoisins; /**< groupes lies aux voisins dans la 
  	topologie cartesienne */	
  vector<MPI_Comm*> m_commVoisins; /**< communicateurs lies aux voisins dans la 
  	topologie cartesienne */
  vector<bool> m_isInCommVoisins; /**< communicateurs dans lequel le processus
  	courant est présent */
  vector<int> m_masterGeoLoc; /**< geolocalisation du master dans les 
  	communicateurs lies aux voisins auquel appartient le proc */	
  int m_nprocs_localComm; /**< nb de processus dans le communicateur lie aux
 	voisins */
  int m_rank_localComm; /**< rang du processus dans le communicateur lie aux
 	voisins */	
  int *m_master_localComm; /**< numero du master dans le communicateur lie 
  	aux voisins */
  static string *m_MPILogString; /**< chaine de log des comm MPI */
  static vector< vector<int> > m_particuleHalozoneToNeighboringProcs; /**< 
  	correspondance entre la geolocalisation d'une particule dans la zone de
	halo et les procs à qui les infos de cette particule doivent etre
	envoyees */
  static vector<int> m_GeoLocReciprocity; /**< correspondance reciproque des
  	geolocalisations (ex: MPIGEO_BOTTOM -> MPIGEO_TOP ) */
  MPI_Group* m_MPI_GROUP_Periodic; /**< groupe des processus contenant des
  	obstacles periodiques */    
  MPI_Comm* m_MPI_COMM_Periodic; /**< communicateur des processus contenant des
  	obstacles periodiques */
  bool m_isInCommPeriodic; /**< le processus courant est il présent dans le
  	comunicateur pour la simulation periodique */ 
  int m_nprocs_per; /**< nombre de processus du communicateur periodique */
  vector<Vecteur> m_MPIperiodes; /**< periodes associees a la topologie MPI */
  multimap<int,Particule*> AccessToClones; /**< facilite l'acces aux clones 
  	par numero */
  //@}


  /** @name Methods  */
  //@{
  /** @brief Definition de particuleHalozoneToNeighboringProcs liee a la
  numerotation dans l'enumeration MPIGeoLocalisation */
  void setParticuleHalozoneToNeighboringProcs(); 
  
  /** @brief Definition de la geolocalisation du master dans les communicateurs
  lies aux voisins auquel appartient le proc */
  void setMasterGeoLocInLocalComm(); 
  
  /** @brief Création & mise à jour des clones sur la base des infos
  communiquees par les autres proc 
  @param time temps physique 
  @param recvsize nombre de particules communiquees
  @param recvbuf_DOUBLE tableau de doubles pour les particules communiquees
  @param NB_DOUBLE_PART nombre de doubles communiques par particule        
  @param particulesClones clones
  @param particules particules actives 
  @param particulesHalozone particules dans la zone de recouvrement
  @param ParticuleClassesReference particules de référence
  @param LC grille de cellules */
  void UpdateOrCreateClones( Scalar time,
 	int const &recvsize, double const* recvbuf_DOUBLE,
	const int& NB_DOUBLE_PART, 
  	list<Particule*>* particulesClones,
	list<Particule*>* particules,
  	list<Particule*> const* particulesHalozone,
	vector<Particule*> const* ParticuleClassesReference,
	LinkedCell* LC );               
  //@}  
};

#endif
