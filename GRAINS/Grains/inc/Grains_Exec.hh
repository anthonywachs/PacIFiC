#ifndef _Grains_Exec__
#define _Grains_Exec__

#include <mpi.h>
#include "Vecteur.H"
#include "Point.H"
#include "VertexBase.H"
#include "IndexArray.H"
#include <string>
#include <fstream>
#include <list>
#include <vector>
#include <set>
using std::string;
using std::ios_base;
using std::list;
using std::vector;
using std::set;
using namespace solid;

class MPIWrapperGrains;
class App;


/** @brief Type de fenetre
@author GRAINS Project - IFPEN - 2015 */
enum FenetreType 
{
  FENETRE_BOX,
  FENETRE_CYLINDER,
  FENETRE_ANNULUS,
  FENETRE_LINE
};


struct Fenetre
{
  /** Type de fenetre: boite ou cylindre */
  FenetreType ftype;
  /** Premier point de la fenetre ou centre du disque inferieur du cylindre */
  Point ptA;
  /** Second point de la fenetre */
  Point ptB;
  /** Rayon du cylindre */
  Scalar radius;
  /** Rayon du cylindre interieur dans le cas d'une fenetre annulaire */
  Scalar radius_int;  
  /** Hauteur du cylindre */
  Scalar hauteur; 
  /** Direction de l'axe du cylindre */
  Direction axisdir; 
};


struct BlocInsertion
{
  /** Fenetre */
  struct Fenetre bloc;
  /** Nb de positions en X */
  size_t NX;
  /** Nb de positions en Y */
  size_t NY;
  /** Nb de positions en Z */
  size_t NZ;    
};


/** @brief Classe d'acces a des parametres globaux
       
    @author A. WACHS - Institut Francais du Petrole - 2010 - Creation */
// ============================================================================
class Grains_Exec
{
public:
  /** @name Methodes d'accès */
  //@{
  /** @brief Acces au wrapper MPI de Grains */
  static MPIWrapperGrains* getComm();

  /**
    @brief Accessor to the list of applications
    @return List of pointer to App
  */
  static list<App*> get_listApp();
  //@}

  
  /** @name Methodes Set */
  //@{
  /**
    @brief Wrapper MPI de Grains
    @param wrapper_ wrapper
  */
  static void setComm(MPIWrapperGrains* wrapper_);

  /**
    @brief Set the list of applications
    @param List of pointer to App
  */
  static void set_listApp( list<App*> allApp_ );
  //@}

  /** @name Gestion des exceptions en parallele */
  //@{  
  static bool m_MPI; /**< Simulation MPI ou sequentielle */  
  static bool m_exception_Contact; /**< Exception sur le contact */
  static bool m_exception_Deplacement; /**< Exception sur le deplacement */
  static bool m_exception_Simulation; /**< Exception de type Simulation */    
  //@} 
  
  /** @name Autres parametres globaux */
  //@{  
  static string m_TIScheme; /**< type de schema d'integration en temps */ 
  static bool m_SphereAsPolyParaview; /**< dans les sorties Paraview, vrai si
  	les spheres sont facetisees avec des polyedres, faux si elles sont
  	post-traitees avec un champ vectoriel */ 
  static int m_MPI_verbose; /**< Ecriture detaillee a l'ecran des operations
  	MPI; 3 levels: 0=rien, 1=particules, 2=particules+contact */
  static string m_ReloadType; /**< Type de reload: "new" pour une nouvelle
  	simulation et "same" pour poursuivre la même simulation */
  static Vecteur m_vgravite; /**< vecteur gravite */
  static Vecteur* m_translationParaviewPostProcessing; /** Dans le cas d'un
  	calcul couple avec le fluide en translation-projection, possibilite de
	translater les sorties Paraview des particules de la distance opposee de
	translation, si cette fonctionnalite est inactive, le pointeur vaut 
	NULL */
  static bool m_withdemcfd; /** Euler-Lagrange coupling ? */
//  static bool m_slipveloutput; /** make output of slip velocity ? */
  static bool m_ContactDissipation; /** With contact dissipation rate output ? */
  static bool m_withlubrication; /** With lubrication correction ? */
  static bool m_ContactforceOutput; /** Contact force on each particle
  	 as text output file */
  static bool m_ContactforceOutput_instantaneous; /** Instantaneous contact force on
          each particle as text output file */
  static bool m_withHydroForce; /** Hydro (drag+potential lift) force on particles
      ? */
  static bool m_withLiftForce; /** Lift Force in Euler-Lagrange coupling ? */
  static bool m_stressTensor; /** Outpout stress tensor? */
  static bool m_particleStressTensor; /** Outpout of individual stress tensor? */
  static bool m_withFluidTemperature; /** Heat transfer between fluid and solid
      in DEM-CFD ? */
  static bool m_withSolidTemperature; /** Heat transfert between particles ? */
  static bool m_withPressureGradient; /** add gradP to total hydro force 
    in Euler-Lagrange coupling ? */
  static bool m_withStochasticDrag; /** the stochastic drag force ? */
  static bool m_withStochasticNusselt; /** the stochastic Nusselt number ? */
  static bool m_addedmass_demcfd; /**< with added mass in demcfd ? */  
  static double m_Prandtl; /** Prandtl number in case of temperature coupling */
  static bool m_periodique; /** Simulation de type periodique ? */
  static bool m_MPIperiodique; /** Simulation de type periodique ou
  	la periodicite est geree par la pattern MPI */
  static bool m_isGrainsCompFeatures; /** Passage dans CompFeatures */
  static bool m_isGrainsPorosity; /** Calcul de la porosite */
  static string m_ReloadDirectory; /** repertoire ou sont stockes tous les
  	fichiers de reload */
  static string m_SaveDirectory; /** repertoire ou sont ecrits tous les
  	fichiers de reload */
  static set<string> m_additionalDataFiles; /** fichiers autres que le fichier
  	de mise en donnees pour decrire le cas test (fichiers polyedres et
	polygones essentiellement) */
  static bool m_writingModeHybrid; /** est ce que le mode d'ecriture des
  	fichiers de reload est hydride i.e. un fichier texte pour entete et 
	obstacles plus un fichier binaire pour les particules de la 
	simulation */
  static string m_GRAINS_HOME; /** nom du repertoire ou se situe le code
  	execute */
  static string m_reloadFile_suffix; /** suffixe du fichier de reload 
  	(A ou B) */
  static bool m_withCohesion; /** with cohesion? */
  static vector<Scalar> m_stressTensorDomain; /** Computation domain of the
    stress tensor */
  //@}   


  /** @name Utilitaires static */
  //@{  
  /** @brief Ecriture d'un double dans un format donné (nombre de caracteres)
  @param figure le double 
  @param size nombre de caractere pour ecrire le double */
  static string doubleToString( const double &figure, const int &size );
  
  /** @brief Ecriture d'un double dans un format donné (precision) 
  @param format le format ios
  @param digits nombre de chiffres significatifs
  @param number le double */
  static string doubleToString( ios_base::fmtflags format, int digits,
      	double const& number ); 
	
  /** @brief Ecriture d'un integer dans un string
  @param figure l'integer */
  static string intToString( const int &figure );
  
  /** @brief Verifie le temps de la derniere ecriture dans un fichier
  @param filename nom du fichier 
  @param current_time temps actuel */
  static void checkTime_outputFile( const string& filename, 
    	double const& current_time ); 
	
  /** @brief Recupere le chemin a partir d'un nom complet de fichier. Exemple:
  Titi/tutu/toto renvoie Titi/tutu
  @param FileName nom du fichier */
  static string extractRoot( string const& FileName ); 
  
  /** @brief Extrait le nom du fichier uniquement a partir d'un nom complet de 
  fichier. Exemple: Titi/tutu/toto renvoie toto
  @param FileName nom du fichier */
  static string extractFileName( string const& FileName );
  
  /** @brief Verifie que tous les fichiers de reload sont bien dans le
  meme repertoire de reload (fichiers polyedres et polygones essentiellement) */
  static void checkAllFilesForReload();   
  
  /** @brief Recupere le nom du fichier de restart a partir du fichier
  RFTable.txt
  @param rootName racine du fichier RFTable.txt 
  @param RFTable_ext extension du fichier RFTable.txt */
  static string restartFileName_AorB( string const& rootName, 
  	string const& RFTable_ext );     	  
  //@}   


  /** @name Garbage collector pour les particules */
  //@{    
  /** @brief Destruction de la memoire (garbage collector) */
  static void GarbageCollector(); 

  /** @brief Ajoute une liste de points de reference utilises par les 
  polytopes
  @param refPB liste de points de reference utilises par les polytopes 
  @param refVB objet contenant le pointeur sur la liste de points de reference
  et utilise pour acceder aux points */
  static void addOnePolytopeRefPointBase( Point* refPB, VertexBase *refVB );
  
  /** @brief Ajoute une description des points sommets pour un type de polytope
  @param idar description des points sommets */
  static void addOnePolytopeNodeNeighbors( IndexArray* idar );
  
  /** @brief Ajoute un tableau d'indice des sommets pour un type de polytope
  @param idar tableau d'indice des sommets */
  static void addOnePolytopeNodeIndex( IndexArray* idar );
  
  /** @brief Ajoute une connectivite des faces pour un type de polyedre
  @param faceCon connectivite des faces */
  static void addOnePolyhedronFaceConnectivity( 
  	vector< vector<int> >* faceCon );
  //@}  
  
  
  /** @name Nombre total de particules dans la simulation */
  //@{  
  /** @brief Nombre total de particules sur tous les procs 
  @param nb_ nombre total de particules sur tous les procs */
  static void setNbreParticulesOnAllProc( const size_t &nb_ ) 
  	{ m_nb_total_particules = nb_; }
  
  /** @brief Nombre total de particules sur tous les procs */
  static size_t nbreParticulesOnAllProc() { return m_nb_total_particules; }
  //@}   
  
  
  /** @name Memoire utilisee */
  //@{  
  /** @brief Memoire utilisee par la simulation sur ce processeur */
  static size_t used_memory( void );
  
  /** @brief Ecriture de la memoire utilisee par la simulation sur ce processeur
  @param os flux de sortie 
  @param memory memoire utilisee */
  static void display_memory( ostream& os, size_t memory );  	        
  //@} 
  
  
private:
  /** @name Parameters */
  //@{
  static MPIWrapperGrains* wrapper; /**< wrapper MPI de Grains */
  static list<App*> m_allApp; /**< Liste des applications utilisees dans la simulation;
      	la premiere application correspond au probleme sec */
  static string *debugString; /**< chaine pour sortie ecran de debug */ 
  static size_t m_nb_total_particules; /**< nombre total de particules dans le 
  	système (actives ou en attente) i.e. sur tous les processeurs */
  static list< pair<Point*,VertexBase *> > m_allPolytopeRefPointBase; /**< 
  	liste des points de reference utilises par les differents polytopes */
  static list<IndexArray*> m_allPolytopeNodeNeighbors; /**< liste 
  	des descriptions des points sommets dans les differents polytopes, les
	pointeurs sont ici des tableaux d'IndexArray */
  static list<IndexArray*> m_allPolytopeNodesIndex; /**< liste 
  	des tableaux d'indice des sommets dans les differents polytopes */
  static list<vector< vector<int> >*> m_allPolyhedronFacesConnectivity; /**< 
  	liste de connectivite des faces dans les differents polyedres */
  //@} 
  

  /** @name Constructors & Destructor */
  //@{
  /** @brief Constructeur */
  Grains_Exec();
  
  /** @brief Destructeur. */
  ~Grains_Exec();
  //@}       
};

#endif
