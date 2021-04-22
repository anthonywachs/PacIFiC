#ifndef _Grains_BuilderFactory
#define _Grains_BuilderFactory

#include "Grains.H"
#include "GrainsCoupledWithFluid.hh"
#include "ReaderXML.hh"


/** @brief Contexte des applications disponibles

    @author GRAINS Project - IFP - 2007 */
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
enum EAPPLI {
  /** @name Enumeration */
  //@{
  /** Simulation 3D */
  DIM_3,
  /** Simulation 2D */
  DIM_2,
  /** Simulation non definie */
  UNDEFINED
  //@}
};


/** @brief Fabrique des use case Grains.

    Construction de l'application en fonction du cas demande.

    @author G.FERRER - Institut Francais du Petrole - 2003 - Creation */
// ============================================================================
class Grains_BuilderFactory
{
public:
  /**@name Methods Static */
  //@{

  /** 
    @brief Creation de l'application avec l'option demandee
    @param root Le noeud "<Grains3D>" ou "<Graind2D>"
    @return L'application definie
    @remark Called by the main class
  */
  static Grains* create( DOMElement* root );
  
  /** @brief Creation de l'application GRAINS_INFLUID avec l'option demandee
  @param root Le noeud principal "<Grains3D>" ou "<Graind2D>"
  @param rhoFluide Densite du fluide
  @param grid_size size of the smallest grid cell
  @return L'application definie */
  static GrainsCoupledWithFluid* createCoupledWithFluid( DOMElement* root, 
	double rhoFluide, double grid_size );  

  /** @brief Contexte de la simulation
  @return Le contexte : 2D ou 3D */
  static EAPPLI getContext();
    
  /** @brief Ajoute le chemin vers les fichiers dtd en fonction de la variable
  shell GRAINS_HOME et du mode d'execution (dimension, avec ou sans 
  fluide, ... ), cree un fichier temporaire contenant ce chemin et renvoie le
  nom de ce fichier, c'est ce dernier qui sera lu par Grains. Cree le repertoire
  local des scripts bash lances par Grains pendant la simulation 
  @param filename nom du fichier 
  @param rank rang du processus 
  @param nprocs nb total de processus */
  static string init( string const& filename, const int& rank, 
  	const int& nprocs );  
  //@}


private:
  /**@name Constructors */
  //@{
  /** @brief Classe Statique : pas de contructeur disponible. */
  Grains_BuilderFactory() {};

  /** @brief Classe Statique : pas de destructeur disponible. */
  ~Grains_BuilderFactory() {};
  //@}


  /** @name Parameters */
  //@{
  static EAPPLI m_context; /**< Contexte de la simulation */
  //@}
};

#endif

