#ifndef _DLMFD_PostProcessingWriter__
#define _DLMFD_PostProcessingWriter__

#include "PostProcessingWriter.hh"
#include <iostream>
#include <fstream>
using std::ofstream;

/** @brief Classe DLMFD_PostProcessingWriter.

    PostProcessing writer for DLMFD solver. For the moment,
    it is created in order to make a text output of time-averaged
    contact force on each particle. 

    @author A.ESTEGHAMATIAN - IFPEN - 2015 - Creation */
//=============================================================================
class DLMFD_PostProcessingWriter : public PostProcessingWriter
{
public:
  /** @name Constructors */
  //@{
  /** Constructeur par defaut 
  @param dn noeud XML 
  @param rank_ rang du processus 
  @param nbranks_ nombre total de processus */
  DLMFD_PostProcessingWriter(DOMNode* dn,int const& rank_,int const& nbranks_);

  /** Destructeur */
  virtual ~DLMFD_PostProcessingWriter();
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
  string getPostProcessingWriterType() const {return "DLMFD";};

  /** @brief Le post-processing writer est il parallele ? */
  bool isParallelWriter() const { return false; }

  //@}
  /** @name Static variables */
  //@{
  static bool b_cumulatedContactForce; /**< if the sum of contact force on each particle
   	is asked as an output. */
  static bool b_cumulatedLubriForce; /**< if the sum of lubri force on each particle
   	is asked as an output. */
  //@}


private:
  /** @name Methods */
  //@{
	
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
  @param File opening mode (here : ios::app) */
  void prepareResultFiles( ios_base::openmode mode ) ;

  //@}
  

private:
  ofstream coordination_number; /**< Flux de nombre de contacts sur chaque
  	particule */ 
  ofstream contact_force_x; /**< Flux de coordonnees X de force de contact */ 
  ofstream contact_force_y; /**< Flux de coordonnees Y de force de contact */    
  ofstream contact_force_z; /**< Flux de coordonnees Z de force de contact */
  ofstream lubri_force_x; /**< Flux de coordonnees X de force de lubri */ 
  ofstream lubri_force_y; /**< Flux de coordonnees Y de force de lubri */    
  ofstream lubri_force_z; /**< Flux de coordonnees Z de force de lubri */
  string simul; /**< racine des fichiers de sortie */

protected:
  /** @name Constructors */
  //@{
  /** @brief Constructeur par defaut */
  DLMFD_PostProcessingWriter() {};  
  //@}  
};

#endif
  
