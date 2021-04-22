#ifndef _STL_PostProcessingWriter__
#define _STL_PostProcessingWriter__

#include "PostProcessingWriter.hh"
#include <iostream>
#include <fstream>
using std::ofstream;

/** @brief Classe STL_PostProcessingWriter.

    Fichier position des particles au format STL pour link avec OpenFoam.
    Seul l'etat final est ecrit (a chaque sauvegarde, le fichier est ecrase). 

    @author A.WACHS - Institut Francais du Petrole - 2011 - Creation */
//=============================================================================
class STL_PostProcessingWriter : public PostProcessingWriter
{
public:
  /** @name Constructors */
  //@{
  /** Constructeur par defaut 
  @param dn noeud XML 
  @param rank_ rang du processus 
  @param nbranks_ nombre total de processus */
  STL_PostProcessingWriter(DOMNode* dn,int const& rank_,int const& nbranks_);

  /** Destructeur */
  virtual ~STL_PostProcessingWriter();
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
  @param LC grille de cellules  */
  void PostProcessing(Scalar const& temps,
  	Scalar const& dt,
  	list<Particule*> const* particules,
	list<Particule*> const* pwait,
	list<Particule*> const* pperiodiques,		
	vector<Particule*> const* ParticuleClassesReference,
	Obstacle *obstacle,
	LinkedCell const* LC);

  /** @brief Clot les ecritures */
  void PostProcessing_end();

  /** @brief Type de post-processing writer */
  string getPostProcessingWriterType() const {return "STL";};
  //@}
  

private:
  string simul; /**< racine des fichiers de sortie */

  /** @name Methods */
  //@{  
  /** @brief Ecriture
  @param temps Temps de sauvegarde  
  @param particules particules actives */
  void one_output(Scalar const& temps,
  	list<Particule*> const* particules); 		
  //@}



protected:
  /** @name Constructors */
  //@{
  /** @brief Constructeur par defaut */
  STL_PostProcessingWriter() {};
  //@}  
};

#endif
  
