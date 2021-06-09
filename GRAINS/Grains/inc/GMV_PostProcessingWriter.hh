#ifndef _GMV_PostProcessingWriter__
#define _GMV_PostProcessingWriter__

#include "PostProcessingWriter.hh"
#include <iostream>
#include <fstream>
using std::ofstream;

/** @brief Classe GMV_PostProcessingWriter.

    Ecriture des resultats pour post-processing par GMV.

    @author A.WACHS - Institut Francais du Petrole - 2010 - Creation */
//=============================================================================
class GMV_PostProcessingWriter : public PostProcessingWriter
{
public:
  /** @name Constructors */
  //@{
  /** @brief Constructeur par defaut 
  @param dn noeud XML 
  @param rank_ rang du processus 
  @param nbranks_ nombre total de processus */
  GMV_PostProcessingWriter(DOMNode* dn,int const& rank_,int const& nbranks_);

  /** @brief Destructeur */
  virtual ~GMV_PostProcessingWriter();
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
  
  /** @brief Numero de cycle initial 
  @param cycle0 numero de cycle initial */
  void setInitialCycleNumber(const int& cycle0) {counter = cycle0;};
  
  /** @brief Type de post-processing writer */
  string getPostProcessingWriterType() const {return "Gmv";};
  //@}


private:
  /** @name Methods */
  //@{
  /** @brief Ecriture
  @param temps Temps de sauvegarde  
  @param particules particules actives
  @param obstacle obstacles */
  void one_output(Scalar const& temps,
  	list<Particule*> const* particules,
  	Obstacle *obstacle);  
  //@}
  
  
private:
  string GMVFilename; /**< racine du nom de fichier de sortie GMV */
  int counter; /**< compteur pour les fichiers de sortie */  
  

protected:
  /** @name Constructors */
  //@{
  /** @brief Constructeur par defaut */
  GMV_PostProcessingWriter() {};
  //@}  
};

#endif
  
