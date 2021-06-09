#ifndef _PostProcessingWriter__
#define _PostProcessingWriter__

#include "Basic.H"
#include <list>
#include <string>
using std::list;
using std::string;
#include "ReaderXML.hh"
#include "LinkedCell.H"

class Particule;
class Obstacle;
class Composant;


/** @brief Classe PostProcessingWriter.

    Ecriture des resultats pour post-processing.

    @author A.WACHS - Institut Francais du Petrole - 2010 - Creation */
//=============================================================================
class PostProcessingWriter
{
public:
  /** @name Constructors */
  //@{
  /** @brief Constructeur avec arguments
  @param rank_ rang du processus 
  @param nbranks_ nombre total de processus */
  PostProcessingWriter(int const& rank_,int const& nbranks_):
  	m_rank(rank_), m_nprocs(nbranks_) {};

  /** @brief Constructeur avec arguments
  @param dn noeud XML 
  @param rank_ rang du processus 
  @param nbranks_ nombre total de processus */
  PostProcessingWriter(DOMNode* dn,int const& rank_,int const& nbranks_):
  	m_rank(rank_), m_nprocs(nbranks_) {};

  /** @brief Destructeur */
  virtual ~PostProcessingWriter() {};
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
  virtual void PostProcessing_start(Scalar const& temps,
  	Scalar const& dt,
  	list<Particule*> const* particules,
	list<Particule*> const* pwait,
	list<Particule*> const* pperiodiques,	
	vector<Particule*> const* ParticuleClassesReference,
	Obstacle* obstacle,
	LinkedCell const* LC,
	vector<Fenetre> const& insert_windows ) = 0;

  /** @brief Ecriture d'evolution 
  @param temps Temps de sauvegarde  
  @param dt pas de temps  
  @param particules particules actives 
  @param pwait particules en attente
  @param pperiodiques particules periodiques  
  @param ParticuleClassesReference classes de reference de particule
  @param obstacle obstacles 
  @param LC grille de cellules  */
  virtual void PostProcessing( Scalar const& temps,
  	Scalar const& dt,
  	list<Particule*> const* particules,
	list<Particule*> const* pwait,
	list<Particule*> const* pperiodiques,		
	vector<Particule*> const* ParticuleClassesReference,
	Obstacle *obstacle,
	LinkedCell const* LC )=0;

  /** @brief Clot les ecritures */
  virtual void PostProcessing_end()=0;
  
  /** @brief Numero de cycle initial 
  @param cycle0 numero de cycle initial */
  virtual void setInitialCycleNumber( const int& cycle0 ) {};
  
  /** @brief Type de post-processing writer */
  virtual string getPostProcessingWriterType() const=0;
  
  /** @brief Le post-processing writer est il parallele ? */
  virtual bool isParallelWriter() const {return false;};
  
  /** @brief Le post-processing writer est il CompFeatures ? */
  virtual bool isCompFeaturesWriter() const { return false; };
  
  /** @brief Ecriture des composants pour une erreur de contact ou de 
  deplacement
  @param filename racine du nom du fichier
  @param errcomposants liste contenant les 2 composants */
  virtual void writeErreurComposantsPostProcessing( string const& filename,
  	list<Composant*> const& errcomposants ) {};  
  //@}
  
  /** @name Static Methods */
  //@{
  static void allocate_PostProcessingWindow( const int& nbRank );
  
  static void set_PostProcessingWindow( const int& rank,
  	const bool& bPPWindow );

  static vector<bool> get_PostProcessingWindow();
  //@}
  

protected:
  /**@name Parameters */
  //@{  
  int m_rank; /**< rang du processus */
  int m_nprocs; /**< nombre total de processus */
  //@}

  /**@name Static Parameters */
  //@{  
  static vector<bool> m_bPPWindow;
  //@}

  /** @name Constructors */
  //@{
  /** @brief Constructeur par defaut */
  PostProcessingWriter();
  //@}
  
};

#endif
  
