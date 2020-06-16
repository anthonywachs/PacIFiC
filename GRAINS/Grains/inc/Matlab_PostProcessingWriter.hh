#ifndef _Matlab_PostProcessingWriter__
#define _Matlab_PostProcessingWriter__

#include "PostProcessingWriter.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
using std::ofstream;
using std::ostringstream;
using std::vector;


/** @brief Classe Matlab_PostProcessingWriter.

    Ecriture des resultats pour post-processing par Matlab.

    @author A.HAMMOUTI - IFPEn - 2013 - Creation */
//=============================================================================
class Matlab_PostProcessingWriter : public PostProcessingWriter
{
public:
  /** @name Constructors */
  //@{
  /** @brief Constructeur par defaut 
  @param dn noeud XML 
  @param rank_ rang du processus 
  @param nbranks_ nombre total de processus */
  Matlab_PostProcessingWriter( DOMNode* dn,int const& rank_,
  	int const& nbranks_ );
  
  /** @brief Constructeur avec arguments 
  @param rank_ rang du processus 
  @param nbranks_ nombre total de processus   
  @param name_ nom des fichiers
  @param root_ racine du nom des fichiers
  @param isBinary ecriture en mode binaire */
  Matlab_PostProcessingWriter( int const& rank_,
  	int const& nbranks_,
  	const string &name_,
	const string &root_,
  	const bool &isBinary );  

  /** @brief Destructeur */
  virtual ~Matlab_PostProcessingWriter();
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
  void PostProcessing( Scalar const& temps,
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
  string getPostProcessingWriterType() const { return "Matlab"; } 
	
  /** @brief Le post-processing writer est il parallele ? */
  bool isParallelWriter() const { return true; }

  //@}
  
  
  /** @name Methods pour ecrire en binaire */
  //@{
  void start_output_binary( int size, int number ) ;
  void write_double_binary( double val ) ;
  int store_int_binary( int val ) ;
  void check_allocated_binary( int size ) ;
  //@}  


private:
  string m_MatlabFilename_root; /**< racine du nom de fichier de sortie 
  	Matlab */
  string m_MatlabFilename; /**< nom de fichier de sortie Paraview */  
  ofstream m_Matlab_Append; /**< Flux de donnees particules */
  int m_MatlabReloadNumber; /**< numero initial des fichiers de sortie */
  bool m_binary; /**< ecriture des fichiers en binaire */
  char * BUFFER ;
  int ALLOCATED ;
  int OFFSET ;
  int CURRENT_LENGTH ;
  
  
  /** @name Methods */
  //@{  
  /** @brief Ecriture
  @param temps Temps de sauvegarde  
  @param dt pas de temps
  @param particules particules actives*/
  void one_output( Scalar const& temps,
  	Scalar const& dt,
	list<Particule*> const* particules ); 
	
  
  /** @brief Efface les fichiers resultats */
  void clearResultFiles() const;
  
  /** @brief Relit un fichier binaire dans le cas d'un restart dans la continuite 
  et transfère son contenu dans le flux correspondant 
  @param filename nom du fichier d'entree
  @param ossflux flux correspondant */
  void readbinFile( string const& filename, ostringstream& ossflux );
  
  /** @brief Recupere le dernier numero de cycle dans le cas d'un restart 
  dans la continuite */
  int getPreviousCycleNumber() const;
  //@}
      

protected:
  /** @name Constructors */
  //@{
  /** @brief Constructeur par defaut */
  Matlab_PostProcessingWriter() {};
  //@}  
};

#endif
  
