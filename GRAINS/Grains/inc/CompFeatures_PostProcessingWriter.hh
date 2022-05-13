#ifndef _CompFeatures_PostProcessingWriter__
#define _CompFeatures_PostProcessingWriter__

#include "PostProcessingWriter.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
using std::ofstream;
using std::ostringstream;
using std::vector;


/** @brief Classe CompFeatures_PostProcessingWriter.

    Ecriture des resultats pour post-processing par
    Paraview de GrainsCompFeatures

    @author D. RAKOTONIRINA - Institut Francais du Petrole - 2014 - Creation */
//=============================================================================
class CompFeatures_PostProcessingWriter : public PostProcessingWriter
{
public:
  /** @name Constructors */
  //@{
  /** @brief Constructeur par defaut 
  @param dn noeud XML 
  @param rank_ rang du processus 
  @param nbranks_ nombre total de processus */
  CompFeatures_PostProcessingWriter( DOMNode* dn,int const& rank_,
  	int const& nbranks_ );
  
  /** @brief Constructeur avec arguments 
  @param rank_ rang du processus 
  @param nbranks_ nombre total de processus   
  @param name_ nom des fichiers
  @param root_ racine du nom des fichiers
  @param isBinary ecriture en mode binaire */
  CompFeatures_PostProcessingWriter( int const& rank_,
  	int const& nbranks_,
  	const string &name_,
	const string &root_,
  	const bool &isBinary );  

  /** @brief Destructeur */
  virtual ~CompFeatures_PostProcessingWriter();
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
  string getPostProcessingWriterType() const { return "Paraview"; } 
	
  /** @brief Numero de cycle initial 
  @param cycle0 numero de cycle initial */
  void setInitialCycleNumber( const int& cycle0 ) 
  	{ m_ParaviewCycleNumber = cycle0; };

  /** @brief Le post-processing writer est il parallele ? */
  bool isCompFeaturesWriter() const { return true; }

  /** @brief Ecriture des composants pour une erreur de contact ou de 
  deplacement
  @param filename racine du nom du fichier
  @param errcomposants liste contenant les 2 composants */
  void writeErreurComposantsPostProcessing( string const& filename,
  	list<Composant*> const& errcomposants );
  //@}
  
  
  /** @name Methods pour ecrire en binaire */
  //@{
  void start_output_binary( int size, int number ) ;
  void write_double_binary( double val ) ;
  void write_int_binary( int val ) ;
  int store_int_binary( int val ) ;
  void flush_binary( std::ofstream& file ) ;
  void check_allocated_binary( int size ) ;
  void compress_segment_binary( int seg ) ;
  //@}  


private:
  string m_ParaviewFilename_root; /**< racine du nom de fichier de sortie 
  	Paraview */
  string m_ParaviewFilename; /**< nom de fichier de sortie Paraview */  
  ostringstream m_Paraview_saveObstacles_pvd; /**< flux pour le fichier pvd des
  	obstacles */
  ostringstream m_Paraview_saveObstaclesPeriodiques_pvd; /**< flux pour le 
  	fichier pvd des obstacles periodiques */	
  vector<ostringstream*> m_Paraview_saveParticules_pvd; /**< flux pour le 
  	fichier pvd des particules */
  ostringstream m_Paraview_saveClonesPeriodiques_pvd; /**< flux pour le fichier
  	pvd des clones periodiques */	
  ostringstream m_Paraview_saveVectors_pvd; /**< flux pour le fichier pvd des
  	vecteurs vitesse de translation des particules */
  ostringstream m_Paraview_saveContactForces_pvd; /**< flux pour le fichier pvd 
  	des vecteurs force de contact entre composants */		
  int m_ParaviewCycleNumber; /**< numero initial des fichiers de sortie */
  bool m_binary; /**< ecriture des fichiers en binaire */
  bool m_postProcessObstacle; /**< post-processing des obstacles */
  char * BUFFER ;
  int ALLOCATED ;
  int OFFSET ;
  int CURRENT_LENGTH ;
  list<string> empty_string_list; 

  
  /** @name Methods */
  //@{  
  /** @brief Ecriture
  @param temps Temps de sauvegarde  
  @param dt pas de temps
  @param particules particules actives
  @param pperiodiques particules periodiques  
  @param ParticuleClassesReference classes de reference de particule
  @param obstacle obstacles 
  @param LC grille de cellules */
  void one_output( Scalar const& temps,
  	Scalar const& dt,
	list<Particule*> const* particules,
	list<Particule*> const* pperiodiques,
	vector<Particule*> const* ParticuleClassesReference,
  	Obstacle *obstacle,
	LinkedCell const* LC ); 

  /** @brief Ecriture
  @param pwait Particules en attente */
  void one_output_features( list<Particule*> const* pwait );

  /** @brief Ecriture des obstacles
  @param allObstacles des obstacles primaires 
  @param obsFilename nom du fichier */
  void writeObstaclesPostProcessing_Paraview(
  	list<MonObstacle*> const &allObstacles,
  	string const &obsFilename );
  
  /** @brief Ecriture des particules
  @param particules particules actives 
  @param partFilename nom du fichier 
  @param forceForAllTag ecrit la particule quel que soit son tag */
  void writeParticulesPostProcessing_Paraview(
  	list<Particule*> const* particules, const string &partFilename,
	bool const& forceForAllTag = false );
	
  /** @brief Ecriture des particules de forme sphérique sous forme d'un vecteur 
  @param particules particules actives 
  @param partFilename nom du fichier 
  @param forceForAllTag ecrit la particule quel que soit son tag */
  void writeSpheresPostProcessing_Paraview(
  	list<Particule*> const* particules, const string &partFilename,
	bool const& forceForAllTag = false );	
	
  /** @brief Ecriture des vecteurs vitesse de translation & rotation des 
  particules
  @param particules particules actives 
  @param partFilename nom du fichier */
  void writeVectorsMotionPostProcessing_Paraview(
  	list<Particule*> const* particules, const string &partFilename );

  /** @brief Ecriture des vecteurs force de contact entre composants
  @param particules particules actives   
  @param LC grille de cellules 
  @param partFilename nom du fichier */
  void writeVectorsForcePostProcessing_Paraview(
  	list<Particule*> const* particules,
  	LinkedCell const* LC, const string &partFilename );
				
  /** @brief Ecriture du fichier pvtu 
  @param filename nom du fichier 
  @param pointVector liste de noms du champs vectoriels sur les points 
  @param pointScalar liste de noms du champs scalaires sur les points   
  @param cellScalar liste de noms du champs scalaires sur les cellules */
  void writePVTU_Paraview( const string &filename,
  	list<string> const* pointVector,
	list<string> const* pointScalar,
	list<string> const* cellScalar );
	
  /** @brief Ecriture du maillage de cellules du LinkedCell
  @param LC grille de cellules 
  @param partFilename nom du fichier */
  void writeLinkedCellPostProcessing_Paraview(
  	LinkedCell const* LC, const string &partFilename );

  /** @brief Ecriture des fenetres d'insertion
  @param insert_windows fenetres d'insertion 
  @param partFilename nom du fichier */
  void writeInsertionPostProcessing_Paraview(
   	vector<Fenetre> const& insert_windows,
  	const string &partFilename );
		
  /** @brief Mise a jour de l'indicateur des obstacles 
  @param temps Temps de sauvegarde  
  @param dt pas de temps
  @param obstacle obstacles */
  void updateObstaclesIndicator( Scalar const& temps,
  	Scalar const& dt, Obstacle *obstacle );
  
  /** @brief Efface les fichiers resultats */
  void clearResultFiles() const;
  
  /** @brief Relit un fichier pvd dans le cas d'un restart dans la continuite 
  et transfère son contenu dans le flux correspondant 
  @param filename nom du fichier d'entree
  @param ossflux flux correspondant */
  void readPVDFile( string const& filename, ostringstream& ossflux );
  
  /** @brief Recupere le dernier numero de cycle dans le cas d'un restart 
  dans la continuite */
  int getPreviousCycleNumber() const;
  //@}
      

protected:
  /** @name Constructors */
  //@{
  /** @brief Constructeur par defaut */
  CompFeatures_PostProcessingWriter() {};
  //@}  
};

#endif
  
