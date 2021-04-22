#ifndef _GrainsCoupledWithFluidMPI
#define _GrainsCoupledWithFluidMPI

#include "GrainsMPI.H"
#include "GrainsCoupledWithFluid.hh"

#include <list>
#include <string>
using namespace std;

#include "ReaderXML.hh"


/** @brief Interface entre la modelisation fluide et la modelisation sec,
    Version MPI
       
    @author A.WACHS - Institut Francais du Petrole - 2009 - Creation */
// ============================================================================
class GrainsCoupledWithFluidMPI : public GrainsCoupledWithFluid, 
	public GrainsMPI
{
public:
  /** @name Constuctors & Destructor */
  //@{
  /** @brief Constructeur
  @param rhoFluide Masse volumique du fluide dans lequel baignent les 
  particules 
  @param grid_size size of the smallest grid cell */
  GrainsCoupledWithFluidMPI( Scalar rhoFluide, double grid_size );

  /** @brief Destructeur */
  virtual ~GrainsCoupledWithFluidMPI();
  //@}


  /** @name Methods Virtual */
  //@{
  /** @brief Construction de la simulation
  @param rootElement Le noeud racine */
  virtual void Chargement( DOMElement* rootElement );

  /** @brief Construction du probleme
  @param rootElement Le noeud racine */
  virtual void Construction( DOMElement* rootElement );

  /** @brief Construction des forces actives
  @param rootElement Le noeud racine */
  virtual void Forces( DOMElement* rootElement );
  
  /** @brief Sauvegarde de l'etat de simulation */
  virtual void Save( const string &ext ) const;

  /** @brief Appel a la simulation granulaire
  @param predict Simulation en etat de prediction
  @param isPredictorCorrector le schéma de couplage avec le fluide est il de
  type prédicteur-correcteur
  @param explicit_added_mass la masse ajoutée est elle traitée de manière
  explicite */
  virtual void Simulation( bool predict = true,
  	bool isPredictorCorrector = false,
	bool explicit_added_mass = false );
	
  /** @brief Mise à jour de la vitesse des particules. 
  @param b_set_velocity_nm1_and_diff mise ou non a jour des vitesses et des
  differences de vitesse du pas de temps flduide precedent */
  virtual void UpdateParticulesVelocities(
  	const bool &b_set_velocity_nm1_and_diff );
  
  /** @brief Ecriture des particules dans un fichier 
  @param filename nom du fichier */
  virtual void WriteParticulesInFluid( const string &filename ) const; 
  
  /** @brief Mise à jour de la vitesse des particules. 
  @param b_set_velocity_nm1_and_diff mise ou non a jour des vitesses et des
  differences de vitesse du pas de temps flduide precedent   
  @param velocities tableau contenant les vitesse de translation & rotation 
  @param b_set_velocity_nm1_and_diff mise a jour de la vitesse au pas de temps
  precedent ainsi que la difference de vitesse explicite */
  virtual void UpdateParticulesVelocities( 
  	const vector<vector<double> > &velocities,
  	const bool &b_set_velocity_nm1_and_diff );
 
  /** @brief Ecriture des particules dans un flux 
  @param is flux d'entrée */
  virtual void WriteParticulesInFluid( istringstream &is ) const;
  
  /** @brief Ecriture de la vitesse et du centre de gravite des particules dans
  un fichier 
  @param filename nom du fichier */
  virtual void WritePVGCInFluid( const string &filename ) const;        

  /** @brief Ecriture de la vitesse et du centre de gravite des particules dans
  un flux 
  @param is flux d'entrée */
  virtual void WritePVGCInFluid( istringstream &is ) const;
  
  /** @brief Sauvegarde par defaut de l'etat initial pour post-processing */
  virtual void InitialPostProcessing();
  
  /** @brief Sauvegarde pour post-processing et restart */
  virtual void doPostProcessing();
  
  /** @brief Operations a effectuer avant appel du destructeur */
  virtual void BeforeDestructor();        	
  //@}


protected:
  /**@name Methods */
  //@{
  /** @brief Insertion d'une particule dans les algorithmes 
  @param mode mode d'insertion */
  virtual bool insertParticule( const PullMode& mode );
  
  /** @brief Positionne les particules en attente pour la simulation a partir
  d'un fichier de positions 
  @param mode mode d'insertion */
  virtual void setPositionParticulesFichier( const PullMode& mode = ORDER );
  
  /** @brief Positionne les particules en attente pour la simulation sous forme
  d'un bloc structure de positions
  @param mode mode d'insertion */
  virtual void setPositionParticulesBloc( const PullMode& mode = ORDER );
  
  /** @brief Lecture des donnees pour les simulations MPI
  @param lx taille du domaine global en X 
  @param ly taille du domaine global en Y   
  @param lz taille du domaine global en Z  
  @param root le noeud XML */
  virtual void readDomainDecomposition( DOMNode* root,
  	const Scalar& lx, const Scalar& ly, const Scalar& lz ); 
  
  /** @brief Renvoie le nom complet du fichier de resultat .result
  @param rootname nom de la racine */
  virtual string fullResultFileName( const string &rootname ) const; 

  /** @brief Lecture des donnees pour les simulations periodiques
  @param rootElement le noeud XML */
  virtual void readPeriodic( DOMElement* rootElement );
  
  /** @brief Definition du LinkedCell 
  @param rayon rayon max des particules */
  virtual void defineLinkedCell( Scalar const& rayon );   
  
  /** @brief Arret du code en cas de probleme */
  virtual void grainsAbort() const;   
  //@}

};

#endif


