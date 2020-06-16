#ifndef _GrainsMPI
#define _GrainsMPI

#include <mpi.h>
#include "MPIWrapperGrains.hh"
#include "Grains.H"

#include <list>
#include <vector>
#include <string>
#include <iostream>
using namespace std;

#include "ReaderXML.hh"

class Voisins;

/** @brief Classe principale de gestion du probleme GrainsMPI.

    Construction du probleme en fonction du fichier de description.
    Gestion de la simulation par décomposition de domaines (MPI).

@author L.GIROLAMI, A.WACHS - Institut Francais du Petrole - 2009 - Creation */
//=============================================================================
class GrainsMPI : virtual public Grains
{
public:
  /**@name Constructors */
  //@{
  /** @brief Constructeur par defaut */
  GrainsMPI();

  /** @brief Destructeur */
  virtual ~GrainsMPI();
  //@}


  /**@name Methods Virtual */
  //@{
  /** @brief Construction du probleme
  @param rootElement Le noeud racine */
  virtual void Construction( DOMElement* rootElement );

  /** @brief Appel a la simulation granulaire.
  @param predict Simulation en etat de prediction
  @param isPredictorCorrector le schéma de couplage avec le fluide est il de
  type prédicteur-correcteur
  @param explicit_added_mass la masse ajoutée est elle traitée de manière
  explicite */
  virtual void Simulation( bool predict = true,
  	bool isPredictorCorrector = false,
	bool explicit_added_mass = false );
  //@}


protected:
  /** @name Parameters */
  //@{
  MPIWrapperGrains *m_wrapper; /**< pour la gestion des communications MPI */
  string m_MPIstrategie; /**< strategie de communication MPI pour échanger entre
  	les sous-domaines */  
  //@}


  /**@name Methods */
  //@{
  /** @brief Insertion d'une particule dans les algorithmes 
  @param mode mode d'insertion */
  virtual bool insertParticule( const PullMode& mode );
  
  /** @brief Positionne les particules en attente pour la simulation 
  @param mode mode d'insertion */
  virtual void setPositionParticulesFichier( const PullMode& mode=ORDER ); 

  /** @brief Positionne les particules en attente pour la simulation sous forme
  d'un bloc structure de positions
  @param mode mode d'insertion */
  virtual void setPositionParticulesBloc( const PullMode& mode = ORDER );  
  
  /** @brief Lecture des donnees pour les simulations MPI
  @param lx taille du domaine global en X 
  @param ly taille du domaine global en Y   
  @param lz taille du domaine global en Z  
  @param root le noeud XML */
  virtual void readDomainDecomposition(DOMNode* root,
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
      
  /** @brief Numero maximum de particules */
  virtual int numeroMaxParticules() const; 
  
  /** @brief Cree, insert et link les nouvelles particules */
  virtual void InsertCreateNewParticules();
  
  /** @brief Tire une classe de particules parmi les nouvelles particules a
  inserer 
  @param mode mode d'insertion 
  @param ParticuleClassesForCreation classes de particules a creer 
  @param random_local random sur tous les proc ou sur le proc */
  Particule* getParticuleClasseForCreation( const PullMode& mode,
  	list< pair<Particule*,int> >& ParticuleClassesForCreation,
	bool const& random_local );
	
  /** @brief Traitement des particules periodiques 
  @param b_perfTiming timing des performances */
  void periodicParticles( bool b_perfTiming );
  
  /** @brief Ecriture de la memoire utilisee par la simulation */
  void display_used_memory() const;  	

  /** Synchronize the PPWindow boolean relative to each sub-domain.
  vector of bool is then correct on every processor*/
  void synchronize_PPWindow();
  //@}

};

#endif
  
