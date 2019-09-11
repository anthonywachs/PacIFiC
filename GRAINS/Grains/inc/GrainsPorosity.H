#ifndef _GrainsPorosity
#define _GrainsPorosity

#include <mpi.h>
#include "MPIWrapperGrains.hh"
#include "Grains.H"
#include "Grains_Exec.hh"
#include "EnsComposant.H"
#include "App.H"
#include "solvercomputingtime.hh"
#include "computingtime.hh"
#include "LinkedCell.H"

#include <list>
#include <vector>
#include <string>
#include <iostream>
using namespace std;

#include "ReaderXML.hh"

/** @brief Classe principale de gestion du probleme GrainsPorosity.

    Construction du probleme en fonction du fichier de description.
    Calcul de la porosite.

    @author D.RAKOTONIRINA - IFP Energies Nouvelles - 2014 - Creation */
//=============================================================================


class GrainsPorosity : virtual public Grains
{
public:
  /** @name Constructors */
  //@{
  /** Constructeur par defaut */
  GrainsPorosity();

  /** Destructeur */
  virtual ~GrainsPorosity();
  //@}



  /** @name Methods Virtual */
  //@{
  /** @brief Construction de la simulation
  @param rootElement Le noeud racine */
  void Chargement( DOMElement* rootElement );

  /** @brief Construction des forces actives
  @param rootElement Le noeud racine */
  void Forces( DOMElement* rootElement );

  /** @brief Appel a la simulation granulaire
  @param predict Simulation en etat de prediction
  @param isPredictorCorrector le schéma de couplage avec le fluide est il de
  type prédicteur-correcteur
  @param explicit_added_mass la masse ajoutée est elle traitée de manière
  explicite */
  virtual void Simulation( bool predict = true,
  	bool isPredictorCorrector = false,
	bool explicit_added_mass = false );
  
  /** @brief Lecture des donnees pour les simulations MPI
  @param lx taille du domaine global en X 
  @param ly taille du domaine global en Y   
  @param lz taille du domaine global en Z  
  @param root le noeud XML */
  void readDomainDecomposition( DOMNode* root, const Scalar& lx,
	const Scalar& ly, const Scalar& lz ); 
  //@}
  
   
protected:
  /** @name Parameters */
  //@{
  
  //@}


  /**@name Methods */
  //@{
      
  //@}
  

private:
  /** @name Parameters */
  //@{
  int m_NX; /**< Decompositon du domaine en X */
  int m_NY; /**< Decompositon du domaine en Y */
  int m_NZ; /**< Decompositon du domaine en Z */
  int m_nptsX; /**< Nombre de points de discretisation suivant x */ 
  int m_nptsY; /**< Nombre de points de discretisation suivant y */ 
  int m_nptsZ; /**< Nombre de points de discretisation suivant z */ 
  int m_lx; /**< Nombre de couches dans la direction x */
  int m_ly; /**< Nombre de couches dans la direction y */
  int m_lz; /**< Nombre de couches dans la direction z */
  string m_ParaviewFilename; /**< nom de fichier de sortie Paraview */
  string m_Structure; /**< Type de structure ( porosite radiale ou
    cartesienne ) */
  ofstream fOut; /**< Output file */
  string m_output_root; /**< Output root */
  string m_output_name; /**< Output file name*/
  double m_radiusIn; /**< Inner radius*/
  double m_radiusOut; /**< Outer radius*/
  Vecteur m_center; /**< Center of the cylinder*/
  int m_lr; /**< Nombre de couches dans la direction radiale */
  MPIWrapperGrains *m_wrapper; /**< pour la gestion des communications MPI */
  string m_MPIstrategie; /**< strategie de communication MPI pour échanger entre
  	les sous-domaines */  
  //@}


  /**@name Methods */
  //@{
  /** @brief Determine si il y a intersection entre deux convexes 
  @return Vrai si les deux convexes sont en contact
  @param a Une Forme A
  @param b Une Forme B */
  bool isIn( const Forme &a, const Forme &b );

  /** @brief Determine si un point se trouve dans une couronne
  r_in  < x < r_out 
  @return Vrai si le point est dans la couronne
  @param x Coordonnes en X
  @param y Coordonnes en Y
  @param z Coordonnes en Z
  @param inRad Rayon interieur
  @param outRad Rayon exterieur
  @param b Une Forme B */
  bool isInRadius( const double &x, const double &y, const double &z,
    const double &inRad, const double &outRad );

  /** @brief Determine s'il y a intersection entre une Forme et
  un Composant
  @return Vrai si Forme et Composant s'intersectent
  @param composant Un Composant
  @param b Une Forme B */
  bool isIn( Composant* composant, const Forme &b ) ;
  
  /** @brief Script pour effacer un fichier */
  void clearOneFile( ofstream&clear_file, string const& filename ) const ;  	  
	
  /** @brief Efface les fichiers resultats */
  void clearResultFiles() const ;

  //@}
};

#endif
  
