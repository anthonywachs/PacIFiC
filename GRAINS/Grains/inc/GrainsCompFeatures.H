#ifndef _GrainsCompFeatures
#define _GrainsCompFeatures

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

/** @brief Classe principale de gestion du probleme GrainsCompFeatures.

    Construction du probleme en fonction du fichier de description.
    Gestion des elements d'inertie dans le cas de particules composites.

    @author D.RAKOTONIRINA - IFP Energies Nouvelles - 2014 - Creation */
//=============================================================================


class GrainsCompFeatures : virtual public Grains
{
public:
  /** @name Constructors */
  //@{
  /** Constructeur par defaut */
  GrainsCompFeatures();

  /** Destructeur */
  virtual ~GrainsCompFeatures();
  //@}



  /** @name Methods Virtual */
  //@{
  /** @brief Construction de la simulation
  @param rootElement Le noeud racine */
  void Chargement( DOMElement* rootElement );

  /** @brief Construction du probleme
  @param rootElement Le noeud racine */
  void Construction( DOMElement* rootElement );

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
  string m_ParaviewFilename; /**< nom de fichier de sortie Paraview */
  //@}


  /**@name Methods */
  //@{
  /** @brief Determine si il y a intersection entre deux convexes 
  @return Vrai si les deux convexes sont en contact
  @param a Une Forme A
  @param b Une Forme B */
  bool isIn( const Forme &a, const Forme &b );

  /** @brief Determine s'il y a intersection entre une Forme et
  un Composant
  @return Vrai si Forme et Composant s'intersectent
  @param composant Un Composant
  @param b Une Forme B */
  bool isIn( Composant* composant, const Forme &b ) ;

  //@}
};

#endif
  
