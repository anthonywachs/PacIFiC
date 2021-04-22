#ifndef _GrainsGeomTest
#define _GrainsGeomTest

#include "Grains.H"
#include "EnsComposant.H"
#include "App.H"
#include "LinkedCell.H"
#include "ReaderXML.hh"
#include <list>
#include <vector>
#include <string>
#include <iostream>
using namespace std;


/** @brief Classe pour tester les intersections de convexes

@author A.WACHS - Institut Francais du Petrole - 2010 - Creation */
//=============================================================================
class GrainsGeomTest : public Grains
{
public:
  /**@name Constructors */
  //@{
  /** @brief Constructeur par defaut */
  GrainsGeomTest();

  /** @brief Destructeur */
  virtual ~GrainsGeomTest();
  //@}


  /**@name Methods Virtual */
  //@{
  /** @brief Construction de la simulation
  @param rootElement Le noeud racine */
  virtual void Chargement(DOMElement* rootElement);

  /** @brief Construction du probleme
  @param rootElement Le noeud racine */
  virtual void Construction(DOMElement* rootElement);

  /** @brief Construction des forces actives
  @param rootElement Le noeud racine */
  virtual void Forces(DOMElement* rootElement);

  /** @brief Appel a la simulation granulaire
  @param predict Simulation en etat de prediction
  @param isPredictorCorrector le schéma de couplage avec le fluide est il de
  type prédicteur-correcteur
  @param explicit_added_mass la masse ajoutée est elle traitée de manière
  explicite */
  virtual void Simulation(bool predict = true,
  	bool isPredictorCorrector = false,
	bool explicit_added_mass = false);
  //@}
  
protected:
  /**@name Methods */
  //@{
  /** @brief Renvoi une transformation de type rotation en fonction des angles
  @param angleX rotation par rapport à l'axe X 
  @param angleY rotation par rapport à l'axe Y   
  @param angleZ rotation par rapport à l'axe Z */
  Matrix getMatrixRotation(const Scalar &angleX,
  	const Scalar &angleY,
	const Scalar &angleZ) const;        
  //@}


};

#endif
  
