#ifndef _GrainsParameters
#define _GrainsParameters

#include "Grains.H"
#include "EnsComposant.H"
#include "App.H"
#include "LinkedCell.H"

#include <list>
#include <vector>
#include <string>
#include <iostream>
using namespace std;

#include "ReaderXML.hh"


/** @brief Classe pour évaluer les paramètres d'une simulation 

@author A.WACHS - Institut Francais du Petrole - 2010 - Creation */
//=============================================================================
class GrainsParameters : public Grains
{
public:
  /**@name Constructors */
  //@{
  /** @brief Constructeur par defaut */
  GrainsParameters();

  /** @brief Destructeur */
  virtual ~GrainsParameters();
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
  /** @name Parameters */
  //@{  
  double m_vRelative; /**< Estimation de la vitesse relative max lors 
  	d'un contact */
  //@}  

};

#endif
  
