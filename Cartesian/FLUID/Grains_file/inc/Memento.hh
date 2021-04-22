#ifndef _ParticuleMemento
#define _ParticuleMemento


#include "Quaternion.H"
#include "Transform.H"
#include "Vecteur.H"

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#include "Particule.H"


/** @brief Memento pour une configuration d'un composant.

Memorise la configuration d'un composant pour permettre de le restaurer
si necessaire. Cette classe n'est accessible que par l'intermediaire d'un
Composant. 

Adaptation "assez libre" du Desing Pattern Memento.
@author GRAINS Project - IFP - 2008 
*/
class ConfigurationMemento
{
public:
  /** @name Constructor & Destructor */
  //@{
  /** @brief Destructor */
  virtual ~ConfigurationMemento() {};
  //@}


private:
  friend class Particule;
  friend class Obstacle;
  friend class CompObstacle;  
  friend class MonObstacle;
  friend class ParticulePeriodique;
  friend class Composant;  
  
  /** @name Constructor */
  //@{
  /** @brief Constructor */
  ConfigurationMemento() {};
  //@}


  /** @name Parameters */
  //@{
  Transform m_position; /**< Transformation entre configuration courante et
  	configuration de base de la particule */
  //@}
};




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#include "CineParticule.H"
#include "LeapFrog_3D.hh"
#include "LeapFrog_2D.H"


class CineParticuleMemento
{
public:
  /** @name Constructor & Destructor */
  //@{
  /** @brief Destructor */
  virtual ~CineParticuleMemento() {};
  //@}
  

private:
  /** @name Classes friend */
  //@{
  friend class CineParticule;
  friend class LeapFrog_3D;
  friend class LeapFrog_2D;
  //@}
  
  /** @name Constructor */
  //@{
  /** @brief Constructor */
  CineParticuleMemento() {};
  //@}


  /** @brief Parameters */
  //@{
  Quaternion m_QuaternionRotation; /**< quaternion de rotation */
  solid::Vecteur m_vitesseT; /**< vitesse de translation */
  Quaternion m_QuaternionVitesseR; /**< quaternion de vitesse de rotation */
  //@}
};




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#include "CineObstacle.H"
#include "LeapFrog_3D.hh"
#include "LeapFrog_2D.H"


class CineObstacleMemento
{
public:
  /** @name Constructor & Destructor */
  //@{
  /** @brief Destructor */
  virtual ~CineObstacleMemento() {};
  //@}
  

private:
  /** @name Classes friend */
  //@{
  friend class CineObstacle;
  //@}
  
  /** @name Constructor */
  //@{
  /** @brief Constructor */
  CineObstacleMemento() {};
  //@}


  /** @brief Parameters */
  //@{
  solid::Vecteur m_vitesseTranslation; /**< vitesse de translation de 
  	l'obstacle */
  solid::Vecteur m_vitesseRotation; /**< vitesse de rotation de l'obstacle */ 
  //@}
};


#endif
