// G.FERRER - Aout.2000 - Creation
// ============================================================================
#ifndef _Erreur__
#define _Erreur__

#include <string>
#include <iostream>
#include <list>
using namespace std;

class Composant;


/** @brief Gestion des erreurs de contact entre formes.

    @author G.FERRER - Institut Francais du Petrole - 2000 - Creation */
// ============================================================================
class ErreurContact
{
public:
  /**@name Constructeurs */
  //@{
  /** @brief Contructeur par defaut. */
  ErreurContact();
  
  /** @brief Destructeur */
  ~ErreurContact();  
  //@}

  /**@name Methods */
  //@{
  /** @brief Message de l'exception : identificateur des composants en contact 
  @param fileOut Flux de sortie */
  void Message( ostream &fileOut ) const;

  /** @brief Affectation des identificateurs des composants en contact
  @param id0_ Premier composant
  @param id1_ Deuxieme composant 
  @param temps_ Temps de simulation */
  void setComposants( Composant* id0_, Composant* id1_,double temps_ );

  /** @brief Affectation du message decrivant l'exception
  @param mes Message de l'exception */
  void setMessage( const string &mes );
  
  /** @brief Renvoie les 2 composants dans une liste pour post-processing */
  list<Composant*> getComposants();
  //@}

private:
  /**@name Parameters */
  //@{
  string message; /**< Message de definition du probleme de contact */
  Composant *id0; /**< Composant en contact 0 */
  Composant *id1; /**< Composant en contact 1 */
  double m_temps; /**< temps de simulation */  
  //@}
};




/** @brief Erreur sur deplacement
    @author GRAINS Project - IFP - 2007 */
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class ErreurDeplacement
{
public:
  /** @name Constructors */
  //@{
  /** @brief Constructeur
  @param id0_ composant
  @param depl Deplacement calcule 
  @param deplMax Deplacement maxi autorise 
  @param temps_ Temps de simulation */
  ErreurDeplacement( Composant* id0_, double depl, double deplMax, 
  	double temps_ );
  
  /** @brief Destructeur */
  ~ErreurDeplacement();
  //@}
  
  /**@name Methods */
  //@{
  /** @brief Message de l'exception : identificateur des composants en contact 
  @param fileOut Flux de sortie */
  void Message( ostream &fileOut ) const;
  
  /** @brief Renvoie le 2 composant dans une liste pour post-processing */
  list<Composant*> getComposant();  
  //@}

private:
  /** @name Parameters */
  //@{
  Composant *id0; /**< Composant */  
  double m_depl; /**< Deplacement calcule */
  double m_deplMax; /**< Deplacement max autorise */
  double m_temps; /**< temps de simulation */
  //@}

  /** @name Contructor */
  //@{
  /** @brief Constructeur par defaut interdit */
  ErreurDeplacement() {};
  //@}
};




/** @brief Gestion des erreurs de simulation.

    @author G.FERRER - Institut Francais du Petrole - 2004 - Creation */
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class ErreurSimulation
{
public:
  /**@name Constructors */
  //@{
  /** @brief Constructeur par defaut */
  ErreurSimulation();

  /** @brief Constructeur avec indication de la methode de levee d'erreur 
  @param str type de methode de levee d'erreur */
  ErreurSimulation( const string &str );

  /** @brief Destructeur */
  ~ErreurSimulation();
  //@}

  /**@name Methods */
  //@{
  /** @brief Message d'erreur pour arret simulation 
  @param fileOut flux de sortie */
  void Message( ostream &fileOut ) const;
  //@}


private:
  /**@name Parameters */
  //@{  
  string methode; /**< Nom de la methode de levee d'exception */
  //@}
};

#endif
