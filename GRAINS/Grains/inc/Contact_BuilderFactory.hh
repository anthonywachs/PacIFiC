#ifndef _Contact_BuilderFactory
#define _Contact_BuilderFactory

#include "ReaderXML.hh"
#include <iostream>
#include <map>
#include <string>
using namespace std;

class ContactLaw;
class Particule;
class Obstacle;


/** @brief Fabrique des Contacts.

Construit les tables de type de contact + parametres associees aux materiaux.
Pour un couple de composants selectionne la loi de contact adequate.

@author Grains Project - IFP - 2007 - Creation */
// ============================================================================
class Contact_BuilderFactory
{
public:
  /** @name Methods Static */
  //@{
  /** @brief Renvoie la loi de contact associee aux materiaux
  @param matA Materiau A
  @param matB Materiau B */
  static ContactLaw* contactForceModel( const string &matA,
	const string &matB );

  /** @brief Construction de la table des contacts
  @param root Le noeud "<Contact>" */
  static void define( DOMNode* root );

  /** @brief Definition des materiaux presents dans la simulation
  @param materiau Nom du materiau 
  @param matForObstacle le materiau s'applique a une obstacle */
  static void defineMaterial( const string &materiau,
  	const bool& matForObstacle );
  
  /** @brief Lecture du fichier de reload des contacts definis
  @param file Flux du fichier principal de reload */
  static void reload( istream &file );
  
  /** @brief Verifie si un materiau ne s'applique qu'a un obstacle en cas de
  reload car cette info n'est pas conservee dans le fichier de reload
  @param ParticuleClassesReference particules de référence */
  static void set_materialsForObstaclesOnly_reload( 
  	vector<Particule*> const* ParticuleClassesReference );  

  /** @brief Ecriture du fichier de reload des contacts definis
  @param file Flux du fichier principal de reload  
  @param contactFile Racine du fichier de reload des contacts 
  @param rank rang du processus */
  static void save( ostream &file, const string& contactFile, const int& rank );

  /** @brief Renvoie le nom du matériau 
  @param num numéro entier dans la map(string,int) */
  static string nomMaterial( int num );

  /** @brief Renvoie le numéro du matériau 
  @param mat nom du matériau dans la map(string,int) */
  static int numeroMaterial( const string &mat );
  
  /** @brief Verifie que toutes les lois de contact existent pour l'ensemble
  des materiaux definis 
  @param matA materiau 1 
  @param matB materiau 2 */
  static bool checkContactLawsExist( string &matA, string &matB ); 
  
  /** @brief Destruit les lois de contact */
  static void eraseAllContactLaws();    
  //@}


private:
  /** @name Enumeration */
  //@{
  enum ContactType {
    /// Contact elastique
    MYCONTACT,
    /// Contact elastique basee sur un coef de restitution
    ERCONTACT,
    /// Contact cohesif
    COHCONTACT,
    /// Contact elastique basee sur un coef de restitution + histoire de
    /// deplacement tangentiel
    ERHCONTACT    
  };
  //@}


  /** @name Structure */
  //@{
  struct ContactFeatures {
    ContactType         name;
    map<string,double> values;
  };
  //@}


  /** @name Methods */
  //@{
  /** @brief Definition des parametres associes au Contact
  @param root Le noeud "<Contact>" 
  @return Les donnees de Contact */
  static pair<Contact_BuilderFactory::ContactFeatures,ContactLaw*> 
  	defineParameters(DOMNode *root) ;
  //@}


  /**@name Constructors */
  //@{
  /** Classe Statique : pas de contructeur disponible. */
  Contact_BuilderFactory() {};

  /** Classe Statique : pas de destructeur disponible. */
  ~Contact_BuilderFactory() {};
  //@}


  /** @name Parameters */
  //@{  
  static map<string,int> m_materials; /**< Valeurs des materiaux dans la 
  	table */
  static map<string,bool> m_materialsForObstaclesOnly; /**< Vrai si le materiau
  	ne s'applique qu'a des obstacles */	  
  static map<int,Contact_BuilderFactory::ContactFeatures> 
  	m_contactParametres; /**< Table des parametres de contact */ 
  static map<int,ContactLaw*> m_contactLaws; /**< Table des lois de contact */
  static int m_value; /**< Valeur pour definir les materiaux */
  //@}
};

#endif
