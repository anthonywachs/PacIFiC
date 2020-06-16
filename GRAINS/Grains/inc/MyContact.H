#ifndef _MyContact
#define _MyContact

#include "ContactLaw.hh"
#include "Basic.H"
#include "Point.H"
#include "PointContact.hh"
#include "Vecteur.H"
#include "ReaderXML.hh"

#include <map>
using namespace std;

class Composant;


/** @brief Loi de contact collision/friction avec comme parametre (entre
autres) le coefficient de dissipation normal muen.

    @author A.WACHS - IFPEN - 2011 - Creation */
// ============================================================================
class MyContact : public ContactLaw
{
public:
  /** @name Constructeurs */
  //@{
  /** @brief Constructeur avec un lot de parametres
  @param parameters Les parametres du contact */
  MyContact(map<string,double>& parameters);

  /** @brief Destructeur. */
  virtual ~MyContact();
  //@}


  /**@name Methods */
  //@{
  /** @brief Nom de la classe de contact. */
  string name() const { return "MyContact"; }

  /** @brief Calcul une estimation du temps de contact et de la penetration
  maximale d'un choc elastique binaire avec dissipation d'energie, puis ecrit
  le resultat dans un flux de sortie
  @param p0_ Premier Composant (Particule)
  @param p1_ Deuxieme Composant (Particule ou Obstacle)
  @param v0 vitesse relative au moment du contact
  @param OUT flux de sortie */
  void computeAndWriteEstimates( Composant* p0_,
	Composant* p1_,
  	const double &v0,
	ostream &OUT ) const;

  /** @brief Calcul des forces & moments de contact
  @param p0_ Premier Composant (Particule)
  @param p1_ Deuxieme Composant (Particule ou Obstacle)
  @param contactInfos Points & Recouvrement entre les deux Composants
  @param LC grille de cellules
  @param dt Increment de temps du pas de calcul
  @param nbContact number of contact */
  bool computeForces( Composant* p0_,
		     Composant* p1_,
		     PointContact &contactInfos,
		     LinkedCell *LC,
		     Scalar dt, int nbContact = 1  ) ;

  /*
    @brief Initialisation des proprietes cohesives
    @param p0_ Premier composant (Particule)
    @param p1_ Deuxieme composant (Particule ou Obstacle)
    @param contactInfos Points & Recouvrement entre les deux Composants
    @param LC grille de cellules
    @param dt Increment de temps du pas de calcul
    @param nbContact number of contact
  bool InitializeCohesiveObjects( Composant* p0_, Composant* p1_,
      const PointContact &contactInfos, LinkedCell *LC,
      Scalar dt, int nbContact = 1 ) { return false; }
*/

  /** @brief Calcul des forces & moments de contact pour le post processing
  @param p0_ Premier Composant (Particule)
  @param p1_ Deuxieme Composant (Particule ou Obstacle)
  @param contactInfos Points & Recouvrement entre les deux Composants
  @param listOfContacts liste de contacts */
  void computeForcesPostProcessing( Composant* p0_,
		     Composant* p1_, Scalar dt,
		     PointContact &contactInfos,
		     list<PointForcePostProcessing>* listOfContacts ) ;
  //@}


  /** @name Static Methods */
  //@{
  /** @brief Definition des parametres pour MyContact
  @param root le noeud "<MyContact>"
  @return valeurs des param�tres */
  static map<string,double> defineParameters(DOMNode* &root);
  //@}


protected:
  /** @name Parameters */
  //@{
  Scalar stiff; /**< Coefficient d'elasticite */
  Scalar muen; /**< Coefficient de frottement direction normale */
  Scalar muet; /**< Coefficient de frottement direction tangente */
  Scalar muec; /**< Coefficient de frottement de coulomb */
  Scalar k_s; /**< Coefficient de rugosite */
  Scalar k_m_s; /**< Coefficient */
  Vecteur m_forceRugosite; /**< Histoire de la force de rugosite du contact */
  //@}


  /**@name Constructeurs */
  //@{
  /** @brief Constructeur par defaut interdit. */
  MyContact();
  //@}
};

#endif
