#ifndef _ERContact
#define _ERContact

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
autres) le coefficient de restitution en.

    @author A.WACHS - IFPEN - 2011 - Creation */
// ============================================================================
class ERContact : public ContactLaw
{
public:
  /** @name Constructeurs */
  //@{
  /** @brief Constructeur avec un lot de parametres
  @param parameters Les parametres du contact */
  ERContact( map<string,double>& parameters );

  /** @brief Destructeur. */
  virtual ~ERContact();
  //@}


  /**@name Methods */
  //@{
  /** @brief Nom de la classe de contact. */
  string name() const { return "ERContact"; }

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
  @param dt Increment de temps du pas de calcul.
  @param nbContact Nombre de points de contact pour les particules
  composites */
  bool computeForces( Composant* p0_,
		     Composant* p1_,
		     PointContact &contactInfos,
		     LinkedCell *LC,
		     Scalar dt, int nbContact = 1 ) ;

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
  /** @brief Definition des parametres pour ERContact
  @param root le noeud "<ERContact>"
  @return valeurs des param�tres */
  static map<string,double> defineParameters( DOMNode* &root );
  //@}


protected:
  /** @name Parameters */
  //@{
  Scalar stiff; /**< Coefficient d'elasticite */
  Scalar en; /**< Coefficient de restitution */
  Scalar muet; /**< Coefficient de frottement direction tangente */
  Scalar muec; /**< Coefficient de frottement de coulomb */
  Scalar k_m_s; /**< Coefficient */
  //@}


  /**@name Constructeurs */
  //@{
  /** @brief Constructeur par defaut interdit. */
  ERContact();
  //@}

  /**@name Methods */
  //@{
  /** @brief Calcul de la penetration maximale sur la base de la solution
  analytique et d'un algo de Newton
  @param theta_ racine( omega0*omega0 - mu*mu )
  @param mu_ coefficient de dissipation
  @param en_ coefficient de restitution
  @param tc_ temps de contact
  @param v0_ vitesse relative */
  double computeDeltaMax( const double &theta_, const double& mu_,
  	const double& en_, const double &tc_, const double &v0_ ) const;

  /**
    @brief Calcul effectif des forces & moments de contact
    @param p0_ Premier Composant (Particule)
    @param p1_ Deuxieme Composant (Particule ou Obstacle)
    @param recouvrement Points & Recouvrement entre les deux Composants
    @param delFN force normale
    @param delFT force tangentielle
    @param delM moment
  */
  void performForcesCalculus(
        Composant* p0_,
	Composant* p1_,
	const PointContact &recouvrement,
	Vecteur &delFN,
	Vecteur &delFT,
	Vecteur &delM );

  /** @brief Compute the contact force for postprocessing
       it also computes the energy dissipation rate
  @param p0_ Premier Composant (Particule)
  @param p1_ Deuxieme Composant (Particule ou Obstacle)
  @param recouvrement Points & Recouvrement entre les deux Composants
  @param delFN force normale
  @param delFT force tangentielle
  @param delM moment */
  void performForcesCalculusPP(
        Composant* p0_,
	Composant* p1_,
	const PointContact &recouvrement,
	Vecteur &delFN,
	Vecteur &delFT,
	Vecteur &delM,
	Scalar &diss  );
  //@}
};

#endif
