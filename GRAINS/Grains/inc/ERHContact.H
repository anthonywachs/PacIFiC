#ifndef _ERHContact
#define _ERHContact

#include "ContactLaw.hh"
#include "Basic.H"
#include "Point.H"
#include "PointContact.hh"
#include "Vecteur.H"
#include "ReaderXML.hh"
#include <map>
using namespace std;

class Composant;


/** @brief xxx.

    @author A.WACHS - UBC - 2017 - Creation */
// ============================================================================
class ERHContact : public ContactLaw
{
public:
  /** @name Constructeurs */
  //@{
  /** @brief Constructeur avec un lot de parametres
  @param parameters Les parametres du contact */
  ERHContact( map<string,double>& parameters );

  /** @brief Destructeur. */
  virtual ~ERHContact();
  //@}


  /**@name Methods */
  //@{
  /** @brief Nom de la classe de contact. */
  string name() const { return "ERHContact"; }

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

  /** @brief Computation of vector tangent to contact according to Costa (2015)
  @param tij The vector where the result is stored
  @param n_t Eta_t, the tangential damping coefficient
  @param u_t The tangential velocity
  @param kdelta The cumulative displacement */
  void computeTangentialVector(Vecteur& tij, double n_t, const Vecteur ut,
        const Vecteur kdelta);

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
  /** @brief Definition des parametres pour ERHContact
  @param root le noeud "<ERHContact>"
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
  Scalar k_m_s; /**< Ancien coefficient de friction de roulement x rayon */
  Scalar eta_r; /**< Viscous rolling damping ratio */
  Scalar mu_r; /**< Rolling resistance coefficient */
  Scalar Jn; /**< dimensionless coefficient used to define
                  k_m_s (Ai et al. 2011) */
  Scalar m_f; /**< binary value: if f=0, no viscous rolling damping when the
                 spring rolling friction is saturated */
  Scalar ks; /**< Coefficient de frottement statique */
  Scalar epsilon; /**< Criterion on velocity norm for cumulative tangential vector */
  bool rolling_friction;
  //@}


  /**@name Constructeurs */
  //@{
  /** @brief Constructeur par defaut interdit. */
  ERHContact();
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
    @brief Computation of contact forces & torques with memory term
    @param p0_ First Composant (Particule)
    @param p1_ Second Composant (Particule or Obstacle)
    @param recouvrement Points & Recouvrement between the two Composants
    @param delFN normal force
    @param delFT tangential force
    @param delM torque
    @param elementary_id0 Id of elementary particle in case of Composite
    @param elementary id1 Id of elementary particle in case of Composite
  */
  void performForcesCalculus(
        Composant* p0_,
	Composant* p1_,
  Scalar dt,
	PointContact &recouvrement,
	Vecteur &delFN,
	Vecteur &delFT,
	Vecteur &delM,
  int elementary_id0,
  int elementary_id1,
  bool forget_contact=false);

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
	PointContact &recouvrement,
	Vecteur &delFN,
	Vecteur &delFT,
	Vecteur &delM,
	Scalar &diss );
  //@}
};

#endif
