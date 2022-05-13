#ifndef _ContactLaw
#define _ContactLaw

#include "Basic.H"
#include "Point.H"
#include "PointContact.hh"
#include "Vecteur.H"
#include "ReaderXML.hh"
#include "LinkedCell.H"
#include <list>
#include <map>
using namespace std;

class Composant;
class LinkedCell;


/** @brief Loi de contact.

    @author A.WACHS - IFPEN - 2011 - Creation */
// ============================================================================
class ContactLaw
{
public:
  /** @name Constructeurs */
  //@{
  /** @brief Constructeur avec un lot de parametres
  @param parameters Les parametres du contact */
  ContactLaw( map<string,double>& parameters );

  /** @brief Destructeur. */
  virtual ~ContactLaw();
  //@}


  /**@name Methods */
  //@{
  /** @brief Nom de la classe de contact. */
  virtual string name() const=0;

  /** @brief Calcul une estimation du temps de contact et de la penetration
  maximale d'un choc elastique binaire avec dissipation d'energie, puis ecrit
  le resultat dans un flux de sortie
  @param p0_ Premier Composant (Particule)
  @param p1_ Deuxieme Composant (Particule ou Obstacle)
  @param v0 vitesse relative au moment du contact
  @param OUT flux de sortie */
  virtual void computeAndWriteEstimates( Composant* p0_,
	Composant* p1_,
  	const double &v0,
	ostream &OUT ) const = 0 ;

  /** @brief Calcul des forces & moments de contact
  @param p0_ Premier Composant (Particule)
  @param p1_ Deuxieme Composant (Particule ou Obstacle)
  @param contactInfos Points & Recouvrement entre les deux Composants
  @param LC grille de cellules
  @param dt Increment de temps du pas de calcul
  @param nbContact nbContact number of contact with the particle */
  virtual bool computeForces( Composant* p0_,
		     Composant* p1_,
		     PointContact &contactInfos,
		     LinkedCell *LC,
		     Scalar dt, int nbContact = 1 ) = 0 ;

  /**
    @brief Initialisation des proprietes cohesives
    @param p0_ Premier composant (Particule)
    @param p1_ Deuxieme composant (Particule ou Obstacle)
    @param contactInfos Points & Recouvrement entre les deux Composants
    @param LC grille de cellules
    @param dt Increment de temps du pas de calcul
    @param nbContact number of contact
  */
  virtual bool InitializeCohesiveObjects( Composant* p0_, Composant* p1_,
      const PointContact &contactInfos, LinkedCell *LC,
      Scalar dt, int nbContact = 1 ) { return false; } ;

  /** @brief Calcul des forces & moments de contact pour le post processing
  @param p0_ Premier Composant (Particule)
  @param p1_ Deuxieme Composant (Particule ou Obstacle)
  @param contactInfos Points & Recouvrement entre les deux Composants
  @param listOfContacts liste de contacts */
  virtual void computeForcesPostProcessing( Composant* p0_,
		     Composant* p1_, Scalar dt,
		     PointContact &contactInfos,
		     list<PointForcePostProcessing>* listOfContacts ) = 0 ;
  //@}


  /**@name Static Methods */
  //@{
  /** @brief Contact entre 2 clones periodiques ?
  @param p0_ Premiere particule (clone periodique)
  @param p1_ Deuxieme Composant (Particule ou clone periodique) */
  static bool is_ClonePer_ClonePer( Composant const* p0_,
		     Composant const* p1_ );

  /** @brief Contact entre 1 clones periodique et une particule possedant un
  clone periodique dans la meme direction (sens oppose) ?
  @param p0_ Premiere particule (clone periodique)
  @param p1_ Deuxieme Composant (Particule ou clone periodique)
  @param LC grille de cellules */
  static bool is_ClonePer_ParticuleWithClonePerSameDirection(
  	Composant const* p0_,
	Composant const* p1_,
	LinkedCell const* LC );

  /** @brief Les directions periodiques font elles un angle de 90 ?
  @param p0_ Premiere particule (clone periodique)
  @param p1_ Deuxieme Composant (clone periodique) */
  static bool are_ClonePerDirections_Perp( Composant const* p0_,
		     Composant const* p1_ );
  //@}


protected:
  /**@name Constructeurs */
  //@{
  /** @brief Constructeur par defaut interdit. */
  ContactLaw() {};
  //@}
};

#endif
