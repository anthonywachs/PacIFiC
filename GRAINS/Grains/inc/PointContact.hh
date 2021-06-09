#ifndef _POINTCONTACT
#define _POINTCONTACT

#include "Point.H"
#include "Vecteur.H"
using namespace solid;

class Composant;


/** @brief Container decrivant le point de contact entre les Composant A & B.

    Le container possede l'ensemble des informations geometriques au point
    de contact. En particulier, il decrit le point de contact dans les trois
    espaces disponibles : "world coordinate", convex A, convex B. On passe
    d'un espace a l'autre par les matrices de transformation.

    @author GRAINS Project - IFP - 2007 - Creation */
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class PointContact
{
public:
  /** @name Constructors */
  //@{
  /** @brief Constructeur par defaut */
  PointContact();

  /** @brief Constructeur avec initialisation
  @param point_  Point de contact dans le WC
  @param recouvrement_ Vecteur recouvrement A-B au point de contact
  @param distance_ Longueur du vecteur A-B
  @param num_iter_ nombre d'iterations de GJK */
  PointContact(const Point &point_,
	       const Vecteur &recouvrement_,
	       double distance_, int num_iter_, int NbCoulombRegime=0);

  /** @brief Constructeur avec initialisation
  @param point_  Point de contact dans le WC
  @param pointA_ Point de contact dans le Convex A
  @param pointB_ Point de contact dans le Convex B
  @param recouvrement_ Vecteur recouvrement A-B au point de contact
  @param distance_ Longueur du vecteur A-B
  @param num_iter_ nombre d'iterations de GJK */
  PointContact(const Point &point_, const Point &pointA_, const Point &pointB_,
	       const Vecteur &recouvrement_, double distance_, int num_iter_,
           int NbCoulombRegime=0);

  /** @brief Constructeur par copie.
  @param pc_ le point contact de reference */
  PointContact(const PointContact &pc_);

  /** @brief Destructeur */
  ~PointContact();
  //@}


  /** @name Methods Get */
  //@{
  /** @brief Point de contact en WC
  @return Les coordonnees du point */
  Point getContact() const;

  /** @brief Point de contact dans l'espace du Convex A
  @return Les coordonnnees du point */
  Point getContactA() const;

  /** @brief Point de contact dans l'espace du Convex B
  @return Les coordonnnees du point */
  Point getContactB() const;

  /** @brief Distance entre les Convex A vers B
  @return La distance */
  double getDistance() const;

  /** @brief Recouvrement entre les Convex A & B
  @return Le vecteur de recouvrement */
  Vecteur getRecouvrement() const;

  /** @brief Nombre d'iterations de GJK
  @return Le nombre d'iterations */
  int getNbIterGJK() const;
  //@}

  int getNbCoulombRegimes() const;

  void incrementNbCoulombRegimes();

  void setNbCoulombRegimes(const int NbCoulombRegimes);

  /** @name Methods Set */
  //@{
  /** @brief Point de contact en WC
  @param point Les coordonnees du point */
  void setContact(const Point &point);

  /** @brief Point de contact dans l'espace du Convex A
  @param point Les coordonnnees du point */
  void setContactA(const Point &point);

  /** @brief Point de contact dans l'espace du Convex B
  @param point Les coordonnnees du point */
  void setContactB(const Point &point);

  /** @brief Distance entre les deux Convex
  @param dist La distance */
  void setDistance(double dist);

  /** @brief Recouvrement entre les Convex A & B
  @param vecteur Le vecteur de recouvrement */
  void setRecouvrement(const Vecteur &vecteur);
  //@}


  /**@name Operators */
  //@{
  /** @brief Operateur d'affectation
  @param rhs point contact � affecter */
  PointContact& operator= (const PointContact& rhs);
  //@}


private:
  /** @name Parameters */
  //@{
  Point m_contact; /**< Point de contact en "World Coordinate" :
      contact = a2w(contactA) & contact = b2w(contactB) */
  Point m_contactA; /**< Point de contact dans l'espace Convex A :
      contactA = w2a(contact) */
  Point m_contactB; /**< Point de contact dans l'espace Convex B :
      contactB = w2b(contact) */
  Vecteur m_recouvrement; /**< Recouvrement du contact en A & B */
  double m_distance; /**< Longueur du vecteur recouvrement A-B */
  int m_nbIterGJK; /**< Nombre d'iterations de GJK */
  int m_nbCoulombRegimes; /**< Nombre de fois où le contact est en régime de Coulomb */
  // int NbStaticRegime;
  //@}
};

static PointContact PointNoContact( OriginePoint, OriginePoint, OriginePoint,
	VecteurNul, 1.e20, 0, 0 );

#endif
