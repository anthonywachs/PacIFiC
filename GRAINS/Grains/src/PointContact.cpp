#include "PointContact.hh"

#include "Composant.H"


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur par defaut
PointContact::PointContact()
{
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur avec initialisation
PointContact::PointContact(const Point &point_,
	const Vecteur &recouvrement_, double distance_, int num_iter_,
	int nb_coulomb_regimes_) :
  m_contact(point_),
  m_recouvrement(recouvrement_),
  m_distance(distance_),
  m_nbIterGJK(num_iter_),
  m_nbCoulombRegimes(nb_coulomb_regimes_)
{
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur avec initialisation
PointContact::PointContact(const Point &point_,
	const Point &pointA_, const Point &pointB_,
	const Vecteur &recouvrement_, double distance_, int num_iter_,
	int nb_coulomb_regimes_) :
  m_contact(point_),
  m_contactA(pointA_),
  m_contactB(pointB_),
  m_recouvrement(recouvrement_),
  m_distance(distance_),
  m_nbIterGJK(num_iter_),
  m_nbCoulombRegimes(nb_coulomb_regimes_)
{
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur par copie
PointContact::PointContact(const PointContact &pc_)
{
  m_contact = pc_.m_contact;
  m_contactA = pc_.m_contactA;
  m_contactB = pc_.m_contactB;
  m_recouvrement = pc_.m_recouvrement;
  m_distance = pc_.m_distance;
  m_nbIterGJK = pc_.m_nbIterGJK;
  m_nbCoulombRegimes = pc_.m_nbCoulombRegimes;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Destructeur
PointContact::~PointContact()
{
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//Point de contact en WC
Point PointContact::getContact() const
{
  return m_contact;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Point de contact dans l'espace du Convex A
Point PointContact::getContactA() const
{
  return m_contactA;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Point de contact dans l'espace du Convex B
Point PointContact::getContactB() const
{
  return m_contactB;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Distance entre les Convex A vers B
double PointContact::getDistance() const
{
  return m_distance;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Recouvrement entre les Convex A & B
Vecteur PointContact::getRecouvrement() const
{
  return m_recouvrement;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Nombre d'iterations de GJK
int PointContact::getNbIterGJK() const
{
  return m_nbIterGJK;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Nombre de fois que les frottements sont dans le régime de Coulomb
int PointContact::getNbCoulombRegimes() const
{
  return m_nbCoulombRegimes;
}



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Point de contact en WC
void PointContact::setContact(const Point &point)
{
  m_contact = point;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Point de contact dans l'espace du Convex A
void PointContact::setContactA(const Point &point)
{
  m_contactA = point;
}



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Nombre de fois que les frottements sont dans le régime de Coulomb
void PointContact::setNbCoulombRegimes(const int NbCoulombRegimes)
{
  m_nbCoulombRegimes = NbCoulombRegimes;
}



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Nombre de fois que les frottements sont dans le régime de Coulomb
void PointContact::incrementNbCoulombRegimes()
{
  m_nbCoulombRegimes++;
  // cout << m_nbCoulombRegimes << endl;
}



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Point de contact dans l'espace du Convex B
void PointContact::setContactB(const Point &point)
{
  m_contactB = point;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Distance entre les deux Convex
void PointContact::setDistance(double dist)
{
  m_distance = dist;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Recouvrement entre les Convex A & B
void PointContact::setRecouvrement(const Vecteur &vecteur)
{
  m_recouvrement = vecteur;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Operateur d'affectation
PointContact& PointContact::operator= (const PointContact& rhs)
{
  if ( this != &rhs )
  {
    m_contact = rhs.m_contact;
    m_contactA = rhs.m_contactA;
    m_contactB = rhs.m_contactB;
    m_recouvrement = rhs.m_recouvrement;
    m_distance = rhs.m_distance;
    m_nbIterGJK = rhs.m_nbIterGJK;
	m_nbCoulombRegimes = rhs.m_nbCoulombRegimes;
  }
  return *this;
}
