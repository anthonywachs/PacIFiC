#include "PointContact.hh"
#include "Component.hh"


// ----------------------------------------------------------------------------
// Default constructor
PointContact::PointContact()
{}




// ----------------------------------------------------------------------------
// Constructor with contact point location in the world reference
// frame, overlap vector, overlap distance and number of iterations of GJK as
// input parameters
PointContact::PointContact( Point3 const& point_, Vector3 const& ov_,
	double distance_, int num_iter_ )
  : m_contact( point_ )
  , m_overlapvector( ov_ )
  , m_overlapdistance( distance_ )
  , m_nbIterGJK( num_iter_ )
{}




// ----------------------------------------------------------------------------
// Constructor with contact point location in the world reference
// frame, the point in the 1st rigid body that realizes the minimal distance in
// the rigid body reference frame, the point in the 2nd rigid body that
// realizes the minimal distance in the rigid body reference frame, overlap
// vector, overlap distance and number of iterations of GJK as input parameters
PointContact::PointContact( Point3 const& point_, Point3 const& pointA_,
    	Point3 const& pointB_, Vector3 const& ov_,
	double distance_, int num_iter_ )
  : m_contact( point_ )
  , m_contactA( pointA_ )
  , m_contactB( pointB_ )
  , m_overlapvector( ov_ )
  , m_overlapdistance( distance_ )
  , m_nbIterGJK( num_iter_ )
{}




// ----------------------------------------------------------------------------
// Copy constructor
PointContact::PointContact( PointContact const& pc_ )
{
  m_contact = pc_.m_contact;
  m_contactA = pc_.m_contactA;
  m_contactB = pc_.m_contactB;
  m_overlapvector = pc_.m_overlapvector;
  m_overlapdistance = pc_.m_overlapdistance;
  m_nbIterGJK = pc_.m_nbIterGJK;
}




// ----------------------------------------------------------------------------
// Destructor
PointContact::~PointContact()
{}




// ----------------------------------------------------------------------------
// Returns the contact point in the world reference frame
Point3 PointContact::getContact() const
{
  return ( m_contact );
}




// ----------------------------------------------------------------------------
// Returns the point in the 1st rigid body that realizes the minimal
// distance in the rigid body reference frame
Point3 PointContact::getContactA() const
{
  return ( m_contactA );
}




// ----------------------------------------------------------------------------
// Returns the point in the 2nd rigid body that realizes the minimal
// distance in the rigid body reference frame
Point3 PointContact::getContactB() const
{
  return ( m_contactB );
}




// ----------------------------------------------------------------------------
// Returns the overlap distance
double PointContact::getOverlapDistance() const
{
  return ( m_overlapdistance );
}




// ----------------------------------------------------------------------------
// Returns the overlap vector
Vector3 PointContact::getOverlapVector() const
{
  return ( m_overlapvector );
}




// ----------------------------------------------------------------------------
// Returns the number of iterations of GJK for convergence
int PointContact::getNbIterGJK() const
{
  return ( m_nbIterGJK );
}




// ----------------------------------------------------------------------------
// Sets the contact point in the world reference frame
void PointContact::setContact( Point3 const& point )
{
  m_contact = point;
}




// ----------------------------------------------------------------------------
// Sets the overlap distance
void PointContact::setOverlapDistance( double dist )
{
  m_overlapdistance = dist;
}




// ----------------------------------------------------------------------------
// Sets the overlap vector
void PointContact::setOverlapVector( Vector3 const& vec )
{
  m_overlapvector = vec;
}




// ----------------------------------------------------------------------------
// Equal operator to another PointContact object
PointContact& PointContact::operator = ( PointContact const& rhs )
{
  if ( this != &rhs )
  {
    m_contact = rhs.m_contact;
    m_contactA = rhs.m_contactA;
    m_contactB = rhs.m_contactB;
    m_overlapvector = rhs.m_overlapvector;
    m_overlapdistance = rhs.m_overlapdistance;
    m_nbIterGJK = rhs.m_nbIterGJK;
  }
  return ( *this );
}
