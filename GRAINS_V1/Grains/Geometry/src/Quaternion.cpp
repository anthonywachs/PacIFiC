#include "GrainsMPIWrapper.hh"
#include "GrainsExec.hh"
#include "Quaternion.hh"
#include "Matrix.hh"
#include <math.h>


size_t Quaternion::m_sizeofQuaternion = solid::Group3::m_sizeofGroup3
	+ sizeof( double );

// ----------------------------------------------------------------------------
// Default constructor
Quaternion::Quaternion()
  : m_w( 1. )
  , m_vqt( 0. )
{}




// ----------------------------------------------------------------------------
// Constructor with 2 scalar as input parameters q and d. Quaternion
// is initialized as [ d, (q,q,q) ]
Quaternion::Quaternion( double q, double d )
  : m_w( d )
  , m_vqt( q )
{}




// ----------------------------------------------------------------------------
// Constructor with a Vector3 vector vec and a scalar d. Quaternion
// is initialized as [ d, vec ]
Quaternion::Quaternion( Vector3 const &vec, double d )
  : m_w( d )
  , m_vqt( vec )
{}




// ----------------------------------------------------------------------------
// Constructor with a vector given by its 3 components (x,y,z) and a
// scalar d. Quaternion is initialized as [ d, (x,y,z) ]
Quaternion::Quaternion( double x, double y, double z, double d )
  : m_w( d )
  , m_vqt( Vector3( x, y, z ) )
{}




// ----------------------------------------------------------------------------
// Copy constructor
Quaternion::Quaternion( Quaternion const& q )
  : m_w( q.m_w )
  , m_vqt( q.m_vqt )
{}




// ----------------------------------------------------------------------------
// Destructor
Quaternion::~Quaternion()
{}




// ----------------------------------------------------------------------------
// Returns the conjugate of the quaternion
Quaternion Quaternion::Conjugate() const
{
  return ( Quaternion( -m_vqt, m_w ) );
}




// ----------------------------------------------------------------------------
// Returns the inverse of the quaternion
Quaternion Quaternion::Inverse() const
{
  return ( Conjugate() * ( 1. / Norm(*this) ) );
}




// ----------------------------------------------------------------------------
// Multiplies the quaternion on the left by a vector lhs, i.e., perform
// [ 0, lhs ] x this and return the product that is a quaternion
Quaternion Quaternion::multLeftVec( Vector3 const& lhs ) const
{
  double tmp = -lhs * m_vqt;
  Vector3 vtmp = ( lhs ^ m_vqt ) + ( lhs * m_w ) ;
  return ( Quaternion( vtmp, tmp ) );
}




// ----------------------------------------------------------------------------
// Multiplies the quaternion on the right by another quaternion rhs,
// i.e., perform this x rhs, and return the vectorial part of this x rhs
Vector3 Quaternion::multToVector3( const Quaternion& rhs ) const
{
  Vector3 vtmp = ( m_vqt ^ rhs.m_vqt ) + ( m_w * rhs.m_vqt )
  	+ ( rhs.m_w * m_vqt );
  return ( vtmp );
}




// ----------------------------------------------------------------------------
// Multiplies the quaternion on the right by the conjugate of another
// quaternion rhs, i.e., perform this x rhs^t, and return the vectorial part of
// this x rhs^t
Vector3 Quaternion::multConjugateToVector3( Quaternion const& rhs ) const
{
  Vector3 vtmp = - ( m_vqt ^ rhs.m_vqt ) - ( m_w * rhs.m_vqt )
  	+ ( rhs.m_w * m_vqt );
  return ( vtmp );
}




// ----------------------------------------------------------------------------
// Returns the value of the scalar part of the quaternion
double Quaternion::getdouble() const
{
  return ( m_w );
}




// ----------------------------------------------------------------------------
// Returns the pointer to the vectorial part of the quaternion
Vector3 const* Quaternion::getVector3() const
{
  return ( &m_vqt );
}




// ----------------------------------------------------------------------------
// Sets the quaternion with a Vector3 vector vec and a scalar d.
// Quaternion is set to [ d, vec ]
void Quaternion::setQuaternion( Vector3 const& vec, double d )
{
  m_w = d;
  m_vqt = vec;
}




// ----------------------------------------------------------------------------
// Sets the quaternion with a vector given by its 3 components (x,y,z)
// and a scalar d. Quaternion is set to [ d, (x,y,z) ]
void Quaternion::setQuaternion( double x, double y, double z, double d )
{
  m_vqt[X] = x;
  m_vqt[Y] = y;
  m_vqt[Z] = z;
  m_w = d;
}




// ----------------------------------------------------------------------------
// Sets the scalar part of the quaternion
void Quaternion::setdouble( double d )
{
  m_w = d;
}




// ----------------------------------------------------------------------------
// Sets the vectorial part of the quaternion
void Quaternion::setVector3( Vector3 const& vec )
{
  m_vqt = vec;
}




// ----------------------------------------------------------------------------
// Sets the quaternion with a rotation matrix
void Quaternion::setQuaternion( Matrix const& rot )
{
  if ( !rot.isRotation() )
    cout << "Matrix is not a rotation in Quaternion::setQuaternion"
    	"( Matrix const& rot )" << endl;
  else
  {
    double den = 0.;

    // Case rotYY > - rotZZ, rotXX > - rotYY and rotXX > - rotZZ
    if ( rot[Y][Y] > - rot[Z][Z] && rot[X][X] > - rot[Y][Y]
    	&& rot[X][X] > - rot[Z][Z] )
    {
      den = pow( 1. + rot[X][X] + rot[Y][Y] + rot[Z][Z], 0.5 );
      m_w = 0.5 * den;
      m_vqt[X] = 0.5 * ( rot[Z][Y] - rot[Y][Z] ) / den;
      m_vqt[Y] = 0.5 * ( rot[X][Z] - rot[Z][X] ) / den;
      m_vqt[Z] = 0.5 * ( rot[Y][X] - rot[X][Y] ) / den;
    }
    // Case rotYY < - rotZZ, rotXX > rotYY and rotXX > rotZZ
    else if ( rot[Y][Y] < - rot[Z][Z] && rot[X][X] > rot[Y][Y]
    	&& rot[X][X] > rot[Z][Z] )
    {
      den = pow( 1. + rot[X][X] - rot[Y][Y] - rot[Z][Z], 0.5 );
      m_w = 0.5 * ( rot[Z][Y] - rot[Y][Z] ) / den;
      m_vqt[X] = 0.5 * den;
      m_vqt[Y] = 0.5 * ( rot[X][Y] + rot[Y][X] ) / den;
      m_vqt[Z] = 0.5 * ( rot[Z][X] + rot[X][Z] ) / den;
    }
    // Case rotYY > rotZZ, rotXX < rotYY and rotXX < - rotZZ
    else if ( rot[Y][Y] > rot[Z][Z] && rot[X][X] < rot[Y][Y]
    	&& rot[X][X] < - rot[Z][Z] )
    {
      den = pow( 1. - rot[X][X] + rot[Y][Y] - rot[Z][Z], 0.5 );
      m_w = 0.5 * ( rot[X][Z] - rot[Z][X] ) / den;
      m_vqt[X] = 0.5 * ( rot[X][Y] + rot[Y][X] ) / den;
      m_vqt[Y] = 0.5 * den;
      m_vqt[Z] = 0.5 * ( rot[Y][Z] + rot[Z][Y] ) / den;
    }
    // Case rotYY < rotZZ, rotXX < - rotYY and rotXX < rotZZ
    else if ( rot[Y][Y] < rot[Z][Z] && rot[X][X] < - rot[Y][Y]
    	&& rot[X][X] < rot[Z][Z] )
    {
      den = pow( 1. - rot[X][X] - rot[Y][Y] + rot[Z][Z], 0.5 );
      m_w = 0.5 * ( rot[Y][X] - rot[X][Y] ) / den;
      m_vqt[X] = 0.5 * ( rot[Z][X] + rot[X][Z] ) / den;
      m_vqt[Y] = 0.5 * ( rot[Y][Z] + rot[Z][Y] ) / den;
      m_vqt[Z] = 0.5 * den;
    }
    else
      cout << "Warning: case not covered in Quaternion::setQuaternion"
    	"( Matrix const& rot )" << endl;
  }
}




// ----------------------------------------------------------------------------
// Build a unit quaternion representing the rotation from u to v. 
// The input vectors need not be normalised. */
// TODO: if the input vectors aren't normalized, normalize them and warn the
// user.
void Quaternion::setRotFromTwoVectors( Vector3 const& u, Vector3 const& v )
{
  double norm_u_norm_v = sqrt( (u*u) * (v*v) );
  double real_part = norm_u_norm_v + u*v;
  Vector3 vect;

  if ( real_part < 1.e-6 * norm_u_norm_v )
  {
    /* If u and v are exactly opposite, rotate 180 degrees
    around an arbitrary orthogonal axis. Axis normalisation
    can happen later, when we normalise the quaternion. */
    real_part = 0.;
    vect = fabs(u[0]) > fabs(u[2]) ? Vector3( -u[1], u[0], 0. )
                                : Vector3( 0., -u[2], u[1] );
  }
  else
  {
    /* Otherwise, build quaternion the standard way. */
    vect = u ^ v;
  }

  Quaternion qq( vect[0], vect[1], vect[2], real_part );
  *this = qq * ( 1. / Norm( qq ) );
}




// ----------------------------------------------------------------------------
// Rotates a vector using the quaternion *this
Vector3 Quaternion::rotateVector( Vector3 const& v ) const
{
  Vector3 v_rotated = ( m_w * m_w - Norm(m_vqt) * Norm(m_vqt) ) * v 
  	+ 2. * ( v * m_vqt ) * m_vqt + 2. * m_w * ( m_vqt ^ v );
  return ( v_rotated );
}




// ----------------------------------------------------------------------------
// Equal operator to another quaternion
Quaternion& Quaternion::operator = ( Quaternion const& rhs )
{
  if ( &rhs != this )
  {
    m_w   = rhs.m_w;
    m_vqt = rhs.m_vqt;
  }
  return ( *this );
}




// ----------------------------------------------------------------------------
// Equal operator to a scalar d, the result is [ d, (d,d,d) ]
Quaternion Quaternion::operator = ( const double d )
{
  m_w = d;
  m_vqt = d;
  return ( *this );
}




// ----------------------------------------------------------------------------
// ith-component accessor: (0,1,2) for the vector components and 3 for
// the scalar
double& Quaternion::operator [] ( int i )
{
  return ( i == 3 ? m_w : m_vqt[i] );
}




// ----------------------------------------------------------------------------
// ith-component accessor: (0,1,2) for the vector components and 3 for
// the scalar
double Quaternion::operator [] ( int i ) const
{
  return ( i == 3 ? m_w :  m_vqt[i] );
}




// ----------------------------------------------------------------------------
// Unitary operator -. Return a quaternion with negative elements
Quaternion Quaternion::operator - ()
{
  return ( Quaternion( - m_vqt, - m_w ) );
}




// ----------------------------------------------------------------------------
// Operator +=
Quaternion& Quaternion::operator += ( Quaternion const& rhs )
{
  m_w += rhs.m_w;
  m_vqt += rhs.m_vqt;
  return ( *this );
}




// ----------------------------------------------------------------------------
// Sum of 2 quaternions
Quaternion Quaternion::operator + ( Quaternion const& rhs ) const
{
  return ( Quaternion( m_vqt + rhs.m_vqt, m_w + rhs.m_w ) );
}




// ----------------------------------------------------------------------------
// Operator -=
Quaternion& Quaternion::operator -= ( Quaternion const& rhs )
{
  m_w -= rhs.m_w;
  m_vqt -= rhs.m_vqt;
  return ( *this );
}




// ----------------------------------------------------------------------------
// Subtraction of 2 quaternions, i.e., compute this - rhs
Quaternion Quaternion::operator - ( Quaternion const& rhs )
{
  return ( Quaternion( m_vqt - rhs.m_vqt, m_w - rhs.m_w ) );
}




// ----------------------------------------------------------------------------
// Unitary operator *= by a scalar
Quaternion& Quaternion::operator *= ( double d )
{
  m_w *= d;
  m_vqt *= d;
  return ( *this );
}




// ----------------------------------------------------------------------------
// Product by a scalar
Quaternion Quaternion::operator * ( double d )
{
  return ( Quaternion( d * m_vqt, d * m_w ) );
}




// ----------------------------------------------------------------------------
// Product by a scalar
Quaternion operator * ( double d, Quaternion const& rhs )
{
  Quaternion result( rhs );
  result *= d;
  return ( result );
}




// ----------------------------------------------------------------------------
// double product this x rhs of 2 quaternions
Quaternion Quaternion::operator * ( Quaternion const& rhs ) const
{
  double tmp = ( m_w * rhs.m_w ) - ( m_vqt * rhs.m_vqt );
  Vector3 vtmp = ( m_vqt ^ rhs.m_vqt ) + ( m_w * rhs.m_vqt )
  	+ ( rhs.m_w * m_vqt );
  return ( Quaternion( vtmp, tmp ) );
}




// ----------------------------------------------------------------------------
// double product on the right of a quaternion by a vector [ 0, rhs ]
Quaternion Quaternion::operator , ( Vector3 const& rhs ) const
{
  double tmp = - m_vqt * rhs;
  Vector3 vtmp = ( m_vqt ^ rhs ) + ( m_w * rhs );
  return ( Quaternion( vtmp, tmp ) );
}




// ----------------------------------------------------------------------------
// double product on the left of a vector [ 0, lhs ] by a quaternion.
// Returns a quaternion.
Quaternion operator , ( Vector3 const& lhs, Quaternion const& q )
{
  return ( q.multLeftVec( lhs ) );
}




// ----------------------------------------------------------------------------
// Comparison operator
bool Quaternion::operator == ( Quaternion const& rhs )
{
  return ( m_w == rhs.m_w
  	&& m_vqt[0] == rhs.m_vqt[0]
	&& m_vqt[1] == rhs.m_vqt[1]
	&& m_vqt[2] == rhs.m_vqt[2] );
}




// ----------------------------------------------------------------------------
// Difference operator
bool Quaternion::operator != ( Quaternion const& rhs )
{
  return ( ! ( *this == rhs ) );
}




// ----------------------------------------------------------------------------
// Output operator
ostream& operator << ( ostream& fileOut, Quaternion const& q )
{
  fileOut << q.m_w << " " << q.m_vqt;
  return ( fileOut );
}




// ----------------------------------------------------------------------------
// Input operator
istream& operator >> ( istream& fileIn, Quaternion& q )
{
  fileIn >> q.m_w >> q.m_vqt;
  return ( fileIn );
}




// ----------------------------------------------------------------------
// Writes the object with a high precision format given by
// FORMAT16DIGITS defined in GrainsExec.hh
void Quaternion::writeQuaternion( ostream &fileOut ) const
{
  fileOut << GrainsExec::doubleToString( ios::scientific, FORMAT16DIGITS,
  	m_w ) << " ";
  m_vqt.writeGroup3( fileOut );
}




// ----------------------------------------------------------------------
// Writes the object in binary format
void Quaternion::writeQuaternion_binary( ostream &fileOut )
{
  fileOut.write( reinterpret_cast<char*>( &m_w ), sizeof( double ) );
  m_vqt.writeGroup3_binary( fileOut );
}




// ----------------------------------------------------------------------
// Reads the object in binary format
void Quaternion::readQuaternion_binary( istream &StreamIN )
{
  StreamIN.read( reinterpret_cast<char*>( &m_w ), sizeof( double ) );
  m_vqt.readGroup3_binary( StreamIN );
}




// ----------------------------------------------------------------------------
// Returns the norm of the quaternion
double Norm( Quaternion const& q )
{
  return ( sqrt( q.m_vqt[X] * q.m_vqt[X] + q.m_vqt[Y] * q.m_vqt[Y]
  	+ q.m_vqt[Z] * q.m_vqt[Z] + q.m_w * q.m_w ) );
}




// ----------------------------------------------------------------------------
// Returns the norm square of the quaternion
double Norm2( Quaternion const& q )
{
  return ( q.m_vqt[X] * q.m_vqt[X] + q.m_vqt[Y] * q.m_vqt[Y]
  	+ q.m_vqt[Z] * q.m_vqt[Z] + q.m_w * q.m_w );
}
