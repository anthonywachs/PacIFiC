#include "Vector3.hh"
#include "Quaternion.hh"
#include <math.h>


namespace solid
{
  // ---------------------------------------------------------------------------
  // Default constructor
  Vector3::Vector3( double def )
    : Group3( def )
  {}




  // ---------------------------------------------------------------------------
  // Constructor with 3 components as inputs
  Vector3::Vector3( double x, double y, double z )
    : Group3( x, y, z )
  {}




  // ---------------------------------------------------------------------------
  // Copy constructor
  Vector3::Vector3( Group3 const& point )
    : Group3( point )
  {}




  // ---------------------------------------------------------------------------
  // Destructor
  Vector3::~Vector3()
  {}




  // ---------------------------------------------------------------------------
  // Determines the direction of lowest absolute component
  int Vector3::closestAxis() const
  {
    double a[2];
    int axis = ( a[X] = fabs( m_comp[X] ) ) < ( a[Y] = fabs( m_comp[Y] ) )
    	? Y : X;
    return ( a[axis] < fabs( m_comp[Z] ) ? Z : axis );
  }




  // ---------------------------------------------------------------------------
  // Unitary nomalization operator
  void Vector3::normalize()
  {
    *this /= Norm( *this );
  }




  // ---------------------------------------------------------------------------
  // Returns a vector corresponding to the normalized vecteur
  Vector3 Vector3::normalized() const
  {
    return ( *this / Norm( *this ) );
  }




  // ---------------------------------------------------------------------------
  // Rotation by an unitary quaternion
  void Vector3::Rotate( Quaternion const& q )
  {
    Quaternion tmp( *this );
    tmp = (q * tmp ) * q.Conjugate();
    *this = *(tmp.getVector3());
  }




  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Equal operator to another Vector3 object
  Vector3& Vector3::operator = ( Vector3 const& g2 )
  {
    if ( &g2 != this )
    {
      m_comp[X] = g2.m_comp[X];
      m_comp[Y] = g2.m_comp[Y];
      m_comp[Z] = g2.m_comp[Z];
    }
    return (*this);
  }




  // ---------------------------------------------------------------------------
  // Cross product this x rhv
  Vector3 Vector3::operator ^ ( Vector3 const& rhv ) const
  {
    return ( Vector3( m_comp[1] * rhv.m_comp[2] - m_comp[2] * rhv.m_comp[1],
	- m_comp[0] * rhv.m_comp[2] + m_comp[2] * rhv.m_comp[0],
	m_comp[0] * rhv.m_comp[1] - m_comp[1] * rhv.m_comp[0] ) );
  }




  // ---------------------------------------------------------------------------
  // Returns whether the vector norm is less than EPSILON2
  // where EPSILON2 is defined in Basic.H
  bool approxZero( Vector3 const& v )
  {
    return ( Norm2( v ) < EPSILON2 );
  }




  // --------------------------------------------------------------------
  // Returns the cosine of the angle between 2 Vector3 objects
  double cos( Vector3 const& v1, Vector3 const& v2 )
  {
    return ( ( v1 * v2 ) / ( Norm( v1 ) * Norm( v2 ) ) );
  }




  // ---------------------------------------------------------------------------
  //Returns the norm of the vector
  double Norm( Vector3 const& v )
  {
    return ( sqrt( v.m_comp[X] * v.m_comp[X]
    	+ v.m_comp[Y] * v.m_comp[Y] + v.m_comp[Z] * v.m_comp[Z] ) );
  }




  // ---------------------------------------------------------------------------
  // Returns the norm square of the vector
  double Norm2( Vector3 const& v )
  {
    return ( v.m_comp[X] * v.m_comp[X] + v.m_comp[Y] * v.m_comp[Y]
    	+ v.m_comp[Z] * v.m_comp[Z] );
  }
}
