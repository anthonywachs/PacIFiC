#include "Vector3.hh"
#include "Quaternion.hh"
#include <math.h>


namespace solid
{
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
  // Rotation by an unitary quaternion
  void Vector3::Rotate( Quaternion const& q )
  {
    Quaternion tmp( *this );
    tmp = (q * tmp ) * q.Conjugate();
    *this = *(tmp.getVector3());
  }




  // --------------------------------------------------------------------
  // Returns the cosine of the angle between 2 Vector3 objects
  double cos( Vector3 const& v1, Vector3 const& v2 )
  {
    return ( ( v1 * v2 ) / ( Norm( v1 ) * Norm( v2 ) ) );
  }
}
