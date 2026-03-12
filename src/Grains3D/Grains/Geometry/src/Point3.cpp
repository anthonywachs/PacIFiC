#include "Point3.hh"
#include "Vector3.hh"
#include <math.h>


namespace solid
{
  // ---------------------------------------------------------------------------
  // Adds the same value to all 3 components, i.e., move the point 
  // by (dist,dist,dist)
  void Point3::Move( double dist )
  {
    for (int i=0; i<3; i++) m_comp[i] += dist;
  }




  // ---------------------------------------------------------------------------
  // Adds values to the 3 components using a 3-component array
  void Point3::Move( double const* dist )
  {
    for (int i=0; i<3; i++) m_comp[i] += dist[i];
  }




  // ---------------------------------------------------------------------------
  // Adds values to the 3 components using 3 scalars
  void Point3::Move( double distX, double distY, double distZ )
  {
    m_comp[0] += distX;
    m_comp[1] += distY;
    m_comp[2] += distZ;
  }
}
