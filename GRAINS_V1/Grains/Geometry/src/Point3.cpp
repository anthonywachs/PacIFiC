#include "Point3.hh"
#include "Vector3.hh"
#include <math.h>


namespace solid
{
  // ---------------------------------------------------------------------------
  // Default constructor
  Point3::Point3( double def )
    : Group3( def )
  {}




  // ---------------------------------------------------------------------------
  // Constructor with 3 components as inputs 
  Point3::Point3( double x, double y, double z )
    : Group3( x, y, z )
  {}
  


  
  // ---------------------------------------------------------------------------
  // Copy constructor
  Point3::Point3( Point3 const& point )
    : Group3( point )
  {}
  
  
  
  
  // ---------------------------------------------------------------------------
  // Copy constructor  
  Point3::Point3( Group3 const& point )
    : Group3( point )
  {}




  // ---------------------------------------------------------------------------
  // Destructor
  Point3::~Point3()
  {}




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




  // ---------------------------------------------------------------------------
  // Distance between 2 points of type Point
  double Point3::DistanceTo( Point3 const& point ) const
  {
    double a = m_comp[0] - point[0];
    double b = m_comp[1] - point[1];
    double c = m_comp[2] - point[2];
    return ( sqrt( a*a + b*b + c*c ) );
  }




  // ---------------------------------------------------------------------------
  // Distance between the point of type Point and a point defined 
  // by a 3-component array
  double Point3::DistanceTo( double const* point) const
  {
    double a = m_comp[0] - point[0];
    double b = m_comp[1] - point[1];
    double c = m_comp[2] - point[2];
    return ( sqrt( a*a + b*b + c*c ) );
  }




  // ---------------------------------------------------------------------------
  // Distance between the point of type Point and a point defined 
  // by 3 scalars
  double Point3::DistanceTo( double x, double y, double z ) const
  {
    double a = m_comp[0] - x;
    double b = m_comp[1] - y;
    double c = m_comp[2] - z;
    return ( sqrt( a*a + b*b + c*c ) );
  }




  // --------------------------------------------------------------------------
  // Equal operator to another Group3 object
  Point3& Point3::operator = ( Point3 const& g2 )
  {
    if ( &g2 != this )
    {      
      m_comp[X] = g2.m_comp[X];
      m_comp[Y] = g2.m_comp[Y];
      m_comp[Z] = g2.m_comp[Z];
    }
    return ( *this );
  }  
}
