#include "Segment.hh"
#include "Transform.hh"
#include "Basic.hh"


// ----------------------------------------------------------------------------
// Constructor with the length as an input parameter
Segment::Segment(double x) 
  : m_halflength( x / 2. )
{}




// ----------------------------------------------------------------------------
// Destructor
Segment::~Segment()
{}




// ----------------------------------------------------------------------------
// Constructor with an input stream
Segment::Segment( istream& fileIn )
{
  readShape( fileIn );
}




// ----------------------------------------------------------------------------
// Returns the convex type
ConvexType Segment::getConvexType() const 
{
  return ( SEGMENT );
}




// ----------------------------------------------------------------------------
// Computes the inertia tensor and the inverse of the inertia
// tensor. Inertia and inverse of inertia are 0 by convention
bool Segment::BuildInertia( double* inertia, double* inertia_1 ) const
{
  inertia[0] = inertia[1] = inertia[2] = inertia[3] = inertia[4]
	= inertia[5] = 0.0;
  inertia_1[0] = inertia_1[1] = inertia_1[2] = inertia_1[3] = inertia_1[4]
	= inertia_1[5] = 0.0;
  return ( true );
}




// ----------------------------------------------------------------------------
// Returns the circumscribed radius of the reference sphere,
// i.e., without applying any transformation
double Segment::computeCircumscribedRadius() const 
{
  return ( m_halflength );
}




// ----------------------------------------------------------------------------
// Returns a clone of the segment
Convex* Segment::clone() const 
{
  return ( new Segment( 2.0 * m_halflength ) );
}




// ----------------------------------------------------------------------------
// Returns the segment length
double Segment::getLength() const 
{
  return ( 2.0 * m_halflength );
}




// ----------------------------------------------------------------------------
// Segment support function, returns the support point P, i.e. the
// point on the surface of the segment that satisfies max(P.v)
Point3 Segment::support( Vector3 const& v ) const
{
  double norm = Norm( v );
  if ( norm > EPSILON ) 
    return ( v[X] > 0.0 ? 
	Point3( m_halflength, 0., 0. ) : Point3( - m_halflength, 0. ,0. ) );
  else
    return ( Point3() );
}




// ----------------------------------------------------------------------------
// Returns the segment volume, here 0 by convention
double Segment::getVolume() const 
{
  return ( 0. );
}




// ----------------------------------------------------------------------------
// Output operator
void Segment::writeShape( ostream& fileOut ) const 
{
  fileOut << "*Segment\n";
  fileOut << 2.0 * m_halflength;
}




// ----------------------------------------------------------------------------
// Input operator
void Segment::readShape( istream& fileIn ) 
{
  fileIn >> m_halflength;
  m_halflength /= 2.0;
}




// ----------------------------------------------------------------------------
// Return the transformation associated to a direction vector and a
// point correponding to the center of mass of the segment
Transform Segment::computeTransform( Vector3 const& v, Point3 const& gc )
{
  Transform tr;

  double normxy = sqrt( v[X] * v[X] + v[Y] * v[Y] ),
	normxz = sqrt( v[X] * v[X] + v[Z] * v[Z] );
  double bx = fabs( v[X] );

  // Angle pr la rotation par rapport à Z
  double angleZ = 0.;
  if ( normxy > 1.e-12 )
  { 
    angleZ = acos( bx / normxy );
    if ( v[X] > 0. )
    {
     if ( v[Y] < 0. ) angleZ *= -1.;
    }
    else
    {
      if ( v[Y] > 0. ) angleZ = PI - angleZ;
      else angleZ += PI;
    }
  }
   
  // Angle pr la rotation par rapport à Y 
  double angleY = 0.;
  if ( normxz > 1.e-12 )
  {
    angleY = acos( bx / normxz );
    if ( v[X] > 0. )
    {
     if ( v[Z] < 0. ) angleY *= -1.;
    }
    else
    {
      if ( v[Z] > 0. ) angleY = PI - angleY;
      else angleY += PI; 
    }
  } 

  // Matrice de rotation
  Matrix rZ( cos(angleZ), - sin(angleZ), 0., 
  	sin(angleZ), cos(angleZ), 0.,
	0., 0., 1. );
  Matrix rY( cos(angleY), 0., sin(angleY),
  	0., 1., 0.,
	- sin(angleY), 0., cos(angleY) );
  Matrix rotation = rY * rZ;
  
  // Set transform
  tr.setBasis( rotation );
  tr.setOrigin( gc );
   
  return ( tr );
}




// ----------------------------------------------------------------------------
// Returns the number of points to write the segment in a Paraview format
int Segment::numberOfPoints_PARAVIEW() const
{
  return ( 2 );
}




// ----------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the segment 
// in a Paraview format
int Segment::numberOfCells_PARAVIEW() const
{
  return ( 1 );  
}




// ----------------------------------------------------------------------------
// Writes a list of points describing the segment in a Paraview format
void Segment::write_polygonsPts_PARAVIEW( ostream& f, 
  	Transform const& transform, Vector3 const* translation ) const
{
  Point3 pp, p;
  
  // Point -
  p[X] = - m_halflength;
  pp = transform( p );
  if ( translation ) pp += *translation;    
  f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
  
  // Point +
  p[X] = m_halflength;
  pp = transform( p );
  if ( translation ) pp += *translation;    
  f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
}




// ----------------------------------------------------------------------------
// Returns a list of points describing the segment in a Paraview format 
list<Point3> Segment::get_polygonsPts_PARAVIEW( Transform const& transform,
  	Vector3 const* translation ) const
{
  list<Point3> ParaviewPoints;
  Point3 pp, p;
  
  // Point -
  p[X] = - m_halflength;
  pp = transform( p );
  if ( translation ) pp += *translation;    
  ParaviewPoints.push_back( pp );
  
  // Point +
  p[X] = m_halflength;
  pp = transform( p );
  if ( translation ) pp += *translation;    
  ParaviewPoints.push_back( pp );
  
  return ( ParaviewPoints ); 
}




// ----------------------------------------------------------------------------
// Writes the segment in a Paraview format
void Segment::write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const
{
  connectivity.push_back( firstpoint_globalnumber );
  connectivity.push_back( firstpoint_globalnumber + 1 ); 
  last_offset += 2;  
  offsets.push_back( last_offset );
  cellstype.push_back( 3 );
  
  firstpoint_globalnumber += 2;
}




// ----------------------------------------------------------------------------
// Returns whether a point lies inside the segment
bool Segment::isIn( Point3 const& pt ) const
{
  return ( pt[X] <= m_halflength && pt[X] >= - m_halflength );
}    
  
