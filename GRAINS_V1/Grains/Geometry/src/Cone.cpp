#include "Cone.hh"

// ----------------------------------------------------------------------
// Constructor with flat base radius and height as input parameters
Cone::Cone( double r, double h ) 
  : m_bottomRadius( r )
  , m_quarterHeight( h / 4. )
  , m_sinAngle( r / sqrt( r * r + h * h ) )  
{} 




// ----------------------------------------------------------------------
// Constructor with an input stream
Cone::Cone( istream& fileIn )
{
  readShape( fileIn );
}




// ----------------------------------------------------------------------
// Destructor
Cone::~Cone() 
{}




// ----------------------------------------------------------------------
// Returns the convex type
ConvexType Cone::getConvexType() const 
{
  return ( CONE );
}




// ----------------------------------------------------------------------
// Cone support function, returns the support point P, i.e. the
// point on the surface of the sphere that satisfies max(P.v)
Point3 Cone::support( Vector3 const& v ) const 
{
  double norm = Norm( v );
  if ( norm > EPSILON ) 
  {
    if ( v[Y] > norm * m_sinAngle ) 
      return ( Point3( 0., 3.0 * m_quarterHeight, 0. ) );
    else 
    {
      double s = sqrt( v[X] * v[X] + v[Z] * v[Z] );
      if ( s > EPSILON ) 
      {
	double d = m_bottomRadius / s;  
	return ( Point3( v[X] * d, - m_quarterHeight, v[Z] * d ) );
      } 
      else
	return ( Point3( 0, - m_quarterHeight, 0 ) );
    }
  } 
  else
    return ( Point3() );
}




// ----------------------------------------------------------------------
// Returns a clone of the cone
Convex* Cone::clone() const 
{
  return ( new Cone( m_bottomRadius, 4.0 * m_quarterHeight ) );
}




// ----------------------------------------------------------------------
// Returns the cone volume
double Cone::getVolume() const
{
  return ( 4.0 * m_quarterHeight * PI * m_bottomRadius * m_bottomRadius / 3.0 );
}




// ----------------------------------------------------------------------
// Computes the inertia tensor and the inverse of the inertia tensor
bool Cone::BuildInertia( double* inertia, double* inertia_1 ) const 
{
  inertia[1] = inertia[2] = inertia[4] = 0.0;
  const double constant = 0.2 * m_quarterHeight * m_bottomRadius
	* m_bottomRadius * PI;
  inertia[0] = inertia[5] = constant *
	( 4. * m_quarterHeight * m_quarterHeight 
	+ m_bottomRadius * m_bottomRadius );
  inertia[3] = 2. * constant * m_bottomRadius * m_bottomRadius;

  inertia_1[1] = inertia_1[2] = inertia_1[4] = 0.0;
  inertia_1[5] = inertia_1[0] = 1.0/inertia[0];
  inertia_1[3] = 1.0/inertia[3];
  return ( true );
}




// ----------------------------------------------------------------------------
// Returns the circumscribed radius of the reference box,
// i.e., without applying any transformation
double Cone::computeCircumscribedRadius() const 
{
  return ( max( sqrt( m_bottomRadius * m_bottomRadius 
  	+ m_quarterHeight * m_quarterHeight ), 3. * m_quarterHeight ) );
}




// ----------------------------------------------------------------------------
// Output operator
void Cone::writeShape( ostream& fileOut ) const 
{
  fileOut << "*Cone\n";
  fileOut << m_bottomRadius     << '\t' << 4.0 * m_quarterHeight << '\n' ;
}




// ----------------------------------------------------------------------
// Input operator
void Cone::readShape( istream& fileIn ) 
{
  fileIn >> m_bottomRadius >> m_quarterHeight;
  m_quarterHeight /= 4.0;
  m_sinAngle = m_bottomRadius / sqrt( m_bottomRadius * m_bottomRadius + 
	m_quarterHeight * m_quarterHeight );
}




// ----------------------------------------------------------------------------
// Returns whether a point lies inside the cone
bool Cone::isIn( Point3 const& pt ) const
{
  // TO DO
  
  return ( false );
}  
