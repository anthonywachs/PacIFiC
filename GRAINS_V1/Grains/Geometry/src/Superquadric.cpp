  #define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <bits/stdc++.h>
#include "Superquadric.hh"
#include "sstream"
#include "BVolume.hh"
#include "OBB.hh"
#include "OBC.hh"

using namespace std;


int Superquadric::visuNodeNbOnPer = 16; // multiple of 4

// ----------------------------------------------------------------------------
// Constructor with input parameters
Superquadric::Superquadric(double a0, double b0, double c0, double n1, double n2)
{
    m_a = a0;
    m_b = b0;
    m_c = c0;
    m_n1 = n1;
    m_n2 = n2;
}




// ----------------------------------------------------------------------------
// Constructor with an input stream
Superquadric::Superquadric( istream& fileIn )
{
  readShape( fileIn );
}




// ----------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
Superquadric::Superquadric( DOMNode* root )
{
  m_a  = ReaderXML::getNodeAttr_Double(root, "a");
  m_b  = ReaderXML::getNodeAttr_Double(root, "b");
  m_c  = ReaderXML::getNodeAttr_Double(root, "c");
  m_n1 = ReaderXML::getNodeAttr_Double(root, "n1");
  m_n2 = ReaderXML::getNodeAttr_Double(root, "n2");
}




// ----------------------------------------------------------------------------
// Destructor
Superquadric::~Superquadric()
{}




// ----------------------------------------------------------------------------
// Returns the convex type
ConvexType Superquadric::getConvexType() const
{
  return ( SUPERQUADRIC );
}




// ----------------------------------------------------------------------------
// Returns whether the convex shape is a superquadric
bool Superquadric::isSuperquadric() const
{
  return ( true );
}




// ----------------------------------------------------------------------------
// Computes the inertia tensor and the inverse of the inertia tensor
bool Superquadric::BuildInertia( double* inertia, double* inertia_1 ) const
{
  const double eps1 = 2.0/m_n1;
  const double eps2 = 2.0/m_n2;
  const double C = 0.4 * m_a * m_b * m_c * eps1 * eps2 ;

  const double prod1 = beta( 1.5*eps2, 0.5*eps2 ) * beta( 2.0*eps1, 0.5*eps1 );
  const double prod2 = m_c*m_c*beta(0.5*eps2, 0.5*eps2) * beta(1.5*eps1, eps1);

  inertia[1] = inertia[2] = inertia[4] = 0.0;
  inertia[0] = C * ( m_b * m_b * prod1 + prod2 );
  inertia[3] = C * ( m_a * m_a * prod1  + prod2 );
  inertia[5] = C * ( m_a * m_a + m_b * m_b ) * prod1;

  inertia_1[1] = inertia_1[2] = inertia_1[4] = 0.0;
  inertia_1[0] = 1.0 / inertia[0];
  inertia_1[3] = 1.0 / inertia[3];
  inertia_1[5] = 1.0 / inertia[5];

  return ( true );
}




// ----------------------------------------------------------------------------
// Returns the circumscribed radius of the reference convex shape,
// i.e., without applying any transformation
double Superquadric::computeCircumscribedRadius() const
{
  double rmax;

  if ( ( m_n1 == 2.0 ) && ( m_n2 == 2.0 ) )
  {
    rmax = max( m_a, max( m_b, m_c ) );
  }

  else if ( m_n1 == 2.0 )
  {
    const double alpha = pow( m_b / m_a, 2.0 / ( m_n2 - 2.0 ) );
    const double xt = 1.0 / pow( 1 + pow( alpha, m_n2 ), 1.0 / m_n2 );

    rmax = max( m_c, sqrt( m_a * xt * m_a * xt
      + alpha * m_b * xt * alpha * m_b * xt ) );
  }

  else if ( m_n2 == 2.0 )
  {
    const double m = max( m_a, m_b );
    const double beta = pow( m_c * m_c / ( m * m ), 1.0 / ( m_n1 - 2.0 ) );
    const double xt = 1.0 / pow( 1.0 + pow( beta, m_n1 ), 1.0 / m_n1 );

    rmax = sqrt( m * xt * m * xt + beta * m_c * xt * beta * m_c * xt ) ;
  }

  else
  {
    const double alpha = pow( m_b / m_a, 2.0 / ( m_n2 - 2.0 ) );
    const double gamma = pow( 1.0 + pow( alpha, m_n2 ), m_n1 / m_n2 - 1.0 );
    const double beta = pow( gamma * m_c*m_c / (m_a*m_a), 1.0 / (m_n1 - 2.0 ) );
    const double xt = 1.0 / pow( pow( 1.0 + pow( alpha, m_n2 ), m_n1 / m_n2 )
      + pow( beta, m_n1 ), 1.0 / m_n1 );

    rmax = sqrt( m_a * xt * m_a * xt + alpha * m_b * xt * alpha * m_b * xt
      + beta * m_c * xt * beta * m_c * xt );
  }

  return ( rmax );
}




// ----------------------------------------------------------------------------
// Returns a clone of the superquadric
Convex* Superquadric::clone() const
{
  return ( new Superquadric( m_a, m_b, m_c, m_n1, m_n2 ) );
}




// ----------------------------------------------------------------------------
// Superquadric support function, returns the support point P, i.e. the
// point on the surface of the superquadric that satisfies max(P.v)
Point3 Superquadric::support( Vector3 const& v ) const
{
  double norm = Norm(v);
  if (norm > EPSILON)
  {
    const double abvx = abs( v[X] );
    const double abvy = abs( v[Y] );
    const double abvz = abs( v[Z] );
    const double signx = copysign( 1.0, v[X] );
    const double signy = copysign( 1.0, v[Y] );
    const double signz = copysign( 1.0, v[Z] );

    Point3 sup;

    if ( abvx == 0. )
    {
      if ( abvy == 0. )
      {
        return ( Point3( 0., 0., signz * m_c ) );
      }

      else
      {
        const double alpha = pow( m_c*abvz / (m_b*abvy), 1.0 / (m_n1 - 1.0) );
        const double yt = 1.0 / pow( 1.0 + pow( alpha, m_n1 ), 1.0 / m_n1 );

        return ( Point3( 0., signy * m_b * yt, signz * alpha * m_c * yt ) );
      }
    }

    else
    {
      const double alpha = pow( m_b*abvy / (m_a*abvx), 1.0 / ( m_n2 - 1.0 ) );
      // const double gamma = pow( 1.0 + pow( alpha, m_n2 ), m_n1 / m_n2 - 1.0 );
      // const double beta = pow( gamma*m_c*abvz / (m_a*abvx), 1.0 / (m_n1-1.0) );
      const double gamma = pow( 1.0 + pow( alpha, m_n2 ), ( m_n1 - m_n2 ) / ( m_n2 * ( m_n1 - 1.0 ) ) );
      const double beta = gamma * pow( m_c*abvz / (m_a*abvx), 1.0 / (m_n1-1.0) );
      const double xt = 1.0 / pow( pow( 1.0+pow(alpha,m_n2) , m_n1 / m_n2 )
        + pow( beta , m_n1 ) , 1.0 / m_n1 ) ;

      return ( Point3( signx * m_a * xt, signy * alpha * m_b * xt, signz * beta
         * m_c * xt) );
    }
  }

  else
  {
    return ( Point3() );
  }
}




// ----------------------------------------------------------------------------
// Returns a vector of points describing the envelope of the sphere
vector<Point3> Superquadric::getEnvelope() const
{
  Point3 point( 0., 0., 0. );
  vector<Point3> envelope( 3, point );
  /**  envelope[0][Y] = - halfHeight;
  envelope[1][Y] = - halfHeight;
  envelope[1][X] = radius;
  envelope[2][Y] = halfHeight; */
  return ( envelope );
}




// ----------------------------------------------------------------------------
// Returns the number of vertices/corners
int Superquadric::getNbCorners() const
{
  return ( 222 );
}




// ----------------------------------------------------------------------------
// Returns a pointer to a 2D array describing the relationship between the face
// indices and the point indices. Returns a null pointer as a convention
vector< vector<int> > const* Superquadric::getFaces() const
{
  vector< vector<int> >* allFaces = NULL;
  return ( allFaces );
}




// ----------------------------------------------------------------------------
// Returns the superquadric volume
double Superquadric::getVolume() const
{
  const double C = 2.0 / 3.0 * m_a * m_b * m_c;
  const double eps1 = 2.0/m_n1;
  const double eps2 = 2.0/m_n2;

  return ( C*eps1*eps2 * beta( eps1, 0.5*eps1 ) * beta( 0.5*eps2, 0.5*eps2 ) );
}




// ----------------------------------------------------------------------------
// Output operator
void Superquadric::writeShape( ostream &fileOut ) const
{
  fileOut << "*Superquadric " << m_a << " " << m_b << " " << m_c << " "
          << m_n1 << " " << m_n2 << " *END";
}




// ----------------------------------------------------------------------------
// Input operator
void Superquadric::readShape( istream &fileIn )
{
  fileIn >> m_a >> m_b >> m_c >> m_n1 >> m_n2;
}




// ----------------------------------------------------------------------------
// Returns the number of points to write the superquadric in a Paraview format
int Superquadric::numberOfPoints_PARAVIEW() const
{
  return ( visuNodeNbOnPer * ( visuNodeNbOnPer - 1 ) + 3 );
}




// ----------------------------------------------------------------------------
// Writes a list of points describing the superquadric in a Paraview format
void Superquadric::write_polygonsPts_PARAVIEW( ostream &f,
  Transform const& transform, Vector3 const* translation ) const
{
  const double eps1 = 2.0 / m_n1;
  const double eps2 = 2.0 / m_n2;
  const double dtheta = PI / visuNodeNbOnPer;
  const double dphi = 2.0 * PI / visuNodeNbOnPer;
  Point3 pp, pptrans;

  pp[X] = pp[Y] = 0.;
  // Center
  pp[Z] = 0.;
  pptrans = transform( pp );
  if ( translation ) pptrans += *translation;
  f << pptrans[X] << " " << pptrans[Y] << " " << pptrans[Z] << endl;

  // Top point
  pp[Z] = m_c;
  pptrans = transform( pp );
  if ( translation ) pptrans += *translation;
  f << pptrans[X] << " " << pptrans[Y] << " " << pptrans[Z] << endl;

  // Bottom point
  pp[Z] = -m_c;
  pptrans = transform( pp ) ;
  if ( translation ) pptrans += *translation;
  f << pptrans[X] << " " << pptrans[Y] << " " << pptrans[Z] << endl;

  // Regular points on the surface
  double cost, sint, costeps1, sinteps1, cosp, sinp;

  for ( int i = 1; i < visuNodeNbOnPer; i++ )
  {
    cost = cos( i * dtheta );
    sint = sin( i * dtheta );

    if ( cost == 0. )
      costeps1 = 0.;
    else if ( cost < 0. )
      costeps1 = -pow( -cost, eps1 );
    else
      costeps1 = pow( cost, eps1 );

    // Theta is always strictly between 0 and pi so sint is strictly positive
    sinteps1 = pow( sint, eps1 );

    for ( int j = 0; j < visuNodeNbOnPer; j++ )
    {
      cosp = cos( j * dphi );
      sinp = sin( j * dphi );

      if ( cosp == 0. )
        pp[X] = 0.;
      else if ( cosp < 0. )
        pp[X] = -m_a * sinteps1 * pow( -cosp, eps2 );
      else
        pp[X] = m_a * sinteps1 * pow( cosp, eps2 );

      if ( sinp == 0. )
        pp[Y] = 0.;
      else if ( sinp < 0. )
        pp[Y] = -m_b * sinteps1 * pow( -sinp, eps2 );
      else
        pp[Y] = m_b * sinteps1 * pow( sinp, eps2 );

      pp[Z] = m_c * costeps1;

      pptrans = transform( pp );
      if ( translation ) pptrans += *translation;
      f << pptrans[X] << " " << pptrans[Y] << " " << pptrans[Z] << endl;
    }
  }
}




// ----------------------------------------------------------------------------
// Returns a list of points describing the superquadric in a Paraview format
list<Point3> Superquadric::get_polygonsPts_PARAVIEW( Transform const& transform,
  Vector3 const* translation ) const
  {
    const double eps1 = 2.0/m_n1 ;
    const double eps2 = 2.0/m_n2 ;
    list<Point3> ParaviewPoints;
    const double dtheta = PI / visuNodeNbOnPer;
    const double dphi = 2.0 * PI / visuNodeNbOnPer;
    Point3 pp, pptrans;

    pp[X] = pp[Y] = 0.;
    // Gravity center
    pp[Z] = 0.;
    pptrans = transform( pp );
    if ( translation ) pptrans += *translation;
    ParaviewPoints.push_back( pptrans );

    // Top point
    pp[Z] = m_c;
    pptrans = transform( pp );
    if ( translation ) pptrans += *translation;
    ParaviewPoints.push_back( pptrans );

    // Bottom point
    pp[Z] = -m_c;
    pptrans = transform( pp ) ;
    if ( translation ) pptrans += *translation;
    ParaviewPoints.push_back( pptrans );

    // Regular points on the surface
    double cost, sint, costeps1, sinteps1, cosp, sinp;

    for ( int i = 1; i < visuNodeNbOnPer; i++ )
    {
      cost = cos( i * dtheta );
      sint = sin( i * dtheta );

      if ( cost == 0. )
        costeps1 = 0.;
      else if ( cost < 0. )
        costeps1 = -pow( -cost, eps1 );
      else
        costeps1 = pow( cost, eps1 );

      // Theta is always strictly between 0 and pi so sint is
      // strictly positive
      sinteps1 = pow( sint, eps1 );

      for ( int j = 0; j < visuNodeNbOnPer; j++ )
      {
        cosp = cos( j * dphi );
        sinp = sin( j * dphi );
        if ( cosp == 0. )
          pp[X] = 0.;
        else if ( cosp < 0. )
          pp[X] = -m_a * sinteps1 * pow( -cosp, eps2 );
        else
          pp[X] = m_a * sinteps1 * pow( cosp, eps2 );

        if ( sinp == 0. )
          pp[Y] = 0.;
        else if ( sinp < 0. )
          pp[Y] = -m_b * sinteps1 * pow( -sinp, eps2 );
        else
          pp[Y] = m_b * sinteps1 * pow( sinp, eps2 );

        pp[Z] = m_c * costeps1 ;

        pptrans = transform( pp ) ;
        if ( translation ) pptrans += *translation;
        ParaviewPoints.push_back( pptrans );
      }
    }

    return ( ParaviewPoints );
}




// ----------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the superquadric in a
// Paraview format
int Superquadric::numberOfCells_PARAVIEW() const
{
  return ( visuNodeNbOnPer * visuNodeNbOnPer );
}




// ----------------------------------------------------------------------------
// Writes the superquadric in a Paraview format
void Superquadric::write_polygonsStr_PARAVIEW( list<int>& connectivity,
      list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
  int& last_offset ) const
{
  // Top cells: tetrahedron
  for ( int j = 0; j < visuNodeNbOnPer - 1; j++ )
  {
    // Center
    connectivity.push_back( firstpoint_globalnumber );
    // Top point
    connectivity.push_back( firstpoint_globalnumber + 1 );
    connectivity.push_back( firstpoint_globalnumber + 3 + j );
    connectivity.push_back( firstpoint_globalnumber + 3 + j + 1);
    last_offset += 4;
    offsets.push_back( last_offset );
    cellstype.push_back( 10 );
  }
  // Center
  connectivity.push_back( firstpoint_globalnumber );
  // Top point
  connectivity.push_back( firstpoint_globalnumber + 1 );
  // Last point of the first latitude (i=1)
  connectivity.push_back( firstpoint_globalnumber + 3 + visuNodeNbOnPer - 1 );
  // First point of the first latitude (i=1)
  connectivity.push_back( firstpoint_globalnumber + 3 );
  last_offset += 4;
  offsets.push_back( last_offset );
  cellstype.push_back( 10 );

  // Regular cells: Pyramid
  int pointact = 3; // Current point (we will browse the grid)

  for ( int i = 1 ; i < visuNodeNbOnPer - 1; ++i )
  {
    for ( int j = 0; j < visuNodeNbOnPer - 1; j++ )
    {
      // Center
      connectivity.push_back( firstpoint_globalnumber );
      // Current point
      connectivity.push_back( firstpoint_globalnumber + pointact );
      // same latitude, +1 longitude
      connectivity.push_back( firstpoint_globalnumber + pointact + 1) ;
      // +1 latitude, +1 longitude
      connectivity.push_back( firstpoint_globalnumber + pointact
        + visuNodeNbOnPer + 1 );
      // +1 latitude, same longitude
      connectivity.push_back( firstpoint_globalnumber + pointact
        + visuNodeNbOnPer );
      last_offset += 5;
      offsets.push_back( last_offset );
      cellstype.push_back( 14 );
      pointact++;
    }
    // Center
    connectivity.push_back( firstpoint_globalnumber );
    // Current point (last of its latitude)
    connectivity.push_back( firstpoint_globalnumber + pointact );
    // First point (j=0) on the same latitude as pointact
    connectivity.push_back( firstpoint_globalnumber + pointact
      - ( visuNodeNbOnPer - 1 ) );
    // First point (j=0) on the latitude under that of pointact
    connectivity.push_back( firstpoint_globalnumber + pointact + 1 );
    // +1 latitude, same longitude
    connectivity.push_back( firstpoint_globalnumber + pointact
      + visuNodeNbOnPer );
    last_offset += 5;
    offsets.push_back( last_offset );
    cellstype.push_back( 14 );
    pointact++;
  }

  // Bottom cells: tetrahedron
  const int Firstptlastlat = 3 + visuNodeNbOnPer * (visuNodeNbOnPer - 2);

  for ( int j = 0; j < visuNodeNbOnPer - 1; j++ )
  {
    connectivity.push_back( firstpoint_globalnumber ); // Center
    connectivity.push_back( firstpoint_globalnumber + 2 ); // Bottom point
    connectivity.push_back( firstpoint_globalnumber + Firstptlastlat + j );
    connectivity.push_back( firstpoint_globalnumber + Firstptlastlat + j + 1 );
    last_offset += 4;
    offsets.push_back( last_offset );
    cellstype.push_back( 10 );
  }
  connectivity.push_back( firstpoint_globalnumber ); // Center
  connectivity.push_back( firstpoint_globalnumber + 2 ); // Bottom point
  // Last point last latitude
  connectivity.push_back( firstpoint_globalnumber + Firstptlastlat
    + visuNodeNbOnPer - 1 );
  // First point last latitude
  connectivity.push_back( firstpoint_globalnumber + Firstptlastlat );
  last_offset += 4;
  offsets.push_back( last_offset );
  cellstype.push_back( 10 );

  firstpoint_globalnumber += visuNodeNbOnPer * ( visuNodeNbOnPer - 1 ) + 3;
}




// ----------------------------------------------------------------------------
// Returns whether a point lies inside the superquadric
bool Superquadric::isIn( Point3 const& pt ) const
{
  return ( pow( pow( pt[X]/m_a, m_n1 ) + pow( pt[Y]/m_b, m_n1 ), m_n1/m_n2 ) +
    pow( pt[Z]/m_c, m_n2 ) <= 1. );
}




// ----------------------------------------------------------------------------
// Returns the bounding volume to superquadrics
BVolume* Superquadric::computeBVolume( unsigned int type ) const
{
  BVolume* bvol = NULL;
  if ( type == 1 ) // OBB
    bvol = new OBB( Vector3( m_a, m_b, m_c ), Matrix() );
  else if ( type == 2 ) // OBC
  {
    double a[2];
    int axis = ( a[X] = fabs( m_a ) ) < ( a[Y] = fabs( m_b ) ) ? Y : X;
    axis = a[axis] < fabs( m_c ) ? Z : axis;
    Vector3 e( 0., 0., 0. );
    e[axis] = 1.;
    double h = 0., r = 0.;
    switch ( axis )
    {
      case 0:
        h = 2. * m_a;
        r = m_b > m_c ? m_b : m_c;
      break;
      case 1:
        h = 2. * m_b;
        r = m_a > m_c ? m_a : m_c;
      break;
      case 2:
        h = 2. * m_c;
        r = m_b > m_a ? m_b : m_a;
      break;
    }

    bvol = new OBC( r, h, e );
  }

  return( bvol );
}




// ----------------------------------------------------------------------------
// Performs advanced comparison of the two superquadrics and returns whether 
// they match
bool Superquadric::equalType_level2( Convex const* other ) const
{
  // We know that other points to a Superquadric, we dynamically cast it to 
  // actual type
  Superquadric const* other_ = dynamic_cast<Superquadric const*>(other);
  
  double lmin = min( computeCircumscribedRadius(),
  	other_->computeCircumscribedRadius() );    

  bool same = ( 
  	fabs( m_a - other_->m_a ) <  LOWEPS * lmin 
	&& fabs( m_b - other_->m_b ) <  LOWEPS * lmin
	&& fabs( m_c - other_->m_c ) <  LOWEPS * lmin	
	&& fabs( m_n1 - other_->m_n1 ) <  LOWEPS 
	&& fabs( m_n2 - other_->m_n2 ) <  LOWEPS );
  
  return ( same );
}
