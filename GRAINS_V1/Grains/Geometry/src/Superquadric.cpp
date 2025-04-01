#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <bits/stdc++.h>
#include "Superquadric.hh"
#include "GrainsExec.hh"
#include "sstream"
using namespace std;

int Superquadric::m_visuNodeNbPerQar = 4; 

// ----------------------------------------------------------------------------
// Constructor with input parameters
Superquadric::Superquadric( double a0, double b0, double c0, double n1, 
	double n2 )
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
  const double eps1 = 2.0 / m_n1;
  const double eps2 = 2.0 / m_n2;
  const double C = 0.4 * m_a * m_b * m_c * eps1 * eps2 ;

  const double prod1 = beta( 1.5 * eps2, 0.5 * eps2 ) 
  	* beta( 2. * eps1, 0.5 * eps1 );
  const double prod2 = m_c * m_c * beta( 0.5 * eps2, 0.5 * eps2 ) 
  	* beta( 1.5 * eps1, eps1 );

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
    double abvx = abs( v[X] ), abvy = abs( v[Y] ), abvz = abs( v[Z] ), 
    	signx = copysign( 1.0, v[X] ), signy = copysign( 1.0, v[Y] ), 
	signz = copysign( 1.0, v[Z] );

    Point3 sup;

    if ( abvx == 0. )
    {
      if ( abvy == 0. )
      {
        return ( Point3( 0., 0., signz * m_c ) );
      }

      else
      {
        double alpha = pow( m_c*abvz / (m_b*abvy), 1.0 / (m_n1 - 1.0) );
        double yt = 1.0 / pow( 1.0 + pow( alpha, m_n1 ), 1.0 / m_n1 );

        return ( Point3( 0., signy * m_b * yt, signz * alpha * m_c * yt ) );
      }
    }

    else
    {
      double alpha = pow( m_b * abvy / ( m_a * abvx ), 
      	1.0 / ( m_n2 - 1.0 ) );
      double gamma = pow( 1.0 + pow( alpha, m_n2 ), 
      	( m_n1 - m_n2 ) / ( m_n2 * ( m_n1 - 1.0 ) ) );
      double beta = gamma * pow( m_c * abvz / ( m_a * abvx ), 
      	1.0 / ( m_n1 - 1.0 ) );
      double xt = 1.0 / pow( pow( 1.0 + pow( alpha, m_n2 ), m_n1 / m_n2 )
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
  const double eps1 = 2.0 / m_n1;
  const double eps2 = 2.0 / m_n2;

  return ( C * eps1 * eps2 * beta( eps1, 0.5 * eps1 ) 
  	* beta( 0.5 * eps2, 0.5 * eps2 ) );
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
  return ( 4 * m_visuNodeNbPerQar * ( 2 * m_visuNodeNbPerQar - 1 ) + 3 );
}




// ----------------------------------------------------------------------------
// Writes a list of points describing the superquadric in a Paraview format
void Superquadric::write_polygonsPts_PARAVIEW( ostream &f,
  Transform const& transform, Vector3 const* translation ) const
{
  double eps1 = 2. / m_n1;
  double eps2 = 2. / m_n2;
  double dtheta = PI / ( 2. * m_visuNodeNbPerQar ), dphi = dtheta;
  Point3 pp, pptrans;

  // Regular points on the surface
  double cost, sint, costeps1, sinteps1, cosp, sinp;

  for ( int i = 1; i < 2*m_visuNodeNbPerQar; i++ )
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

    for ( int j = 0; j < 4*m_visuNodeNbPerQar; j++ )
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
  
  pp[X] = 0.;
  pp[Y] = 0.;
  // Bottom point
  pp[Z] = - m_c;
  pptrans = transform( pp );
  if ( translation ) pptrans += *translation;
  f << pptrans[X] << " " << pptrans[Y] << " " << pptrans[Z] << endl;
	
  // Top point
  pp[Z] = m_c;
  pptrans = transform( pp );
  if ( translation ) pptrans += *translation;
  f << pptrans[X] << " " << pptrans[Y] << " " << pptrans[Z] << endl;
	
  // Gravity center
  pp[Z] = 0.;
  pptrans = transform( pp );
  if ( translation ) pptrans += *translation;
  f << pptrans[X] << " " << pptrans[Y] << " " << pptrans[Z] << endl;  
}




// ----------------------------------------------------------------------------
// Returns a list of points describing the superquadric in a Paraview format
list<Point3> Superquadric::get_polygonsPts_PARAVIEW( Transform const& transform,
  Vector3 const* translation ) const
{
  double eps1 = 2. / m_n1 ;
  double eps2 = 2. / m_n2 ;
  double dtheta = PI / ( 2. * m_visuNodeNbPerQar ), dphi = dtheta;
  list<Point3> ParaviewPoints;
  Point3 pp, pptrans;

  // Regular points on the surface
  double cost, sint, costeps1, sinteps1, cosp, sinp;

  for ( int i = 1; i < 2*m_visuNodeNbPerQar; i++ )
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

    for ( int j = 0; j < 4*m_visuNodeNbPerQar; j++ )
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
  
  pp[X] = 0.;
  pp[Y] = 0.;
  // Bottom point
  pp[Z] = - m_c;
  pptrans = transform( pp );
  if ( translation ) pptrans += *translation;
  ParaviewPoints.push_back( pptrans );
	
  // Top point
  pp[Z] = m_c;
  pptrans = transform( pp );
  if ( translation ) pptrans += *translation;
  ParaviewPoints.push_back( pptrans );
	
  // Gravity center
  pp[Z] = 0.;
  pptrans = transform( pp );  
  if ( translation ) pptrans += *translation;
  ParaviewPoints.push_back( pptrans );  

  return ( ParaviewPoints );
}




// ----------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the superquadric in a
// Paraview format
int Superquadric::numberOfCells_PARAVIEW() const
{
  return ( 8 * m_visuNodeNbPerQar * m_visuNodeNbPerQar );
}




// ----------------------------------------------------------------------------
// Writes the superquadric in a Paraview format
void Superquadric::write_polygonsStr_PARAVIEW( list<int>& connectivity,
      list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
  int& last_offset ) const
{
  int i, k, ptsPerlevel = 4 * m_visuNodeNbPerQar,
	Bottom_number = ptsPerlevel * ( 2 * m_visuNodeNbPerQar - 1 ),
	Top_number = ptsPerlevel * ( 2 * m_visuNodeNbPerQar - 1 ) + 1,  
  	GC_number = ptsPerlevel * ( 2 * m_visuNodeNbPerQar - 1 ) + 2;
  
  // Regular cells: Pyramid
  for ( k = 0; k < 2*m_visuNodeNbPerQar-2 ; ++k ) 
  {  
    for ( i = 0; i < ptsPerlevel-1 ; ++i )
    {
      connectivity.push_back( firstpoint_globalnumber + k * ptsPerlevel + i ); 
      connectivity.push_back( firstpoint_globalnumber + k * ptsPerlevel + i 
      	+ 1);
      connectivity.push_back( firstpoint_globalnumber + ( k + 1) * ptsPerlevel
      	+ i + 1 );	
      connectivity.push_back( firstpoint_globalnumber + ( k + 1 ) * ptsPerlevel
	+ i );
      connectivity.push_back( firstpoint_globalnumber + GC_number );
      last_offset += 5;
      offsets.push_back( last_offset );
      cellstype.push_back( 14 );		
    }
    connectivity.push_back( firstpoint_globalnumber + k * ptsPerlevel + 
    	ptsPerlevel - 1 );
    connectivity.push_back( firstpoint_globalnumber + k * ptsPerlevel );
    connectivity.push_back( firstpoint_globalnumber + ( k + 1 ) * ptsPerlevel );
    connectivity.push_back( firstpoint_globalnumber + ( k + 1 ) * ptsPerlevel 
    	+ ptsPerlevel - 1 );
    connectivity.push_back( firstpoint_globalnumber + GC_number );
    last_offset += 5;
    offsets.push_back( last_offset );
    cellstype.push_back( 14 );    
  }  

  // Top cells: tetraedron
  for ( i = 0; i < ptsPerlevel-1 ; ++i )
  {
    connectivity.push_back( firstpoint_globalnumber + i ); 
    connectivity.push_back( firstpoint_globalnumber + i + 1 );
    connectivity.push_back( firstpoint_globalnumber + Top_number );	
    connectivity.push_back( firstpoint_globalnumber + GC_number );
    last_offset += 4;
    offsets.push_back( last_offset );
    cellstype.push_back( 10 );   
  }
  connectivity.push_back( firstpoint_globalnumber + ptsPerlevel - 1 );
  connectivity.push_back( firstpoint_globalnumber );
  connectivity.push_back( firstpoint_globalnumber + Top_number );	
  connectivity.push_back( firstpoint_globalnumber + GC_number );
  last_offset += 4;
  offsets.push_back( last_offset );
  cellstype.push_back( 10 );  
  
  // Bottom cells: tetraedron  
  for ( i = 0; i < ptsPerlevel-1 ; ++i )
  {
    connectivity.push_back( firstpoint_globalnumber
    	+ ( 2 * m_visuNodeNbPerQar - 2 ) * ptsPerlevel + i );
    connectivity.push_back( firstpoint_globalnumber
    	+ ( 2 * m_visuNodeNbPerQar - 2 ) * ptsPerlevel + i + 1 );
    connectivity.push_back( firstpoint_globalnumber + Bottom_number );	
    connectivity.push_back( firstpoint_globalnumber + GC_number );
    last_offset += 4;
    offsets.push_back( last_offset );
    cellstype.push_back( 10 ); 
  }
  connectivity.push_back( firstpoint_globalnumber
  	+ ( 2 * m_visuNodeNbPerQar - 1 ) * ptsPerlevel - 1 );
  connectivity.push_back( firstpoint_globalnumber
  	+ ( 2 * m_visuNodeNbPerQar - 2 ) * ptsPerlevel );
  connectivity.push_back( firstpoint_globalnumber + Bottom_number );	
  connectivity.push_back( firstpoint_globalnumber + GC_number );
  last_offset += 4;
  offsets.push_back( last_offset );
  cellstype.push_back( 10 ); 

  firstpoint_globalnumber += 4 * m_visuNodeNbPerQar * 
  	( 2 * m_visuNodeNbPerQar - 1 ) + 3 ;
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
    // cylinder-like SQ -- in this case, bigger from m_a and m_b is radius
    // and height is 2. * m_c
    if ( m_n2 == 2. )
    {
      Vector3 e( 0., 0., 1. );
      double r = m_a > m_b ? m_a : m_b;
      double h = 2. * m_c;
      bvol = new OBC( r, h, e );
    }
    // cube-like SQ -- in this case, the OBC is the same as OBC for a box
    // if ( m_n1 == m_n2 )
    else
    {
      double xy = fabs( m_a - m_b );
      double xz = fabs( m_a - m_c );
      double yz = fabs( m_b - m_c );
      // pick from xy and xz, store to xy
      int zAxis = xy < xz ? Z : Y;
      xy = xy < xz ? xy : xz;
      // pick from xy and yz
      zAxis = xy < yz ? zAxis : X;
          
      double h;
      if ( zAxis == X ) h = 2. * m_a;
      else if ( zAxis == Y ) h = 2. * m_b;
      else if ( zAxis == Z ) h = 2. * m_c;

      Vector3 e( 0., 0., 0. );
      e[ zAxis ] = 1.;
      double r = sqrt( m_a * m_a + m_b * m_b + m_c * m_c - h * h / 4. );

      bvol = new OBC( r, h, e );
    }
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




// ----------------------------------------------------------------------------
// Sets the number of point per quarter of the equator line for Paraview 
// post-processing, i.e., controls the number of facets in the superquadric
// reconstruction in Paraview
void Superquadric::SetvisuNodeNbPerQar( int nbpts )
{
  m_visuNodeNbPerQar = nbpts;
}




// ----------------------------------------------------------------------------
// Writes the sphere in an OBJ format
void Superquadric::write_convex_OBJ( ostream& f, Transform  const& transform,
    	size_t& firstpoint_number ) const
{
  double eps1 = 2. / m_n1;
  double eps2 = 2. / m_n2;
  double dtheta = PI / ( 2. * m_visuNodeNbPerQar ), dphi = dtheta;
  Point3 p, pp;
  int i, k, ptsPerlevel = 4 * m_visuNodeNbPerQar,
	Bottom_number = ptsPerlevel * ( 2 * m_visuNodeNbPerQar - 1 ),
	Top_number = ptsPerlevel * ( 2 * m_visuNodeNbPerQar - 1 ) + 1;
  double cost, sint, costeps1, sinteps1, cosp, sinp;

  // Vertices
  for (i = 1; i < 2*m_visuNodeNbPerQar; i++)
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

    for ( int j = 0; j < 4*m_visuNodeNbPerQar; j++ )
    {
      cosp = cos( j * dphi );
      sinp = sin( j * dphi );

      if ( cosp == 0. )
        p[X] = 0.;
      else if ( cosp < 0. )
        p[X] = -m_a * sinteps1 * pow( -cosp, eps2 );
      else
        p[X] = m_a * sinteps1 * pow( cosp, eps2 );

      if ( sinp == 0. )
        p[Y] = 0.;
      else if ( sinp < 0. )
        p[Y] = -m_b * sinteps1 * pow( -sinp, eps2 );
      else
        p[Y] = m_b * sinteps1 * pow( sinp, eps2 );

      p[Z] = m_c * costeps1;

      pp = transform( p );
      f << "v " << GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
			pp[X] ) << " " << 
	GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
			pp[Y] ) << " " <<
	GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
			pp[Z] ) << " " << endl;	
    }
  }
  
  p[X] = 0.;
  p[Y] = 0.;
  // Bottom point
  p[Z] = - m_c;
  pp = transform( p );
  f << "v " << GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
			pp[X] ) << " " << 
	GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
			pp[Y] ) << " " <<
	GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
			pp[Z] ) << " " << endl;
	
  // Top point
  p[Z] = m_c;
  pp = transform( p );
  f << "v " << GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
			pp[X] ) << " " << 
	GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
			pp[Y] ) << " " <<
	GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
			pp[Z] ) << " " << endl;

  // Faces  
  // Square faces
  for ( k = 0; k < 2*m_visuNodeNbPerQar-2 ; ++k ) 
  {  
    for ( i = 0; i < ptsPerlevel-1 ; ++i )
    {
      f << "f " << firstpoint_number + k * ptsPerlevel + i << " " 
      	<< firstpoint_number + k * ptsPerlevel + i + 1 << " "
      	<< firstpoint_number + ( k + 1) * ptsPerlevel + i + 1 << " "	
      	<< firstpoint_number + ( k + 1 ) * ptsPerlevel + i << endl;		
    }
    f << "f " << firstpoint_number + k * ptsPerlevel + ptsPerlevel - 1 
    	<< " " << firstpoint_number + k * ptsPerlevel << " " 
	<< firstpoint_number + ( k + 1 ) * ptsPerlevel << " "
	<< firstpoint_number + ( k + 1 ) * ptsPerlevel + ptsPerlevel - 1 
	<< endl;   
  }  

  // Top triangular faces
  for ( i = 0; i < ptsPerlevel-1 ; ++i )
  {
    f << "f " << firstpoint_number + i << " " 
    	<< firstpoint_number + i + 1 << " "
    	<< firstpoint_number + Top_number << endl;
  }
  f << "f " << firstpoint_number + ptsPerlevel - 1 << " "
  	<< firstpoint_number << " "
	<< firstpoint_number + Top_number << endl;
  
  // Bottom triangular faces
  for ( i = 0; i < ptsPerlevel-1 ; ++i )
  {
    f << "f " << firstpoint_number
    		+ ( 2 * m_visuNodeNbPerQar - 2 ) * ptsPerlevel + i << " " 
    	<<  firstpoint_number
    		+ ( 2 * m_visuNodeNbPerQar - 2 ) * ptsPerlevel + i + 1 << " " 
    	<< firstpoint_number + Bottom_number << endl;
  }
  f << "f " << firstpoint_number
  		+ ( 2 * m_visuNodeNbPerQar - 1 ) * ptsPerlevel - 1 << " " 
  	<< firstpoint_number
  		+ ( 2 * m_visuNodeNbPerQar - 2 ) * ptsPerlevel  << " "
	<< firstpoint_number + Bottom_number << endl;

  firstpoint_number += 4 * m_visuNodeNbPerQar * 
  	( 2 * m_visuNodeNbPerQar - 1 ) + 2 ;  
}
