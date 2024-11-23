#include "GrainsExec.hh"
#include "GrainsBuilderFactory.hh"
#include "SpheroCylindricalPrism.hh"
#include "BBox.hh"
#include <sstream>

int SpheroCylindricalPrism::m_visuNodeNbPerHalf = 16;


// ----------------------------------------------------------------------------
// Constructor with edge length as input parameters
SpheroCylindricalPrism::SpheroCylindricalPrism( double r, double l, double h )
  : m_radius( r )
  , m_length( l )
  , m_height( h )
{}




// ----------------------------------------------------------------------------
// Constructor with an input stream
SpheroCylindricalPrism::SpheroCylindricalPrism( istream& fileIn )
{
  readShape( fileIn );
}




// ----------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
SpheroCylindricalPrism::SpheroCylindricalPrism( DOMNode* root )
{
  m_radius = ReaderXML::getNodeAttr_Double( root, "Radius" );
  m_length = ReaderXML::getNodeAttr_Double( root, "Length" );
  m_height = ReaderXML::getNodeAttr_Double( root, "Height" );
}




// ----------------------------------------------------------------------------
// Destructor
SpheroCylindricalPrism::~SpheroCylindricalPrism()
{}




// ----------------------------------------------------------------------
// Returns the convex type
ConvexType SpheroCylindricalPrism::getConvexType() const
{
  return ( SPHEROCYLINDRICALPRISM );
}




// ----------------------------------------------------------------------------
// Returns a clone of the spherocylindrical prism
Convex* SpheroCylindricalPrism::clone() const
{
  return ( new SpheroCylindricalPrism( m_radius, m_length, m_height ) );
}




// ----------------------------------------------------------------------
// Computes the inertia tensor and the inverse of the inertia tensor
bool SpheroCylindricalPrism::BuildInertia( double* inertia, double* inertia_1 ) 
	const
{  
  inertia[1] = inertia[2] = inertia[4] = 0.;
  inertia_1[1] = inertia_1[2] = inertia_1[4] = 0.0;  
  
  inertia[0] = 2. * m_height * ( m_length + 2. * m_radius ) * m_radius
    * ( 0.25 * m_height * m_height + m_radius * m_radius ) / 3.0;
  inertia[3] = 2. * m_height * ( m_length + 2. * m_radius ) * m_radius
    * ( ( 0.5 * m_length + m_radius ) * ( 0.5 * m_length + m_radius ) 
    + m_radius * m_radius ) / 3.0;
  inertia[5] = 2. * m_height * ( m_length + 2. * m_radius ) * m_radius
    * ( 0.25 * m_height * m_height 
    + ( 0.5 * m_length + m_radius ) * ( 0.5 * m_length + m_radius ) ) / 3.0;
	
  inertia_1[0] = 1.0 / inertia[0];
  inertia_1[3] = 1.0 / inertia[3];
  inertia_1[5] = 1.0 / inertia[5];
  
  return ( true );
}




// ----------------------------------------------------------------------------
// Returns the circumscribed radius of the reference spherocylindrical prism,
// i.e., without applying any transformation
double SpheroCylindricalPrism::computeCircumscribedRadius() const
{
  return ( sqrt( m_height * m_height 
  	+ ( m_length + m_radius ) * ( m_length + m_radius ) ) );
}




// ----------------------------------------------------------------------------
// Returns a vector of points describing the envelope of the
// cylinder. Here simply returns 3 points as follows: center of bottom circular
// face of the elementary cylinder, an arbitrary point on the lateral surface 
// of the elementary cylinder and center of top circular face of the elementary 
// cylinder
vector<Point3> SpheroCylindricalPrism::getEnvelope() const
{
  Point3 point;
  vector<Point3> envelope( 1, point );

  // TO DO

  return ( envelope );
}




// ----------------------------------------------------------------------------
// Returns a pointer to a 2D array describing the relationship between the face
// indices and the point indices. Returns a null pointer as a convention
vector< vector<int> > const* SpheroCylindricalPrism::getFaces() const
{
  vector< vector<int> >* allFaces = NULL;
  return ( allFaces );
}




// ----------------------------------------------------------------------------
// Returns the number of vertices/corners or a code corresponding to
// a specific convex shape. Here returns the code 4444
int SpheroCylindricalPrism::getNbCorners() const
{
  return ( 4444 );
}




// ----------------------------------------------------------------------------
// Returns the spherocylindrical prism volume
double SpheroCylindricalPrism::getVolume() const
{
  return ( m_height * m_radius * ( PI * m_radius + 2. * m_length ) );
}




// ----------------------------------------------------------------------------
// Output operator
void SpheroCylindricalPrism::writeShape( ostream& fileOut ) const
{
  fileOut << "*SpheroCyl " << m_radius << " " << m_length 
  	<< " " << m_height << " *END";
}




// ----------------------------------------------------------------------------
// Input operator
void SpheroCylindricalPrism::readShape( istream& fileIn )
{
  fileIn >> m_radius >> m_length >> m_height;
}




// ----------------------------------------------------------------------------
// spherocylindrical prism support function, returns the support point P, i.e. 
// the point on the surface of the spherocylindrical prism that satisfies 
// max(P.v)
Point3 SpheroCylindricalPrism::support( Vector3 const& v ) const
{
  double norm = Norm( v );
  if ( norm > EPSILON )
  {
    double s = sqrt( v[X] * v[X] + v[Z] * v[Z] );
    if ( s > EPSILON )
    {
      double d = m_radius / s;
      if ( fabs( v[Y] ) < EPSILON )
      {
        if ( fabs( v[X] ) < EPSILON )
	  return ( Point3( 0., 0., v[Z] * d ) ); 
	else
          return ( Point3( v[X] * d 
	  	+ ( v[X] > 0. ? 0.5 : -0.5 ) * m_length, 0., v[Z] * d ) ); 
      }     
      else 
      {
        if ( fabs( v[X] ) < EPSILON )
	  return ( Point3( 0., 
		( v[Y] > 0. ? 0.5 : -0.5 ) * m_height, v[Z] * d ) );
	else     
          return ( Point3( v[X] * d 
	  	+ ( v[X] > 0. ? 0.5 : -0.5 ) * m_length, 
		( v[Y] > 0. ? 0.5 : -0.5 ) * m_height, v[Z] * d ) );
      }
    }
    else
      return ( Point3( 0., ( v[Y] > 0. ? 0.5 : -0.5 ) * m_height, 0. ) );
  }
  else
    return ( Point3() );    
}




// ----------------------------------------------------------------------------
// Returns the number of points to write the spherocylindrical prism in a  
// Paraview format
int SpheroCylindricalPrism::numberOfPoints_PARAVIEW() const
{
  return ( 4 * ( m_visuNodeNbPerHalf + 2 ) );
}




// ----------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the spherocylindrical 
// prism in a Paraview format
int SpheroCylindricalPrism::numberOfCells_PARAVIEW() const
{
  return ( 2 * m_visuNodeNbPerHalf + 1 );
}




// ----------------------------------------------------------------------------
// Writes a list of points describing the sphere in a Paraview format
void SpheroCylindricalPrism::write_polygonsPts_PARAVIEW( ostream& f,
  	Transform const& transform, Vector3 const* translation ) const
{	 	
  Point3 pp, p;
  double dtheta = PI / m_visuNodeNbPerHalf, tstartleft = 1.5 * PI,
  	tstartright = - 0.5 * PI;

  // Left half cylinder
  // Lower disk rim
  p[Y] = - 0.5 * m_height;
  for (int i=0;i<m_visuNodeNbPerHalf+1;++i)
  {
    p[X] = m_radius * cos ( tstartleft - i * dtheta ) - 0.5 * m_length;
    p[Z] = m_radius * sin ( tstartleft - i * dtheta );
    pp = transform( p );
    if ( translation ) pp += *translation;
    f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
  }

  // Upper disk rim
  p[Y] = 0.5 * m_height;
  for (int i=0;i<m_visuNodeNbPerHalf+1;++i)
  {
    p[X] = m_radius * cos ( tstartleft - i * dtheta ) - 0.5 * m_length;
    p[Z] = m_radius * sin ( tstartleft - i * dtheta );
    pp = transform( p );
    if ( translation ) pp += *translation;
    f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
  }

  // Lower disk center
  p[X] = - 0.5 * m_length;
  p[Y] = - 0.5 * m_height;
  p[Z] = 0.;
  pp = transform( p );
  if ( translation ) pp += *translation;
  f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;

  // Upper disk center
  p[Y] = 0.5 * m_height;
  pp = transform( p );
  if ( translation ) pp += *translation;
  f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
  
  
  // Right half cylinder
  // Lower disk rim
  p[Y] = - 0.5 * m_height;
  for (int i=0;i<m_visuNodeNbPerHalf+1;++i)
  {
    p[X] = m_radius * cos ( tstartright + i * dtheta ) + 0.5 * m_length;
    p[Z] = m_radius * sin ( tstartright + i * dtheta );
    pp = transform( p );
    if ( translation ) pp += *translation;
    f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
  }

  // Upper disk rim
  p[Y] = 0.5 * m_height;
  for (int i=0;i<m_visuNodeNbPerHalf+1;++i)
  {
    p[X] = m_radius * cos ( tstartright + i * dtheta ) + 0.5 * m_length;
    p[Z] = m_radius * sin ( tstartright + i * dtheta );
    pp = transform( p );
    if ( translation ) pp += *translation;
    f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
  }

  // Lower disk center
  p[X] = 0.5 * m_length;
  p[Y] = - 0.5 * m_height;
  p[Z] = 0.;
  pp = transform( p );
  if ( translation ) pp += *translation;
  f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;

  // Upper disk center
  p[Y] = 0.5 * m_height;
  pp = transform( p );
  if ( translation ) pp += *translation;
  f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
  
  
  // Box: no additional point needed  
}




// ----------------------------------------------------------------------------
// Returns a list of points describing the spherocylindrical prism in a Paraview 
// format
list<Point3> SpheroCylindricalPrism::get_polygonsPts_PARAVIEW( 
	Transform const& transform,
	Vector3 const* translation ) const
{
  list<Point3> ParaviewPoints;
  Point3 pp, p;
  double dtheta = PI / m_visuNodeNbPerHalf, tstartleft = 1.5 * PI,
  	tstartright = - 0.5 * PI;

  // Left half cylinder
  // Lower disk rim
  p[Y] = - 0.5 * m_height;
  for (int i=0;i<m_visuNodeNbPerHalf+1;++i)
  {
    p[X] = m_radius * cos ( tstartleft - i * dtheta ) - 0.5 * m_length;
    p[Z] = m_radius * sin ( tstartleft - i * dtheta );
    pp = transform( p );
    if ( translation ) pp += *translation;
    ParaviewPoints.push_back( pp );
  }

  // Upper disk rim
  p[Y] = 0.5 * m_height;
  for (int i=0;i<m_visuNodeNbPerHalf+1;++i)
  {
    p[X] = m_radius * cos ( tstartleft - i * dtheta ) - 0.5 * m_length;
    p[Z] = m_radius * sin ( tstartleft - i * dtheta );
    pp = transform( p );
    if ( translation ) pp += *translation;
    ParaviewPoints.push_back( pp );
  }

  // Lower disk center
  p[X] = - 0.5 * m_length;
  p[Y] = - 0.5 * m_height;
  p[Z] = 0.;
  pp = transform( p );
  if ( translation ) pp += *translation;
  ParaviewPoints.push_back( pp );

  // Upper disk center
  p[Y] = 0.5 * m_height;
  pp = transform( p );
  if ( translation ) pp += *translation;
  ParaviewPoints.push_back( pp );
  
  
  // Right half cylinder
  // Lower disk rim
  p[Y] = - 0.5 * m_height;
  for (int i=0;i<m_visuNodeNbPerHalf+1;++i)
  {
    p[X] = m_radius * cos ( tstartright + i * dtheta ) + 0.5 * m_length;
    p[Z] = m_radius * sin ( tstartright + i * dtheta );
    pp = transform( p );
    if ( translation ) pp += *translation;
    ParaviewPoints.push_back( pp );
  }

  // Upper disk rim
  p[Y] = 0.5 * m_height;
  for (int i=0;i<m_visuNodeNbPerHalf+1;++i)
  {
    p[X] = m_radius * cos ( tstartright + i * dtheta ) + 0.5 * m_length;
    p[Z] = m_radius * sin ( tstartright + i * dtheta );
    pp = transform( p );
    if ( translation ) pp += *translation;
    ParaviewPoints.push_back( pp );
  }

  // Lower disk center
  p[X] = 0.5 * m_length;
  p[Y] = - 0.5 * m_height;
  p[Z] = 0.;
  pp = transform( p );
  if ( translation ) pp += *translation;
  ParaviewPoints.push_back( pp );

  // Upper disk center
  p[Y] = 0.5 * m_height;
  pp = transform( p );
  if ( translation ) pp += *translation;
  ParaviewPoints.push_back( pp );
  
  
  // Box: no additional point needed      
      
  return ( ParaviewPoints );
}




// ----------------------------------------------------------------------------
// Writes the spherocylindrical prism in a Paraview format
void SpheroCylindricalPrism::write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const
{
  // Left half cylinder
  for (int i=0;i<m_visuNodeNbPerHalf;++i)
  {
    connectivity.push_back( firstpoint_globalnumber + i );
    connectivity.push_back( firstpoint_globalnumber + i + 1 );
    connectivity.push_back( firstpoint_globalnumber + 
    	2 * ( m_visuNodeNbPerHalf + 1 ) );
    connectivity.push_back( firstpoint_globalnumber + i + 
    	m_visuNodeNbPerHalf + 1 );
    connectivity.push_back( firstpoint_globalnumber + i + 
    	m_visuNodeNbPerHalf + 2 );
    connectivity.push_back( firstpoint_globalnumber +
    	2 * ( m_visuNodeNbPerHalf + 1 ) + 1 );
    last_offset += 6;
    offsets.push_back( last_offset );
    cellstype.push_back( 13 );
  }
  
  // Right half cylinder 
  int rightstart = firstpoint_globalnumber +  2 * ( m_visuNodeNbPerHalf + 1 ) 
  	+ 2;
  for (int i=0;i<m_visuNodeNbPerHalf;++i)
  {
    connectivity.push_back( rightstart + i );
    connectivity.push_back( rightstart + i + 1 );
    connectivity.push_back( rightstart + 2 * ( m_visuNodeNbPerHalf + 1 ) );
    connectivity.push_back( rightstart + i + m_visuNodeNbPerHalf + 1 );
    connectivity.push_back( rightstart + i + m_visuNodeNbPerHalf + 2 );
    connectivity.push_back( rightstart + 2 * ( m_visuNodeNbPerHalf + 1 ) + 1 );
    last_offset += 6;
    offsets.push_back( last_offset );
    cellstype.push_back( 13 );
  }
  
  // Box
  connectivity.push_back( firstpoint_globalnumber );
  connectivity.push_back( firstpoint_globalnumber + m_visuNodeNbPerHalf );
  connectivity.push_back( firstpoint_globalnumber + 
  	2 * m_visuNodeNbPerHalf + 1 );
  connectivity.push_back( firstpoint_globalnumber + m_visuNodeNbPerHalf + 1 );
  connectivity.push_back( firstpoint_globalnumber + 2 * m_visuNodeNbPerHalf + 
  	4 );
  connectivity.push_back( firstpoint_globalnumber + 3 * m_visuNodeNbPerHalf + 
  	4 );
  connectivity.push_back( firstpoint_globalnumber + 4 * m_visuNodeNbPerHalf + 
  	5 );	
  connectivity.push_back( firstpoint_globalnumber + 3 * m_visuNodeNbPerHalf + 
  	5 );	
  last_offset += 8;
  offsets.push_back( last_offset );
  cellstype.push_back( 12 );
  
  firstpoint_globalnumber += 4 * ( m_visuNodeNbPerHalf + 2 ) ;    	
}




// ----------------------------------------------------------------------------
// Returns whether a point lies inside the spherocylindrical prism
bool SpheroCylindricalPrism::isIn( Point3 const& pt ) const
{
  bool isIn = false;
  double halfHeight = 0.5 * m_height, halfLength = 0.5 * m_length,
  	r2 = m_radius * m_radius;
  
  if ( pt[Y] >= - halfHeight && pt[Y] <= halfHeight )
  { 
    if ( pt[X] >= - halfLength - m_radius && pt[Y] < - halfLength )
      isIn = ( pt[X] + halfLength ) * ( pt[X] + halfLength ) 
	+ pt[Z] * pt[Z] <= r2;
    else if ( pt[X] >= - halfLength && pt[X] <= halfLength )
      isIn = fabs( pt[Z] ) <= m_radius;
    else if ( pt[X] > halfLength && pt[X] <= halfLength + m_radius )
      isIn = ( pt[X] - halfLength ) * ( pt[X] - halfLength ) 
	+ pt[Z] * pt[Z] <= r2;
  }
	
  return ( isIn );
}




// ----------------------------------------------------------------------------
// Returns the bounding volume to spherocylindrical prism
BVolume* SpheroCylindricalPrism::computeBVolume( unsigned int type ) const
{
  BVolume* bvol = NULL;
  if ( type == 1 ) // OBB
  {
    Vector3 const& extent = Vector3( 0.5 * m_length + m_radius, 0.5 * m_height, 
    	m_radius );
    bvol = new OBB( extent, Matrix() );
  }

  else if ( type == 2 ) // OBC
  {
    Vector3 const& e = Vector3( 0., 1., 0. );
    bvol = new OBC( 0.5 * m_length + m_radius, m_height, e );
  }
  
  return( bvol );
}



// ----------------------------------------------------------------------------
// Performs advanced comparison of the two spherocylindrical prisms and returns 
// whether they match
bool SpheroCylindricalPrism::equalType_level2( Convex const* other ) const
{
  // We know that other points to a spherocylindrical prism, we dynamically 
  // cast it to actual type
  SpheroCylindricalPrism const* other_ = 
  	dynamic_cast<SpheroCylindricalPrism const*>(other); 

  double lmin = min( computeCircumscribedRadius(),
  	other_->computeCircumscribedRadius() ); 
  
  bool same = ( 
  	fabs( m_radius - other_->m_radius ) <  LOWEPS * lmin 
	&& fabs( m_length - other_->m_length ) <  LOWEPS * lmin 
	&& fabs( m_height - other_->m_height ) <  LOWEPS * lmin );
  
  return ( same );
} 
