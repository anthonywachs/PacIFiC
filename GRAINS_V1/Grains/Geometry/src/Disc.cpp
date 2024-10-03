#include "Disc.hh"
#include <math.h>

int Disc::m_visuNodeNb = 20;

// ----------------------------------------------------------------------------
// Constructor with radius as input parameter
Disc::Disc( double r ) 
  : m_radius(r)
{}




// ----------------------------------------------------------------------------
// Constructor with an input stream
Disc::Disc( istream& fileIn )
{
  readShape( fileIn );
}




// ----------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
Disc::Disc( DOMNode* root )
{
  m_radius = ReaderXML::getNodeAttr_Double( root, "Radius" );
}




// ----------------------------------------------------------------------------
// Destructor
Disc::~Disc() 
{}




// ----------------------------------------------------------------------
// Returns the convex type
ConvexType Disc::getConvexType() const 
{
  return ( DISC2D );
}




// ----------------------------------------------------------------------------
// Computes the inertia tensor and the inverse of the inertia tensor
bool Disc::BuildInertia( double* inertia, double* inertia_1 ) const
{
  // Active 2D plane is XY -> rotation around Z
  inertia[1] = inertia[2] = inertia[4]= 0.;
  inertia[0] = inertia[3] = PI * pow( m_radius, 4. ) / 4.;
  inertia[5] = PI * pow( m_radius, 4. ) / 2.;
  
  inertia_1[1] = inertia_1[2] = inertia_1[4] = 0.0;
  inertia_1[0] = 1.0 / inertia[0];
  inertia_1[3] = 1.0 / inertia[3];
  inertia_1[5] = 1.0 / inertia[5];

  return ( true );
}




// ----------------------------------------------------------------------------
// Returns the circumscribed radius of the reference disc,
// i.e., without applying any transformation
double Disc::computeCircumscribedRadius() const 
{
  return ( m_radius );
}




// ----------------------------------------------------------------------------
// Returns a clone of the disc
Convex* Disc::clone() const
{
  return ( new Disc( m_radius ) );
}




// ----------------------------------------------------------------------------
// Disc support function, returns the support point P, i.e. the
// point on the surface of the disc that satisfies max(P.v)
Point3 Disc::support( Vector3 const& v ) const 
{
  double s = Norm( v );
  if ( s > EPSILON ) 
  {
    double r = m_radius / s;
    return ( Point3( v[X] * r, v[Y] * r, v[Z] * r ) );
  } 
  else 
    return ( Point3() );
}




// ----------------------------------------------------------------------------
// Returns a vector of points describing the envelope of the disc.
// Here simply returns the point (0,0,0) as a convention
vector<Point3> Disc::getEnvelope() const
{
  vector<Point3> envelope;
  Point3 point( 0.,0.,0. );
  envelope.push_back( point );
  return ( envelope );
}




// ----------------------------------------------------------------------------
// Returns the number of vertices/corners or a code corresponding to
// a specific convex shape. Here 1 as a convention
int Disc::getNbCorners() const
{
  return ( 1 );
}




// ----------------------------------------------------------------------------
// Returns the disc surface area
double Disc::getVolume()const 
{
  return ( PI * m_radius * m_radius );
}




// ----------------------------------------------------------------------------
// Output operator
void Disc::writeShape( ostream& fileOut ) const 
{
  fileOut << "*Disc " << m_radius << " *END";   
}




// ----------------------------------------------------------------------------
// Input operator
void Disc::readShape( istream& fileIn ) 
{
  fileIn >> m_radius;
}




// ----------------------------------------------------------------------------
// Returns the number of points to write the disc in a Paraview format
int Disc::numberOfPoints_PARAVIEW() const 
{ 
  return ( m_visuNodeNb );
}




// ----------------------------------------------------------------------------
// Writes a list of points describing the disc in a Paraview format
void Disc::write_polygonsPts_PARAVIEW( ostream& f, 
  	Transform const& transform, Vector3 const* translation ) const
{
  double angle = 2. * PI / m_visuNodeNb ;
  int i;
  Point3 pp, pptrans;
  
  // Regular points on the surface
  for ( i = 0; i < m_visuNodeNb ; ++i )
  {
    pp[X] = m_radius * cos( i * angle );
    pp[Y] = m_radius * sin( i * angle );
    pptrans = transform( pp );
    if ( translation ) pptrans += *translation;
    f << pptrans[X] << " " << pptrans[Y] << " " << pptrans[Z] << endl;
  }
}




// ----------------------------------------------------------------------------
// Returns a list of points describing the disc in a Paraview format
list<Point3> Disc::get_polygonsPts_PARAVIEW( Transform const& transform,
  	Vector3 const* translation ) const
{
  list<Point3> ParaviewPoints;
  double angle = 2. * PI / m_visuNodeNb ;
  int i;
  Point3 pp, pptrans;
  
  // Regular points on the surface
  for ( i = 0; i < m_visuNodeNb ; ++i )
  {
    pp[X] = m_radius * cos( i * angle );
    pp[Y] = m_radius * sin( i * angle );
    pptrans = transform( pp );
    if ( translation ) pptrans += *translation;
    ParaviewPoints.push_back( pptrans );
  } 
  
  return ( ParaviewPoints );
}




// ----------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the sphere 
// in a Paraview format
int Disc::numberOfCells_PARAVIEW() const
{
  return ( 1 ); 
}




// ----------------------------------------------------------------------------
// Writes the disc in a Paraview format
void Disc::write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const
{
  int ncorners = m_visuNodeNb;

  int count = firstpoint_globalnumber;
  for (int i=0;i<ncorners;++i)
  {
    connectivity.push_back( count );
    ++count;
  }  
  last_offset += ncorners;    
  offsets.push_back( last_offset );
  cellstype.push_back( 7 );
  
  firstpoint_globalnumber += ncorners;
}




// ----------------------------------------------------------------------------
// Sets the number of points over the perimeter of the disc for 
// Paraview post-processing
void Disc::SetvisuNodeNb( int nbpts )
{ 
  m_visuNodeNb = nbpts; 
}




// ----------------------------------------------------------------------------
// Returns whether a point lies inside the disc
bool Disc::isIn( Point3 const& pt ) const
{
  return ( pt[X] * pt[X] + pt[Y] * pt[Y] <= m_radius * m_radius );
}  




// ----------------------------------------------------------------------------
// Returns the bounding volume to disc
BVolume* Disc::computeBVolume( unsigned int type ) const
{
  BVolume* bvol = NULL;
  if ( type == 1 ) // OBB
    bvol = new OBB( Vector3( m_radius, m_radius, 0. ), Matrix() );
  else if ( type == 2 ) // OBC
    bvol = new OBC( m_radius, 0., Vector3( 0., 0., 1. ) );

  return( bvol );
}




// ----------------------------------------------------------------------------
// Performs advanced comparison of the two discs and returns whether 
// they match
bool Disc::equalType_level2( Convex const* other ) const
{
  // We know that other points to a Disc, we dynamically cast it to actual 
  // type
  Disc const* other_ = dynamic_cast<Disc const*>(other);
  
  double lmin = min( m_radius, other_->m_radius );

  bool same = ( fabs( m_radius - other_->m_radius ) <  LOWEPS * lmin );
  
  return ( same );
} 
