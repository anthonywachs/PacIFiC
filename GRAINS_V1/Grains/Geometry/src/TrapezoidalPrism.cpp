#include "GrainsExec.hh"
#include "GrainsBuilderFactory.hh"
#include "TrapezoidalPrism.hh"
#include "BBox.hh"
#include <sstream>

// Static attribute
vector< vector<int> > TrapezoidalPrism::m_allFaces;



// ----------------------------------------------------------------------------
// Constructor with edge length as input parameters
TrapezoidalPrism::TrapezoidalPrism( double width, double a, double b, 
	double h )
{
  m_extent = new double[5];
  m_extent[0] = 0.5 * width;
  m_extent[1] = 0.5 * a;  
  m_extent[2] = 0.5 * b;  
  m_extent[3] = ( 2. * b + a ) * h / ( 3. * ( a + b ) );
  m_extent[4] = ( b + 2. * a ) * h / ( 3. * ( a + b ) ); 
  setCornersFaces();
  double l = 0.5 * ( b - a );
  m_sinAngle = l / sqrt( l * l + h * h );
}




// ----------------------------------------------------------------------------
// Constructor with an input stream
TrapezoidalPrism::TrapezoidalPrism( istream& fileIn )
{
  m_extent = new double[5];
  readShape( fileIn );
  setCornersFaces();
  double l = m_extent[2] - m_extent[1], h = m_extent[3] + m_extent[4];  
  m_sinAngle = l / sqrt( l * l + h * h );  
}




// ----------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
TrapezoidalPrism::TrapezoidalPrism( DOMNode* root )
{
  double width = ReaderXML::getNodeAttr_Double( root, "W" );
  double a = ReaderXML::getNodeAttr_Double( root, "A" );
  double b = ReaderXML::getNodeAttr_Double( root, "B" );
  if ( a > b )
    cout << "Warning: In TrapezoidalPrism, B must be larger than A" << endl; 
  double h = ReaderXML::getNodeAttr_Double( root, "H" );
  m_extent = new double[5];
  m_extent[0] = 0.5 * width;
  m_extent[1] = 0.5 * a;  
  m_extent[2] = 0.5 * b;  
  m_extent[3] = ( 2. * b + a ) * h / ( 3. * ( a + b ) );
  m_extent[4] = ( b + 2. * a ) * h / ( 3. * ( a + b ) );    
  setCornersFaces(); 
  double l = 0.5 * ( b - a );
  m_sinAngle = l / sqrt( l * l + h * h );
}




// ----------------------------------------------------------------------------
// Destructor
TrapezoidalPrism::~TrapezoidalPrism()
{
  m_corners.clear();
  delete [] m_extent;
}




// ----------------------------------------------------------------------
// Returns the convex type
ConvexType TrapezoidalPrism::getConvexType() const
{
  return ( TRAPEZOIDALPRISM );
}




// ----------------------------------------------------------------------------
// Sets the corner/vertex coordinates and the face/vertices numbering
void TrapezoidalPrism::setCornersFaces()
{
  m_corners.reserve( 8 );
  Point3 sommet;
  sommet.setValue( - m_extent[0], - m_extent[1], - m_extent[3] );
  m_corners.push_back(sommet);
  sommet.setValue( m_extent[0], - m_extent[1], - m_extent[3] );
  m_corners.push_back(sommet);
  sommet.setValue( m_extent[0], - m_extent[2], m_extent[4] );
  m_corners.push_back(sommet);
  sommet.setValue( - m_extent[0], - m_extent[2], m_extent[4] );
  m_corners.push_back(sommet);
  sommet.setValue( - m_extent[0], m_extent[1], - m_extent[3] );
  m_corners.push_back(sommet);
  sommet.setValue( m_extent[0], m_extent[1], - m_extent[3] );
  m_corners.push_back(sommet);
  sommet.setValue( m_extent[0], m_extent[2], m_extent[4] );
  m_corners.push_back(sommet);
  sommet.setValue( - m_extent[0], m_extent[2], m_extent[4] );
  m_corners.push_back(sommet);

  if ( m_allFaces.empty() )
  {
    vector<int> oneFace( 4, 0 );
    m_allFaces.reserve( 6 );
    for (int i=0;i<6;++i) m_allFaces.push_back( oneFace );

    m_allFaces[0][0]=4;
    m_allFaces[0][1]=5;
    m_allFaces[0][2]=6;
    m_allFaces[0][3]=7;

    m_allFaces[1][0]=5;
    m_allFaces[1][1]=1;
    m_allFaces[1][2]=2;
    m_allFaces[1][3]=6;

    m_allFaces[2][0]=1;
    m_allFaces[2][1]=0;
    m_allFaces[2][2]=3;
    m_allFaces[2][3]=2;

    m_allFaces[3][0]=0;
    m_allFaces[3][1]=4;
    m_allFaces[3][2]=7;
    m_allFaces[3][3]=3;

    m_allFaces[4][0]=4;
    m_allFaces[4][1]=0;
    m_allFaces[4][2]=1;
    m_allFaces[4][3]=5;

    m_allFaces[5][0]=3;
    m_allFaces[5][1]=7;
    m_allFaces[5][2]=6;
    m_allFaces[5][3]=2;
  }
}




// ----------------------------------------------------------------------------
// Returns a clone of the trapezoidal prism
Convex* TrapezoidalPrism::clone() const
{
  return ( new TrapezoidalPrism( 2. * m_extent[0], 2. * m_extent[1], 
  	2. * m_extent[2], m_extent[3] + m_extent[4] ) );
}




// ----------------------------------------------------------------------
// Computes the inertia tensor and the inverse of the inertia tensor
bool TrapezoidalPrism::BuildInertia( double* inertia, double* inertia_1 ) const
{  
  double vol = 0.;
  inertia[0] = inertia[1] = inertia[2] = inertia[3] = inertia[4] = 
  	inertia[5] = 0;

  GrainsExec::computeVolumeInertiaContrib( m_corners[0], 
      	m_corners[1], m_corners[2], vol, inertia ); 
  GrainsExec::computeVolumeInertiaContrib( m_corners[0], 
      	m_corners[2], m_corners[3], vol, inertia ); 
	
  GrainsExec::computeVolumeInertiaContrib( m_corners[4], 
      	m_corners[5], m_corners[6], vol, inertia ); 
  GrainsExec::computeVolumeInertiaContrib( m_corners[4], 
      	m_corners[6], m_corners[7], vol, inertia );
	
  GrainsExec::computeVolumeInertiaContrib( m_corners[0], 
      	m_corners[1], m_corners[5], vol, inertia ); 
  GrainsExec::computeVolumeInertiaContrib( m_corners[0], 
      	m_corners[5], m_corners[4], vol, inertia );		

  GrainsExec::computeVolumeInertiaContrib( m_corners[3], 
      	m_corners[7], m_corners[6], vol, inertia ); 
  GrainsExec::computeVolumeInertiaContrib( m_corners[3], 
      	m_corners[6], m_corners[2], vol, inertia );

  GrainsExec::computeVolumeInertiaContrib( m_corners[0], 
      	m_corners[4], m_corners[7], vol, inertia ); 
  GrainsExec::computeVolumeInertiaContrib( m_corners[0], 
      	m_corners[7], m_corners[3], vol, inertia );

  GrainsExec::computeVolumeInertiaContrib( m_corners[1], 
      	m_corners[2], m_corners[6], vol, inertia ); 
  GrainsExec::computeVolumeInertiaContrib( m_corners[1], 
      	m_corners[6], m_corners[5], vol, inertia );  
	
  inertia[1] = inertia[2] = inertia[4] = 0.0;
  inertia_1[1] = inertia_1[2] = inertia_1[4] = 0.0;
  inertia_1[0] = 1.0 / inertia[0];
  inertia_1[3] = 1.0 / inertia[3];
  inertia_1[5] = 1.0 / inertia[5];  

  return ( true );
}




// ----------------------------------------------------------------------------
// Returns the edge half-lengths of the trapezoidal prism
double const* TrapezoidalPrism::getExtent() const
{
  return ( m_extent );
}




// ----------------------------------------------------------------------------
// Returns the circumscribed radius of the reference trapezoidal prism,
// i.e., without applying any transformation
double TrapezoidalPrism::computeCircumscribedRadius() const
{
  return ( max( Norm( m_corners[2] ), Norm( m_corners[1] ) ) );
}




// ----------------------------------------------------------------------------
// Returns a vector of points describing the envelope of the trapezoidal prism
vector<Point3> TrapezoidalPrism::getEnvelope() const
{
  return ( m_corners );
}




// ----------------------------------------------------------------------------
// Returns a pointer to a 2D array describing the relationship
// between the face indices and the point indices
vector<vector<int> > const* TrapezoidalPrism::getFaces() const
{
  return ( &m_allFaces );
}




// ----------------------------------------------------------------------------
// Returns the number of vertices/corners or a code corresponding to
// a specific convex shape. Here returns the code 1111
int TrapezoidalPrism::getNbCorners() const
{
  return ( 1111 );
}




// ----------------------------------------------------------------------------
// Returns the trapezoidal prism volume
double TrapezoidalPrism::getVolume() const
{
  return ( 2. * m_extent[0] * ( m_extent[1] + m_extent[2] ) * 
  	( m_extent[3] + m_extent[4] ) );
}




// ----------------------------------------------------------------------------
// Output operator
void TrapezoidalPrism::writeShape( ostream& fileOut ) const
{
  fileOut << "*TrapezoidalPrism " << m_extent[0] << " " << m_extent[1] 
  	<< " " << m_extent[2] << " " << m_extent[3] << " " << m_extent[4] 
	<< " *END";
}




// ----------------------------------------------------------------------------
// Input operator
void TrapezoidalPrism::readShape( istream& fileIn )
{
  fileIn >> m_extent[0] >> m_extent[1] >> m_extent[2] >> m_extent[3] 
  	>> m_extent[4];
}




// ----------------------------------------------------------------------------
// Trapezoidal prism support function, returns the support point P, i.e. the
// point on the surface of the trapezoidal prism that satisfies max(P.v)
Point3 TrapezoidalPrism::support( Vector3 const& v ) const
{
  double norm = Norm( v );
  if ( norm > EPSILON ) 
  {
    double s = sqrt( v[Y] * v[Y] + v[Z] * v[Z] );
    return ( Point3( v[X] < 0. ? - m_extent[0] : m_extent[0],
	v[Y] < 0. ? 
	( v[Z] < - s * m_sinAngle ? - m_extent[1] : - m_extent[2] ) : 
	( v[Z] < - s * m_sinAngle ? m_extent[1] : m_extent[2] ),
	v[Z] < - s * m_sinAngle ? - m_extent[3] : m_extent[4] ) );    
  }   
  else
    return ( Point3() );
}




// ----------------------------------------------------------------------------
// Returns the number of points to write the trapezoidal prism in a Paraview 
// format
int TrapezoidalPrism::numberOfPoints_PARAVIEW() const
{
  return ( 8 );
}




// ----------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the trapezoidal prism
// in a Paraview format
int TrapezoidalPrism::numberOfCells_PARAVIEW() const
{
  return ( 1 );
}




// ----------------------------------------------------------------------------
// Writes a list of points describing the sphere in a Paraview format
void TrapezoidalPrism::write_polygonsPts_PARAVIEW( ostream& f,
  	Transform const& transform, Vector3 const* translation ) const
{
  Point3 pp;
  for (int i=0;i<8;++i)
  {
    pp = transform( m_corners[i] );
    if ( translation ) pp += *translation;
    f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
  }
}




// ----------------------------------------------------------------------------
// Returns a list of points describing the trapezoidal prism in a Paraview 
// format
list<Point3> TrapezoidalPrism::get_polygonsPts_PARAVIEW( 
	Transform const& transform,
	Vector3 const* translation ) const
{
  list<Point3> ParaviewPoints;
  Point3 pp;
  for (int i=0;i<8;++i)
  {
    pp = transform( m_corners[i] );
    if ( translation ) pp += *translation;
    ParaviewPoints.push_back( pp );
  }

  return ( ParaviewPoints );
}




// ----------------------------------------------------------------------------
// Writes the trapezoidal prism in a Paraview format
void TrapezoidalPrism::write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const
{
  int count=firstpoint_globalnumber;
  for (int i=0;i<8;++i)
  {
    connectivity.push_back( count );
    ++count;
  }
  last_offset += 8;
  offsets.push_back( last_offset );
  cellstype.push_back( 12 );

  firstpoint_globalnumber += 8;
}




// ----------------------------------------------------------------------------
// Returns point[0] in face number i
Point3 TrapezoidalPrism::getFirstPointFace( int i ) const
{
  return ( m_corners[m_allFaces[i][0]] );
}




// ----------------------------------------------------------------------------
// Returns whether a point lies inside the trapezoidal prism
bool TrapezoidalPrism::isIn( Point3 const& pt ) const
{
  double slope = ( m_extent[2] - m_extent[1] ) / ( m_extent[3] + m_extent[4] );
  return ( pt[X] <= m_extent[0] && pt[X] >= - m_extent[0]
	&& pt[Y] <= slope * ( pt[Z] + m_extent[3] ) + m_extent[1] 
	&& pt[Y] >= - slope * ( pt[Z] + m_extent[3] ) - m_extent[1]
	&& pt[Z] <= m_extent[4] && pt[Z] >= - m_extent[3] );
}




// ----------------------------------------------------------------------------
// Returns the bounding volume to trapezoidal prism
BVolume* TrapezoidalPrism::computeBVolume( unsigned int type ) const
{
  BVolume* bvol = NULL;
  if ( type == 1 ) // OBB
  {
    Vector3 ext( m_extent[0], m_extent[2], m_extent[4] );
    bvol = new OBB( ext, Matrix() );
  }
  else if ( type == 2 ) // OBC
  {
    double a[2];
    int axis = ( a[X] = fabs( m_extent[0] ) ) < ( a[Y] = fabs( m_extent[2] ) )
      ? Y : X;
    int i = a[axis] < fabs( m_extent[4] ) ? Z : axis;
    Vector3 e( 0., 0., 0. );
    e[i] = 1.;
    double h = 2. * m_extent[i];
    double r = sqrt( m_extent[0] * m_extent[0]
    	+ m_extent[2] * m_extent[2] +  m_extent[4] * m_extent[4]
	- m_extent[i] * m_extent[i] );

    bvol = new OBC( r, h, e );
  }

  return( bvol );
}



// ----------------------------------------------------------------------------
// Performs advanced comparison of the two trapezoidal prisms and returns 
// whether they match
bool TrapezoidalPrism::equalType_level2( Convex const* other ) const
{
  // We know that other points to a trapezoidal prism, we dynamically cast it 
  // to actual type
  TrapezoidalPrism const* other_ = 
  	dynamic_cast<TrapezoidalPrism const*>(other); 

  bool same = true;  
  double lmin = min( computeCircumscribedRadius(),
  	other_->computeCircumscribedRadius() ); 
  
  for (size_t j=0;j<8 && same;++j)
    same = ( m_corners[j].DistanceTo( (other_->m_corners)[j] ) 
      	< LOWEPS * lmin );
  
  return ( same );
} 




// ----------------------------------------------------------------------------
// Writes the trapezoidal prism in an OBJ format
void TrapezoidalPrism::write_convex_OBJ( ostream& f, 
	Transform  const& transform, size_t& firstpoint_number ) const
{
  Point3 pp;

  // Vertices  
  for (int i=0;i<8;++i)
  {
    pp = transform( m_corners[i] );
    f << "v " << GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
		pp[X] ) << " " << 
	GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
		pp[Y] ) << " " <<
	GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
		pp[Z] ) << " " << endl;	
  }
  
  // Faces
  f << "f " << firstpoint_number + 4 << " "
  	<< firstpoint_number + 5 << " "
	<< firstpoint_number + 6 << " "
  	<< firstpoint_number + 7 << endl;
  f << "f " << firstpoint_number + 5 << " "
  	<< firstpoint_number + 1 << " "
	<< firstpoint_number + 2 << " "
  	<< firstpoint_number + 6 << endl;	
  f << "f " << firstpoint_number + 1 << " "
  	<< firstpoint_number << " "
	<< firstpoint_number + 3 << " "
  	<< firstpoint_number + 2 << endl;	
  f << "f " << firstpoint_number << " "
  	<< firstpoint_number + 4 << " "
	<< firstpoint_number + 7 << " "
  	<< firstpoint_number + 3 << endl;
  f << "f " << firstpoint_number + 4 << " "
  	<< firstpoint_number << " "
	<< firstpoint_number + 1 << " "
  	<< firstpoint_number + 5 << endl;	
  f << "f " << firstpoint_number + 3 << " "
  	<< firstpoint_number + 7 << " "
	<< firstpoint_number + 6 << " "
  	<< firstpoint_number + 2 << endl;
	
  firstpoint_number += 8;
}
