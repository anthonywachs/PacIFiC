#include "Rectangle.hh"
#include "GrainsExec.hh"


// NOTE: while the Rectangle is a bounded 2D plane in a 3d space and therefore
// has 0 width, we represent it as a thin plate in Paraview using an arbitrary
// width of 2*LOWEPS. This facilitates the 3D visualization of Rectangle in 
// Paraview


// ----------------------------------------------------------------------------
// Constructor with the dimensions of the wall
Rectangle::Rectangle( double LX, double LY, double pw, 
	RectangleVisuExpansion pwt )
  : m_LX( LX )
  , m_LY( LY )
  , m_ParaviewWidth( pw )
  , m_ParaviewWidth_type( pwt )
{
  setCorners();
}



// ----------------------------------------------------------------------------
// Constructor with an input stream
Rectangle::Rectangle( istream& fileIn )
{
  readShape( fileIn );
  setCorners();  
}



// ----------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
Rectangle::Rectangle( DOMNode* root )
{
  m_LX =  ReaderXML::getNodeAttr_Double( root, "LX" );
  m_LY =  ReaderXML::getNodeAttr_Double( root, "LY" );
  if ( ReaderXML::hasNodeAttr( root, "ParaviewWidth" ) )
    m_ParaviewWidth = ReaderXML::getNodeAttr_Double( root, "ParaviewWidth" );
  else m_ParaviewWidth = 2 * LOWEPS;
  m_ParaviewWidth_type = RVE_CENTERED;
  if ( ReaderXML::hasNodeAttr( root, "ParaviewExpand" ) )
  {
    string sexp = ReaderXML::getNodeAttr_String( root, "ParaviewExpand" );
    if ( sexp == "Z+" ) m_ParaviewWidth_type = RVE_ZPLUS;
    else if ( sexp == "Z-" ) m_ParaviewWidth_type = RVE_ZMINUS;
  }
  setCorners();
}



// ----------------------------------------------------------------------------
// Destructor
Rectangle::~Rectangle()
{
  m_corners.clear();
}



// ----------------------------------------------------------------------------
// Sets the corner/vertex coordinates
void Rectangle::setCorners()
{
  m_corners.reserve( 4 );
  Point3 pp;
  pp.setValue( -m_LX/2., -m_LY/2., 0. );
  m_corners.push_back( pp );
  pp.setValue( m_LX/2., -m_LY/2., 0. );
  m_corners.push_back( pp );
  pp.setValue( m_LX/2., m_LY/2., 0. );
  m_corners.push_back( pp );
  pp.setValue( -m_LX/2., m_LY/2., 0. );
  m_corners.push_back( pp );
}




// ----------------------------------------------------------------------------
// Returns the convex type
ConvexType Rectangle::getConvexType() const
{
  return ( RECTANGLE2D );
}




// ----------------------------------------------------------------------------
// Computes the inertia tensor and the inverse of the inertia tensor
bool Rectangle::BuildInertia( double* inertia, double* inertia_1 ) const
{
  // Active 2D plane is XY -> rotation around Z
  inertia[1] = inertia[2] = inertia[4]= 0.0;
  inertia[0] = m_LX * pow( m_LY, 3. ) / 12.;
  inertia[3] = pow( m_LX, 3. ) * m_LY / 12.;
  inertia[5] = m_LX * m_LY * ( m_LX * m_LX + m_LY * m_LY ) / 12.;

  inertia_1[1] = inertia_1[2] = inertia_1[4] = 0.0;
  inertia_1[0] = 1.0 / inertia[0];
  inertia_1[3] = 1.0 / inertia[3];
  inertia_1[5] = 1.0 / inertia[5];

  return ( true );
}




// ----------------------------------------------------------------------------
// Returns the circumscribed radius of the reference rectangle,
// i.e., without applying any transformation
double Rectangle::computeCircumscribedRadius() const
{
  return ( sqrt( m_LX * m_LX + m_LY * m_LY ) / 2. );
}




// ----------------------------------------------------------------------------
// Returns a clone of the rectangle
Convex* Rectangle::clone() const
{
  return ( new Rectangle( m_LX, m_LY ) );
}




// ----------------------------------------------------------------------------
// Rectangle support function, returns the support point P, i.e. the
// point on the surface of the rectangle that satisfies max(P.v)
Point3 Rectangle::support( Vector3 const& v ) const
{
  double s = Norm( v );
  if ( s > EPSILON )
  {
    return ( Point3( v[X] < 0. ? -m_LX / 2. : m_LX / 2.,
                     v[Y] < 0. ? -m_LY / 2. : m_LY / 2.,
                     0. ) );
  }
  else
    return ( Point3() );
}




// ----------------------------------------------------------------------------
// Returns a vector of points describing the surface of the rectangle.
// Here simply returns the point (0,0,0) as a convention
vector<Point3> Rectangle::getEnvelope() const
{
  vector<Point3> surface;
  Point3 point( 0., 0., 0. );
  surface.push_back( point );
  return ( surface );
}




// ----------------------------------------------------------------------------
// Returns the number of vertices/corners or a code corresponding to
// a specific convex shape. Here 444 as a convention
int Rectangle::getNbCorners() const
{
  return ( 444 );
}




// ----------------------------------------------------------------------------
// Returns the rectangle surface area
double Rectangle::getVolume()const
{
  return ( m_LX * m_LY );
}




// ----------------------------------------------------------------------------
// Output operator
void Rectangle::writeShape( ostream& fileOut ) const
{
  fileOut << "*Rectangle " << m_LX << " " << m_LY << " " <<	
  	m_ParaviewWidth << " " << static_cast<int>(m_ParaviewWidth_type) 
	<< " *END";
}




// ----------------------------------------------------------------------------
// Input operator
void Rectangle::readShape( istream& fileIn )
{
  int nn;
  fileIn >> m_LX >> m_LY >> m_ParaviewWidth >> nn;
  switch( nn )
  {
    case 0: m_ParaviewWidth_type = RVE_CENTERED; break;
    case 1: m_ParaviewWidth_type = RVE_ZPLUS; break;    
    case 2: m_ParaviewWidth_type = RVE_ZMINUS; break;    
  }
}




// ----------------------------------------------------------------------------
// Returns the number of points to write the rectangle in a Paraview format
int Rectangle::numberOfPoints_PARAVIEW() const
{
  return ( 8 );
}




// ----------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the rectangle
// in a Paraview format
int Rectangle::numberOfCells_PARAVIEW() const
{
  return ( 1 );
}




// ----------------------------------------------------------------------------
// Returns a pointer to a 2D array describing the relationship between the face
// indices and the point indices. Returns a null pointer as a convention
vector< vector<int> > const* Rectangle::getFaces() const
{
  vector< vector<int> >* allFaces = NULL;
  return ( allFaces );
}



// ----------------------------------------------------------------------------
// Writes a list of points describing the rectangle in a Paraview format
void Rectangle::write_polygonsPts_PARAVIEW( ostream& f,
  	Transform const& transform, Vector3 const* translation ) const
{
  Point3 pp;
  for ( int i = 0; i < 4; ++i )
  {
    pp = m_corners[i];
    switch( m_ParaviewWidth_type )
    {
      case RVE_CENTERED: pp[Z] = - 0.5 * m_ParaviewWidth; break;
      case RVE_ZPLUS: pp[Z] = 0.; break;
      case RVE_ZMINUS: pp[Z] = - m_ParaviewWidth; break;      
    }
    pp = transform( pp );
    if ( translation ) pp += *translation;
    f << pp[X] << " " << pp[Y] << " " << 
	GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS, pp[Z] )
	<< endl;
  }
  for ( int i = 0; i < 4; ++i )
  {
    pp = m_corners[i];
    switch( m_ParaviewWidth_type )
    {
      case RVE_CENTERED: pp[Z] = 0.5 * m_ParaviewWidth; break;
      case RVE_ZPLUS: pp[Z] = m_ParaviewWidth; break;
      case RVE_ZMINUS: pp[Z] = 0.; break;      
    }
    pp = transform( pp );
    if ( translation ) pp += *translation;
    f << pp[X] << " " << pp[Y] << " " << 
	GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS, pp[Z] )
	<< endl;
  }  
}




// ----------------------------------------------------------------------------
// Returns a list of points describing the rectangle in a Paraview format
list<Point3> Rectangle::get_polygonsPts_PARAVIEW( Transform const& transform,
	Vector3 const* translation ) const
{
  list<Point3> ParaviewPoints;
  Point3 pp;
  for ( int i = 0; i < 4; ++i )
  {
    pp = m_corners[i];
    switch( m_ParaviewWidth_type )
    {
      case RVE_CENTERED: pp[Z] = - 0.5 * m_ParaviewWidth; break;
      case RVE_ZPLUS: pp[Z] = 0.; break;
      case RVE_ZMINUS: pp[Z] = - m_ParaviewWidth; break;      
    }
    pp = transform( pp );
    if ( translation ) pp += *translation;
    ParaviewPoints.push_back( pp );
  }
  for ( int i = 0; i < 4; ++i )
  {
    pp = m_corners[i];
    switch( m_ParaviewWidth_type )
    {
      case RVE_CENTERED: pp[Z] = 0.5 * m_ParaviewWidth; break;
      case RVE_ZPLUS: pp[Z] = m_ParaviewWidth; break;
      case RVE_ZMINUS: pp[Z] = 0.; break;      
    }
    pp = transform( pp );
    if ( translation ) pp += *translation;
    ParaviewPoints.push_back( pp );
  }  
  return ( ParaviewPoints );
}




// ----------------------------------------------------------------------------
// Writes the rectangle in a Paraview format
void Rectangle::write_polygonsStr_PARAVIEW( list<int>& connectivity,
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
// Returns whether the convex shape is a rectangle
bool Rectangle::isRectangle() const
{
  return ( true );
}




// ----------------------------------------------------------------------------
// Returns whether a point lies in the rectangle
bool Rectangle::isIn( Point3 const& pt ) const
{
  return ( fabs( pt[Z] ) < EPSILON &&
           fabs( pt[Y] ) < m_LY/2. &&
           fabs( pt[X] ) < m_LX/2. );
}




// ----------------------------------------------------------------------------
// Performs advanced comparison of the two rectangles and returns whether 
// they match
bool Rectangle::equalType_level2( Convex const* other ) const
{
  // We know that other points to a Rectangle, we dynamically cast it to 
  // actual type
  Rectangle const* other_ = dynamic_cast<Rectangle const*>(other); 
  
  double lmin = min( computeCircumscribedRadius(),
  	other_->computeCircumscribedRadius() ); 
  
  bool same = ( fabs( m_LX - other_->m_LX ) <  LOWEPS * lmin 
	&& fabs( m_LY - other_->m_LY ) <  LOWEPS * lmin );  

  for (size_t j=0;j<4 && same;++j)
      same = ( m_corners[j].DistanceTo( (other_->m_corners)[j] ) 
      	< LOWEPS * lmin );
  
  return ( same );
} 




// ----------------------------------------------------------------------------
// Returns the bounding volume to rectangle
BVolume* Rectangle::computeBVolume( unsigned int type ) const
{
  BVolume* bvol = NULL;
  if ( type == 1 ) // OBB
    bvol = new OBB( Vector3( m_LX, m_LY, EPSILON ), Matrix() );
  else if ( type == 2 ) // OBC
  {
    bvol = new OBC( sqrt( m_LX * m_LX + m_LY * m_LY ) / 2., 
                    EPSILON,
                    Vector3( 0., 0., 1. ) );
  }
  return( bvol );
}
