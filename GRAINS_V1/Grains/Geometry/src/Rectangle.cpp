#include "Rectangle.hh"



// ----------------------------------------------------------------------------
// Constructor with the dimensions of the wall
Rectangle::Rectangle( double LX, double LY )
  : m_LX( LX )
  , m_LY( LY )
{
  setCorners();
}



// ----------------------------------------------------------------------------
// Constructor with an input stream
Rectangle::Rectangle( istream& fileIn )
{
  readShape( fileIn );
}



// ----------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
Rectangle::Rectangle( DOMNode* root )
{
  m_LX =  ReaderXML::getNodeAttr_Double( root, "LX" );
  m_LY =  ReaderXML::getNodeAttr_Double( root, "LY" );

  setCorners();
}



// ----------------------------------------------------------------------------
// Destructor
Rectangle::~Rectangle()
{
  m_corners.clear();
}



// ----------------------------------------------------------------------------
// Construit les sommets du pav� et les faces
void Rectangle::setCorners()
{
  m_corners.reserve( 4 );
  Point3 pp;
  pp.setValue( -m_LX/2., -m_LY/2., 0. );
  m_corners.push_back( pp );
  pp.setValue( m_LX/2., -m_LY/2., 0. );
  m_corners.push_back( pp );
  pp.setValue( -m_LX/2., m_LY/2., 0. );
  m_corners.push_back( pp );
  pp.setValue( m_LX/2., m_LY/2., 0. );
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
    return ( Point3( v[X] < 0. ? -m_LX/2. : m_LX/2.,
                     v[Y] < 0. ? -m_LY/2. : m_LY/2.,
                     0. ) );
  }
  else
    return ( Point3() );
}




// ----------------------------------------------------------------------------
// Returns a vector of points describing the envelope of the rectangle.
// Here simply returns the point (0,0,0) as a convention
vector<Point3> Rectangle::getEnvelope() const
{
  vector<Point3> envelope;
  Point3 point( 0.,0.,0. );
  envelope.push_back( point );
  return ( envelope );
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
  fileOut << "*Rectangle " << m_LX << " " << m_LY << " END*";
}




// ----------------------------------------------------------------------------
// Input operator
void Rectangle::readShape( istream& fileIn )
{
  cerr << "Program Error :\n" << "Rectangle::readShape non accessible.\n";
  exit( 3 );
}




// ----------------------------------------------------------------------------
// Returns the number of points to write the rectangle in a Paraview format
int Rectangle::numberOfPoints_PARAVIEW() const
{
  return ( 4 );
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
    pp = transform( m_corners[i] );
    if ( translation ) pp += *translation;
    f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
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
    pp = transform( m_corners[i] );
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
  int count = firstpoint_globalnumber;
  for ( int i = 0; i < 4; ++i )
  {
    connectivity.push_back( count );
    ++count;
  }
  last_offset += 4;
  offsets.push_back( last_offset );
  cellstype.push_back( 8 );

  firstpoint_globalnumber += 4;
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
