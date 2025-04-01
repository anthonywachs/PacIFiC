#include "GrainsExec.hh"
#include "GrainsBuilderFactory.hh"
#include "Box.hh"
#include "BBox.hh"
#include <sstream>

// Static attribute
vector< vector<int> > Box::m_allFaces;


// ----------------------------------------------------------------------------
// Constructor with a vector containing the edge half-lengths as
// input parameters
Box::Box( Vector3 const& extent_ )
  : m_extent( extent_ )
  , m_corners2D_XY( NULL )
{
  if ( GrainsBuilderFactory::getContext() == DIM_2 ) m_extent[Z] = 0.;
  setCornersFaces();
}




// ----------------------------------------------------------------------------
// Constructor with edge length as input parameters
Box::Box( double x, double y, double z )
  : m_extent( fabs(x)/2., fabs(y)/2., fabs(z)/2. )
  , m_corners2D_XY( NULL )
{
  if ( GrainsBuilderFactory::getContext() == DIM_2 ) m_extent[Z] = 0.;
  setCornersFaces();
}




// ----------------------------------------------------------------------------
// Constructor with an input stream
Box::Box( istream& fileIn )
  : m_corners2D_XY( NULL )
{
  readShape( fileIn );
  setCornersFaces();
}




// ----------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
Box::Box( DOMNode* root )
  : m_corners2D_XY( NULL )
{
  m_extent[X] = ReaderXML::getNodeAttr_Double( root, "LX" ) / 2.;
  m_extent[Y] = ReaderXML::getNodeAttr_Double( root, "LY" ) / 2.;
  m_extent[Z] = ReaderXML::getNodeAttr_Double( root, "LZ" ) / 2.;
  if ( GrainsBuilderFactory::getContext() == DIM_2 ) m_extent[Z] = 0.;
  setCornersFaces(); 
}




// ----------------------------------------------------------------------------
// Destructor
Box::~Box()
{
  m_corners.clear();
  if ( m_corners2D_XY )
  {
    m_corners2D_XY->clear();
    delete m_corners2D_XY;
  }
}




// ----------------------------------------------------------------------
// Returns the convex type
ConvexType Box::getConvexType() const
{
  return ( BOX );
}




// ----------------------------------------------------------------------------
// Sets the corner/vertex coordinates and the face/vertices numbering
void Box::setCornersFaces()
{
  m_corners.reserve( 8 );
  Point3 sommet;
  sommet.setValue( - m_extent[X], - m_extent[Y], - m_extent[Z] );
  m_corners.push_back(sommet);
  sommet.setValue( m_extent[X], - m_extent[Y], - m_extent[Z] );
  m_corners.push_back(sommet);
  sommet.setValue( m_extent[X], - m_extent[Y], m_extent[Z] );
  m_corners.push_back(sommet);
  sommet.setValue( - m_extent[X], - m_extent[Y], m_extent[Z] );
  m_corners.push_back(sommet);
  sommet.setValue( - m_extent[X], m_extent[Y], - m_extent[Z] );
  m_corners.push_back(sommet);
  sommet.setValue( m_extent[X], m_extent[Y], - m_extent[Z] );
  m_corners.push_back(sommet);
  sommet.setValue( m_extent[X], m_extent[Y], m_extent[Z] );
  m_corners.push_back(sommet);
  sommet.setValue( - m_extent[X], m_extent[Y], m_extent[Z] );
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

  // In case of a 2D simulation
  if ( GrainsBuilderFactory::getContext() == DIM_2 )
  {
    m_corners2D_XY = new vector<Point3>( 4, sommet );
    (*m_corners2D_XY)[0] = m_corners[0] ;
    (*m_corners2D_XY)[0][Z] = 0.;
    (*m_corners2D_XY)[1] = m_corners[1] ;
    (*m_corners2D_XY)[1][Z] = 0.;
    (*m_corners2D_XY)[2] = m_corners[5] ;
    (*m_corners2D_XY)[2][Z] = 0.;
    (*m_corners2D_XY)[3] = m_corners[7] ;
    (*m_corners2D_XY)[3][Z] = 0.;
  }
}




// ----------------------------------------------------------------------------
// Returns a clone of the box
Convex* Box::clone() const
{
  return ( new Box( 2.0 * m_extent[X], 2.0 * m_extent[Y], 
  	2.0 * m_extent[Z] ) );
}




// ----------------------------------------------------------------------
// Computes the inertia tensor and the inverse of the inertia tensor
bool Box::BuildInertia( double* inertia, double* inertia_1 ) const
{
  inertia[1] = inertia[2] = inertia[4] = 0.0;
  if ( GrainsBuilderFactory::getContext() == DIM_2 )
  {
    inertia[0] = 4.0 * m_extent[X] * m_extent[Y] * m_extent[Y] 
    	* m_extent[Y]  / 3.0;
    inertia[3] = 4.0 * m_extent[Y] * m_extent[X] * m_extent[X] 
    	* m_extent[X]  / 3.0;
    inertia[5] = inertia[0] + inertia[3];
  }
  else
  {
    inertia[0] = 8.0 * m_extent[X] * m_extent[Y] * m_extent[Z]
    * ( m_extent[Y] * m_extent[Y] + m_extent[Z] * m_extent[Z] ) / 3.0;
    inertia[3] = 8.0 * m_extent[X] * m_extent[Y] * m_extent[Z]
    * ( m_extent[X] * m_extent[X] + m_extent[Z] * m_extent[Z] ) / 3.0;
    inertia[5] = 8.0 * m_extent[X] * m_extent[Y] * m_extent[Z]
    * ( m_extent[Y] * m_extent[Y] + m_extent[X] * m_extent[X] ) / 3.0;
  }  

  inertia_1[1] = inertia_1[2] = inertia_1[4] = 0.0;
  inertia_1[0] = 1.0 / inertia[0];
  inertia_1[3] = 1.0 / inertia[3];
  inertia_1[5] = 1.0 / inertia[5];

  return ( true );
}




// ----------------------------------------------------------------------------
// Returns the edge half-lengths of the box
Vector3 const& Box::getExtent() const
{
  return ( m_extent );
}




// ----------------------------------------------------------------------------
// Returns the circumscribed radius of the reference box,
// i.e., without applying any transformation
double Box::computeCircumscribedRadius() const
{
  return ( Norm( m_extent ) );
}




// ----------------------------------------------------------------------------
// Returns a vector of points describing the envelope of the box
vector<Point3> Box::getEnvelope() const
{
  if ( GrainsBuilderFactory::getContext() == DIM_2 )
    return ( *m_corners2D_XY );
  else
    return ( m_corners );
}




// ----------------------------------------------------------------------------
// Returns a pointer to a 2D array describing the relationship
// between the face indices and the point indices
vector<vector<int> > const* Box::getFaces() const
{
  return ( &m_allFaces );
}




// ----------------------------------------------------------------------------
// Returns the number of vertices/corners or a code corresponding to
// a specific convex shape. Here returns the code 444 for a 2D box and 666 for
// a 3D box
int Box::getNbCorners() const
{
  if ( GrainsBuilderFactory::getContext() == DIM_2 ) return ( 444 );
  else return ( 666 );
}




// ----------------------------------------------------------------------------
// Returns the box volume
double Box::getVolume() const
{
  if ( GrainsBuilderFactory::getContext() == DIM_2 )
    return ( 4.0 * m_extent[X] * m_extent[Y] ) ;
  else
    return ( 8.0 * m_extent[X] * m_extent[Y] * m_extent[Z] );
}




// ----------------------------------------------------------------------------
// Output operator
void Box::writeShape( ostream& fileOut ) const
{
  fileOut << "*Box " << m_extent << " *END";
}




// ----------------------------------------------------------------------------
// Input operator
void Box::readShape( istream& fileIn )
{
  fileIn >> m_extent;
}




// ----------------------------------------------------------------------------
// Box support function, returns the support point P, i.e. the
// point on the surface of the box that satisfies max(P.v)
Point3 Box::support( Vector3 const& v ) const
{
  double norm = Norm( v );
  if ( norm < EPSILON )
    return ( Point3() );
  else
    return ( Point3( v[X] < 0. ? -m_extent[X] : m_extent[X],
	v[Y] < 0. ? -m_extent[Y] : m_extent[Y],
	v[Z] < 0. ? -m_extent[Z] : m_extent[Z]) );
}




// ----------------------------------------------------------------------------
// Returns the number of points to write the box in a Paraview format
int Box::numberOfPoints_PARAVIEW() const
{
  return ( 8 );
}




// ----------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the box
// in a Paraview format
int Box::numberOfCells_PARAVIEW() const
{
  return ( 1 );
}




// ----------------------------------------------------------------------------
// Writes a list of points describing the box in a Paraview format
void Box::write_polygonsPts_PARAVIEW( ostream& f,
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
// Returns a list of points describing the box in a Paraview format
list<Point3> Box::get_polygonsPts_PARAVIEW( Transform const& transform,
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
// Writes the box in a Paraview format
void Box::write_polygonsStr_PARAVIEW( list<int>& connectivity,
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
// Writes the box in an OBJ format
void Box::write_convex_OBJ( ostream& f, Transform  const& transform,
    	size_t& firstpoint_number ) const
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




// ----------------------------------------------------------------------------
// Returns point[0] in face number i
Point3 Box::getFirstPointFace( int i ) const
{
  return ( m_corners[m_allFaces[i][0]] );
}




// ----------------------------------------------------------------------------
// Returns the contact point in the box reference frame and the
// overlapping distance between the box and a sphere. If no contact, returned
// contact point is world reference frame origin (0,0,0)
Point3 Box::IntersectionPointSPHERE( Point3 const& SphereCenter,
  	double const& SphereRadius, double& overlap,
	bool warningSphereCenterInBox ) const
{
  try{
  Point3 contactPoint;
  double gx = SphereCenter[X], gy = SphereCenter[Y], gz = SphereCenter[Z];
  Vector3 distance;
  double normDistance = 0.;
  overlap = 1.;

  if ( gx > m_extent[X] )
  {
    if ( gy > m_extent[Y] )
    {
      if ( gz > m_extent[Z] )
      {
        // Distance au corner 6
	contactPoint = ContactCornerSPHERE( SphereCenter, SphereRadius, 6,
		overlap );
      }
      else if ( gz < - m_extent[Z] )
      {
        // Distance au corner 5
        contactPoint = ContactCornerSPHERE( SphereCenter, SphereRadius, 5,
		overlap );
      }
      else
      {
        // Distance a l'arete 5-6
        contactPoint = ContactEdgeSPHERE( SphereCenter, SphereRadius, 5, Z,
		overlap );
      }
    }
    else if ( gy < - m_extent[Y] )
    {
      if ( gz > m_extent[Z] )
      {
        // Distance au corner 2
	contactPoint = ContactCornerSPHERE( SphereCenter, SphereRadius, 2,
		overlap );
      }
      else if ( gz < - m_extent[Z] )
      {
        // Distance au corner 1
	contactPoint = ContactCornerSPHERE( SphereCenter, SphereRadius, 1,
		overlap );
      }
      else
      {
        // Distance a l'arete 1-2
	contactPoint = ContactEdgeSPHERE( SphereCenter, SphereRadius, 1, Z,
		overlap );
      }
    }
    else
    {
      if ( gz > m_extent[Z] )
      {
        // Distance a l'arete 2-6
	contactPoint = ContactEdgeSPHERE( SphereCenter, SphereRadius, 6, Y,
		overlap );
      }
      else if ( gz < - m_extent[Z] )
      {
        // Distance a l'arete 1-5
	contactPoint = ContactEdgeSPHERE( SphereCenter, SphereRadius, 1, Y,
		overlap );
      }
      else
      {
        // Distance a la face 1
	normDistance = gx - m_extent[X];
	if ( normDistance < SphereRadius )
        {
          overlap = normDistance - SphereRadius;
          contactPoint.setValue( m_extent[X] + 0.5 * overlap, SphereCenter[Y],
	  	SphereCenter[Z] );
        }
      }
    }
  }
  else if ( gx < - m_extent[X] )
  {
    if ( gy > m_extent[Y] )
    {
      if ( gz > m_extent[Z] )
      {
        // Distance au corner 7
	contactPoint = ContactCornerSPHERE( SphereCenter, SphereRadius, 7,
		overlap );
      }
      else if ( gz < - m_extent[Z] )
      {
        // Distance au corner 4
	contactPoint = ContactCornerSPHERE( SphereCenter, SphereRadius, 4,
		overlap );
      }
      else
      {
        // Distance a l'arete 4-7
	contactPoint = ContactEdgeSPHERE( SphereCenter, SphereRadius, 4, Z,
		overlap );
      }
    }
    else if ( gy < -m_extent[Y] )
    {
      if ( gz > m_extent[Z] )
      {
        // Distance au corner 3
	contactPoint = ContactCornerSPHERE( SphereCenter, SphereRadius, 3,
		overlap );
      }
      else if ( gz < - m_extent[Z] )
      {
        // Distance au corner 0
	contactPoint = ContactCornerSPHERE( SphereCenter, SphereRadius, 0,
		overlap );
      }
      else
      {
        // Distance a l'arete 0-3
	contactPoint = ContactEdgeSPHERE( SphereCenter, SphereRadius, 0, Z,
		overlap );
      }
    }
    else
    {
      if ( gz > m_extent[Z] )
      {
        // Distance a l'arete 3-7
	contactPoint = ContactEdgeSPHERE( SphereCenter, SphereRadius, 3, Y,
		overlap );
      }
      else if ( gz < - m_extent[Z] )
      {
        // Distance a l'arete 0-4
	contactPoint = ContactEdgeSPHERE( SphereCenter, SphereRadius, 0, Y,
		overlap );
      }
      else
      {
        // Distance a la face 3
	normDistance = - m_extent[X] - gx;
	if ( normDistance < SphereRadius )
        {
          overlap = normDistance - SphereRadius;
          contactPoint.setValue( -m_extent[X] - 0.5 * overlap, SphereCenter[Y],
	  	SphereCenter[Z] );
        }
      }
    }
  }
  else
  {
    if ( gy > m_extent[Y] )
    {
      if ( gz > m_extent[Z] )
      {
        // Distance a l'arete 6-7
	contactPoint = ContactEdgeSPHERE( SphereCenter, SphereRadius, 6, X,
		overlap );
      }
      else if ( gz < - m_extent[Z] )
      {
        // Distance a l'arete 4-5
	contactPoint = ContactEdgeSPHERE( SphereCenter, SphereRadius, 4, X,
		overlap );
      }
      else
      {
        // Distance a la face 0
	normDistance = gy - m_extent[Y];
	if ( normDistance < SphereRadius )
        {
          overlap = normDistance - SphereRadius;
          contactPoint.setValue( SphereCenter[X], m_extent[Y] + 0.5 * overlap,
	  	SphereCenter[Z] );
        }
      }
    }
    else if ( gy < - m_extent[Y] )
    {
      if ( gz > m_extent[Z] )
      {
        // Distance a l'arete 3-2
	contactPoint = ContactEdgeSPHERE( SphereCenter, SphereRadius, 3, X,
		overlap );
      }
      else if ( gz < - m_extent[Z] )
      {
        // Distance a l'arete 0-1
	contactPoint = ContactEdgeSPHERE( SphereCenter, SphereRadius, 0, X,
		overlap );
      }
      else
      {
        // Distance a la face 2
	normDistance = - m_extent[Y] - gy;
	if ( normDistance < SphereRadius )
        {
          overlap = normDistance - SphereRadius;
          contactPoint.setValue( SphereCenter[X], - m_extent[Y] 
	  	- 0.5 * overlap, SphereCenter[Z] );
        }
      }
    }
    else
    {
      if ( gz > m_extent[Z] )
      {
        // Distance a la face 5
	normDistance = gz - m_extent[Z];
	if ( normDistance < SphereRadius )
        {
          overlap = normDistance - SphereRadius;
          contactPoint.setValue( SphereCenter[X], SphereCenter[Y],
	  	m_extent[Z] + 0.5 * overlap );
        }
      }
      else if ( gz < - m_extent[Z] )
      {
        // Distance a la face 4
	normDistance = - m_extent[Z] - gz;
	if ( normDistance < SphereRadius )
        {
          overlap = normDistance - SphereRadius;
          contactPoint.setValue( SphereCenter[X], SphereCenter[Y],
	  	- m_extent[Z] - 0.5 * overlap );
        }
      }
      else
      {
        if ( warningSphereCenterInBox )
	{
	  cout << "Warning: sphere center in box in "
	  	"Box::IntersectionPointSPHERE" << endl;
	  GrainsExec::m_exception_Contact = true;
          throw ContactError();
	}
	else overlap = -10000.;
      }
    }
  }

  return ( contactPoint );
  }
  catch (ContactError&){
    throw ContactError();
  }
}




// ----------------------------------------------------------------------------
// Returns the contact point in the box reference frame and the
// overlapping distance between a sphere and a box corner. If no contact,
// returned contact point is world reference frame origin (0,0,0)
Point3 Box::ContactCornerSPHERE( Point3 const& SphereCenter,
  	double const& SphereRadius,
	int cornerNumber,
	double& overlap ) const
{
  Point3 contactPoint;
  Vector3 distance = SphereCenter - m_corners[cornerNumber];
  double normDistance = Norm(distance);
  if ( normDistance < SphereRadius )
  {
    overlap = normDistance - SphereRadius;
    contactPoint = SphereCenter - ( 1. + 0.5 * overlap / normDistance )
    	* distance;
  }

  return ( contactPoint );
}




// ----------------------------------------------------------------------------
// Returns the contact point in the box reference frame and the
// overlapping distance between a sphere and a box edge. If no contact,
// returned contact point is world reference frame origin (0,0,0)
Point3 Box::ContactEdgeSPHERE( Point3 const& SphereCenter,
  	double const& SphereRadius,
	int cornerNumber,
	int projectionDirection,
	double& overlap ) const
{
  Point3 contactPoint, SphereCenterProjected( SphereCenter );
  SphereCenterProjected[projectionDirection] =
  	m_corners[cornerNumber][projectionDirection];
  Vector3 distance = SphereCenterProjected - m_corners[cornerNumber];
  double normDistance = Norm(distance);
  if ( normDistance < SphereRadius )
  {
    overlap = normDistance - SphereRadius;
    contactPoint = SphereCenter - ( 1. + 0.5 * overlap / normDistance )
    	* distance;
  }

  return ( contactPoint );
}




// ----------------------------------------------------------------------------
// This method sends back the projection of sphere center on obstacle's faces,
// edges or corners (based on local position of the sphere with respect to the
// obstacle).
// Pay attention that only in the case of projection on faces, the gap (normal
// distance between sphere and the wall) is  modified. In case of projection on
// edges and corners, the gap is set to zero (can be modified later) */
Point3 Box::ProjectedPointSPHERE( const Point3& SphereCenter,
  	const double& SphereRadius,
	double& gap ) const
{
  Point3 ProjectedPoint;
  double gx = SphereCenter[X], gy = SphereCenter[Y], gz = SphereCenter[Z];
  double normDistance = 0.;
  gap = 0.;

  if ( gx > m_extent[X] )
  {
    if ( gy > m_extent[Y] )
    {
      if ( gz > m_extent[Z] )
      {
        // Distance au corner 6
	ProjectedPoint = m_corners[6];
      }
      else if ( gz < - m_extent[Z] )
      {
        // Distance au corner 5
        ProjectedPoint = m_corners[5];
      }
      else
      {
        // Distance a l'arete 5-6
        ProjectedPoint = SphereCenter;
        ProjectedPoint[Z] = m_corners[5][Z];
      }
    }
    else if ( gy < - m_extent[Y] )
    {
      if ( gz > m_extent[Z] )
      {
        // Distance au corner 2
	ProjectedPoint = m_corners[2];
      }
      else if ( gz < - m_extent[Z] )
      {
        // Distance au corner 1
	ProjectedPoint = m_corners[1];
      }
      else
      {
        // Distance a l'arete 1-2
	ProjectedPoint = SphereCenter;
        ProjectedPoint[Z] = m_corners[1][Z];
      }
    }
    else
    {
      if ( gz > m_extent[Z] )
      {
        // Distance a l'arete 2-6
	ProjectedPoint = SphereCenter;
        ProjectedPoint[Y] = m_corners[6][Y];
      }
      else if ( gz < - m_extent[Z] )
      {
        // Distance a l'arete 1-5
	ProjectedPoint = SphereCenter;
        ProjectedPoint[Y] = m_corners[1][Y];
      }
      else
      {
        // Distance a la face 1
	normDistance = gx - m_extent[X];
	if ( normDistance > SphereRadius )
        {
          gap = normDistance - SphereRadius;
          ProjectedPoint.setValue( m_extent[X], SphereCenter[Y],
	  	SphereCenter[Z] );
        }
      }
    }
  }
  else if ( gx < - m_extent[X] )
  {
    if ( gy > m_extent[Y] )
    {
      if ( gz > m_extent[Z] )
      {
        // Distance au corner 7
	ProjectedPoint = m_corners[7];
      }
      else if ( gz < - m_extent[Z] )
      {
        // Distance au corner 4
	ProjectedPoint = m_corners[4];
      }
      else
      {
        // Distance a l'arete 4-7
	ProjectedPoint = SphereCenter;
        ProjectedPoint[Z] = m_corners[4][Z];
      }
    }
    else if ( gy < - m_extent[Y] )
    {
      if ( gz > m_extent[Z] )
      {
        // Distance au corner 3
	ProjectedPoint = m_corners[3];
      }
      else if ( gz < - m_extent[Z] )
      {
        // Distance au corner 0
	ProjectedPoint = m_corners[0];
      }
      else
      {
        // Distance a l'arete 0-3
	ProjectedPoint = SphereCenter;
        ProjectedPoint[Z] = m_corners[0][Z];
      }
    }
    else
    {
      if ( gz > m_extent[Z] )
      {
        // Distance a l'arete 3-7
	ProjectedPoint = SphereCenter;
        ProjectedPoint[Y] = m_corners[3][Y];
      }
      else if ( gz < - m_extent[Z] )
      {
        // Distance a l'arete 0-4
	ProjectedPoint = SphereCenter;
        ProjectedPoint[Y] = m_corners[0][Y];
      }
      else
      {
        // Distance a la face 3
	normDistance = - m_extent[X] - gx;
	if ( normDistance > SphereRadius )
        {
          gap = normDistance - SphereRadius;
          ProjectedPoint.setValue( - m_extent[X], SphereCenter[Y],
	  	SphereCenter[Z] );
        }
      }
    }
  }
  else
  {
    if ( gy > m_extent[Y] )
    {
      if ( gz > m_extent[Z] )
      {
        // Distance a l'arete 6-7
	ProjectedPoint = SphereCenter;
        ProjectedPoint[X] = m_corners[6][X];
      }
      else if ( gz < - m_extent[Z] )
      {
        // Distance a l'arete 4-5
	ProjectedPoint = SphereCenter;
        ProjectedPoint[X] = m_corners[4][X];
      }
      else
      {
        // Distance a la face 0
	normDistance = gy - m_extent[Y];
	if ( normDistance > SphereRadius )
        {
          gap = normDistance - SphereRadius;
          ProjectedPoint.setValue( SphereCenter[X], m_extent[Y],
	  	SphereCenter[Z] );
        }
      }
    }
    else if ( gy < - m_extent[Y] )
    {
      if ( gz > m_extent[Z] )
      {
        // Distance a l'arete 3-2
	ProjectedPoint = SphereCenter;
        ProjectedPoint[X] = m_corners[3][X];
      }
      else if ( gz < - m_extent[Z] )
      {
        // Distance a l'arete 0-1
	ProjectedPoint = SphereCenter;
        ProjectedPoint[X] = m_corners[0][X];
      }
      else
      {
        // Distance a la face 2
	normDistance = - m_extent[Y] - gy;
	if ( normDistance > SphereRadius )
        {
          gap = normDistance - SphereRadius;
          ProjectedPoint.setValue( SphereCenter[X], - m_extent[Y],
	  	SphereCenter[Z] );
        }
      }
    }
    else
    {
      if ( gz > m_extent[Z] )
      {
        // Distance a la face 5
	normDistance = gz - m_extent[Z];
	if ( normDistance > SphereRadius )
        {
          gap = normDistance - SphereRadius;
          ProjectedPoint.setValue( SphereCenter[X], SphereCenter[Y],
	  	m_extent[Z] );
        }
      }
      else if ( gz < - m_extent[Z] )
      {
        // Distance a la face 4
	normDistance = - m_extent[Z] - gz;
	if ( normDistance > SphereRadius )
        {
          gap = normDistance - SphereRadius;
          ProjectedPoint.setValue( SphereCenter[X], SphereCenter[Y],
	  	- m_extent[Z] );
        }
      }
      else
      {
        cout << "Warning: sphere center in box in Box::ProjectedPointSPHERE"
		<< endl;
      }
    }
  }

  return ( ProjectedPoint );

}




// ----------------------------------------------------------------------------
// Returns whether a point lies inside the box
bool Box::isIn( Point3 const& pt ) const
{
  return ( pt[X] <= m_extent[X] && pt[X] >= - m_extent[X]
	&& pt[Y] <= m_extent[Y] && pt[Y] >= - m_extent[Y]
	&& pt[Z] <= m_extent[Z] && pt[Z] >= - m_extent[Z] );
}




// ----------------------------------------------------------------------------
// Returns the bounding volume to box
BVolume* Box::computeBVolume( unsigned int type ) const
{
  BVolume* bvol = NULL;
  if ( type == 1 ) // OBB
    bvol = new OBB( m_extent, Matrix() );
  else if ( type == 2 ) // OBC
  {
    double xy = fabs( m_extent[X] - m_extent[Y] );
    double xz = fabs( m_extent[X] - m_extent[Z] );
    double yz = fabs( m_extent[Y] - m_extent[Z] );
    // pick from xy and xz, store to xy
    int zAxis = xy < xz ? Z : Y;
    xy = xy < xz ? xy : xz;
    // pick from xy and yz
    zAxis = xy < yz ? zAxis : X;

    Vector3 e( 0., 0., 0. );
    e[ zAxis ] = 1.;
    double h = 2. * m_extent[zAxis];
    double r = sqrt( m_extent[X] * m_extent[X] + 
                     m_extent[Y] * m_extent[Y] +
                     m_extent[Z] * m_extent[Z] -
                     m_extent[zAxis] * m_extent[zAxis] );

    bvol = new OBC( r, h, e );
  }

  return( bvol );
}



// ----------------------------------------------------------------------------
// Performs advanced comparison of the two boxes and returns whether 
// they match
bool Box::equalType_level2( Convex const* other ) const
{
  // We know that other points to a Box, we dynamically cast it to actual type
  Box const* other_ = dynamic_cast<Box const*>(other); 

  size_t dim = ( GrainsBuilderFactory::getContext() == DIM_2 ? 2: 3 ), j;
  bool same = true;  
  double lmin = min( computeCircumscribedRadius(),
  	other_->computeCircumscribedRadius() ); 
  
  if ( dim == 2 )
    for (j=0;j<4 && same;++j)
      same = ( (*m_corners2D_XY)[j].DistanceTo( (*other_->m_corners2D_XY)[j] ) 
      	< LOWEPS * lmin );
  else
    for (j=0;j<8 && same;++j)
      same = ( m_corners[j].DistanceTo( (other_->m_corners)[j] ) 
      	< LOWEPS * lmin );
  
  return ( same );
} 
