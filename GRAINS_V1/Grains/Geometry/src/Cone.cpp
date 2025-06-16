#include "Cone.hh"
#include "GrainsExec.hh"

int Cone::m_visuNodeNbOnPer = 32;

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




// -----------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
Cone::Cone( DOMNode* root )
{
  m_bottomRadius = ReaderXML::getNodeAttr_Double( root, "Radius" );
  m_quarterHeight = ReaderXML::getNodeAttr_Double( root, "Height") / 4.;
  m_sinAngle = m_bottomRadius / sqrt( m_bottomRadius * m_bottomRadius 
  	+ 16. * m_quarterHeight * m_quarterHeight ); 
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




// ----------------------------------------------------------------------------
// Returns a vector of points describing the surface of the
// cone. TO DO
vector<Point3> Cone::getEnvelope() const
{
  Point3 point(0.,0.,0.);
  vector<Point3> surface( 3, point );
  surface[0][Y] = - m_quarterHeight;
  surface[1][Y] = - m_quarterHeight;
  surface[1][X] = m_bottomRadius;
  surface[2][Y] = 3. * m_quarterHeight;
  return ( surface );
}




// ----------------------------------------------------------------------------
// Returns the number of vertices/corners or a code corresponding to
// a specific convex shape. Here returns the code 888
int Cone::getNbCorners() const
{
  return ( 888 );
}




// ----------------------------------------------------------------------------
// Returns a pointer to a 2D array describing the relationship between the face
// indices and the point indices. Returns a null pointer as a convention
vector< vector<int> > const* Cone::getFaces() const
{
  vector< vector<int> >* allFaces = NULL;
  return ( allFaces );
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
  fileOut << "*Cone " << m_bottomRadius << " " << 4.0 * m_quarterHeight 
  	<< " *END";
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
// Returns the number of points to write the cone in a Paraview format
int Cone::numberOfPoints_PARAVIEW() const
{
  return ( m_visuNodeNbOnPer + 2 );
}




// ----------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the cone in a
// Paraview format
int Cone::numberOfCells_PARAVIEW() const
{
  return ( m_visuNodeNbOnPer );
}




// ----------------------------------------------------------------------------
// Writes a list of points describing the cone in a Paraview format
void Cone::write_polygonsPts_PARAVIEW( ostream& f,
  	Transform const& transform, Vector3 const* translation ) const
{
  Point3 pp, p;
  double dtheta = 2.* PI / m_visuNodeNbOnPer;

  // Disk rim
  p[Y] = - m_quarterHeight;
  for (int i=0;i<m_visuNodeNbOnPer;++i)
  {
    p[X] = m_bottomRadius * cos ( i * dtheta );
    p[Z] = m_bottomRadius * sin ( i * dtheta );
    pp = transform( p );
    if ( translation ) pp += *translation;
    f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
  }

  // Disk center
  p[X] = 0.;
  p[Y] = - m_quarterHeight;
  p[Z] = 0.;
  pp = transform( p );
  if ( translation ) pp += *translation;
  f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;

  // Upper tip
  p[X] = 0.;
  p[Y] = 3. * m_quarterHeight;
  p[Z] = 0.;
  pp = transform( p );
  if ( translation ) pp += *translation;
  f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
}




// ----------------------------------------------------------------------------
// Returns a list of points describing the cylinder in a Paraview format
list<Point3> Cone::get_polygonsPts_PARAVIEW( Transform const& transform,
  	Vector3 const* translation ) const
{
  list<Point3> ParaviewPoints;
  Point3 pp,p;
  double dtheta = 2.* PI / m_visuNodeNbOnPer;

  // Disk rim
  p[Y] = - m_quarterHeight;
  for (int i=0;i<m_visuNodeNbOnPer;++i)
  {
    p[X] = m_bottomRadius * cos ( i * dtheta );
    p[Z] = m_bottomRadius * sin ( i * dtheta );
    pp = transform( p );
    if ( translation ) pp += *translation;
    ParaviewPoints.push_back( pp );
  }

  // Disk center
  p[X] = 0.;
  p[Y] = - m_quarterHeight;
  p[Z] = 0.;
  pp = transform( p );
  if ( translation ) pp += *translation;
  ParaviewPoints.push_back( pp );  

  // Upper tip
  p[X] = 0.;
  p[Y] = 3. * m_quarterHeight;
  p[Z] = 0.;
  pp = transform( p );
  if ( translation ) pp += *translation;
  ParaviewPoints.push_back( pp );

  return ( ParaviewPoints );
}




// ----------------------------------------------------------------------------
// Writes the cylinder in a Paraview format
void Cone::write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const
{
  for (int i=0;i<m_visuNodeNbOnPer-1;++i)
  {
    connectivity.push_back( firstpoint_globalnumber + i );
    connectivity.push_back( firstpoint_globalnumber + i + 1 );
    connectivity.push_back( firstpoint_globalnumber + m_visuNodeNbOnPer );
    connectivity.push_back( firstpoint_globalnumber + m_visuNodeNbOnPer + 1 );
    last_offset += 4;
    offsets.push_back( last_offset );
    cellstype.push_back( 10 );
  }
  connectivity.push_back( firstpoint_globalnumber + m_visuNodeNbOnPer - 1 );
  connectivity.push_back( firstpoint_globalnumber );
  connectivity.push_back( firstpoint_globalnumber + m_visuNodeNbOnPer );
  connectivity.push_back( firstpoint_globalnumber + m_visuNodeNbOnPer + 1 );
  last_offset += 4;
  offsets.push_back( last_offset );
  cellstype.push_back( 10 );

  firstpoint_globalnumber += m_visuNodeNbOnPer + 2;
}




// ----------------------------------------------------------------------------
// Returns whether a point lies inside the cone
bool Cone::isIn( Point3 const& pt ) const
{
  return ( pt[Y] >= - m_quarterHeight && pt[Y] <= 3. * m_quarterHeight
  	&& sqrt( pt[X] * pt[X] + pt[Z] * pt[Z] ) <= 
		( 3. * m_quarterHeight - pt[Y] ) * m_sinAngle 
			/ sqrt( 1. - m_sinAngle * m_sinAngle ) );
}  




// ----------------------------------------------------------------------------
// Returns the bounding volume to cone
BVolume* Cone::computeBVolume( unsigned int type ) const
{
  BVolume* bvol = NULL;
  if ( type == 1 ) // OBB
    bvol = new OBB( Vector3( m_bottomRadius, 
	3. * m_quarterHeight, m_bottomRadius ), Matrix() );
  else if ( type == 2 ) // OBC
    bvol = new OBC( m_bottomRadius, 6. * m_quarterHeight, Vector3(0., 1., 0.) );

  return( bvol );
}




// ----------------------------------------------------------------------------
// Performs advanced comparison of the two cones and returns whether 
// they match
bool Cone::equalType_level2( Convex const* other ) const
{
  // We know that other points to a Cone, we dynamically cast it to actual type
  Cone const* other_ = dynamic_cast<Cone const*>(other);
  
  double lmin = min( computeCircumscribedRadius(),
  	other_->computeCircumscribedRadius() );    

  bool same = ( 
  	fabs( m_bottomRadius - other_->m_bottomRadius ) <  LOWEPS * lmin 
	&& fabs( m_quarterHeight - other_->m_quarterHeight ) <  LOWEPS * lmin
	&& fabs( m_sinAngle - other_->m_sinAngle ) <  LOWEPS );
  
  return ( same );
}



 
// ----------------------------------------------------------------------------
// Sets the number of point over the cone perimeter for Paraview 
// post-processing, i.e., controls the number of facets in the cone 
// reconstruction in Paraview
void Cone::SetvisuNodeNbOverPer( int nbpts )
{
  m_visuNodeNbOnPer = nbpts;
}




// ----------------------------------------------------------------------------
// Writes the cone in an OBJ format
void Cone::write_convex_OBJ( ostream& f, Transform  const& transform,
    	size_t& firstpoint_number ) const
{
  Point3 pp, p;
  double dtheta = 2.* PI / m_visuNodeNbOnPer;

  // Vertices  
  // Disk rim
  p[Y] = - m_quarterHeight;
  for (int i=0;i<m_visuNodeNbOnPer;++i)
  {
    p[X] = m_bottomRadius * cos ( i * dtheta );
    p[Z] = m_bottomRadius * sin ( i * dtheta );
    pp = transform( p );
    f << "v " << GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
		pp[X] ) << " " << 
	GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
		pp[Y] ) << " " <<
	GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
		pp[Z] ) << " " << endl;	
  }

  // Disk center
  p[X] = 0.;
  p[Y] = - m_quarterHeight;
  p[Z] = 0.;
  pp = transform( p );
  f << "v " << GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
			pp[X] ) << " " << 
	GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
			pp[Y] ) << " " <<
	GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
			pp[Z] ) << " " << endl;

  // Upper tip
  p[X] = 0.;
  p[Y] = 3. * m_quarterHeight;
  p[Z] = 0.;
  pp = transform( p );
  f << "v " << GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
			pp[X] ) << " " << 
	GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
			pp[Y] ) << " " <<
	GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
			pp[Z] ) << " " << endl;
  
  // Faces 
  // Triangular lateral faces 
  for (int i=0;i<m_visuNodeNbOnPer-1;++i)
  {
    f << "f "<< firstpoint_number + i << " "
    	<< firstpoint_number + i + 1 << " "
    	<< firstpoint_number + m_visuNodeNbOnPer + 1 << endl;
  }
  f << "f " << firstpoint_number + m_visuNodeNbOnPer - 1 << " "
  	<< firstpoint_number << " "
  	<< firstpoint_number + m_visuNodeNbOnPer + 1 << endl;

  // Triangular bottom faces 
  for (int i=0;i<m_visuNodeNbOnPer-1;++i)
  {
    f << "f "<< firstpoint_number + i << " "
    	<< firstpoint_number + i + 1 << " "
    	<< firstpoint_number + m_visuNodeNbOnPer << endl;
  }
  f << "f " << firstpoint_number + m_visuNodeNbOnPer - 1 << " "
  	<< firstpoint_number << " "
  	<< firstpoint_number + m_visuNodeNbOnPer << endl;

  firstpoint_number += m_visuNodeNbOnPer + 2;    			  
}
