#include "Cylinder.hh"

int Cylinder::m_visuNodeNbOnPer = 36;


// ----------------------------------------------------------------------------
// Constructor with radius and height as input parameters
Cylinder::Cylinder( double r, double h ) 
  : m_radius( r )
  , m_halfHeight( h / 2. ) 
{}




// ----------------------------------------------------------------------------
// Constructor with an input stream
Cylinder::Cylinder( istream &fileIn ) 
{
  readShape( fileIn );
}




// -----------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
Cylinder::Cylinder( DOMNode* root )
{ 
  m_radius = ReaderXML::getNodeAttr_Double( root, "Radius" );
  m_halfHeight = ReaderXML::getNodeAttr_Double( root, "Height") / 2.;
}


  

// ----------------------------------------------------------------------------
// Destructor
Cylinder::~Cylinder() 
{}




// ----------------------------------------------------------------------
// Returns the convex type
ConvexType Cylinder::getConvexType() const 
{
  return ( CYLINDER );
}




// ----------------------------------------------------------------------------
// Computes the inertia tensor and the inverse of the inertia tensor
bool Cylinder::BuildInertia( double* inertia, double* inertia_1 ) const 
{
  inertia[1] = inertia[2] = inertia[4] = 0.0;
  double constant = 0.5 * m_halfHeight * m_radius * m_radius * PI;
  inertia[0] = inertia[5] = constant 
  	* ( 4.0 * m_halfHeight * m_halfHeight / 3.0 + m_radius * m_radius );
  inertia[3] = 2.0 * constant * m_radius * m_radius;

  inertia_1[1] = inertia_1[2] = inertia_1[4] = 0.0;
  inertia_1[5] = inertia_1[0] = 1.0 / inertia[0];
  inertia_1[3] = 1.0 / inertia[3];
  return true;
}




// ----------------------------------------------------------------------------
// Returns the circumscribed radius of the reference sphere,
// i.e., without applying any transformation
double Cylinder::computeCircumscribedRadius() const 
{
  return ( sqrt( m_radius * m_radius + m_halfHeight * m_halfHeight ) );
}




// ----------------------------------------------------------------------------
// Returns a clone of the cylinder
Convex* Cylinder::clone() const
{
  return ( new Cylinder( m_radius, 2. * m_halfHeight ) );
}




// ----------------------------------------------------------------------------
// Returns the cylinder volume
double Cylinder::getVolume() const
{
  return ( 2 * m_halfHeight * PI * m_radius * m_radius );
}




// ----------------------------------------------------------------------------
// Cylinder support function, returns the support point P, i.e. the
// point on the surface of the cylinder that satisfies max(P.v)
Point3 Cylinder::support( Vector3 const& v ) const 
{
  double norm = Norm( v );
  if ( norm > EPSILON ) 
  {
    double s = sqrt( v[X] * v[X] + v[Z] * v[Z] );
    if ( s > EPSILON ) 
    {
      double d = m_radius / s;  
      return ( Point3( v[X] * d, v[Y] < 0. ? - m_halfHeight : m_halfHeight, 
      	v[Z] * d ) );
    } 
    else
      return ( Point3( 0., v[Y] < 0. ? - m_halfHeight : m_halfHeight, 0. ) );
  } 
  else 
    return ( Point3() );
}




// ----------------------------------------------------------------------------
// Returns a vector of points describing the envelope of the
// cylinder. Here simply returns 3 points as follows: center of bottom circular
// face, center of top circular face and an arbitrary point on the lateral
// surface of the cylinder
vector<Point3> Cylinder::getEnvelope() const
{
  Point3 point(0.,0.,0.);
  vector<Point3> enveloppe(3,point);
  enveloppe[0][Y] = - m_halfHeight;
  enveloppe[1][Y] = - m_halfHeight;
  enveloppe[1][X] = m_radius;  
  enveloppe[2][Y] = m_halfHeight;    
  return ( enveloppe );
}




// ----------------------------------------------------------------------------
// Returns the number of vertices/corners or a code corresponding to
// a specific convex shape. Here returns the code 777
int Cylinder::getNbCorners() const
{
  return ( 777 );
}




// ----------------------------------------------------------------------------
// Returns a pointer to a 2D array describing the relationship between the face
// indices and the point indices. Returns a null pointer as a convention
vector< vector<int> > const* Cylinder::getFaces() const
{
  vector< vector<int> >* allFaces = NULL;
  return ( allFaces );
}




// ----------------------------------------------------------------------------
// Output operator
void Cylinder::writeShape( ostream& fileOut ) const 
{
  fileOut << "*Cylinder " << m_radius << " " << 2.0 * m_halfHeight << " *END";
}



// ----------------------------------------------------------------------------
// Input operator
void Cylinder::readShape( istream& fileIn ) 
{
  fileIn >> m_radius >> m_halfHeight;
  m_halfHeight /= 2.0;
}




// ----------------------------------------------------------------------------
// Returns the number of points to write the cylinder in a Paraview format
int Cylinder::numberOfPoints_PARAVIEW() const
{
  return ( 2 * m_visuNodeNbOnPer + 2 );
}




// ----------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the cylinder in a 
// Paraview format
int Cylinder::numberOfCells_PARAVIEW() const
{
  return ( m_visuNodeNbOnPer );  
}




// ----------------------------------------------------------------------------
// Writes a list of points describing the cylinder in a Paraview format
void Cylinder::write_polygonsPts_PARAVIEW( ostream& f, 
  	Transform const& transform, Vector3 const* translation ) const
{
  Point3 pp,p;
  double dtheta = 2.* PI / m_visuNodeNbOnPer;
  
  // Couronne inferieure
  p[Y] = - m_halfHeight;
  for (int i=0;i<m_visuNodeNbOnPer;++i)
  {
    p[X] = m_radius * cos ( i * dtheta );
    p[Z] = m_radius * sin ( i * dtheta ); 
    pp = transform( p );
    if ( translation ) pp += *translation;    
    f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
  }
  
  // Couronne superieure
  p[Y] = m_halfHeight;
  for (int i=0;i<m_visuNodeNbOnPer;++i)
  {
    p[X] = m_radius * cos ( i * dtheta );
    p[Z] = m_radius * sin ( i * dtheta ); 
    pp = transform( p );
    if ( translation ) pp += *translation;    
    f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
  }
  
  // Centre inferieur
  p[X] = 0.;
  p[Y] = - m_halfHeight;
  p[Z] = 0.;
  pp = transform( p );

  if ( translation ) pp += *translation;  
  f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
  
  // Centre superieur
  p[Y] = m_halfHeight;
  pp = transform( p );
  if ( translation ) pp += *translation;    
  f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
}




// ----------------------------------------------------------------------------
// Returns a list of points describing the cylinder in a Paraview format 
list<Point3> Cylinder::get_polygonsPts_PARAVIEW( Transform const& transform,
  	Vector3 const* translation ) const
{
  list<Point3> ParaviewPoints;
  Point3 pp,p;
  double dtheta = 2.* PI / m_visuNodeNbOnPer;
  
  // Couronne inferieure
  p[Y] = - m_halfHeight;  
  for (int i=0;i<m_visuNodeNbOnPer;++i)
  {
    p[X] = m_radius * cos ( i * dtheta );
    p[Z] = m_radius * sin ( i * dtheta ); 
    pp = transform( p );
    if ( translation ) pp += *translation;
    ParaviewPoints.push_back( pp );
  }
  
  // Couronne superieure
  p[Y] = m_halfHeight;
  for (int i=0;i<m_visuNodeNbOnPer;++i)
  {
    p[X] = m_radius * cos ( i * dtheta );
    p[Z] = m_radius * sin ( i * dtheta ); 
    pp = transform( p );
    if ( translation ) pp += *translation;
    ParaviewPoints.push_back( pp );
  }
  
  // Centre inferieur
  p[X] = 0.;
  p[Y] = - m_halfHeight;
  p[Z] = 0.;
  pp = transform( p );
  if ( translation ) pp += *translation;
  ParaviewPoints.push_back( pp );
  
  // Centre superieur
  p[Y] = m_halfHeight;
  pp = transform( p );
  if ( translation ) pp += *translation;
  ParaviewPoints.push_back( pp );
  
  return ( ParaviewPoints ); 
}




// ----------------------------------------------------------------------------
// Writes the cylinder in a Paraview format
void Cylinder::write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const
{
  for (int i=0;i<m_visuNodeNbOnPer-1;++i)
  {
    connectivity.push_back( firstpoint_globalnumber + i );
    connectivity.push_back( firstpoint_globalnumber + i + 1 );    
    connectivity.push_back( firstpoint_globalnumber + 2 * m_visuNodeNbOnPer );
    connectivity.push_back( firstpoint_globalnumber + i + m_visuNodeNbOnPer);
    connectivity.push_back( firstpoint_globalnumber + i + m_visuNodeNbOnPer 
    	+ 1 );
    connectivity.push_back( firstpoint_globalnumber + 2 * m_visuNodeNbOnPer 
    	+ 1 );        
    last_offset += 6;     
    offsets.push_back( last_offset );
    cellstype.push_back( 13 );
  }
  connectivity.push_back( firstpoint_globalnumber + m_visuNodeNbOnPer - 1 );
  connectivity.push_back( firstpoint_globalnumber );    
  connectivity.push_back( firstpoint_globalnumber + 2 * m_visuNodeNbOnPer );    
  connectivity.push_back( firstpoint_globalnumber + 2 * m_visuNodeNbOnPer - 1 );
  connectivity.push_back( firstpoint_globalnumber + m_visuNodeNbOnPer );
  connectivity.push_back( firstpoint_globalnumber + 2 * m_visuNodeNbOnPer + 1 );
  last_offset += 6;    
  offsets.push_back( last_offset );
  cellstype.push_back( 13 );  
  
  firstpoint_globalnumber += 2 * m_visuNodeNbOnPer + 2;
}




// ----------------------------------------------------------------------------
// Returns whether a point lies inside the cylinder
bool Cylinder::isIn( Point3 const& pt ) const
{
  return ( pt[Y] > - m_halfHeight && pt[Y] < m_halfHeight 
  	&& sqrt( pt[X] * pt[X] + pt[Z] * pt[Z] ) < m_radius );
}  
