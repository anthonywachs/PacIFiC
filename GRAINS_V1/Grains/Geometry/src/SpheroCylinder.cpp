#include "GrainsExec.hh"
#include "GrainsBuilderFactory.hh"
#include "SpheroCylinder.hh"
#include "BBox.hh"
#include <sstream>

int SpheroCylinder::m_visuNodeNbPerQar = 8;


// ----------------------------------------------------------------------------
// Constructor with edge length as input parameters
SpheroCylinder::SpheroCylinder( double r, double h )
  : m_radius( r )
  , m_height( h )
{}




// ----------------------------------------------------------------------------
// Constructor with an input stream
SpheroCylinder::SpheroCylinder( istream& fileIn )
{
  readShape( fileIn );
}




// ----------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
SpheroCylinder::SpheroCylinder( DOMNode* root )
{
  m_radius = ReaderXML::getNodeAttr_Double( root, "Radius" );
  m_height = ReaderXML::getNodeAttr_Double( root, "Height" );
}




// ----------------------------------------------------------------------------
// Destructor
SpheroCylinder::~SpheroCylinder()
{}




// ----------------------------------------------------------------------
// Returns the convex type
ConvexType SpheroCylinder::getConvexType() const
{
  return ( SPHEROCYLINDER );
}




// ----------------------------------------------------------------------------
// Returns a clone of the SpheroCylinder
Convex* SpheroCylinder::clone() const
{
  return ( new SpheroCylinder( m_radius, m_height ) );
}




// ----------------------------------------------------------------------
// Computes the inertia tensor and the inverse of the inertia tensor
bool SpheroCylinder::BuildInertia( double* inertia, double* inertia_1 ) const
{  
  inertia[1] = inertia[2] = inertia[4] = 0.;
  inertia_1[1] = inertia_1[2] = inertia_1[4] = 0.0;  
  
  inertia[0] = inertia[5] = 
  	( 1. / 12. ) * PI * m_height * m_radius * m_radius
  	 	* ( 3. * m_radius * m_radius + m_height * m_height 
			+ 4. * m_radius * m_height ) 
	+ ( 8. / 15. ) * PI * pow( m_radius, 5. )
	+ ( 1. / 2. ) * PI * m_height * pow( m_radius, 4. );
  inertia[3] = ( 1. / 30. ) * PI * pow( m_radius, 4. )
  	* ( 16. * m_radius + 15. * m_height );
	
  inertia_1[0] = 1.0 / inertia[0];
  inertia_1[3] = 1.0 / inertia[3];
  inertia_1[5] = 1.0 / inertia[5];
  
  return ( true );
}




// ----------------------------------------------------------------------------
// Returns the circumscribed radius of the reference SpheroCylinder,
// i.e., without applying any transformation
double SpheroCylinder::computeCircumscribedRadius() const
{
  return ( 0.5 * m_height + m_radius );
}




// ----------------------------------------------------------------------------
// Returns a vector of points describing the surface of the
// cylinder. Here simply returns 3 points as follows: center of bottom circular
// face of the elementary cylinder, an arbitrary point on the lateral surface 
// of the elementary cylinder and center of top circular face of the elementary 
// cylinder
vector<Point3> SpheroCylinder::getEnvelope() const
{
  Point3 point( 0., 0., 0. );
  vector<Point3> surface( 3, point );
  surface[0][Y] = - 0.5 * m_height;
  surface[1][Y] = - 0.5 * m_height;
  surface[1][X] = m_radius;
  surface[2][Y] = 0.5 * m_height;
  return ( surface );
}




// ----------------------------------------------------------------------------
// Returns a pointer to a 2D array describing the relationship between the face
// indices and the point indices. Returns a null pointer as a convention
vector< vector<int> > const* SpheroCylinder::getFaces() const
{
  vector< vector<int> >* allFaces = NULL;
  return ( allFaces );
}




// ----------------------------------------------------------------------------
// Returns the number of vertices/corners or a code corresponding to
// a specific convex shape. Here returns the code 3333
int SpheroCylinder::getNbCorners() const
{
  return ( 3333 );
}




// ----------------------------------------------------------------------------
// Returns the SpheroCylinder volume
double SpheroCylinder::getVolume() const
{
  return ( PI * m_radius * m_radius * 
  	( m_height + ( 4. / 3. ) * m_radius ) );
}




// ----------------------------------------------------------------------------
// Output operator
void SpheroCylinder::writeShape( ostream& fileOut ) const
{
  fileOut << "*SpheroCylinder " << m_radius << " " << m_height << " *END";
}




// ----------------------------------------------------------------------------
// Input operator
void SpheroCylinder::readShape( istream& fileIn )
{
  fileIn >> m_radius >> m_height;
}




// ----------------------------------------------------------------------------
// SpheroCylinder support function, returns the support point P, i.e. the
// point on the surface of the SpheroCylinder that satisfies max(P.v)
Point3 SpheroCylinder::support( Vector3 const& v ) const
{
  double s = Norm( v );
  if ( s > EPSILON ) 
  {
    double r = m_radius / s;
    if ( fabs( v[Y] ) < EPSILON )
      return ( Point3( v[X] * r, 0., v[Z] * r ) );      
    else
      return ( Point3( v[X] * r, v[Y] * r 
    	+ ( v[Y] > 0. ? 0.5 : -0.5 ) * m_height, v[Z] * r ) );
  } 
  else 
    return ( Point3() );    
}




// ----------------------------------------------------------------------------
// Returns the number of points to write the SpheroCylinder in a Paraview 
// format
int SpheroCylinder::numberOfPoints_PARAVIEW() const
{
  return ( 2 * ( 4 * m_visuNodeNbPerQar * m_visuNodeNbPerQar + 2 ) );
}




// ----------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the SpheroCylinder
// in a Paraview format
int SpheroCylinder::numberOfCells_PARAVIEW() const
{
  return ( 2 * ( 4 * m_visuNodeNbPerQar * m_visuNodeNbPerQar ) 
  	+ 4 * m_visuNodeNbPerQar );
}




// ----------------------------------------------------------------------------
// Writes a list of points describing the sphere in a Paraview format
void SpheroCylinder::write_polygonsPts_PARAVIEW( ostream& f,
  	Transform const& transform, Vector3 const* translation ) const
{	 	
  double angle = PI / ( 2. * m_visuNodeNbPerQar ) ;
  double angleY = 0., local_radius = 0.;
  int k, i, ptsPerlevel = 4 * m_visuNodeNbPerQar;
  Point3 pp, pptrans;
  
  // Top sphere 
  // Regular points on the surface
  for ( k = m_visuNodeNbPerQar - 1; k < 2*m_visuNodeNbPerQar-1 ; ++k ) 
  {  
    angleY = - PI / 2. + ( k + 1 ) * angle;
    local_radius = m_radius * cos( angleY );
    pp[Y] = m_radius * sin( angleY ) + 0.5 * m_height;
    for ( i = 0; i < ptsPerlevel ; ++i )
    {
      pp[X] = local_radius * cos( i * angle );
      pp[Z] = local_radius * sin( i * angle );
      pptrans = transform( pp );
      if ( translation ) pptrans += *translation;
      f << pptrans[X] << " " << pptrans[Y] << " " << pptrans[Z] << endl;
    }
  }
   
  pp[X] = 0.;
  pp[Z] = 0.;	
  // Top point
  pp[Y] = m_radius + 0.5 * m_height;
  pptrans = transform( pp );
  if ( translation ) pptrans += *translation;
  f << pptrans[X] << " " << pptrans[Y] << " " << pptrans[Z] << endl;
	
  // Gravity center
  pp[Y] = 0.5 * m_height;
  pptrans = transform( pp );
  if ( translation ) pptrans += *translation;
  f << pptrans[X] << " " << pptrans[Y] << " " << pptrans[Z] << endl; 
  
  
  // Bottom sphere	 	  
  // Regular points on the surface
  for ( k = m_visuNodeNbPerQar - 1; k >=0 ; --k ) 
  {  
    angleY = - PI / 2. + ( k + 1 ) * angle;
    local_radius = m_radius * cos( angleY );
    pp[Y] = m_radius * sin( angleY ) - 0.5 * m_height;
    for ( i = 0; i < ptsPerlevel ; ++i )
    {
      pp[X] = local_radius * cos( i * angle );
      pp[Z] = local_radius * sin( i * angle );
      pptrans = transform( pp );
      if ( translation ) pptrans += *translation;
      f << pptrans[X] << " " << pptrans[Y] << " " << pptrans[Z] << endl;
    }
  }
   
  pp[X] = 0.;
  pp[Z] = 0.;
  // Bottom point
  pp[Y] = - m_radius - 0.5 * m_height;
  pptrans = transform( pp );
  if ( translation ) pptrans += *translation;
  f << pptrans[X] << " " << pptrans[Y] << " " << pptrans[Z] << endl;	
	
  // Gravity center
  pp[Y] = - 0.5 * m_height;
  pptrans = transform( pp );
  if ( translation ) pptrans += *translation;
  f << pptrans[X] << " " << pptrans[Y] << " " << pptrans[Z] << endl;
  
  
  // Cylinder: no additional point needed
}




// ----------------------------------------------------------------------------
// Returns a list of points describing the SpheroCylinder in a Paraview 
// format
list<Point3> SpheroCylinder::get_polygonsPts_PARAVIEW( 
	Transform const& transform,
	Vector3 const* translation ) const
{
  list<Point3> ParaviewPoints;
  double angle = PI / ( 2. * m_visuNodeNbPerQar ) ;
  double angleY = 0., local_radius = 0.;
  int k, i, ptsPerlevel =  4 * m_visuNodeNbPerQar;
  Point3 pp, pptrans;

  // Top sphere
  // Regular points on the surface
  for ( k = m_visuNodeNbPerQar - 1; k < 2*m_visuNodeNbPerQar-1 ; ++k ) 
  {  
    angleY = -PI / 2. + ( k + 1 ) * angle;
    local_radius = m_radius * cos( angleY );
    pp[Y] = m_radius * sin( angleY ) + 0.5 * m_height;
    for ( i = 0; i < ptsPerlevel ; ++i )
    {
      pp[X] = local_radius * cos( i * angle );
      pp[Z] = local_radius * sin( i * angle );
      pptrans = transform( pp );
      if ( translation ) pptrans += *translation;
      ParaviewPoints.push_back( pptrans );
    }
  }
   
  pp[X] = 0.;
  pp[Z] = 0.;	
  // Top point
  pp[Y] = m_radius + 0.5 * m_height;
  pptrans = transform( pp );
  if ( translation ) pptrans += *translation;
  ParaviewPoints.push_back( pptrans );
	
  // Gravity center
  pp[Y] = 0.5 * m_height;
  pptrans = transform( pp );  
  if ( translation ) pptrans += *translation;
  ParaviewPoints.push_back( pptrans );
  
  
  // Bottom sphere  
  // Regular points on the surface
  for ( k = m_visuNodeNbPerQar - 1; k >=0 ; --k ) 
  {  
    angleY = -PI / 2. + ( k + 1 ) * angle;
    local_radius = m_radius * cos( angleY );
    pp[Y] = m_radius * sin( angleY ) - 0.5 * m_height;
    for ( i = 0; i < ptsPerlevel ; ++i )
    {
      pp[X] = local_radius * cos( i * angle );
      pp[Z] = local_radius * sin( i * angle );
      pptrans = transform( pp );
      if ( translation ) pptrans += *translation;
      ParaviewPoints.push_back( pptrans );
    }
  }
   
  pp[X] = 0.;
  pp[Z] = 0.;
  // Bottom point
  pp[Y] = - m_radius - 0.5 * m_height;
  pptrans = transform( pp );
  if ( translation ) pptrans += *translation;
  ParaviewPoints.push_back( pptrans );
	
  // Gravity center
  pp[Y] = - 0.5 * m_height;
  pptrans = transform( pp );  
  if ( translation ) pptrans += *translation;
  ParaviewPoints.push_back( pptrans );
  
  
  // Cylinder: no additional point needed     
      
  return ( ParaviewPoints );
}




// ----------------------------------------------------------------------------
// Writes the SpheroCylinder in a Paraview format
void SpheroCylinder::write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const
{
  int i, k, ptsPerlevel = 4 * m_visuNodeNbPerQar;
  
  // Top sphere
  int Top_number = firstpoint_globalnumber + ptsPerlevel * m_visuNodeNbPerQar,  
  	Top_GC_number = firstpoint_globalnumber 
		+ ptsPerlevel * m_visuNodeNbPerQar + 1;
  // Regular cells: Pyramid
  for ( k = 0; k < m_visuNodeNbPerQar-1 ; ++k ) 
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
      connectivity.push_back( Top_GC_number );
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
    connectivity.push_back( Top_GC_number );
    last_offset += 5;
    offsets.push_back( last_offset );
    cellstype.push_back( 14 );    
  }   
  
  // Top cells: tetraedron  
  for ( i = 0; i < ptsPerlevel-1 ; ++i )
  {
    connectivity.push_back( firstpoint_globalnumber
    	+ ( m_visuNodeNbPerQar - 1 ) * ptsPerlevel + i );
    connectivity.push_back( firstpoint_globalnumber
    	+ ( m_visuNodeNbPerQar - 1 ) * ptsPerlevel + i + 1 );
    connectivity.push_back( Top_number );	
    connectivity.push_back( Top_GC_number );
    last_offset += 4;
    offsets.push_back( last_offset );
    cellstype.push_back( 10 ); 
  }
  connectivity.push_back( firstpoint_globalnumber
  	+  m_visuNodeNbPerQar * ptsPerlevel - 1 );
  connectivity.push_back( firstpoint_globalnumber
  	+ ( m_visuNodeNbPerQar - 1 ) * ptsPerlevel );
  connectivity.push_back( Top_number );	
  connectivity.push_back( Top_GC_number );
  last_offset += 4;
  offsets.push_back( last_offset );
  cellstype.push_back( 10 ); 
	

  // Bottom sphere
  int shift = firstpoint_globalnumber + ptsPerlevel * m_visuNodeNbPerQar + 2 ;  
  int Bottom_number = shift + ptsPerlevel * m_visuNodeNbPerQar,  
  	Bottom_GC_number = shift + ptsPerlevel * m_visuNodeNbPerQar + 1; 
  // Regular cells: Pyramid
  for ( k = 0; k < m_visuNodeNbPerQar-1 ; ++k ) 
  {  
    for ( i = 0; i < ptsPerlevel-1 ; ++i )
    {
      connectivity.push_back( shift + k * ptsPerlevel + i ); 
      connectivity.push_back( shift + k * ptsPerlevel + i + 1);
      connectivity.push_back( shift + ( k + 1) * ptsPerlevel + i + 1 );	
      connectivity.push_back( shift + ( k + 1 ) * ptsPerlevel + i );
      connectivity.push_back( Bottom_GC_number );
      last_offset += 5;
      offsets.push_back( last_offset );
      cellstype.push_back( 14 );		
    }
    connectivity.push_back( shift + k * ptsPerlevel + ptsPerlevel - 1 );
    connectivity.push_back( shift + k * ptsPerlevel );
    connectivity.push_back( shift + ( k + 1 ) * ptsPerlevel );
    connectivity.push_back( shift + ( k + 1 ) * ptsPerlevel + ptsPerlevel - 1 );
    connectivity.push_back( Bottom_GC_number );
    last_offset += 5;
    offsets.push_back( last_offset );
    cellstype.push_back( 14 );    
  }

  // Bottom cells: tetraedron
  for ( i = 0; i < ptsPerlevel-1 ; ++i )
  {
    connectivity.push_back( shift + ( m_visuNodeNbPerQar - 1 ) * ptsPerlevel 
    	+ i );
    connectivity.push_back( shift + ( m_visuNodeNbPerQar - 1 ) * ptsPerlevel 
    	+ i + 1 );
    connectivity.push_back( Bottom_number );	
    connectivity.push_back( Bottom_GC_number );
    last_offset += 4;
    offsets.push_back( last_offset );
    cellstype.push_back( 10 ); 
  }
  connectivity.push_back( shift +  m_visuNodeNbPerQar * ptsPerlevel - 1 );
  connectivity.push_back( shift + ( m_visuNodeNbPerQar - 1 ) * ptsPerlevel );
  connectivity.push_back( Bottom_number );	
  connectivity.push_back( Bottom_GC_number );
  last_offset += 4;
  offsets.push_back( last_offset );
  cellstype.push_back( 10 ); 
  
  
  // Cylinder	
  for ( i = 0; i < ptsPerlevel-1 ; ++i )
  {
    connectivity.push_back( firstpoint_globalnumber + i );
    connectivity.push_back( firstpoint_globalnumber + i + 1 );
    connectivity.push_back( Top_GC_number );
    connectivity.push_back( shift + i );
    connectivity.push_back( shift + i + 1 );
    connectivity.push_back( Bottom_GC_number );
    last_offset += 6;
    offsets.push_back( last_offset );
    cellstype.push_back( 13 );
  }
  connectivity.push_back( firstpoint_globalnumber + ptsPerlevel - 1 );
  connectivity.push_back( firstpoint_globalnumber );
  connectivity.push_back( Top_GC_number );
  connectivity.push_back( shift + ptsPerlevel - 1 );
  connectivity.push_back( shift );
  connectivity.push_back( Bottom_GC_number );
  last_offset += 6;
  offsets.push_back( last_offset );
  cellstype.push_back( 13 );


  firstpoint_globalnumber += 2 * ( ptsPerlevel * m_visuNodeNbPerQar + 2 ) ;
}




// ----------------------------------------------------------------------------
// Returns whether a point lies inside the SpheroCylinder
bool SpheroCylinder::isIn( Point3 const& pt ) const
{
  bool isIn = false;
  double halfHeight = 0.5 * m_height, r2 = m_radius * m_radius;
  
  if ( pt[Y] >= - halfHeight - m_radius && pt[Y] < - halfHeight )
    isIn = pt[X] * pt[X] 
    	+ ( pt[Y] + halfHeight ) * ( pt[Y] + halfHeight )
	+ pt[Z] * pt[Z] <= r2;
  else if ( pt[Y] >= - halfHeight && pt[Y] <= halfHeight )
    isIn = pt[X] * pt[X] + pt[Z] * pt[Z] <= r2;
  else if ( pt[Y] > halfHeight && pt[Y] <= halfHeight + m_radius )
    isIn = pt[X] * pt[X] 
    	+ ( pt[Y] - halfHeight ) * ( pt[Y] - halfHeight )
	+ pt[Z] * pt[Z] <= r2;
	
  return ( isIn );
}




// ----------------------------------------------------------------------------
// Returns the bounding volume to SpheroCylinder
BVolume* SpheroCylinder::computeBVolume( unsigned int type ) const
{
  BVolume* bvol = NULL;
  if ( type == 1 ) // OBB
  {
    Vector3 const& extent = Vector3( m_radius, m_radius + 0.5 * m_height, 
    	m_radius );
    bvol = new OBB( extent, Matrix() );
  }

  else if ( type == 2 ) // OBC
  {
    Vector3 const& e = Vector3( 0., 1., 0. );
    bvol = new OBC( m_radius, m_radius + m_height, e );
  }
  
  return( bvol );
}



// ----------------------------------------------------------------------------
// Performs advanced comparison of the two SpheroCylinders and returns 
// whether they match
bool SpheroCylinder::equalType_level2( Convex const* other ) const
{
  // We know that other points to a SpheroCylinder, we dynamically cast it 
  // to actual type
  SpheroCylinder const* other_ = 
  	dynamic_cast<SpheroCylinder const*>(other); 

  double lmin = min( computeCircumscribedRadius(),
  	other_->computeCircumscribedRadius() ); 
  
  bool same = ( 
  	fabs( m_radius - other_->m_radius ) <  LOWEPS * lmin 
	&& fabs( m_height - other_->m_height ) <  LOWEPS * lmin );
  
  return ( same );
} 




// ----------------------------------------------------------------------------
// Sets the number of point per quarter of the elementary spherocylinder 
// perimeter for Paraview post-processing, i.e., controls the number of facets 
// in the spherocylinder reconstruction in Paraview
void SpheroCylinder::SetvisuNodeNbPerQar( int nbpts )
{
  m_visuNodeNbPerQar = nbpts;
}




// ----------------------------------------------------------------------------
// Writes the SpheroCylinder in an OBJ format
void SpheroCylinder::write_convex_OBJ( ostream& f, Transform  const& transform,
    	size_t& firstpoint_number ) const
{
  double angle = PI / ( 2. * m_visuNodeNbPerQar ) ;
  double angleY = 0., local_radius = 0.;
  int k, i, ptsPerlevel = 4 * m_visuNodeNbPerQar,
  	Top_number = int(firstpoint_number) + ptsPerlevel * m_visuNodeNbPerQar;
  Point3 p, pp;
  
  // Top sphere 
  // Regular points on the surface
  for ( k = m_visuNodeNbPerQar - 1; k < 2*m_visuNodeNbPerQar-1 ; ++k ) 
  {  
    angleY = - PI / 2. + ( k + 1 ) * angle;
    local_radius = m_radius * cos( angleY );
    p[Y] = m_radius * sin( angleY ) + 0.5 * m_height;
    for ( i = 0; i < ptsPerlevel ; ++i )
    {
      p[X] = local_radius * cos( i * angle );
      p[Z] = local_radius * sin( i * angle );
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
  p[Z] = 0.;	
  // Top point
  p[Y] = m_radius + 0.5 * m_height;
  pp = transform( p );
  f << "v " << GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
			pp[X] ) << " " << 
	GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
			pp[Y] ) << " " <<
	GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
			pp[Z] ) << " " << endl;  
  
  // Bottom sphere	 	  
  // Regular points on the surface
  for ( k = m_visuNodeNbPerQar - 1; k >=0 ; --k ) 
  {  
    angleY = - PI / 2. + ( k + 1 ) * angle;
    local_radius = m_radius * cos( angleY );
    p[Y] = m_radius * sin( angleY ) - 0.5 * m_height;
    for ( i = 0; i < ptsPerlevel ; ++i )
    {
      p[X] = local_radius * cos( i * angle );
      p[Z] = local_radius * sin( i * angle );
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
  p[Z] = 0.;
  // Bottom point
  p[Y] = - m_radius - 0.5 * m_height;
  pp = transform( p );
  f << "v " << GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
			pp[X] ) << " " << 
	GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
			pp[Y] ) << " " <<
	GrainsExec::doubleToString( ios::scientific, FORMAT10DIGITS,
			pp[Z] ) << " " << endl;    
  
  // Cylinder: no additional point needed 
  
  // Faces  
  // Top sphere
  // Square faces
  for ( k = 0; k < m_visuNodeNbPerQar-1 ; ++k ) 
  {  
    for ( i = 0; i < ptsPerlevel-1 ; ++i )
    {
      f << "f " << firstpoint_number + k * ptsPerlevel + i << " "
      	<< firstpoint_number + k * ptsPerlevel + i + 1 << " " 
      	<< firstpoint_number + ( k + 1) * ptsPerlevel + i + 1 << " " 
      	<< firstpoint_number + ( k + 1 ) * ptsPerlevel + i << endl;		
    }
    f << "f " << firstpoint_number + k * ptsPerlevel + ptsPerlevel - 1 << " "
    	<< firstpoint_number + k * ptsPerlevel << " "
    	<< firstpoint_number + ( k + 1 ) * ptsPerlevel << " "
    	<< firstpoint_number + ( k + 1 ) * ptsPerlevel 
		+ ptsPerlevel - 1 << endl; 
  }   
  
  // Triangular faces  
  for ( i = 0; i < ptsPerlevel-1 ; ++i )
  {
    f << "f " << firstpoint_number
    		+ ( m_visuNodeNbPerQar - 1 ) * ptsPerlevel + i << " " 
    	<< firstpoint_number
    		+ ( m_visuNodeNbPerQar - 1 ) * ptsPerlevel + i + 1 << " " 
    	<< Top_number << endl;
  }
  f << "f " << firstpoint_number +  m_visuNodeNbPerQar * ptsPerlevel - 1 
  	<< " " << firstpoint_number 
		+ ( m_visuNodeNbPerQar - 1 ) * ptsPerlevel << " " 
  	<< Top_number  << endl; 
	
  // Bottom sphere
  int shift = int(firstpoint_number) + ptsPerlevel * m_visuNodeNbPerQar + 1 ;  
  int Bottom_number = shift + ptsPerlevel * m_visuNodeNbPerQar; 
  // Square faces
  for ( k = 0; k < m_visuNodeNbPerQar-1 ; ++k ) 
  {  
    for ( i = 0; i < ptsPerlevel-1 ; ++i )
    {
      f << "f " << shift + k * ptsPerlevel + i << " " 
      	<< shift + k * ptsPerlevel + i + 1 << " " 
      	<< shift + ( k + 1) * ptsPerlevel + i + 1 << " " 	
      	<< shift + ( k + 1 ) * ptsPerlevel + i << endl;
    }
    f << "f " << shift + k * ptsPerlevel + ptsPerlevel - 1 << " "
    	<< shift + k * ptsPerlevel << " "
    	<< shift + ( k + 1 ) * ptsPerlevel << " "
    	<< shift + ( k + 1 ) * ptsPerlevel + ptsPerlevel - 1 << endl;    
  }

  // Triangular faces 
  for ( i = 0; i < ptsPerlevel-1 ; ++i )
  {
    f << "f " << shift + ( m_visuNodeNbPerQar - 1 ) * ptsPerlevel + i << " "
    	<< shift + ( m_visuNodeNbPerQar - 1 ) * ptsPerlevel + i + 1 << " "
    	<< Bottom_number << endl; 
  }
  f << "f " << shift +  m_visuNodeNbPerQar * ptsPerlevel - 1 << " "
  	<< shift + ( m_visuNodeNbPerQar - 1 ) * ptsPerlevel << " "
  	<< Bottom_number << endl;	
	
  // Cylinder	
  for ( i = 0; i < ptsPerlevel-1 ; ++i )
  {
    f << "f " << firstpoint_number + i << " "
    	<< firstpoint_number + i + 1 << " "
    	<< shift + i + 1 << " "
    	<< shift + i << endl;
  }
  f << "f " << firstpoint_number + ptsPerlevel - 1 << " "
  	<< firstpoint_number << " "
  	<< shift << " "
  	<< shift + ptsPerlevel - 1 << endl;
	
  firstpoint_number += size_t(2 * ( ptsPerlevel * m_visuNodeNbPerQar + 1 )) ;
}
