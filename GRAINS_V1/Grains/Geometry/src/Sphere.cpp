#include "Sphere.hh"
#include "BVolume.hh"
#include "OBB.hh"
#include "OBC.hh"
#include "sstream"
#include <math.h>


int Sphere::m_visuNodeNbPerQar = 8;

// ----------------------------------------------------------------------------
// Constructor with radius as input parameter
Sphere::Sphere( double r ) 
  : m_radius( r )
{}




// ----------------------------------------------------------------------------
// Constructor with an input stream
Sphere::Sphere( istream& fileIn )
{
  readShape( fileIn );
}




// ----------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
Sphere::Sphere( DOMNode* root )
{
  m_radius = ReaderXML::getNodeAttr_Double(root, "Radius");
}




// ----------------------------------------------------------------------------
// Destructor
Sphere::~Sphere() 
{}




// ----------------------------------------------------------------------------
// Returns the convex type
ConvexType Sphere::getConvexType() const 
{
  return ( SPHERE );
}




// ----------------------------------------------------------------------------
// Returns whether the convex shape is a sphere
bool Sphere::isSphere() const 
{
  return ( true );
}




// ----------------------------------------------------------------------------
// Computes the inertia tensor and the inverse of the inertia tensor
bool Sphere::BuildInertia( double* inertia, double* inertia_1 ) const
{
  inertia[1] = inertia[2] = inertia[4] = 0.0;
  inertia[5] = inertia[3] = inertia[0] = 8.0 * PI * pow( m_radius, 5. ) / 15.0;
  
  inertia_1[1] = inertia_1[2] = inertia_1[4] = 0.0;
  inertia_1[5] = inertia_1[3] = inertia_1[0] = 1.0 / inertia[0];

  return ( true );
}




// ----------------------------------------------------------------------------
// Returns the circumscribed radius of the reference convex shape,
// i.e., without applying any transformation
double Sphere::computeCircumscribedRadius() const 
{
  return ( m_radius );
}




// ----------------------------------------------------------------------------
// Returns a clone of the sphere
Convex* Sphere::clone() const
{
  return ( new Sphere( m_radius ) );
}




// ----------------------------------------------------------------------------
// Sphere support function, returns the support point P, i.e. the
// point on the surface of the sphere that satisfies max(P.v)
Point3 Sphere::support( Vector3 const& v ) const 
{
  double s = Norm( v );
  if (s > EPSILON) 
  {
    double r = m_radius / s;
    return ( Point3(v[X] * r, v[Y] * r, v[Z] * r) );
  } 
  else 
    return ( Point3() );    
}




// ----------------------------------------------------------------------------
// Returns a vector of points describing the envelope of the sphere
vector<Point3> Sphere::getEnvelope() const
{
  vector<Point3> envelope;
  Point3 point( 0., 0., 0. );
  envelope.push_back( point );

  return ( envelope );
}




// ----------------------------------------------------------------------------
// Returns the number of vertices/corners, here 1 as a convention
int Sphere::getNbCorners() const
{
  return ( 1 );
}




// ----------------------------------------------------------------------------
// Returns a pointer to a 2D array describing the relationship between the face 
// indices and the point indices. Returns a null pointer as a convention
vector< vector<int> > const* Sphere::getFaces() const
{
  vector< vector<int> > const* allFaces = NULL;
  return ( allFaces );
}




// ----------------------------------------------------------------------------
// Returns the sphere volume
double Sphere::getVolume() const
{
  return ( 4.0 * PI * pow( m_radius, 3. ) / 3.0 );
}




// ----------------------------------------------------------------------------
// Output operator
void Sphere::writeShape( ostream &fileOut ) const 
{
  fileOut << "*Sphere " << m_radius << " *END"; 
}




// ----------------------------------------------------------------------------
// Input operator
void Sphere::readShape( istream &fileIn ) 
{
  fileIn >> m_radius;
}




// ----------------------------------------------------------------------------
// Returns the number of points to write the sphere in a Paraview format
int Sphere::numberOfPoints_PARAVIEW() const 
{ 
  return ( 4 * m_visuNodeNbPerQar * ( 2 * m_visuNodeNbPerQar - 1 ) + 3 );
}




// ----------------------------------------------------------------------------
// Writes a list of points describing the sphere in a Paraview format
void Sphere::write_polygonsPts_PARAVIEW( ostream& f, 
  	Transform const& transform, Vector3 const* translation ) const
{
  double angle = PI / ( 2. * m_visuNodeNbPerQar ) ;
  double angleZ = 0., local_radius = 0.;
  int k, i, ptsPerlevel = 4 * m_visuNodeNbPerQar;
  Point3 pp, pptrans;
  
  // Regular points on the surface
  for ( k = 0; k < 2*m_visuNodeNbPerQar-1 ; ++k ) 
  {  
    angleZ = - PI / 2. + ( k + 1 ) * angle;
    local_radius = m_radius * cos( angleZ );
    pp[Z] = m_radius * sin( angleZ );
    for ( i = 0; i < ptsPerlevel ; ++i )
    {
      pp[X] = local_radius * cos( i * angle );
      pp[Y] = local_radius * sin( i * angle );
      pptrans = transform( pp );
      if ( translation ) pptrans += *translation;
      f << pptrans[X] << " " << pptrans[Y] << " " << pptrans[Z] << endl;
    }
  }
   
  pp[X] = 0.;
  pp[Y] = 0.;
  // Bottom point
  pp[Z] = - m_radius;
  pptrans = transform( pp );
  if ( translation ) pptrans += *translation;
  f << pptrans[X] << " " << pptrans[Y] << " " << pptrans[Z] << endl;
	
  // Top point
  pp[Z] = m_radius;
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
// Writes the sphere in a STL format
void Sphere::write_convex_STL( ostream& f, Transform const& transform ) const
{
  int visuNodeNbPerQar_STL = 16;
  double angle = PI / ( 2. * visuNodeNbPerQar_STL ) ;
  double angleZ = 0., local_radius = 0.;
  int k, i, ptsPerlevel = 4 * visuNodeNbPerQar_STL, 
  	nbLevels = 2 * visuNodeNbPerQar_STL - 1;
  Point3 pp, pptrans, ppTop, ppBottom, GC;
  vector<Point3> work( ptsPerlevel, pp );
  vector< vector<Point3> > STLPoints( nbLevels, work );  
  
  // Regular points on the surface
  for ( k = 0; k < nbLevels ; ++k ) 
  {  
    angleZ = -PI / 2. + ( k + 1 ) * angle;
    local_radius = m_radius * cos( angleZ );
    pp[Z] = m_radius * sin( angleZ );
    for ( i = 0; i < ptsPerlevel ; ++i )
    {
      pp[X] = local_radius * cos( i * angle );
      pp[Y] = local_radius * sin( i * angle );
      pptrans = transform( pp );
      STLPoints[k][i] = pptrans;
    }
  }
   
  pp[X] = 0.;
  pp[Y] = 0.;
  // Bottom point
  pp[Z] = - m_radius;
  ppBottom = transform( pp );
	
  // Top point
  pp[Z] = m_radius;
  ppTop = transform( pp );
	
  // Gravity center
  pp[Z] = 0.;
  GC = transform( pp );  
  
  
  // Writing facets
  // Regular facets
  for ( k = 0; k < nbLevels-1 ; ++k ) 
  {
    for ( i = 0; i < ptsPerlevel-1 ; ++i )
    {
      // Bottom triangle
      write_STLfacet_sphere( f, GC, STLPoints[k][i], STLPoints[k][i+1],
      	STLPoints[k+1][i] );
            
      // Top triangle
      write_STLfacet_sphere( f, GC, STLPoints[k][i+1], STLPoints[k+1][i+1],
      	STLPoints[k+1][i] );
    }
    
    // Last bottom triangle
    write_STLfacet_sphere( f, GC, STLPoints[k][ptsPerlevel-1], STLPoints[k][0],
      	STLPoints[k+1][ptsPerlevel-1] ); 
	
    // Last top triangle
    write_STLfacet_sphere( f, GC, STLPoints[k][0], STLPoints[k+1][0],
      	STLPoints[k+1][ptsPerlevel-1] );	 
  }

  // Bottom facets
  for ( i = 0; i < ptsPerlevel-1 ; ++i )
    write_STLfacet_sphere( f, GC, ppBottom, STLPoints[0][i], 
    	STLPoints[0][i+1] );
  write_STLfacet_sphere( f, GC, ppBottom, STLPoints[0][ptsPerlevel-1],
  	STLPoints[0][0] );  

  // Top facets
  for ( i = 0; i < ptsPerlevel-1 ; ++i )
    write_STLfacet_sphere( f, GC, ppTop, STLPoints[nbLevels-1][i],
    	STLPoints[nbLevels-1][i+1] );
  write_STLfacet_sphere( f, GC, ppTop, STLPoints[nbLevels-1][ptsPerlevel-1],
    	STLPoints[nbLevels-1][0] );	  
}




// ----------------------------------------------------------------------------
// Writes a triangular facet in the STL format based on its 3
// vertices and the center of mass coordinates
void Sphere::write_STLfacet_sphere( ostream &f, Point3 const& GC,
  	Point3 const& pp1,
  	Point3 const& pp2,	
  	Point3 const& pp3 ) const
{  
  Point3 triangleGC = pp1 + ( 2. / 3. ) * ( 0.5 * ( pp2 + pp3 ) - pp1 );
  Vector3 outward_normal = triangleGC - GC;
  outward_normal.normalize();
  f << "  facet normal " << outward_normal[X]
      	<< " " << outward_normal[Y]
      	<< " " << outward_normal[Z] << endl;
  f << "    outer loop" << endl;
  f << "      vertex " << pp1[X] << " " << pp1[Y] << " " << pp1[Z] << endl;
  f << "      vertex " << pp2[X] << " " << pp2[Y] << " " << pp2[Z] << endl;
  f << "      vertex " << pp3[X] << " " << pp3[Y] << " " << pp3[Z] << endl;
  f << "    endloop" << endl;
  f << "  endfacet" << endl;     
}




// ----------------------------------------------------------------------------
// Returns a list of points describing the sphere in a Paraview format
list<Point3> Sphere::get_polygonsPts_PARAVIEW( Transform const& transform,
	Vector3 const* translation ) const
{
  list<Point3> ParaviewPoints;
  double angle = PI / ( 2. * m_visuNodeNbPerQar ) ;
  double angleZ = 0., local_radius = 0.;
  int k, i, ptsPerlevel =  4 * m_visuNodeNbPerQar;
  Point3 pp, pptrans;
  
  // Regular points on the surface
  for ( k = 0; k < 2*m_visuNodeNbPerQar-1 ; ++k ) 
  {  
    angleZ = -PI / 2. + ( k + 1 ) * angle;
    local_radius = m_radius * cos( angleZ );
    pp[Z] = m_radius * sin( angleZ );
    for ( i = 0; i < ptsPerlevel ; ++i )
    {
      pp[X] = local_radius * cos( i * angle );
      pp[Y] = local_radius * sin( i * angle );
      pptrans = transform( pp );
      if ( translation ) pptrans += *translation;
      ParaviewPoints.push_back( pptrans );
    }
  }
   
  pp[X] = 0.;
  pp[Y] = 0.;
  // Bottom point
  pp[Z] = - m_radius;
  pptrans = transform( pp );
  if ( translation ) pptrans += *translation;
  ParaviewPoints.push_back( pptrans );
	
  // Top point
  pp[Z] = m_radius;
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
// Returns the number of elementary polytopes to write the sphere in a 
// Paraview format
int Sphere::numberOfCells_PARAVIEW() const
{
  return ( 8 * m_visuNodeNbPerQar * m_visuNodeNbPerQar ); 
}




// ----------------------------------------------------------------------------
// Writes the sphere in a Paraview format
void Sphere::write_polygonsStr_PARAVIEW( list<int>& connectivity,
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

  // Bottom cells: tetraedron
  for ( i = 0; i < ptsPerlevel-1 ; ++i )
  {
    connectivity.push_back( firstpoint_globalnumber + i ); 
    connectivity.push_back( firstpoint_globalnumber + i + 1 );
    connectivity.push_back( firstpoint_globalnumber + Bottom_number );	
    connectivity.push_back( firstpoint_globalnumber + GC_number );
    last_offset += 4;
    offsets.push_back( last_offset );
    cellstype.push_back( 10 );   
  }
  connectivity.push_back( firstpoint_globalnumber + ptsPerlevel - 1 );
  connectivity.push_back( firstpoint_globalnumber );
  connectivity.push_back( firstpoint_globalnumber + Bottom_number );	
  connectivity.push_back( firstpoint_globalnumber + GC_number );
  last_offset += 4;
  offsets.push_back( last_offset );
  cellstype.push_back( 10 );  
  
  // Top cells: tetraedron  
  for ( i = 0; i < ptsPerlevel-1 ; ++i )
  {
    connectivity.push_back( firstpoint_globalnumber
    	+ ( 2 * m_visuNodeNbPerQar - 2 ) * ptsPerlevel + i );
    connectivity.push_back( firstpoint_globalnumber
    	+ ( 2 * m_visuNodeNbPerQar - 2 ) * ptsPerlevel + i + 1 );
    connectivity.push_back( firstpoint_globalnumber + Top_number );	
    connectivity.push_back( firstpoint_globalnumber + GC_number );
    last_offset += 4;
    offsets.push_back( last_offset );
    cellstype.push_back( 10 ); 
  }
  connectivity.push_back( firstpoint_globalnumber
  	+ ( 2 * m_visuNodeNbPerQar - 1 ) * ptsPerlevel - 1 );
  connectivity.push_back( firstpoint_globalnumber
  	+ ( 2 * m_visuNodeNbPerQar - 2 ) * ptsPerlevel );
  connectivity.push_back( firstpoint_globalnumber + Top_number );	
  connectivity.push_back( firstpoint_globalnumber + GC_number );
  last_offset += 4;
  offsets.push_back( last_offset );
  cellstype.push_back( 10 ); 

  firstpoint_globalnumber += 4 * m_visuNodeNbPerQar * 
  	( 2 * m_visuNodeNbPerQar - 1 ) + 3 ;
}




// ----------------------------------------------------------------------------
// Returns an orientation vector describing the convex shape angular
// position
Vector3 Sphere::computeOrientationVector( Transform const* transform ) const
{
  Point3 pp( 0., m_radius, 0. );
  Point3 pptrans = (*transform)( pp ); 

  return ( pptrans - *transform->getOrigin() );
} 




// ----------------------------------------------------------------------------
// Sets the number of point per quarter of the equator line for 
// Paraview post-processing, i.e., controls the number of facets in the sphere
// reconstruction in Paraview
void Sphere::SetvisuNodeNbPerQar( int nbpts )
{
  m_visuNodeNbPerQar = nbpts;
}




// ----------------------------------------------------------------------------
// Returns whether a point lies inside the sphere
bool Sphere::isIn( Point3 const& pt ) const
{
  return ( pt[X] * pt[X] + pt[Y] * pt[Y] +  pt[Z] * pt[Z] 
  	<= m_radius * m_radius );
}



 
// ----------------------------------------------------------------------------
// Returns the bounding volume to box
BVolume* Sphere::computeBVolume( unsigned int type ) const
{
  BVolume* bvol = NULL;
  if ( type == 1 ) // OBB
  {
    Vector3 ext( m_radius, m_radius, m_radius );    
    bvol = new OBB( ext, Matrix() );
  }
  else if ( type == 2 ) // OBC
  {
    Vector3 e( 0., 0., 1. );
    bvol = new OBC( m_radius, 2. * m_radius, e );
  }

  return( bvol );
}



 
// ----------------------------------------------------------------------------
// Returns the bounding volume to sphere
BVolume* Sphere::computeBVolume( unsigned int type ) const
{
  BVolume* bvol = NULL;
  if ( type == 1 ) // OBB
    bvol = new OBB( Vector3( m_radius, m_radius, m_radius ), Matrix() );
  else if ( type == 2 ) // OBC
    bvol = new OBC( m_radius, 2. * m_radius, Vector3( 1., 0., 0. ) );

  return( bvol );
}




// ----------------------------------------------------------------------------
// Performs advanced comparison of the two spheres and returns whether 
// they match
bool Sphere::equalType_level2( Convex const* other ) const
{
  // We know that other points to a Sphere, we dynamically cast it to actual 
  // type
  Sphere const* other_ = dynamic_cast<Sphere const*>(other);
  
  double lmin = min( m_radius, other_->m_radius );

  bool same = ( fabs( m_radius - other_->m_radius ) <  LOWEPS * lmin );
  
  return ( same );
} 
