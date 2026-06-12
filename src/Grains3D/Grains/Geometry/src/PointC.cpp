#include "OBB.hh"
#include "OBC.hh"
#include "PointC.hh"



// ----------------------------------------------------------------------------
// Default constructor
PointC::PointC() 
{}




// ----------------------------------------------------------------------------
// Constructor with an input stream
PointC::PointC( istream& fileIn )
{
  readShape( fileIn );
}




// ----------------------------------------------------------------------------
// Destructor
PointC::~PointC() 
{}




// ----------------------------------------------------------------------
// Returns the convex type
ConvexType PointC::getConvexType() const 
{
  return ( POINT );
}




// ----------------------------------------------------------------------------
// Computes the inertia tensor and the inverse of the inertia tensor
bool PointC::BuildInertia( double* inertia, double* inertia_1 ) const
{
  inertia[0] = inertia[1] = inertia[2] = inertia[3] = inertia[4]
	= inertia[5] = 0.;
  inertia_1[0] = inertia_1[1] = inertia_1[2] = inertia_1[3] = inertia_1[4]
	= inertia_1[5] = 0.;
  return ( true );
}




// ----------------------------------------------------------------------------
// Returns the circumscribed radius of the reference disc,
// i.e., without applying any transformation
double PointC::computeCircumscribedRadius() const 
{
  return ( 0. );
}




// ----------------------------------------------------------------------------
// Returns a clone of the point
Convex* PointC::clone() const 
{
  return ( new PointC() );
}


// ----------------------------------------------------------------------------
// Returns the point volume, 0 by convention
double PointC::getVolume() const
{
  return ( 0. );
}




// ----------------------------------------------------------------------------
// Point support function, returns the support point P, i.e. the
// point on the surface of the Point that satisfies max(P.v)
Point3 PointC::support( Vector3 const& v ) const 
{
  return ( Point3(0., 0., 0.) );
}




// ----------------------------------------------------------------------------
// Output operator
void PointC::writeShape( ostream& fileOut ) const 
{
  fileOut << "*PointC *END"; ;
}



  
// ----------------------------------------------------------------------------
// Input operator
void PointC::readShape( istream& fileIn ) 
{}




// ----------------------------------------------------------------------------
// Returns the number of points to write the PointC in a Paraview format
int PointC::numberOfPoints_PARAVIEW() const
{
  return ( 1 );
}




// ----------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the PointC 
// in a Paraview format
int PointC::numberOfCells_PARAVIEW() const
{
  return ( 1 );  
}




// ----------------------------------------------------------------------------
// Writes a list of points describing the PointC in a Paraview format
void PointC::write_polygonsPts_PARAVIEW( ostream& f, 
  	Transform const& transform, Vector3 const* translation ) const
{
  Point3 pp, p;
  
  pp = transform( p );
  if ( translation ) pp += *translation;    
  f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
}




// ----------------------------------------------------------------------------
// Returns a list of points describing the PointC in a Paraview format 
list<Point3> PointC::get_polygonsPts_PARAVIEW( Transform const& transform,
  	Vector3 const* translation ) const
{
  list<Point3> ParaviewPoints;
  Point3 pp, p;
  
  pp = transform( p );
  if ( translation ) pp += *translation;    
  ParaviewPoints.push_back( pp );
  
  return ( ParaviewPoints ); 
}




// ----------------------------------------------------------------------------
// Writes the PointC in a Paraview format
void PointC::write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const
{
  connectivity.push_back( firstpoint_globalnumber );
  last_offset += 1;  
  offsets.push_back( last_offset );
  cellstype.push_back( 1 );
  
  firstpoint_globalnumber += 1;
}




// ----------------------------------------------------------------------------
// Returns whether a point lies inside the point (returns false
// by convention)
bool PointC::isIn( Point3 const& pt ) const
{
  return ( false );
}  




// ----------------------------------------------------------------------------
// Performs advanced comparison of the two points and returns whether 
// they match
bool PointC::equalType_level2( Convex const* other ) const
{
  return ( true );
} 




// ----------------------------------------------------------------------------
// Returns the bounding volume to point
BVolume* PointC::computeBVolume( unsigned int type ) const
{
  BVolume* bvol = NULL;
  if ( type == 1 ) // OBB
    bvol = new OBB( Vector3( 0., 0., 0. ), Matrix() );
  else if ( type == 2 ) // OBC
  {
    bvol = new OBC( 0., 
                    0.,
                    Vector3( 0., 0., 1. ) );
  }
  return( bvol );
}
