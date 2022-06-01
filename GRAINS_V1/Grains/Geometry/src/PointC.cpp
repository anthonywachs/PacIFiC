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
// Returns whether a point lies inside the point (returns false
// by convention)
bool PointC::isIn( Point3 const& pt ) const
{
  return ( false );
}  
