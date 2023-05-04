#include "BBox.hh"
#include <iostream>
using namespace std;


// --------------------------------------------------------------------
// Default constructor
BBox::BBox()
{}




// --------------------------------------------------------------------
// Constructors with 2 corners as inputs, the 1st point with the
// lowest coordinates and the 2nd point with the largest coordinates
BBox::BBox( Point3 const& min, Point3 const& max )
{
  setValue( min, max );
}




// --------------------------------------------------------------------
// Copy constructor
BBox::BBox( BBox const& bbox_ )
{
  m_center = bbox_.m_center;
  m_extent = bbox_.m_extent;
}




// --------------------------------------------------------------------
// Destructeur
BBox::~BBox()
{}




// --------------------------------------------------------------------
// Operateur d'affectation
// A.WACHS - Aout.2009 - Creation
BBox& BBox::operator= (const BBox& rhs)
{
  if ( &rhs != this )
  {
    m_center = rhs.m_center;
    m_extent = rhs.m_extent;
  }
  return ( *this );
}




// ----------------------------------------------------------------------------
// Sets the bounding box to the intersection of 2 bounding boxes. If
// there is no intersection, leaves the bounding box unchanged
void BBox::closest( BBox const& a, BBox const& b )
{
  if ( intersect(a,b) )
  {
    Point3 lower( max( a.getLower(X), b.getLower(X) ),
	max( a.getLower(Y), b.getLower(Y) ),
	max( a.getLower(Z), b.getLower(Z) ) );
    Point3 upper( min( a.getUpper(X), b.getUpper(X) ),
	min( a.getUpper(Y), b.getUpper(Y) ),
	min( a.getUpper(Z), b.getUpper(Z) ) );
    setValue( lower, upper );
  }
}




// --------------------------------------------------------------------
// Sets the bounding box to the union of the 2 bounding boxes a and b
void BBox::enclose( BBox const& a, BBox const& b )
{
  Point3 lower( min(a.getLower(X), b.getLower(X) ),
	min( a.getLower(Y), b.getLower(Y) ),
	min( a.getLower(Z), b.getLower(Z) ) );
  Point3 upper( max( a.getUpper(X), b.getUpper(X) ),
	max( a.getUpper(Y), b.getUpper(Y) ),
	max( a.getUpper(Z), b.getUpper(Z) ) );
  setValue( lower, upper );
}




// --------------------------------------------------------------------
// Returns the bounding box center
Point3 const& BBox::getCenter() const
{
  return ( m_center );
}




// --------------------------------------------------------------------
// Returns the bounding box center
Vector3 const& BBox::getExtent() const
{
  return ( m_extent );
}




// --------------------------------------------------------------------
// Returns the ith minimum coordinate
double BBox::getLower( int i ) const
{
  return ( m_center[i] - m_extent[i] );
}




// --------------------------------------------------------------------
// Returns the ith maximum coordinate
double BBox::getUpper( int i ) const
{
  return ( m_center[i] + m_extent[i] );
}




// --------------------------------------------------------------------
// Extends the bounding box to a point p if p is outside the box
void BBox::include( Point3 const& p )
{
  Point3 lower( min( getLower(X), p[X] ),
	min( getLower(Y), p[Y] ),
	min( getLower(Z), p[Z] ) );
  Point3 upper( max( getUpper(X), p[X] ),
	max( getUpper(Y), p[Y] ),
	max( getUpper(Z), p[Z] ) );
  setValue( lower, upper );
}




// ----------------------------------------------------------------------
// Sets the bounding box to the union of itself and another bounding box
void BBox::include( BBox const& b )
{
  enclose( *this, b );
}




// --------------------------------------------------------------------
// Returns whether a cubic box defined by its center and its half
// edge length intersects
bool BBox::InZone( Point3 const* p, double halfEdgeLength ) const
{
  return ( fabs( m_center[X] - (*p)[X] ) <= m_extent[X] + halfEdgeLength &&
    fabs( m_center[Y] - (*p)[Y] ) <= m_extent[Y] + halfEdgeLength &&
    fabs( m_center[Z] - (*p)[Z] ) <= m_extent[Z] + halfEdgeLength );
}




// --------------------------------------------------------------------
// Returns whether a box defined by its center and its half edge lengths
// intersects
bool BBox::InZone( Point3 const* p, double halfEdgeLength_X,
    	double halfEdgeLength_Y, double halfEdgeLength_Z ) const
{
  return ( fabs( m_center[X] - (*p)[X] ) <= m_extent[X] + halfEdgeLength_X &&
    fabs( m_center[Y] - (*p)[Y] ) <= m_extent[Y] + halfEdgeLength_Y &&
    fabs( m_center[Z] - (*p)[Z] ) <= m_extent[Z] + halfEdgeLength_Z );
}




// --------------------------------------------------------------------
// Returns the direction of longest edge
int BBox::longestAxis() const
{
  return ( m_extent.closestAxis() );
}




// --------------------------------------------------------------------
// Sets the bounding box center
void BBox::setCenter( Point3 const& p )
{
  m_center = p;
}




// ----------------------------------------------------------------------------
// Sets the bounding box to an empty bounding box. This is done by
// assigning minus infinity extensions
void BBox::setEmpty()
{
  m_center.setValue( 0, 0, 0 );
  m_extent.setValue( -INFINITY, -INFINITY, -INFINITY );
}




// --------------------------------------------------------------------
// Sets the bounding boxes half lengths
void BBox::setExtent( Vector3 const& v )
{
  m_extent = v;
}




// --------------------------------------------------------------------
// Sets the box dimensions using the point with the
// lowest coordinates and the point with the largest coordinates
void BBox::setValue( Point3 const& min, Point3 const& max )
{
  m_extent = (max - min) / 2.;
  m_center = min + m_extent;
}




// ----------------------------------------------------------------------------
// Returns the largest half length of the bounding box
double BBox::size() const
{
  return ( max( max( m_extent[X], m_extent[Y] ), m_extent[Z] ) );
}




// ----------------------------------------------------------------------------
// Returns whether the 2 bounding boxes a and b intersect
bool intersect( BBox const& a, BBox const& b )
{
  return (
    fabs( a.m_center[X] - b.m_center[X] ) <= a.m_extent[X] + b.m_extent[X] &&
    fabs( a.m_center[Y] - b.m_center[Y] ) <= a.m_extent[Y] + b.m_extent[Y] &&
    fabs( a.m_center[Z] - b.m_center[Z] ) <= a.m_extent[Z] + b.m_extent[Z] );
}





// ----------------------------------------------------------------------------
// Debugging method
void BBox::debug( char const* s ) const
{
  cerr << s;
  cerr << "BBox Center " << m_center << "     Extent " << m_extent;
}




// ----------------------------------------------------------------------------
// Output operator
ostream& operator << ( ostream& f, BBox const& B )
{
  f << "BBox: Center = " << B.m_center;
  f << "      Extent = " << B.m_extent;
  return ( f );
}




// ----------------------------------------------------------------------------
// Returns whether the bounding box fully contains the other bounding box
bool BBox::fullyContain( BBox const& a )
{
  return (
    a.m_center[X] - a.m_extent[X] >= m_center[X] - m_extent[X] &&
    a.m_center[X] + a.m_extent[X] <= m_center[X] + m_extent[X] &&
    a.m_center[Y] - a.m_extent[Y] >= m_center[Y] - m_extent[Y] &&
    a.m_center[Y] + a.m_extent[Y] <= m_center[Y] + m_extent[Y] &&
    a.m_center[Z] - a.m_extent[Z] >= m_center[Z] - m_extent[Z] &&
    a.m_center[Z] + a.m_extent[Z] <= m_center[Z] + m_extent[Z] ) ;
}