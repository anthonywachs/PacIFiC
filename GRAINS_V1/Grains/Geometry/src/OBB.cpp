#include "OBB.hh"


// --------------------------------------------------------------------
// Default constructor
OBB::OBB()
{}




// --------------------------------------------------------------------
// Constructor with half-lengths extent and initial orientation ori
OBB::OBB( Vector3 const& extent, Matrix const& ori )
{
  m_extent = extent;
  m_initOrientation = ori;
}




// --------------------------------------------------------------------
// Copy constructor
OBB::OBB( OBB const& obb_ )
{
  m_extent = obb_.m_extent;
  m_initOrientation = obb_.m_initOrientation;
}




// --------------------------------------------------------------------
// Destructor
OBB::~OBB()
{}




// --------------------------------------------------------------------
// Returns the OBC radius
BVolumeType OBB::getBVolumeType() const
{
  return ( typeOBB );
}




// --------------------------------------------------------------------
// Returns a clone of the OBB
BVolume* OBB::clone() const
{
  return ( new OBB( m_extent, m_initOrientation ) );
}




// --------------------------------------------------------------------
// Returns the OBB half-lengths
Vector3 const& OBB::getExtent() const
{
  return ( m_extent );
}




// --------------------------------------------------------------------
// Returns the OBB initial orientation
Matrix const& OBB::getInitOrientation() const
{
  return ( m_initOrientation );
}




// --------------------------------------------------------------------
// Sets the OBB radius
void OBB::setExtent( Vector3 const& extent )
{
  m_extent = extent;
}




// --------------------------------------------------------------------
// Sets the OBB initial orientation
void OBB::setInitOrientation( Matrix const& ori )
{
  m_initOrientation = ori;
}




// ----------------------------------------------------------------------------
// Output operator
void OBB::writeShape ( ostream& fileOut ) const
{
  fileOut << "*OBB " << m_extent << " " << m_initOrientation;
}




// ----------------------------------------------------------------------------
#define TESTCASE1( i ) \
  ( fabs( cen[i] ) > \
  ( a[i] + b[0]*fabs(ori[0][i]) + b[1]*fabs(ori[1][i]) + b[2]*fabs(ori[2][i]) ) )

#define TESTCASE2( i ) \
  ( fabs( cen[0]*ori[i][0] + cen[1]*ori[i][1] + cen[2]*ori[i][2] ) > \
  ( b[i] + a[0]*fabs(ori[i][0]) + a[1]*fabs(ori[i][1]) + a[2]*fabs(ori[i][2]) ) )

// #define TESTCASE3(i, j) \
//   ( fabs( cen[(j+2)%3]*ori[i][(j+1)%3] - cen[(j+1)%3]*ori[i][(j+2)%3] ) > \
//   ( a[(i+1)%3]*fabs(ori[(i+2)%3][j]) + a[(i+2)%3]*fabs(ori[(i+1)%3][j]) + \
//     b[(j+1)%3]*fabs(ori[i][(j+2)%3]) + b[(j+2)%3]*fabs(ori[i][(j+1)%3]) ) )

#define TESTCASE3(i, j) \
  ( fabs( cen[(j+2)%3]*ori[i][(j+1)%3] - cen[(j+1)%3]*ori[i][(j+2)%3] ) > \
  ( a[(i+1)%3]*oriAbs[(i+2)%3][j] + a[(i+2)%3]*oriAbs[(i+1)%3][j] + \
    b[(j+1)%3]*oriAbs[i][(j+2)%3] + b[(j+2)%3]*oriAbs[i][(j+1)%3] ) )

// Returns whether two OBBs are in contact
bool isContactBVolume( OBB const& obbA,
                       OBB const& obbB,
                       Transform const& a2w,
                       Transform const& b2w )
{
  Vector3 const& a = obbA.getExtent();
  Vector3 const& b = obbB.getExtent();
  Point3 const& cen = transpose( a2w.getBasis() ) * 
                      ( *( b2w.getOrigin() ) - *( a2w.getOrigin() ) );
  Matrix const& ori = transpose( b2w.getBasis() ) * a2w.getBasis();
  Matrix const oriAbs( fabs(ori[0][0]), fabs(ori[0][1]), fabs(ori[0][2]),
                       fabs(ori[1][0]), fabs(ori[1][1]), fabs(ori[1][2]),
                       fabs(ori[2][0]), fabs(ori[2][1]), fabs(ori[2][2]) );

  // CASE 1: ( three of them )
  if TESTCASE1( 0 ) return ( false );
  if TESTCASE1( 1 ) return ( false );
  if TESTCASE1( 2 ) return ( false );

  // CASE 2: ( three of them )
  if TESTCASE2( 0 ) return ( false );
  if TESTCASE2( 1 ) return ( false );
  if TESTCASE2( 2 ) return ( false );

  // CASE 3: ( nine of them )
  if TESTCASE3( 0, 0 ) return ( false );
  if TESTCASE3( 1, 0 ) return ( false );
  if TESTCASE3( 2, 0 ) return ( false );
  if TESTCASE3( 0, 1 ) return ( false );
  if TESTCASE3( 1, 1 ) return ( false );
  if TESTCASE3( 2, 1 ) return ( false );
  if TESTCASE3( 0, 2 ) return ( false );
  if TESTCASE3( 1, 2 ) return ( false );
  if TESTCASE3( 2, 2 ) return ( false );
  
  return ( true );
}
#undef TESTCASE1
#undef TESTCASE2
#undef TESTCASE3
