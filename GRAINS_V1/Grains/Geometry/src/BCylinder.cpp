#include "BCylinder.hh"
#include <iostream>
using namespace std;


// --------------------------------------------------------------------
// Default constructor
BCylinder::BCylinder()
{}




// --------------------------------------------------------------------
// Constructor with radius r, height h, and axis e
BCylinder::BCylinder( double r, double h, Vector3 const& v )
{
  m_radius = r;
  m_height = h;
  m_axis = v;
}




// --------------------------------------------------------------------
// Copy constructor
BCylinder::BCylinder( BCylinder const& bcylinder_ )
{
  m_radius = bcylinder_.m_radius;
  m_height = bcylinder_.m_height;
  m_axis = bcylinder_.m_axis;
}




// --------------------------------------------------------------------
// Destructor
BCylinder::~BCylinder()
{}




// --------------------------------------------------------------------
// Returns the bounding cylinder radius
double BCylinder::getRadius() const
{
  return ( m_radius );
}




// --------------------------------------------------------------------
// Returns the bounding cylinder height
double BCylinder::getHeight() const
{
  return ( m_height );
}




// --------------------------------------------------------------------
// Returns the bounding cylinder axis
Vector3 const& BCylinder::getAxis() const
{
  return ( m_axis );
}




// --------------------------------------------------------------------
// Sets the bounding cylinder radius
void BCylinder::setRadius( double r )
{
  m_radius = r;
}




// --------------------------------------------------------------------
// Sets the bounding cylinder height
void BCylinder::setHeight( double h )
{
  m_height = h;
}




// --------------------------------------------------------------------
// Sets the bounding cylinder axis
void BCylinder::setAxis( Vector3 const& v )
{
  m_axis = v;
}




// ----------------------------------------------------------------------------
// Output operator
ostream& operator << ( ostream& f, BCylinder const& B )
{
  f << "BCylinder: Radius = " << B.getRadius();
  f << "      Height = " << B.getHeight();
  f << "      Axis = " << B.getAxis();
  return ( f );
}




// ----------------------------------------------------------------------------
// Returns the contact point of two cylinders
PointContact intersect( BCylinder const& a, BCylinder const& b,
                        Transform const& a2w, Transform const& b2w )
{
  // Variables
  double rA = a.getRadius();
  double rB = b.getRadius();
  double hA = a.getHeight();
  double hB = b.getHeight();
  Vector3 e_A = a2w.getBasis() * a.getAxis();
  Vector3 e_B = b2w.getBasis() * b.getAxis();

  Vector3 zAxis( 0., 0., 1. );
  PointContact ptCont = PointNoContact;

  // Relative positions - B w.r.t. A
  Matrix rotMatA = getRotationMatrix( e_A, zAxis );
  Point3 x_B2A = rotMatA * ( *( b2w.getOrigin() ) - *( a2w.getOrigin() ) );
  Vector3 e_B2A =  ( rotMatA * e_B ).normalized();
  e_B2A.round( EPSILON );
  // Relative positions - A w.r.t. B
  Matrix rotMatB = getRotationMatrix( e_B, zAxis );
  Point3 x_A2B = rotMatB * ( *( a2w.getOrigin() ) - *( b2w.getOrigin() ) );
  Vector3 e_A2B = ( rotMatB * e_A ).normalized();
  e_A2B.round( EPSILON );

  int iter = 0;
  while( ptCont.getOverlapDistance() >= 0 && iter < 10 )
  {
    iter++;
    BCylinderContactWrapper( rA, hA, rB, hB, e_B2A, x_B2A, e_A2B, x_A2B, iter,
                             ptCont );
  }

  switch ( iter )
  {
    case 1: case 2: case 4: case 6: case 7: case 9:
      ptCont.setContact( *a2w.getOrigin() +
                        Point3( transpose( rotMatA ) * ptCont.getContact() ) );
      ptCont.setOverlapVector( transpose( rotMatA ) * ptCont.getOverlapVector() );
      break;
    case 3: case 5: case 8:
      ptCont.setContact( *b2w.getOrigin() +
                        Point3( transpose( rotMatB ) * ptCont.getContact() ) );
      ptCont.setOverlapVector( transpose( rotMatB ) * ptCont.getOverlapVector() );
      break;
    default:
      break;
  }

  return ( ptCont );
}




// ----------------------------------------------------------------------------
// Wrapper function for cylinders contacts
void BCylinderContactWrapper( double rA, double hA, double rB, double hB,
                              Vector3 const& e_B2A, Point3 const& x_B2A,
                              Vector3 const& e_A2B, Point3 const& x_A2B,
                              int method, PointContact& ptCont )
{
  switch ( method )
  {
    case 1:
      F2FB2BParContact( rA, hA, rB, hB, e_B2A, x_B2A, ptCont );
      break;
    case 2:
      F2BContact( rA, hA, rB, hB, e_B2A, x_B2A, ptCont );
      break;
    case 3:
      F2BContact( rB, hB, rA, hA, e_A2B, x_A2B, ptCont );
      ptCont.setOverlapVector( - ptCont.getOverlapVector() );
      break;
    case 4:
      F2EContact( rA, hA, rB, hB, e_B2A, x_B2A, ptCont );
      break;
    case 5:
      F2EContact( rB, hB, rA, hA, e_A2B, x_A2B, ptCont );
      ptCont.setOverlapVector( - ptCont.getOverlapVector() );
      break;
    case 6:
      B2BSkewContact( rA, hA, rB, hB, e_B2A, x_B2A, ptCont );
      break;
    case 7:
      B2EContact( rA, hA, rB, hB, e_B2A, x_B2A, ptCont );
      break;
    case 8:
      B2EContact( rB, hB, rA, hA, e_A2B, x_A2B, ptCont );
      ptCont.setOverlapVector( - ptCont.getOverlapVector() );
      break;
    case 9:
      E2EContact( rA, hA, rB, hB, e_B2A, x_B2A, ptCont );
      break;
    default:
      break;
  }
}




// ----------------------------------------------------------------------------
#define NormXY( A ) ( A[X]*A[X] + A[Y]*A[Y] )
#define DotXY( A, B ) ( A[X]*B[X] + A[Y]*B[Y] )
// Returns whether the cylinders are in contact
bool isContact( BCylinder const& a, BCylinder const& b,
                Transform const& a2w, Transform const& b2w )
{
  // Variables
  double rA = a.getRadius();
  double rB = b.getRadius();
  double hA = a.getHeight();
  double hB = b.getHeight();
  Vector3 e_A = a2w.getBasis() * a.getAxis();
  Vector3 e_B = b2w.getBasis() * b.getAxis();

  Vector3 zAxis( 0., 0., 1. );
  bool contCond;

  // Relative positions - B w.r.t. A
  Matrix rotMatA = getRotationMatrix( e_A, zAxis );
  Point3 x_B2A = rotMatA * ( *( b2w.getOrigin() ) - *( a2w.getOrigin() ) );
  Vector3 e_B2A =  ( rotMatA * e_B ).normalized();
  e_B2A.round( EPSILON );
  // Relative positions - A w.r.t. B
  Matrix rotMatB = getRotationMatrix( e_B, zAxis );
  Point3 x_A2B = rotMatB * ( *( a2w.getOrigin() ) - *( b2w.getOrigin() ) );
  Vector3 e_A2B = ( rotMatB * e_A ).normalized();
  e_A2B.round( EPSILON );

  // Contact scenarios
  // Face-Face / Band-Band
  {
    contCond = ( fabs( e_B2A[Z] - 1. ) < EPSILON ) &&
               ( fabs( x_B2A[Z] ) < .5 * ( hA + hB ) ) &&
               ( NormXY( x_B2A ) < pow( rA + rB, 2 ) );
    if ( contCond == true ) return ( true );
  }
  // Face-Edge
  {
    Vector3 r = ( e_B2A ^ ( e_B2A ^ zAxis ) ).normalized();
    r = sgn( x_B2A[Z] ) * ( rB * r - .5 * hB * sgn( e_B2A[Z] ) * e_B2A );
    Point3 ptE = x_B2A + r;
    contCond = ( fabs( ptE[Z] ) < .5 * hA ) &&
               ( NormXY( ptE ) < pow( rA, 2 ) ) &&
               ( fabs( x_B2A[Z] ) > .5 * hA );
    if ( contCond == true ) return ( true );
  }
  // Edge-Face
  {
    Vector3 r = ( e_A2B ^ ( e_A2B ^ zAxis ) ).normalized();
    r = sgn( x_A2B[Z] ) * ( rA * r - .5 * hA * sgn( e_A2B[Z] ) * e_A2B );
    Point3 ptE = x_A2B + r;
    contCond = ( fabs( ptE[Z] ) < .5 * hB ) &&
               ( NormXY( ptE ) < pow( rB, 2 ) ) &&
               ( fabs( x_A2B[Z] ) > .5 * hB );
    if ( contCond == true ) return ( true );
  }
  // Band-Band
  {
    Vector3 r = ( zAxis ^ e_B2A ).normalized();
    double d = fabs( x_B2A * r );
    double lBStar = ( e_B2A[Z] * x_B2A[Z] - e_B2A * x_B2A )
                    / ( 1 - pow( e_B2A[Z], 2 ) );
    double lAStar = zAxis * ( x_B2A + lBStar * e_B2A );
    contCond = ( d < rA + rB ) &&
               ( fabs( lAStar ) < .5 * hA ) &&
               ( fabs( lBStar ) < .5 * hB );
    if ( contCond == true ) return ( true );
  }
  // Edge-Edge A2B
  {
    Point3 c1 = x_B2A + .5 * hB * e_B2A;
    Point3 c2 = x_B2A - .5 * hB * e_B2A;
    double d1 = pow( sqrt( NormXY( c1 ) ) - rA, 2) +
                pow( fabs( c1[Z] ) - .5 * hA, 2);
    double d2 = pow( sqrt( NormXY( c2 ) ) - rA, 2) +
                pow( fabs( c2[Z] ) - .5 * hA, 2);
    Point3 ptCenter = d2 < d1 ? c2 : c1; // decide on edge of B
    Point3 ptA, ptB;
    edgePointsOnCyl( rA, hA, rB, e_B2A, ptCenter, ptA, ptB );
    contCond = ( ( ptA != OriginePoint ) &&
                 ( fabs( ptA[Z] ) < .5 * hA ) );
    if ( contCond == true ) return ( true );
  }
  // Edge-Edge B2A
  {
    Point3 c1 = x_A2B + .5 * hA * e_A2B;
    Point3 c2 = x_A2B - .5 * hA * e_A2B;
    double d1 = pow( sqrt( NormXY( c1 ) ) - rB, 2) +
                pow( fabs( c1[Z] ) - .5 * hB, 2);
    double d2 = pow( sqrt( NormXY( c2 ) ) - rB, 2) +
                pow( fabs( c2[Z] ) - .5 * hB, 2);
    Point3 ptCenter = d2 < d1 ? c2 : c1; // decide on edge of A
    Point3 ptA, ptB;
    edgePointsOnCyl( rB, hB, rA, e_A2B, ptCenter, ptA, ptB );
    contCond = ( ( ptA != OriginePoint ) &&
                 ( fabs( ptA[Z] ) < .5 * hB ) );
    if ( contCond == true ) return ( true );
  }
  return ( false );

  // // Contact scenarios
  // // Face-Face / Band-Band
  // {
  //   contCond = ( fabs( e_B2A[Z] - 1. ) < EPSILON ) &&
  //              ( fabs( x_B2A[Z] ) < .5 * ( hA + hB ) ) &&
  //              ( NormXY( x_B2A ) < pow( rA + rB, 2 ) );
  //   if ( contCond == true ) return ( true );
  // }
  // // Face-Band
  // {
  //   double S = fabs( x_B2A * e_B2A );
  //   double t = sqrt( NormXY( x_B2A ) - pow( S, 2 ) );
  //   double tStar = S < .5 * hB ?
  //                 rA : sqrt( pow( rA, 2 ) - pow( S - .5 * hB, 2 ) );
  //   contCond = ( fabs( e_B2A[Z] ) < EPSILON ) &&
  //              ( fabs( x_B2A[Z] ) < .5 * hA + rB ) &&
  //              ( fabs( x_B2A[Z] ) > .5 * hA ) &&
  //              ( S < .5 * hB + rA ) &&
  //              ( t < tStar );
  //   if ( contCond == true ) return ( true );
  // }
  // // Band-Face
  // {
  //   double S = fabs( x_A2B * e_A2B );
  //   double t = sqrt( NormXY( x_A2B ) - pow( S, 2 ) );
  //   double tStar = S < .5 * hA ?
  //                 rB : sqrt( pow( rB, 2 ) - pow( S - .5 * hA, 2 ) );
  //   contCond = ( fabs( e_A2B[Z] ) < EPSILON ) &&
  //              ( fabs( x_A2B[Z] ) < .5 * hB + rA ) &&
  //              ( fabs( x_A2B[Z] ) > .5 * hB ) &&
  //              ( S < .5 * hA + rB ) &&
  //              ( t < tStar );
  //   if ( contCond == true ) return ( true );
  // }
  // // Face-Edge
  // {
  //   Vector3 r = ( e_B2A ^ ( e_B2A ^ zAxis ) ).normalized();
  //   r = sgn( x_B2A[Z] ) * ( rB * r - .5 * hB * sgn( e_B2A[Z] ) * e_B2A );
  //   Point3 ptE = x_B2A + r;
  //   contCond = ( fabs( ptE[Z] ) < .5 * hA ) &&
  //              ( NormXY( ptE ) < pow( rA, 2 ) ) &&
  //              ( .5 * hA - fabs( ptE[Z] ) < rA - sqrt( NormXY( ptE ) ) ) &&
  //              ( fabs( x_B2A[Z] ) > .5 * hA );
  //   if ( contCond == true ) return ( true );
  // }
  // // Edge-Face
  // {
  //   Vector3 r = ( e_A2B ^ ( e_A2B ^ zAxis ) ).normalized();
  //   r = sgn( x_A2B[Z] ) * ( rA * r - .5 * hA * sgn( e_A2B[Z] ) * e_A2B );
  //   Point3 ptE = x_A2B + r;
  //   contCond = ( fabs( ptE[Z] ) < .5 * hB ) &&
  //              ( NormXY( ptE ) < pow( rB, 2 ) ) &&
  //              ( .5 * hB - fabs( ptE[Z] ) < rB - sqrt( NormXY( ptE ) ) ) &&
  //              ( fabs( x_A2B[Z] ) > .5 * hB );
  //   if ( contCond == true ) return ( true );
  // }
  // // Band-Band
  // {
  //   Vector3 r = ( zAxis ^ e_B2A ).normalized();
  //   double d = fabs( x_B2A * r );
  //   double lBStar = ( e_B2A[Z] * x_B2A[Z] - e_B2A * x_B2A )
  //                   / ( 1 - pow( e_B2A[Z], 2 ) );
  //   double lAStar = zAxis * ( x_B2A + lBStar * e_B2A );
  //   contCond = ( d < rA + rB ) &&
  //              ( fabs( lAStar ) < .5 * hA ) &&
  //              ( fabs( lBStar ) < .5 * hB );
  //   if ( contCond == true ) return ( true );
  // }
  // // Band-Edge
  // {
  //   Point3 c1 = x_B2A + .5 * hB * e_B2A;
  //   Point3 c2 = x_B2A - .5 * hB * e_B2A;
  //   Point3 ptCenter = NormXY( c1 ) < NormXY( c2 ) ? c1 : c2; // decide the edge
  //   Point3 ptA;
  //   edgePointClose2Z( rA, rB, e_B2A, ptCenter, ptA );
  //   contCond = ( fabs( ptA[Z] ) < .5 * hA ) &&
  //              ( NormXY( ptA ) < pow( rA, 2 ) ) &&
  //              ( rA - sqrt( NormXY( ptA ) ) < .5 * hA - fabs( ptA[Z] ) );
  //   if ( contCond == true ) return ( true );
  // }
  // // Edge-Band
  // {
  //   Point3 c1 = x_A2B + .5 * hA * e_A2B;
  //   Point3 c2 = x_A2B - .5 * hA * e_A2B;
  //   Point3 ptCenter = NormXY( c1 ) < NormXY( c2 ) ? c1 : c2; // decide the edge
  //   Point3 ptA;
  //   edgePointClose2Z( rB, rA, e_A2B, ptCenter, ptA );
  //   contCond = ( fabs( ptA[Z] ) < .5 * hB ) &&
  //              ( NormXY( ptA ) < pow( rB, 2 ) ) &&
  //              ( rB - sqrt( NormXY( ptA ) ) < .5 * hB - fabs( ptA[Z] ) );
  //   if ( contCond == true ) return ( true );
  // }
  // // Edge-Edge
  // {
  //   Point3 c1 = x_B2A + .5 * hB * e_B2A;
  //   Point3 c2 = x_B2A - .5 * hB * e_B2A;
  //   double d1 = pow( sqrt( NormXY( c1 ) ) - rA, 2) +
  //               pow( fabs( c1[Z] ) - .5 * hA, 2);
  //   double d2 = pow( sqrt( NormXY( c2 ) ) - rA, 2) +
  //               pow( fabs( c2[Z] ) - .5 * hA, 2);
  //   Point3 ptCenter = d2 < d1 ? c2 : c1; // decide on edge of B
  //   Point3 ptA, ptB;
  //   edgePointsOnCyl( rA, hA, rB, e_B2A, ptCenter, ptA, ptB );
  //   contCond = ( ( ptA != OriginePoint ) &&
  //                ( fabs( ptA[Z] ) < .5 * hA ) &&
  //                ( ptB != OriginePoint ) &&
  //                ( NormXY( ptB ) < pow( rA, 2 ) ) );
  //   if ( contCond == true ) return ( true );
  // }
  // return ( false );
}




// ----------------------------------------------------------------------------
// Returns the contact point of two cylinders in the world of the first cylinder
// if the contact is either Face-Face or Band-Band (Parallel)
void F2FB2BParContact( double rA, double hA, double rB, double hB,
                       Vector3 const& e, Point3 const& x, PointContact& ptCont )
{
  // Contact variables
  double overlap;
  Vector3 contVec; // Contact vector directed from b to a
  Point3 contPt;
  bool contCond; // TRUE if contact occurs in this scenario

  // Contact
  contCond = ( fabs( e[Z] - 1. ) < 1.e-5 ) &&
             ( fabs( x[Z] ) < .5 * ( hA + hB ) ) &&
             ( NormXY( x ) < pow( rA + rB, 2 ) );

  if ( contCond )
  {
    double axialOverlap = .5 * ( hA + hB ) - fabs( x[Z] );
    double radialOverlap = ( rA + rB ) - sqrt( NormXY( x ) );
    if ( axialOverlap < radialOverlap ) // Face-Face contact
    {
      // amount of overlap
      overlap = axialOverlap;
      // contact vector
      contVec = Vector3( 0., 0., -sgn( x[Z] ) );
    }
    else // Band-Band (parallel) contact
    {
      // amount of overlap
      overlap = radialOverlap;
      // contact vector
      contVec = - x;
      contVec[Z] = 0.;
      contVec = contVec.normalized();
    }
    // contact point
    Vector3 r = Vector3( x[X], x[Y], 0.).normalized();
    contPt = (rA - .5 * radialOverlap ) * r; // assigning X & Y components
    contPt[Z] = x[Z] - .5 * sgn( x[Z] ) * ( hB - axialOverlap );
    // output
    ptCont.setContact( contPt );
    ptCont.setOverlapVector( overlap * contVec );
    ptCont.setOverlapDistance( - overlap );
  }
}




// ----------------------------------------------------------------------------
// Returns the contact point of two cylinders in the world of the first cylinder
// if the contact is Face-Band
void F2BContact( double rA, double hA, double rB, double hB,
                 Vector3 const& e, Point3 const& x, PointContact& ptCont )
{
  // Contact variables
  double overlap;
  Vector3 contVec; // Contact vector directed from b to a
  Point3 contPt;
  bool contCond; // TRUE if contact occurs in this scenario

  // Some variables
  double S = fabs( x * e );
  double t = sqrt( NormXY( x ) - pow( S, 2 ) );
  double tStar = S < .5 * hB ? rA : sqrt( pow( rA, 2 ) - pow( S - .5 * hB, 2 ) );

  // Contact
  contCond = ( fabs( e[Z] ) < 1.e-5 ) &&
             ( fabs( x[Z] ) < .5 * hA + rB ) &&
             ( S < .5 * hB + rA ) &&
             ( t < tStar );

  if ( contCond )
  {
    // additional condition - assuring Face of A is in contact with Band of B
    if ( NormXY( x ) > rA * rA &&
         ( .5 * hB + rA - S ) < .5 * hA + rB - fabs( x[Z] ) )
      return;
    // amount of overlap
    overlap = .5 * hA + rB - fabs( x[Z] );
    // contact vector
    contVec = Vector3( 0., 0., - sgn( x[Z] ));
    // contact point
    Point3 ptA = x + rB * contVec + .5 * hB * e;
    Point3 ptB = x + rB * contVec - .5 * hB * e;
    double kappa[2];
    solveQuadratic( NormXY( ptA ) + NormXY( ptB ) - 2 * DotXY( ptA, ptB ),
                    - 2 * NormXY( ptB ) + 2 * DotXY( ptA, ptB ),
                    NormXY( ptB ) - pow( rA, 2 ),
                    kappa );
    if ( NormXY( ptA ) > pow( rA, 2 ) && NormXY( ptB ) < pow( rA, 2 ) )
    {
      double k = fabs( kappa[0] - .5 ) < .5 ? kappa[0] : kappa[1];
      ptA = k * ptA + ( 1. - k ) * ptB;
    }
    if ( NormXY( ptA ) < pow( rA, 2 ) && NormXY( ptB ) > pow( rA, 2 ) )
    {
      double k = fabs( kappa[0] - .5 ) < .5 ? kappa[0] : kappa[1];
      ptB = k * ptA + ( 1. - k ) * ptB;
    }
    if ( NormXY( ptA ) > pow( rA, 2 ) && NormXY( ptB ) > pow( rA, 2 ) )
    {
      Point3 temp = kappa[0] * ptA + ( 1. - kappa[0] ) * ptB;
      ptB = kappa[1] * ptA + ( 1. - kappa[1] ) * ptB;
      ptA = temp;
    }
    contPt = .5 * ( ptA + ptB ) + .5 * overlap * contVec;
    // output
    ptCont.setContact( contPt );
    ptCont.setOverlapVector( overlap * contVec );
    ptCont.setOverlapDistance( - overlap );
  }
}




// ----------------------------------------------------------------------------
// Returns the contact point of two cylinders in the world of the first cylinder
// if the contact is Face-Edge
void F2EContact( double rA, double hA, double rB, double hB,
                 Vector3 const& e, Point3 const& x, PointContact& ptCont )
{
  // Contact variables
  double overlap;
  Vector3 contVec; // Contact vector directed from b to a
  Point3 contPt;
  bool contCond; // TRUE if contact occurs in this scenario

  // Some variables
  Vector3 r = ( e ^ ( e ^ Vector3( 0., 0., 1. ) ) ).normalized();
  r = sgn( x[Z] ) * ( rB * r - .5 * hB * sgn( e[Z] ) * e );
  Point3 ptE = x + r;

  // Contact
  contCond = ( fabs( ptE[Z] ) < .5 * hA ) &&
             ( NormXY( ptE ) < pow( rA, 2 ) ) &&
             // ( .5 * hA - fabs( ptE[Z] ) < rA - sqrt( NormXY( ptE ) ) ) &&
             ( fabs( x[Z] ) > .5 * hA );

  if ( contCond )
  {
    // amount of overlap
    overlap = .5 * hA - fabs( ptE[Z] );
    // contact vector
    contVec = Vector3( 0., 0., - sgn( x[Z] ));
    // contact point
    contPt = ptE - .5 * overlap * contVec;
    // output
    ptCont.setContact( contPt );
    ptCont.setOverlapVector( overlap * contVec );
    ptCont.setOverlapDistance( - overlap );
  }
}




// ----------------------------------------------------------------------------
// Returns the contact point of two cylinders in the world of the first cylinder
// if the contact is Band-Band (Skewed)
void B2BSkewContact( double rA, double hA, double rB, double hB,
                     Vector3 const& e, Point3 const& x, PointContact& ptCont )
{
  // Contact variables
  double overlap;
  Vector3 contVec; // Contact vector directed from b to a
  Point3 contPt;
  bool contCond; // TRUE if contact occurs in this scenario

  // Some variables
  Vector3 r = ( Vector3( 0., 0., 1. ) ^ e ).normalized();
  double d = fabs( x * r );
  double lBStar = ( e[Z] * x[Z] - e * x ) / ( 1 - pow( e[Z], 2 ) );
  double lAStar = Vector3( 0., 0., 1. ) * ( x + lBStar * e );

  // Contact
  contCond = ( d < rA + rB ) &&
             ( fabs( lAStar ) < .5 * hA ) &&
             ( fabs( lBStar ) < .5 * hB );

  if ( contCond )
  {
    // amount of overlap
    overlap = rA + rB - d;
    // contact vector
    Point3 ptP = lAStar * Vector3( 0., 0., 1. );
    Point3 ptQ = x + lBStar * e;
    contVec = Vector3( ptP - ptQ ).normalized();
    // contact point
    contPt = .5 * ( ptP + ptQ );
    // output
    ptCont.setContact( contPt );
    ptCont.setOverlapVector( overlap * contVec );
    ptCont.setOverlapDistance( - overlap );
  }
}




// ----------------------------------------------------------------------------
// Returns the contact point of two cylinders in the world of the first cylinder
// if the contact is Band-Edge
void B2EContact( double rA, double hA, double rB, double hB,
                 Vector3 const& e, Point3 const& x, PointContact& ptCont )
{
  // Contact variables
  double overlap;
  Vector3 contVec; // Contact vector directed from b to a
  Point3 contPt;
  bool contCond; // TRUE if contact occurs in this scenario

  // Some variables
  Point3 c1 = x + .5 * hB * e;
  Point3 c2 = x - .5 * hB * e;
  Point3 ptCenter = NormXY( c1 ) < NormXY( c2 ) ? c1 : c2; // decide on the edge
  // Vector3 majorAxis = ( e ^ Vector3( 0., 0., 1. ) ).normalized();
  // Vector3 minorAxis = ( majorAxis ^ Vector3( 0., 0., 1. ) ).normalized();
  // Vector3 minorAxisInPlane = ( majorAxis ^ e ).normalized();
  // double majorRadius = rB;
  // double minorRadius = fabs( minorAxisInPlane * minorAxis ) * rB;
  // double theta = atan2( - majorRadius * ( majorAxis * ptCenter ),
  //                       - minorRadius * ( minorAxis * ptCenter ) );
  // Point3 ptA = ptCenter
  //             + rB * cos( theta ) * minorAxisInPlane
  //             + rB * sin( theta ) * majorAxis;
  Point3 ptA;
  edgePointClose2Z( rA, rB, e, ptCenter, ptA );

  // Contact
  contCond = ( ptA != OriginePoint ) &&
             ( fabs( ptA[Z] ) < .5 * hA ) &&
             ( NormXY( ptA ) < pow( rA, 2 ) ) &&
             ( rA - sqrt( NormXY( ptA ) ) < .5 * hA - fabs( ptA[Z] ) );

  if ( contCond )
  {
    // amount of overlap
    overlap = rA - sqrt( NormXY( ptA ) );
    // contact vector
    contVec = ptA;
    contVec[Z] = 0.;
    contVec = contVec.normalized();
    contVec = - contVec;
    // contact point
    contPt = ptA - .5 * overlap * contVec;
    // output
    ptCont.setContact( contPt );
    ptCont.setOverlapVector( overlap * contVec );
    ptCont.setOverlapDistance( - overlap );
  }
}




// ----------------------------------------------------------------------------
// Returns the contact point of two cylinders in the world of the first cylinder
// if the contact is Edge-Edge
void E2EContact( double rA, double hA, double rB, double hB,
                 Vector3 const& e, Point3 const& x,
                 PointContact& ptCont )
{
  // Contact variables
  double overlap;
  Vector3 contVec; // Contact vector directed from b to a
  Point3 contPt;
  bool contCond; // TRUE if contact occurs in this scenario

  // Some variables
  Point3 c_t = x + .5 * hB * e;
  Point3 c_b = x - .5 * hB * e;
  double d_t = pow( sqrt( NormXY( c_t ) ) - rA, 2) +
               pow( fabs( c_t[Z] ) - .5 * hA, 2);
  double d_b = pow( sqrt( NormXY( c_b ) ) - rA, 2) +
               pow( fabs( c_b[Z] ) - .5 * hA, 2);
  Point3 ptCenter = d_b < d_t ? c_b : c_t; // decide on edge of B
  Point3 ptA, ptB;
  edgePointsOnCyl( rA, hA, rB, e, ptCenter, ptA, ptB );

  // contact
  contCond = ( ptA != OriginePoint ) &&
             ( fabs( ptA[Z] ) < .5 * hA ) &&
             ( ptB != OriginePoint ) &&
             ( NormXY( ptB ) < pow( rA, 2 ) );

  if ( contCond )
  {
    // Finding contacting edge points of cylinder A
    Matrix rotMatA2B = getRotationMatrix( e, Vector3( 0., 0., 1. ) );
    Point3 x_A2B = rotMatA2B * ( - x );
    Vector3 e_A2B = ( rotMatA2B * Vector3( 0., 0., 1. ) ).normalized();
    c_t = x_A2B + .5 * hA * e_A2B;
    c_b = x_A2B - .5 * hA * e_A2B;
    d_t = pow( sqrt( NormXY( c_t ) ) - rB, 2) +
          pow( fabs( c_t[Z] ) - .5 * hB, 2);
    d_b = pow( sqrt( NormXY( c_b ) ) - rB, 2) +
          pow( fabs( c_b[Z] ) - .5 * hB, 2);
    ptCenter = d_b < d_t ? c_b : c_t;
    Point3 ptC, ptD;
    edgePointsOnCyl( rB, hB, rA, e_A2B, ptCenter, ptC, ptD );
    ptC = x + Point3( transpose( rotMatA2B ) * ( ptC ) );
    ptD = x + Point3( transpose( rotMatA2B ) * ( ptD ) );

    // Distance between two lines
    Vector3 vec_rA = .5 * ( ptC + ptD );
    Vector3 vec_eA = ( Vector3( ptD - ptC ) ).normalized();
    Vector3 vec_rB = .5 * ( ptA + ptB );
    Vector3 vec_eB = ( Vector3( ptB - ptA ) ).normalized();
    Vector3 vec_rAB = vec_rA - vec_rB;
    double muStar = - ( vec_rAB * vec_eA - ( vec_rAB * vec_eB ) * ( vec_eA * vec_eB ) ) / ( 1. - pow( vec_eA * vec_eB, 2 ) );
    double lambdaStar = ( vec_rAB * vec_eB - ( vec_rAB * vec_eA ) * ( vec_eA * vec_eB ) ) / ( 1. - pow( vec_eA * vec_eB, 2 ) );
    Point3 ptPa = vec_rA + muStar * vec_eA;
    Point3 ptPb = vec_rB + lambdaStar * vec_eB;
    // amount of overlap
    overlap = Norm( vec_rAB + muStar * vec_eA - lambdaStar * vec_eB );
    // contact vector
    contVec = ( Vector3( ptPb - ptPa ) ).normalized();
    // contact point
    contPt = .5 * ( ptPa + ptPb );
    // output
    ptCont.setContact( contPt );
    ptCont.setOverlapVector( overlap * contVec );
    ptCont.setOverlapDistance( - overlap );
  }
}




// ----------------------------------------------------------------------------
// Returns the band (ptA) and Face (ptB) points of the given circle
// intersecting with a cylinder oriented along Z axis and centered at the
// origin with radius rA and height hA
void edgePointClose2Z( double rA, double r, Vector3 const& e,
                       Point3 const& c, Point3& ptA )
{
  Vector3 u = ( e ^ Vector3( 0., 0., 1. ) ).normalized();
  Vector3 v = ( u ^ e ).normalized();

  // Finding the band point
  double A, B, C, D, E; // coefficients of degree-four polynomial
  double RVxVyUxUy, UxCxUyCy, CxVxCyVy; // misc variables
  RVxVyUxUy = r * ( NormXY( v ) - NormXY( u ) );
  UxCxUyCy = DotXY( u, c );
  CxVxCyVy = DotXY( c, v );

  A = pow(RVxVyUxUy, 2);
  B = 2. * CxVxCyVy * RVxVyUxUy;
  C = pow(UxCxUyCy, 2) + pow(CxVxCyVy, 2) - pow(RVxVyUxUy, 2);
  D = - 2. * CxVxCyVy * RVxVyUxUy;
  E = - pow(CxVxCyVy, 2);
  // std::cout << "\n\n" << A/A << " " << B/A << " " << C/A << " " << D/A << " " << E/A << '\n';
  double sol[4];
  solveQuartic( A, B, C, D, E, sol );
  // std::cout << sol[0] << " " << sol[1] << " " << sol[2] << " " << sol[3] << '\n';
  double sn, cs;
  for ( int i = 0; i < 4; i++ )
  {
    sn = sol[i];
    if ( fabs( sn ) <= 1 )
    {
      cs = ( UxCxUyCy * sn ) / ( CxVxCyVy + RVxVyUxUy * sn );
      if ( fabs( cs ) > 1. )
        cs = sgn( cs ) * sqrt( 1. - sn * sn );
      ptA = c + r * cs * u + r * sn * v;
      if ( NormXY( ptA ) < pow( rA, 2 ) )
        break;
    }
  }
}




// ----------------------------------------------------------------------------
// Returns the band (ptA) and Face (ptB) points of the given circle
// intersecting with a cylinder oriented along Z axis and centered at the
// origin with radius rA and height hA
void edgePointsOnCyl( double rA, double hA, double rB, Vector3 const& e,
                      Point3 const& c, Point3& ptA, Point3& ptB )
{
  Vector3 u = ( e ^ Vector3( 0., 0., 1. ) ).normalized();
  Vector3 v = ( u ^ e ).normalized();

  // Finding the band point
  double A, B, C, D, E; // coefficients of degree-four polynomial
  double VxVyUxUy, CxUxCyUy, CxVxCyVy, RCxCyUxUy; // misc variables
  VxVyUxUy = rB * rB * ( NormXY( v ) - NormXY( u ) );
  CxUxCyUy = rB * DotXY( c, u );
  CxVxCyVy = rB * DotXY( c, v );
  RCxCyUxUy = rA * rA - NormXY( c ) - rB * rB * NormXY( u );

  A = pow(VxVyUxUy, 2);
  B = 4. * VxVyUxUy * CxVxCyVy;
  C = 4. * pow(CxUxCyUy, 2) + 4. * pow(CxVxCyVy, 2) - 2. * VxVyUxUy * RCxCyUxUy;
  D = - 4. * CxVxCyVy * RCxCyUxUy;
  E = pow(RCxCyUxUy, 2) - 4. * pow(CxUxCyUy, 2);

  double sol[4];
  solveQuartic( A, B, C, D, E, sol );
  double sn, cs;
  for ( int i = 0; i < 4; i++ )
  {
    sn = sol[i];
    if ( fabs( sn ) <= 1 )
    {
      cs = ( RCxCyUxUy - 2. * CxVxCyVy * sn - VxVyUxUy * sn * sn )
                                                            / ( 2. * CxUxCyUy );
      if ( fabs( cs ) > 1. )
        cs = sgn( cs ) * sqrt( 1. - sn * sn );
      ptA = c + rB * cs * u + rB * sn * v;
      if ( fabs( ptA[Z] ) < .5 * hA )
        break;
    }
  }

  // Finding the face point
  sn = ( .5 * sgn( ptA[Z] ) * hA - c[Z] ) / ( rB * v[Z] );
  for ( int i = 0; i < 2; i++ )
  {
    cs = pow( -1., i) * sqrt( 1. - pow( sn, 2 ) );
    ptB = c + rB * cs * u + rB * sn * v;
    if ( NormXY( ptB ) < pow( rA, 2 ) )
      break;
  }
}



// ----------------------------------------------------------------------------
// Returns the solutions to a quartic equation
void solveQuartic( double a, double b, double c, double d, double e,
                   double sol[4] )
{
  // For more stability
  e /= a; d /= a; c /= a; b /= a; a /= a;
  // Deprressed quartic: y^4 + p*y^2 + q*y + r = 0
  double p = ( 8.*c - 3.*b*b ) / 8.;
  double q = ( b*b*b - 4.*b*c + 8.*d ) / 8.;
  double r = ( -3.*b*b*b*b + 256.*e - 64.*b*d + 16.*b*b*c ) / 256.;

  // Solve
  double m = solveCubicReal( 8, 8*p, 2*p*p - 8*r, -q*q );
  if ( fabs( m ) < EPSILON2 )
  {
    double temp[2];
    solveQuadratic( 1., p, r, temp );
    sol[0] = sqrt( temp[0] ) - b/4.;
    sol[1] = - sqrt( temp[0] ) - b/4.;
    sol[2] = sqrt( temp[1] ) - b/4.;
    sol[3] = - sqrt( temp[1] ) - b/4.;
  }
  else
  {
    sol[0] = ( sqrt(2*m) + sqrt( - ( 2*p + 2*m + sqrt(2/m)*q ) ) ) / 2. - b/4.;
    sol[1] = ( sqrt(2*m) - sqrt( - ( 2*p + 2*m + sqrt(2/m)*q ) ) ) / 2. - b/4.;
    sol[2] = ( -sqrt(2*m) + sqrt( - ( 2*p + 2*m - sqrt(2/m)*q ) ) ) / 2. - b/4.;
    sol[3] = ( -sqrt(2*m) - sqrt( - ( 2*p + 2*m - sqrt(2/m)*q ) ) ) / 2. - b/4.;
  }
}




// ----------------------------------------------------------------------------
// Returns a REAL solution to a cubic equation
double solveCubicReal( double a, double b, double c, double d )
{
  // For more stability
  d = d/a; c = c/a; b = b/a; a = a/a;
  // Deprressed cubic: y^3 + p*y + q = 0
  double p = ( 3.*c - b*b ) / 3.;
  double q = ( 2*b*b*b - 9.*b*c + 27.*d ) / 27.;

  double del = q*q/4. + p*p*p/27.;
  double u;
  if ( del > 0 )
      u = pow( -q/2. + sqrt( del ), 1./3. );
  else
      u = pow( sqrt( q*q/4. + abs(del) ), 1./3. );
  return ( u - p / 3. / u - b / 3. );
}




// ----------------------------------------------------------------------------
// Returns solutions to a quadratic equation
void solveQuadratic( double a, double b, double c, double sol[2] )
{
  double delta = b * b - 4 * a * c;
  sol[0] = ( - b + sqrt( delta ) ) / ( 2 * a );
  sol[1] = ( - b - sqrt( delta ) ) / ( 2 * a );
}




// ----------------------------------------------------------------------------
// Sign function
template <typename T>
int sgn(T val)
{
    return ( ( T(0) < val ) - ( val < T(0) ) );
}
#undef NormXY
#undef DotXY
