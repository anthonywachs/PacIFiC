#include "BCylinder.hh"
#include <iostream>

using namespace std;

double tol = 1.e-6; // Tolerance used in this class

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
// Sign function
template < typename T >
inline int sgn( T val )
{
    return ( ( T(0) < val ) - ( val < T(0) ) );
}




// ----------------------------------------------------------------------------
// Returns the norm of a Point3 object in the xy-plane
inline double normXY( Point3 x )
{
  return ( x[X]*x[X] + x[Y]*x[Y] );
}




// ----------------------------------------------------------------------------
// Returns the dot product of two Vector3 objects in the xy-plane
inline double dotXY( Vector3 x, Vector3 y )
{
  return ( x[X]*y[X] + x[Y]*y[Y] );
}




// ----------------------------------------------------------------------------
// Returns solutions to the quadratic equation ax^2 + bx + c
inline void solveQuadratic( double const a, double const b, double const c,
                            double sol[2] )
{
  double delta = b * b - 4 * a * c;
  sol[0] = ( - b + sqrt( delta ) ) / ( 2 * a );
  sol[1] = ( - b - sqrt( delta ) ) / ( 2 * a );
}




// ----------------------------------------------------------------------------
// Returns the solutions to to the quartic equation x^4 + bx^3 + cx^2 + dx + e
inline void solveQuartic( double const b, double const c, double const d,
                          double const e, double sol[4], int& nbRoots )
{
  // reseting the number of roots
  nbRoots = 0;
  // Deprressed quartic: y^4 + p*y^2 + q*y + r = 0
  double const b2 = b*b;
  double const p = c - 3.*b2/8.;
  double const q = b2*b/8. - b*c/2. + d;
  double const r = -3.*b2*b2/256. + e - b*d/4. + b2*c/16.;
  double const p2 = p*p;

  // Solve
  if ( fabs( q ) < EPSILON )
  {
    // finding solutions to the quadratic equation x^2 + px + r = 0.
    double const del = p2 / 4. - r; // this is actually del/4.!
    if ( del < 0. )
      return;
    else
    {
      double const m1 = - p / 2. + sqrt( del );
      double const m2 = - p / 2. - sqrt( del );
      if ( m1 > 0. )
      {
        sol[ nbRoots++ ] = sqrt( m1 ) - b / 4.;
        sol[ nbRoots++ ] = - sqrt( m1 ) - b / 4.;
      }
      if ( m2 > 0. )
      {
        sol[ nbRoots++ ] = sqrt( m2 ) - b / 4.;
        sol[ nbRoots++ ] = - sqrt( m2 ) - b / 4.;
      }
    }
  }
  else
  {
    // finding a real root to cubic equation x^3 + px^2 + (p*p/4. - r)x - q*q/8.
    double const u = -p2/36. - r/3.; // this is actually p/3.!
    double const v = -p2*p/216. + r*p/6. - q*q/16.; // this is actually v/2.!

    double const del = u*u*u + v*v;
    double m = 0.;
    if ( del < 0 )
      m = 2. * sqrt( -u ) * cos( acos( v / sqrt( -u ) / u ) / 3. ) - p / 3.;
    else
    {
      m = cbrt( -v + sqrt( del ) );
      m = m - u / m - p / 3.;
    }

    // roots
    if ( m < 0. )
      return;
    else
    {
      double const sqrt_mhalf = sqrt( m / 2. );
      double const first_var = - p / 2. - m / 2. - q / sqrt_mhalf / 4.;
      double const second_var = first_var + q / sqrt_mhalf / 2.;

      if ( first_var > 0. )
      {
        sol[ nbRoots++ ] = sqrt_mhalf + sqrt( first_var ) - b / 4.;
        sol[ nbRoots++ ] = sqrt_mhalf - sqrt( first_var ) - b / 4.;
      }
      if ( second_var > 0. )
      {
        sol[ nbRoots++ ] = - sqrt_mhalf + sqrt( second_var ) - b / 4.;
        sol[ nbRoots++ ] = - sqrt_mhalf - sqrt( second_var ) - b / 4.;
      }
    }
  }
}




// ----------------------------------------------------------------------------
// Rotation matrix that transforms v to Vector3(0., 0., 1.)
inline void rotateVec2VecZ( Vector3 const& v, Matrix& rotMat )
{
  double const c = 1. + v[Z];
  if ( fabs( c ) < tol )
    rotMat.setValue( -1.,  0.,  0.,
                      0., -1.,  0.,
                      0.,  0., -1. );
  else
  {
    double const vx2 = v[X]*v[X]/c;
    double const vy2 = v[Y]*v[Y]/c;
    double const vxvy = v[X]*v[Y]/c;
    rotMat.setValue( 1. - vx2, -vxvy, -v[X],
                     -vxvy, 1. - vy2, -v[Y],
                     v[X], v[Y], 1. - vx2 - vy2 );
  }

  rotMat.round( EPSILON2 );
}




// ----------------------------------------------------------------------------
// Returns the contact point of two cylinders
PointContact intersect( BCylinder const& a, BCylinder const& b,
                        Transform const& a2w, Transform const& b2w )
{
  // Variables
  double const rA = a.getRadius();
  double const rB = b.getRadius();
  double const hA = a.getHeight();
  double const hB = b.getHeight();
  Vector3 const& e_A = a2w.getBasis() * a.getAxis();
  Vector3 const& e_B = b2w.getBasis() * b.getAxis();

  PointContact ptCont = PointNoContact;

  // Relative positions - B w.r.t. A
  Matrix rotMatA;
  rotateVec2VecZ( e_A, rotMatA );
  Point3 const& x_B2A = rotMatA * ( *(b2w.getOrigin()) - *(a2w.getOrigin()) );
  Vector3 const& e_B2A = ( rotMatA * e_B ).normalized();
  // e_B2A.round( EPSILON );
  // Relative positions - A w.r.t. B
  Matrix rotMatB;
  rotateVec2VecZ( e_B, rotMatB );
  Point3 const& x_A2B = rotMatB * ( *(a2w.getOrigin()) - *(b2w.getOrigin()) );
  Vector3 const& e_A2B = ( rotMatB * e_A ).normalized();
  // e_A2B.round( EPSILON );

  int counter = 0;
  while( ptCont.getOverlapDistance() >= 0 && counter < 10 )
  {
    counter++;
    switch ( counter )
    {
      case 1:
        F2FB2BParContact( rA, hA, rB, hB, e_B2A, x_B2A, ptCont );
        break;
      case 2:
        F2BContact( rA, hA, rB, hB, e_B2A, x_B2A, ptCont );
        break;
      case 3:
        F2BContact( rB, hB, rA, hA, e_A2B, x_A2B, ptCont );
        break;
      case 4:
        F2EContact( rA, hA, rB, hB, e_B2A, x_B2A, ptCont );
        break;
      case 5:
        F2EContact( rB, hB, rA, hA, e_A2B, x_A2B, ptCont );
        break;
      case 6:
        B2BSkewContact( rA, hA, rB, hB, e_B2A, x_B2A, ptCont );
        break;
      case 7:
        B2EContact( rA, hA, rB, hB, e_B2A, x_B2A, ptCont );
        break;
      case 8:
        B2EContact( rB, hB, rA, hA, e_A2B, x_A2B, ptCont );
        break;
      case 9:
        E2EContact( rA, hA, rB, hB, e_B2A, x_B2A, ptCont );
        break;
      default:
        break;
    }
  }

  switch ( counter )
  {
    case 1: case 2: case 4: case 6: case 7: case 9:
      ptCont.setContact( *a2w.getOrigin() +
                        Point3( transpose(rotMatA) * ptCont.getContact() ) );
      ptCont.setOverlapVector( transpose(rotMatA) * ptCont.getOverlapVector() );
      break;
    case 3: case 5: case 8:
      ptCont.setContact( *b2w.getOrigin() +
                        Point3( transpose(rotMatB) * ptCont.getContact() ) );
      ptCont.setOverlapVector( - ptCont.getOverlapVector() );
      ptCont.setOverlapVector( transpose(rotMatB) * ptCont.getOverlapVector() );
      break;
    default:
      break;
  }

  return ( ptCont );
}




// ----------------------------------------------------------------------------
// Returns whether the cylinders are in contact
bool isContact( BCylinder const& a, BCylinder const& b,
                Transform const& a2w, Transform const& b2w )
{
  // Variables
  double const rA = a.getRadius();
  double const rB = b.getRadius();
  double const hA = a.getHeight();
  double const hB = b.getHeight();
  Vector3 const& e_A = a2w.getBasis() * a.getAxis();
  Vector3 const& e_B = b2w.getBasis() * b.getAxis();

  // Relative positions - B w.r.t. A
  Matrix rotMatA;
  rotateVec2VecZ( e_A, rotMatA );
  Point3 const& x_B2A = rotMatA * ( *(b2w.getOrigin()) - *(a2w.getOrigin()) );
  Vector3 const& e_B2A = ( rotMatA * e_B ).normalized();
  // e_B2A.round( EPSILON );
  // Relative positions - A w.r.t. B
  Matrix rotMatB;
  rotateVec2VecZ( e_B, rotMatB );
  Point3 const& x_A2B = rotMatB * ( *(a2w.getOrigin()) - *(b2w.getOrigin()) );
  Vector3 const& e_A2B = ( rotMatB * e_A ).normalized();
  // e_A2B.round( EPSILON );


  // General variables to use later
  // Vector3 const& u1 = ( e_B2A ^ Vector3( 0., 0., 1. ) ).normalized();
  Vector3 const& u1 = ( Vector3( e_B2A[Y], -e_B2A[X], 0. ) ).normalized();
  Vector3 const& v1 = ( u1 ^ e_B2A ).normalized();

  // Vector3 const& u2 = ( e_A2B ^ Vector3( 0., 0., 1. ) ).normalized();
  Vector3 const& u2 = ( Vector3( e_A2B[Y], -e_A2B[X], 0. ) ).normalized();
  Vector3 const& v2 = ( u2 ^ e_A2B ).normalized();


  // Contact scenarios
  // Face-Face / Band-Band (Parallel)
  if ( ( fabs( fabs( e_B2A[Z] ) - 1. ) < tol ) &&
       ( fabs( x_B2A[Z] ) < .5 * ( hA + hB ) ) &&
       ( normXY( x_B2A ) < ( rA + rB )*( rA + rB ) ) )
    return ( true );
  // Face-Edge
  if ( fabs( x_B2A[Z] ) > .5 * hA )
  {
    // Vector3 r = ( e_B2A ^ ( e_B2A ^ zAxis ) ).normalized();
    // r = sgn( x_B2A[Z] ) * ( rB * r - .5 * hB * sgn( e_B2A[Z] ) * e_B2A );
    Vector3 r = sgn( -x_B2A[Z] ) * ( rB * v1 + .5 * hB * sgn( e_B2A[Z] ) * e_B2A );
    Point3 ptE = x_B2A + r;
    if ( ( fabs( ptE[Z] ) < .5 * hA ) && ( normXY( ptE ) < rA * rA ) )
      return ( true );
  }
  // Edge-Face
  if ( fabs( x_A2B[Z] ) > .5 * hB )
  {
    // Vector3 r = ( e_A2B ^ ( e_A2B ^ zAxis ) ).normalized();
    // r = sgn( x_A2B[Z] ) * ( rA * r - .5 * hA * sgn( e_A2B[Z] ) * e_A2B );
    Vector3 r = sgn( -x_A2B[Z] ) * ( rA * v2 + .5 * hA * sgn( e_A2B[Z] ) * e_A2B );
    Point3 ptE = x_A2B + r;
    if ( ( fabs( ptE[Z] ) < .5 * hB ) && ( normXY( ptE ) < rB * rB ) )
      return ( true );
  }
  // Band-Band (Skewed)
  // Vector3 r = ( zAxis ^ e_B2A ).normalized();
  // double d = fabs( x_B2A * r );
  double d = fabs( x_B2A * ( -u1 ) );
  if ( d < rA + rB )
  {
    double lBStar = ( e_B2A[Z] * x_B2A[Z] - e_B2A * x_B2A )
                    / ( 1 - e_B2A[Z]*e_B2A[Z] );
    if ( fabs( lBStar ) < .5 * hB )
    {
      double lAStar = x_B2A[Z] + lBStar * e_B2A[Z];
      if ( fabs( lAStar ) < .5 * hA )
        return ( true );
    }
  }
  // Band-Edge
  {
    Point3 c1 = x_B2A + .5 * hB * e_B2A;
    Point3 c2 = x_B2A - .5 * hB * e_B2A;
    Point3 const& ptCenter = normXY( c1 ) < normXY( c2 ) ? c1 : c2;

    // additional condition to avoid solving the polynomial
    if ( ( fabs( ptCenter[Z] ) < 0.5 * hA + rB ) &&
         ( fabs( sqrt( normXY( ptCenter ) ) - rA ) < rB ) )
    {
      // Misc variables
      double const p = rB * ( normXY( v1 ) - normXY( u1 ) );
      double const q = dotXY( u1, ptCenter ) / p;
      double const r = dotXY( ptCenter, v1 ) / p;

      double sint[4];
      int nbRoots = 0;
      solveQuartic( 2.*r, q*q + r*r - 1., -2.*r, -r*r, sint, nbRoots );

      Point3 ptA;
      double cost;
      for ( int i = 0; i < nbRoots; i++ )
      {
        if ( fabs( sint[i] ) <= 1. )
        {
          cost = ( q * sint[i] ) / ( r + sint[i] );
          if ( fabs( cost ) > 1. )
            cost = sgn( cost ) * sqrt( 1. - sint[i]*sint[i] );
          ptA = ptCenter + rB * cost * u1 + rB * sint[i] * v1;
          if ( normXY( ptA ) < rA * rA && fabs( ptA[Z] ) < .5 * hA )
            return ( true );
        }
      }
    }
  }
  // Edge-Band
  {
    Point3 c1 = x_A2B + .5 * hA * e_A2B;
    Point3 c2 = x_A2B - .5 * hA * e_A2B;
    Point3 const& ptCenter = normXY( c1 ) < normXY( c2 ) ? c1 : c2;

    // additional condition to avoid solving the polynomial
    if ( ( fabs( ptCenter[Z] ) < 0.5 * hB + rA ) &&
         ( fabs( sqrt( normXY( ptCenter ) ) - rB ) < rA ) )
    {
      // Misc variables
      double const p = rA * ( normXY( v2 ) - normXY( u2 ) );
      double const q = dotXY( u2, ptCenter ) / p;
      double const r = dotXY( ptCenter, v2 ) / p;

      double sint[4];
      int nbRoots = 0;
      solveQuartic( 2.*r, q*q + r*r - 1., -2.*r, -r*r, sint, nbRoots );

      Point3 ptA;
      double cost;
      for ( int i = 0; i < nbRoots; i++ )
      {
        if ( fabs( sint[i] ) <= 1. )
        {
          cost = ( q * sint[i] ) / ( r + sint[i] );
          if ( fabs( cost ) > 1. )
            cost = sgn( cost ) * sqrt( 1. - sint[i]*sint[i] );
          ptA = ptCenter + rA * cost * u2 + rA * sint[i] * v2;
          if ( normXY( ptA ) < rB * rB && fabs( ptA[Z] ) < .5 * hB )
            return ( true );
        }
      }
    }
  }
  // Edge-Edge
  {
    Point3 c1 = x_B2A + .5 * hB * e_B2A;
    Point3 c2 = x_B2A - .5 * hB * e_B2A;
    double d1 = pow( sqrt( normXY( c1 ) ) - rA, 2 ) +
                pow( fabs( c1[Z] ) - .5 * hA, 2 );
    double d2 = pow( sqrt( normXY( c2 ) ) - rA, 2 ) +
                pow( fabs( c2[Z] ) - .5 * hA, 2 );
    Point3 const& ptCenter = d1 < d2 ? c1 : c2; // decide on edge of B

    if ( ( d1 < d2 ? d1 : d2 ) <= rB * rB )
    {
      // Misc variables
      double const p = rB * rB * ( normXY( v1 ) - normXY( u1 ) ) / 2.;
      double const q = rB * dotXY( ptCenter, u1 ) / p;
      double const r = rB * dotXY( ptCenter, v1 ) / p;
      double const s = ( rA*rA - normXY( ptCenter ) - rB*rB*normXY( u1 ) ) / p;

      double sint[4];
      int nbRoots;
      solveQuartic( 2.*r, q*q + r*r - s, -r*s, s*s/4. - q*q, sint, nbRoots );

      Point3 ptA;
      double cost;
      for ( int i = 0; i < nbRoots; i++ )
      {
        if ( fabs( sint[i] ) <= 1 )
        {
          cost = ( s/2. - r*sint[i] - sint[i]*sint[i] ) / q;
          ptA = ptCenter + rB * cost * u1 + rB * sint[i] * v1;
          if ( fabs( ptA[Z] ) < .5 * hA )
          return ( true );
        }
      }
    }
  }

  return ( false );
}




// ----------------------------------------------------------------------------
// LOW-LEVEL ROUTINES
// ----------------------------------------------------------------------------
// Returns the contact point of two cylinders in the world of the first cylinder
// if the contact is either Face-Face or Band-Band (Parallel)
void F2FB2BParContact( double rA, double hA, double rB, double hB,
                       Vector3 const& e, Point3 const& x, PointContact& ptCont )
{
  // Contact happens if
  bool contCond = false; // TRUE if contact occurs in this scenario
  contCond = ( fabs( fabs( e[Z] ) - 1. ) < tol ) &&
             ( fabs( x[Z] ) < .5 * ( hA + hB ) ) &&
             ( normXY( x ) < ( rA + rB ) * ( rA + rB ) );

  if ( contCond )
  {
    // Contact variables
    double overlap;
    Vector3 contVec; // Contact vector directed from b to a
    Point3 contPt;

    double axialOverlap = .5 * ( hA + hB ) - fabs( x[Z] );
    double radialOverlap = ( rA + rB ) - sqrt( normXY( x ) );
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
    Vector3 r = Vector3( x[X], x[Y], 0. ).normalized();
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
  // Contact happens if
  bool contCond = false; // TRUE if contact occurs in this scenario
  if ( fabs( e[Z] ) < tol && fabs( x[Z] ) < .5 * hA + rB )
  {
    double S = fabs( x * e );
    if ( S < .5 * hB + rA )
    {
      double t = sqrt( normXY( x ) - S * S );
      double tStar = S < .5*hB ? rA : sqrt( rA*rA - (S - .5*hB) * (S - .5*hB) );
      // Last condition - assuring Face A is in contact with Band B
      contCond = ( t < tStar ) &&
          ( !( S > rA && ( .5 * hB + rA - S ) < .5 * hA + rB - fabs( x[Z] ) ) );
    }
  }

  if ( contCond )
  {
    // Contact variables
    double overlap;
    Vector3 contVec; // Contact vector directed from b to a
    Point3 contPt;

    // amount of overlap
    overlap = .5 * hA + rB - fabs( x[Z] );
    // contact vector
    contVec = Vector3( 0., 0., - sgn( x[Z] ));
    // contact point
    Point3 ptA = x + rB * contVec + .5 * hB * e;
    Point3 ptB = x + rB * contVec - .5 * hB * e;
    double kappa[2];
    solveQuadratic( normXY( ptA ) + normXY( ptB ) - 2 * dotXY( ptA, ptB ),
                    - 2 * normXY( ptB ) + 2 * dotXY( ptA, ptB ),
                    normXY( ptB ) - pow( rA, 2 ),
                    kappa );
    if ( normXY( ptA ) > rA * rA && normXY( ptB ) < rA * rA )
    {
      double k = fabs( kappa[0] - .5 ) < .5 ? kappa[0] : kappa[1];
      ptA = k * ptA + ( 1. - k ) * ptB;
    }
    if ( normXY( ptA ) < rA * rA && normXY( ptB ) > rA * rA )
    {
      double k = fabs( kappa[0] - .5 ) < .5 ? kappa[0] : kappa[1];
      ptB = k * ptA + ( 1. - k ) * ptB;
    }
    if ( normXY( ptA ) > rA * rA && normXY( ptB ) > rA * rA )
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
  // Contact happens if
  bool contCond = false; // TRUE if contact occurs in this scenario
  Point3 ptE;
  if ( fabs( x[Z] ) > .5 * hA )
  {
    // Vector3 r = ( e ^ ( e ^ zAxis ) ).normalized();
    Vector3 r = ( e ^ Vector3( e[Y], -e[X], 0. ) ).normalized();
    r = sgn( x[Z] ) * ( rB * r - .5 * hB * sgn( e[Z] ) * e );
    ptE = x + r;
    contCond = ( fabs( ptE[Z] ) < .5 * hA ) &&
               ( normXY( ptE ) < rA * rA ) &&
               ( .5 * hA - fabs( ptE[Z] ) < rA - sqrt( normXY( ptE ) ) );
  }

  if ( contCond )
  {
    // Contact variables
    double overlap;
    Vector3 contVec; // Contact vector directed from b to a
    Point3 contPt;

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
  // Contact happens if
  bool contCond = false; // TRUE if contact occurs in this scenario
  double lAStar, lBStar;
  // Vector3 r = ( Vector3( 0., 0., 1. ) ^ e ).normalized();
  Vector3 r = ( Vector3( -e[Y], e[X], 0. ) ).normalized();
  double d = fabs( x * r );
  if ( d < rA + rB )
  {
    lBStar = ( e[Z] * x[Z] - e * x ) / ( 1 - e[Z] * e[Z] );
    if ( fabs( lBStar ) < .5 * hB )
    {
      lAStar = x[Z] + lBStar * e[Z];
      contCond = ( fabs( lAStar ) < .5 * hA );
    }
  }

  if ( contCond )
  {
    // Contact variables
    double overlap;
    Vector3 contVec; // Contact vector directed from b to a
    Point3 contPt;

    // amount of overlap
    overlap = rA + rB - d;
    // contact vector
    Point3 ptP = Vector3( 0., 0., lAStar );
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
  // Contact happens if
  bool contCond = false; // TRUE if contact occurs in this scenario
  // Some variables
  Point3 const& c1 = x + .5 * hB * e;
  Point3 const& c2 = x - .5 * hB * e;
  Point3 const& ptCenter = normXY( c1 ) < normXY( c2 ) ? c1 : c2;
  // additional condition to avoid solving the polynomial
  if ( ( fabs( ptCenter[Z] ) > 0.5 * hA + rB ) ||
       ( fabs( sqrt( normXY( ptCenter ) ) - rA ) > rB ) )
    return;

  // Vector3 const& u = ( e ^ Vector3( 0., 0., 1. ) ).normalized();
  Vector3 const& u = ( Vector3( e[Y], -e[X], 0. ) ).normalized();
  Vector3 const& v = ( u ^ e ).normalized();

  // Misc variables
  double const a = rB * ( normXY( v ) - normXY( u ) );
  double const b = dotXY( u, ptCenter ) / a;
  double const c = dotXY( ptCenter, v ) / a;

  double sint[4];
  int nbRoots = 0;
  solveQuartic( 2.*c, b*b + c*c - 1., -2.*c, -c*c, sint, nbRoots );

  Point3 ptA;
  double cost;
  for ( int i = 0; i < nbRoots; i++ )
  {
    if ( fabs( sint[i] ) <= 1. )
    {
      cost = ( b * sint[i] ) / ( c + sint[i] );
      if ( fabs( cost ) > 1. )
        cost = sgn( cost ) * sqrt( 1. - sint[i]*sint[i] );
      ptA = ptCenter + rB * cost * u + rB * sint[i] * v;
      if ( normXY( ptA ) < rA * rA && fabs( ptA[Z] ) < .5 * hA )
      {
        contCond = ( rA - sqrt( normXY( ptA ) ) < .5 * hA - fabs( ptA[Z] ) );
        break;
      }
    }
  }


  // Contact
  if ( contCond )
  {
    // Contact variables
    double overlap;
    Vector3 contVec; // Contact vector directed from b to a
    Point3 contPt;

    // amount of overlap
    overlap = rA - sqrt( normXY( ptA ) );
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
// Returns the contact point of two cylinders in the global world if the contact
// is Edge-Edge
void E2EContact( double rA, double hA, double rB, double hB,
                 Vector3 const& e, Point3 const& x, PointContact& ptCont )
{
  // Contact happens if
  bool contCond = false; // TRUE if contact occurs in this scenario
  // Some variables
  Point3 c_t = x + .5 * hB * e;
  Point3 c_b = x - .5 * hB * e;
  double d_t = pow( sqrt( normXY( c_t ) ) - rA, 2 ) +
               pow( fabs( c_t[Z] ) - .5 * hA, 2 );
  double d_b = pow( sqrt( normXY( c_b ) ) - rA, 2 ) +
               pow( fabs( c_b[Z] ) - .5 * hA, 2 );
  Point3 const& ptCenter1 = d_b < d_t ? c_b : c_t; // decide on edge of B
  // additional condition to avoid solving the polynomial
  if ( d_t > rB * rB && d_b > rB * rB )
    return;

  // Vector3 const& u1 = ( e1 ^ Vector3( 0., 0., 1. ) ).normalized();
  Vector3 const& u1 = ( Vector3( e[Y], -e[X], 0. ) ).normalized();
  Vector3 const& v1 = ( u1 ^ e ).normalized();

  // Misc variables
  double const a1 = rB * rB * ( normXY( v1 ) - normXY( u1 ) ) / 2.;
  double const b1 = rB * dotXY( ptCenter1, u1 ) / a1;
  double const c1 = rB * dotXY( ptCenter1, v1 ) / a1;
  double const d1 = ( rA*rA - normXY( ptCenter1 ) - rB*rB*normXY( u1 ) ) / a1;

  double sint1[4];
  int nbRoots1;
  solveQuartic( 2.*c1, b1*b1 + c1*c1 - d1, -c1*d1, d1*d1/4. - b1*b1,
                sint1, nbRoots1 );

  Point3 ptA, ptB;
  double cost1;
  for ( int i = 0; i < nbRoots1; i++ )
  {
    if ( fabs( sint1[i] ) <= 1. )
    {
      cost1 = ( d1/2. - c1*sint1[i] - sint1[i]*sint1[i] ) / b1;
      // if ( fabs( cs ) > 1. )
      //   cs = sgn( cs ) * sqrt( 1. - sol[i]*sol[i] );
      ptA = ptCenter1 + rB * cost1 * u1 + rB * sint1[i] * v1;
      if ( fabs( ptA[Z] ) < .5 * hA )
      {
        // Finding the face point
        double sina1 = ( .5 * sgn( ptA[Z] ) * hA - ptCenter1[Z] ) / ( rB*v1[Z] );
        double cosa1 = sqrt( 1. - sina1 * sina1 );
        ptB = ptCenter1 + rB * cosa1 * u1 + rB * sina1 * v1;
        if ( normXY( ptB ) < rA * rA )
        {
          contCond = true;
          break;
        }
        ptB = ptCenter1 + rB * (-cosa1) * u1 + rB * sina1 * v1;
        if ( normXY( ptB ) < rA * rA )
        {
          contCond = true;
          break;
        }
      }
    }
  }

  // Contact
  if ( contCond )
  {
    // Contact variables
    double overlap;
    Vector3 contVec; // Contact vector directed from b to a
    Point3 contPt;

    // Finding contacting edge points of cylinder A
    Matrix rotMatA2B;
    rotateVec2VecZ( e, rotMatA2B );
    Point3 const& x2 = rotMatA2B * ( - x );
    Vector3 const& e2 = ( rotMatA2B * Vector3( 0., 0., 1. ) ).normalized();
    c_t = x2 + .5 * hA * e2;
    c_b = x2 - .5 * hA * e2;
    d_t = pow( sqrt( normXY( c_t ) ) - rB, 2) +
          pow( fabs( c_t[Z] ) - .5 * hB, 2);
    d_b = pow( sqrt( normXY( c_b ) ) - rB, 2) +
          pow( fabs( c_b[Z] ) - .5 * hB, 2);
    Point3 const& ptCenter2 = d_b < d_t ? c_b : c_t;

    // Vector3 const& u2 = ( e2 ^ Vector3( 0., 0., 1. ) ).normalized();
    Vector3 const& u2 = ( Vector3( e2[Y], -e2[X], 0. ) ).normalized();
    Vector3 const& v2 = ( u2 ^ e2 ).normalized();

    // Misc variables
    double const a2 = rA * rA * ( normXY( v2 ) - normXY( u2 ) ) / 2.;
    double const b2 = rA * dotXY( ptCenter2, u2 ) / a2;
    double const c2 = rA * dotXY( ptCenter2, v2 ) / a2;
    double const d2 = ( rB*rB - normXY( ptCenter2 ) - rA*rA*normXY( u2 ) ) / a2;

    double sint2[4];
    int nbRoots2;
    solveQuartic( 2.*c2, b2*b2 + c2*c2 - d2, -c2*d2, d2*d2/4. - b2*b2,
                  sint2, nbRoots2 );

    Point3 ptC, ptD;
    double cost2;
    for ( int i = 0; i < nbRoots2; i++ )
    {
      if ( fabs( sint2[i] ) <= 1. )
      {
        cost2 = ( d2/2. - c2*sint2[i] - sint2[i]*sint2[i] ) / b2;
        if ( fabs( cost2 ) > 1. )
          cost2 = sgn( cost2 ) * sqrt( 1. - sint2[i]*sint2[i] );
        ptC = ptCenter2 + rA * cost2 * u2 + rA * sint2[i] * v2;
        if ( fabs( ptC[Z] ) < .5 * hB )
          break;
      }
    }
    // Finding the face point
    double sina2 = ( .5 * sgn( ptC[Z] ) * hB - ptCenter2[Z] ) / ( rA*v2[Z] );
    double cosa2 = sqrt( 1. - sina2 * sina2 );
    ptD = ptCenter2 + rA * cosa2 * u2 + rA * sina2 * v2;
    if ( normXY( ptD ) > rB * rB )
      ptD = ptCenter2 + rA * (-cosa2) * u2 + rA * sina2 * v2;

    // Contact points in the coordinate system of cylinder A
    ptC = x + Point3( transpose( rotMatA2B ) * ( ptC ) );
    ptD = x + Point3( transpose( rotMatA2B ) * ( ptD ) );

    // Distance between two lines
    Vector3 vec_rA = .5 * ( ptC + ptD );
    Vector3 vec_eA = ( Vector3( ptD - ptC ) ).normalized();
    Vector3 vec_rB = .5 * ( ptA + ptB );
    Vector3 vec_eB = ( Vector3( ptB - ptA ) ).normalized();
    Vector3 vec_rAB = vec_rA - vec_rB;
    double a = vec_rAB * vec_eA;
    double b = vec_rAB * vec_eB;
    double c = vec_eA * vec_eB;
    double muStar = ( b * c - a ) / ( 1. - c * c );
    double lambdaStar = ( b - a * c ) / ( 1. - c * c );
    Point3 ptPa = vec_rA + muStar * vec_eA;
    Point3 ptPb = vec_rB + lambdaStar * vec_eB;
    // contact vector
    contVec = Vector3( ptPb - ptPa );
    // amount of overlap
    overlap = Norm( contVec );
    // contact point
    contPt = .5 * ( ptPa + ptPb );
    // output
    ptCont.setContact( contPt );
    ptCont.setOverlapVector( contVec );
    ptCont.setOverlapDistance( - overlap );
  }
}
// #undef normXY
// #undef dotXY
