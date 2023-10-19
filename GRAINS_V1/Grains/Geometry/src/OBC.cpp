#include "OBC.hh"

static double tol = EPSILON; // Tolerance used in this class
// --------------------------------------------------------------------
// Default constructor
OBC::OBC()
{}




// --------------------------------------------------------------------
// Constructor with radius r, height h, and initial orientation ori
OBC::OBC( double r, double h, Vector3 const& ori )
{
  m_radius = r;
  m_height = h;
  m_initOrientation = ori;
}




// --------------------------------------------------------------------
// Copy constructor
OBC::OBC( OBC const& obc_ )
{
  m_radius = obc_.m_radius;
  m_height = obc_.m_height;
  m_initOrientation = obc_.m_initOrientation;
}




// --------------------------------------------------------------------
// Destructor
OBC::~OBC()
{}




// --------------------------------------------------------------------
// Returns the OBC radius
BVolumeType OBC::getBVolumeType() const
{
  return ( typeOBC );
}




// --------------------------------------------------------------------
// Returns a clone of the OBB
BVolume* OBC::clone() const
{
  return ( new OBC( m_radius, m_height, m_initOrientation ) );
}




// --------------------------------------------------------------------
// Returns the OBC radius
double OBC::getRadius() const
{
  return ( m_radius );
}




// --------------------------------------------------------------------
// Returns the OBC height
double OBC::getHeight() const
{
  return ( m_height );
}




// --------------------------------------------------------------------
// Returns the OBC initial orientation
Vector3 const& OBC::getInitOrientation() const
{
  return ( m_initOrientation );
}




// --------------------------------------------------------------------
// Sets the OBC radius
void OBC::setRadius( double r )
{
  m_radius = r;
}




// --------------------------------------------------------------------
// Sets the OBC height
void OBC::setHeight( double h )
{
  m_height = h;
}




// --------------------------------------------------------------------
// Sets the OBC initial orientation
void OBC::setInitOrientation( Vector3 const& ori )
{
  m_initOrientation = ori;
}




// ----------------------------------------------------------------------------
// Output operator
void OBC::writeShape ( ostream& fileOut ) const
{
  fileOut << "*OBC " << m_radius << " " << m_height << " " << m_initOrientation;
}




// ----------------------------------------------------------------------------
// Sign function
template < typename T >
inline int sgn( T const val )
{
    return ( ( T(0) < val ) - ( val < T(0) ) );
}




// ----------------------------------------------------------------------------
// Returns the norm of a Point3 object in the xy-plane
inline double normXY( Point3 const& x )
{
  return ( x[X]*x[X] + x[Y]*x[Y] );
}




// ----------------------------------------------------------------------------
// Returns the dot product of two Vector3 objects in the xy-plane
inline double dotXY( Vector3 const& x, Vector3 const& y )
{
  return ( x[X]*y[X] + x[Y]*y[Y] );
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
  if ( fabs( c ) < EPSILON2 )
    rotMat.setValue( -1.,  0.,  0.,
                      0., -1.,  0.,
                      0.,  0., -1. );
  else if ( fabs( c - 2. ) < EPSILON2 )
    rotMat.setValue( 1.,  0.,  0.,
                     0., 1.,  0.,
                     0.,  0., 1. );
  else
  {
    double const vx2 = v[X]*v[X]/c;
    double const vy2 = v[Y]*v[Y]/c;
    double const vxvy = v[X]*v[Y]/c;
    rotMat.setValue( 1. - vx2, -vxvy, -v[X],
                     -vxvy, 1. - vy2, -v[Y],
                     v[X], v[Y], 1. - vx2 - vy2 );
  }
}




// ----------------------------------------------------------------------------
// Returns whether two OBCs are in contact
bool isContactBVolume( OBC const& obcA,
                       OBC const& obcB,
                       Transform const& a2w,
                       Transform const& b2w )
{
  // Variables
  double const rA = obcA.getRadius();
  double const rB = obcB.getRadius();
  double const hA = obcA.getHeight() / 2.; // half-height
  double const hB = obcB.getHeight() / 2.; // half-height
  Vector3 const& e_A = a2w.getBasis() * obcA.getInitOrientation();
  Vector3 const& e_B = b2w.getBasis() * obcB.getInitOrientation();

  // Relative positions - B w.r.t. A
  Matrix rotMatA;
  rotateVec2VecZ( e_A, rotMatA );
  Point3 const& x_B2A = rotMatA * ( *(b2w.getOrigin()) - *(a2w.getOrigin()) );
  Vector3 const& e_B2A = ( rotMatA * e_B ).normalized();

  // General variables to use later
  Vector3 const& u1 = ( Vector3( e_B2A[Y], -e_B2A[X], 0. ) ).normalized();
  Vector3 const& v1 = ( u1 ^ e_B2A ).normalized();

  //* Contact Scenarios: A primary, B secondary *//

  // Face A - Edge B
  if ( fabs( x_B2A[Z] ) > hA )
  {
    // Vector3 r = ( e_B2A ^ ( e_B2A ^ zAxis ) ).normalized();
    // r = sgn( x_B2A[Z] ) * ( rB * r - hB * sgn( e_B2A[Z] ) * e_B2A );
    Vector3 r = sgn( -x_B2A[Z] ) * ( rB * v1 + hB * sgn( e_B2A[Z] ) * e_B2A );
    Point3 const& ptE = x_B2A + r;
    if ( ( fabs( ptE[Z] ) < hA ) && ( normXY( ptE ) < rA * rA ) )
      return ( true );
  }

  // Edge A - Edge B
  {
    // To find the correct edge on the secondary cylinder
    Point3 c1 = x_B2A + hB * e_B2A;
    Point3 c2 = x_B2A - hB * e_B2A;
    double distBand1 = pow( sqrt( normXY( c1 ) ) - rA, 2 );
    double distBand2 = pow( sqrt( normXY( c2 ) ) - rA, 2 );
    double distEdge1 = distBand1 + pow( fabs( c1[Z] ) - hA, 2 );
    double distEdge2 = distBand2 + pow( fabs( c2[Z] ) - hA, 2 );
    if ( ( distBand1 - distBand2 ) * ( distEdge1 - distEdge2 ) < 0. &&
         ( distEdge1 < distEdge2 ? distEdge1 : distEdge2 ) >= rB * rB )
    {
      distEdge1 = distBand1;
      distEdge2 = distBand2;
    }
    Point3 const& ptCenter = distEdge1 < distEdge2 ? c1 : c2;

    if ( ( distBand1 < distBand2 ? distBand1 : distBand2 ) < rB * rB  &&
          fabs( ptCenter[Z] ) < hA + rB )
    {
      // Misc variables
      // double const p = rB * rB * ( normXY( v1 ) - normXY( u1 ) ) / 2.;
      double const p = rB * rB * ( normXY( v1 ) - 1. ) / 2.;
      double const q = rB * dotXY( u1, ptCenter ) / p;
      double const r = rB * dotXY( v1, ptCenter ) / p;
      // double const s = ( rA*rA - normXY( ptCenter ) - rB*rB*normXY(u1) ) / p;
      double const s = ( rA * rA - normXY( ptCenter ) - rB * rB ) / p;

      double sint[4];
      int nbRoots;
      solveQuartic( 2.*r, q*q + r*r - s, -r*s, s*s/4. - q*q, sint, nbRoots );

      for ( int i = 0; i < nbRoots; i++ )
      {
        if ( fabs( sint[i] ) <= 1. )
        {
          double cost = ( s/2. - r*sint[i] - sint[i]*sint[i] ) / q;
          double zVal = ptCenter[Z] + rB * cost * u1[Z] + rB * sint[i] * v1[Z];
          if ( fabs( zVal ) < hA )
            return ( true );
        }
      }
    }
  }

  // Special Case: Face A - Face B
  if ( ( fabs( fabs( e_B2A[Z] ) - 1. ) < tol ) &&
       ( fabs( x_B2A[Z] ) < hA + hB ) &&
       ( normXY( x_B2A ) < ( rA + rB ) * ( rA + rB ) ) )
    return ( true );

  // Special Case: Band A - Band B
  if ( fabs( x_B2A * ( - u1 ) ) < rA + rB )
  {
    double lBStar = ( e_B2A[Z] * x_B2A[Z] - e_B2A * x_B2A )
    / ( 1. - e_B2A[Z] * e_B2A[Z] );
    if ( fabs( lBStar ) < hB )
    {
      double lAStar = x_B2A[Z] + lBStar * e_B2A[Z];
      if ( fabs( lAStar ) < hA )
        return ( true );
    }
  }


  // Relative positions - A w.r.t. B
  Matrix rotMatB;
  rotateVec2VecZ( e_B, rotMatB );
  Point3 const& x_A2B = rotMatB * ( *(a2w.getOrigin()) - *(b2w.getOrigin()) );
  Vector3 const& e_A2B = ( rotMatB * e_A ).normalized();

  // Vector3 const& u2 = ( e_A2B ^ Vector3( 0., 0., 1. ) ).normalized();
  Vector3 const& u2 = ( Vector3( e_A2B[Y], -e_A2B[X], 0. ) ).normalized();
  Vector3 const& v2 = ( u2 ^ e_A2B ).normalized();

  //* Contact Scenarios: B secondary, A pirmary *//

  // Edge B - Face A
  if ( fabs( x_A2B[Z] ) > hB )
  {
    // Vector3 r = ( e_A2B ^ ( e_A2B ^ zAxis ) ).normalized();
    // r = sgn( x_A2B[Z] ) * ( rA * r - hA * sgn( e_A2B[Z] ) * e_A2B );
    Vector3 r = sgn( -x_A2B[Z] ) * ( rA * v2 + hA * sgn( e_A2B[Z] ) * e_A2B );
    Point3 ptE = x_A2B + r;
    if ( ( fabs( ptE[Z] ) < hB ) && ( normXY( ptE ) < rB * rB ) )
      return ( true );
  }

  // Edge B - Band A
  {
    Point3 c1 = x_A2B + hA * e_A2B;
    Point3 c2 = x_A2B - hA * e_A2B;
    Point3 const& ptCenter = normXY( c1 ) < normXY( c2 ) ? c1 : c2;

    // additional condition to avoid solving the polynomial
    if ( ( fabs( ptCenter[Z] ) < hB + rA ) &&
         ( sqrt( normXY( ptCenter ) ) < rA + rB ) )
    {
      // Misc variables
      // double const p = rA * ( normXY( v2 ) - normXY( u2 ) );
      double const p = rA * ( normXY( v2 ) - 1. );
      double const q = dotXY( u2, ptCenter ) / p;
      double const r = dotXY( v2, ptCenter ) / p;

      double sint[4];
      int nbRoots = 0;
      solveQuartic( 2.*r, q*q + r*r - 1., -2.*r, -r*r, sint, nbRoots );

      for ( int i = 0; i < nbRoots; i++ )
      {
        if ( sint[i] * r >= 0. ) // the min distance, not the max!
        {
          double cs = sgn( q ) * sqrt( 1. - sint[i] * sint[i] );
          Point3 const& ptA = ptCenter + rA * cs * u2 + rA * sint[i] * v2;
          if ( ( fabs( ptA[Z] ) < hB ) && ( normXY( ptA ) < rB * rB ) )
            return ( true );
        }
      }
    }
  }

  return ( false );
}




// // ----------------------------------------------------------------------------
// // Returns whether two OBCs are in contact
// bool isContactBVolume( OBC const& obcA,
//                        OBC const& obcB,
//                        Transform const& a2w,
//                        Transform const& b2w )
// {
//   // Variables
//   double const rA = obcA.getRadius();
//   double const rB = obcB.getRadius();
//   double const hA = obcA.getHeight() / 2.; // half-height
//   double const hB = obcB.getHeight() / 2.; // half-height
//   Vector3 const& e_A = a2w.getBasis() * obcA.getInitOrientation();
//   Vector3 const& e_B = b2w.getBasis() * obcB.getInitOrientation();

//   // Relative positions - B w.r.t. A
//   Matrix rotMat;
//   rotateVec2VecZ( e_A, rotMat );
//   Point3 const& x = rotMat * ( *(b2w.getOrigin()) - *(a2w.getOrigin()) );
//   Vector3 const& e = ( rotMat * e_B ).normalized();

//   // General variables to use later
//   Vector3 const& u = ( Vector3( e[Y], -e[X], 0. ) ).normalized();
//   Vector3 const& v = ( u ^ e ).normalized();

//   //* Contact Scenarios: A primary, B secondary *//

//   // Face A - Face B
//   if ( ( fabs( fabs( e[Z] ) - 1. ) < tol ) &&
//        ( fabs( x[Z] ) < hA + hB ) &&
//        ( normXY( x ) < ( rA + rB ) * ( rA + rB ) ) )
//     return ( true );

//   // Band A - Band B
//   if ( fabs( x * ( - u1 ) ) < rA + rB )
//   {
//     double lBStar = ( e[Z] * x[Z] - e * x ) / ( 1. - e[Z] * e[Z] );
//     if ( fabs( lBStar ) < hB )
//     {
//       double lAStar = x[Z] + lBStar * e[Z];
//       if ( fabs( lAStar ) < hA )
//         return ( true );
//     }
//   }

//   // Face A - Edge B
//   if ( fabs( x[Z] ) > hA )
//   {
//     // Vector3 r = ( e ^ ( e ^ zAxis ) ).normalized();
//     // r = sgn( x[Z] ) * ( rB * r - hB * sgn( e[Z] ) * e );
//     Vector3 r = sgn( -x[Z] ) * ( rB * v + hB * sgn( e[Z] ) * e );
//     Point3 const& ptE = x + r;
//     if ( ( fabs( ptE[Z] ) < hA ) && ( normXY( ptE ) < rA * rA ) )
//       return ( true );
//   }

//   // Face B - Edge A
//   if ( fabs( x[Z] ) > hA )
//   {
//     // Vector3 r = ( e ^ ( e ^ zAxis ) ).normalized();
//     // r = sgn( x[Z] ) * ( rB * r - hB * sgn( e[Z] ) * e );
//     Vector3 r = sgn( -x[Z] ) * ( rB * v + hB * sgn( e[Z] ) * e );
//     Point3 const& ptE = x + r;
//     if ( ( fabs( ptE[Z] ) < hA ) && ( normXY( ptE ) < rA * rA ) )
//       return ( true );
//   }



//   // Edge A - Edge B
//   {
//     // To find the correct edge on the secondary cylinder
//     Point3 c1 = x_B2A + hB * e_B2A;
//     Point3 c2 = x_B2A - hB * e_B2A;
//     double distBand1 = pow( sqrt( normXY( c1 ) ) - rA, 2 );
//     double distBand2 = pow( sqrt( normXY( c2 ) ) - rA, 2 );
//     double distEdge1 = distBand1 + pow( fabs( c1[Z] ) - hA, 2 );
//     double distEdge2 = distBand2 + pow( fabs( c2[Z] ) - hA, 2 );
//     if ( ( distBand1 - distBand2 ) * ( distEdge1 - distEdge2 ) < 0. &&
//          ( distEdge1 < distEdge2 ? distEdge1 : distEdge2 ) >= rB * rB )
//     {
//       distEdge1 = distBand1;
//       distEdge2 = distBand2;
//     }
//     Point3 const& ptCenter = distEdge1 < distEdge2 ? c1 : c2;

//     if ( ( distBand1 < distBand2 ? distBand1 : distBand2 ) < rB * rB  &&
//           fabs( ptCenter[Z] ) < hA + rB )
//     {
//       // Misc variables
//       // double const p = rB * rB * ( normXY( v1 ) - normXY( u1 ) ) / 2.;
//       double const p = rB * rB * ( normXY( v1 ) - 1. ) / 2.;
//       double const q = rB * dotXY( u1, ptCenter ) / p;
//       double const r = rB * dotXY( v1, ptCenter ) / p;
//       // double const s = ( rA*rA - normXY( ptCenter ) - rB*rB*normXY(u1) ) / p;
//       double const s = ( rA * rA - normXY( ptCenter ) - rB * rB ) / p;

//       double sint[4];
//       int nbRoots;
//       solveQuartic( 2.*r, q*q + r*r - s, -r*s, s*s/4. - q*q, sint, nbRoots );

//       for ( int i = 0; i < nbRoots; i++ )
//       {
//         if ( fabs( sint[i] ) <= 1. )
//         {
//           double cost = ( s/2. - r*sint[i] - sint[i]*sint[i] ) / q;
//           double zVal = ptCenter[Z] + rB * cost * u1[Z] + rB * sint[i] * v1[Z];
//           if ( fabs( zVal ) < hA )
//             return ( true );
//         }
//       }
//     }
//   }

  


//   // Relative positions - A w.r.t. B
//   Matrix rotMatB;
//   rotateVec2VecZ( e_B, rotMatB );
//   Point3 const& x_A2B = rotMatB * ( *(a2w.getOrigin()) - *(b2w.getOrigin()) );
//   Vector3 const& e_A2B = ( rotMatB * e_A ).normalized();

//   // Vector3 const& u2 = ( e_A2B ^ Vector3( 0., 0., 1. ) ).normalized();
//   Vector3 const& u2 = ( Vector3( e_A2B[Y], -e_A2B[X], 0. ) ).normalized();
//   Vector3 const& v2 = ( u2 ^ e_A2B ).normalized();

//   //* Contact Scenarios: B secondary, A pirmary *//

//   // Edge B - Face A
//   if ( fabs( x_A2B[Z] ) > hB )
//   {
//     // Vector3 r = ( e_A2B ^ ( e_A2B ^ zAxis ) ).normalized();
//     // r = sgn( x_A2B[Z] ) * ( rA * r - hA * sgn( e_A2B[Z] ) * e_A2B );
//     Vector3 r = sgn( -x_A2B[Z] ) * ( rA * v2 + hA * sgn( e_A2B[Z] ) * e_A2B );
//     Point3 ptE = x_A2B + r;
//     if ( ( fabs( ptE[Z] ) < hB ) && ( normXY( ptE ) < rB * rB ) )
//       return ( true );
//   }

//   // Edge B - Band A
//   {
//     Point3 c1 = x_A2B + hA * e_A2B;
//     Point3 c2 = x_A2B - hA * e_A2B;
//     Point3 const& ptCenter = normXY( c1 ) < normXY( c2 ) ? c1 : c2;

//     // additional condition to avoid solving the polynomial
//     if ( ( fabs( ptCenter[Z] ) < hB + rA ) &&
//          ( sqrt( normXY( ptCenter ) ) < rA + rB ) )
//     {
//       // Misc variables
//       // double const p = rA * ( normXY( v2 ) - normXY( u2 ) );
//       double const p = rA * ( normXY( v2 ) - 1. );
//       double const q = dotXY( u2, ptCenter ) / p;
//       double const r = dotXY( v2, ptCenter ) / p;

//       double sint[4];
//       int nbRoots = 0;
//       solveQuartic( 2.*r, q*q + r*r - 1., -2.*r, -r*r, sint, nbRoots );

//       for ( int i = 0; i < nbRoots; i++ )
//       {
//         if ( sint[i] * r >= 0. ) // the min distance, not the max!
//         {
//           double cs = sgn( q ) * sqrt( 1. - sint[i] * sint[i] );
//           Point3 const& ptA = ptCenter + rA * cs * u2 + rA * sint[i] * v2;
//           if ( ( fabs( ptA[Z] ) < hB ) && ( normXY( ptA ) < rB * rB ) )
//             return ( true );
//         }
//       }
//     }
//   }

//   return ( false );
// }