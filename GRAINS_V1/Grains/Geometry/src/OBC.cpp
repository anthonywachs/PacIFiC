#include "OBC.hh"


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
static inline int sgn( T const val )
{
    return ( ( T(0) < val ) - ( val < T(0) ) );
}




// ----------------------------------------------------------------------------
// Returns the solutions to to the quartic equation x^4 + bx^3 + cx^2 + dx + e
static inline void solveQuartic( double const b, double const c, double const d,
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
// Returns whether two OBCs are in contact
bool isContactBVolume( OBC const& obcA,
                       OBC const& obcB,
                       Transform const& a2w,
                       Transform const& b2w )
{
  // Variables
  double const r1 = obcA.getRadius();
  double const r2 = obcB.getRadius();
  double const h1 = obcA.getHeight() / 2.; // half-height
  double const h2 = obcB.getHeight() / 2.; // half-height
  Vector3 const eZ = a2w.getBasis() * obcA.getInitOrientation(); // e1 = eZ
  Vector3 e2 = b2w.getBasis() * obcB.getInitOrientation();

  // Variables to represent the 2nd cyl in the 1st cyl local coordinates
  Point3 xPt( *( b2w.getOrigin() ) - *( a2w.getOrigin() ) );
  Vector3 eX( ( e2 ^ eZ ).normalized() );
  Vector3 eY( ( eZ ^ eX ).normalized() );
  // Vector3 eZ( e1 );
  const double x = eX * xPt;
  const double y = eY * xPt;
  const double z = eZ * xPt;
  const double ey = eY * e2;
  const double ez = eZ * e2;


  /* Step one: Shortest distance -- Projection onto XZ */
  if ( fabs( x ) >= r1 + r2 )
    return ( false );
  else 
  {
    const double s2 = y / ey;
    const double s1 = z + s2 * ez;
    if ( fabs( s1 ) < h1 && fabs( s2 ) < h2 )
      return ( true );
  }
  
  /* Step two: Projection onto YZ and check rectangles intersection */
  // We project the rectangles onto four axes; local coords of the
  // first and second rectangles.
  {
    const double fey = fabs( ey ) + 1.e-5;
    const double fez = fabs( ez ) + 1.e-5;
    // onto y-axis
    if ( fabs( y ) > r1 + h2 * fey + r2 * fez )
      return ( false ); 
    // onto z-axis
    if ( fabs( z ) > h1 + h2 * fez + r2 * fey )
      return ( false );
    // onto e_normal
    if ( fabs( y * ez - z * ey ) > r2 + h1 * fey + r1 * fez )
      return ( false );
    // onto e
    if ( fabs( y * ey + z * ez ) > h2 + h1 * fez + r1 * fey )
      return ( false );
  }

  /* Step three: Projection onto XY */
  // We check the overlap of an ellipse w/ a circle
  if ( fabs( ez ) < 0.99 && fabs( ez ) > 0.01 )
  {
    // First, secondary is an ellipse
    {
      const double topy = y + h2 * ey;
      const double bottomy = y - h2 * ey;
      if ( std::signbit( topy ) == std::signbit( bottomy ) )
      {
        const double x0 = x;
        const double y0 = abs( topy ) > abs( bottomy ) ? bottomy : topy;
        // check if origin is not in the ellipse -- and cyls are not par/per
        if ( x0 * x0 + y0 * y0 / ez / ez > r2 * r2 )
        {
          const double ey_sq = ey * ey;
          const double A = r2 * r2 * ey_sq * ey_sq;
          const double B = 2. * y0 * ez / r2 / ey_sq;
          const double C = x0 * x0 / A;
          const double D = y0 * y0 * ez * ez / A;

          double sint[4];
          int nbRoots;
          solveQuartic( -B, C + D - 1., B, -D, sint, nbRoots );

          for ( int i = 0; i < nbRoots; i++ )
          {
            if ( std::signbit( sint[i] ) != std::signbit( B ) )
            {
              const double cost = - sgn( x0 ) * sqrt( 1. - sint[i] * sint[i] );
              const double ptX = x0 + r2 * cost;
              const double ptY = y0 + r2 * ez * sint[i];
              if ( ptX * ptX + ptY * ptY > r1 * r1 )
                return ( false );
            }
          }
        }
      }
    } 
    // Next, primary is an ellipse
    {
      Vector3 eY2( ( eX ^ e2 ).normalized() );
      // don't update x2, it is the same as x with sign flipped
      const double y2 = -( eY2 * xPt );
      const double ey2 = eY2 * eZ;
      // don't update ez, it is the same as before
      const double topy = y2 + h1 * ey2;
      const double bottomy = y2 - h1 * ey2;
      if ( std::signbit( topy ) == std::signbit( bottomy ) )
      {
        const double x0 = -x;
        const double y0 = abs( topy ) > abs( bottomy ) ? bottomy : topy;
        if ( x0 * x0 + y0 * y0 / ez / ez > r1 * r1 )
        {
          const double ey_sq = ey2 * ey2;
          const double A = r1 * r1 * ey_sq * ey_sq;
          const double B = 2. * y0 * ez / r1 / ey_sq;
          const double C = x0 * x0 / A;
          const double D = y0 * y0 * ez * ez / A;

          double sint[4];
          int nbRoots;
          solveQuartic( -B, C + D - 1., B, -D, sint, nbRoots );

          for ( int i = 0; i < nbRoots; i++ )
          {
            if ( std::signbit( sint[i] ) != std::signbit( B ) )
            {
              const double cost = - sgn( x0 ) * sqrt( 1. - sint[i] * sint[i] );
              const double ptX = x0 + r1 * cost;
              const double ptY = y0 + r1 * ez * sint[i];
              if ( ptX * ptX + ptY * ptY > r2 * r2 )
                return ( false );
            }
          }
        }
      }
    }
  }
  
  return ( true );
}
