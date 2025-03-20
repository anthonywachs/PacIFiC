#include "Convex.hh"
#include "Basic.hh"
#include "BBox.hh"
#include "Sphere.hh"
#include "Vector3.hh"
#include "Transform.hh"
#include "BVolume.hh"
#include "GrainsExec.hh"


// ----------------------------------------------------------------------------
// Default constructor (forbidden)
Convex::Convex() : Shape()
{}




// ----------------------------------------------------------------------------
// Destructor
Convex::~Convex()
{}




// ----------------------------------------------------------------------------
// Returns the convex shape bounding box
BBox Convex::bbox( Transform const& t ) const
{
  Point3 const* ori = t.getOrigin() ;

  Point3 min( (*ori)[X] + t.getBasis()[X] * support(-t.getBasis()[X]),
	    (*ori)[Y] + t.getBasis()[Y] * support(-t.getBasis()[Y]),
	    (*ori)[Z] + t.getBasis()[Z] * support(-t.getBasis()[Z]) );
  Point3 max( (*ori)[X] + t.getBasis()[X] * support(t.getBasis()[X]),
	    (*ori)[Y] + t.getBasis()[Y] * support(t.getBasis()[Y]),
	    (*ori)[Z] + t.getBasis()[Z] * support(t.getBasis()[Z]) );

  return ( BBox( min, max ) );
}




// ----------------------------------------------------------------------------
// Returns the convex shape bounding volume
BVolume* Convex::computeBVolume( unsigned int type ) const
{
  if ( type != 0 )
  {
    cout << "Warning for this Convex (" 
         << this->getConvexType() 
         << ") method Convex::computeBVolume() "
         << "is not yet implemented!\n"
         << "Changing the bounding volume collision detection OFF!\n";
    GrainsExec::m_colDetBoundingVolume = 0;
  }
  return( nullptr );
}




// ----------------------------------------------------------------------------
// Returns a vector of points describing the envelope of the convex
vector<Point3> Convex::getEnvelope() const
{
  cout << "Warning for this Convex the method Convex::getEnvelope() "
       << "is not yet implemented !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);

  vector<Point3> envelope;
  return ( envelope );
}




// ----------------------------------------------------------------------------
// Returns a pointer to a 2D array describing the relationship
// between the face indices and the point indices
vector<vector<int> > const* Convex::getFaces() const
{
  cout << "Warning for this Convex the method Convex::getFaces() "
       << "is not yet implemented !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);

  vector<vector<int> >* allFaces = NULL;

  return ( allFaces );
}




// ----------------------------------------------------------------------------
// Returns the number of vertices/corners or a code corresponding to
// a specific convex shape having to vertices/corners as, e.g., a cylinder
int Convex::getNbCorners() const
{
  cout << "Warning for this Convex the method Convex::getNbCorners() "
       << "is not yet implemented !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);

  return ( 0 );
}




// ----------------------------------------------------------------------------
// Returns the number of points to write the convex shape in a Paraview format
int Convex::numberOfPoints_PARAVIEW() const
{
  cout << "Warning for this Convex the method Convex::numberOfPoints_PARAVIEW()"
       << " is not yet implemented !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);

  return ( 0 );
}




// ----------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the convex
// shape in a Paraview format
int Convex::numberOfCells_PARAVIEW() const
{
  cout << "Warning for this Convex the method Convex::numberOfCells_PARAVIEW()"
       << " is not yet implemented !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);

  return ( 0 );
}




// ----------------------------------------------------------------------------
// Writes the points describing the convex shape in a Paraview format
void Convex::write_polygonsPts_PARAVIEW( ostream& f,
  	Transform const& transform, Vector3 const* translation ) const
{
  cout << "Warning for this Convex the method "
       << "Convex::write_polygonsPts_PARAVIEW() is not yet implemented !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Writes the convex shape in a STL format
void Convex::write_convex_STL( ostream& f, Transform  const& transform ) const
{
  cout << "Warning for this Convex the method "
       << "Convex::write_convex_STL() is not yet implemented !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Writes a list of points describing the convex shape in a Paraview format
list<Point3> Convex::get_polygonsPts_PARAVIEW( Transform const& transform,
  	Vector3 const* translation )
	const
{
  cout << "Warning for this Convex the method "
       << "Convex::get_polygonsPts_PARAVIEW() is not yet implemented !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}



// ----------------------------------------------------------------------------
// Writes the convex shape in a Paraview format
void Convex::write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const
{
  cout << "Warning for this Convex the method "
       << "Convex::write_polygonsStr_PARAVIEW() is not yet implemented !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);
}




// ----------------------------------------------------------------------------
// Returns an orientation vector describing the convex shape angular position
Vector3 Convex::computeOrientationVector( Transform const* transform ) const
{
  return ( Vector3Null );
}




/* ========================================================================== */
/*                             Low-Level Methods                              */
/* ========================================================================== */
static inline void computeDet( unsigned int const bits,
	unsigned int const last,
	unsigned int const last_bit,
	unsigned int const all_bits,
	Vector3 const y[4],
	double dp[4][4],
	double det[16][4] )
{
  for ( unsigned int i = 0, bit = 1; i < 4; ++i, bit <<= 1 )
    if ( bits & bit ) 
      dp[i][last] = dp[last][i] = y[i] * y[last];
  dp[last][last] = y[last] * y[last];

  det[last_bit][last] = 1.;
  for ( unsigned int j = 0, sj = 1; j < 4; ++j, sj <<= 1 )
  {
    if ( bits & sj )
    {
      unsigned int s2 = sj | last_bit;
      det[s2][j] = dp[last][last] - dp[last][j];
      det[s2][last] = dp[j][j] - dp[j][last];
      for ( unsigned int k = 0, sk = 1; k < j; ++k, sk <<= 1 )
      {
        if ( bits & sk )
	{
	  int s3 = sk | s2;
	  det[s3][k] = det[s2][j] * (dp[j][j] - dp[j][k]) +
	  det[s2][last] * (dp[last][j] - dp[last][k]);
	  det[s3][j] = det[sk|last_bit][k] * (dp[k][k] - dp[k][j]) +
	  det[sk|last_bit][last] * (dp[last][k] - dp[last][j]);
	  det[s3][last] = det[sk|sj][k] * (dp[k][k] - dp[k][last]) +
	  det[sk|sj][j] * (dp[j][k] - dp[j][last]);
	}
      }
    }
  }

  if ( all_bits == 15 )
  {
    det[15][0] = det[14][1] * (dp[1][1] - dp[1][0]) +
                     det[14][2] * (dp[2][1] - dp[2][0]) +
                     det[14][3] * (dp[3][1] - dp[3][0]);
    det[15][1] = det[13][0] * (dp[0][0] - dp[0][1]) +
                     det[13][2] * (dp[2][0] - dp[2][1]) +
                     det[13][3] * (dp[3][0] - dp[3][1]);
    det[15][2] = det[11][0] * (dp[0][0] - dp[0][2]) +
                     det[11][1] * (dp[1][0] - dp[1][2]) +
                     det[11][3] * (dp[3][0] - dp[3][2]);
    det[15][3] = det[7][0] * (dp[0][0] - dp[0][3]) +
                     det[7][1] * (dp[1][0] - dp[1][3]) +
                     det[7][2] * (dp[2][0] - dp[2][3]);
  }
}




// ----------------------------------------------------------------------------- 
static inline bool valid( unsigned int const s,
	unsigned int const all_bits,
	double const det[16][4] )
{
  for ( unsigned int i = 0, bit = 1; i < 4; ++i, bit <<= 1 )
  {
    if ( all_bits & bit )
    {
      if ( s & bit )
      {
        if ( det[s][i] <= EPSILON3 )
          return ( false );
      }
      else if ( det[s|bit][i] > 0 )
        return ( false );
    }
  }
  return ( true );
}




// -----------------------------------------------------------------------------
static inline void computeVector( unsigned int const bits_,
	Vector3 const y[4],
	double const det[16][4],
	Vector3& v )
{
  double sum = 0.;
  v.setValue( 0., 0., 0. );
  for ( unsigned int i = 0, bit = 1; i < 4; ++i, bit <<= 1 )
  {
    if ( bits_ & bit )
    {
      sum += det[bits_][i];
      v += det[bits_][i] * y[i];
    }
  }
  v *= 1. / sum;
}




// -----------------------------------------------------------------------------
static inline void computePoints( unsigned int const bits_,
	double const det[16][4],
	Point3 const p[4],
	Point3 const q[4],
	Point3& p1,
	Point3& p2 )
{
  double sum = 0.;
  p1.setValue( 0., 0., 0. );
  p2.setValue( 0., 0., 0. );
  for ( unsigned int i = 0, bit = 1; i < 4; ++i, bit <<= 1 )
  {
    if ( bits_ & bit )
    {
      sum += det[bits_][i];
      p1 += det[bits_][i] * p[i];
      p2 += det[bits_][i] * q[i];
    }
  }
  double s = 1. / sum;
  p1 *= s;
  p2 *= s;
}




// -----------------------------------------------------------------------------
static inline bool proper( unsigned int const s,
	double const det[16][4] )
{
  for ( unsigned int i = 0, bit = 1; i < 4; ++i, bit <<= 1 )
    if ( ( s & bit ) && det[s][i] <= EPSILON3 )
      return ( false );
  return ( true );
}




// -----------------------------------------------------------------------------
static inline bool closest( unsigned int& bits,
	unsigned int const last,
	unsigned int const last_bit,
	unsigned int const all_bits,
	Vector3 const y[4],
	double dp[4][4],
	double det[16][4],
	Vector3& v )
{
  unsigned int s;
  computeDet( bits, last, last_bit, all_bits, y, dp, det );
  for ( s = bits; s; --s )
  {
    if ( ( s & bits ) == s )
    {
      if ( valid( s | last_bit, all_bits, det ) )
      {
        bits = s | last_bit;
	computeVector( bits, y, det, v );
	return( true );
      }
    }
  }
  if ( valid( last_bit, all_bits, det ) )
  {
    bits = last_bit;
    v = y[last];
    return( true );
  }
  // Original GJK calls the backup procedure at this point.
  double min_dist2 = INFINITY;
  for ( s = all_bits; s; --s )
  {
    if ( ( s & all_bits ) == s )
    {
      if ( proper( s, det ) )
      {
        Vector3 u;
	computeVector( s, y, det, u );
	double dist2 = Norm2( u );
	if ( dist2 < min_dist2 )
	{
	  min_dist2 = dist2;
	  bits = s;
	  v = u;
	}
      }
    }
  }
  return ( false );
}




// -----------------------------------------------------------------------------
// The next function is used for detecting degenerate cases that cause
// termination problems due to rounding errors.
static inline bool degenerate( unsigned int const all_bits,
	Vector3 const y[4],
	Vector3 const& w )
{
  for ( unsigned int i = 0, bit = 1; i < 4; ++i, bit <<= 1 )
    if ( ( all_bits & bit ) && y[i] == w )
      return ( true );
  return ( false );
}




// -----------------------------------------------------------------------------
// For num_iterations > 1000
void catch_me()
{
  cerr << "closest_points : Out on iteration > 1000\n";
}




/* ========================================================================== */
/*                            High-Level Methods                              */
/* ========================================================================== */
// Returns whether 2 convex shapes intersect using the GJK algorithm
bool intersect( Convex const& a, Convex const& b,
	Transform const& a2w, Transform const& b2w, Vector3& v )
{
  unsigned int bits = 0;           // identifies current simplex
  unsigned int last = 0;           // identifies last found support point
  unsigned int last_bit = 0;       // last_bit = 1<<last
  unsigned int all_bits = 0;       // all_bits = bits|last_bit
  Vector3 y[4];                    // support points of A-B in world
  double det[16][4] = { 0. };      // cached sub-determinants
  double dp[4][4] = { 0. };

  Vector3 w;
  double prod;
    
  do 
  {        
    last = 0;
    last_bit = 1;
    while( bits & last_bit )
    {
      ++last;
      last_bit <<= 1;
    }
    w = a2w( a.support( ( -v ) * a2w.getBasis() ) ) -
		b2w( b.support(    v   * b2w.getBasis() ) );
    prod = v * w;
    if ( prod > 0. || fabs( prod ) < EPSILON2 )
      return ( false );
    if ( degenerate( all_bits, y, w ) )
      return ( false );      
    y[last] = w;
    all_bits = bits | last_bit;
    if ( !closest( bits, last, last_bit, all_bits, y, dp, det, v ) )
      return ( false );
  } while ( bits < 15 && !approxZero( v ) );
  return ( true );
}




// -----------------------------------------------------------------------------
// Returns whether 2 convex shapes intersect using the GJK algorithm - relative
// transformation
bool intersect( Convex const& a, Convex const& b,
	Transform const& b2a, Vector3& v )
{
  unsigned int bits = 0;           // identifies current simplex
  unsigned int last = 0;           // identifies last found support point
  unsigned int last_bit = 0;       // last_bit = 1<<last
  unsigned int all_bits = 0;       // all_bits = bits|last_bit
  Vector3 y[4];                    // support points of A-B in world
  double det[16][4] = {0.};       // cached sub-determinants
  double dp[4][4] = {0.};

  Vector3 w;
  double prod;
    
  do 
  {        
    last = 0;
    last_bit = 1;
    while( bits & last_bit )
    {
      ++last;
      last_bit <<= 1;
    }
    w = a.support( -v ) - b2a( b.support( v * b2a.getBasis() ) );
    prod = v * w;
    if ( prod > 0. || fabs( prod ) < EPSILON2 )
      return ( false );
    if ( degenerate( all_bits, y, w ) )
      return ( false );      
    y[last] = w;
    all_bits = bits | last_bit;
    if ( !closest( bits, last, last_bit, all_bits, y, dp, det, v ) )
      return ( false );
  } while ( bits < 15 && !approxZero( v ) );
  return ( true );
}




// -----------------------------------------------------------------------------
// Returns the minimal distance between 2 convex shapes and a point per convex
// shape that represents the tips of the minimal distance segment
double closest_points( Convex const& a, Convex const& b, 
	Transform const& a2w, Transform const& b2w, 
	Point3& pa, Point3& pb, int& nbIter )
{
  // GJK variables
  unsigned int bits = 0;           // identifies current simplex
  unsigned int last = 0;           // identifies last found support point
  unsigned int last_bit = 0;       // last_bit = 1<<last
  unsigned int all_bits = 0;       // all_bits = bits|last_bit
  Point3 p[4];                     // support points of A in local
  Point3 q[4];                     // support points of B in local
  Vector3 y[4];                    // support points of A-B in world
  double mu = 0.;                  // optimality gap
  int numIterations = 0;           // No. iterations
  double det[16][4] = { 0. };      // cached sub-determinants
  double dp[4][4] = { 0. };

  // Misc variables, e.g. tolerance, ...
  double relError = GrainsExec::m_colDetTolerance; // rel error for opt gap
  double absError = 1.e-4 * relError; // abs error for optimality gap
  bool acceleration = GrainsExec::m_colDetAcceleration; // isAcceleration?
  double momentum = 0., oneMinusMomentum = 1.; // in case we use acceleration

  // Initializing vectors
  Vector3 v( a2w( a.support( Vector3Null ) ) 
  	- b2w( b.support( Vector3Null ) ) );
  Vector3 w( v );
  Vector3 d( v );
  double dist = Norm( v );

  while ( bits < 15 && dist > EPSILON2 && numIterations < 1000 )
  {
    ++numIterations;
    // Updating the bits, ...
    last = 0;
    last_bit = 1;
    while (bits & last_bit) { ++last; last_bit <<= 1; }

    // Finding the suitable direction using either Nesterov or original
    // The number 8 is hard-coded. Emprically, it shows the best convergence
    // for superquadrics. For the rest of shapes, we really do not need to 
    // use Nesterov as the improvemenet is marginal.
    if ( acceleration && numIterations % 8 != 0 )
    {
      momentum = numIterations / ( numIterations + 2. );
      oneMinusMomentum = 1. - momentum;
      d = momentum * d + momentum * oneMinusMomentum * v +
		oneMinusMomentum * oneMinusMomentum * w;
    }
    else
      d = v;

    p[last] = a.support( ( -d ) * a2w.getBasis() );
    q[last] = b.support( (  d ) * b2w.getBasis() );
    w = a2w( p[last] ) - b2w( q[last] );
        
    // termination criteria -- optimiality gap
    mu = dist - v * w / dist;
    if ( mu < dist * relError || mu < absError )
    {
      if ( acceleration )
      {
        if ( Norm( d - v ) < EPSILON )
          break;
        acceleration = false;
	p[last] = a.support( ( -v ) * a2w.getBasis() );
	q[last] = b.support( (  v ) * b2w.getBasis() );
	w = a2w( p[last] ) - b2w( q[last] );
      }
      else
        break;
    }
    // termination criteria -- degenerate case
    if ( degenerate( all_bits, y, w ) )
      break;
        
    // if not terminated, get ready for the next iteration
    y[last] = w;
    all_bits = bits | last_bit;
    if ( !closest( bits, last, last_bit, all_bits, y, dp, det, v ) )
      break;
    dist = Norm( v );
  }
  // compute witness points
  computePoints( bits, det, p, q, pa, pb );
  // if has not converged after lot of iterations
  if ( numIterations > 1000 ) 
    catch_me();
  else // otherwise, report the No. iterations
    nbIter = numIterations;
  return ( dist );
}




// -----------------------------------------------------------------------------
// Returns the minimal distance between 2 convex shapes and a point per convex
// shape that represents the tips of the minimal distance segment -- initial
// guess
double closest_points( Convex const& a, Convex const& b, 
	Transform const& a2w, Transform const& b2w, Vector3& v, 
	Point3& pa, Point3& pb, int& nbIter )
{
  // GJK variables
  unsigned int bits = 0;           // identifies current simplex
  unsigned int last = 0;           // identifies last found support point
  unsigned int last_bit = 0;       // last_bit = 1<<last
  unsigned int all_bits = 0;       // all_bits = bits|last_bit
  Point3 p[4];                     // support points of A in local
  Point3 q[4];                     // support points of B in local
  Vector3 y[4];                    // support points of A-B in world
  double det[16][4] = { 0. };      // cached sub-determinants
  double dp[4][4] = { 0. };
  double dist = 0.;                // distance

  // Misc variables, e.g. tolerance, ...
  double relError = GrainsExec::m_colDetTolerance; // rel error for opt gap
  double absError = 1.e-4 * relError; // abs error for optimality gap
  bool acceleration = GrainsExec::m_colDetAcceleration; // isAcceleration?
  double momentum = 0., oneMinusMomentum = 1.; // in case we use acceleration
  int numIterations = 0;          // No. iterations
  double mu = 0.;                 // optimality gap

  // Initializing vectors
  Vector3 w;
  if ( v != Vector3Null )
    w = a2w( a.support( ( -v ) * a2w.getBasis() ) ) - 
		b2w( b.support( (  v ) * b2w.getBasis() ) );
  else // if we don't have a better guess
  {
    v = a2w( a.support( Vector3Null ) ) - b2w( b.support( Vector3Null ) );
    w = v;
  }
  Vector3 d( v );
  dist = 1.;

  while ( bits < 15 && dist > EPSILON2 && numIterations < 1000 )
  {
    ++numIterations;
    last = 0;
    last_bit = 1;
    while (bits & last_bit) { ++last; last_bit <<= 1; }

    // Finding the suitable direction using either Nesterov or original
    // The number 8 is hard-coded. Emprically, it shows the best convergence
    // for superquadrics. For the rest of shapes, we really do not need to 
    // use Nesterov as the improvemenet is marginal.
    if ( acceleration && numIterations % 8 != 0 )
    {
      momentum = numIterations / ( numIterations + 2. );
      oneMinusMomentum = 1. - momentum;
      d = momentum * d + momentum * oneMinusMomentum * v +
		oneMinusMomentum * oneMinusMomentum * w;
    }
    else
      d = v;

    p[last] = a.support( ( -d ) * a2w.getBasis() );
    q[last] = b.support( (  d ) * b2w.getBasis() );
    w = a2w( p[last] ) - b2w( q[last] );
        
    // termination criteria -- optimiality gap
    mu = dist - v * w / dist;
    if ( mu < dist * relError || mu < absError )
    {
      if ( acceleration )
      {
        if ( Norm( d - v ) < EPSILON )
	  break;
        acceleration = false;
        p[last] = a.support( ( -v ) * a2w.getBasis() );
        q[last] = b.support( (  v ) * b2w.getBasis() );
        w = a2w( p[last] ) - b2w( q[last] );
      }
      else
        break;
    }
    // termination criteria -- degenerate case
    if ( degenerate( all_bits, y, w ) )
      break;
        
    // if not terminated, get ready for the next iteration
    y[last] = w;
    all_bits = bits | last_bit;
    if ( !closest( bits, last, last_bit, all_bits, y, dp, det, v ) )
      break;
    dist = Norm( v );
  }
  // compute witness points
  computePoints( bits, det, p, q, pa, pb );
  // if has not converged after lot of iterations
  if ( numIterations > 1000 ) 
    catch_me();
  else // otherwise, report the No. iterations
    nbIter = numIterations;
  return ( dist );
}




// ---------------------------------------------------------------------
// Output operator
ostream& operator << ( ostream& fileOut, Convex const& convex )
{
  convex.writeShape( fileOut );

  return ( fileOut );
}




// ----------------------------------------------------------------------------
// Input operator
istream& operator >> ( istream& fileIn, Convex& convex )
{
  convex.readShape( fileIn );

  return ( fileIn );
}




// ----------------------------------------------------------------------------
// Returns shape type
ShapeType Convex::getType() const
{
  return ( CONVEX );
}




// ----------------------------------------------------------------------------
// Returns whether the convex shape is a sphere
bool Convex::isSphere() const
{
  return ( false );
}




// ----------------------------------------------------------------------------
// Returns whether two convexes are of the same type 
bool Convex::equalType( Convex const* other, bool const& level2 ) const
{
  bool same = ( getConvexType() == other->getConvexType() );

  if ( same && level2 ) same = equalType_level2( other );
  
  return ( same );
}
