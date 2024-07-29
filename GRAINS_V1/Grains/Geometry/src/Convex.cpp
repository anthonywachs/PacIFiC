#include "Convex.hh"
#include "Basic.hh"
#include "BBox.hh"
#include "Sphere.hh"
#include "Vector3.hh"
#include "Transform.hh"
#include "BVolume.hh"
#include "GrainsExec.hh"


// Tolerances are commented out, as we can dynamically set them through the XML
// file 
// double rel_error = EPSILON;   // relative error in the computed distance
// double abs_error = EPSILON2;  // absolute error if the distance is almost zero
int numIterations = 0;

static Point3 p[4];         // support points of object A in local coordinates
static Point3 q[4];         // support points of object B in local coordinates
static Vector3 y[4];       // support points of A - B in world coordinates

static int bits;           // identifies current simplex
static int last;           // identifies last found support point
static int last_bit;       // last_bit = 1<<last
static int all_bits;       // all_bits = bits|last_bit

static double det[16][4];  // cached sub-determinants


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
  cout << "Warning for this Convex (" 
       << this->getConvexType() 
       << ") method Convex::computeBVolume() "
       << "is not yet implemented !\n"
       << "Need for an assistance ! Stop running !\n";
  exit(10);

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




// ----------------------------------------------------------------------------
// SOLID Software low-level routines
// ----------------------------------------------------------------------------
void compute_det()
{
  static double dp[4][4];

  for (int i = 0, bit = 1; i < 4; ++i, bit <<=1)
    if (bits & bit) dp[i][last] = dp[last][i] = y[i] * y[last];
  dp[last][last] = y[last] * y[last];

  det[last_bit][last] = 1;
  for (int j = 0, sj = 1; j < 4; ++j, sj <<= 1) {
    if (bits & sj) {
      int s2 = sj|last_bit;
      det[s2][j] = dp[last][last] - dp[last][j];
      det[s2][last] = dp[j][j] - dp[j][last];
      for (int k = 0, sk = 1; k < j; ++k, sk <<= 1) {
	if (bits & sk) {
	  int s3 = sk|s2;
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
  if (all_bits == 15) {
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




// ----------------------------------------------------------------------------
inline bool valid( int s )
{
  for (int i = 0, bit = 1; i < 4; ++i, bit <<= 1) {
    if (all_bits & bit) {
      if (s & bit) {
	if (det[s][i] <= EPSILON3)
	  return false;
      } else if (det[s|bit][i] > 0.)
	return false;
    }
  }

  return ( true );
}




// ----------------------------------------------------------------------------
inline void compute_vector( int bits_, Vector3& v )
{
  double sum = 0.;
  v.setValue(0., 0., 0.);
  for (int i = 0, bit = 1; i < 4; ++i, bit <<= 1) {
    if (bits_ & bit) {
      sum += det[bits_][i];
      v += y[i] * det[bits_][i];
    }
  }
  v *= 1 / sum;
}




// ----------------------------------------------------------------------------
inline void compute_points( int bits_, Point3& p1, Point3& p2 )
{
  double sum = 0;
  p1.setValue(0, 0, 0);
  p2.setValue(0, 0, 0);
  for (int i = 0, bit = 1; i < 4; ++i, bit <<= 1) {
    if (bits_ & bit) {
      sum += det[bits_][i];
      p1 += p[i] * det[bits_][i];
      p2 += q[i] * det[bits_][i];
    }
  }
  double s = 1 / sum;
  p1 *= s;
  p2 *= s;
}




// ----------------------------------------------------------------------------
inline bool proper( int s )
{
  for (int i = 0, bit = 1; i < 4; ++i, bit <<= 1)
    if ((s & bit) && det[s][i] <= EPSILON3) return ( false );

  return ( true );
}




// ----------------------------------------------------------------------------
inline bool closest( Vector3& v )
{
  int s;
  compute_det();
  for (s = bits; s; --s) {
    if ((s & bits) == s) {
      if (valid(s|last_bit)) {
	bits = s|last_bit;
 	compute_vector(bits, v);
	return ( true );
      }
    }
  }
  if (valid(last_bit)) {
    bits = last_bit;
    v = y[last];
    return ( true );
  }
  // Original GJK calls the backup procedure at this point.
  double min_dist2 = INFINITY;
  for (s = all_bits; s; --s) {
    if ((s & all_bits) == s) {
      if (proper(s)) {
	Vector3 u;
 	compute_vector(s, u);
	double dist2 = Norm2(u);
	if (dist2 < min_dist2) {
	  min_dist2 = dist2;
	  bits = s;
	  v = u;
	}
      }
    }
  }

  return ( false );
}




// ----------------------------------------------------------------------------
// The next function is used for detecting degenerate cases that cause
// termination problems due to rounding errors.
inline bool degenerate(const Vector3& w)
{
  for (int i = 0, bit = 1; i < 4; ++i, bit <<= 1)
    if ((all_bits & bit) && y[i] == w) return ( true );

  return ( false );
}




// ----------------------------------------------------------------------------
// For numIterations > 1000
void catch_me()
{
  cerr << "closest_points : Out on iteration > 1000\n";
}
// ----------------------------------------------------------------------------
// END of SOLID Software low-level routines
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// Returns whether 2 convex shapes intersect
bool intersect( Convex const& a, Convex const& b, Transform const& a2w,
	Transform const& b2w, Vector3& v  )
{
  Vector3 w;

  bits = 0;
  all_bits = 0;

  double prod;

  do {
    last = 0;
    last_bit = 1;
    while (bits & last_bit) { ++last; last_bit <<= 1; }
    w = a2w(a.support((-v) * a2w.getBasis())) -
      b2w(b.support(v * b2w.getBasis()));
    prod = v * w;
    if (prod > 0. || fabs(prod) < EPSILON2) {
      return false;
    }
    if (degenerate(w)) {
      return false;
    }
    y[last] = w;
    all_bits = bits|last_bit;
    if (!closest(v)) {
      return false;
    }
  }
  while (bits < 15 && !approxZero(v));

  return ( true );
}




// ----------------------------------------------------------------------------
// Returns whether 2 convex shapes intersect
bool intersect( Convex const& a, Convex const& b, Transform const& b2a,
	Vector3& v )
{
  Vector3 w;

  bits = 0;
  all_bits = 0;
  double prod;

  do {
    last = 0;
    last_bit = 1;
    while (bits & last_bit) { ++last; last_bit <<= 1; }
    w = a.support(-v) - b2a(b.support(v * b2a.getBasis()));
    prod = v * w;
    if (prod > 0. || fabs(prod) < EPSILON2) {
      return false;
    }
    if (degenerate(w)) {
      return false;
    }
    y[last] = w;
    all_bits = bits|last_bit;
    if (!closest(v)) {
      return false;
    }
  }
  while (bits < 15 && !approxZero(v));

  return ( true );
}




// ----------------------------------------------------------------------------
// Returns whether 2 convex shapes intersect and if they intersect
// returns an intersection point per convex shape in each convex reference frame
bool common_point( Convex const& a, Convex const& b, Transform const& a2w,
	Transform const& b2w, Vector3& v, Point3& pa, Point3& pb )
{
  Vector3 w;

  bits = 0;
  all_bits = 0;
  double prod;

  do {
    last = 0;
    last_bit = 1;
    while (bits & last_bit) { ++last; last_bit <<= 1; }
    p[last] = a.support((-v) * a2w.getBasis());
    q[last] = b.support(v * b2w.getBasis());
    w = a2w(p[last]) - b2w(q[last]);
    prod = v * w;
    if (prod > 0. || fabs(prod) < EPSILON2) {
      return false;
    }
    if (degenerate(w)) {
      return false;
    }
    y[last] = w;
    all_bits = bits|last_bit;
    if (!closest(v)) {
      return false;
    }
  }
  while (bits < 15 && !approxZero(v) ) ;
  compute_points(bits, pa, pb);

  return ( true );
}




// ----------------------------------------------------------------------------
// Returns whether 2 convex shapes intersect and if they intersect
// returns an intersection point per convex shape in each convex reference frame
bool common_point( Convex const& a, Convex const& b, Transform const& b2a,
	Vector3& v, Point3& pa, Point3& pb )
{
  Vector3 w;

  bits = 0;
  all_bits = 0;
  double prod;

  do {
    last = 0;
    last_bit = 1;
    while (bits & last_bit) { ++last; last_bit <<= 1; }
    p[last] = a.support(-v);
    q[last] = b.support(v * b2a.getBasis());
    w = p[last] - b2a(q[last]);
    prod = v * w;
    if (prod > 0. || fabs(prod) < EPSILON2) {
      return false;
    }
    if (degenerate(w)) {
      return false;
    }
    y[last] = w;
    all_bits = bits|last_bit;
    if (!closest(v)) {
      return false;
    }
  }
  while (bits < 15 && !approxZero(v) );
  compute_points(bits, pa, pb);

  return ( true );
}









// ----------------------------------------------------------------------------
// Returns the minimal distance between 2 convex shapes and a point per
// convex shape that represents the tips of the minimal distance segment
double closest_points( Convex const& a, Convex const& b, Transform const& a2w,
	Transform const& b2w, Point3& pa, Point3& pb, int& nbIter )
{
  Vector3 v = a2w(a.support(Vector3Null)) - b2w(b.support(Vector3Null));
  Vector3 w = v;
  Vector3 d = v;
  double dist = Norm(v);

  bits = 0;
  all_bits = 0;
  numIterations = 0;
  double mu = 0.;

  double momentum, oneMinusMomentum;
  bool acceleration = GrainsExec::m_colDetAcceleration;
  double relError = GrainsExec::m_colDetTolerance;
  double absError = 1.e-4 * relError;

  while (bits < 15 && dist > EPSILON2 && numIterations < 1000)
  {
    ++numIterations;
    last = 0;
    last_bit = 1;
    while (bits & last_bit) { ++last; last_bit <<= 1; }

    // Finding the suitable direction using either Nesterov or original
    // The number 8 is hard-coded. Emprically, it shows the best convergence for
    // superquadrics. For the rest of shapes, we really do not need to use
    // Nesterov as the improvemenet is marginal.
    if ( acceleration && numIterations % 8 != 0 )
    {
      momentum = numIterations / ( numIterations + 2. );
      oneMinusMomentum = 1. - momentum;
      d = momentum * d + 
          momentum * oneMinusMomentum * v +
          oneMinusMomentum * oneMinusMomentum * w;
    }
    else
      d = v;

    p[last] = a.support( ( -d ) * a2w.getBasis() );
    q[last] = b.support( (  d ) * b2w.getBasis() );
    w = a2w(p[last]) - b2w(q[last]);

    // termination criteria
    mu = dist - v * w / dist;
    if ( mu < dist * relError ||
         mu < absError )
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

    if (degenerate(w))
      break;

    
    y[last] = w;
    all_bits = bits|last_bit;

    if ( !closest( v ) )
      break;
    dist = Norm(v);
  }
  compute_points(bits, pa, pb);
  if (numIterations > 1000) 
    catch_me();
  else 
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
