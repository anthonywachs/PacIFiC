#include "Cylinder.hh"
#include "PointContact.hh"
#include "BVolume.hh"
#include "OBB.hh"
#include "OBC.hh"
#include "Matrix.hh"

int Cylinder::m_visuNodeNbOnPer = 32;
static double tol = EPSILON; // Tolerance used in this class

// ----------------------------------------------------------------------------
// Constructor with radius and height as input parameters
Cylinder::Cylinder( double r, double h )
  : m_radius( r )
  , m_halfHeight( h / 2. )
{}




// ----------------------------------------------------------------------------
// Constructor with an input stream
Cylinder::Cylinder( istream &fileIn )
{
  readShape( fileIn );
}




// -----------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
Cylinder::Cylinder( DOMNode* root )
{
  m_radius = ReaderXML::getNodeAttr_Double( root, "Radius" );
  m_halfHeight = ReaderXML::getNodeAttr_Double( root, "Height") / 2.;
}




// ----------------------------------------------------------------------------
// Destructor
Cylinder::~Cylinder()
{}




// ----------------------------------------------------------------------
// Returns the convex type
ConvexType Cylinder::getConvexType() const
{
  return ( CYLINDER );
}




// ----------------------------------------------------------------------------
// Computes the inertia tensor and the inverse of the inertia tensor
bool Cylinder::BuildInertia( double* inertia, double* inertia_1 ) const
{
  inertia[1] = inertia[2] = inertia[4] = 0.0;
  double constant = 0.5 * m_halfHeight * m_radius * m_radius * PI;
  inertia[0] = inertia[5] = constant
  	* ( 4.0 * m_halfHeight * m_halfHeight / 3.0 + m_radius * m_radius );
  inertia[3] = 2.0 * constant * m_radius * m_radius;

  inertia_1[1] = inertia_1[2] = inertia_1[4] = 0.0;
  inertia_1[5] = inertia_1[0] = 1.0 / inertia[0];
  inertia_1[3] = 1.0 / inertia[3];
  return true;
}




// ----------------------------------------------------------------------------
// Returns the circumscribed radius of the reference sphere,
// i.e., without applying any transformation
double Cylinder::computeCircumscribedRadius() const
{
  return ( sqrt( m_radius * m_radius + m_halfHeight * m_halfHeight ) );
}




// ----------------------------------------------------------------------------
// Returns a clone of the cylinder
Convex* Cylinder::clone() const
{
  return ( new Cylinder( m_radius, 2. * m_halfHeight ) );
}




// ----------------------------------------------------------------------------
// Returns the cylinder volume
double Cylinder::getVolume() const
{
  return ( 2 * m_halfHeight * PI * m_radius * m_radius );
}




// ----------------------------------------------------------------------------
// Cylinder support function, returns the support point P, i.e. the
// point on the surface of the cylinder that satisfies max(P.v)
Point3 Cylinder::support( Vector3 const& v ) const
{
  double norm = Norm( v );
  if ( norm > EPSILON )
  {
    double s = sqrt( v[X] * v[X] + v[Z] * v[Z] );
    if ( s > EPSILON )
    {
      double d = m_radius / s;
      return ( Point3( v[X] * d, v[Y] < 0. ? - m_halfHeight : m_halfHeight,
      	v[Z] * d ) );
    }
    else
      return ( Point3( 0., v[Y] < 0. ? - m_halfHeight : m_halfHeight, 0. ) );
  }
  else
    return ( Point3() );
}




// ----------------------------------------------------------------------------
// Returns a vector of points describing the envelope of the
// cylinder. Here simply returns 3 points as follows: center of bottom circular
// face, center of top circular face and an arbitrary point on the lateral
// surface of the cylinder
vector<Point3> Cylinder::getEnvelope() const
{
  Point3 point(0.,0.,0.);
  vector<Point3> enveloppe(3,point);
  enveloppe[0][Y] = - m_halfHeight;
  enveloppe[1][Y] = - m_halfHeight;
  enveloppe[1][X] = m_radius;
  enveloppe[2][Y] = m_halfHeight;
  return ( enveloppe );
}




// ----------------------------------------------------------------------------
// Returns the number of vertices/corners or a code corresponding to
// a specific convex shape. Here returns the code 777
int Cylinder::getNbCorners() const
{
  return ( 777 );
}




// ----------------------------------------------------------------------------
// Returns a pointer to a 2D array describing the relationship between the face
// indices and the point indices. Returns a null pointer as a convention
vector< vector<int> > const* Cylinder::getFaces() const
{
  vector< vector<int> >* allFaces = NULL;
  return ( allFaces );
}




// ----------------------------------------------------------------------------
// Output operator
void Cylinder::writeShape( ostream& fileOut ) const
{
  fileOut << "*Cylinder " << m_radius << " " << 2.0 * m_halfHeight << " *END";
}



// ----------------------------------------------------------------------------
// Input operator
void Cylinder::readShape( istream& fileIn )
{
  fileIn >> m_radius >> m_halfHeight;
  m_halfHeight /= 2.0;
}




// ----------------------------------------------------------------------------
// Returns the number of points to write the cylinder in a Paraview format
int Cylinder::numberOfPoints_PARAVIEW() const
{
  return ( 2 * m_visuNodeNbOnPer + 2 );
}




// ----------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the cylinder in a
// Paraview format
int Cylinder::numberOfCells_PARAVIEW() const
{
  return ( m_visuNodeNbOnPer );
}




// ----------------------------------------------------------------------------
// Writes a list of points describing the cylinder in a Paraview format
void Cylinder::write_polygonsPts_PARAVIEW( ostream& f,
  	Transform const& transform, Vector3 const* translation ) const
{
  Point3 pp, p;
  double dtheta = 2.* PI / m_visuNodeNbOnPer;

  // Lower disk rim
  p[Y] = - m_halfHeight;
  for (int i=0;i<m_visuNodeNbOnPer;++i)
  {
    p[X] = m_radius * cos ( i * dtheta );
    p[Z] = m_radius * sin ( i * dtheta );
    pp = transform( p );
    if ( translation ) pp += *translation;
    f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
  }

  // Upper disk rim
  p[Y] = m_halfHeight;
  for (int i=0;i<m_visuNodeNbOnPer;++i)
  {
    p[X] = m_radius * cos ( i * dtheta );
    p[Z] = m_radius * sin ( i * dtheta );
    pp = transform( p );
    if ( translation ) pp += *translation;
    f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
  }

  // Lower disk center
  p[X] = 0.;
  p[Y] = - m_halfHeight;
  p[Z] = 0.;
  pp = transform( p );
  if ( translation ) pp += *translation;
  f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;

  // Upper disk center
  p[Y] = m_halfHeight;
  pp = transform( p );
  if ( translation ) pp += *translation;
  f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
}




// ----------------------------------------------------------------------------
// Returns a list of points describing the cylinder in a Paraview format
list<Point3> Cylinder::get_polygonsPts_PARAVIEW( Transform const& transform,
  	Vector3 const* translation ) const
{
  list<Point3> ParaviewPoints;
  Point3 pp,p;
  double dtheta = 2.* PI / m_visuNodeNbOnPer;

  // Lower disk rim
  p[Y] = - m_halfHeight;
  for (int i=0;i<m_visuNodeNbOnPer;++i)
  {
    p[X] = m_radius * cos ( i * dtheta );
    p[Z] = m_radius * sin ( i * dtheta );
    pp = transform( p );
    if ( translation ) pp += *translation;
    ParaviewPoints.push_back( pp );
  }

  // Upper disk rim
  p[Y] = m_halfHeight;
  for (int i=0;i<m_visuNodeNbOnPer;++i)
  {
    p[X] = m_radius * cos ( i * dtheta );
    p[Z] = m_radius * sin ( i * dtheta );
    pp = transform( p );
    if ( translation ) pp += *translation;
    ParaviewPoints.push_back( pp );
  }

  // Lower disk center
  p[X] = 0.;
  p[Y] = - m_halfHeight;
  p[Z] = 0.;
  pp = transform( p );
  if ( translation ) pp += *translation;
  ParaviewPoints.push_back( pp );

  // Upper disk center
  p[Y] = m_halfHeight;
  pp = transform( p );
  if ( translation ) pp += *translation;
  ParaviewPoints.push_back( pp );

  return ( ParaviewPoints );
}




// ----------------------------------------------------------------------------
// Writes the cylinder in a Paraview format
void Cylinder::write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const
{
  for (int i=0;i<m_visuNodeNbOnPer-1;++i)
  {
    connectivity.push_back( firstpoint_globalnumber + i );
    connectivity.push_back( firstpoint_globalnumber + i + 1 );
    connectivity.push_back( firstpoint_globalnumber + 2 * m_visuNodeNbOnPer );
    connectivity.push_back( firstpoint_globalnumber + i + m_visuNodeNbOnPer);
    connectivity.push_back( firstpoint_globalnumber + i + m_visuNodeNbOnPer
    	+ 1 );
    connectivity.push_back( firstpoint_globalnumber + 2 * m_visuNodeNbOnPer
    	+ 1 );
    last_offset += 6;
    offsets.push_back( last_offset );
    cellstype.push_back( 13 );
  }
  connectivity.push_back( firstpoint_globalnumber + m_visuNodeNbOnPer - 1 );
  connectivity.push_back( firstpoint_globalnumber );
  connectivity.push_back( firstpoint_globalnumber + 2 * m_visuNodeNbOnPer );
  connectivity.push_back( firstpoint_globalnumber + 2 * m_visuNodeNbOnPer - 1 );
  connectivity.push_back( firstpoint_globalnumber + m_visuNodeNbOnPer );
  connectivity.push_back( firstpoint_globalnumber + 2 * m_visuNodeNbOnPer + 1 );
  last_offset += 6;
  offsets.push_back( last_offset );
  cellstype.push_back( 13 );

  firstpoint_globalnumber += 2 * m_visuNodeNbOnPer + 2;
}




// ----------------------------------------------------------------------------
// Returns whether a point lies inside the cylinder
bool Cylinder::isIn( Point3 const& pt ) const
{
  return ( pt[Y] >= - m_halfHeight && pt[Y] <= m_halfHeight
  	&& sqrt( pt[X] * pt[X] + pt[Z] * pt[Z] ) <= m_radius );
}




// ----------------------------------------------------------------------------
// Returns the bounding volume to box
BVolume* Cylinder::computeBVolume( unsigned int type ) const
{
  BVolume* bvol = NULL;
  if ( type == 1 ) // OBB
  {
    Vector3 const& extent = Vector3( m_radius, m_halfHeight, m_radius );
    bvol = new OBB( extent, Matrix() );
  }

  else if ( type == 2 ) // OBC
  {
    Vector3 const& e = Vector3( 0., 1., 0. );
    bvol = new OBC( m_radius, 2*m_halfHeight, e );
  }
  
  return( bvol );
}




// ----------------------------------------------------------------------------
// LOW-LEVEL ROUTINES FOR CYLINDERS CONTACTS
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
// Returns the contact point of two cylinders in the world of the first cylinder
// if the contact is either Face-Face or Band-Band (Parallel)
inline void F2FB2BParContact( double rA, double hA, double rB, double hB,
                       Vector3 const& e, Point3 const& x, PointContact& ptCont )
{
  // Contact happens if
  bool contCond = false; // TRUE if contact occurs in this scenario
  contCond = ( fabs( fabs( e[Z] ) - 1. ) < tol ) &&
             ( fabs( x[Z] ) < hA + hB ) &&
             ( normXY( x ) < ( rA + rB ) * ( rA + rB ) );

  if ( contCond )
  {
    // Contact variables
    double overlap;
    Vector3 contVec; // Contact vector directed from b to a
    Point3 contPt;

    double axialOverlap = hA + hB - fabs( x[Z] );
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
    contPt = ( rA - .5 * radialOverlap ) * r; // assigning X & Y components
    contPt[Z] = x[Z] - sgn( x[Z] ) * ( hB - .5 * axialOverlap );
    // output
    ptCont.setContact( contPt );
    ptCont.setOverlapVector( overlap * contVec );
    ptCont.setOverlapDistance( - overlap );
  }
}




// ----------------------------------------------------------------------------
// Returns the contact point of two cylinders in the world of the first cylinder
// if the contact is Face-Band
inline void F2BContact(double rA, double hA, double rB, double hB,
                       Vector3 const& e, Point3 const& x, PointContact& ptCont)
{
  // Contact happens if
  bool contCond = false; // TRUE if contact occurs in this scenario
  if ( fabs( e[Z] ) < tol && fabs( x[Z] ) < hA + rB )
  {
    double S = fabs( x * e );
    if ( S < hB + rA )
    {
      double t = sqrt( normXY( x ) - S * S );
      double tStar = S < hB ? rA : sqrt( rA * rA - ( S - hB ) * ( S - hB ) );
      // Last condition - assuring Face A is in contact with Band B
      contCond = ( t < tStar ) &&
          ( !( S > rA && ( hB + rA - S ) < hA + rB - fabs( x[Z] ) ) );
    }
  }

  if ( contCond )
  {
    // Contact variables
    double overlap;
    Vector3 contVec; // Contact vector directed from b to a
    Point3 contPt;

    // amount of overlap
    overlap = hA + rB - fabs( x[Z] );
    // contact vector
    contVec = Vector3( 0., 0., - sgn( x[Z] ));
    // contact point
    Point3 ptA = x + rB * contVec + hB * e;
    Point3 ptB = x + rB * contVec - hB * e;
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
inline void F2EContact(double rA, double hA, double rB, double hB,
                       Vector3 const& e, Point3 const& x, PointContact& ptCont)
{
  // Contact happens if
  bool contCond = false; // TRUE if contact occurs in this scenario
  Point3 ptE;
  if ( fabs( x[Z] ) > hA )
  {
    // Vector3 r = ( e ^ ( e ^ zAxis ) ).normalized();
    Vector3 r = ( e ^ Vector3( e[Y], -e[X], 0. ) ).normalized();
    r = sgn( x[Z] ) * ( rB * r - hB * sgn( e[Z] ) * e );
    ptE = x + r;
    contCond = ( fabs( ptE[Z] ) < hA ) &&
               ( normXY( ptE ) < rA * rA ) &&
               ( hA - fabs( ptE[Z] ) < rA - sqrt( normXY( ptE ) ) );
  }

  if ( contCond )
  {
    // Contact variables
    double overlap;
    Vector3 contVec; // Contact vector directed from b to a
    Point3 contPt;

    // amount of overlap
    overlap = hA - fabs( ptE[Z] );
    // contact vector
    contVec = Vector3( 0., 0., - sgn( x[Z] ) );
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
inline void B2BSkewContact( double rA, double hA, double rB, double hB,
                            Vector3 const& e, Point3 const& x,
                            PointContact& ptCont )
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
    if ( fabs( lBStar ) < hB )
    {
      lAStar = x[Z] + lBStar * e[Z];
      contCond = ( fabs( lAStar ) < hA );
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
inline void B2EContact(double rA, double hA, double rB, double hB,
                       Vector3 const& e, Point3 const& x, PointContact& ptCont)
{
  // Contact happens if
  bool contCond = false; // TRUE if contact occurs in this scenario
  // Some variables
  Point3 const& c1 = x + hB * e;
  Point3 const& c2 = x - hB * e;
  Point3 const& ptCenter = normXY( c1 ) < normXY( c2 ) ? c1 : c2;
  // additional condition to avoid solving the polynomial
  if ( ( fabs( ptCenter[Z] ) > hA + rB ) ||
       ( sqrt( normXY( ptCenter ) ) > rA + rB ) )
    return;

  // Vector3 const& u = ( e ^ Vector3( 0., 0., 1. ) ).normalized();
  Vector3 const& u = ( Vector3( e[Y], -e[X], 0. ) ).normalized();
  Vector3 const& v = ( u ^ e ).normalized();

  // Misc variables
  // double const a = rB * ( normXY( v ) - normXY( u ) );
  double const a = rB * ( normXY( v ) - 1. );
  double const b = dotXY( u, ptCenter ) / a;
  double const c = dotXY( v, ptCenter ) / a;

  // double sint[4];
  // int nbRoots = 0;
  // solveQuartic( 2.*c, b*b + c*c - 1., -2.*c, -c*c, sint, nbRoots );
  //
  // Point3 ptA;
  // double cost;
  // for ( int i = 0; i < nbRoots; i++ )
  // {
  //   if ( fabs( sint[i] ) <= 1. )
  //   {
  //     cost = ( b * sint[i] ) / ( c + sint[i] );
  //     if ( fabs( cost ) > 1. )
  //       cost = sgn( cost ) * sqrt( 1. - sint[i]*sint[i] );
  //     ptA = ptCenter + rB * cost * u + rB * sint[i] * v;
  //     if ( normXY( ptA ) < rA * rA && fabs( ptA[Z] ) < hA )
  //     {
  //       contCond = ( rA - sqrt( normXY( ptA ) ) < hA - fabs( ptA[Z] ) );
  //       break;
  //     }
  //   }
  // }

  double sint[4];
  int nbRoots = 0;
  solveQuartic( 2.*c, b*b + c*c - 1., -2.*c, -c*c, sint, nbRoots );

  double sn = 1., cs = 0.;
  for ( int i = 0; i < nbRoots; i++ )
  {
    sn = sint[i];
    if ( sn * c >= 0. )
    {
      cs = sgn( b ) * sqrt( 1. - sn * sn );
      break;
    }
  }

  Point3 const& ptA = ptCenter + rB * cs * u + rB * sn * v;
  contCond = ( fabs( ptA[Z] ) < hA ) &&
             ( normXY( ptA ) < rA * rA ) &&
             ( rA - sqrt( normXY( ptA ) ) < hA - fabs( ptA[Z] ) );

  // // Gradient descent for finging the angle corresponding to minimal distance
  // double theta = PI/4.;
  // double delta = 1.;
  // double sn, cs;
  // int iter;
  // for ( iter = 0; iter < 100 && abs( delta ) > tol; iter++ )
  // {
  //   sn = sin(theta);
  //   cs = cos(theta);
  //   delta = 2. * ( ( c + sn ) * cs - b * sn );
  //   theta -= 0.05 * delta;
  // }
  // sn = sin(theta);
  // cs = cos(theta);
  // Point3 const& ptA = ptCenter + rB * cs * u + rB * sn * v;
  // contCond = ( fabs( ptA[Z] ) < hA ) &&
  //            ( normXY( ptA ) < rA * rA ) &&
  //            ( rA - sqrt( normXY( ptA ) ) < hA - fabs( ptA[Z] ) );

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
inline void E2EContact(double rA, double hA, double rB, double hB,
                       Vector3 const& e, Point3 const& x, PointContact& ptCont)
{
  // Contact happens if
  bool contCond = false; // TRUE if contact occurs in this scenario
  // Some variables
  Point3 c1 = x + hB * e;
  Point3 c2 = x - hB * e;
  double d1 = pow( sqrt( normXY(c1) ) - rA, 2 ) + pow( fabs( c1[Z] ) - hA, 2 );
  double d2 = pow( sqrt( normXY(c2) ) - rA, 2 ) + pow( fabs( c2[Z] ) - hA, 2 );
  Point3 const& ptCenter1 = d1 < d2 ? c1 : c2; // decide on edge of B
  // additional condition to avoid solving the polynomial
  if ( ( d1 < d2 ? d1 : d2 ) > rB * rB )
    return;

  // Vector3 const& u1 = ( e1 ^ Vector3( 0., 0., 1. ) ).normalized();
  Vector3 const& u1 = ( Vector3( e[Y], -e[X], 0. ) ).normalized();
  Vector3 const& v1 = ( u1 ^ e ).normalized();

  // Misc variables
  // double const p1 = rB * rB * ( normXY( v1 ) - normXY( u1 ) ) / 2.;
  double const p1 = rB * rB * ( normXY( v1 ) - 1. ) / 2.;
  double const q1 = rB * dotXY( u1, ptCenter1 ) / p1;
  double const r1 = rB * dotXY( v1, ptCenter1 ) / p1;
  // double const s1 = ( rA*rA - normXY( ptCenter1 ) - rB*rB*normXY(u1) ) / p1;
  double const s1 = ( rA * rA - normXY( ptCenter1 ) - rB * rB ) / p1;

  double sol1[4];
  int nbRoots;
  solveQuartic( 2.*r1, q1*q1 + r1*r1 - s1, -r1*s1, s1*s1/4. - q1*q1,
                sol1, nbRoots );

  double sint1 = 1., cost1 = 0.;
  for ( int i = 0; i < nbRoots; i++ )
  {
    sint1 = sol1[i];
    cost1 = ( s1 / 2. - ( r1 + sint1 ) * sint1 ) / q1;
    if ( fabs( ptCenter1[Z] + rB * cost1 * u1[Z] + rB * sint1 * v1[Z] ) < hA )
    {
      contCond = true;
      break;
    }
  }

  // Contact
  if ( contCond )
  {
    // Contact variables
    double overlap;
    Vector3 contVec; // Contact vector directed from b to a
    Point3 contPt;

    // The band point
    Point3 const& ptA = ptCenter1 + rB * cost1 * u1 + rB * sint1 * v1;

    // The face point
    sint1 = ( sgn( ptA[Z] ) * hA - ptCenter1[Z] ) / ( rB * v1[Z] );
    if ( fabs( sint1 ) > 1. )
      return;
    cost1 = sqrt( 1. - sint1 * sint1 );
    Point3 ptB = ptCenter1 + rB * cost1 * u1 + rB * sint1 * v1;
    if ( normXY( ptB ) > rA * rA )
      ptB = ptCenter1 + rB * (-cost1) * u1 + rB * sint1 * v1;
    // Theoretically, we don't need to check the following. Numerically,
    // it should be checked.
    if ( normXY( ptB ) > rA * rA )
      return;

    // Finding contacting edge points of cylinder A
    Matrix rotMatA2B;
    rotateVec2VecZ( e, rotMatA2B );
    Point3 const& x2 = rotMatA2B * ( - x );
    Vector3 const& e2 = ( rotMatA2B * Vector3( 0., 0., 1. ) ).normalized();
    c1 = x2 + hA * e2;
    c2 = x2 - hA * e2;
    d1 = pow( sqrt( normXY( c1 ) ) - rB, 2) + pow( fabs( c1[Z] ) - hB, 2);
    d2 = pow( sqrt( normXY( c2 ) ) - rB, 2) + pow( fabs( c2[Z] ) - hB, 2);
    Point3 const& ptCenter2 = d1 < d2 ? c1 : c2;

    // Vector3 const& u2 = ( e2 ^ Vector3( 0., 0., 1. ) ).normalized();
    Vector3 const& u2 = ( Vector3( e2[Y], -e2[X], 0. ) ).normalized();
    Vector3 const& v2 = ( u2 ^ e2 ).normalized();

    // Misc variables
    // double const p2 = rA * rA * ( normXY( v2 ) - normXY( u2 ) ) / 2.;
    double const p2 = rA * rA * ( normXY( v2 ) - 1. ) / 2.;
    double const q2 = rA * dotXY( u2, ptCenter2 ) / p2;
    double const r2 = rA * dotXY( v2, ptCenter2 ) / p2;
    // double const s2 = ( rB*rB - normXY(ptCenter2) - rA*rA*normXY(u2) ) / p2;
    double const s2 = ( rB * rB - normXY( ptCenter2 ) - rA * rA ) / p2;

    double sol2[4];
    solveQuartic( 2.*r2, q2*q2 + r2*r2 - s2, -r2*s2, s2*s2/4. - q2*q2,
                  sol2, nbRoots );

    double sint2 = 1., cost2 = 0.;
    for ( int i = 0; i < nbRoots; i++ )
    {
      sint2 = sol2[i];
      cost2 = ( s2 / 2. - ( r2 + sint2 ) * sint2 ) / q2;
      if ( fabs( ptCenter2[Z] + rA * cost2 * u2[Z] + rA * sint2 * v2[Z] ) < hB )
        break;
    }
    Point3 ptC = ptCenter2 + rA * cost2 * u2 + rA * sint2 * v2;
    // Theoretically, we don't need to check the following. Numerically,
    //  it should be checked.
    if ( fabs( ptC[Z] ) > hB )
      return;

    // Finding the face point
    sint2 = ( sgn( ptC[Z] ) * hB - ptCenter2[Z] ) / ( rA * v2[Z] );
    if ( fabs( sint2 ) > 1. )
      return;
    cost2 = sqrt( 1. - sint2 * sint2 );
    Point3 ptD = ptCenter2 + rA * cost2 * u2 + rA * sint2 * v2;
    if ( normXY( ptD ) > rB * rB )
      ptD = ptCenter2 + rA * (-cost2) * u2 + rA * sint2 * v2;
    // Theoretically, we don't need to check the following. Numerically,
    //  it should be checked.
    if ( normXY( ptD ) > rB * rB )
      return;

    // Contact points in the coordinate system of cylinder A
    ptC = x + Point3( transpose( rotMatA2B ) * ( ptC ) );
    ptD = x + Point3( transpose( rotMatA2B ) * ( ptD ) );

    // Theoretically, we don't need to check the following. Numerically,
    //  it should be checked.
    if ( fabs( ptC[Z] ) > hA || normXY( ptC ) > rA * rA ||
         fabs( ptD[Z] ) > hA || normXY( ptD ) > rA * rA )
      return;

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




// ----------------------------------------------------------------------------
// END OF LOW-LEVEL ROUTINES FOR CYLINDERS CONTACTS
// ----------------------------------------------------------------------------
// Returns the contact point of two cylinders
PointContact intersect( Cylinder const& a,
                        Cylinder const& b,
                        Transform const& a2w,
                        Transform const& b2w )
{
  // Variables
  double const rA = a.m_radius;
  double const rB = b.m_radius;
  double const hA = b.m_halfHeight;
  double const hB = b.m_halfHeight;
  Vector3 const& eA = a2w.getBasis() * Vector3( 0., 1., 0. );
  Vector3 const& eB = b2w.getBasis() * Vector3( 0., 1., 0. );

  PointContact ptCont = PointNoContact;

  // Relative positions - B w.r.t. A
  Matrix rotMatA;
  rotateVec2VecZ( eA, rotMatA );
  Point3 const& x_B2A = rotMatA * ( *(b2w.getOrigin()) - *(a2w.getOrigin()) );
  Vector3 const& e_B2A = ( rotMatA * eB ).normalized();
  
  // Relative positions - A w.r.t. B
  Matrix rotMatB;
  rotateVec2VecZ( eB, rotMatB );
  Point3 const& x_A2B = rotMatB * ( *(a2w.getOrigin()) - *(b2w.getOrigin()) );
  Vector3 const& e_A2B = ( rotMatB * eA ).normalized();
  

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