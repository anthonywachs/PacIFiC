#include "GrainsBuilderFactory.hh"
#include "RigidBody.hh"
#include "ConvexBuilderFactory.hh"
#include "Torsor.hh"
#include "Vector3.hh"
#include <math.h>
#include <stdlib.h>


// ----------------------------------------------------------------------------
// Default constructor
RigidBody::RigidBody() 
  : m_convex( NULL )
  , m_circumscribedRadius( 0.0 )
{}




// ----------------------------------------------------------------------------
// Constructor with a convex and a transformation as input parameters
RigidBody::RigidBody( Convex* convex_, Transform const& position_ ) 
  : m_transform( position_ )
  , m_convex( convex_ )
{
  m_circumscribedRadius = m_convex->computeCircumscribedRadius();
}




// ----------------------------------------------------------------------------
// Copy constructor
RigidBody::RigidBody( RigidBody const& form ) 
  : m_transform( form.m_transform )
  , m_convex( NULL )
  , m_circumscribedRadius( form.m_circumscribedRadius )
{
  if ( form.m_convex ) m_convex = form.m_convex->clone();
}




// ----------------------------------------------------------------------------
// Destructor
RigidBody::~RigidBody()
{
  delete m_convex;
}




// ----------------------------------------------------------------------------
// Returns the bounding box of the rigid body in its current configuration
BBox RigidBody::BoxRigidBody() const 
{
  return ( m_convex->bbox( m_transform ) );
}




// ----------------------------------------------------------------------------
// Computes the inertia tensor and the inverse of the inertia tensor
bool RigidBody::BuildInertia( double *inertia, double *inertia_1 ) const
{
  bool status = m_convex->BuildInertia( inertia, inertia_1 );

  Matrix m( inertia[0], inertia[1], inertia[2],
	inertia[1], inertia[3], inertia[4],
	inertia[2], inertia[4], inertia[5] );
  Matrix m_1( inertia_1[0], inertia_1[1], inertia_1[2],
	inertia_1[1], inertia_1[3], inertia_1[4],
	inertia_1[2], inertia_1[4], inertia_1[5] );

  Matrix base = m_transform.getBasis();
  m   = base * m   * base.transpose();
  m_1 = base * m_1 * base.transpose();

  inertia[0] = m[0][0];
  inertia[1] = m[1][0];
  inertia[2] = m[2][0];
  inertia[3] = m[1][1];
  inertia[4] = m[2][1];
  inertia[5] = m[2][2];

  inertia_1[0] = m_1[0][0];
  inertia_1[1] = m_1[1][0];
  inertia_1[2] = m_1[2][0];
  inertia_1[3] = m_1[1][1];
  inertia_1[4] = m_1[2][1];
  inertia_1[5] = m_1[2][2];

  return ( status );
}




// ----------------------------------------------------------------------------
// Returns the distance to another rigid body
double RigidBody::DistanceTo( RigidBody const& neighbor ) const 
{
  Point3 A, B;
  int nbIter;
  return ( closest_points( *m_convex, *neighbor.m_convex, m_transform, 
  	neighbor.m_transform, A, B, nbIter ) );
}




// ----------------------------------------------------------------------------
// Returns a pointer to the rigid body' center of mass position
Point3 const* RigidBody::getCentre() const
{
  return ( m_transform.getOrigin() );
}




// ----------------------------------------------------------------------------
// Copies the rigid body' center of mass position in a 3-element 1D array 
void RigidBody::getCentre( double *pos ) const
{
  Point3 const* ori = m_transform.getOrigin() ;  
  pos[X] =  (*ori)[X];
  pos[Y] =  (*ori)[Y];
  pos[Z] =  (*ori)[Z];
}




// ----------------------------------------------------------------------------
// Returns a pointer to the convex shape
Convex const* RigidBody::getConvex() const
{
  return ( m_convex );
}




// ----------------------------------------------------------------------------
// Returns a pointer to the rigid body's transformation
Transform const* RigidBody::getTransform() const
{
  return ( &m_transform );
}




// ----------------------------------------------------------------------------
// Returns a pointer to the rigid body's transformation
Transform* RigidBody::getTransform()
{
  return ( &m_transform );
}




// ----------------------------------------------------------------------------
// Returns the rigid body's circumscribed radius
double RigidBody::getCircumscribedRadius() const
{
  return ( m_circumscribedRadius );
}




// ----------------------------------------------------------------------------
// Returns the rigid body volume
double RigidBody::getVolume() const
{
  return ( m_convex->getVolume() );
}




// ----------------------------------------------------------------------------
// Returns whether 2 rigid bodies intersect 
bool intersect( RigidBody const& a, RigidBody const& b )
{
  Vector3 v = *(b.m_transform.getOrigin()) - *(a.m_transform.getOrigin());
  if ( Norm(v) < EPSILON ) return true;  
  return ( intersect( *a.m_convex, *b.m_convex, a.m_transform, b.m_transform, 
  	v ) );
}




// ----------------------------------------------------------------------------
// Returns whether there is geometric contact with another rigid body
bool RigidBody::isContact( RigidBody const& neighbor ) const 
{
  bool contact;

  Vector3 v = *(neighbor.m_transform.getOrigin()) - *(m_transform.getOrigin());
  
  if ( Norm(v) < EPSILON ) return true;
  
  contact = intersect( *m_convex, *(neighbor.m_convex), 
	m_transform, neighbor.m_transform, v );    

  return ( contact );
}




// ----------------------------------------------------------------------------
// Returns whether the rigid body is close to another rigid body in
// the sense of whether their respective bounding boxes overlap
bool RigidBody::isClose( RigidBody const& neighbor ) const
{
  bool contact;
  BBox box0 = (*this).BoxRigidBody();
  BBox box1 = neighbor.BoxRigidBody();
  contact = intersect( box0, box1 );

  return ( contact );
}




// ----------------------------------------------------------------------------
// Applies a rotation defined by a quaternion to the rigid body
void RigidBody::Rotate( Quaternion const& q ) 
{
  m_transform.composeLeftByRotation( q );
}




// ----------------------------------------------------------------------------
// Sets the origin of the rigid body's transformation
void RigidBody::setOrigin( double const* pos ) 
{
  m_transform.setOrigin( pos );
}




// ----------------------------------------------------------------------------
// Sets the origin of the rigid body's transformation
void RigidBody::setOrigin( double gx, double gy, double gz ) 
{
  m_transform.setOrigin( gx, gy, gz );
}




// ----------------------------------------------------------------------------
// Sets the origin of the rigid body's transformation
void RigidBody::setOrigin( Point3 const& pos ) 
{
  m_transform.setOrigin( pos[X], pos[Y], pos[Z] );
}




// ----------------------------------------------------------------------------
// Sets the rigid body's transformation with an 1D array of 16 
// values (see class Transform for details)
void RigidBody::setTransform( double const* pos ) 
{
  m_transform.setValue( pos );
}




// ----------------------------------------------------------------------------
// Sets the rigid body's circumscribed radius
void RigidBody::setCircumscribedRadius( double r ) 
{
  m_circumscribedRadius = r;
}




// ----------------------------------------------------------------------------
// Sets the rigid body's transformation with a transformation
void RigidBody::setTransform( Transform const& transform_ )
{
  m_transform = transform_;
}




// ----------------------------------------------------------------------------
// Applies a transformation trot to the right, i.e., this = this o trot 
void RigidBody::composeRightByTransform( Transform const& trot )
{
  m_transform.composeRightByTransform( trot );
}  




// ----------------------------------------------------------------------------
// Applies a transformation trot to the left, i.e., this = trot o 
// this, which means first this then trot
void RigidBody::composeLeftByTransform( Transform const& trot )
{
  m_transform.composeLeftByTransform( trot );
} 




// ----------------------------------------------------------------------------
// Applies a rotation defined by a transformation trot to the left, 
// i.e., this = trot o this, which means first this then trot. This composition
// leaves the origin unchanged but does not check that trot is indeed a 
// rotation
void RigidBody::composeLeftByRotation( Transform const& trot )
{
  m_transform.composeLeftByRotation( trot );
} 




// ----------------------------------------------------------------------------
// Applies a translation to the left, i.e., this = translation o this, which 
// means first this then translation. 
void RigidBody::composeLeftByTranslation( Vector3 const& v )
{
  m_transform.composeLeftByTranslation( v );
} 




// ----------------------------------------------------------------------------
// Reads the rigid body's transformation from an input stream
void RigidBody::readPosition( istream& fileIn )
{
  fileIn >> m_transform;
  m_circumscribedRadius = m_convex->computeCircumscribedRadius();
}




// ----------------------------------------------------------------------------
// Reads the rigid body's transformation from an input stream with
// the 2014 reload format 
void RigidBody::readPosition2014( istream& fileIn )
{
  m_transform.readTransform2014( fileIn );
  m_circumscribedRadius = m_convex->computeCircumscribedRadius();
}




// ----------------------------------------------------------------------------
// Reads the rigid body's transformation in binary format from an 
// input stream with the 2014 reload format
void RigidBody::readPosition2014_binary( istream& fileIn )
{
  m_transform.readTransform2014_binary( fileIn );
  m_circumscribedRadius = m_convex->computeCircumscribedRadius();
}




// ----------------------------------------------------------------------------
// Writes the rigid body's transformation in an output stream
void RigidBody::writePosition( ostream& fileOut ) const
{
  m_transform.writeTransform( fileOut ); 
}




// ----------------------------------------------------------------------------
// Writes the rigid body's "static" data, i.e., the convex geometric
// description only (without any transformation)
void RigidBody::writeStatic( ostream& fileOut ) const
{
  fileOut << *m_convex;
}




// ----------------------------------------------------------------------------
// Writes the geometric features of the rigid body in its current
// position in an output stream in a format suitable to the coupling with a 
// fluid solver. Note: this method works for discs, polygons, polyhedrons, 
// spheres and 3D cylinders
void RigidBody::writePositionInFluid( ostream& fluid )
{
  Point3 pointEnvelop;
  vector<Point3> allPoints = m_convex->getEnvelope();
  vector<Point3>::iterator point;
  
  fluid << " " << allPoints.size() << endl;

  // 2D case
  if ( GrainsBuilderFactory::getContext() == DIM_2 )
  {
    // Points describing the shape
    for (point=allPoints.begin(); point!=allPoints.end(); point++) 
    {
      pointEnvelop = m_transform(*point);
      fluid << GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
		pointEnvelop[X] ) << " " << 
	GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
		pointEnvelop[Y] ) << endl;     
    }
  }
  
  // 3D case 
  else if ( GrainsBuilderFactory::getContext() == DIM_3 )
  {
    // Points describing the shape
    for (point=allPoints.begin(); point!=allPoints.end(); point++) 
    {
      pointEnvelop = m_transform(*point);
      fluid << GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
		pointEnvelop[X] ) << " " <<  
	GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
		pointEnvelop[Y] ) << " " <<  
	GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
		pointEnvelop[Z] ) << endl; 
    }

    // Faces describing the shape
    vector< vector<int> > const* allFaces  = m_convex->getFaces();
    vector< vector<int> >::const_iterator face;
    if ( allFaces )
    {
      fluid << allFaces->size() << endl;
      for (face=allFaces->begin(); face!=allFaces->end(); face++) 
      {
        vector<int>::const_iterator index;
        fluid << (*face).size() << " ";
        for (index=(*face).begin(); index!=(*face).end(); index++)
          fluid << (*index) << " ";
        fluid << endl;
      }
    }
    else fluid << "0" << endl;  
  }
  
  // Problem with dimension
  else
  {
    cout << "!!! Warning: Physical dimension undefined (DIM_2 or DIM_3)" 
    	<< endl;
    exit(0);
  }
}




// ----------------------------------------------------------------------------
// Copies the rigid body transformation in a 1D array
void RigidBody::copyTransform( double* vit, int i ) const
{
  m_transform.copyTransform( vit, i );
}




// ----------------------------------------------------------------------------
// Copies the rigid body transformation in a 1D array composed on 
// the left by a translation (useful for periodic particles in parallel)
void RigidBody::copyTransform( double* vit, int i, Vector3 const& vec ) const
{
  m_transform.copyTransform( vit, i, vec );
}




// ----------------------------------------------------------------------------
// Writes a list of points describing the rigid body's convex shape 
// in a Paraview format 
void RigidBody::write_polygonsPts_PARAVIEW( ostream& f, 
	Vector3 const* translation ) const
{
  m_convex->write_polygonsPts_PARAVIEW( f, m_transform, translation );
}




// ----------------------------------------------------------------------------
// Writes the rigid body's convex shape in a STL format
void RigidBody::write_convex_STL( ostream& f ) const
{
  m_convex->write_convex_STL( f, m_transform );
}




// ----------------------------------------------------------------------------
// Returns a list of points describing the rigid body's convex shape
// in a Paraview format
list<Point3> RigidBody::get_polygonsPts_PARAVIEW( Vector3 const* translation ) 
	const
{
  return ( m_convex->get_polygonsPts_PARAVIEW( m_transform, translation ) );
}




// ----------------------------------------------------------------------------
// Returns whether a point lies inside the rigid body
bool RigidBody::isIn( Point3 const& pt ) const
{
  Transform invT; 
  invT.setToInverseTransform( m_transform );
  return ( m_convex->isIn( invT( pt ) ) );
}
