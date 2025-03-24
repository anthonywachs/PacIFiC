#include "GrainsMPIWrapper.hh"
#include "RigidBodyWithCrust.hh"
#include "Box.hh"
#include "Cylinder.hh"
#include "ConvexBuilderFactory.hh"
#include "PointContact.hh"
#include "PointC.hh"
#include "GrainsExec.hh"
#include "Particle.hh"
#include "BBox.hh"
#include "GJK_SV.hh"

#include <fstream>
#include <sstream>
#include <iostream>

// ----------------------------------------------------------------------------
// Default constructor
RigidBodyWithCrust::RigidBodyWithCrust()
  : RigidBody()
  , m_crustThickness( 0.0 )
  , m_scaling( NULL )
  , m_transformWithCrust( NULL )
  , m_transformWithCrust_computed( false )
{}




// ----------------------------------------------------------------------------
// Copy constructor
RigidBodyWithCrust::RigidBodyWithCrust( RigidBodyWithCrust const& rbwc )
  : RigidBody( rbwc )
  , m_crustThickness( rbwc.m_crustThickness )
  , m_scaling( NULL )
  , m_transformWithCrust( NULL )
  , m_transformWithCrust_computed( false )
{
  if ( rbwc.m_scaling ) m_scaling = new Vector3( *rbwc.m_scaling ) ;
  if ( rbwc.m_transformWithCrust )
    m_transformWithCrust = new Transform( *rbwc.m_transformWithCrust );
}




// ----------------------------------------------------------------------------
// Constructor with input parameters, used by composite component 
// whose own shape is not defined ( composite = true ) or to create component 
// from scratch in the code ( composite = false )
RigidBodyWithCrust::RigidBodyWithCrust( Convex* convex_,
	Transform const& position_,
	bool composite, double const& crust_thickness )
  : RigidBody( convex_, position_ )
  , m_crustThickness( crust_thickness )
  , m_scaling( NULL )
  , m_transformWithCrust( NULL )
  , m_transformWithCrust_computed( false )
{
  if ( !composite )
  {
    m_circumscribedRadius = m_convex->computeCircumscribedRadius();
    m_scaling = new Vector3;
    BBox box = m_convex->bbox( TransformIdentity );
    Vector3 const& extent = box.getExtent();

    // Scaling factor from bounding box
    (*m_scaling)[X] = extent[X] < EPSILON ?
  	1. : ( extent[X] - m_crustThickness ) / extent[X];
    (*m_scaling)[Y] = extent[Y] < EPSILON ?
  	1. : ( extent[Y] - m_crustThickness ) / extent[Y];
    (*m_scaling)[Z] = extent[Z] < EPSILON ?
  	1. : ( extent[Z] - m_crustThickness ) / extent[Z];

    // Transformation with crust
    m_transformWithCrust = new Transform( position_ );
    m_transformWithCrust_computed = false ;  
  }
  GrainsExec::setMinCrustThickness( m_crustThickness );
}




// ----------------------------------------------------------------------------
// Constructor from an input stream and a convex type
// The type is used in the case of elementary particles only, for standard
// particles type is an empty string
RigidBodyWithCrust::RigidBodyWithCrust( istream& fileIn, string type )
  : RigidBody()
{
  string cle;
  if ( !type.empty() ) cle = type;
  else fileIn >> cle;

  // Read the rigid body shape
  m_convex = ConvexBuilderFactory::create( cle, fileIn );
  m_boundingVolume = m_convex->
	computeBVolume( GrainsExec::m_colDetBoundingVolume );
  fileIn >> cle;
  assert( cle == "*END" );

  // Read the crust thickness
  fileIn >> cle >> m_crustThickness;
  GrainsExec::setMinCrustThickness( m_crustThickness );
  
  // Circumscribed radius and bounding box
  m_circumscribedRadius = m_convex->computeCircumscribedRadius();
  m_volume = m_convex->getVolume();
  m_scaling = new Vector3;
  BBox box = m_convex->bbox( TransformIdentity );
  Vector3 const& extent = box.getExtent();

  // Scaling factor from bounding box
  (*m_scaling)[X] = extent[X] < EPSILON ?
  	1. : ( extent[X] - m_crustThickness ) / extent[X];
  (*m_scaling)[Y] = extent[Y] < EPSILON ?
  	1. : ( extent[Y] - m_crustThickness ) / extent[Y];
  (*m_scaling)[Z] = extent[Z] < EPSILON ?
  	1. : ( extent[Z] - m_crustThickness ) / extent[Z];

  // Transformation with crust
  m_transformWithCrust = new Transform();
  m_transformWithCrust_computed = false ;
}




// ----------------------------------------------------------------------------
// Constructor with an XML node as an input parameter
RigidBodyWithCrust::RigidBodyWithCrust( DOMNode* root )
{
  // Convex
  DOMNode* forme = ReaderXML::getNode( root, "Convex" );
  m_convex = ConvexBuilderFactory::create( forme );
  m_crustThickness = ReaderXML::getNodeAttr_Double( forme, "CrustThickness" );
  m_boundingVolume = m_convex->
	computeBVolume( GrainsExec::m_colDetBoundingVolume );
  GrainsExec::setMinCrustThickness( m_crustThickness );

  // Transformation
  m_transform.load( root );

  // Circumscribed radius, bounding box
  m_circumscribedRadius = m_convex->computeCircumscribedRadius();
  m_volume = m_convex->getVolume();
  m_scaling = new Vector3;
  BBox box = m_convex->bbox( TransformIdentity );
  const Vector3& extent = box.getExtent();

  // Scaling factor from bounding box
  (*m_scaling)[X] = extent[X] < EPSILON ?
        1. : ( extent[X] - m_crustThickness ) / extent[X];
  (*m_scaling)[Y] = extent[Y] < EPSILON ?
        1. : ( extent[Y] - m_crustThickness ) / extent[Y];
  (*m_scaling)[Z] = extent[Z] < EPSILON ?
        1. : ( extent[Z] - m_crustThickness ) / extent[Z];

  // Transformation with crust
  m_transformWithCrust = new Transform();
  m_transformWithCrust_computed = false ;
}




// ----------------------------------------------------------------------------
// Destructor
RigidBodyWithCrust::~RigidBodyWithCrust()
{
  if ( m_scaling ) delete m_scaling;
  if ( m_transformWithCrust ) delete m_transformWithCrust;
}




// ----------------------------------------------------------------------------
// Returns the bounding box extended by the crust thickness
BBox RigidBodyWithCrust::BoxRigidBody() const
{
  BBox box = m_convex->bbox( m_transform );

  Vector3 extent = box.getExtent(), vec_addVdW( m_crustThickness );
  extent += vec_addVdW;
  box.setExtent( extent );

  return ( box );
}




// ----------------------------------------------------------------------------
// Returns the features of the contact: contact point location,
// overlap vector (vector joining the points on each rigid body surface that
// realize the minimal distance between the shrunk rigid bodies, divided by the
// minimal distance between the shrunk rigid bodies and multiplied by the
// sum of the crust thicknesses minus the minimal distance between the shrunk
// rigid bodies, i.e., minus the overlap), overlap distance = minimal distance
// between the shrunk rigid bodies minus the sum of the crust thicknesses.
// Note: contact exists if overlap distance is negative, i.e., minimal distance
// between the shrunk rigid bodies < sum of the crust thicknesses
PointContact RigidBodyWithCrust::ClosestPoint( RigidBodyWithCrust &neighbor )
{
  try 
  {
    Convex const* convexA = m_convex;
    Convex const* convexB = neighbor.m_convex;

    // Comment on the direction of the overlap vector
    // Assuming A and B are the centers of the 2 convex bodies
    // overlap_vector = overlap * Vector3(A to B)
    // If contact, overlap is negative and overlap_vector is from B to A
    // If no contact, overlap is positive and we do not care about the direction
    // of overlap_vector

    switch( convexA->getConvexType() )
    {
      case SPHERE:
        switch( convexB->getConvexType() )
        {
          case SPHERE:
	    return ( ClosestPointSPHERESPHERE( *this, neighbor ) );
	  
	  case BOX:
	    return ( ClosestPointSPHEREBOX( *this, neighbor ) );
          
	  default:
            break;
        }
        break;

      case DISC2D:
        switch( convexB->getConvexType() )	
        {
          case DISC2D:
            return ( ClosestPointSPHERESPHERE( *this, neighbor ) );
          
	  case BOX:
            return ( ClosestPointSPHEREBOX( *this, neighbor ) );
          
	  default:
            break;
        }
        break;

      case BOX:
        switch( convexB->getConvexType() )	
        {
          case SPHERE:
            return ( ClosestPointSPHEREBOX( *this, neighbor ) );
          
	  case DISC2D:
            return ( ClosestPointSPHEREBOX( *this, neighbor ) );
          
	  default:
            break;
        }
        break;
	
      // case CYLINDER:
      //   switch( convexB->getConvexType() )	
      //   {
      //     case CYLINDER:
      //       return( ClosestPointCYLINDERS( *this, neighbor ) );
      //     default:
      //       break;
      //   }
      //   break;

      default:
        break;
    }

    // Bounding Volume Check
    ++GrainsExec::m_nb_GJK_narrow_collision_detections;
    Vector3 gcagcb = *m_transform.getOrigin() -
    	*neighbor.m_transform.getOrigin();
    bool BVCheck = ( Norm(gcagcb) < m_circumscribedRadius + 
	neighbor.m_circumscribedRadius );
    if ( GrainsExec::m_colDetBoundingVolume && BVCheck )
      BVCheck = isContactBVolume( *this, neighbor );

    // General case for any pair of convex rigid bodies
    if ( BVCheck )
    {
      // In case one rigid body is a rectangle
      if ( convexA->getConvexType() == RECTANGLE2D ||
           convexB->getConvexType() == RECTANGLE2D )
      {
        --GrainsExec::m_nb_GJK_narrow_collision_detections;
	return ( ClosestPointRECTANGLE( *this, neighbor, false ) );
      }
      ++GrainsExec::m_nb_GJK_calls;

      // Distance between the 2 rigid bodies shrunk by their crust thickness
      Transform const* a2w = this->getTransformWithCrust();
      Transform const* b2w = neighbor.getTransformWithCrust();
      Point3 pointA, pointB;
      int nbIterGJK = 0;
      double distance = 0;
      // Choose the appropriate GJK version according to the input XML
      if ( GrainsExec::m_colDetGJK_SV ) // Signed-Volume
        distance = closest_points_GJK_SV( *m_convex, 
		*(neighbor.m_convex), *a2w, *b2w, pointA, pointB, nbIterGJK );
      else // default: Johnson
        distance = closest_points( *m_convex, 
		*(neighbor.m_convex), *a2w, *b2w, pointA, pointB, nbIterGJK );

      if ( distance < EPSILON )
      {
        cout << "ERR RigidBodyWithCrust::ClosestPoint on Processor "
             << (GrainsExec::m_MPI ? GrainsExec::getComm()->get_rank() : 0 ) 
             << endl;
        throw ContactError();
      }

      // Points A and B are in their respective local coordinate systems
      // Thus we transform them into the world coordinate system
      pointA = (*a2w)( pointA );
      pointB = (*b2w)( pointB );

      // Comment on the ba vector
      // pointA is the point realizing the shortest distance in rigid body A
      // pointB is the point realizing the shortest distance in rigid body B
      // thus pointA - pointB = ba is directed from B to A
      Vector3 ba = pointA - pointB;

      // Contact point definition as the mid point between pointA and pointB
      Point3 contact = pointA / 2.0 + pointB / 2.0;

      // Computation of the actual overlap vector
      // If contact, crustA + crustB - distance > 0, the overlap vector is
      // directed from B to A
      // If no contact, crustA + crustB - distance < 0 and we do not care about
      // the direction of the overlap vector
      Vector3 overlap_vector = ba / distance;
      overlap_vector.round();
      overlap_vector *= m_crustThickness + neighbor.m_crustThickness - distance;

      // Computation of the actual overlap distance = distance - crustA - crustB
      // If actual overlap distance < 0 => contact
      // otherwise no contact
      distance -= m_crustThickness + neighbor.m_crustThickness;

      return ( PointContact( contact, overlap_vector, distance, nbIterGJK ) );
    }
    else
      return ( PointNoContact );
  }
  catch ( ContactError const& ) { throw; }
}




// ----------------------------------------------------------------------------
// Returns the same as above, but we feed the initial direction for the GJK alg
// using the previous contact between the same two rigid bodies
PointContact RigidBodyWithCrust::ClosestPoint( RigidBodyWithCrust &neighbor,
	Vector3& initialDirection )
{
  try 
  {
    Convex const* convexA = m_convex;
    Convex const* convexB = neighbor.m_convex;

    // Comment on the direction of the overlap vector
    // Assuming A and B are the centers of the 2 convex bodies
    // overlap_vector = overlap * Vector3(A to B)
    // If contact, overlap is negative and overlap_vector is from B to A
    // If no contact, overlap is positive and we do not care about the direction
    // of overlap_vector

    switch( convexA->getConvexType() )
    {
      case SPHERE:
        switch( convexB->getConvexType() )
        {
          case SPHERE:
	    return ( ClosestPointSPHERESPHERE( *this, neighbor ) );

	  case BOX:
            return ( ClosestPointSPHEREBOX( *this, neighbor ) );
          
	  default:
            break;
        }
        break;

      case DISC2D:
        switch( convexB->getConvexType() )	
        {
          case DISC2D:
            return ( ClosestPointSPHERESPHERE( *this, neighbor ) );

          case BOX:
            return ( ClosestPointSPHEREBOX( *this, neighbor ) );

          default:
            break;
        }
        break;

      case BOX:
        switch( convexB->getConvexType() )	
        {
          case SPHERE:
            return ( ClosestPointSPHEREBOX( *this, neighbor ) );

          case DISC2D:
            return ( ClosestPointSPHEREBOX( *this, neighbor ) );

          default:
            break;
        }
        break;
	
      // case CYLINDER:
      //   switch( convexB->getConvexType() )	
      //   {
      //     case CYLINDER:
      //       return( ClosestPointCYLINDERS( *this, neighbor ) );
      //     default:
      //       break;
      //   }
      //   break;

      default:
        break;
    }

    // Bounding Volume Check
    ++GrainsExec::m_nb_GJK_narrow_collision_detections;
    Vector3 gcagcb = *m_transform.getOrigin() -
	*neighbor.m_transform.getOrigin();
    bool BVCheck = ( Norm(gcagcb) < m_circumscribedRadius + 
	neighbor.m_circumscribedRadius );
    if ( GrainsExec::m_colDetBoundingVolume && BVCheck )
      BVCheck = isContactBVolume( *this, neighbor );
            
    if ( BVCheck )
    {
      // In case one rigid body is a rectangle
      if ( convexA->getConvexType() == RECTANGLE2D ||
           convexB->getConvexType() == RECTANGLE2D )
      {
        --GrainsExec::m_nb_GJK_narrow_collision_detections;
	return ( ClosestPointRECTANGLE( *this, neighbor, false ) );
      }
      ++GrainsExec::m_nb_GJK_calls;
      
      // Distance between the 2 rigid bodies shrunk by their crust thickness
      Transform const* a2w = this->getTransformWithCrust();
      Transform const* b2w = neighbor.getTransformWithCrust();
      Point3 pointA, pointB;
      int nbIterGJK = 0;
      double distance = 0.;

      // Choose the appropriate GJK version according to the input XML
      if ( GrainsExec::m_colDetGJK_SV ) // Signed-Volume
        distance = closest_points_GJK_SV( *m_convex, 
		*(neighbor.m_convex), *a2w, *b2w, initialDirection, 
		pointA, pointB, nbIterGJK );
      else // default: Johnson
        distance = closest_points( *m_convex, 
		*(neighbor.m_convex), *a2w, *b2w, initialDirection, 
		pointA, pointB, nbIterGJK );

      if ( distance < EPSILON )
      {
        cout << "ERR RigidBodyWithCrust::ClosestPoint on Processor "
     	       << (GrainsExec::m_MPI ? GrainsExec::getComm()->get_rank() : 0 ) 
	           << endl;
        throw ContactError();
      }

      // Points A and B are in their respective local coordinate systems
      // Thus we transform them into the world coordinate system
      pointA = (*a2w)( pointA );
      pointB = (*b2w)( pointB );

      // Comment on the ba vector
      // pointA is the point realizing the shortest distance in rigid body A
      // pointB is the point realizing the shortest distance in rigid body B
      // thus pointA - pointB = ba is directed from B to A
      Vector3 ba = pointA - pointB;

      // Contact point definition as the mid point between pointA and pointB
      Point3 contact = pointA / 2.0 + pointB / 2.0;

      // Computation of the actual overlap vector
      // If contact, crustA + crustB - distance > 0, the overlap vector is
      // directed from B to A
      // If no contact, crustA + crustB - distance < 0 and we do not care 
      // about the direction of the overlap vector
      Vector3 overlap_vector = ba / distance;
      overlap_vector.round();
      overlap_vector *= m_crustThickness + neighbor.m_crustThickness - distance;

      // Computation of the actual overlap distance = distance - crustA 
      // - crustB
      // If actual overlap distance < 0 => contact
      // otherwise no contact
      distance -= m_crustThickness + neighbor.m_crustThickness;

      return ( PointContact( contact, overlap_vector, distance, nbIterGJK ) );
    }
    else
    {
      initialDirection = Vector3Null;
      return ( PointNoContact );
    }
  }
  catch ( ContactError const& ) { throw; }
}




// ----------------------------------------------------------------------------
// Returns the features of the contact when the overlap computed by
// ClosestPoint is too large, the method artificially increases the size of the
// crust thickness for this particular contact detection by a factor > 1 and
// imposes the overlap distance to minus the sum of the crust thicknesses. This
// is useful when a few contacts involve a slightly large overlap due for
// instance (i) to the user's contact parameters not being properly set, (ii)
// an unexpectedly high colliding velocity between 2 rigid bodies that
// constitutes a rare by physically meaningful event or (iii) when a contact
// has not been detected by GJK for a few time steps preceding the call to
// this method that constitutes an extremely rare event.
PointContact RigidBodyWithCrust::ClosestPoint_ErreurHandling(
	RigidBodyWithCrust const& neighbor, double const& factor, int const& id,
	int const& id_neighbor )
{
  try 
  {
    // General case for any pair of convex rigid bodies
    Point3 pointA, pointB;
    int nbIterGJK = 0;
    Transform a2w = this->getTransformWithCrust( factor, 0.5 );
    Transform b2w = neighbor.getTransformWithCrust( factor, 0.5 );
    double distance = closest_points( *m_convex, *(neighbor.m_convex), a2w, b2w,
	pointA, pointB, nbIterGJK );

    if ( distance < EPSILON )
    {
      cout << "ERR RigidBodyWithCrust::ClosestPoint_ErreurHandling on "
      	<< " Processor " << (GrainsExec::m_MPI ?
	GrainsExec::getComm()->get_rank() : 0 )
	<< " between components " << id << " and " << id_neighbor
	<< endl;
      throw ContactError();
    }
    else
    {
      cout << "Handling contact error on Processor "
	   << (GrainsExec::m_MPI ? GrainsExec::getComm()->get_rank() : 0 )
	   << " between components " << id << " and " << id_neighbor
	   << endl;
    }

    // Points A and B are in their respective local coordinate systems
    // Thus we transform them into the world coordinate system
    pointA = a2w( pointA );
    pointB = b2w( pointB );

    // Comment on the ba vector
    // pointA is the point realizing the shortest distance in rigid body A
    // pointB is the point realizing the shortest distance in rigid body B
    // thus pointA - pointB = ba is directed from B to A
    Vector3 ba = pointA - pointB;
    ba.normalize();
 
    // Contact point definition as the mid point between pointA and pointB
    Point3 contact = pointA / 2.0 + pointB / 2.0;

    // Imposed overlap distance equal to the sum of the crust thicknesses
    double imposed_overlap_distance = 1. * ( m_crustThickness
  	+ neighbor.m_crustThickness );

    return ( PointContact( contact, imposed_overlap_distance * ba,
  	- imposed_overlap_distance, nbIterGJK ) );
  }
  catch ( ContactError const& ) { throw; }
}




// ----------------------------------------------------------------------------
// Returns the features of the contact when the 2 rigid bodies are
// spheres, i.e., a SPHERE-SPHERE contact
PointContact ClosestPointSPHERESPHERE( RigidBodyWithCrust const& rbA,
	RigidBodyWithCrust const& rbB )
{
  try 
  {
    // Comment on the direction of the overlap vector
    // Assuming A and B are the centers of the 2 convex bodies
    // overlap_vector = overlap * Vector3(A to B)
    // If contact, overlap is negative and overlap_vector is from B to A
    // If no contact, overlap is positive and we do not care about the direction
    // of overlap_vector

    Point3 const* pointA  = rbA.getTransform()->getOrigin();
    Point3 const* pointB  = rbB.getTransform()->getOrigin();

    Vector3 vecteurAB = *pointB - *pointA;
    double  rayonA    = rbA.getCircumscribedRadius();
    double  rayonB    = rbB.getCircumscribedRadius();

    double  distance  = Norm( vecteurAB ) - ( rayonA + rayonB );
    if ( distance > 0. )
      return ( PointNoContact );
    else
    {
      double rdwA = rbA.getCrustThickness();
      double rdwB = rbB.getCrustThickness();
      if ( - distance >= rdwA + rdwB )
      {
        cout << "ERR RigidBodyWithCrust::ClosestPointSPHERE on Processor "
		<< (GrainsExec::m_MPI ? GrainsExec::getComm()->get_rank() : 0 )
		<< ": " << - distance << " & " << rdwA + rdwB << "\n";
        cout << "Position: " << *pointA << " " << *pointB << endl;
        throw ContactError();
      }

      Point3 contact  = *pointA + ( rayonA + 0.5 * distance ) *
    	vecteurAB / Norm( vecteurAB );
      Vector3 overlap_vector = distance * ( vecteurAB / Norm( vecteurAB ) );

      return ( PointContact( contact, overlap_vector, distance, 1 ) );
    }
  }
  catch ( ContactError const& ) { throw; }
}




// ----------------------------------------------------------------------------
// Returns the features of the contact when the 1 rigid body is a sphere
// and the other rigid body is a box, i.e., a SPHERE-BOX contact
PointContact ClosestPointSPHEREBOX( RigidBodyWithCrust const& rbA,
	RigidBodyWithCrust const& rbB )
{
  try 
  {
    // Comment on the direction of the overlap vector
    // Assuming A and B are the centers of the 2 convex bodies
    // overlap_vector = overlap * Vector3(A to B)
    // If contact, overlap is negative and overlap_vector is from B to A
    // If no contact, overlap is positive and we do not care about the direction
    // of overlap_vector

    Convex const* convexA = rbA.getConvex();
    Convex const* convexB = rbB.getConvex();
    double rdwA = rbA.getCrustThickness();
    double rdwB = rbB.getCrustThickness();
    double overlap=0.;
    Point3 contactPoint, contact;

    if ( convexA->getConvexType() == SPHERE
  	|| convexA->getConvexType() == DISC2D )
    {
      Box const* convexBoxB = (Box const*)(convexB);
      Point3 const* pointA  = rbA.getTransform()->getOrigin();
      double rayonA = rbA.getCircumscribedRadius();
      Transform const* transfB = rbB.getTransform();
      Transform w2b;
      w2b.setToInverseTransform( *transfB );
      contactPoint = convexBoxB->IntersectionPointSPHERE( w2b(*pointA), rayonA,
    	overlap );	
      if ( overlap < 0. )
      {
        if ( - overlap >=  rdwA + rdwB )
        {
          cout << "ERR RigidBodyWithCrust::ClosestPointSPHEREBOX on Processor "
      	  	<< ( GrainsExec::m_MPI ? GrainsExec::getComm()->get_rank() : 0 )
		<< ": " << - overlap << " & " << rdwA + rdwB << endl;
	  throw ContactError();
        }
        contact = (*transfB)( contactPoint );
        Vector3 AB = contact - *pointA;
        Vector3 overlap_vector( ( overlap / Norm(AB) ) * AB );
        return ( PointContact( contact, overlap_vector, overlap, 1 ) );
      }
      else return ( PointNoContact );
    }
    else
    {
      Box const* convexBoxA = (Box const*)(convexA);
      Point3 const* pointB  = rbB.getTransform()->getOrigin();
      double rayonB = rbB.getCircumscribedRadius();
      Transform const* transfA = rbA.getTransform();
      Transform w2a;
      w2a.setToInverseTransform( *transfA );
      contactPoint = convexBoxA->IntersectionPointSPHERE( w2a(*pointB), rayonB,
    	overlap );
      if ( overlap < 0. )
      {
        if ( - overlap >= rdwA + rdwB )
        {
          cout << "ERR RigidBodyWithCrust::ClosestPointSPHEREBOX on Processor "
      		<< ( GrainsExec::m_MPI ? GrainsExec::getComm()->get_rank() : 0 )
	     	 << ": " << - overlap << " & " << rdwA + rdwB << endl;
	  throw ContactError();
        }
        contact = (*transfA)( contactPoint );
        Vector3 AB = *pointB - contact;
        Vector3 overlap_vector( ( overlap / Norm(AB) ) * AB );
        return ( PointContact( contact, overlap_vector, overlap, 1 ) );
      }
      else return ( PointNoContact );
    }
  }
  catch ( ContactError const& ) { throw; }
}



// ----------------------------------------------------------------------------
// Returns the crust thickness
double RigidBodyWithCrust::getCrustThickness() const
{
  return ( m_crustThickness );
}




// ----------------------------------------------------------------------------
// Sets the crust thickness
void RigidBodyWithCrust::setCrustThickness( double cthickness_ )
{
  m_crustThickness = cthickness_;
  GrainsExec::setMinCrustThickness( m_crustThickness );  
}




// ----------------------------------------------------------------------------
// Returns a pointer to the rigid body's transformation with the
// scaling by the crust thickness to shrink the rigid body
Transform const* RigidBodyWithCrust::getTransformWithCrust()
{
  if ( !m_transformWithCrust_computed )
  {
    *m_transformWithCrust = m_transform;
    m_transformWithCrust->composeWithScaling( (*m_scaling)[X], (*m_scaling)[Y],
  	(*m_scaling)[Z] );
    m_transformWithCrust_computed = true;
  }

  return ( m_transformWithCrust );
}




// ----------------------------------------------------------------------------
// Returns the rigid body's transformation with the
// scaling by the crust thickness multiplied by a factor > 1 to shrink the
// rigid body. The shrinkage is capped by a minimum scaling usually set to a
// value between 0.5 and 1 when this method is called
Transform RigidBodyWithCrust::getTransformWithCrust( double const& factor,
	double const& min_scaling ) const
{
  BBox box = m_convex->bbox( TransformIdentity );
  Vector3 const& extent = box.getExtent();

  double scaleX, scaleY, scaleZ;
  scaleX = extent[X] < EPSILON ? 1. :
  	( extent[X] - factor * m_crustThickness ) / extent[X];
  scaleY = ( extent[Y] < EPSILON ) ? 1. :
  	( extent[Y] - factor * m_crustThickness ) / extent[Y];
  scaleZ = ( extent[Z]<EPSILON ) ? 1. :
  	( extent[Z] - factor * m_crustThickness ) / extent[Z];

  if ( min_scaling != 1. )
  {
    scaleX = max( scaleX, min_scaling );
    scaleY = max( scaleY, min_scaling );
    scaleZ = max( scaleZ, min_scaling );
  }

  Transform vdw = m_transform;
  vdw.composeWithScaling( scaleX, scaleY, scaleZ );

  return ( vdw );
}




// ----------------------------------------------------------------------------
// Returns whether the rigid body is close to another rigid body in
// the sense of whether their respective bounding boxes minus their crust
// thickness overlap
bool RigidBodyWithCrust::isClose( RigidBodyWithCrust const& neighbor ) const
{
  bool contact;

  BBox boxA = (*this).BoxRigidBody();
  BBox boxB = neighbor.BoxRigidBody();
  Vector3 const& extentA = boxA.getExtent();
  Vector3 const& extentB = boxB.getExtent();
  Point3 const& pointA  = boxA.getCenter();
  Point3 const& pointB  = boxB.getCenter();
  Vector3 ab = pointA - pointB;

  double x, y, z;
  double ctA = (*this).m_crustThickness;
  double ctB = neighbor.m_crustThickness;
  x = fabs( ab[X] ) - ( extentA[X]
  	+ extentB[X] * ( extentA[X] + ctA ) / extentA[X]
  	* ( extentB[X] + ctB ) / extentB[X] );
  y = fabs( ab[Y] ) - ( extentA[Y]
  	+ extentB[Y] * ( extentA[Y] + ctA ) / extentA[Y]
  	* ( extentB[Y] + ctB ) / extentB[Y] );
  z = fabs( ab[Z] ) - ( extentA[Z]
  	+ extentB[Z] * ( extentA[Z] + ctA ) / extentA[Z]
  	* ( extentB[Z] + ctB ) / extentB[Z] );

  contact = true;
  if ( x>0. || y>0. || z>0. ) contact = false;

  return ( contact );
}




// ----------------------------------------------------------------------------
// Returns whether there is geometric contact with another rigid
// body in the sense of ClosestPoint, i.e., if minimal distance
// between the shrunk rigid bodies < sum of the crust thicknesses
bool RigidBodyWithCrust::isContact( RigidBodyWithCrust& neighbor )
{
  bool contact = false;
  PointContact pc;
  Convex const* convexA = m_convex;
  Convex const* convexB = neighbor.m_convex;

  switch( convexA->getConvexType() )
  {
    case SPHERE:
      switch( convexB->getConvexType() )	
      {
	case SPHERE:
	  return ( isContactSPHERESPHERE( *this,  neighbor ) );
	  break;
	    
	case BOX:
	  return ( isContactSPHEREBOX( *this, neighbor ) );
	  break;
	    
	default:
	  break;
      }  	      
      break;
	
    case DISC2D:
      switch( convexB->getConvexType() )	
      {
	case DISC2D:
	  return ( isContactSPHERESPHERE( *this, neighbor ) );
	  break;
	    
	case BOX:
	  return ( isContactSPHEREBOX( *this, neighbor ) );
	  break;
	    
	default:
	  break;
      }  	      
      break;        

    case BOX:
      switch( convexB->getConvexType() )	
      {
	case SPHERE:
	  return ( isContactSPHEREBOX( *this, neighbor ) );
	  break;
	    
	case DISC2D:
	  return ( isContactSPHEREBOX( *this, neighbor ) );
	  break;
	    
	default:
	  break;  	      
      }  	      
      break;  
	
    default:
      break;    
  }

  // Bounding Volume Check
  Vector3 gcagcb = *m_transform.getOrigin() -
    	*neighbor.m_transform.getOrigin();
  bool BVCheck = ( Norm(gcagcb) < m_circumscribedRadius + 
	neighbor.m_circumscribedRadius );
  if ( GrainsExec::m_colDetBoundingVolume && BVCheck )
    BVCheck = isContactBVolume( *this, neighbor );
  
  // Check bounding volume overlap only
  if ( GrainsExec::m_InsertionWithBVonly )
    contact = BVCheck;
  // General case with GJK
  // Comment: GJK has consistently shown accuracy issues when 2 particles
  // overlap a lot. Instead of returning a distance of zero to machine 
  // precision, it returns a small number that scales with the size of the 
  // particle. Consequently, some particles are mistakenly inserted in the 
  // simulation. This requires a fix in the future  
  else
  { 
    ++GrainsExec::m_nb_GJK_narrow_collision_detections;

    if ( BVCheck )
    {
      // In case one rigid body is a rectangle
      if ( convexA->getConvexType() == RECTANGLE2D ||
	convexB->getConvexType() == RECTANGLE2D )
      {
        --GrainsExec::m_nb_GJK_narrow_collision_detections;
        pc = ClosestPointRECTANGLE( *this, neighbor, false );
        if ( pc.getOverlapDistance() < 0. ) contact = true;	
      }
      else
      {            
        ++GrainsExec::m_nb_GJK_calls;
      
        Point3 pointA, pointB;
        int nbIterGJK = 0;
        Transform const* a2w = this->getTransformWithCrust();
        Transform const* b2w = neighbor.getTransformWithCrust();  
   
        double distanceMin = (*this).m_crustThickness 
		+ neighbor.m_crustThickness - EPSILON;

        double distance = closest_points( *m_convex, *(neighbor.m_convex), 
		*a2w, *b2w, pointA, pointB, nbIterGJK );

        if ( distance < distanceMin ) contact = true;
      }
    }
  }

  return ( contact );
}




// ----------------------------------------------------------------------------
// Returns whether there is geometric contact with another rigid body
// in the sense of ClosestPoint when the 2 rigid bodies are spheres, i.e., a
// SPHERE-SPHERE contact
bool isContactSPHERESPHERE( RigidBodyWithCrust const& rbA,
	RigidBodyWithCrust const& rbB )
{
  bool contact = false;

  Point3 const* pointA  = rbA.getTransform()->getOrigin();
  Point3 const* pointB  = rbB.getTransform()->getOrigin();

  Vector3 vecteurAB = *pointB - *pointA;
  double  rayonA    = rbA.getCircumscribedRadius();
  double  rayonB    = rbB.getCircumscribedRadius();

  double  distance  = Norm( vecteurAB ) - ( rayonA + rayonB );
  if ( distance < 0. ) contact = true;

  return ( contact );
}




// ----------------------------------------------------------------------------
// Returns whether there is geometric contact with another rigid body
// in the sense of ClosestPoint when 1 rigid body is a sphere and the other
// rigid body is a box, i.e., a SPHERE-BOX contact
bool isContactSPHEREBOX( RigidBodyWithCrust const& rbA,
	RigidBodyWithCrust const& rbB )
{
  bool contact = false;

  const Convex* convexA = rbA.getConvex();
  const Convex* convexB = rbB.getConvex();
  double overlap = 0.;
  Point3 contactPoint;

  if ( convexA->getConvexType() == SPHERE
  	|| convexA->getConvexType() == DISC2D )
  {
    Box const* convexBoxB = (Box const*)(convexB);
    Point3 const* pointA  = rbA.getTransform()->getOrigin();
    double rayonA = rbA.getCircumscribedRadius();
    const Transform* transfB = rbB.getTransform();
    Transform w2b;
    w2b.setToInverseTransform( *transfB );
    contactPoint = convexBoxB->IntersectionPointSPHERE( w2b(*pointA), rayonA,
    	overlap, false );
    if ( overlap < 0. ) contact = true;
  }
  else
  {
    Box const* convexBoxA = (Box const*)(convexA);
    Point3 const* pointB  = rbB.getTransform()->getOrigin();
    double rayonB = rbB.getCircumscribedRadius();
    const Transform* transfA = rbA.getTransform();
    Transform w2a;
    w2a.setToInverseTransform( *transfA );
    contactPoint = convexBoxA->IntersectionPointSPHERE( w2a(*pointB), rayonB,
    	overlap, false );
    if ( overlap < 0. ) contact = true;
  }

  return ( contact );
}




// ----------------------------------------------------------------------------
// Writes the rigid body's "static" data, i.e., the convex geometric
// description only (without any transformation)
void RigidBodyWithCrust::writeStatic( ostream& fileOut ) const
{
  fileOut << *m_convex << endl;
  fileOut << "*CrustThickness " << m_crustThickness;
}




// ----------------------------------------------------------------------------
// Returns whether the rigid body is close to another rigid body in
// the sense of whether their respective bounding boxes plus their crust
// thickness overlap. Slightly different from isClose, needs further
// clarification.
bool intersect( RigidBodyWithCrust const& a, RigidBodyWithCrust const& b )
{
  return ( intersect( a.BoxRigidBody(), b.BoxRigidBody() ) );
}




// ----------------------------------------------------------------------------
// Sets the boolean that tells that the rigid body's transformation
// with the scaling by the crust thickness to shrink the rigid bodies has
// already been computed to false
void RigidBodyWithCrust::initialize_transformWithCrust_to_notComputed()
{
  m_transformWithCrust_computed = false;
}




// ----------------------------------------------------------------------------
// Returns the features of the contact when the 2 rigid bodies are cylinders
PointContact ClosestPointCYLINDERS( RigidBodyWithCrust const& rbA,
                                    RigidBodyWithCrust const& rbB )
{
  try 
  {
    Cylinder const* cylA = (Cylinder const*)( rbA.getConvex() );
    Cylinder const* cylB = (Cylinder const*)( rbB.getConvex() );
    Transform const* a2w = rbA.getTransform();
    Transform const* b2w = rbB.getTransform();
    Vector3 gcagcb = *( a2w->getOrigin() ) - *( b2w->getOrigin() );
    if ( Norm(gcagcb) < rbA.getCircumscribedRadius() +
                      rbB.getCircumscribedRadius() )
      return( intersect( *cylA, *cylB, *a2w, *b2w ) );
    else
      return ( PointNoContact );
  }
  catch ( ContactError const& ) { throw ContactError(); }
}




// ----------------------------------------------------------------------------
// Returns whether there is a contact between the bounding volumes of two rigid
// bodies
bool isContactBVolume( RigidBodyWithCrust const& rbA, 
	RigidBodyWithCrust const& rbB )
{
  return ( isContact( rbA.getBVolume(), rbB.getBVolume(),*rbA.getTransform(),
	*rbB.getTransform() ) );
}




// ----------------------------------------------------------------------------
// Returns the features of the contact when the 1 rigid body is a rectangle
PointContact ClosestPointRECTANGLE( RigidBodyWithCrust const& rbA,
  RigidBodyWithCrust const& rbB, bool const& checkoverlap,
  bool checkCGInRec )
{
  // Note: the rectangle has conceptually no width and therefore no crust 
  // thickness, so instead of using the sum of the crust thicknesses of the 
  // two rigid bodies as the maximum allowed overlap, we use twice the crust 
  // thickness of the rigid body that is not the rectangle.
  // If the overlap is larger than this maximum allowed overlap, we output an 
  // error message and saturate the overlap to the maximum allowed overlap but 
  // do not throw any exception as there is nothing to do at the computing 
  // level.

  try 
  {
    Convex const* convexA = rbA.getConvex();
    Convex const* convexB = rbB.getConvex();
    Transform const* a2w = rbA.getTransform();
    Transform const* b2w = rbB.getTransform();
    double overlap = 0.;
    bool proj = false;

    if ( convexA->getConvexType() == RECTANGLE2D )
    {
      Point3 const* rPt = a2w->getOrigin(); // rectangle center
      Point3 cPt = *( b2w->getOrigin() ); // center of the convex body
      Vector3 rNorm = a2w->getBasis() * Vector3( 0., 0., 1. ); // rect normal
      rNorm.normalized();
      rNorm = copysign( 1., rNorm * ( cPt - *rPt ) ) * rNorm;
      Point3 pointA = (*b2w) 
                      ( convexB->support( ( -rNorm ) * b2w->getBasis() ) );
      if ( ( rNorm * ( pointA - *rPt ) ) < 0. )
      {
        Point3 pointB = ( ( *rPt - pointA ) * rNorm ) * rNorm + pointA;

        // The projection point lies in the rectangle?
        Transform invTransform;
        invTransform.setToInverseTransform( *a2w );      
        if ( convexA->isIn( ( invTransform )( pointB ) ) ) proj = true;
        // The projection of the center of mass lies in the rectangle?	
	else if ( checkCGInRec )
	{
	  Point3 pointC = ( ( *rPt - cPt ) * rNorm ) * rNorm + cPt;     
          if ( convexA->isIn( ( invTransform )( pointC ) ) ) proj = true;
	}
	
	if ( proj )
        {
          Point3 contact = pointA / 2.0 + pointB / 2.0;
          Vector3 overlap_vector = pointA - pointB;
          overlap = - Norm( overlap_vector );

	        if ( checkoverlap )
          {
            double rdwB = rbB.getCrustThickness();	  
            if ( - overlap >= 2. * rdwB )
            {
              cout << "ERR RigidBodyWithCrust::ClosestPointRECTANGLE on "
	      	<< "Processor "
		<< ( GrainsExec::m_MPI ? GrainsExec::getComm()->get_rank() : 0 )
		<< ": " << - overlap << " & " << 2. * rdwB << endl;
              overlap = - 2. * rdwB;
            }
          }	  

          return ( PointContact( contact, overlap_vector, overlap, 0 ) );
        }
      }
    }
    else
    {
      Point3 const* rPt = b2w->getOrigin(); // rectangle center
      Point3 cPt = *( a2w->getOrigin() ); // center of the convex body
      Vector3 rNorm = b2w->getBasis() * Vector3( 0., 0., 1. ); // rect normal
      rNorm.normalized();      
      rNorm = copysign( 1., rNorm * ( cPt - *rPt ) ) * rNorm;
      Point3 pointA = (*a2w)
                      ( convexA->support( ( -rNorm ) * a2w->getBasis() ) );

      if ( ( rNorm * ( pointA - *rPt ) ) < 0. )
      {
        Point3 pointB = ( ( *rPt - pointA ) * rNorm ) * rNorm + pointA;

        // The projection point lies in the rectangle?
        Transform invTransform;
        invTransform.setToInverseTransform( *b2w );
        if ( convexB->isIn( ( invTransform )( pointB ) ) ) proj = true;
	else if ( checkCGInRec )
	{
	  Point3 pointC = ( ( *rPt - cPt ) * rNorm ) * rNorm + cPt;     
          if ( convexB->isIn( ( invTransform )( pointC ) ) ) proj = true;
	}
		
	if ( proj )
        {
          Point3 contact = pointA / 2.0 + pointB / 2.0;
          Vector3 overlap_vector = pointB - pointA;
          overlap = - Norm( overlap_vector );
	  
          if ( checkoverlap )
          {
            double rdwA = rbA.getCrustThickness();
            if ( - overlap >= 2. * rdwA )
            {
              cout << "ERR RigidBodyWithCrust::ClosestPointRECTANGLE on "
	      	"Processor "
		<< ( GrainsExec::m_MPI ? GrainsExec::getComm()->get_rank() : 0 )
		<< ": " << - overlap << " & " << 2. * rdwA << endl;
              overlap = - 2. * rdwA;		
            }
          }

          return ( PointContact( contact, overlap_vector, overlap, 0 ) );
        }
      }
    }

    return ( PointNoContact );
  }
  catch ( ContactError const& ) { throw; }
}
