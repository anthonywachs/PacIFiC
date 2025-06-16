#include "GrainsMPIWrapper.hh"
#include "CompositeParticle.hh"
#include "ContactBuilderFactory.hh"
#include "Memento.hh"
#include "KinematicsBuilderFactory.hh"
#include "GrainsBuilderFactory.hh"
#include "PointC.hh"
#include "RigidBody.hh"
#include "GrainsExec.hh"
#include "ContactForceModel.hh"
#include <iterator>
#include <algorithm>


// ----------------------------------------------------------------------------
// Constructor with autonumbering as input parameter
CompositeParticle::CompositeParticle( bool const& autonumbering ):
  Particle( autonumbering )
{}





// ----------------------------------------------------------------------------
// Constructor with an XML node as an input parameter. This constructor is
// expected to be used for reference composite particles
CompositeParticle::CompositeParticle( DOMNode* root, int const& pc )
  : Particle( false )
{
  // Geometric type
  m_GeomType = pc;

  // The composite particle does not have a shape per se, its shape is made
  // of the glued elementary particles. Hence we defines its shape by a point
  // corresponding to its center of mass position (same as in CompositeObstacle)
  m_geoRBWC = new RigidBodyWithCrust( new PointC(), Transform(), true,
  	EPSILON );

  // Create kinematics
  m_kinematics = KinematicsBuilderFactory::create(
  	m_geoRBWC->getConvex() );

  // Particle density
  if ( ReaderXML::hasNodeAttr( root, "Density" ) )
    m_density = ReaderXML::getNodeAttr_Double( root, "Density" );

  // Particle volume, mass and weight
  DOMNode* nVolume = ReaderXML::getNode( root, "Volume" );
  double vol = ReaderXML::getNodeAttr_Double( nVolume, "Value" );
  m_mass = m_density * vol;
  computeWeight();

  // Material
  DOMNode* material_ = ReaderXML::getNode( root, "Material" );
  if ( material_ )
  {
    m_materialName = ReaderXML::getNodeValue_String( material_ );
    ContactBuilderFactory::defineMaterial( m_materialName, false );
  }

  // Angular position of the composite particle
  m_geoRBWC->getTransform()->load( root );

  // Moment of inertia tensor of the composite
  DOMNode* nMomentOfInertiaTensor =
  	ReaderXML::getNode( root, "MomentOfInertiaTensor" );
  DOMNode* nIxx =
  	ReaderXML::getNode( nMomentOfInertiaTensor, "Ixx" );
  m_inertia[0] = ReaderXML::getNodeAttr_Double( nIxx, "Value" );
  DOMNode* nIxy =
  	ReaderXML::getNode( nMomentOfInertiaTensor, "Ixy" );
  m_inertia[1] = ReaderXML::getNodeAttr_Double( nIxy, "Value" );
  DOMNode* nIxz =
  	ReaderXML::getNode( nMomentOfInertiaTensor, "Ixz" );
  m_inertia[2] = ReaderXML::getNodeAttr_Double( nIxz, "Value" );
  DOMNode* nIyy =
  	ReaderXML::getNode( nMomentOfInertiaTensor, "Iyy" );
  m_inertia[3] = ReaderXML::getNodeAttr_Double( nIyy, "Value" );
  DOMNode* nIyz =
  	ReaderXML::getNode( nMomentOfInertiaTensor, "Iyz" );
  m_inertia[4] = ReaderXML::getNodeAttr_Double( nIyz, "Value" );
  DOMNode* nIzz =
  	ReaderXML::getNode( nMomentOfInertiaTensor, "Izz" );
  m_inertia[5] = ReaderXML::getNodeAttr_Double( nIzz, "Value" );
  BuildInertia();

  // Elementary particles
  DOMNode* nElementaryParticles =
  	ReaderXML::getNode( root, "ElementaryParticles" );
  DOMNodeList* allElemParticles = ReaderXML::getNodes( nElementaryParticles );

  // Number of elementary particles
  m_nbElemPart = allElemParticles->getLength();

  // Allocate containers that scale with the number of elementary particles
  Particle* ppp = NULL;
  Matrix ttt;
  m_elementaryParticles.reserve( m_nbElemPart );
  m_InitialRelativePositions.reserve( m_nbElemPart );
  m_InitialRotationMatrices.reserve( m_nbElemPart );
  for ( size_t j=0; j<m_nbElemPart; ++j )
  {
    m_elementaryParticles.push_back( ppp );
    m_InitialRelativePositions.push_back( Vector3Null );
    m_InitialRotationMatrices.push_back( ttt );
  }

  // Read the elementary particles
  for ( XMLSize_t i=0; i<m_nbElemPart; ++i )
  {
    DOMNode* nElemParticle = allElemParticles->item( i );

    // Construction of the elementary particle
    m_elementaryParticles[i] = new Particle( nElemParticle, pc );

    // Set its material as the composite material
    m_elementaryParticles[i]->setMaterial( m_materialName );

    // Set its density as the composite density
    m_elementaryParticles[i]->setDensity( m_density );

    // Initial relative position
    DOMNode* nRelPos  = ReaderXML::getNode( nElemParticle,
	"RelativePosition" );
    m_InitialRelativePositions[i][X] =
    	ReaderXML::getNodeAttr_Double( nRelPos, "X" );
    m_InitialRelativePositions[i][Y] =
    	ReaderXML::getNodeAttr_Double( nRelPos, "Y" );
    m_InitialRelativePositions[i][Z] =
    	ReaderXML::getNodeAttr_Double( nRelPos, "Z" );
    m_elementaryParticles[i]->setPosition( m_InitialRelativePositions[i] );

    // Set the composite as the master particle
    m_elementaryParticles[i]->setMasterParticle( this );
    
    // Set its ID number (always > 0 for particles )
    m_elementaryParticles[i]->setID( int(i) + 1 );
  }

  // Set crust thickness of the composite as the minimum of the crust
  // thickness of elementary particles
  double minCT = 1.e20;
  for ( size_t j=0; j<m_nbElemPart; ++j )
    minCT = min( minCT, m_elementaryParticles[j]->getCrustThickness() );
  m_geoRBWC->setCrustThickness( minCT );

  // Save the initial rotation matrix of elementary particles
  for ( size_t j=0; j<m_nbElemPart; ++j )
    m_InitialRotationMatrices[j] = m_elementaryParticles[j]->getRigidBody()
	->getTransform()->getBasis();

  // Apply the initial rotation of the composite particle to its elementary
  // particles
  for ( size_t i=0; i<m_elementaryParticles.size(); ++i )
    m_elementaryParticles[i]->getRigidBody()
	->composeLeftByTransform( *(m_geoRBWC->getTransform()) );

  // Compute and set the the circumscribed radius
  setCircumscribedRadius();
  
  // Compute and set the non-spherical bounding volume
  if ( GrainsExec::m_colDetBoundingVolume ) createBoundingVolume();

  // In case part of the particle acceleration computed explicity
  if ( Particle::m_splitExplicitAcceleration ) createVelocityInfosNm1();
}





// ----------------------------------------------------------------------------
// Constructor with input parameters
CompositeParticle::CompositeParticle( int const& id_,
	Particle const* ParticleRef,
	double const& vx, double const& vy, double const& vz,
	double const& qrotationx, double const& qrotationy,
	double const& qrotationz, double const& qrotations,
	double const& rx, double const& ry, double const& rz,
	const double m[12],
	ParticleActivity const& activ,
	int const& tag_,
	int const& coordination_number_ )
  : Particle( id_, ParticleRef, vx, vy, vz,
	qrotationx, qrotationy, qrotationz, qrotations,
	rx, ry, rz, m, activ, tag_, coordination_number_ )
{
  // We know that ParticleRef points to a CompositeParticle, such that
  // we can dynamic cast it to actual type and use -> instead of using
  // get methods through virtual typing
  CompositeParticle const* CompParticleRef =
  	dynamic_cast<CompositeParticle const*>(ParticleRef);

  // Creates and sets the elementary particles
  createSetElementaryParticles( CompParticleRef );
}




// ----------------------------------------------------------------------------
// Constructor with input parameters
CompositeParticle::CompositeParticle( int const& id_,
	Particle const* ParticleRef,
	Vector3 const& vtrans,
	Quaternion const& qrot,
	Vector3 const& vrot,
	Transform const& config,
	ParticleActivity const& activ,
     	map< std::tuple<int,int,int>,
     	std::tuple<bool, Vector3, Vector3, Vector3> > const* contactMap )
  : Particle( id_, ParticleRef, vtrans, qrot, vrot, config, activ, contactMap )
{
  // We know that ParticleRef points to a CompositeParticle, such that
  // we can dynamic cast it to actual type and use -> instead of using
  // get methods through virtual typing
  CompositeParticle const* CompParticleRef =
  	dynamic_cast<CompositeParticle const*>(ParticleRef);

  // Creates and sets the elementary particles
  createSetElementaryParticles( CompParticleRef );
}




// ----------------------------------------------------------------------------
// Destructor
CompositeParticle::~CompositeParticle()
{
  for (size_t i=0;i<m_elementaryParticles.size();++i)
    delete m_elementaryParticles[i];
  m_elementaryParticles.clear();
}




// ----------------------------------------------------------------------------
// Copy constructor (the torsor is initialized to 0)
CompositeParticle::CompositeParticle( CompositeParticle const& other, 
    	bool const& autonumbering )
  : Particle( other, autonumbering )
{
  // Elementary particles
  m_nbElemPart = other.m_nbElemPart;
  m_elementaryParticles.reserve( m_nbElemPart );
  Particle* ppp = NULL;
  for ( size_t i=0; i<m_nbElemPart; ++i )
    m_elementaryParticles.push_back( ppp );
  for ( size_t i=0; i<m_nbElemPart; ++i )
  {
    m_elementaryParticles[i] = new Particle(
      *((other.m_elementaryParticles)[i]), false );    
    m_elementaryParticles[i]->setMasterParticle( this );
    m_elementaryParticles[i]->setID( int(i) + 1 );
  }

  // Initial relative positions
  m_InitialRelativePositions.reserve( m_nbElemPart );
  for ( size_t i=0; i<m_nbElemPart; ++i )
    m_InitialRelativePositions.push_back(
    	(other.m_InitialRelativePositions)[i] );

  // Initial rotation matrices
  m_InitialRotationMatrices.reserve( m_nbElemPart );
  for ( size_t i=0; i<m_nbElemPart; ++i )
    m_InitialRotationMatrices.push_back(
    	(other.m_InitialRotationMatrices)[i] );
}




// ----------------------------------------------------------------------------
// Creates a clone of the particle. This method calls the standard
// copy constructor and is used for new particles to be inserted in the
// simulation. Numbering is automatic, total number of components is
// incremented by 1 and activity is set to WAIT. The calling object is
// expected to be a reference particle
Particle* CompositeParticle::createCloneCopy( bool const& autonumbering ) const
{
  Particle* particle = new CompositeParticle( *this, autonumbering );

  return ( particle );
}




// ----------------------------------------------------------------------------
// Creates a clone of the composite particle. This method calls the
// constructor CompositeParticle( int const& id_, Particle const* ParticleRef,
// Vector3 const& vtrans, Quaternion const& qrot, Vector3 const& vrot,
// Transform const& config, ParticleActivity const& activ ) and is used for
// periodic clone composite particles to be inserted in the simulation.
// Autonumbering is set to false and numbering is set with the parameter id_
Particle* CompositeParticle::createCloneCopy( int const& id_,
    	Particle const* ParticleRef, Vector3 const& vtrans,
	Quaternion const& qrot,	Vector3 const& vrot,
	Transform const& config, ParticleActivity const& activ,
	map< std::tuple<int,int,int>,
     	std::tuple<bool, Vector3, Vector3, Vector3> > const* contactMap ) const
{
  Particle* particle = new CompositeParticle( id_, ParticleRef, vtrans,
	qrot, vrot, config, activ, contactMap );

  return ( particle );
}




// ----------------------------------------------------------------------------
// Computes and sets the moment of inertia tensor and its inverse
void CompositeParticle::BuildInertia()
{
  for ( size_t i=0; i<6; ++i ) m_inertia[i] *= m_density;

  double determinant = m_inertia[0] * m_inertia[3] * m_inertia[5]
    - m_inertia[0] * m_inertia[4] * m_inertia[4]
    - m_inertia[5] * m_inertia[1] * m_inertia[1]
    - m_inertia[3] * m_inertia[2] * m_inertia[2]
    + 2. * m_inertia[1] * m_inertia[2] * m_inertia[4];

  m_inertia_1[0] = ( m_inertia[3] * m_inertia[5] - m_inertia[4] * m_inertia[4] )
  	/ determinant;
  m_inertia_1[1] = ( m_inertia[2] * m_inertia[4] - m_inertia[1] * m_inertia[5] )
  	/ determinant;
  m_inertia_1[2] = ( m_inertia[1] * m_inertia[4] - m_inertia[2] * m_inertia[3] )
  	/ determinant;
  m_inertia_1[3] = ( m_inertia[0] * m_inertia[5] - m_inertia[2] * m_inertia[2] )
  	/ determinant;
  m_inertia_1[4] = ( m_inertia[1] * m_inertia[2] - m_inertia[0] * m_inertia[4] )
  	/ determinant;
  m_inertia_1[5] = ( m_inertia[0] * m_inertia[3] - m_inertia[1] * m_inertia[1] )
  	/ determinant;
}




// ----------------------------------------------------------------------------
// Computes and sets the circumscribed radius
void CompositeParticle::setCircumscribedRadius()
{
  double crad = 0.;

  for ( size_t i=0; i<m_nbElemPart; ++i )
    crad = max( crad, Norm( m_InitialRelativePositions[i] )
	      + m_elementaryParticles[i]->getCircumscribedRadius() );

  m_geoRBWC->setCircumscribedRadius( crad );
}




// ----------------------------------------------------------------------------
// Returns whether the component is a composite particle
bool CompositeParticle::isCompositeParticle() const
{
  return ( true );
}




// ----------------------------------------------------------------------------
// Returns the number of points to write the composite particle in a
// Paraview format
int CompositeParticle::numberOfPoints_PARAVIEW() const
{
  int nbpts = 0 ;
  for ( size_t i=0; i<m_nbElemPart; ++i )
    nbpts += m_elementaryParticles[i]->getRigidBody()->getConvex()
	->numberOfPoints_PARAVIEW();

  return ( nbpts );
}




// ----------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the
// composite particle shape in a Paraview format
int CompositeParticle::numberOfCells_PARAVIEW() const
{
  int nbcells = 0 ;
  for ( size_t i=0; i<m_nbElemPart; ++i )
    nbcells += m_elementaryParticles[i]->getRigidBody()->getConvex()
	->numberOfCells_PARAVIEW();

  return ( nbcells );
}




// ----------------------------------------------------------------------------
// Writes the points describing the composite particle in a Paraview format
void CompositeParticle::write_polygonsPts_PARAVIEW( ostream& f,
	Vector3 const* translation )const
{
  for ( size_t i=0; i<m_nbElemPart; ++i )
    m_elementaryParticles[i]->write_polygonsPts_PARAVIEW( f, translation ) ;
}




// ----------------------------------------------------------------------------
// Returns a list of points describing the component in a Paraview format
list<Point3> CompositeParticle::get_polygonsPts_PARAVIEW(
	Vector3 const* translation ) const
{
  list<Point3> ppp, ParaviewPoints;
  list<Point3>::const_iterator itpp;

  for ( size_t i=0; i<m_nbElemPart; ++i )
  {
    ppp = m_elementaryParticles[i]->get_polygonsPts_PARAVIEW( translation );
    for ( itpp=ppp.begin(); itpp!=ppp.end(); ++itpp )
      ParaviewPoints.push_back( *itpp );
  }

  return ( ParaviewPoints );
}




// ----------------------------------------------------------------------------
// Writes the composite particle in a Paraview format
void CompositeParticle::write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const
{
  for ( size_t i=0; i<m_nbElemPart; ++i )
    m_elementaryParticles[i]->write_polygonsStr_PARAVIEW( connectivity,
	    offsets, cellstype, firstpoint_globalnumber, last_offset );
}





// ----------------------------------------------------------------------------
// Returns the number of elementary particles
size_t CompositeParticle::getNbElemPart() const
{
  return ( m_nbElemPart );
}




// ----------------------------------------------------------------------------
// Returns the number of corners of the rigib body shape and a code
// describing the rigid body shape
int CompositeParticle::getNbCorners() const
{
  return ( 999 );
}




// ----------------------------------------------------------------------------
// Sets the angular velocity
void CompositeParticle::setAngularVelocity( Vector3 const& vrot )
{
  Particle::setAngularVelocity( vrot );

  // Note: my impression is that we never use the elementary particle
  // velocity anywhere, so it might be valuable not to update it
  // and save computing time

  for ( size_t i=0; i<m_nbElemPart; ++i )
    m_elementaryParticles[i]->setAngularVelocity( vrot );
}




// ----------------------------------------------------------------------------
// Sets the angular velocity
void CompositeParticle::setAngularVelocity( double const& omx, 
	double const& omy, double const& omz )
{
  Particle::setAngularVelocity( omx, omy, omz );

  // Note: my impression is that we never use the elementary particle
  // velocity anywhere, so it might be valuable not to update it
  // and save computing time

  for ( size_t i=0; i<m_nbElemPart; ++i )
    m_elementaryParticles[i]->setAngularVelocity( omx, omy, omz );
}




// ----------------------------------------------------------------------------
// Sets the translation velocity
void CompositeParticle::setTranslationalVelocity( Vector3 const& vtrans )
{
  Particle::setTranslationalVelocity( vtrans );

  // Note: my impression is that we never use the elementary particle
  // velocity anywhere, so it might be valuable not to update it
  // and save computing time

  // Translational velocity of elementary particles satisfies rigid-body motion
  // so U_elempart = U_composite + om_composite ^ leverarm_at_present_time
  Matrix rota = m_geoRBWC->getTransform()->getBasis() ;
  Vector3 vrot = *(m_kinematics->getAngularVelocity());
  for ( size_t i=0; i<m_nbElemPart; ++i )
  {
    m_elementaryParticles[i]->setTranslationalVelocity(
    	vtrans + ( vrot ^ ( rota * m_InitialRelativePositions[i] ) ) );
  }
}




// ----------------------------------------------------------------------------
// Sets the translation velocity
void CompositeParticle::setTranslationalVelocity( double const& vx, 
	double const& vy, double const& vz )
{
  Particle::setTranslationalVelocity( vx, vy, vz );

  // Note: my impression is that we never use the elementary particle
  // velocity anywhere, so it might be valuable not to update it
  // and save computing time

  // Translational velocity of elementary particles satisfies rigid-body motion
  // so U_elempart = U_composite + om_composite ^ leverarm_at_present_time
  Matrix rota = m_geoRBWC->getTransform()->getBasis() ;
  Vector3 vrot = *(m_kinematics->getAngularVelocity());  
  Vector3 vrotcrossl;
  for ( size_t i=0; i<m_nbElemPart; ++i )
  {
    vrotcrossl = vrot ^ ( rota * m_InitialRelativePositions[i] );
    m_elementaryParticles[i]->setTranslationalVelocity(
    	vx + vrotcrossl[X], vy + vrotcrossl[Y], vz + vrotcrossl[Z] );
  }
}




// ----------------------------------------------------------------------------
// Sets the boolean that tells that the rigid body's transformation
// with the scaling by the crust thickness to shrink the rigid bodies has
// already been computed to false
void CompositeParticle::initialize_transformWithCrust_to_notComputed()
{
  for ( size_t i=0; i<m_nbElemPart; ++i )
    m_elementaryParticles[i]->initialize_transformWithCrust_to_notComputed();
}




// ----------------------------------------------------------------------------
// Returns the volume of the composite particle
double CompositeParticle::getVolume() const
{
  return ( m_mass / m_density );
}




// ----------------------------------------------------------------------------
// Returns the radius of the sphere of volume equivalent to that of the 
// composite particle 
double CompositeParticle::getEquivalentSphereRadius() const
{
  return ( pow( ( 0.75 / PI ) * m_mass / m_density, 1. / 3. ) );
}




// ----------------------------------------------------------------------------
// Sets the elementary particle position given the position of the composite
// particle
void CompositeParticle::setElementaryParticlesPosition()
{
  // Note: we know the composite particle transformation at this time and we 
  // use it to update the elementary particle transformations. However, we do 
  // not have access to the composite particle relative tranformation from t 
  // to t+dt and cannot simply left compose the elementary particle 
  // transformation by the composite particle relative tranformation. Instead 
  // we work with the elementary particle initial transformation defined
  // by their initial relative position with respect to the composite particle
  // center of mass required to be at (0,0,0) and their initial rotation matrix
  // Therefore there is an additional component to the translation coming from
  // the rotation of the composite particle and the elementary particle initial
  // relative position
  Matrix rota = m_geoRBWC->getTransform()->getBasis() ;
  for ( size_t i=0; i<m_nbElemPart; ++i )
  {
    // Angular position = composite particle rotation matrix * initial rotation
    // matrix
    m_elementaryParticles[i]->getRigidBody()->getTransform()->setBasis(
   	rota * m_InitialRotationMatrices[i] );

    // Center of mass position = composite particle center of mass position +
    // composite particle rotation matrix * initial relative position
    m_elementaryParticles[i]->setPosition( *(m_geoRBWC->getCentre())
      + rota * m_InitialRelativePositions[i] );
  }
}




// ----------------------------------------------------------------------------
// Sets the origin of the composite particle's transformation
void CompositeParticle::setPosition( Point3 const& centre )
{
  Particle::setPosition( centre );
  Matrix rota = m_geoRBWC->getTransform()->getBasis() ;
  for ( size_t i=0; i<m_nbElemPart; ++i )
    m_elementaryParticles[i]->setPosition( centre
      + rota * m_InitialRelativePositions[i] );
}




// ----------------------------------------------------------------------------
// Sets the composite particle's transformation with an 1D array of 12
// values (see class Transform for details)
void CompositeParticle::setPosition( double const* pos )
{
  Particle::setPosition( pos );
  setElementaryParticlesPosition();
}




// ----------------------------------------------------------------------------
// Sets the composite particle's transformation with a transformation
void CompositeParticle::setTransform( Transform const& transform_ )
{
  Particle::setTransform( transform_ );
  setElementaryParticlesPosition();
}




// ----------------------------------------------------------------------------
// Returns whether there is geometric contact with another component.
// Note: the other component must not be of the derived type CompositeObstacle
bool CompositeParticle::isContact( Component const* voisin ) const
{
  bool contact = false;

  for ( size_t i=0; i<m_nbElemPart && !contact; ++i )
    if ( voisin->isCompositeParticle() || voisin->isSTLObstacle() )
      contact = voisin->isContact( m_elementaryParticles[i] );
    else
      contact = m_elementaryParticles[i]->isContact( voisin );

  return ( contact );
}




// ----------------------------------------------------------------------------
// Returns whether there is geometric contact with another component accounting
// for crust thickness. Note: the other component must not be of the derived
// type CompositeObstacle
bool CompositeParticle::isContactWithCrust( Component const* voisin ) const
{
  bool contact = false;

  for ( size_t i=0; i<m_nbElemPart && !contact; ++i )
    if ( voisin->isCompositeParticle() || voisin->isSTLObstacle() )
      contact = voisin->isContactWithCrust( m_elementaryParticles[i] );
    else
      contact = m_elementaryParticles[i]->isContactWithCrust( voisin );

  return ( contact );
}




// ----------------------------------------------------------------------------
// Returns whether there is geometric proximity with another
// component in the sense of whether their respective bounding boxes overlap.
// Note: the other component must not be of the derived type CompositeObstacle
bool CompositeParticle::isClose( Component const* voisin ) const
{
  bool contact = false;

  for ( size_t i=0; i<m_nbElemPart && !contact; ++i )
    if ( voisin->isCompositeParticle() || voisin->isSTLObstacle() )
      contact = voisin->isClose( m_elementaryParticles[i] );
    else
      contact = m_elementaryParticles[i]->isClose( voisin );

  return ( contact );
}




// ----------------------------------------------------------------------------
// Returns whether there is geometric proximity with another
// component in the sense of whether their respective bounding boxes minus
// their crust thickness overlap. Note: the other component must not be of the
// derived type CompositeObstacle
bool CompositeParticle::isCloseWithCrust( Component const* voisin ) const
{
  bool contact = false;

  for ( size_t i=0; i<m_nbElemPart && !contact; ++i )
    if ( voisin->isCompositeParticle() || voisin->isSTLObstacle() )
      contact = voisin->isCloseWithCrust( m_elementaryParticles[i] );
    else
      contact = m_elementaryParticles[i]->isCloseWithCrust( voisin );

  return ( contact );
}




// ----------------------------------------------------------------------------
// Checks whether 2 components are in contact. If contact, computes
// contact force & torque and adds to each component torsor
void CompositeParticle::InterAction( Component* voisin,
	double dt, double const& time, LinkedCell* LC )
{
  try 
  {
    list<ContactInfos*> listContactInfos;

    // Search all contact points between the composite particle and the 
    // component and store them in the list listContactInfos
    for ( size_t i=0; i<m_nbElemPart; ++i )
    {
      if ( voisin->isCompositeParticle() || voisin->isSTLObstacle() )
        voisin->SearchContact( m_elementaryParticles[i], dt,
	  time, LC, listContactInfos );
      else m_elementaryParticles[i]->SearchContact( voisin, dt,
	  time, LC, listContactInfos );
    }

    // Loop over all contact points and compute the contact force & torque
    int nbContact = int(listContactInfos.size());
    for ( list<ContactInfos*>::iterator il=listContactInfos.begin();
      il!=listContactInfos.end(); il++ )
    {
      LC->addToContactsFeatures( time, (*il)->ContactPoint );

      if ( ContactBuilderFactory::contactForceModel(
		(*il)->p0->getMaterial(), (*il)->p1->getMaterial() )
      		->computeForces( (*il)->p0, (*il)->p1, (*il)->ContactPoint,
		LC, dt, nbContact ) )
      {
        (*il)->p0->getMasterComponent()->addToCoordinationNumber( 1 );
        (*il)->p1->getMasterComponent()->addToCoordinationNumber( 1 );
      }
      delete *il;
    }

    // Free the list
    listContactInfos.clear();
  }
  catch ( ContactError const& ) { throw; }
}




// ----------------------------------------------------------------------------
// Searches and stores all contact points between a composite particle and
// a component
void CompositeParticle::SearchContact( Component* voisin, double dt,
      double const& time, LinkedCell *LC,
      list<ContactInfos*>& listContact )
{
  try 
  {
    for ( size_t i=0; i<m_nbElemPart; ++i )
      m_elementaryParticles[i]->SearchContact( voisin, dt, time, LC,
	listContact );
  }
  catch ( ContactError const& ) { throw; }
}




// ----------------------------------------------------------------------------
// Compose the component transformation on the left by another
// transformation: this = t o this (this first followed by t)
void CompositeParticle::composePositionLeftByTransform( Transform const& t )
{
  m_geoRBWC->composeLeftByTransform( t );
  setElementaryParticlesPosition();
}




// ----------------------------------------------------------------------------
// Compose the component transformation on the right by another
// transformation: this = this o t (t first followed by this)
void CompositeParticle::composePositionRightByTransform( Transform const& t )
{
  m_geoRBWC->composeRightByTransform( t );
  setElementaryParticlesPosition();
}




// ----------------------------------------------------------------------------
// Solves the Newton's law and move particle to their new position
void CompositeParticle::Move( double time, 
	double const& dt_particle_vel, 
    	double const& dt_particle_disp )
{
  // Move the composite particle
  Particle::Move( time, dt_particle_vel, dt_particle_disp );

  // Update elementary particles' position
  setElementaryParticlesPosition();

  // Update elementary particles' velocity
  // Note: my impression is that we never use the elementary particle
  // velocity anywhere, so it might be valuable not to update it
  // and save computing time
  Matrix rota = m_geoRBWC->getTransform()->getBasis() ;
  Vector3 vrot = *(m_kinematics->getAngularVelocity());
  Vector3 vtrans =  *(m_kinematics->getTranslationalVelocity());
  for ( size_t i=0; i<m_nbElemPart; ++i )
  {
    m_elementaryParticles[i]->setAngularVelocity( vrot );
    m_elementaryParticles[i]->setTranslationalVelocity(
   	vtrans + ( vrot ^ ( rota * m_InitialRelativePositions[i] ) ) );
  }
}




// ----------------------------------------------------------------------------
// Translates the composite particle
void CompositeParticle::Translate( Vector3 const& translation )
{
  Component::Translate( translation );
  for ( size_t i=0; i<m_nbElemPart; ++i )
    m_elementaryParticles[i]->Translate( translation );
}




// ----------------------------------------------------------------------------
// Saves additional features of a (in practice reference) composite particle
// for reload
void CompositeParticle::writeAdditionalFeatures( ostream& fileSave ) const
{
  fileSave << endl << "<ElementaryParticles> " << m_nbElemPart << endl;
  for ( size_t i=0; i<m_nbElemPart; ++i )
  {
    fileSave << "<ElementaryParticle>" << endl;
    m_elementaryParticles[i]->writeStatic( fileSave );
    fileSave << "*RelativePosition" << endl;
    m_InitialRelativePositions[i].writeGroup3( fileSave );
    fileSave << endl << "*InitialRotationMatrix" << endl;
    m_InitialRotationMatrices[i].writeMatrix( fileSave );
    fileSave << endl;
    fileSave << "</ElementaryParticle>" << endl;
  }
  fileSave << "</ElementaryParticles>";
}




// ----------------------------------------------------------------------------
// Reads additional features of a (in practice reference) composite particle
// data from a stream
void CompositeParticle::readAdditionalFeatures( istream& fileIn )
{
  string buffer;
  ParticleKinematics* pkine = NULL;

  // Read the buffer "<ElementaryParticles>" and the number of elementary
  // particles
  fileIn >> buffer >> m_nbElemPart;

  // Allocate containers that scale with the number of elementary particles
  Particle* ppp = NULL;
  Matrix ttt;
  m_elementaryParticles.reserve( m_nbElemPart );
  m_InitialRelativePositions.reserve( m_nbElemPart );
  m_InitialRotationMatrices.reserve( m_nbElemPart );
  for ( size_t j=0; j<m_nbElemPart; ++j )
  {
    m_elementaryParticles.push_back( ppp );
    m_InitialRelativePositions.push_back( Vector3Null );
    m_InitialRotationMatrices.push_back( ttt );
  }

  // Read elementary particles
  for ( size_t j=0; j<m_nbElemPart; ++j )
  {
    // Read the buffer "<ElementaryParticle>"
    fileIn >> buffer;

    // Construct an empty particle
    m_elementaryParticles[j] = new Particle( false );

    // Read static data from stream
    m_elementaryParticles[j]->read( fileIn, true );

    // Create the particle kinematics
    pkine = KinematicsBuilderFactory::create(
    	m_elementaryParticles[j]->getRigidBody()->getConvex() );
    m_elementaryParticles[j]->setKinematics( pkine );

    // Read relative positions
    fileIn >> buffer;
    fileIn >> m_InitialRelativePositions[j];

    // Read initial rotation matrices
    fileIn >> buffer;
    fileIn >> m_InitialRotationMatrices[j];

    // Read the buffer "</ElementaryParticle>"
    fileIn >> buffer;

    // Set its density as the composite density
    m_elementaryParticles[j]->setDensity( m_density );

    // Set the composite as the master particle
    m_elementaryParticles[j]->setMasterParticle( this );

    // Set the activity as the composite activity
    m_elementaryParticles[j]->setActivity( m_activity );
    
    // Set its ID number
    m_elementaryParticles[j]->setID( int(j) + 1 );    
  }

  // Read the buffer "</ElementaryParticles>"
  fileIn >> buffer;

  // Set elementary particle positions
  setElementaryParticlesPosition();

  // Compute and set the the circumscribed radius
  setCircumscribedRadius();
    
  // Compute and set the non-spherical bounding volume
  if ( GrainsExec::m_colDetBoundingVolume ) createBoundingVolume();  

  // Set elementary particle kinematics
  // Note: again not sure we do anything with the elementary particle kinematics
  Matrix rota = m_geoRBWC->getTransform()->getBasis() ;
  Vector3 vrot = *(m_kinematics->getAngularVelocity());
  Vector3 vtrans = *(m_kinematics->getTranslationalVelocity());
  Quaternion qrot = *(m_kinematics->getQuaternionRotation());
  for ( size_t i=0; i<m_nbElemPart; ++i )
  {
    m_elementaryParticles[i]->setTranslationalVelocity( vtrans
      + ( vrot ^ ( rota * m_InitialRelativePositions[i] ) ) );
    m_elementaryParticles[i]->setAngularVelocity( vrot );
    m_elementaryParticles[i]->setQuaternionRotation( qrot );
  }

  // Set crust thickness of the composite as the minimum of the crust
  // thickness of elementary particles
  double minCT = 1.e20;
  for ( size_t j=0; j<m_nbElemPart; ++j )
    minCT = min( minCT, m_elementaryParticles[j]->getCrustThickness() );
  m_geoRBWC->setCrustThickness( minCT );
}




// ----------------------------------------------------------------------------
// Reads composite particle data from a stream. Usage: for standard composite
// particles in the 2014 reload format
void CompositeParticle::read2014( istream& fileIn, vector<Particle*> const*
  	referenceParticles )
{
  Particle::read2014( fileIn, referenceParticles );

  // We know that (*referenceParticles)[m_GeomType] points to a
  // CompositeParticle, such that we can dynamic cast it to actual type and
  // use -> instead of using get methods through virtual typing
  CompositeParticle const* CompParticleRef =
  	dynamic_cast<CompositeParticle const*>(
		(*referenceParticles)[m_GeomType]);

  // Creates and sets the elementary particles
  createSetElementaryParticles( CompParticleRef );
}




// ----------------------------------------------------------------------------
// Reads composite particle data from a stream in a binary form.
// Usage: for standard composite particles in the 2014 reload format
void CompositeParticle::read2014_binary( istream& fileIn,
	vector<Particle*> const* referenceParticles )
{
  Particle::read2014_binary( fileIn, referenceParticles );

  // We know that (*referenceParticles)[m_GeomType] points to a
  // CompositeParticle, such that we can dynamic cast it to actual type and
  // use -> instead of using get methods through virtual typing
  CompositeParticle const* CompParticleRef =
  	dynamic_cast<CompositeParticle const*>(
		(*referenceParticles)[m_GeomType]);

  // Creates and sets the elementary particles
  createSetElementaryParticles( CompParticleRef );
}




// ----------------------------------------------------------------------------
// Creates and sets the elementary particles
void CompositeParticle::createSetElementaryParticles(
	CompositeParticle const* CompParticleRef )
{
  // Elementary particles
  m_nbElemPart = CompParticleRef->m_nbElemPart;
  m_elementaryParticles.reserve( m_nbElemPart );
  Particle* ppp = NULL;
  for ( size_t i=0; i<m_nbElemPart; ++i )
    m_elementaryParticles.push_back( ppp );
  for ( size_t i=0; i<m_nbElemPart; ++i )
  {
    m_elementaryParticles[i] = new Particle(
      *((CompParticleRef->m_elementaryParticles)[i]), false );
    m_elementaryParticles[i]->setMasterParticle( this );
    m_elementaryParticles[i]->setID( int(i) + 1 );
  }

  // Initial relative positions
  m_InitialRelativePositions.reserve( m_nbElemPart );
  for ( size_t i=0; i<m_nbElemPart; ++i )
    m_InitialRelativePositions.push_back(
    	(CompParticleRef->m_InitialRelativePositions)[i] );

  // Initial rotation matrices
  m_InitialRotationMatrices.reserve( m_nbElemPart );
  for ( size_t i=0; i<m_nbElemPart; ++i )
    m_InitialRotationMatrices.push_back(
    	(CompParticleRef->m_InitialRotationMatrices)[i] );

  // Set elementary particle positions
  setElementaryParticlesPosition();

  // Set elementary particle kinematics
  // Note: again not sure we do anything with the elementary particle kinematics
  Matrix rota = m_geoRBWC->getTransform()->getBasis() ;
  Vector3 vrot = *(m_kinematics->getAngularVelocity());
  Vector3 vtrans =  *(m_kinematics->getTranslationalVelocity());
  Quaternion qrot = *(m_kinematics->getQuaternionRotation());
  for ( size_t i=0; i<m_nbElemPart; ++i )
  {
    m_elementaryParticles[i]->setTranslationalVelocity( vtrans
      + ( vrot ^ ( rota * m_InitialRelativePositions[i] ) ) );
    m_elementaryParticles[i]->setAngularVelocity( vrot );
    m_elementaryParticles[i]->setQuaternionRotation( qrot );
  }
}



// ----------------------------------------------------------------------------
// Outputs information to be transferred to the fluid
void CompositeParticle::writePositionInFluid( ostream& fluid )
{
  // TO DO
}




// ----------------------------------------------------------------------------
// Returns whether a point lies inside the composite particle
bool CompositeParticle::isIn( Point3 const& pt ) const
{
  bool bisIn = false;

  for ( size_t i=0; i<m_nbElemPart && !bisIn; ++i )
    bisIn = m_elementaryParticles[i]->isIn( pt );

  return ( bisIn );
}




// ----------------------------------------------------------------------------
// Returns the bounding box of the composite particle
BBox CompositeParticle::BoundingBox() const
{
  BBox compBBox = m_elementaryParticles[0]->BoundingBox();
  for ( size_t i=1; i<m_nbElemPart; ++i )
    compBBox.enclose( compBBox, m_elementaryParticles[i]->BoundingBox() );

  return ( compBBox );
}




// ----------------------------------------------------------------------------
// Returns a pointer to the vector of initial relative positions */
vector<Vector3> const* CompositeParticle::getInitialRelativePositions() const
{
  return ( &m_InitialRelativePositions );
}




// ----------------------------------------------------------------------------
// Sets the rotation quaternion
void CompositeParticle::setQuaternionRotation( double const& vecteur0,
	double const& vecteur1,
	double const& vecteur2,
	double const& scalaire )
{
  Particle::setQuaternionRotation( vecteur0, vecteur1, vecteur2, scalaire );
  for ( size_t i=0; i<m_nbElemPart; ++i )
    m_elementaryParticles[i]->setQuaternionRotation( vecteur0, vecteur1, 
    	vecteur2, scalaire );  
}	
	
	
// ----------------------------------------------------------------------------
// Sets the rotation quaternion
void CompositeParticle::setQuaternionRotation( Quaternion const& qrot )
{
  Particle::setQuaternionRotation( qrot );
  for ( size_t i=0; i<m_nbElemPart; ++i )
    m_elementaryParticles[i]->setQuaternionRotation( qrot );  
}




// ----------------------------------------------------------------------------
// Returns whether two particles are of the same type
bool CompositeParticle::equalType( Particle const* other ) const
{
  bool same = false;
  
  // Check if particle/composite particle matches
  same = other->isCompositeParticle();

  // Check if density matches
  if ( same ) same = fabs( m_density - other->getDensity() ) < LOWEPS;

  // Check if material matches
  if ( same ) same = ( m_materialName == other->getMaterial() );
  
  // Check if geometric features match
  if ( same )
  {
    // We know that ParticleRef points to a CompositeParticle, such that
    // we can dynamically cast it to actual type
    CompositeParticle const* other_ =
  	dynamic_cast<CompositeParticle const*>(other);
	
    // Check if number of elementary particles matches
    same = ( m_nbElemPart == other_->m_nbElemPart ); 
    
    // Check if elementary particles match
    if ( same )
      for ( size_t i=0; i<m_nbElemPart && same; ++i )
        same = m_elementaryParticles[i]->equalType( 
		other_->m_elementaryParticles[i] ); 
		
    // Check if relative positions of elementary particles match
    double lmin = m_geoRBWC->getCircumscribedRadius();
    if ( same )
      for ( size_t i=0; i<m_nbElemPart && same; ++i )
        same = ( Norm( m_InitialRelativePositions[i]
		- other_->m_InitialRelativePositions[i] ) < LOWEPS * lmin ) ;
		
    // Check if initial rotation matrices of elementary particles match
    if ( same )
      for ( size_t i=0; i<m_nbElemPart && same; ++i )
        for ( int j=0; j<int(m_nbElemPart) && same; ++j )  
        same = ( Norm( m_InitialRotationMatrices[i][j]
		- other_->m_InitialRotationMatrices[i][j] ) < LOWEPS * lmin ) ;
  }
  
  return ( same );  
}




// ----------------------------------------------------------------------------
// Creates the non-spherical bounding volume */
void CompositeParticle::createBoundingVolume()
{
  BVolume* bvol = NULL;
  size_t ndim = 3;
  if ( GrainsBuilderFactory::getContext() == DIM_2 ) ndim = 2;
  Convex const* convexA = NULL;
  Transform const* a2w = NULL;
  RigidBodyWithCrust const* rbA = NULL;
  Vector3 u;
  Point3 pointA;
   
  if ( GrainsExec::m_colDetBoundingVolume == 1 ) // OBB
  {
    Vector3 maxcoord( -1.20 );
    for ( size_t i=0; i<m_nbElemPart; ++i )
    {
      rbA = m_elementaryParticles[i]->getRigidBody();
      convexA = rbA->getConvex();
      a2w = rbA->getTransform();
      for ( size_t m=0; m<ndim; ++m )
      {
        u.reset();
	
	// Positive coordinate
	u[m] = 1.;
	pointA = (*a2w)( convexA->support( ( u ) * a2w->getBasis() ) );
	maxcoord[m] = max( maxcoord[m], pointA[m] );
	
	// Negative coordinate
	u[m] = - 1.;
	pointA = (*a2w)( convexA->support( ( u ) * a2w->getBasis() ) );
	maxcoord[m] = max( maxcoord[m], fabs( pointA[m] ) );		
      }
    }

    bvol = new OBB( maxcoord, Matrix() );
  }
  else if ( GrainsExec::m_colDetBoundingVolume == 2 ) // OBC
  {
    size_t const nn = 720;
    double const deltatheta = 2. * PI / double(nn);
    double theta = 0., minVol = 1.e20;
    Vector3 radius, height( -1.20 );
    Vector3 initOri[3];
    size_t coord0, coord1, axisDir = 0;
    for ( size_t m=0; m<3; ++m ) initOri[m][m] = 1.;
    
    for ( size_t i=0; i<m_nbElemPart; ++i )
    {    
      rbA = m_elementaryParticles[i]->getRigidBody();
      convexA = rbA->getConvex();
      a2w = rbA->getTransform();
      for ( size_t m=0; m<3; ++m )
      {      
        u.reset();
      
        // Height
	// Positive coordinate
	u[m] = initOri[m][m];
	pointA = (*a2w)( convexA->support( ( u ) * a2w->getBasis() ) );
	height[m] = max( height[m], 2. * pointA[m] );
	
	// Negative coordinate
	u[m] = - initOri[m][m];
	pointA = (*a2w)( convexA->support( ( u ) * a2w->getBasis() ) );
	height[m] = max( height[m], fabs( 2. * pointA[m] ) );
	
	
	// Radius
	u[m] = 0.;
	coord0 = m == 0 ? 1 : 0;
	coord1 = m == 2 ? 1 : 2;	
	for ( size_t k=0; k<nn; ++k )
	{
	  theta = double(k) * deltatheta;
	  u[coord0] = cos( theta );
	  u[coord1] = sin( theta );
	  pointA = (*a2w)( convexA->support( ( u ) * a2w->getBasis() ) );
	  radius[m] = max( radius[m], Norm( pointA ) );  
	}
      }      
    }
    
    // Select the direction that heads to the smallest volume cylinder
    for ( size_t m=0; m<3; ++m )  
      if ( PI * height[m] * pow( radius[m], 2. ) < minVol ) 
      {  
        minVol = PI * height[m] * pow( radius[m], 2. );
	axisDir = m;
      }
      
    bvol = new OBC( radius[axisDir], height[axisDir], initOri[axisDir] );      
  }
  
  m_geoRBWC->setBoundingVolume( bvol ); 
}




// ----------------------------------------------------------------------------
// Writes the particle in an OBJ format
void CompositeParticle::write_OBJ( ostream& f, size_t const& group_number,
    	size_t& firstpoint_number ) const
{
  for ( size_t i=0; i<m_nbElemPart; ++i )
  {
    m_elementaryParticles[i]->write_OBJ( f, i, firstpoint_number );
    f << endl;
  } 
}
