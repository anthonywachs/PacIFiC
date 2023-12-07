#include "GrainsMPIWrapper.hh"
#include "SpheroCylinder.hh"
#include "ContactBuilderFactory.hh"
#include "Memento.hh"
#include "KinematicsBuilderFactory.hh"
#include "GrainsBuilderFactory.hh"
#include "PointC.hh"
#include "RigidBody.hh"
#include "GrainsExec.hh"
#include "ContactForceModel.hh"
#include "Cylinder.hh"
#include "Sphere.hh"
#include <iterator>
#include <algorithm>


// ----------------------------------------------------------------------------
// Constructor with autonumbering as input parameter
SpheroCylinder::SpheroCylinder( bool const& autonumbering )
  : CompositeParticle( autonumbering )
{
  m_specific_composite_shape = "SpheroCylinder";
  m_height = 0.;
  m_radius = 0.;
}





// ----------------------------------------------------------------------------
// Constructor with an XML node as an input parameter. This constructor is
// expected to be used for reference composite particles
SpheroCylinder::SpheroCylinder( DOMNode* root,
	bool const& autonumbering, int const& pc )
  : CompositeParticle( autonumbering )
{
  m_specific_composite_shape = "SpheroCylinder";

  // Geometric type
  m_GeomType = pc;

  // The composite particle does not have a shape per se, its shape is made
  // of the glued elementary particles. Hence we defines its shape by a point
  // corresponding to its center of mass position (same as in CompositeObstacle)
  m_geoRBWC = new RigidBodyWithCrust( new PointC(), Transform() );

  // Create kinematics
  m_kinematics = KinematicsBuilderFactory::create(
  	m_geoRBWC->getConvex() );

  // Particle density
  if ( ReaderXML::hasNodeAttr( root, "Density" ) )
    m_density = ReaderXML::getNodeAttr_Double( root, "Density" );

  // Height (of the cylinder only), radius, mass, weight and crust thickness
  DOMNode* nGeometry = ReaderXML::getNode( root, "Geometry" );
  m_height = ReaderXML::getNodeAttr_Double( nGeometry, "Height" );
  m_radius = ReaderXML::getNodeAttr_Double( nGeometry, "Radius" ); 
  m_mass = m_density * PI * m_radius * m_radius * 
  	( m_height + ( 4. / 3. ) * m_radius );
  computeWeight();
  double crust_thickness = 
  	ReaderXML::getNodeAttr_Double( nGeometry, "CrustThickness" );
  m_geoRBWC->setCrustThickness( crust_thickness ); 

  // Material
  DOMNode* material_ = ReaderXML::getNode( root, "Material" );
  if ( material_ )
  {
    m_materialName = ReaderXML::getNodeValue_String( material_ );
    ContactBuilderFactory::defineMaterial( m_materialName, false );
  }

  // Angular position of the composite particle
  m_geoRBWC->getTransform()->load( root );

  // Moment of inertia tensor of the spherocylinder
  m_inertia[1] = m_inertia[2] = m_inertia[4] = 0.;
  m_inertia[0] = m_inertia[5] = 
  	( 1. / 12. ) * PI * m_height * m_radius * m_radius
  	 	* ( 3. * m_radius * m_radius + m_height * m_height 
			+ 4. * m_radius * m_height ) 
	+ ( 8. / 15. ) * PI * pow( m_radius, 5. );
  m_inertia[3] = ( 1. / 30. ) * PI * pow( m_radius, 4. )
  	* ( 16. * m_radius + 15. * m_height );  
  cout << m_mass / m_density << " " << m_inertia[0] << " " << m_inertia[3]
  	<< endl;
  BuildInertia();


  // Number of elementary particles
  m_nbElemPart = 3;

  // Allocate containers that scale with the number of elementary particles
  Particle* ppp = NULL;
  Matrix ttt;
  m_elementaryParticles.reserve( m_nbElemPart );
  m_InitialRelativePositions.reserve( m_nbElemPart );
  m_InitialRotationMatrices.reserve( m_nbElemPart );
  for ( size_t j=0; j<m_nbElemPart; ++j )
  {
    m_elementaryParticles.push_back( ppp );
    m_InitialRelativePositions.push_back( Vector3Nul );
    m_InitialRotationMatrices.push_back( ttt );
  }

  // Cylinder
  Cylinder* cyl = new Cylinder( m_radius, m_height );
  RigidBodyWithCrust* geoRBWC_cyl = new RigidBodyWithCrust( cyl, Transform(),
  	false, crust_thickness ); 
  m_elementaryParticles[0] = new Particle( geoRBWC_cyl, m_density,
      	m_materialName, false, pc );
  m_InitialRelativePositions[0][X] = 0.;
  m_InitialRelativePositions[0][Y] = 0.;  
  m_InitialRelativePositions[0][Z] = 0.;  
  m_elementaryParticles[0]->setPosition( m_InitialRelativePositions[0] );  
  m_elementaryParticles[0]->setMasterParticle( this );
  m_InitialRotationMatrices[0] = m_elementaryParticles[0]->getRigidBody()
	->getTransform()->getBasis();
  m_elementaryParticles[0]->getRigidBody()
 	->composeLeftByTransform( *(m_geoRBWC->getTransform()) );

  // Top sphere
  Sphere* spht = new Sphere( m_radius );
  RigidBodyWithCrust* geoRBWC_spht = 
  	new RigidBodyWithCrust( spht, Transform(),
  	false, crust_thickness );
  m_elementaryParticles[1] = new Particle( geoRBWC_spht, m_density,
      	m_materialName, false, pc );
  m_InitialRelativePositions[1][X] = 0.;
  m_InitialRelativePositions[1][Y] = m_height / 2.;  
  m_InitialRelativePositions[1][Z] = 0.;  
  m_elementaryParticles[1]->setPosition( m_InitialRelativePositions[1] );  
  m_elementaryParticles[1]->setMasterParticle( this );
  m_InitialRotationMatrices[1] = m_elementaryParticles[1]->getRigidBody()
	->getTransform()->getBasis();
  m_elementaryParticles[1]->getRigidBody()
 	->composeLeftByTransform( *(m_geoRBWC->getTransform()) );

  // Bottom sphere
  Sphere* sphb = new Sphere( m_radius );
  RigidBodyWithCrust* geoRBWC_sphb = 
  	new RigidBodyWithCrust( sphb, Transform(),
  	false, crust_thickness );
  m_elementaryParticles[2] = new Particle( geoRBWC_sphb, m_density,
      	m_materialName, false, pc );
  m_InitialRelativePositions[2][X] = 0.;
  m_InitialRelativePositions[2][Y] = - m_height / 2.;  
  m_InitialRelativePositions[2][Z] = 0.;  
  m_elementaryParticles[2]->setPosition( m_InitialRelativePositions[2] );  
  m_elementaryParticles[2]->setMasterParticle( this );
  m_InitialRotationMatrices[2] = m_elementaryParticles[1]->getRigidBody()
	->getTransform()->getBasis();
  m_elementaryParticles[2]->getRigidBody()
 	->composeLeftByTransform( *(m_geoRBWC->getTransform()) );

  // Compute and set the the circumscribed radius
  setCircumscribedRadius();

  // In case part of the particle acceleration computed explicity
  if ( Particle::m_splitExplicitAcceleration ) createVelocityInfosNm1();
}





// ----------------------------------------------------------------------------
// Constructor with input parameters
SpheroCylinder::SpheroCylinder( int const& id_,
	Particle const* ParticleRef,
	double const& vx, double const& vy, double const& vz,
	double const& qrotationx, double const& qrotationy,
	double const& qrotationz, double const& qrotations,
	double const& rx, double const& ry, double const& rz,
	const double m[12],
	ParticleActivity const& activ,
	int const& tag_,
	int const& coordination_number_ ,
 	bool const& updatePosition )
  : CompositeParticle( id_, ParticleRef,
	vx, vy, vz, qrotationx, qrotationy, qrotationz, qrotations,
	rx, ry, rz, m, activ, tag_, coordination_number_ , updatePosition )	
{
  m_specific_composite_shape = "SpheroCylinder";
}




// ----------------------------------------------------------------------------
// Constructor with input parameters
SpheroCylinder::SpheroCylinder( int const& id_,
	Particle const* ParticleRef,
	Vector3 const& vtrans,
	Quaternion const& qrot,
	Vector3 const& vrot,
	Transform const& config,
	ParticleActivity const& activ )
  : CompositeParticle( id_, ParticleRef, vtrans, qrot, vrot, config, activ )
{
  m_specific_composite_shape = "SpheroCylinder";
  SpheroCylinder const* SpheroCylRef =
  	dynamic_cast<SpheroCylinder const*>(ParticleRef);
  m_height = SpheroCylRef->m_height;
  m_radius = SpheroCylRef->m_radius;
}




// ----------------------------------------------------------------------------
// Destructor
SpheroCylinder::~SpheroCylinder()
{}




// ----------------------------------------------------------------------------
// Copy constructor (the torsor is initialized to 0)
SpheroCylinder::SpheroCylinder( SpheroCylinder const& other )
  : CompositeParticle( other )
{
  m_height = other.m_height;
  m_radius = other.m_radius;
}




// ----------------------------------------------------------------------------
// Creates a clone of the particle. This method calls the standard
// copy constructor and is used for new particles to be inserted in the
// simulation. Numbering is automatic, total number of components is
// incremented by 1 and activity is set to WAIT. The calling object is
// expected to be a reference particle
Particle* SpheroCylinder::createCloneCopy() const
{
  Particle* particle = new SpheroCylinder( *this );

  return ( particle );
}




// ----------------------------------------------------------------------------
// Creates a clone of the composite particle. This method calls the
// constructor SpheroCylinder( int const& id_, Particle const* ParticleRef,
// Vector3 const& vtrans, Quaternion const& qrot, Vector3 const& vrot,
// Transform const& config, ParticleActivity const& activ ) and is used for
// periodic clone composite particles to be inserted in the simulation.
// Numbering is set with the parameter id_ and total number of components left
// unchanged.
Particle* SpheroCylinder::createCloneCopy( int const& id_,
    	Particle const* ParticleRef, Vector3 const& vtrans,
	Quaternion const& qrot,	Vector3 const& vrot,
	Transform const& config, ParticleActivity const& activ ) const
{
  Particle* particle = new SpheroCylinder( id_, ParticleRef, vtrans,
	qrot, vrot, config, activ );

  return ( particle );
}




// ----------------------------------------------------------------------------
// Returns the number of points to write the composite particle in a
// Paraview format
int SpheroCylinder::numberOfPoints_PARAVIEW() const
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
int SpheroCylinder::numberOfCells_PARAVIEW() const
{
  int nbcells = 0 ;
  for ( size_t i=0; i<m_nbElemPart; ++i )
    nbcells += m_elementaryParticles[i]->getRigidBody()->getConvex()
	->numberOfCells_PARAVIEW();

  return ( nbcells );
}




// ----------------------------------------------------------------------------
// Writes the points describing the composite particle in a Paraview format
void SpheroCylinder::write_polygonsPts_PARAVIEW( ostream& f,
	Vector3 const* translation )const
{
  for ( size_t i=0; i<m_nbElemPart; ++i )
    m_elementaryParticles[i]->write_polygonsPts_PARAVIEW( f, translation ) ;
}




// ----------------------------------------------------------------------------
// Returns a list of points describing the component in a Paraview format
list<Point3> SpheroCylinder::get_polygonsPts_PARAVIEW(
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
// Writes the points describing the composite particle in a
// Paraview format with a transformation that may be different than the current
// transformation of the particle
void SpheroCylinder::write_polygonsPts_PARAVIEW( ostream& f,
	Transform const& transform, Vector3 const* translation ) const
{
  for ( size_t i=0; i<m_nbElemPart; ++i )
    m_elementaryParticles[i]->write_polygonsPts_PARAVIEW( f,
	transform, translation );
}




// ----------------------------------------------------------------------------
// Writes the composite particle in a Paraview format
void SpheroCylinder::write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const
{
  for ( size_t i=0; i<m_nbElemPart; ++i )
    m_elementaryParticles[i]->write_polygonsStr_PARAVIEW( connectivity,
	    offsets, cellstype, firstpoint_globalnumber, last_offset );
}






// ----------------------------------------------------------------------------
// Returns the number of corners of the rigib body shape and a code
// describing the rigid body shape
int SpheroCylinder::getNbCorners() const
{
  return ( 101 );
}





// ----------------------------------------------------------------------------
// Outputs information to be transferred to the fluid
void SpheroCylinder::writePositionInFluid( ostream& fluid )
{
  // TO DO
}
