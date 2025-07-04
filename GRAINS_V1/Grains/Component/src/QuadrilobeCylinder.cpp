#include "GrainsMPIWrapper.hh"
#include "QuadrilobeCylinder.hh"
#include "ContactBuilderFactory.hh"
#include "Memento.hh"
#include "KinematicsBuilderFactory.hh"
#include "GrainsBuilderFactory.hh"
#include "PointC.hh"
#include "RigidBody.hh"
#include "GrainsExec.hh"
#include "ContactForceModel.hh"
#include "SpheroCylindricalPrism.hh"
#include <iterator>
#include <algorithm>

int QuadrilobeCylinder::m_visuNodeNbPerHalf = 16;


// ----------------------------------------------------------------------------
// Constructor with autonumbering as input parameter
QuadrilobeCylinder::QuadrilobeCylinder( bool const& autonumbering )
  : CompositeParticle( autonumbering )
{
  m_specific_composite_shape = "QuadrilobeCylinder";
  m_radius = 0.;
  m_height = 0.;
  m_armLength = 0.;
}





// ----------------------------------------------------------------------------
// Constructor with an XML node as an input parameter. This constructor is
// expected to be used for reference composite particles
QuadrilobeCylinder::QuadrilobeCylinder( DOMNode* root, int const& pc )
  : CompositeParticle( false )
{
  m_specific_composite_shape = "QuadrilobeCylinder";

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

  // Height (of the cylinder only), radius, mass, weight and crust thickness
  DOMNode* nGeometry = ReaderXML::getNode( root, "Geometry" );
  m_radius = ReaderXML::getNodeAttr_Double( nGeometry, "Radius" ); 
  m_height = ReaderXML::getNodeAttr_Double( nGeometry, "Height" );
  m_armLength = ReaderXML::getNodeAttr_Double( nGeometry, "ArmLength" ); 
  m_mass = m_density * ( m_height * pow( m_radius, 2. ) * ( 2. * PI
  	- 4. ) + 8. * m_radius * m_armLength * m_height );
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
  
  // Moment of inertia tensor of the QuadrilobeCylinder
  m_inertia[1] = m_inertia[2] = m_inertia[4] = 0.;
  double h2 = m_height * m_height, r2 = m_radius * m_radius;
  double a1 = 4 * r2 * m_height,
  	a2 = ( m_armLength - m_radius ) * 2. * m_radius * m_height,
	a3 = 0.5 * PI * r2 * m_height;  
  m_inertia[0] = a1 * ( 4. * r2 + h2 ) / 12.
  	+ 0.5 * a2 * ( ( pow( m_armLength - m_radius, 2. ) + h2 ) / 3.   
		+ pow( m_armLength +  m_radius, 2. ) 
		+ ( 4. * r2 + h2 ) / 3. )
	+ 2. * a3 * ( ( 0.25 - 16. / ( 9. * PI * PI ) ) * r2
		+ h2 / 6. + pow( m_armLength + 4. 
		* m_radius / ( 3. * PI ), 2. ) + 0.25 * r2 ) ;
  m_inertia[3] = 2. * a1 * r2 / 3.
  	+ a2 * ( ( pow( m_armLength - m_radius, 2. ) + 4. * r2 ) / 3. 
		+ pow( m_armLength + m_radius, 2.) )
	+ 4. * a3 * ( ( 0.5 - 16. / ( 9. * PI * PI ) ) * r2 
		+ pow( m_armLength + 4. * m_radius / ( 3. * PI ), 2. ) ); 
  m_inertia[5] = m_inertia[0];  
  BuildInertia();


  // Number of elementary particles
  m_nbElemPart = 2;

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

  // First spherocylindrical prism
  SpheroCylindricalPrism* scp1 = new SpheroCylindricalPrism( m_radius, 
  	2. * m_armLength, m_height );
  RigidBodyWithCrust* geoRBWC_scp1 = new RigidBodyWithCrust( scp1, Transform(),
  	false, crust_thickness ); 
  m_elementaryParticles[0] = new Particle( geoRBWC_scp1, m_density,
      	m_materialName, pc );
  m_InitialRelativePositions[0][X] = 0.;
  m_InitialRelativePositions[0][Y] = 0.;  
  m_InitialRelativePositions[0][Z] = 0.;  
  m_elementaryParticles[0]->setPosition( m_InitialRelativePositions[0] );  
  m_elementaryParticles[0]->setMasterParticle( this );
  m_InitialRotationMatrices[0] = m_elementaryParticles[0]->getRigidBody()
	->getTransform()->getBasis();
  m_elementaryParticles[0]->getRigidBody()
 	->composeLeftByTransform( *(m_geoRBWC->getTransform()) );

  // Second spherocylindrical prism
  SpheroCylindricalPrism* scp2 = new SpheroCylindricalPrism( m_radius, 
  	2. * m_armLength, m_height );
  RigidBodyWithCrust* geoRBWC_scp2 = new RigidBodyWithCrust( scp2, Transform(),
  	false, crust_thickness ); 
  m_elementaryParticles[1] = new Particle( geoRBWC_scp2, m_density,
      	m_materialName, pc );
  m_InitialRelativePositions[1][X] = 0.;
  m_InitialRelativePositions[1][Y] = 0.;  
  m_InitialRelativePositions[1][Z] = 0.;  
  m_elementaryParticles[1]->setPosition( m_InitialRelativePositions[1] );  
  m_elementaryParticles[1]->setMasterParticle( this );
  m_InitialRotationMatrices[1].setValue( 
  	cos( 0.5 * PI ), 0., - sin( 0.5 * PI ), 
	0., 1., 0.,
	sin( 0.5 * PI ), 0., cos( 0.5 * PI ) );  
  m_elementaryParticles[1]->getRigidBody()->getTransform()->setBasis(
  	m_InitialRotationMatrices[1] );
  m_elementaryParticles[1]->getRigidBody()
 	->composeLeftByTransform( *(m_geoRBWC->getTransform()) );

  // Set the the circumscribed radius
  setCircumscribedRadius();

  // Compute and set the non-spherical bounding volume
  if ( GrainsExec::m_colDetBoundingVolume ) createBoundingVolume();

  // In case part of the particle acceleration computed explicity
  if ( Particle::m_splitExplicitAcceleration ) createVelocityInfosNm1();
}





// ----------------------------------------------------------------------------
// Constructor with input parameters
QuadrilobeCylinder::QuadrilobeCylinder( int const& id_,
	Particle const* ParticleRef,
	double const& vx, double const& vy, double const& vz,
	double const& qrotationx, double const& qrotationy,
	double const& qrotationz, double const& qrotations,
	double const& rx, double const& ry, double const& rz,
	const double m[12],
	ParticleActivity const& activ,
	int const& tag_,
	int const& coordination_number_ )
  : CompositeParticle( id_, ParticleRef,
	vx, vy, vz, qrotationx, qrotationy, qrotationz, qrotations,
	rx, ry, rz, m, activ, tag_, coordination_number_ )	
{
  m_specific_composite_shape = "QuadrilobeCylinder";
  QuadrilobeCylinder const* QuadrilobeCylinderRef =
  	dynamic_cast<QuadrilobeCylinder const*>(ParticleRef);
  m_radius = QuadrilobeCylinderRef->m_radius; 
  m_armLength = QuadrilobeCylinderRef->m_armLength;    
  m_height = QuadrilobeCylinderRef->m_height;
}




// ----------------------------------------------------------------------------
// Constructor with input parameters
QuadrilobeCylinder::QuadrilobeCylinder( int const& id_,
	Particle const* ParticleRef,
	Vector3 const& vtrans,
	Quaternion const& qrot,
	Vector3 const& vrot,
	Transform const& config,
	ParticleActivity const& activ,
     	map< std::tuple<int,int,int>,
     	std::tuple<bool, Vector3, Vector3, Vector3> > const* contactMap )
  : CompositeParticle( id_, ParticleRef, vtrans, qrot, vrot, config, activ,
  	contactMap )
{
  m_specific_composite_shape = "QuadrilobeCylinder";
  QuadrilobeCylinder const* QuadrilobeCylinderRef =
  	dynamic_cast<QuadrilobeCylinder const*>(ParticleRef);
  m_radius = QuadrilobeCylinderRef->m_radius; 
  m_armLength = QuadrilobeCylinderRef->m_armLength;    
  m_height = QuadrilobeCylinderRef->m_height;
}




// ----------------------------------------------------------------------------
// Destructor
QuadrilobeCylinder::~QuadrilobeCylinder()
{}




// ----------------------------------------------------------------------------
// Copy constructor (the torsor is initialized to 0)
QuadrilobeCylinder::QuadrilobeCylinder( QuadrilobeCylinder const& other, 
    	bool const& autonumbering )
  : CompositeParticle( other, autonumbering )
{
  m_radius = other.m_radius;
  m_armLength = other.m_armLength;
  m_height = other.m_height;
}




// ----------------------------------------------------------------------------
// Creates a clone of the particle. This method calls the standard
// copy constructor and is used for new particles to be inserted in the
// simulation. Numbering is automatic, total number of components is
// incremented by 1 and activity is set to WAIT. The calling object is
// expected to be a reference particle
Particle* QuadrilobeCylinder::createCloneCopy( bool const& autonumbering ) const
{
  Particle* particle = new QuadrilobeCylinder( *this, autonumbering );

  return ( particle );
}




// ----------------------------------------------------------------------------
// Creates a clone of the composite particle. This method calls the
// constructor QuadrilobeCylinder( int const& id_, Particle const* ParticleRef,
// Vector3 const& vtrans, Quaternion const& qrot, Vector3 const& vrot,
// Transform const& config, ParticleActivity const& activ ) and is used for
// periodic clone composite particles to be inserted in the simulation.
// Autonumbering is set to false and numbering is set with the parameter id_
Particle* QuadrilobeCylinder::createCloneCopy( int const& id_,
    	Particle const* ParticleRef, Vector3 const& vtrans,
	Quaternion const& qrot,	Vector3 const& vrot,
	Transform const& config, ParticleActivity const& activ,
	map< std::tuple<int,int,int>,
     	std::tuple<bool, Vector3, Vector3, Vector3> > const* contactMap ) const
{
  Particle* particle = new QuadrilobeCylinder( id_, ParticleRef, vtrans,
	qrot, vrot, config, activ, contactMap );

  return ( particle );
}




// ----------------------------------------------------------------------------
// Returns the number of corners of the rigib body shape and a code
// describing the rigid body shape
int QuadrilobeCylinder::getNbCorners() const
{
  return ( 1001 );
}





// ----------------------------------------------------------------------------
// Outputs information to be transferred to the fluid
// Same format data as a regular cylinder: center of bottom circular
// face of the elementary cylinder, an arbitrary point on the lateral surface 
// of the elementary cylinder and center of top circular face of the elementary
// cylinder    
void QuadrilobeCylinder::writePositionInFluid( ostream& fluid )
{
  // TO DO
}




// ----------------------------------------------------------------------------
// Reads composite particle data from a stream. Usage: for standard composite
// particles in the 2014 reload format
void QuadrilobeCylinder::read2014( istream& fileIn, vector<Particle*> const*
  	referenceParticles )
{
  CompositeParticle::read2014( fileIn, referenceParticles );

  // We know that (*referenceParticles)[m_GeomType] points to a
  // CompositeParticle, such that we can dynamic cast it to actual type and
  // use -> instead of using get methods through virtual typing
  QuadrilobeCylinder const* CompParticleRef =
  	dynamic_cast<QuadrilobeCylinder const*>(
		(*referenceParticles)[m_GeomType]);

  m_radius = CompParticleRef->m_radius;
  m_armLength = CompParticleRef->m_armLength;
  m_height = CompParticleRef->m_height;
}




// ----------------------------------------------------------------------------
// Reads composite particle data from a stream in a binary form.
// Usage: for standard composite particles in the 2014 reload format
void QuadrilobeCylinder::read2014_binary( istream& fileIn,
	vector<Particle*> const* referenceParticles )
{
  CompositeParticle::read2014_binary( fileIn, referenceParticles );

  // We know that (*referenceParticles)[m_GeomType] points to a
  // CompositeParticle, such that we can dynamic cast it to actual type and
  // use -> instead of using get methods through virtual typing
  QuadrilobeCylinder const* CompParticleRef =
  	dynamic_cast<QuadrilobeCylinder const*>(
		(*referenceParticles)[m_GeomType]);

  m_radius = CompParticleRef->m_radius;
  m_armLength = CompParticleRef->m_armLength;
  m_height = CompParticleRef->m_height;
}




// ----------------------------------------------------------------------------
// Computes and sets the circumscribed radius
void QuadrilobeCylinder::setCircumscribedRadius()
{
  m_geoRBWC->setCircumscribedRadius( sqrt( 0.25 * m_height * m_height 
  	+ ( m_armLength + m_radius ) * ( m_armLength + m_radius ) ) );
}




// ----------------------------------------------------------------------------
// Saves additional features of a (in practice reference) composite particle
// for reload
void QuadrilobeCylinder::writeAdditionalFeatures( ostream& fileSave ) const
{
  fileSave << endl << "*RadiusArmLengthHeight " << m_radius << " " << 
  	m_armLength << " " << m_height;  
  CompositeParticle::writeAdditionalFeatures( fileSave );
}




// ----------------------------------------------------------------------------
// Reads additional features of a (in practice reference) composite particle
// data from a stream
void QuadrilobeCylinder::readAdditionalFeatures( istream& fileIn )
{
  string buffer;
  fileIn >> buffer >> m_radius >> m_armLength >> m_height;
  CompositeParticle::readAdditionalFeatures( fileIn );
}
