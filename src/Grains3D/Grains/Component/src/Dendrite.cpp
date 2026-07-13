#include "GrainsMPIWrapper.hh"
#include "Dendrite.hh"
#include "ContactBuilderFactory.hh"
#include "Memento.hh"
#include "KinematicsBuilderFactory.hh"
#include "GrainsBuilderFactory.hh"
#include "PointC.hh"
#include "RigidBody.hh"
#include "GrainsExec.hh"
#include "ContactForceModel.hh"
#include "Box.hh"
#include <iterator>
#include <algorithm>


// ----------------------------------------------------------------------------
// Constructor with autonumbering as input parameter
Dendrite::Dendrite( bool const& autonumbering )
  : CompositeParticle( autonumbering )
{
  m_specific_composite_shape = "Dendrite";
  m_armWidth = 0;
  m_armLength = 0;
  m_depth = 0;
}




// ----------------------------------------------------------------------------
// Constructor with an XML node as an input parameter. This constructor is
// expected to be used for reference composite particles
Dendrite::Dendrite( DOMNode* root, int const& pc )
  : CompositeParticle( false )
{
  m_specific_composite_shape = "Dendrite";

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

  m_armWidth = ReaderXML::getNodeAttr_Double( nGeometry, "ArmWidth" );
  m_armLength = ReaderXML::getNodeAttr_Double( nGeometry, "ArmLength" );
  m_depth = ReaderXML::getNodeAttr_Double( nGeometry, "Depth" );

  m_mass = m_density * ( ( 6. * m_armWidth * m_armLength * m_depth ) 
  	+ 1.5 * ( m_armWidth * m_armWidth * m_depth ) );
  computeWeight();
  double crust_thickness =
  	ReaderXML::getNodeAttr_Double( nGeometry, "CrustThickness" );
  m_geoRBWC->setCrustThickness( crust_thickness );

  // Material
  DOMNode* nmaterial = ReaderXML::getNode( root, "Material" );
  if ( nmaterial )
  {
    m_materialName = ReaderXML::getNodeValue_String( nmaterial );
    ContactBuilderFactory::defineMaterial( m_materialName, false );
  }

  // Angular position of the composite particle
  m_geoRBWC->getTransform()->load( root );

  // Moment of inertia tensor of the Dendrite
  double global_ixx = 0., global_iyy = 0., global_izz = 0.;
  vector<double> hex_inertia = Dendrite::hex_moment( m_armWidth, m_depth );
  vector<double> straight_arm_inertia = Dendrite::straight_arm_moment(
  	m_armWidth, m_armLength, m_depth );
  vector<double> angled_arm_inertia = Dendrite::angled_arm_moment( m_armWidth, 
  	m_armLength, m_depth );
  global_ixx += hex_inertia[0] + 2. * straight_arm_inertia[0] + 
  	4. * angled_arm_inertia[0];
  global_iyy += hex_inertia[1] + 2. * straight_arm_inertia[1] + 
  	4. * angled_arm_inertia[1];
  global_izz += hex_inertia[2] + 2. * straight_arm_inertia[2] + 
  	4. * angled_arm_inertia[2];
  m_inertia[1] = m_inertia[2] = m_inertia[4] = 0.;
  m_inertia[0] = global_ixx;
  m_inertia[3] = global_iyy;
  m_inertia[5] = global_izz;
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
    m_InitialRelativePositions.push_back( Vector3Null );
    m_InitialRotationMatrices.push_back( ttt );
  }

  // Create the 3 rectangular spines
  double spine_length = 2. * m_armLength + 2. * m_armWidth*sin( M_PI / 3. );
  double ang[] = { PI / 6., PI / 2., 5. * PI / 6. };

  Box* spine_1 = new Box( m_armWidth, spine_length, m_depth );
  Box* spine_2 = new Box( m_armWidth, spine_length, m_depth );
  Box* spine_3 = new Box( m_armWidth, spine_length, m_depth );

  RigidBodyWithCrust* geoRBWC_box1 = new RigidBodyWithCrust( spine_1, 
  	Transform(), false, crust_thickness );
  m_elementaryParticles[0] = new Particle( geoRBWC_box1, m_density, 
  	m_materialName, pc );
  RigidBodyWithCrust* geoRBWC_box2 = new RigidBodyWithCrust( spine_2, 
  	Transform(), false, crust_thickness );
  m_elementaryParticles[1] = new Particle( geoRBWC_box2, m_density,
  	 m_materialName, pc );
  RigidBodyWithCrust* geoRBWC_box3 = new RigidBodyWithCrust( spine_3, 
  	Transform(), false, crust_thickness );
  m_elementaryParticles[2] = new Particle( geoRBWC_box3, m_density, 
  	m_materialName, pc );

  for (size_t m=0; m<m_nbElemPart;++m) 
  {
    m_InitialRelativePositions[m][X] = 0.;
    m_InitialRelativePositions[m][Y] = 0.;
    m_InitialRelativePositions[m][Z] = 0.;
    m_elementaryParticles[m]->setPosition( m_InitialRelativePositions[m] );
    m_elementaryParticles[m]->setMasterParticle( this );

    m_InitialRotationMatrices[m].setValue(
        cos(ang[m]), -sin(ang[m]), 0,
        sin(ang[m]), cos(ang[m]), 0,
        0, 0, 1 );
    m_elementaryParticles[m]->getRigidBody()->getTransform()->setBasis( 
    	m_InitialRotationMatrices[m] );
    m_elementaryParticles[m]->getRigidBody()->composeLeftByTransform( 
    	*(m_geoRBWC->getTransform()) );
  }

  // Set the the circumscribed radius
  CompositeParticle::setCircumscribedRadius();

  // Compute and set the non-spherical bounding volume
  if ( GrainsExec::m_colDetBoundingVolume ) createBoundingVolume();

  // In case part of the particle acceleration computed explicity
  if ( Particle::m_splitExplicitAcceleration ) createVelocityInfosNm1();
}




// ----------------------------------------------------------------------------
// Constructor with input parameters
Dendrite::Dendrite( int const& id_,
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
  m_specific_composite_shape = "Dendrite";
  Dendrite const* DendriteRef =
  	dynamic_cast<Dendrite const*>(ParticleRef);
  m_armWidth = DendriteRef->m_armWidth;
  m_armLength = DendriteRef->m_armLength;
  m_depth = DendriteRef->m_depth;
}




// ----------------------------------------------------------------------------
// Constructor with input parameters
Dendrite::Dendrite( int const& id_,
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
  m_specific_composite_shape = "Dendrite";
  Dendrite const* DendriteRef =
  	dynamic_cast<Dendrite const*>(ParticleRef);
  m_armWidth = DendriteRef->m_armWidth;
  m_armLength = DendriteRef->m_armLength;
  m_depth = DendriteRef->m_depth;
}




// ----------------------------------------------------------------------------
// Destructor
Dendrite::~Dendrite()
{}




// ----------------------------------------------------------------------------
// Copy constructor (the torsor is initialized to 0)
Dendrite::Dendrite( Dendrite const& other,
    	bool const& autonumbering )
  : CompositeParticle( other, autonumbering )
{
  m_armWidth = other.m_armWidth;
  m_armLength = other.m_armLength;
  m_depth = other.m_depth;
}




// ----------------------------------------------------------------------------
// Creates a clone of the particle. This method calls the standard
// copy constructor and is used for new particles to be inserted in the
// simulation. Numbering is automatic, total number of components is
// incremented by 1 and activity is set to WAIT. The calling object is
// expected to be a reference particle
Particle* Dendrite::createCloneCopy( bool const& autonumbering ) const
{
  Particle* particle = new Dendrite( *this, autonumbering );

  return ( particle );
}




// ----------------------------------------------------------------------------
// Creates a clone of the composite particle. This method calls the
// constructor Dendrite( int const& id_, Particle const* ParticleRef,
// Vector3 const& vtrans, Quaternion const& qrot, Vector3 const& vrot,
// Transform const& config, ParticleActivity const& activ ) and is used for
// periodic clone composite particles to be inserted in the simulation.
// Autonumbering is set to false and numbering is set with the parameter id_
Particle* Dendrite::createCloneCopy( int const& id_,
    	Particle const* ParticleRef, Vector3 const& vtrans,
	Quaternion const& qrot,	Vector3 const& vrot,
	Transform const& config, ParticleActivity const& activ,
	map< std::tuple<int,int,int>,
     	std::tuple<bool, Vector3, Vector3, Vector3> > const* contactMap ) const
{
  Particle* particle = new Dendrite( id_, ParticleRef, vtrans,
	qrot, vrot, config, activ, contactMap );

  return ( particle );
}




// ----------------------------------------------------------------------------
// Returns a code describing the rigid body shape
int Dendrite::getShapeCode() const
{
  return ( 1000003 );
}




// ----------------------------------------------------------------------------
// Outputs information to be transferred to the fluid
// Same format data as a regular cylinder: center of bottom circular
// face of the elementary cylinder, an arbitrary point on the lateral surface
// of the elementary cylinder and center of top circular face of the elementary
// cylinder
void Dendrite::writePositionInFluid( ostream& fluid )
{
  fluid << " " << m_armWidth << " " << m_armLength << " " << m_depth << endl;
}




// ----------------------------------------------------------------------------
// Reads composite particle data from a stream. Usage: for standard composite
// particles in the 2014 reload format
void Dendrite::read2014( istream& fileIn, vector<Particle*> const*
  	referenceParticles )
{
  CompositeParticle::read2014( fileIn, referenceParticles );

  // We know that (*referenceParticles)[m_GeomType] points to a
  // CompositeParticle, such that we can dynamic cast it to actual type and
  // use -> instead of using get methods through virtual typing
  Dendrite const* CompParticleRef =
  	dynamic_cast<Dendrite const*>(
		(*referenceParticles)[m_GeomType]);

  m_armWidth = CompParticleRef->m_armWidth;
  m_armLength = CompParticleRef->m_armLength;
  m_depth = CompParticleRef->m_depth;
}




// ----------------------------------------------------------------------------
// Reads composite particle data from a stream in a binary form.
// Usage: for standard composite particles in the 2014 reload format
void Dendrite::read2014_binary( istream& fileIn,
	vector<Particle*> const* referenceParticles )
{
  CompositeParticle::read2014_binary( fileIn, referenceParticles );

  // We know that (*referenceParticles)[m_GeomType] points to a
  // CompositeParticle, such that we can dynamic cast it to actual type and
  // use -> instead of using get methods through virtual typing
  Dendrite const* CompParticleRef =
  	dynamic_cast<Dendrite const*>(
		(*referenceParticles)[m_GeomType]);

  m_armWidth = CompParticleRef->m_armWidth;
  m_armLength = CompParticleRef->m_armLength;
  m_depth = CompParticleRef->m_depth;
}




// ----------------------------------------------------------------------------
// Saves additional features of a (in practice reference) composite particle
// for reload
void Dendrite::writeAdditionalFeatures( ostream& fileSave ) const
{
  fileSave << endl << "*ArmLengthWidthLengthDepth " << m_armWidth << " " << 
  	m_armLength << " " << m_depth;  
  CompositeParticle::writeAdditionalFeatures( fileSave );
}




// ----------------------------------------------------------------------------
// Reads additional features of a (in practice reference) composite particle
// data from a stream
void Dendrite::readAdditionalFeatures( istream& fileIn )
{
  string buffer;
  fileIn >> buffer >> m_armWidth >> m_armLength >> m_depth;
  CompositeParticle::readAdditionalFeatures( fileIn );
}




// ----------------------------------------------------------------------------
// Inertia helper functions
vector<double> Dendrite::hex_moment( double sideLength, double depth ) 
{
  double ixx = (2./3.) * pow(sideLength, 2) * sin(M_PI / 3.) * depth 
  	* (1./4. * pow(depth, 2) + pow(sideLength, 2) * pow(sin(M_PI / 3.), 2));
  ixx += 8./12. * depth * sin(M_PI/3.) * (1./4.*pow(sin(M_PI/3.), 2) 
  	* pow(sideLength, 4) + 1./8. * pow(depth, 2) * pow(sideLength, 2));

  double iyy = 2 * depth * pow(sideLength, 2) * sin(M_PI/3.) * cos(M_PI/3.) 
  	* ((2./3.) * pow(sideLength * cos(M_PI/3.), 2) + (1./6.)*pow(depth,2));
  iyy += 8 * depth * sin(M_PI/3.) * ((7./24.-15./64.) * pow(sideLength, 4) 
  	+ pow(sideLength*depth, 2) * (1./24.-1./32.));

  double izz = ixx + iyy - sin(M_PI/3.)*cos(M_PI/3.) * pow(sideLength,2) 
  	* pow(depth,3);

  return {ixx, iyy, izz};
}

vector<double> Dendrite::straight_arm_moment( double width, double height, 
	double depth ) 
{
  double mass_branch = (width * height * depth); // not including density 
  double dist_com = height / 2. + (width * sin(M_PI / 3.)); // dist to center 
  	// of mass

  double local_ixx = (1./12.) * mass_branch * (pow(height, 2) + pow(depth, 2));
  double local_iyy = (1./12.) * mass_branch * (pow(depth, 2) + pow(width, 2));
  double local_izz = (1./12.) * mass_branch * (pow(height, 2) + pow(width, 2));

  double ixx = local_ixx;
  double iyy = local_iyy;
  double izz = local_izz;

  ixx += mass_branch * pow(dist_com, 2);
  izz += mass_branch * pow(dist_com, 2);

  return {ixx, iyy, izz};
}

vector<double> Dendrite::angled_arm_moment( double width, double height, 
	double depth ) 
{
  double mass_branch = (width * height * depth); // not including density 
  double dist_com = height / 2. + (width * sin(M_PI / 3.)); // dist to center 
  	// of mass

  double ixx = (1./12.) * mass_branch * (pow(height, 2) + pow(depth, 2));
  double iyy = (1./12.) * mass_branch * (pow(depth, 2) + pow(width, 2));
  double izz = (1./12.) * mass_branch * (pow(height, 2) + pow(width, 2));

  double ang = M_PI / 3.;
  double dist_x = dist_com * cos(ang);
  double dist_y = dist_com * sin(ang);

  ixx += mass_branch * pow(dist_com, 2);
  izz += mass_branch * pow(dist_com, 2);

  // return inertia tensor with applied change of basis
  return {iyy * pow(sin(ang), 2) + ixx * pow(cos(ang), 2), 
  	ixx * pow(sin(ang), 2) + iyy * pow(cos(ang), 2), izz};
}
