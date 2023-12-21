#include "GrainsMPIWrapper.hh"
#include "Particle.hh"
#include "KinematicsBuilderFactory.hh"
#include "ContactBuilderFactory.hh"
#include "GrainsBuilderFactory.hh"
#include "Obstacle.hh"
#include "SimpleObstacle.hh"
#include "Memento.hh"
#include "Sphere.hh"
#include "Disc.hh"
#include "LinkedCell.hh"
#include "GrainsExec.hh"
#include "ContactForceModel.hh"
#include <algorithm>
#include <sstream>
#include <string>
using namespace std;


// Initialisation des attributs static
double Particle::m_fluidDensity = 0.;
double Particle::m_fluidViscosity = 0.;
bool Particle::m_fluidCorrectedAcceleration = true ;
bool Particle::m_splitExplicitAcceleration = false;


// ----------------------------------------------------------------------------
// Constructor with autonumbering as input parameter
Particle::Particle( bool const& autonumbering )
  : Component( autonumbering )
  , m_masterParticle( this )
  , m_kinematics( NULL )
  , m_density( 2500. )
  , m_activity( WAIT )
  , m_VelocityInfosNm1( NULL )
  , m_tag( 0 )
  , m_GeoLoc( GEOPOS_NONE )
  , m_cellule( NULL )
  , m_tag_nm1( 0 )
  , m_GeoLoc_nm1( GEOPOS_NONE )
  , m_cellule_nm1( NULL )
  , m_GeomType( 0 )
  , m_coordination_number( 0 )
  , m_specific_composite_shape( "none" )
{
  // Initialize inertia
  for (int i=0; i<6; i++)
  {
    m_inertia[i] = 0.0;
    m_inertia_1[i] = 0.0;
  }
}




// ----------------------------------------------------------------------------
// Copy constructor (the torsor is initialized to 0)
Particle::Particle( Particle const& other, bool const& autonumbering )
  : Component( other, autonumbering )
  , m_masterParticle( this )
  , m_density( other.m_density )
  , m_activity( WAIT )
  , m_VelocityInfosNm1( NULL )
  , m_tag( other.m_tag )
  , m_GeoLoc( other.m_GeoLoc )
  , m_cellule( other.m_cellule )
  , m_tag_nm1( other.m_tag )
  , m_GeoLoc_nm1( other.m_GeoLoc )
  , m_cellule_nm1( other.m_cellule_nm1 )
  , m_GeomType( other.m_GeomType )
  , m_coordination_number( other.m_coordination_number )
  , m_weight( other.m_weight )
  , m_specific_composite_shape( other.m_specific_composite_shape )  
{
  m_kinematics = other.m_kinematics->clone();
  copy( &other.m_inertia[0], &other.m_inertia[6], &m_inertia[0] );
  copy( &other.m_inertia_1[0], &other.m_inertia_1[6], &m_inertia_1[0] );

  if ( other.m_VelocityInfosNm1 )
  {
    m_VelocityInfosNm1 = new struct VelocityInfosNm1;
    m_VelocityInfosNm1->TranslationalVelocity_nm1 = other.m_VelocityInfosNm1->
    	TranslationalVelocity_nm1;
    m_VelocityInfosNm1->TranslationalVelocity_difference =
    	other.m_VelocityInfosNm1->TranslationalVelocity_difference;
    m_VelocityInfosNm1->RotationalVelocity_nm1 = other.m_VelocityInfosNm1->
    	RotationalVelocity_nm1;
    m_VelocityInfosNm1->RotationalVelocity_difference =
    	other.m_VelocityInfosNm1->RotationalVelocity_difference;
  }
}




// ----------------------------------------------------------------------------
// Constructor with an XML node as an input parameter. This constructor is
// expected to be used for reference particles
Particle::Particle( DOMNode* root, int const& pc )
  : Component( false )
  , m_masterParticle( this )
  , m_kinematics( NULL )
  , m_density( 2500. )
  , m_activity( WAIT )
  , m_VelocityInfosNm1( NULL )
  , m_tag( 0 )
  , m_GeoLoc( GEOPOS_NONE )
  , m_cellule( NULL )
  , m_tag_nm1( 0 )
  , m_GeoLoc_nm1( GEOPOS_NONE )
  , m_cellule_nm1( NULL )
  , m_GeomType( pc )
  , m_coordination_number( 0 )
  , m_specific_composite_shape( "none" )
{
  for (int i=0; i<6; i++)
  {
    m_inertia[i] = 0.0;
    m_inertia_1[i] = 0.0;
  }

  m_geoRBWC = new RigidBodyWithCrust( root );
  m_kinematics = KinematicsBuilderFactory::create(
  	m_geoRBWC->getConvex() );

  // Particle density
  if ( ReaderXML::hasNodeAttr( root, "Density" ) )
    m_density = ReaderXML::getNodeAttr_Double( root, "Density" );

  // Material
  DOMNode* material_ = ReaderXML::getNode( root, "Material" );
  if ( material_ )
  {
    m_materialName = ReaderXML::getNodeValue_String( material_ );
    ContactBuilderFactory::defineMaterial( m_materialName, false );
  }

  // Mass and inertia
  m_mass = m_density * m_geoRBWC->getVolume();
  m_geoRBWC->BuildInertia( m_inertia, m_inertia_1 );
  for (int i=0; i<6; i++)
  {
    m_inertia[i] *= m_density;
    m_inertia_1[i] /= m_density;
  }

  // Weight
  computeWeight();

  // In case part of the particle acceleration is computed explicity
  if ( Particle::m_splitExplicitAcceleration ) createVelocityInfosNm1();
}




// ----------------------------------------------------------------------------
// Constructor with input parameters. This constructor is
// expected to be used for reference particles
Particle::Particle( RigidBodyWithCrust* georbwc, double const& density,
    	string const& mat, int const& pc )
  : Component( false )
  , m_masterParticle( this )
  , m_kinematics( NULL )
  , m_density( 2500. )
  , m_activity( WAIT )
  , m_VelocityInfosNm1( NULL )
  , m_tag( 0 )
  , m_GeoLoc( GEOPOS_NONE )
  , m_cellule( NULL )
  , m_tag_nm1( 0 )
  , m_GeoLoc_nm1( GEOPOS_NONE )
  , m_cellule_nm1( NULL )
  , m_GeomType( pc )
  , m_coordination_number( 0 )
  , m_specific_composite_shape( "none" )
{
  for (int i=0; i<6; i++)
  {
    m_inertia[i] = 0.0;
    m_inertia_1[i] = 0.0;
  }

  m_geoRBWC = georbwc;
  m_kinematics = KinematicsBuilderFactory::create(
  	m_geoRBWC->getConvex() );

  // Particle density
  m_density = density;

  // Material
  m_materialName = mat;
  ContactBuilderFactory::defineMaterial( m_materialName, false );

  // Mass and inertia
  m_mass = m_density * m_geoRBWC->getVolume();
  m_geoRBWC->BuildInertia( m_inertia, m_inertia_1 );
  for (int i=0; i<6; i++)
  {
    m_inertia[i] *= m_density;
    m_inertia_1[i] /= m_density;
  }

  // Weight
  computeWeight();

  // In case part of the particle acceleration is computed explicity
  if ( Particle::m_splitExplicitAcceleration ) createVelocityInfosNm1();
}




// ----------------------------------------------------------------------------
// Constructor with input parameters
Particle::Particle( int const& id_, Particle const* ParticleRef,
	double const& vx, double const& vy, double const& vz,
	double const& qrotationx, double const& qrotationy,
	double const& qrotationz, double const& qrotations,
	double const& rx, double const& ry, double const& rz,
	const double m[12],
	ParticleActivity const& activ,
	int const& tag_,
	int const& coordination_number_ )
  : Component( false )
  , m_masterParticle( this )
  , m_kinematics( NULL )
  , m_activity( activ )
  , m_VelocityInfosNm1( NULL )
  , m_tag( tag_ )
  , m_GeoLoc( GEOPOS_NONE )
  , m_cellule_nm1( NULL )
  , m_coordination_number( coordination_number_ )
  , m_specific_composite_shape( "none" )
{
  // Initialize inertia
  for (int i=0; i<6; i++)
  {
    m_inertia[i] = 0.0;
    m_inertia_1[i] = 0.0;
  }

  // ID number
  m_id = id_;

  // Particle class
  m_GeomType = ParticleRef->m_GeomType;

  // Rigid body shape
  m_geoRBWC = new RigidBodyWithCrust( *ParticleRef->m_geoRBWC );

  // Transform
  m_geoRBWC->setTransform( m );

  // Kinematics
  m_kinematics = ParticleRef->m_kinematics->clone();
  m_kinematics->setTranslationalVelocity( Vector3( vx, vy, vz ) );
  m_kinematics->setQuaternionRotation( qrotationx, qrotationy,
	qrotationz, qrotations );
  m_kinematics->setAngularVelocity( Vector3( rx, ry, rz ) );

  // Material
  m_materialName = ParticleRef->m_materialName;

  // Mass and inertia
  m_density = ParticleRef->m_density;
  m_mass = ParticleRef->m_mass;
  for (int i=0; i<6; i++)
  {
    m_inertia[i] = ParticleRef->m_inertia[i];
    m_inertia_1[i] = ParticleRef->m_inertia_1[i];
  }

  // Weight
  computeWeight();

  // In case part of the particle acceleration is computed explicity
  if ( Particle::m_splitExplicitAcceleration ) createVelocityInfosNm1();
}




// ----------------------------------------------------------------------------
// Constructor with input parameters
Particle::Particle( int const& id_, Particle const* ParticleRef,
	Vector3 const& vtrans, Quaternion const& qrot,
	Vector3 const& vrot, Transform const& config,
	ParticleActivity const& activ,
     	map< std::tuple<int,int,int>,
     	std::tuple<bool, Vector3, Vector3, Vector3> > const* contactMap )
  : Component( false )
  , m_masterParticle( this )
  , m_kinematics( NULL )
  , m_activity( activ )
  , m_VelocityInfosNm1( NULL )
  , m_tag( 0 )
  , m_GeoLoc( GEOPOS_NONE )
  , m_cellule( NULL )
  , m_tag_nm1( 0 )
  , m_GeoLoc_nm1( GEOPOS_NONE )
  , m_cellule_nm1( NULL )
  , m_coordination_number( 0 )
  , m_specific_composite_shape( "none" )  
{
  // ID number
  m_id = id_;

  // Particle class
  m_GeomType = ParticleRef->m_GeomType;

  // Rigid body shape
  m_geoRBWC = new RigidBodyWithCrust( *(ParticleRef->m_geoRBWC) );

  // Transform
  m_geoRBWC->setTransform( config );

  // Kinematics
  m_kinematics = ParticleRef->m_kinematics->clone();
  m_kinematics->setTranslationalVelocity( vtrans );
  m_kinematics->setQuaternionRotation( qrot );
  m_kinematics->setAngularVelocity( vrot );

  // Material
  m_materialName = ParticleRef->m_materialName;

  // Mass and inertia
  m_density = ParticleRef->m_density;
  m_mass = ParticleRef->m_mass;
  for (int i=0; i<6; i++)
  {
    m_inertia[i] = ParticleRef->m_inertia[i];
    m_inertia_1[i] = ParticleRef->m_inertia_1[i];
  }
  
  // Contact map
  if ( contactMap->size() ) m_contactMap = *contactMap;   

  // Weight
  computeWeight();

  // In case part of the particle acceleration is computed explicity
  if ( Particle::m_splitExplicitAcceleration ) createVelocityInfosNm1();
}




// ----------------------------------------------------------------------------
// Creates a clone of the particle. This method calls the standard
// copy constructor and is used for new particles to be inserted in the
// simulation. Activity is set to WAIT. The calling object is
// expected to be a reference particle
Particle* Particle::createCloneCopy( bool const& autonumbering ) const
{
  Particle* particle = new Particle( *this, autonumbering );

  return ( particle );
}




// ----------------------------------------------------------------------------
// Creates a clone of the particle. This method calls the
// constructor Particle( int const& id_, Particle const* ParticleRef, Vector3
// const& vtrans, Quaternion const& qrot, Vector3 const& vrot,	Transform
// const& config, ParticleActivity const& activ ) and is used for periodic
// clone particles to be inserted in the simulation. Autonumbering
// is set to false and numbering is set with the parameter id_
Particle* Particle::createCloneCopy( int const& id_,
    	Particle const* ParticleRef, Vector3 const& vtrans,
	Quaternion const& qrot,	Vector3 const& vrot,
	Transform const& config, ParticleActivity const& activ,
     	map< std::tuple<int,int,int>,
     	std::tuple<bool, Vector3, Vector3, Vector3> > const* contactMap ) const
{
  Particle* particle = new Particle( id_, ParticleRef, vtrans,
	qrot, vrot, config, activ, contactMap );

  return ( particle );
}




// ----------------------------------------------------------------------------
// Computes particle net weight
void Particle::computeWeight()
{
  m_weight = m_mass
  	* ( 1. - Particle::m_fluidDensity / m_density )
  	* GrainsExec::m_vgravity ;
}




// ----------------------------------------------------------------------------
// Destructor
Particle::~Particle()
{
  delete m_kinematics;
  if ( m_VelocityInfosNm1 ) delete m_VelocityInfosNm1;
}




// ----------------------------------------------------------------------------
// Solves the Newton's law and move particle to their new position
void Particle::Move( double time, 
	double const& dt_particle_vel, 
    	double const& dt_particle_disp )
{
  try {
  // Time integration of Newton's law and kinematic equations
  double depl = m_kinematics->Move( this, dt_particle_vel, 
    	dt_particle_disp );

  // Check that translational displacement is smaller than crust thickness
  double crust = m_geoRBWC->getCrustThickness();
  if ( depl > crust )
  {
    cout << endl << "Processor = " <<
    	(GrainsExec::m_MPI ? GrainsExec::getComm()->get_rank_active() : 0 )
	<< " has thrown an DisplacementError exception" <<  endl;
    GrainsExec::m_exception_Displacement = true;
    DisplacementError erreur( this, depl, crust, time );
    throw erreur;
  }

  }
  catch (const DisplacementError&) {
    throw DisplacementError();
  }
}




// ----------------------------------------------------------------------------
// Computes acceleration
void Particle::computeAcceleration( double time )
{
  m_kinematics->computeAcceleration( m_torsor, this );
}




// ----------------------------------------------------------------------------
// Advances velocity over dt_particle_vel
void Particle::advanceVelocity( double time, double const& dt_particle_vel )
{
  m_kinematics->advanceVelocity( dt_particle_vel );
}




// ----------------------------------------------------------------------------
// Contact between a particle and a component. If contact exists,
// computes the contact force and torque and adds to each component
void Particle::InterAction( Component* voisin,
	double dt, double const& time, LinkedCell* LC )
{
  try {
  PointContact closestPoint;
  double delta=0.;

  try {
    closestPoint = m_geoRBWC->ClosestPoint( *(voisin->getRigidBody()) );
  }
  catch ( ContactError &erreur )
  {
    try {
      closestPoint = m_geoRBWC->ClosestPoint_ErreurHandling(
      	*(voisin->getRigidBody()), 10., m_id, voisin->getID() );
    }
    catch (ContactError &erreur_level2)
    {
      cout << endl << "Processor = "
	<< ( GrainsExec::m_MPI ? GrainsExec::getComm()->get_rank_active() : 0 )
	<< " has thrown an ContactError exception" <<  endl;
      erreur_level2.setMessage( "Particle::InterAction : choc de croute a t="
      	+ GrainsExec::doubleToString( time, TIMEFORMAT ) );
      erreur_level2.setComponents( this, voisin, time );
      GrainsExec::m_exception_Contact = true;
      throw(erreur_level2);
    }
  }

  LC->addToContactsFeatures( time, closestPoint );
  delta = closestPoint.getOverlapDistance();

  if( delta < 0. )
  {
    if ( ContactBuilderFactory::contactForceModel(
		m_materialName, voisin->getMaterial() )
    		->computeForces( this, voisin, closestPoint, LC, dt ) )
    {
      // Note: in this method this and voisin cannot be a CompositeParticle
      // thus there is no need to call getMasterComponent() before
      // addToCoordinationNumber( 1 )
      this->addToCoordinationNumber( 1 );
      voisin->addToCoordinationNumber( 1 );
    }
  }

  }
  catch (const ContactError&) {
    throw ContactError();
  }
}




// ----------------------------------------------------------------------------
// Searches and stores all contact points between a composite particle and a
// component.
void Particle::SearchContact( Component* voisin, double dt,
      double const& time, LinkedCell *LC, list<ContactInfos*>& listContact )
{
  try{
  PointContact closestPoint;

  try {
    closestPoint = m_geoRBWC->ClosestPoint( *(voisin->getRigidBody()) );
  }
  catch ( ContactError &erreur )
  {
    try {
      closestPoint = m_geoRBWC->ClosestPoint_ErreurHandling(
      	*(voisin->getRigidBody()), 10., m_id, voisin->getID() );
    }
    catch (ContactError &erreur_level2)
    {
      cout << endl << "Processor = "
	<< ( GrainsExec::m_MPI ? GrainsExec::getComm()->get_rank_active() : 0 )
	<< " has thrown an ContactError exception" <<  endl;
      erreur_level2.setMessage( "Particle::InterAction : choc de croute a t="
      	+ GrainsExec::doubleToString( time, TIMEFORMAT ) );
      erreur_level2.setComponents( this, voisin, time );
      GrainsExec::m_exception_Contact = true;
      throw(erreur_level2);
    }
  }

  ContactInfos* result = NULL ;
  if ( closestPoint.getOverlapDistance() < 0. )
  {
    result = new struct ContactInfos;
    result->ContactPoint = closestPoint;
    result->p0 = this;
    result->p1 = voisin;
    listContact.push_back( result );
  }

  }
  catch (const ContactError&) {
    throw ContactError();
  }
}




// ----------------------------------------------------------------------------
// Returns a pointer to the particle kinematics
ParticleKinematics const* Particle::getKinematics() const
{
  return ( m_kinematics );
}




// ----------------------------------------------------------------------------
// Returns the fluid density
double Particle::getFluidDensity()
{
  return ( Particle::m_fluidDensity );
}




// ----------------------------------------------------------------------------
// Returns the viscosity of the surrounding fluid
double Particle::getFluidViscosity()
{
  return Particle::m_fluidViscosity ;
}




// ----------------------------------------------------------------------------
// Returns particle inertia tensor
double const* Particle::getInertiaTensorBodyFixed() const
{
  return ( m_inertia );
}




// ----------------------------------------------------------------------------
// Inertie inverse de la particle.
double const* Particle::getInverseInertiaTensorBodyFixed() const
{
  return ( m_inertia_1 );
}




// ----------------------------------------------------------------------------
// Returns particle density
double Particle::getDensity() const
{
  return ( m_density );
}




// ----------------------------------------------------------------------------
// Returns the radius of the sphere of same volume
double Particle::getEquivalentSphereRadius() const
{
  return ( pow( ( 0.75 / PI ) * m_mass / m_density, 1. / 3. ) );
}




// ----------------------------------------------------------------------------
// Returns the velocity at a point in space based on the
// translational and angular velocity of the particle. This method assumes
// that the point belongs to the particle but this assumption is not verified.
Vector3 Particle::getVelocityAtPoint( Point3 const& pt ) const
{
  Vector3 levier = pt - *m_geoRBWC->getCentre();
  return ( m_kinematics->Velocity( levier ) );
}




// ----------------------------------------------------------------------------
// Returns angular velocity
Vector3 const* Particle::getAngularVelocity() const
{
  return ( m_kinematics->getAngularVelocity() );
}




// ----------------------------------------------------------------------------
// Returns translational velocity
Vector3 const* Particle::getTranslationalVelocity() const
{
  return ( m_kinematics->getTranslationalVelocity() );
}




// ----------------------------------------------------------------------------
// Returns total force exerted on the particle
Vector3 const* Particle::getForce() const
{
  return ( m_torsor.getForce() );
}




// ----------------------------------------------------------------------------
// Returns the rotation quaternion
Quaternion const* Particle::getQuaternionRotation() const
{
  return ( m_kinematics->getQuaternionRotation() );
}




// ----------------------------------------------------------------------------
// Copies rotation quaternion in a 1D array
void Particle::copyQuaternionRotation( double* vit, int i ) const
{
  Quaternion const* qr = m_kinematics->getQuaternionRotation();
  Vector3 const* vqr = qr->getVector3();
  for (int j=0 ;j<3; j++) vit[i+j] = (*vqr)[j];
  vit[i+3] = qr->getdouble();
}




// ----------------------------------------------------------------------------
// Copies force and torque in a 1D array
void Particle::copyForceTorque( double* fm, int i ) const
{
  m_torsor.copyForceTorque( fm, i );
}




// ----------------------------------------------------------------------------
// Copies angular velocity in a 1D array
void Particle::copyAngularVelocity( double* vit, int i ) const
{
  Vector3 const* vr = m_kinematics->getAngularVelocity();
  for (int j=0 ;j<3; j++) vit[i+j] = (*vr)[j];
}




// ----------------------------------------------------------------------------
// Copies translational velocity in a 1D array
void Particle::copyTranslationalVelocity( double* vit, int i ) const
{
  Vector3 const* vt = m_kinematics->getTranslationalVelocity();
  for (int j=0 ;j<3; j++) vit[i+j] = (*vt)[j];
}




// ----------------------------------------------------------------------------
// Copies kinematics at time t-2dt (translational velocity, angular
// velocity, variation of translational velocity, variation of angular
// velocity) in a 1D array
void Particle::copyKinematicsNm2( double* vit, int i ) const
{
  m_kinematics->copyKinematicsNm2( vit, i );
}




// ----------------------------------------------------------------------------
// Adds a force whose point of application is different from the
// reference point of the torsor (additional torque contribution)
void Particle::addForce( Point3 const& point, Vector3 const& f_ )
{
  m_torsor.addForce( point, f_ );
}




// ----------------------------------------------------------------------------
// Resets kinematics and transformation to 0
void Particle::reset()
{
  Point3 zero;
  m_geoRBWC->setOrigin( (double *) &zero );
  m_kinematics->reset();
}




// ----------------------------------------------------------------------------
// Resets kinematics to 0
void Particle::resetKinematics()
{
  m_kinematics->reset();
}




// ----------------------------------------------------------------------------
// Saves particle state
void Particle::saveState()
{
  if (!m_memento)
    m_memento = new ConfigurationMemento();
  m_memento->m_position = *m_geoRBWC->getTransform();
  m_kinematics->saveState();
}




// ----------------------------------------------------------------------------
// Creates and returns particle state
pair<ConfigurationMemento*,ParticleKinematicsMemento*> Particle::createState()
{
  ConfigurationMemento* Pmemento_ = new ConfigurationMemento();
  Pmemento_->m_position = *m_geoRBWC->getTransform();
  ParticleKinematicsMemento* Cmemento_ = m_kinematics->createState();
  pair<ConfigurationMemento*,ParticleKinematicsMemento*> ppp( Pmemento_,
  	Cmemento_ );

  return ( ppp );
}




// ----------------------------------------------------------------------------
// Restores particle state
void Particle::restoreState()
{
  m_geoRBWC->setTransform( m_memento->m_position );
  m_kinematics->restoreState();
}




// ----------------------------------------------------------------------------
// Sets the fluid density
void Particle::setFluidDensity( double rho )
{
  m_fluidDensity = rho;
}




// ----------------------------------------------------------------------------
// Sets the viscosity of the surrounding fluid
void Particle::setFluidViscosity( double mu )
{
  m_fluidViscosity = mu;
}




// ----------------------------------------------------------------------------
// Sets the angular velocity
void Particle::setAngularVelocity( Vector3 const& vrot )
{
  m_kinematics->setAngularVelocity( vrot );
}




// ----------------------------------------------------------------------------
// Sets the translation velocity
void Particle::setTranslationalVelocity( Vector3 const& vtrans )
{
  m_kinematics->setTranslationalVelocity( vtrans );
}




// ----------------------------------------------------------------------------
// Sets kinematics at time t-2dt as a 1D array of 12 scalars
void Particle::setKinematicsNm2( double const* tab )
{
  m_kinematics->setKinematicsNm2( tab );
}




// ----------------------------------------------------------------------------
// Reads a (in practice reference) particle data from a stream
void Particle::read( istream& fileIn, bool elemPart )
{
  string buffer;

  // ID number
  fileIn >> buffer >> m_id;

  // Material name
  fileIn >> buffer >> m_materialName;

  // Construction of the rigid body with crust
  // The flow is RigidBodyWithCrust( fileIn )
  // => ConvexBuilderFactory::create( cle, fileIn )
  // => xxx::create( fileIn )
  // => constructeur de xxx
  m_geoRBWC = new RigidBodyWithCrust( fileIn );

  if ( !elemPart )
  {
    // Geometric type of particle
    fileIn >> buffer >> m_GeomType;

    // Particle tag
    fileIn >> buffer >> m_tag;

    // Particle mass
    fileIn >> buffer >> m_mass ;

    // Particle density
    fileIn >> buffer >> m_density ;

    // Moment of inertia tensor and its inverse
    fileIn >> buffer
	   >> m_inertia[0] >> m_inertia[1] >> m_inertia[2]
	   >> m_inertia[3] >> m_inertia[4] >> m_inertia[5];
    fileIn >> buffer
	   >> m_inertia_1[0] >> m_inertia_1[1] >> m_inertia_1[2]
	   >> m_inertia_1[3] >> m_inertia_1[4] >> m_inertia_1[5];

    // Particle position
    fileIn >> buffer;
    m_geoRBWC->readPosition( fileIn );

    // Particle activity
    bool actif;
    fileIn >> buffer >> actif;
    m_activity = ( actif == true ) ? COMPUTE : WAIT;

    // Particle kinematics
    if ( m_kinematics ) delete m_kinematics;
    m_kinematics = KinematicsBuilderFactory::read( fileIn,
  	m_geoRBWC->getConvex() );

    // Read additional features
    readAdditionalFeatures( fileIn );

    // In case part of the particle acceleration is computed explicity
    if ( Particle::m_splitExplicitAcceleration ) createVelocityInfosNm1();

    // Compute particle weight
    computeWeight();
  }
}




// ----------------------------------------------------------------------------
// Reads additional features of a (in practice reference) particle data from
// a stream
void Particle::readAdditionalFeatures( istream& fileIn )
{}




// ----------------------------------------------------------------------------
// Reads particle data from a stream. Usage: for standard particles
// in the 2014 reload format
void Particle::read2014( istream& fileIn, vector<Particle*> const*
  	referenceParticles )
{
  // ID number
  // Note: the geometric type of particle m_GeomType is read and set prior
  // to constructing the particle
  fileIn >> m_id;

  // Create the rigid body using the copy constructor of RigidBodyWithCrust
  m_geoRBWC = new RigidBodyWithCrust(
    *(*referenceParticles)[m_GeomType]->getRigidBody() );

  // Material name, mass and density
  m_materialName =
  	(*referenceParticles)[m_GeomType]->m_materialName ;
  m_mass = (*referenceParticles)[m_GeomType]->m_mass ;
  m_density = (*referenceParticles)[m_GeomType]->m_density ;

  // Moment of inertia tensor and its inverse
  for (size_t i=0;i<6;++i)
  {
    m_inertia[i] =
    	(*referenceParticles)[m_GeomType]->m_inertia[i] ;
    m_inertia_1[i] =
    	(*referenceParticles)[m_GeomType]->m_inertia_1[i] ;
  }

  // Particle tag
  fileIn >> m_tag;

  // Particle activity
  bool actif;
  fileIn >> actif;
  m_activity = ( actif == true ) ? COMPUTE : WAIT;

  // Particle position
  m_geoRBWC->readPosition2014( fileIn );

  // Particle kinematics
  if ( m_kinematics ) delete m_kinematics;
  m_kinematics = KinematicsBuilderFactory::create(
  	m_geoRBWC->getConvex() );
  m_kinematics->readParticleKinematics2014( fileIn );
  
  // Read contact map
  readContactMap2014( fileIn );

  // In case part of the particle acceleration is computed explicity
  if ( Particle::m_splitExplicitAcceleration ) createVelocityInfosNm1();

  // Compute particle weight
  computeWeight();
}




// ----------------------------------------------------------------------------
// Reads particle data from a stream in a binary form. Usage: for
// standard particles in the 2014 reload format
void Particle::read2014_binary( istream& fileIn, vector<Particle*> const*
  	referenceParticles )
{
  // ID number
  // Note: the geometric type of particle m_GeomType is read and set prior
  // to constructing the particle
  // Note: the m_GeomType is read and set prior to constructing the particle
  fileIn.read( reinterpret_cast<char*>( &m_id ), sizeof(int) );

  // Create the rigid body using the copy constructor of RigidBodyWithCrust
  m_geoRBWC = new RigidBodyWithCrust(
    *(*referenceParticles)[m_GeomType]->getRigidBody() );

  // Material name, mass and density
  m_materialName =
  	(*referenceParticles)[m_GeomType]->m_materialName ;
  m_mass = (*referenceParticles)[m_GeomType]->m_mass ;
  m_density = (*referenceParticles)[m_GeomType]->m_density ;

  // Moment of inertia tensor and its inverse
  for (size_t i=0;i<6;++i)
  {
    m_inertia[i] =
    	(*referenceParticles)[m_GeomType]->m_inertia[i] ;
    m_inertia_1[i] =
    	(*referenceParticles)[m_GeomType]->m_inertia_1[i] ;
  }

  // Particle tag
  fileIn.read( reinterpret_cast<char*>( &m_tag ), sizeof(int) );

  // Particle activity
  unsigned int iact = 0;
  fileIn.read( reinterpret_cast<char*>( &iact ), sizeof(unsigned int) );
  if ( iact == 1 ) m_activity = COMPUTE;
  else m_activity = WAIT;

  // Particle position
  m_geoRBWC->readPosition2014_binary( fileIn );

  // Particle kinematics
  if ( m_kinematics ) delete m_kinematics;
  m_kinematics = KinematicsBuilderFactory::create(
  	m_geoRBWC->getConvex() );
  m_kinematics->readParticleKinematics2014_binary( fileIn );

  // Read contact map
  readContactMap2014_binary( fileIn );

  // In case part of the particle acceleration is computed explicity
  if ( Particle::m_splitExplicitAcceleration ) createVelocityInfosNm1();

  // Compute particle weight
  computeWeight();
}




// ----------------------------------------------------------------------------
// Saves a (in practice reference) particle for reload
void Particle::write( ostream& fileSave ) const
{
  fileSave << endl << ( isCompositeParticle() ? "<CompositeParticle>" :
  	"<Particle>" ) << endl;
  if ( isCompositeParticle() ) fileSave << "*SpecificShape " 
  	<< m_specific_composite_shape << endl;
  writeStatic( fileSave );
  fileSave << "*Type " << m_GeomType << endl;
  fileSave << "*Tag " << m_tag << endl;
  fileSave << "*Mass "
	   << GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
  		m_mass ) << endl;
  fileSave << "*Density "
	   << GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
  		m_density ) << endl;
  fileSave << "*MomentOfInertiaTensor" << endl
	<< GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
  		m_inertia[0] ) << " "
	<< GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
  		m_inertia[1] ) << " "
	<< GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
  		m_inertia[2] ) << " "
	<< GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
  		m_inertia[3] ) << " "
	<< GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
  		m_inertia[4] ) << " "
	<< GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
  		m_inertia[5] ) << endl;
  fileSave << "*InverseOfMOIT" << endl
	<< GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
  		m_inertia_1[0] ) << " "
	<< GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
  		m_inertia_1[1] ) << " "
	<< GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
  		m_inertia_1[2] ) << " "
	<< GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
  		m_inertia_1[3] ) << " "
	<< GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
  		m_inertia_1[4] ) << " "
	<< GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
  		m_inertia_1[5] ) << endl;

  Component::writePosition( fileSave );
  fileSave << endl;
  bool b_activ = ( m_activity == COMPUTE ) ? true : false;
  fileSave << "*Activity " << b_activ << endl;
  m_kinematics->writeParticleKinematics( fileSave );
  writeAdditionalFeatures( fileSave );
  fileSave << endl << ( isCompositeParticle() ? "</CompositeParticle>" :
  	"</Particle>" ) << endl;
}




// ----------------------------------------------------------------------------
// Saves additional features of a (in practice reference) particle for reload
void Particle::writeAdditionalFeatures( ostream& fileSave ) const
{}




// ----------------------------------------------------------------------------
// Writes the particle's "static" data
void Particle::writeStatic( ostream& fileOut ) const
{
  Component::writeStatic( fileOut );
}




// ----------------------------------------------------------------------------
// Saves particle for reload in 2014 format
void Particle::write2014( ostream& fileSave ) const
{
  fileSave << m_GeomType << " " << m_id << " " << m_tag << " " <<
  	( m_activity == COMPUTE ? true : false ) << " ";
  m_geoRBWC->getTransform()->writeTransform2014( fileSave );
  fileSave << " ";
  m_kinematics->writeParticleKinematics2014( fileSave );
  writeContactMemory2014( fileSave );
  fileSave << endl;
}




// ----------------------------------------------------------------------------
// Saves particle for reload in 2014 and binary format
void Particle::write2014_binary( ostream& fileSave )
{
  fileSave.write( reinterpret_cast<char*>( &m_GeomType ), sizeof(int) );
  fileSave.write( reinterpret_cast<char*>( &m_id ), sizeof(int) );
  fileSave.write( reinterpret_cast<char*>( &m_tag ), sizeof(int) );
  unsigned int iact = m_activity == COMPUTE ? true : false ;
  fileSave.write( reinterpret_cast<char*>( &iact ), sizeof(unsigned int) );
  m_geoRBWC->getTransform()->writeTransform2014_binary( fileSave );
  m_kinematics->writeParticleKinematics2014_binary( fileSave );
  writeContactMemory2014_binary( fileSave );
}




// ----------------------------------------------------------------------------
// Saves the identity of the particle (often id number and address)
void Particle::writeIdentity( ostream& file ) const
{
  file << m_id;
}




// ----------------------------------------------------------------------------
// Returns translational velocity difference at the previous discrete time
Vector3 Particle::getTranslationalVelocityDifferencePreviousTime() const
{
  return ( m_VelocityInfosNm1->TranslationalVelocity_difference );
}




// ----------------------------------------------------------------------------
// Returns angular velocity difference at the previous discrete time
Vector3 Particle::getRotationalVelocityDifferencePreviousTime() const
{
  return ( m_VelocityInfosNm1->RotationalVelocity_difference );
}




// ----------------------------------------------------------------------------
// Updates velocity difference and velocity at previous discrete time
void Particle::setVelocityAndVelocityDifferencePreviousTime()
{
  m_VelocityInfosNm1->TranslationalVelocity_difference =
  	*m_kinematics->getTranslationalVelocity()
  	- m_VelocityInfosNm1->TranslationalVelocity_nm1;
  m_VelocityInfosNm1->TranslationalVelocity_nm1 =
  	*m_kinematics->getTranslationalVelocity();
  m_VelocityInfosNm1->RotationalVelocity_difference =
  	*m_kinematics->getAngularVelocity()
  	- m_VelocityInfosNm1->RotationalVelocity_nm1;
  m_VelocityInfosNm1->RotationalVelocity_nm1 =
  	*m_kinematics->getAngularVelocity();
}




// ----------------------------------------------------------------------------
// Sets the velocity at the previous discrete time when a simulation
// is restarted (as Grains3D does not save this information)
void Particle::setVelocityPreviousTimeRestart(
  	double const& vx, double const& vy, double const& vz,
  	double const& omx, double const& omy, double const& omz )
{
  m_VelocityInfosNm1->TranslationalVelocity_nm1[X] = vx ;
  m_VelocityInfosNm1->TranslationalVelocity_nm1[Y] = vy ;
  m_VelocityInfosNm1->TranslationalVelocity_nm1[Z] = vz ;
  m_VelocityInfosNm1->RotationalVelocity_nm1[X] = omx ;
  m_VelocityInfosNm1->RotationalVelocity_nm1[Y] = omy ;
  m_VelocityInfosNm1->RotationalVelocity_nm1[Z] = omz ;
}




// ----------------------------------------------------------------------------
// Defines whether the particle acceleration (i.e. change of
// momentum) is corrected by the fluid density in case of immersed rigid
// bodies
void Particle::setFluidCorrectedAcceleration( bool correct )
{
  Particle::m_fluidCorrectedAcceleration = correct;
}




// ----------------------------------------------------------------------------
// Returns whether the particle acceleration (i.e. change of
// momentum) is corrected by the fluid density in case of immersed rigid
// bodies
bool Particle::getFluidCorrectedAcceleration()
{
  return ( Particle::m_fluidCorrectedAcceleration );
}




// ----------------------------------------------------------------------------
// Output operator
ostream& operator << ( ostream& f, Particle const& P )
{
  f << "ID = " << P.m_id << endl;
  f << "Type = " << P.m_GeomType << endl;
  f << "Position = " << *P.getPosition();
  f << "Tag = " << P.m_tag;
  return ( f );
}




// ----------------------------------------------------------------------------
// Sets the rotation quaternion
void Particle::setQuaternionRotation(
	double const& vecteur0,
	double const& vecteur1,
	double const& vecteur2,
	double const& scalaire )
{
  m_kinematics->setQuaternionRotation( vecteur0, vecteur1, vecteur2,
    scalaire );
}




// ----------------------------------------------------------------------------
// Sets the rotation quaternion
void Particle::setQuaternionRotation( Quaternion const& qrot )
{
  m_kinematics->setQuaternionRotation( qrot );
}




// ----------------------------------------------------------------------------
// Updates geographic localisation in the LinkedCell. Note that
// this method uses the cell from the previous time
void Particle::updateGeoPosition()
{
  m_GeoLoc = m_cellule_nm1->getGeoPosition();
}




// ----------------------------------------------------------------------------
// Returns particle geographic localisation at the current discrete time
GeoPosition Particle::getGeoPosition() const
{
  return ( m_GeoLoc );
}




// ----------------------------------------------------------------------------
// Returns particle geographic localisation at the previous discrete time
GeoPosition Particle::getGeoPositionNm1() const
{
  return ( m_GeoLoc_nm1 );
}




// ----------------------------------------------------------------------------
// Creates the VelocityInfosNm1 structure
void Particle::createVelocityInfosNm1()
{
  if ( !m_VelocityInfosNm1 )
    m_VelocityInfosNm1 = new struct VelocityInfosNm1;
}




// ----------------------------------------------------------------------------
// Returns an orientation vector to describe the angular position of the
// particle
Vector3 Particle::computeOrientationVector() const
{
  return( m_geoRBWC->getConvex()->computeOrientationVector(
  	m_geoRBWC->getTransform() ) );
}




// ----------------------------------------------------------------------------
// Returns particle coordination number
int Particle::getCoordinationNumber() const
{
  return ( m_coordination_number );
}




// ----------------------------------------------------------------------------
// Initializes the particle torsor
void Particle::InitializeForce( bool const& withWeight )
{
  if ( withWeight )
    m_torsor.setToBodyForce( *m_geoRBWC->getCentre(), m_weight );
  else m_torsor.setToBodyForce( *m_geoRBWC->getCentre(), Vector3Nul );
  m_coordination_number = 0 ;
}




// ----------------------------------------------------------------------------
// Writes the particle in a Paraview format
void Particle::write_polygonsStr_PARAVIEW( list<int> &connectivity,
    	list<int> &offsets, list<int> &cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const
{
  m_geoRBWC->getConvex()->write_polygonsStr_PARAVIEW( connectivity,
	offsets, cellstype, firstpoint_globalnumber, last_offset );
}




// ----------------------------------------------------------------------------
// Returns the number of points to write the particle in a Paraview format
int Particle::numberOfPoints_PARAVIEW() const
{
  return ( m_geoRBWC->getConvex()->numberOfPoints_PARAVIEW() );
}




// ----------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the particle
// shape in a Paraview format
int Particle::numberOfCells_PARAVIEW() const
{
  return ( m_geoRBWC->getConvex()->numberOfCells_PARAVIEW() );
}




// ----------------------------------------------------------------------------
// Returns a list of points describing the particle in a Paraview format
list<Point3> Particle::get_polygonsPts_PARAVIEW( Vector3 const* translation )
	const
{
  return ( m_geoRBWC->get_polygonsPts_PARAVIEW( translation ) );
}




// ----------------------------------------------------------------------------
// Writes the points describing the particle in a Paraview format
void Particle::write_polygonsPts_PARAVIEW( ostream &f,
  	Vector3 const* translation ) const
{
  m_geoRBWC->write_polygonsPts_PARAVIEW( f, translation );
}




// ----------------------------------------------------------------------------
// Sets the boolean that tells that the rigid body's transformation
// with the scaling by the crust thickness to shrink the rigid bodies has
// already been computed to false
void Particle::initialize_transformWithCrust_to_notComputed()
{
  m_geoRBWC->initialize_transformWithCrust_to_notComputed();
}




// ----------------------------------------------------------------------------
// Sets the velocity with a random motion
void Particle::setRandomMotion( double const& coefTrans,
	double const& coefRot )
{
  if ( coefTrans )
    if ( m_tag != 2 )
    {
      Vector3 rvel;
      rvel[X] = coefTrans * (
	    	2. * (double(random()) / double(INT_MAX)) - 1. ) ;
      rvel[Y] = coefTrans * (
	    	2. * (double(random()) / double(INT_MAX)) - 1. ) ;
      if ( GrainsBuilderFactory::getContext() == DIM_3 )
        rvel[Z] = coefTrans * (
	      	2. * (double(random()) / double(INT_MAX)) - 1. ) ;
      setTranslationalVelocity( rvel );
    }

  if ( coefRot )
    if ( m_tag != 2 )
    {
      Vector3 rvel;
      if ( GrainsBuilderFactory::getContext() == DIM_3 )
      {      
        rvel[X] = coefTrans * (
	    	2. * (double(random()) / double(INT_MAX)) - 1. ) ;
        rvel[Y] = coefTrans * (
	    	2. * (double(random()) / double(INT_MAX)) - 1. ) ;
      }
      rvel[Z] = coefTrans * (
	      	2. * (double(random()) / double(INT_MAX)) - 1. ) ;
      setAngularVelocity( rvel );
    }
}




// ----------------------------------------------------------------------------
// Returns the number of corners of the rigib body shape and a code
// describing the rigid body shape
int Particle::getNbCorners() const
{
  return ( m_geoRBWC->getConvex()->getNbCorners() );
}




// ----------------------------------------------------------------------------
// Outputs information to be transferred to the fluid
void Particle::writePositionInFluid( ostream& fileOut )
{
  m_geoRBWC->writePositionInFluid( fileOut );
}




// ----------------------------------------------------------------------------
// Sets particle activity
void Particle::setActivity( ParticleActivity activity )
{
  m_activity = activity;
}




// ----------------------------------------------------------------------------
// Sets and returns the particle tag
int Particle::setTag( int tag_ )
{
  m_tag = tag_;
  return ( m_tag );
}




// ----------------------------------------------------------------------------
// Sets the geographic location of the particle in the LinkedCell
void Particle::setGeoPosition( GeoPosition const& geoloc_ )
{
  m_GeoLoc = geoloc_;
}




// ----------------------------------------------------------------------------
// Sets the cell the particle belonged to at the previous discrete time
void Particle::setCellNm1( Cell* cel )
{
  m_cellule_nm1 = cel;
}




// ----------------------------------------------------------------------------
// Sets the particle geometric type
void Particle::setGeometricType( int const& pc )
{
  m_GeomType = pc;
}




// ----------------------------------------------------------------------------
// Returns particle activity
ParticleActivity Particle::getActivity() const
{
  return ( m_activity );
}




// ----------------------------------------------------------------------------
// Returns particle tag at the current discrete time
int Particle::getTag() const
{
  return ( m_tag );
}




// ----------------------------------------------------------------------------
// Returns particle tag at the previous discrete time
int Particle::getTagNm1() const
{
  return ( m_tag_nm1 );
}




// ----------------------------------------------------------------------------
// Returns the cell the particle belonged to at the previous discrete time
Cell* Particle::getCellNm1() const
{
  return ( m_cellule_nm1 );
}




// ----------------------------------------------------------------------------
// Returns the cell the particle belongs to at the current discrete time
Cell* Particle::getCell() const
{
  return ( m_cellule );
}




// ----------------------------------------------------------------------------
// Returns particle class
int Particle::getGeometricType() const
{
  return ( m_GeomType );
}




// ----------------------------------------------------------------------------
// Sets the cell the particle belonged to, the particle tag and
// the geographic location of the particle at the current time
void Particle::setCellTagGeoPosition( Cell* cell_, int tag_,
    	GeoPosition const& geoloc_ )
{
  m_cellule = cell_;
  m_tag = tag_;
  m_GeoLoc = geoloc_;
}




// ----------------------------------------------------------------------------
// Sets the cell the particle belonged to, the particle tag and
// the geographic location of the particle at the previous time
void Particle::setCellTagGeoPosition_nm1( Cell* cell_, int tag_,
    	GeoPosition const& geoloc_ )
{
  m_cellule_nm1 = cell_;
  m_tag_nm1 = tag_;
  m_GeoLoc_nm1 = geoloc_;
}




// ----------------------------------------------------------------------------
// Copies the cell the particle belonged to, the particle tag and
// the geographic location of the particle from current time to previous time
void Particle::copyCellTagGeoPosition_n_to_nm1()
{
  m_cellule_nm1 = m_cellule;
  m_tag_nm1 = m_tag;
  m_GeoLoc_nm1 = m_GeoLoc;
}




// ----------------------------------------------------------------------------
// Increments the coordination number by nc
void Particle::addToCoordinationNumber( int const& nc )
{
  m_coordination_number += nc;
}




// ----------------------------------------------------------------------------
// Returns a pointer to the reference component of the component:
// this in general and the CompositeParticle for an elementary particle
Component* Particle::getMasterComponent()
{
  return ( m_masterParticle );
}




// ----------------------------------------------------------------------------
// Sets the pointer to the master particle of the particle
// (in general the CompositeParticle for an elementary particle)
void Particle::setMasterParticle( Particle* master_ )
{
  m_masterParticle = master_ ;
}




// ----------------------------------------------------------------------------
// Returns whether the particle is an elementary particle */
bool Particle::isElementaryParticle() const
{
  return ( m_masterParticle != this );
}




// ----------------------------------------------------------------------------
// Sets the particle density
void Particle::setDensity( double const& density_)
{
  m_density = density_;
}




// ----------------------------------------------------------------------------
// Compose the component transformation on the left by another
// transformation: this = t o this (this first followed by t)
void Particle::composePositionLeftByTransform( Transform const& t )
{
  m_geoRBWC->composeLeftByTransform( t );
}




// ----------------------------------------------------------------------------
// Compose the component transformation on the right by another
// transformation: this = this o t (t first followed by this)
void Particle::composePositionRightByTransform( Transform const& t )
{
  m_geoRBWC->composeRightByTransform( t );
}




// ----------------------------------------------------------------------------
// Sets the particle's transformation with a transformation
void Particle::setTransform( Transform const& transform_ )
{
   m_geoRBWC->setTransform( transform_ );
}




// ----------------------------------------------------------------------------
// Sets the pointer to the particle's kinematics
void Particle::setKinematics( ParticleKinematics* pkine )
{
  m_kinematics = pkine;
}




// ----------------------------------------------------------------------------
// Computes particle inertia tensor in the space fixed coordinate frame
void Particle::computeInertiaTensorSpaceFixed( vector<double>& inertia ) const
{
  // Get rotation matrix from ref position to current position
  Matrix mr = m_geoRBWC->getTransform()->getBasis();

  // Compute the rotation matrix transposed
  Matrix mrt = mr.transpose();

  // Transfer inertia vector into a matrix
  Matrix inertiaBodyFixed( m_inertia[0], m_inertia[1], m_inertia[2],
  	m_inertia[1], m_inertia[3], m_inertia[4],
	m_inertia[2], m_inertia[4], m_inertia[5] );

  // Compute the inertia matrix in the space fixed frame
  // I_space = Mr * I_body * Mr^t
  Matrix inertiaSpaceFixed = mr * ( inertiaBodyFixed * mrt );

  // Transfer the inertia matrix into the inertia vector
  inertia[0] = inertiaSpaceFixed[X][X];
  inertia[1] = inertiaSpaceFixed[X][Y];
  inertia[2] = inertiaSpaceFixed[X][Z];
  inertia[3] = inertiaSpaceFixed[Y][Y];
  inertia[4] = inertiaSpaceFixed[Y][Z];
  inertia[5] = inertiaSpaceFixed[Z][Z];
}




// ----------------------------------------------------------------------------
// Returns the specific composite shape name: "none" for standard particles
// and non-specific composite particle and class name for specific composite
// particles 
string Particle::getSpecificCompositeShapeName() const
{
  return ( m_specific_composite_shape );
}
