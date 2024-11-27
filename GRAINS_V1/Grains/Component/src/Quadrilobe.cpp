#include "GrainsMPIWrapper.hh"
#include "Quadrilobe.hh"
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

int Quadrilobe::m_visuNodeNbPerHalf = 16;


// ----------------------------------------------------------------------------
// Constructor with autonumbering as input parameter
Quadrilobe::Quadrilobe( bool const& autonumbering )
  : CompositeParticle( autonumbering )
{
  m_specific_composite_shape = "Quadrilobe";
  m_radius = 0.;
  m_height = 0.;
  m_length = 0.;
}





// ----------------------------------------------------------------------------
// Constructor with an XML node as an input parameter. This constructor is
// expected to be used for reference composite particles
Quadrilobe::Quadrilobe( DOMNode* root, int const& pc )
  : CompositeParticle( false )
{
  m_specific_composite_shape = "Quadrilobe";

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
  m_length = ReaderXML::getNodeAttr_Double( nGeometry, "Length" ); 
  m_mass = m_density * ( m_height * pow( m_radius, 2. ) * ( 2. * PI
  	- 4. ) + 4. * m_radius * m_length * m_height );
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

  // Moment of inertia tensor of the Quadrilobe
  m_inertia[1] = m_inertia[2] = m_inertia[4] = 0.;
  double h2 = m_height * m_height, r2 = m_radius * m_radius, 
  	halfl = 0.5 * m_length;
  double a1 = 4 * r2 * m_height,
  	a2 = ( halfl - m_radius ) * 2. * m_radius * m_height,
	a3 = 0.5 * PI * r2 * m_height;  
  m_inertia[0] = a1 * ( 4. * r2 + h2 ) / 12.
  	+ 0.5 * a2 * ( ( pow( halfl - m_radius, 2. ) + h2 ) / 3.   
		+ pow( halfl +  m_radius, 2. ) 
		+ ( 4. * r2 + h2 ) / 3. )
	+ 2. * a3 * ( ( 0.25 - 16. / ( 9. * PI * PI ) ) * r2
		+ h2 / 6. + pow( halfl + 4. * m_radius / ( 3. * PI ), 2. ) 
		+ 0.25 * r2 ) ;
  m_inertia[3] = 2. * a1 * r2 / 3.
  	+ a2 * ( ( pow( halfl - m_radius, 2. ) + 4. * r2 ) / 3. 
		+ pow( halfl + m_radius, 2.) )
	+ 4. * a3 * ( ( 0.5 - 16. / ( 9. * PI * PI ) ) * r2 
		+ pow( halfl + 4. * m_radius / ( 3. * PI ), 2. ) ); 
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
  	m_length, m_height );
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
  	m_length, m_height );
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

  // In case part of the particle acceleration computed explicity
  if ( Particle::m_splitExplicitAcceleration ) createVelocityInfosNm1();
}





// ----------------------------------------------------------------------------
// Constructor with input parameters
Quadrilobe::Quadrilobe( int const& id_,
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
  m_specific_composite_shape = "Quadrilobe";
  Quadrilobe const* QuadrilobeRef =
  	dynamic_cast<Quadrilobe const*>(ParticleRef);
  m_radius = QuadrilobeRef->m_radius; 
  m_length = QuadrilobeRef->m_length;    
  m_height = QuadrilobeRef->m_height;
}




// ----------------------------------------------------------------------------
// Constructor with input parameters
Quadrilobe::Quadrilobe( int const& id_,
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
  m_specific_composite_shape = "Quadrilobe";
  Quadrilobe const* QuadrilobeRef =
  	dynamic_cast<Quadrilobe const*>(ParticleRef);
  m_radius = QuadrilobeRef->m_radius; 
  m_length = QuadrilobeRef->m_length;    
  m_height = QuadrilobeRef->m_height;
}




// ----------------------------------------------------------------------------
// Destructor
Quadrilobe::~Quadrilobe()
{}




// ----------------------------------------------------------------------------
// Copy constructor (the torsor is initialized to 0)
Quadrilobe::Quadrilobe( Quadrilobe const& other, 
    	bool const& autonumbering )
  : CompositeParticle( other, autonumbering )
{
  m_radius = other.m_radius;
  m_length = other.m_length;
  m_height = other.m_height;
}




// ----------------------------------------------------------------------------
// Creates a clone of the particle. This method calls the standard
// copy constructor and is used for new particles to be inserted in the
// simulation. Numbering is automatic, total number of components is
// incremented by 1 and activity is set to WAIT. The calling object is
// expected to be a reference particle
Particle* Quadrilobe::createCloneCopy( bool const& autonumbering ) const
{
  Particle* particle = new Quadrilobe( *this, autonumbering );

  return ( particle );
}




// ----------------------------------------------------------------------------
// Creates a clone of the composite particle. This method calls the
// constructor Quadrilobe( int const& id_, Particle const* ParticleRef,
// Vector3 const& vtrans, Quaternion const& qrot, Vector3 const& vrot,
// Transform const& config, ParticleActivity const& activ ) and is used for
// periodic clone composite particles to be inserted in the simulation.
// Autonumbering is set to false and numbering is set with the parameter id_
Particle* Quadrilobe::createCloneCopy( int const& id_,
    	Particle const* ParticleRef, Vector3 const& vtrans,
	Quaternion const& qrot,	Vector3 const& vrot,
	Transform const& config, ParticleActivity const& activ,
	map< std::tuple<int,int,int>,
     	std::tuple<bool, Vector3, Vector3, Vector3> > const* contactMap ) const
{
  Particle* particle = new Quadrilobe( id_, ParticleRef, vtrans,
	qrot, vrot, config, activ, contactMap );

  return ( particle );
}




// ----------------------------------------------------------------------------
// Returns the number of points to write the composite particle in a
// Paraview format
int Quadrilobe::numberOfPoints_PARAVIEW() const
{
 return ( 8 * ( m_visuNodeNbPerHalf + 2 ) );  
}




// ----------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the
// composite particle shape in a Paraview format
int Quadrilobe::numberOfCells_PARAVIEW() const
{  
  return ( 2 * ( 2 * m_visuNodeNbPerHalf + 1 ) );   
}




// ----------------------------------------------------------------------------
// Writes the points describing the composite particle in a Paraview format
void Quadrilobe::write_polygonsPts_PARAVIEW( ostream& f,
	Vector3 const* translation )const
{ 
  Point3 pp, p;
  double dtheta = PI / m_visuNodeNbPerHalf, tstartleft = 1.5 * PI,
  	tstartright = - 0.5 * PI;

  for (size_t m=0;m<2;++m)
  {
    Transform const* transform = m_elementaryParticles[m]->getRigidBody()
  	->getTransform();
	  
    // Left half cylinder
    // Lower disk rim
    p[Y] = - 0.5 * m_height;
    for (int i=0;i<m_visuNodeNbPerHalf+1;++i)
    {
      p[X] = m_radius * cos ( tstartleft - i * dtheta ) - 0.5 * m_length;
      p[Z] = m_radius * sin ( tstartleft - i * dtheta );
      pp = (*transform)( p );
      if ( translation ) pp += *translation;
      f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
    }

    // Upper disk rim
    p[Y] = 0.5 * m_height;
    for (int i=0;i<m_visuNodeNbPerHalf+1;++i)
    {
      p[X] = m_radius * cos ( tstartleft - i * dtheta ) - 0.5 * m_length;
      p[Z] = m_radius * sin ( tstartleft - i * dtheta );
      pp = (*transform)( p );
      if ( translation ) pp += *translation;
      f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
    }

    // Lower disk center
    p[X] = - 0.5 * m_length;
    p[Y] = - 0.5 * m_height;
    p[Z] = 0.;
    pp = (*transform)( p );
    if ( translation ) pp += *translation;
    f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;

    // Upper disk center
    p[Y] = 0.5 * m_height;
    pp = (*transform)( p );
    if ( translation ) pp += *translation;
    f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
  
  
    // Right half cylinder
    // Lower disk rim
    p[Y] = - 0.5 * m_height;
    for (int i=0;i<m_visuNodeNbPerHalf+1;++i)
    {
      p[X] = m_radius * cos ( tstartright + i * dtheta ) + 0.5 * m_length;
      p[Z] = m_radius * sin ( tstartright + i * dtheta );
      pp = (*transform)( p );
      if ( translation ) pp += *translation;
      f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
    }

    // Upper disk rim
    p[Y] = 0.5 * m_height;
    for (int i=0;i<m_visuNodeNbPerHalf+1;++i)
    {
      p[X] = m_radius * cos ( tstartright + i * dtheta ) + 0.5 * m_length;
      p[Z] = m_radius * sin ( tstartright + i * dtheta );
      pp = (*transform)( p );
      if ( translation ) pp += *translation;
      f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
    }

    // Lower disk center
    p[X] = 0.5 * m_length;
    p[Y] = - 0.5 * m_height;
    p[Z] = 0.;
    pp = (*transform)( p );
    if ( translation ) pp += *translation;
    f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;

    // Upper disk center
    p[Y] = 0.5 * m_height;
    pp = (*transform)( p );
    if ( translation ) pp += *translation;
    f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
  }          
}




// ----------------------------------------------------------------------------
// Returns a list of points describing the component in a Paraview format
list<Point3> Quadrilobe::get_polygonsPts_PARAVIEW(
	Vector3 const* translation ) const
{
  list<Point3> ParaviewPoints;
  Point3 pp, p;
  double dtheta = PI / m_visuNodeNbPerHalf, tstartleft = 1.5 * PI,
  	tstartright = - 0.5 * PI;

  for (size_t m=0;m<2;++m)
  {
    Transform const* transform = m_elementaryParticles[m]->getRigidBody()
  	->getTransform();
	  
    // Left half cylinder
    // Lower disk rim
    p[Y] = - 0.5 * m_height;
    for (int i=0;i<m_visuNodeNbPerHalf+1;++i)
    {
      p[X] = m_radius * cos ( tstartleft - i * dtheta ) - 0.5 * m_length;
      p[Z] = m_radius * sin ( tstartleft - i * dtheta );
      pp = (*transform)( p );
      if ( translation ) pp += *translation;
      ParaviewPoints.push_back( pp );
    }

    // Upper disk rim
    p[Y] = 0.5 * m_height;
    for (int i=0;i<m_visuNodeNbPerHalf+1;++i)
    {
      p[X] = m_radius * cos ( tstartleft - i * dtheta ) - 0.5 * m_length;
      p[Z] = m_radius * sin ( tstartleft - i * dtheta );
      pp = (*transform)( p );
      if ( translation ) pp += *translation;
      ParaviewPoints.push_back( pp );
    }

    // Lower disk center
    p[X] = - 0.5 * m_length;
    p[Y] = - 0.5 * m_height;
    p[Z] = 0.;
    pp = (*transform)( p );
    if ( translation ) pp += *translation;
    ParaviewPoints.push_back( pp );

    // Upper disk center
    p[Y] = 0.5 * m_height;
    pp = (*transform)( p );
    if ( translation ) pp += *translation;
    ParaviewPoints.push_back( pp );
  
  
    // Right half cylinder
    // Lower disk rim
    p[Y] = - 0.5 * m_height;
    for (int i=0;i<m_visuNodeNbPerHalf+1;++i)
    {
      p[X] = m_radius * cos ( tstartright + i * dtheta ) + 0.5 * m_length;
      p[Z] = m_radius * sin ( tstartright + i * dtheta );
      pp = (*transform)( p );
      if ( translation ) pp += *translation;
      ParaviewPoints.push_back( pp );
    }

    // Upper disk rim
    p[Y] = 0.5 * m_height;
    for (int i=0;i<m_visuNodeNbPerHalf+1;++i)
    {
      p[X] = m_radius * cos ( tstartright + i * dtheta ) + 0.5 * m_length;
      p[Z] = m_radius * sin ( tstartright + i * dtheta );
      pp = (*transform)( p );
      if ( translation ) pp += *translation;
      ParaviewPoints.push_back( pp );
    }

    // Lower disk center
    p[X] = 0.5 * m_length;
    p[Y] = - 0.5 * m_height;
    p[Z] = 0.;
    pp = (*transform)( p );
    if ( translation ) pp += *translation;
    ParaviewPoints.push_back( pp );

    // Upper disk center
    p[Y] = 0.5 * m_height;
    pp = (*transform)( p );
    if ( translation ) pp += *translation;
    ParaviewPoints.push_back( pp );
  }
  
  return ( ParaviewPoints );  
}




// ----------------------------------------------------------------------------
// Writes the composite particle in a Paraview format
void Quadrilobe::write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const
{
  // First spherocylindrical prism 
  // Left half cylinder
  for (int i=0;i<m_visuNodeNbPerHalf;++i)
  {
    connectivity.push_back( firstpoint_globalnumber + i );
    connectivity.push_back( firstpoint_globalnumber + i + 1 );
    connectivity.push_back( firstpoint_globalnumber + 
    	2 * ( m_visuNodeNbPerHalf + 1 ) );
    connectivity.push_back( firstpoint_globalnumber + i + 
    	m_visuNodeNbPerHalf + 1 );
    connectivity.push_back( firstpoint_globalnumber + i + 
    	m_visuNodeNbPerHalf + 2 );
    connectivity.push_back( firstpoint_globalnumber +
    	2 * ( m_visuNodeNbPerHalf + 1 ) + 1 );
    last_offset += 6;
    offsets.push_back( last_offset );
    cellstype.push_back( 13 );
  }
  
  // Right half cylinder 
  int rightstart = firstpoint_globalnumber +  2 * ( m_visuNodeNbPerHalf + 1 ) 
  	+ 2;
  for (int i=0;i<m_visuNodeNbPerHalf;++i)
  {
    connectivity.push_back( rightstart + i );
    connectivity.push_back( rightstart + i + 1 );
    connectivity.push_back( rightstart + 2 * ( m_visuNodeNbPerHalf + 1 ) );
    connectivity.push_back( rightstart + i + m_visuNodeNbPerHalf + 1 );
    connectivity.push_back( rightstart + i + m_visuNodeNbPerHalf + 2 );
    connectivity.push_back( rightstart + 2 * ( m_visuNodeNbPerHalf + 1 ) + 1 );
    last_offset += 6;
    offsets.push_back( last_offset );
    cellstype.push_back( 13 );
  }
  
  // Box
  connectivity.push_back( firstpoint_globalnumber );
  connectivity.push_back( firstpoint_globalnumber + m_visuNodeNbPerHalf );
  connectivity.push_back( firstpoint_globalnumber + 
  	2 * m_visuNodeNbPerHalf + 1 );
  connectivity.push_back( firstpoint_globalnumber + m_visuNodeNbPerHalf + 1 );
  connectivity.push_back( firstpoint_globalnumber + 2 * m_visuNodeNbPerHalf + 
  	4 );
  connectivity.push_back( firstpoint_globalnumber + 3 * m_visuNodeNbPerHalf + 
  	4 );
  connectivity.push_back( firstpoint_globalnumber + 4 * m_visuNodeNbPerHalf + 
  	5 );	
  connectivity.push_back( firstpoint_globalnumber + 3 * m_visuNodeNbPerHalf + 
  	5 );	
  last_offset += 8;
  offsets.push_back( last_offset );
  cellstype.push_back( 12 );
  
  firstpoint_globalnumber += 4 * ( m_visuNodeNbPerHalf + 2 ) ;
  
  
  // Second spherocylindrical prism 
  // Left half cylinder
  for (int i=0;i<m_visuNodeNbPerHalf;++i)
  {
    connectivity.push_back( firstpoint_globalnumber + i );
    connectivity.push_back( firstpoint_globalnumber + i + 1 );
    connectivity.push_back( firstpoint_globalnumber + 
    	2 * ( m_visuNodeNbPerHalf + 1 ) );
    connectivity.push_back( firstpoint_globalnumber + i + 
    	m_visuNodeNbPerHalf + 1 );
    connectivity.push_back( firstpoint_globalnumber + i + 
    	m_visuNodeNbPerHalf + 2 );
    connectivity.push_back( firstpoint_globalnumber +
    	2 * ( m_visuNodeNbPerHalf + 1 ) + 1 );
    last_offset += 6;
    offsets.push_back( last_offset );
    cellstype.push_back( 13 );
  }
  
  // Right half cylinder 
  rightstart = firstpoint_globalnumber +  2 * ( m_visuNodeNbPerHalf + 1 ) + 2;
  for (int i=0;i<m_visuNodeNbPerHalf;++i)
  {
    connectivity.push_back( rightstart + i );
    connectivity.push_back( rightstart + i + 1 );
    connectivity.push_back( rightstart + 2 * ( m_visuNodeNbPerHalf + 1 ) );
    connectivity.push_back( rightstart + i + m_visuNodeNbPerHalf + 1 );
    connectivity.push_back( rightstart + i + m_visuNodeNbPerHalf + 2 );
    connectivity.push_back( rightstart + 2 * ( m_visuNodeNbPerHalf + 1 ) + 1 );
    last_offset += 6;
    offsets.push_back( last_offset );
    cellstype.push_back( 13 );
  }
  
  // Box
  connectivity.push_back( firstpoint_globalnumber );
  connectivity.push_back( firstpoint_globalnumber + m_visuNodeNbPerHalf );
  connectivity.push_back( firstpoint_globalnumber + 
  	2 * m_visuNodeNbPerHalf + 1 );
  connectivity.push_back( firstpoint_globalnumber + m_visuNodeNbPerHalf + 1 );
  connectivity.push_back( firstpoint_globalnumber + 2 * m_visuNodeNbPerHalf + 
  	4 );
  connectivity.push_back( firstpoint_globalnumber + 3 * m_visuNodeNbPerHalf + 
  	4 );
  connectivity.push_back( firstpoint_globalnumber + 4 * m_visuNodeNbPerHalf + 
  	5 );	
  connectivity.push_back( firstpoint_globalnumber + 3 * m_visuNodeNbPerHalf + 
  	5 );	
  last_offset += 8;
  offsets.push_back( last_offset );
  cellstype.push_back( 12 );
  
  firstpoint_globalnumber += 4 * ( m_visuNodeNbPerHalf + 2 ) ;  
}






// ----------------------------------------------------------------------------
// Returns the number of corners of the rigib body shape and a code
// describing the rigid body shape
int Quadrilobe::getNbCorners() const
{
  return ( 1001 );
}





// ----------------------------------------------------------------------------
// Outputs information to be transferred to the fluid
// Same format data as a regular cylinder: center of bottom circular
// face of the elementary cylinder, an arbitrary point on the lateral surface 
// of the elementary cylinder and center of top circular face of the elementary
// cylinder    
void Quadrilobe::writePositionInFluid( ostream& fluid )
{
  // TO DO
}




// ----------------------------------------------------------------------------
// Reads composite particle data from a stream. Usage: for standard composite
// particles in the 2014 reload format
void Quadrilobe::read2014( istream& fileIn, vector<Particle*> const*
  	referenceParticles )
{
  CompositeParticle::read2014( fileIn, referenceParticles );

  // We know that (*referenceParticles)[m_GeomType] points to a
  // CompositeParticle, such that we can dynamic cast it to actual type and
  // use -> instead of using get methods through virtual typing
  Quadrilobe const* CompParticleRef =
  	dynamic_cast<Quadrilobe const*>(
		(*referenceParticles)[m_GeomType]);

  m_radius = CompParticleRef->m_radius;
  m_length = CompParticleRef->m_length;
  m_height = CompParticleRef->m_height;
}




// ----------------------------------------------------------------------------
// Reads composite particle data from a stream in a binary form.
// Usage: for standard composite particles in the 2014 reload format
void Quadrilobe::read2014_binary( istream& fileIn,
	vector<Particle*> const* referenceParticles )
{
  CompositeParticle::read2014_binary( fileIn, referenceParticles );

  // We know that (*referenceParticles)[m_GeomType] points to a
  // CompositeParticle, such that we can dynamic cast it to actual type and
  // use -> instead of using get methods through virtual typing
  Quadrilobe const* CompParticleRef =
  	dynamic_cast<Quadrilobe const*>(
		(*referenceParticles)[m_GeomType]);

  m_radius = CompParticleRef->m_radius;
  m_length = CompParticleRef->m_length;
  m_height = CompParticleRef->m_height;
}




// ----------------------------------------------------------------------------
// Computes and sets the circumscribed radius
void Quadrilobe::setCircumscribedRadius()
{
  m_geoRBWC->setCircumscribedRadius( sqrt( 0.25 * m_height * m_height 
  	+ ( 0.5 * m_length + m_radius ) * ( 0.5 * m_length + m_radius ) ) );
}




// ----------------------------------------------------------------------------
// Saves additional features of a (in practice reference) composite particle
// for reload
void Quadrilobe::writeAdditionalFeatures( ostream& fileSave ) const
{
  fileSave << endl << "*RadiusLengthHeight " << m_radius << " " << m_length
  	<< " " << m_height;  
  CompositeParticle::writeAdditionalFeatures( fileSave );
}




// ----------------------------------------------------------------------------
// Reads additional features of a (in practice reference) composite particle
// data from a stream
void Quadrilobe::readAdditionalFeatures( istream& fileIn )
{
  string buffer;
  fileIn >> buffer >> m_radius >> m_length >> m_height;
  CompositeParticle::readAdditionalFeatures( fileIn );
}
