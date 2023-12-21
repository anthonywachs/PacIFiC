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

int SpheroCylinder::m_visuNodeNbPerQar = 8;


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
SpheroCylinder::SpheroCylinder( DOMNode* root, int const& pc )
  : CompositeParticle( false )
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
	+ ( 8. / 15. ) * PI * pow( m_radius, 5. )
	+ ( 1. / 2. ) * PI * m_height * pow( m_radius, 4. );
  m_inertia[3] = ( 1. / 30. ) * PI * pow( m_radius, 4. )
  	* ( 16. * m_radius + 15. * m_height );  
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

  // Top sphere
  Sphere* spht = new Sphere( m_radius );
  RigidBodyWithCrust* geoRBWC_spht = 
  	new RigidBodyWithCrust( spht, Transform(),
  	false, crust_thickness );
  m_elementaryParticles[1] = new Particle( geoRBWC_spht, m_density,
      	m_materialName, pc );
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
      	m_materialName, pc );
  m_InitialRelativePositions[2][X] = 0.;
  m_InitialRelativePositions[2][Y] = - m_height / 2.;  
  m_InitialRelativePositions[2][Z] = 0.;  
  m_elementaryParticles[2]->setPosition( m_InitialRelativePositions[2] );  
  m_elementaryParticles[2]->setMasterParticle( this );
  m_InitialRotationMatrices[2] = m_elementaryParticles[1]->getRigidBody()
	->getTransform()->getBasis();
  m_elementaryParticles[2]->getRigidBody()
 	->composeLeftByTransform( *(m_geoRBWC->getTransform()) );

  // Set the the circumscribed radius
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
  SpheroCylinder const* SpheroCylRef =
  	dynamic_cast<SpheroCylinder const*>(ParticleRef);
  m_height = SpheroCylRef->m_height;
  m_radius = SpheroCylRef->m_radius;  
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
SpheroCylinder::SpheroCylinder( SpheroCylinder const& other, 
    	bool const& autonumbering )
  : CompositeParticle( other, autonumbering )
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
Particle* SpheroCylinder::createCloneCopy( bool const& autonumbering ) const
{
  Particle* particle = new SpheroCylinder( *this, autonumbering );

  return ( particle );
}




// ----------------------------------------------------------------------------
// Creates a clone of the composite particle. This method calls the
// constructor SpheroCylinder( int const& id_, Particle const* ParticleRef,
// Vector3 const& vtrans, Quaternion const& qrot, Vector3 const& vrot,
// Transform const& config, ParticleActivity const& activ ) and is used for
// periodic clone composite particles to be inserted in the simulation.
// Autonumbering is set to false and numbering is set with the parameter id_
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
 return ( 2 * ( 4 * m_visuNodeNbPerQar * m_visuNodeNbPerQar + 2 ) );  
}




// ----------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the
// composite particle shape in a Paraview format
int SpheroCylinder::numberOfCells_PARAVIEW() const
{  
  return ( 2 * ( 4 * m_visuNodeNbPerQar * m_visuNodeNbPerQar ) 
  	+ 4 * m_visuNodeNbPerQar );   
}




// ----------------------------------------------------------------------------
// Writes the points describing the composite particle in a Paraview format
void SpheroCylinder::write_polygonsPts_PARAVIEW( ostream& f,
	Vector3 const* translation )const
{ 
  // Top sphere
  Transform const* transform = m_elementaryParticles[1]->getRigidBody()
  	->getTransform();  	 	
  double angle = PI / ( 2. * m_visuNodeNbPerQar ) ;
  double angleY = 0., local_radius = 0.;
  int k, i, ptsPerlevel = 4 * m_visuNodeNbPerQar;
  Point3 pp, pptrans;
  
  // Regular points on the surface
  for ( k = m_visuNodeNbPerQar - 1; k < 2*m_visuNodeNbPerQar-1 ; ++k ) 
  {  
    angleY = - PI / 2. + ( k + 1 ) * angle;
    local_radius = m_radius * cos( angleY );
    pp[Y] = m_radius * sin( angleY );
    for ( i = 0; i < ptsPerlevel ; ++i )
    {
      pp[X] = local_radius * cos( i * angle );
      pp[Z] = local_radius * sin( i * angle );
      pptrans = (*transform)( pp );
      if ( translation ) pptrans += *translation;
      f << pptrans[X] << " " << pptrans[Y] << " " << pptrans[Z] << endl;
    }
  }
   
  pp[X] = 0.;
  pp[Z] = 0.;	
  // Top point
  pp[Y] = m_radius;
  pptrans = (*transform)( pp );
  if ( translation ) pptrans += *translation;
  f << pptrans[X] << " " << pptrans[Y] << " " << pptrans[Z] << endl;
	
  // Gravity center
  pp[Y] = 0.;
  pptrans = (*transform)( pp );
  if ( translation ) pptrans += *translation;
  f << pptrans[X] << " " << pptrans[Y] << " " << pptrans[Z] << endl; 
  
  
  // Bottom sphere
  transform = m_elementaryParticles[2]->getRigidBody()
  	->getTransform();  	 	
  
  // Regular points on the surface
  for ( k = m_visuNodeNbPerQar - 1; k >=0 ; --k ) 
  {  
    angleY = - PI / 2. + ( k + 1 ) * angle;
    local_radius = m_radius * cos( angleY );
    pp[Y] = m_radius * sin( angleY );
    for ( i = 0; i < ptsPerlevel ; ++i )
    {
      pp[X] = local_radius * cos( i * angle );
      pp[Z] = local_radius * sin( i * angle );
      pptrans = (*transform)( pp );
      if ( translation ) pptrans += *translation;
      f << pptrans[X] << " " << pptrans[Y] << " " << pptrans[Z] << endl;
    }
  }
   
  pp[X] = 0.;
  pp[Z] = 0.;
  // Bottom point
  pp[Y] = - m_radius;
  pptrans = (*transform)( pp );
  if ( translation ) pptrans += *translation;
  f << pptrans[X] << " " << pptrans[Y] << " " << pptrans[Z] << endl;	
	
  // Gravity center
  pp[Y] = 0.;
  pptrans = (*transform)( pp );
  if ( translation ) pptrans += *translation;
  f << pptrans[X] << " " << pptrans[Y] << " " << pptrans[Z] << endl;
  
  
  // Cylinder: no additional point needed      
}




// ----------------------------------------------------------------------------
// Returns a list of points describing the component in a Paraview format
list<Point3> SpheroCylinder::get_polygonsPts_PARAVIEW(
	Vector3 const* translation ) const
{
  // Top sphere
  Transform const* transform = m_elementaryParticles[1]->getRigidBody()
  	->getTransform();
  list<Point3> ParaviewPoints;
  double angle = PI / ( 2. * m_visuNodeNbPerQar ) ;
  double angleY = 0., local_radius = 0.;
  int k, i, ptsPerlevel =  4 * m_visuNodeNbPerQar;
  Point3 pp, pptrans;
  
  // Regular points on the surface
  for ( k = m_visuNodeNbPerQar - 1; k < 2*m_visuNodeNbPerQar-1 ; ++k ) 
  {  
    angleY = -PI / 2. + ( k + 1 ) * angle;
    local_radius = m_radius * cos( angleY );
    pp[Y] = m_radius * sin( angleY );
    for ( i = 0; i < ptsPerlevel ; ++i )
    {
      pp[X] = local_radius * cos( i * angle );
      pp[Z] = local_radius * sin( i * angle );
      pptrans = (*transform)( pp );
      if ( translation ) pptrans += *translation;
      ParaviewPoints.push_back( pptrans );
    }
  }
   
  pp[X] = 0.;
  pp[Z] = 0.;	
  // Top point
  pp[Y] = m_radius;
  pptrans = (*transform)( pp );
  if ( translation ) pptrans += *translation;
  ParaviewPoints.push_back( pptrans );
	
  // Gravity center
  pp[Y] = 0.;
  pptrans = (*transform)( pp );  
  if ( translation ) pptrans += *translation;
  ParaviewPoints.push_back( pptrans );
  
  
  // Bottom sphere
  transform = m_elementaryParticles[2]->getRigidBody()
  	->getTransform();
  
  // Regular points on the surface
  for ( k = m_visuNodeNbPerQar - 1; k >=0 ; --k ) 
  {  
    angleY = -PI / 2. + ( k + 1 ) * angle;
    local_radius = m_radius * cos( angleY );
    pp[Y] = m_radius * sin( angleY );
    for ( i = 0; i < ptsPerlevel ; ++i )
    {
      pp[X] = local_radius * cos( i * angle );
      pp[Z] = local_radius * sin( i * angle );
      pptrans = (*transform)( pp );
      if ( translation ) pptrans += *translation;
      ParaviewPoints.push_back( pptrans );
    }
  }
   
  pp[X] = 0.;
  pp[Z] = 0.;
  // Bottom point
  pp[Y] = - m_radius;
  pptrans = (*transform)( pp );
  if ( translation ) pptrans += *translation;
  ParaviewPoints.push_back( pptrans );
	
  // Gravity center
  pp[Y] = 0.;
  pptrans = (*transform)( pp );  
  if ( translation ) pptrans += *translation;
  ParaviewPoints.push_back( pptrans );
  
  
  // Cylinder: no additional point needed     
      
  return ( ParaviewPoints );
}




// ----------------------------------------------------------------------------
// Writes the composite particle in a Paraview format
void SpheroCylinder::write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const
{
  int i, k, ptsPerlevel = 4 * m_visuNodeNbPerQar;
  
  // Top sphere
  int Top_number = firstpoint_globalnumber + ptsPerlevel * m_visuNodeNbPerQar,  
  	Top_GC_number = firstpoint_globalnumber 
		+ ptsPerlevel * m_visuNodeNbPerQar + 1;
  // Regular cells: Pyramid
  for ( k = 0; k < m_visuNodeNbPerQar-1 ; ++k ) 
  {  
    for ( i = 0; i < ptsPerlevel-1 ; ++i )
    {
      connectivity.push_back( firstpoint_globalnumber + k * ptsPerlevel + i ); 
      connectivity.push_back( firstpoint_globalnumber + k * ptsPerlevel + i 
      	+ 1);
      connectivity.push_back( firstpoint_globalnumber + ( k + 1) * ptsPerlevel
      	+ i + 1 );	
      connectivity.push_back( firstpoint_globalnumber + ( k + 1 ) * ptsPerlevel
	+ i );
      connectivity.push_back( Top_GC_number );
      last_offset += 5;
      offsets.push_back( last_offset );
      cellstype.push_back( 14 );		
    }
    connectivity.push_back( firstpoint_globalnumber + k * ptsPerlevel + 
    	ptsPerlevel - 1 );
    connectivity.push_back( firstpoint_globalnumber + k * ptsPerlevel );
    connectivity.push_back( firstpoint_globalnumber + ( k + 1 ) * ptsPerlevel );
    connectivity.push_back( firstpoint_globalnumber + ( k + 1 ) * ptsPerlevel 
    	+ ptsPerlevel - 1 );
    connectivity.push_back( Top_GC_number );
    last_offset += 5;
    offsets.push_back( last_offset );
    cellstype.push_back( 14 );    
  }   
  
  // Top cells: tetraedron  
  for ( i = 0; i < ptsPerlevel-1 ; ++i )
  {
    connectivity.push_back( firstpoint_globalnumber
    	+ ( m_visuNodeNbPerQar - 1 ) * ptsPerlevel + i );
    connectivity.push_back( firstpoint_globalnumber
    	+ ( m_visuNodeNbPerQar - 1 ) * ptsPerlevel + i + 1 );
    connectivity.push_back( Top_number );	
    connectivity.push_back( Top_GC_number );
    last_offset += 4;
    offsets.push_back( last_offset );
    cellstype.push_back( 10 ); 
  }
  connectivity.push_back( firstpoint_globalnumber
  	+  m_visuNodeNbPerQar * ptsPerlevel - 1 );
  connectivity.push_back( firstpoint_globalnumber
  	+ ( m_visuNodeNbPerQar - 1 ) * ptsPerlevel );
  connectivity.push_back( Top_number );	
  connectivity.push_back( Top_GC_number );
  last_offset += 4;
  offsets.push_back( last_offset );
  cellstype.push_back( 10 ); 
	

  // Bottom sphere
  int shift = firstpoint_globalnumber + ptsPerlevel * m_visuNodeNbPerQar + 2 ;  
  int Bottom_number = shift + ptsPerlevel * m_visuNodeNbPerQar,  
  	Bottom_GC_number = shift + ptsPerlevel * m_visuNodeNbPerQar + 1; 
  // Regular cells: Pyramid
  for ( k = 0; k < m_visuNodeNbPerQar-1 ; ++k ) 
  {  
    for ( i = 0; i < ptsPerlevel-1 ; ++i )
    {
      connectivity.push_back( shift + k * ptsPerlevel + i ); 
      connectivity.push_back( shift + k * ptsPerlevel + i + 1);
      connectivity.push_back( shift + ( k + 1) * ptsPerlevel + i + 1 );	
      connectivity.push_back( shift + ( k + 1 ) * ptsPerlevel + i );
      connectivity.push_back( Bottom_GC_number );
      last_offset += 5;
      offsets.push_back( last_offset );
      cellstype.push_back( 14 );		
    }
    connectivity.push_back( shift + k * ptsPerlevel + ptsPerlevel - 1 );
    connectivity.push_back( shift + k * ptsPerlevel );
    connectivity.push_back( shift + ( k + 1 ) * ptsPerlevel );
    connectivity.push_back( shift + ( k + 1 ) * ptsPerlevel + ptsPerlevel - 1 );
    connectivity.push_back( Bottom_GC_number );
    last_offset += 5;
    offsets.push_back( last_offset );
    cellstype.push_back( 14 );    
  }

  // Bottom cells: tetraedron
  for ( i = 0; i < ptsPerlevel-1 ; ++i )
  {
    connectivity.push_back( shift + ( m_visuNodeNbPerQar - 1 ) * ptsPerlevel 
    	+ i );
    connectivity.push_back( shift + ( m_visuNodeNbPerQar - 1 ) * ptsPerlevel 
    	+ i + 1 );
    connectivity.push_back( Bottom_number );	
    connectivity.push_back( Bottom_GC_number );
    last_offset += 4;
    offsets.push_back( last_offset );
    cellstype.push_back( 10 ); 
  }
  connectivity.push_back( shift +  m_visuNodeNbPerQar * ptsPerlevel - 1 );
  connectivity.push_back( shift + ( m_visuNodeNbPerQar - 1 ) * ptsPerlevel );
  connectivity.push_back( Bottom_number );	
  connectivity.push_back( Bottom_GC_number );
  last_offset += 4;
  offsets.push_back( last_offset );
  cellstype.push_back( 10 ); 
  
  
  // Cylinder	
  for ( i = 0; i < ptsPerlevel-1 ; ++i )
  {
    connectivity.push_back( firstpoint_globalnumber + i );
    connectivity.push_back( firstpoint_globalnumber + i + 1 );
    connectivity.push_back( Top_GC_number );
    connectivity.push_back( shift + i );
    connectivity.push_back( shift + i + 1 );
    connectivity.push_back( Bottom_GC_number );
    last_offset += 6;
    offsets.push_back( last_offset );
    cellstype.push_back( 13 );
  }
  connectivity.push_back( firstpoint_globalnumber + ptsPerlevel - 1 );
  connectivity.push_back( firstpoint_globalnumber );
  connectivity.push_back( Top_GC_number );
  connectivity.push_back( shift + ptsPerlevel - 1 );
  connectivity.push_back( shift );
  connectivity.push_back( Bottom_GC_number );
  last_offset += 6;
  offsets.push_back( last_offset );
  cellstype.push_back( 13 );


  firstpoint_globalnumber += 2 * ( ptsPerlevel * m_visuNodeNbPerQar + 2 ) ;
}






// ----------------------------------------------------------------------------
// Returns the number of corners of the rigib body shape and a code
// describing the rigid body shape
int SpheroCylinder::getNbCorners() const
{
  return ( 1001 );
}





// ----------------------------------------------------------------------------
// Outputs information to be transferred to the fluid
// Same format data as a regular cylinder: center of bottom circular
// face of the elementary cylinder, an arbitrary point on the lateral surface 
// of the elementary cylinder and center of top circular face of the elementary
// cylinder    
void SpheroCylinder::writePositionInFluid( ostream& fluid )
{
  Point3 pointEnvelope;
  vector<Point3> allPoints( 3, pointEnvelope );
  vector<Point3>::iterator point;
  Transform const* transform = m_geoRBWC->getTransform();

  fluid << " 3" << endl;

  allPoints[0][Y] = - m_height / 2.;
  allPoints[1][Y] = - m_height / 2.;
  allPoints[1][X] = m_radius;
  allPoints[2][Y] = m_height / 2.;

  for (point=allPoints.begin(); point!=allPoints.end(); point++)
  {
    pointEnvelope = (*transform)(*point);
    fluid << GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
		pointEnvelope[X] ) << " " <<
	GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
		pointEnvelope[Y] ) << " " <<
	GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
		pointEnvelope[Z] ) << endl;	
  }
  
  // No faces
  fluid << "0" << endl;
}




// ----------------------------------------------------------------------------
// Reads composite particle data from a stream. Usage: for standard composite
// particles in the 2014 reload format
void SpheroCylinder::read2014( istream& fileIn, vector<Particle*> const*
  	referenceParticles )
{
  CompositeParticle::read2014( fileIn, referenceParticles );

  // We know that (*referenceParticles)[m_GeomType] points to a
  // CompositeParticle, such that we can dynamic cast it to actual type and
  // use -> instead of using get methods through virtual typing
  SpheroCylinder const* CompParticleRef =
  	dynamic_cast<SpheroCylinder const*>(
		(*referenceParticles)[m_GeomType]);

  m_height = CompParticleRef->m_height;
  m_radius = CompParticleRef->m_radius;  
}




// ----------------------------------------------------------------------------
// Reads composite particle data from a stream in a binary form.
// Usage: for standard composite particles in the 2014 reload format
void SpheroCylinder::read2014_binary( istream& fileIn,
	vector<Particle*> const* referenceParticles )
{
  CompositeParticle::read2014_binary( fileIn, referenceParticles );

  // We know that (*referenceParticles)[m_GeomType] points to a
  // CompositeParticle, such that we can dynamic cast it to actual type and
  // use -> instead of using get methods through virtual typing
  SpheroCylinder const* CompParticleRef =
  	dynamic_cast<SpheroCylinder const*>(
		(*referenceParticles)[m_GeomType]);

  m_height = CompParticleRef->m_height;
  m_radius = CompParticleRef->m_radius;
}




// ----------------------------------------------------------------------------
// Computes and sets the circumscribed radius
void SpheroCylinder::setCircumscribedRadius()
{
  m_geoRBWC->setCircumscribedRadius( m_height / 2. + m_radius );
}




// ----------------------------------------------------------------------------
// Saves additional features of a (in practice reference) composite particle
// for reload
void SpheroCylinder::writeAdditionalFeatures( ostream& fileSave ) const
{
  fileSave << endl << "*HeightAndRadius " << m_height << " " << m_radius;  
  CompositeParticle::writeAdditionalFeatures( fileSave );
}




// ----------------------------------------------------------------------------
// Reads additional features of a (in practice reference) composite particle
// data from a stream
void SpheroCylinder::readAdditionalFeatures( istream& fileIn )
{
  string buffer;
  fileIn >> buffer >> m_height >> m_radius;
  CompositeParticle::readAdditionalFeatures( fileIn );
}
