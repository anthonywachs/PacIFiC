#include "GrainsMPIWrapper.hh"
#include "TrilobeCylinder.hh"
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

int TrilobeCylinder::m_visuNodeNbPerHalf = 16;
bool TrilobeCylinder::m_minParaview = false;


// ----------------------------------------------------------------------------
// Constructor with autonumbering as input parameter
TrilobeCylinder::TrilobeCylinder( bool const& autonumbering )
  : CompositeParticle( autonumbering )
{
  m_specific_composite_shape = "TrilobeCylinder";
  m_radius = 0.;
  m_height = 0.;
  m_armLength = 0.;
}





// ----------------------------------------------------------------------------
// Constructor with an XML node as an input parameter. This constructor is
// expected to be used for reference composite particles
TrilobeCylinder::TrilobeCylinder( DOMNode* root, int const& pc )
  : CompositeParticle( false )
{
  m_specific_composite_shape = "TrilobeCylinder";

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
  m_mass = m_density * ( 6. * m_armLength * m_radius
  	+ ( 3. * PI / 2. - sqrt( 3. ) ) * m_radius * m_radius ) * m_height;	
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
  
  // Paraview
  DOMNode* npar = ReaderXML::getNode( root, "Paraview" );
  if ( npar )
  {
    string spar = ReaderXML::getNodeAttr_String( npar, "Minimize" );
    if ( spar == "True" ) m_minParaview = true;
  }  

  // Angular position of the composite particle
  m_geoRBWC->getTransform()->load( root );
  
  // Moment of inertia tensor of the TrilobeCylinder
  m_inertia[1] = m_inertia[2] = m_inertia[4] = 0.;
  double h2 = m_height * m_height, r2 = m_radius * m_radius,
  	lmr3 = m_armLength - m_radius / sqrt( 3. ),
  	lpr3 = m_armLength + m_radius / sqrt( 3. );
  double a1 = sqrt( 3. ) * r2 * m_height,
  	a2 = 2. * lmr3 * m_radius * m_height,
	a3 = 0.5 * PI * r2 * m_height;
  double Ixpxp = ( a2 * ( 4. * r2 + h2 ) + a3 * ( 3. * r2 + h2 ) ) / 12.,
	Izpzp = a2 * ( lmr3 * lmr3 + h2 + 3. * lpr3 * lpr3 ) / 12.
	+ a3 * ( ( 0.25 - 16. / ( 9. * PI * PI ) ) * r2 + h2 / 12. 
		+ pow( m_armLength + 4. * m_radius / ( 3. * PI ), 2. ) );
  m_inertia[0] = a1 * r2 / 6. + 1.5 * Ixpxp + 1.5 * Izpzp ;
  m_inertia[3] = a1 * h2 / 12. + 3. * ( a2 * ( 4. * r2 + lmr3 * lmr3 
  	+ 3. * lpr3 * lpr3 ) / 12. 
 	+ a3 * ( ( 0.5 - 16. / ( 9. * PI * PI ) ) * r2  
		+ pow( m_armLength + 4. * m_radius / ( 3. * PI ), 2. ) ) );
  m_inertia[5] = m_inertia[0]; 
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

  // Create the 3 spherocylindrical prisms
  double ang[] = { 0., - 2. * PI / 3., 2. * PI / 3.}; 
  for (size_t m=0;m<m_nbElemPart;++m)
  {
    SpheroCylindricalPrism* scp = new SpheroCylindricalPrism( m_radius, 
  	m_armLength, m_height );
    RigidBodyWithCrust* geoRBWC_scp = new RigidBodyWithCrust( scp, Transform(),
  	false, crust_thickness ); 
    m_elementaryParticles[m] = new Particle( geoRBWC_scp, m_density,
      	m_materialName, pc );
    m_InitialRelativePositions[m][X] = 0.5 * m_armLength * cos( ang[m] );
    m_InitialRelativePositions[m][Y] = 0.;  
    m_InitialRelativePositions[m][Z] = - 0.5 * m_armLength * sin( ang[m] );   
    m_elementaryParticles[m]->setPosition( m_InitialRelativePositions[m] );  
    m_elementaryParticles[m]->setMasterParticle( this );
    m_InitialRotationMatrices[m].setValue( 
  	cos( ang[m] ), 0., sin( ang[m] ), 
	0., 1., 0.,
	- sin( ang[m] ), 0., cos( ang[m] ) );  
    m_elementaryParticles[m]->getRigidBody()->getTransform()->setBasis(
  	m_InitialRotationMatrices[m] );
    m_elementaryParticles[m]->getRigidBody()
 	->composeLeftByTransform( *(m_geoRBWC->getTransform()) );
  }

  // Set the the circumscribed radius
  setCircumscribedRadius();
  
  // Compute and set the non-spherical bounding volume
  if ( GrainsExec::m_colDetBoundingVolume ) createBoundingVolume();  

  // In case part of the particle acceleration computed explicity
  if ( Particle::m_splitExplicitAcceleration ) createVelocityInfosNm1();
}





// ----------------------------------------------------------------------------
// Constructor with input parameters
TrilobeCylinder::TrilobeCylinder( int const& id_,
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
  m_specific_composite_shape = "TrilobeCylinder";
  TrilobeCylinder const* TrilobeCylinderRef =
  	dynamic_cast<TrilobeCylinder const*>(ParticleRef);
  m_radius = TrilobeCylinderRef->m_radius; 
  m_armLength = TrilobeCylinderRef->m_armLength;    
  m_height = TrilobeCylinderRef->m_height;
}




// ----------------------------------------------------------------------------
// Constructor with input parameters
TrilobeCylinder::TrilobeCylinder( int const& id_,
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
  m_specific_composite_shape = "TrilobeCylinder";
  TrilobeCylinder const* TrilobeCylinderRef =
  	dynamic_cast<TrilobeCylinder const*>(ParticleRef);
  m_radius = TrilobeCylinderRef->m_radius; 
  m_armLength = TrilobeCylinderRef->m_armLength;    
  m_height = TrilobeCylinderRef->m_height;
}




// ----------------------------------------------------------------------------
// Destructor
TrilobeCylinder::~TrilobeCylinder()
{}




// ----------------------------------------------------------------------------
// Copy constructor (the torsor is initialized to 0)
TrilobeCylinder::TrilobeCylinder( TrilobeCylinder const& other, 
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
Particle* TrilobeCylinder::createCloneCopy( bool const& autonumbering ) const
{
  Particle* particle = new TrilobeCylinder( *this, autonumbering );

  return ( particle );
}




// ----------------------------------------------------------------------------
// Creates a clone of the composite particle. This method calls the
// constructor TrilobeCylinder( int const& id_, Particle const* ParticleRef,
// Vector3 const& vtrans, Quaternion const& qrot, Vector3 const& vrot,
// Transform const& config, ParticleActivity const& activ ) and is used for
// periodic clone composite particles to be inserted in the simulation.
// Autonumbering is set to false and numbering is set with the parameter id_
Particle* TrilobeCylinder::createCloneCopy( int const& id_,
    	Particle const* ParticleRef, Vector3 const& vtrans,
	Quaternion const& qrot,	Vector3 const& vrot,
	Transform const& config, ParticleActivity const& activ,
	map< std::tuple<int,int,int>,
     	std::tuple<bool, Vector3, Vector3, Vector3> > const* contactMap ) const
{
  Particle* particle = new TrilobeCylinder( id_, ParticleRef, vtrans,
	qrot, vrot, config, activ, contactMap );

  return ( particle );
}




// ----------------------------------------------------------------------------
// Returns the number of points to write the composite particle in a
// Paraview format
int TrilobeCylinder::numberOfPoints_PARAVIEW() const
{
  if ( m_minParaview )
    return ( int(m_nbElemPart) * ( 2 * m_visuNodeNbPerHalf + 8 ) ); 
  else
    return ( CompositeParticle::numberOfPoints_PARAVIEW() ); 
}




// ----------------------------------------------------------------------------
// Returns the number of elementary polytopes to write the
// composite particle shape in a Paraview format
int TrilobeCylinder::numberOfCells_PARAVIEW() const
{  
  if ( m_minParaview )
    return ( int(m_nbElemPart) * ( m_visuNodeNbPerHalf + 1 ) );
  else
    return ( CompositeParticle::numberOfCells_PARAVIEW() );       
}




// ----------------------------------------------------------------------------
// Writes the points describing the composite particle in a Paraview format
void TrilobeCylinder::write_polygonsPts_PARAVIEW( ostream& f,
	Vector3 const* translation )const
{ 
  if ( m_minParaview )
  {  
    Point3 pp, p;
    double dtheta = PI / m_visuNodeNbPerHalf, tstartright = - 0.5 * PI;

    for (size_t m=0;m<m_nbElemPart;++m)
    {
      Transform const* transform = m_elementaryParticles[m]->getRigidBody()
  	->getTransform();
 
      // Half cylinder
      // Lower disk rim
      p[Y] = - 0.5 * m_height;
      for (int i=0;i<m_visuNodeNbPerHalf+1;++i)
      {
        p[X] = m_radius * cos ( tstartright + i * dtheta ) + 0.5 * m_armLength;
        p[Z] = m_radius * sin ( tstartright + i * dtheta );
        pp = (*transform)( p );
        if ( translation ) pp += *translation;
        f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
      }

      // Upper disk rim
      p[Y] = 0.5 * m_height;
      for (int i=0;i<m_visuNodeNbPerHalf+1;++i)
      {
        p[X] = m_radius * cos ( tstartright + i * dtheta ) + 0.5 * m_armLength;
        p[Z] = m_radius * sin ( tstartright + i * dtheta );
        pp = (*transform)( p );
        if ( translation ) pp += *translation;
        f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
      }

      // Lower disk center
      p[X] = 0.5 * m_armLength;
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
    
      // Lower corners
      p[X] = - 0.5 * m_armLength;
      p[Y] = - 0.5 * m_height;
      p[Z] = - m_radius;
      pp = (*transform)( p );
      if ( translation ) pp += *translation;
      f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
    
      p[Z] = m_radius;
      pp = (*transform)( p );
      if ( translation ) pp += *translation;
      f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
    
      // Upper corners
      p[X] = - 0.5 * m_armLength;
      p[Y] = 0.5 * m_height;
      p[Z] = - m_radius;
      pp = (*transform)( p );
      if ( translation ) pp += *translation;
      f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;
    
      p[Z] = m_radius;
      pp = (*transform)( p );
      if ( translation ) pp += *translation;
      f << pp[X] << " " << pp[Y] << " " << pp[Z] << endl;                   
    } 
  }
  else
    CompositeParticle::write_polygonsPts_PARAVIEW( f, translation );         
}




// ----------------------------------------------------------------------------
// Returns a list of points describing the component in a Paraview format
list<Point3> TrilobeCylinder::get_polygonsPts_PARAVIEW(
	Vector3 const* translation ) const
{
  list<Point3> ParaviewPoints;
  
  if ( m_minParaview )
  { 
    Point3 pp, p;
    double dtheta = PI / m_visuNodeNbPerHalf, tstartright = - 0.5 * PI;

    for (size_t m=0;m<m_nbElemPart;++m)
    {
      Transform const* transform = m_elementaryParticles[m]->getRigidBody()
  	->getTransform();
 
      // Half cylinder
      // Lower disk rim
      p[Y] = - 0.5 * m_height;
      for (int i=0;i<m_visuNodeNbPerHalf+1;++i)
      {
        p[X] = m_radius * cos ( tstartright + i * dtheta ) + 0.5 * m_armLength;
        p[Z] = m_radius * sin ( tstartright + i * dtheta );
        pp = (*transform)( p );
        if ( translation ) pp += *translation;
        ParaviewPoints.push_back( pp );
      }

      // Upper disk rim
      p[Y] = 0.5 * m_height;
      for (int i=0;i<m_visuNodeNbPerHalf+1;++i)
      {
        p[X] = m_radius * cos ( tstartright + i * dtheta ) + 0.5 * m_armLength;
        p[Z] = m_radius * sin ( tstartright + i * dtheta );
        pp = (*transform)( p );
        if ( translation ) pp += *translation;
        ParaviewPoints.push_back( pp );
      }

      // Lower disk center
      p[X] = 0.5 * m_armLength;
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
    
      // Lower corners
      p[X] = - 0.5 * m_armLength;
      p[Y] = - 0.5 * m_height;
      p[Z] = - m_radius;
      pp = (*transform)( p );
      if ( translation ) pp += *translation;
      ParaviewPoints.push_back( pp );
    
      p[Z] = m_radius;
      pp = (*transform)( p );
      if ( translation ) pp += *translation;
      ParaviewPoints.push_back( pp );
    
      // Upper corners
      p[X] = - 0.5 * m_armLength;
      p[Y] = 0.5 * m_height;
      p[Z] = - m_radius;
      pp = (*transform)( p );
      if ( translation ) pp += *translation;
      ParaviewPoints.push_back( pp );
    
      p[Z] = m_radius;
      pp = (*transform)( p );
      if ( translation ) pp += *translation;
      ParaviewPoints.push_back( pp );
    }
  }
  else
    ParaviewPoints = CompositeParticle::get_polygonsPts_PARAVIEW( translation );
  
  return ( ParaviewPoints );
    
}




// ----------------------------------------------------------------------------
// Writes the composite particle in a Paraview format
void TrilobeCylinder::write_polygonsStr_PARAVIEW( list<int>& connectivity,
    	list<int>& offsets, list<int>& cellstype, int& firstpoint_globalnumber,
	int& last_offset ) const
{
  if ( m_minParaview )
  { 
    for (size_t m=0;m<m_nbElemPart;++m)
    {  
      // Half cylinder 
      for (int i=0;i<m_visuNodeNbPerHalf;++i)
      {
        connectivity.push_back( firstpoint_globalnumber + i );
        connectivity.push_back( firstpoint_globalnumber + i + 1 );
        connectivity.push_back( firstpoint_globalnumber 
      	+ 2 * ( m_visuNodeNbPerHalf + 1 ) );
        connectivity.push_back( firstpoint_globalnumber + i 
      	+ m_visuNodeNbPerHalf + 1 );
        connectivity.push_back( firstpoint_globalnumber + i 
      	+ m_visuNodeNbPerHalf + 2 );
        connectivity.push_back( firstpoint_globalnumber + 
      	2 * ( m_visuNodeNbPerHalf + 1 ) + 1 );
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
      connectivity.push_back( firstpoint_globalnumber + 2 * m_visuNodeNbPerHalf + 
  	5 );
      connectivity.push_back( firstpoint_globalnumber + 2 * m_visuNodeNbPerHalf + 
  	7 );	
      connectivity.push_back( firstpoint_globalnumber + 2 * m_visuNodeNbPerHalf + 
  	6 );	
      last_offset += 8;
      offsets.push_back( last_offset );
      cellstype.push_back( 12 );
  
      firstpoint_globalnumber += 2 * m_visuNodeNbPerHalf + 8;
    }
  }
  else
    CompositeParticle::write_polygonsStr_PARAVIEW( connectivity,
    	offsets, cellstype, firstpoint_globalnumber, last_offset ); 
}






// ----------------------------------------------------------------------------
// Returns the number of corners of the rigib body shape and a code
// describing the rigid body shape
int TrilobeCylinder::getNbCorners() const
{
  return ( 1001 );
}





// ----------------------------------------------------------------------------
// Outputs information to be transferred to the fluid
// Same format data as a regular cylinder: center of bottom circular
// face of the elementary cylinder, an arbitrary point on the lateral surface 
// of the elementary cylinder and center of top circular face of the elementary
// cylinder    
void TrilobeCylinder::writePositionInFluid( ostream& fluid )
{
  // TO DO
}




// ----------------------------------------------------------------------------
// Reads composite particle data from a stream. Usage: for standard composite
// particles in the 2014 reload format
void TrilobeCylinder::read2014( istream& fileIn, vector<Particle*> const*
  	referenceParticles )
{
  CompositeParticle::read2014( fileIn, referenceParticles );

  // We know that (*referenceParticles)[m_GeomType] points to a
  // CompositeParticle, such that we can dynamic cast it to actual type and
  // use -> instead of using get methods through virtual typing
  TrilobeCylinder const* CompParticleRef =
  	dynamic_cast<TrilobeCylinder const*>(
		(*referenceParticles)[m_GeomType]);

  m_radius = CompParticleRef->m_radius;
  m_armLength = CompParticleRef->m_armLength;
  m_height = CompParticleRef->m_height;
}




// ----------------------------------------------------------------------------
// Reads composite particle data from a stream in a binary form.
// Usage: for standard composite particles in the 2014 reload format
void TrilobeCylinder::read2014_binary( istream& fileIn,
	vector<Particle*> const* referenceParticles )
{
  CompositeParticle::read2014_binary( fileIn, referenceParticles );

  // We know that (*referenceParticles)[m_GeomType] points to a
  // CompositeParticle, such that we can dynamic cast it to actual type and
  // use -> instead of using get methods through virtual typing
  TrilobeCylinder const* CompParticleRef =
  	dynamic_cast<TrilobeCylinder const*>(
		(*referenceParticles)[m_GeomType]);

  m_radius = CompParticleRef->m_radius;
  m_armLength = CompParticleRef->m_armLength;
  m_height = CompParticleRef->m_height;
}




// ----------------------------------------------------------------------------
// Computes and sets the circumscribed radius
void TrilobeCylinder::setCircumscribedRadius()
{
  m_geoRBWC->setCircumscribedRadius( sqrt( 0.25 * m_height * m_height 
  	+ ( m_armLength + m_radius ) * ( m_armLength + m_radius ) ) );
}




// ----------------------------------------------------------------------------
// Saves additional features of a (in practice reference) composite particle
// for reload
void TrilobeCylinder::writeAdditionalFeatures( ostream& fileSave ) const
{
  fileSave << endl << "*RadiusArmLengthHeight " << m_radius << " " << 
  	m_armLength << " " << m_height;  
  CompositeParticle::writeAdditionalFeatures( fileSave );
}




// ----------------------------------------------------------------------------
// Reads additional features of a (in practice reference) composite particle
// data from a stream
void TrilobeCylinder::readAdditionalFeatures( istream& fileIn )
{
  string buffer;
  fileIn >> buffer >> m_radius >> m_armLength >> m_height;
  CompositeParticle::readAdditionalFeatures( fileIn );
}
