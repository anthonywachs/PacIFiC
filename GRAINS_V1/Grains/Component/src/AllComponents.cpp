#include "GrainsMPIWrapper.hh"
#include "GrainsExec.hh"
#include "AllComponents.hh"
#include "App.hh"
#include "AppCollision.hh"
#include "RigidBodyWithCrust.hh"
#include "CompositeObstacle.hh"
#include "ObstacleBuilderFactory.hh"
#include "PostProcessingWriter.hh"
#include "ParaviewPostProcessingWriter.hh"
#include "GrainsBuilderFactory.hh"
#include <math.h>
#include <stdlib.h>



// ----------------------------------------------------------------------------
// Default constructor
AllComponents::AllComponents()
  : m_wait( NULL )
  , m_nb_particles( 0 )
  , m_total_nb_particles( 0 )
  , m_nb_active_particles( 0 )
  , m_total_nb_active_particles( 0 )
  , m_nb_active_particles_on_proc( 0 )
  , m_total_nb_active_particles_on_all_procs( 0 )
  , m_nb_inactive_particles( 0 )
  , m_total_nb_physical_particles( 0 )
  , m_obstacle( NULL )
  , m_outputTorsorObstacles_counter( 1 )
  , m_outputTorsorObstacles_frequency( 0 )
{
  m_obstacle = new CompositeObstacle( "__AllObstacles___" );
}




// ----------------------------------------------------------------------------
// Destructor
AllComponents::~AllComponents()
{
  list<Particle*>::iterator particle;

  for (particle=m_ActiveParticles.begin();
  	particle!=m_ActiveParticles.end(); particle++)
    delete *particle;

  for (particle=m_InactiveParticles.begin(); 
  	particle!=m_InactiveParticles.end();particle++)  delete *particle;

  m_ActiveParticles.clear();
  m_InactiveParticles.clear();
  m_PeriodicCloneParticles.clear();

  m_ParticlesInHalozone.clear();
  m_CloneParticles.clear();

  vector<Particle*>::iterator ivp;
  for (ivp=m_ReferenceParticles.begin();
  	ivp!=m_ReferenceParticles.end(); ivp++)
    delete *ivp;
  m_ReferenceParticles.clear();

  delete m_obstacle;

  list<PostProcessingWriter*>::iterator pp;
  for (pp=m_postProcessors.begin();pp!=m_postProcessors.end();pp++)
    delete *pp;
  m_postProcessors.clear();
 
  // Note: the loadings are destroyed by the destructors of the classes
  // ObstacleKinematicsVelocity and ObstacleKinematicsForce
  // Hence we are free the lists here but do not destroy the pointed objects
  m_AllImposedVelocitiesOnObstacles.clear();
  m_AllImposedForcesOnObstacles.clear();
}




// ----------------------------------------------------------------------------
// Updates particle activity
void AllComponents::UpdateParticleActivity()
{
  list<Particle*>::iterator particle;
  for (particle=m_ActiveParticles.begin();
  	particle!=m_ActiveParticles.end(); )
  {
    switch ( (*particle)->getActivity() )
    {
      case COMPUTE :
        particle++;
        break;

      case CLEARandWAIT:
        (*particle)->reset();
        m_InactiveParticles.push_back(*particle);
        particle = m_ActiveParticles.erase( particle );
        break;

      default:
        break;
    }
  }
}




// ----------------------------------------------------------------------------
// Adds a particle
void AllComponents::AddParticle( Particle* particle )
{
  switch ( particle->getActivity() )
  {
    case WAIT:
      m_InactiveParticles.push_back( particle );
      break;

    case COMPUTE:
      m_ActiveParticles.push_back( particle );
      break;

    default:
      break;
  }
}




// ----------------------------------------------------------------------------
// Adds a reference particle
void AllComponents::AddReferenceParticle( Particle* particle )
{
  m_ReferenceParticles.reserve( m_ReferenceParticles.size() + 1 );
  m_ReferenceParticles.push_back( particle );
}




// ----------------------------------------------------------------------------
// Adds an obstacle
void AllComponents::AddObstacle( Obstacle* obstacle_ )
{
  if ( !m_obstacle ) m_obstacle = obstacle_;
  else m_obstacle->append( obstacle_ );
}




// ----------------------------------------------------------------------------
// Associates the imposed velocity to the obstacle
void AllComponents::LinkImposedMotion( ObstacleImposedVelocity* impvel )
{
  m_obstacle->LinkImposedMotion( impvel );
  m_AllImposedVelocitiesOnObstacles.push_back( impvel );
}




// ----------------------------------------------------------------------------
// Associates the imposed velocity to the obstacle
void AllComponents::LinkImposedMotion( ObstacleImposedForce* load )
{
  m_obstacle->LinkImposedMotion( load );
  m_AllImposedForcesOnObstacles.push_back( load );
}




// ----------------------------------------------------------------------------
// Moves all components
list<SimpleObstacle*> AllComponents::Move( double time,
	double const& dt_particle_vel, 
    	double const& dt_particle_disp,
	double const& dt_obstacle )
{
  try{
  // Particles displacement
  list<Particle*>::iterator particle;
  for (particle=m_ActiveParticles.begin();
      particle!=m_ActiveParticles.end(); particle++)
    if ( (*particle)->getTag() != 2 )
      (*particle)->Move( time, dt_particle_vel, dt_particle_disp );

  // Obstacles displacement
  list<SimpleObstacle*> displacedObstacles;
  if ( !m_AllImposedVelocitiesOnObstacles.empty()
  	|| !m_AllImposedForcesOnObstacles.empty() )
  {
    m_obstacle->resetKinematics();
    displacedObstacles = m_obstacle->Move( time, dt_obstacle, false, false );

    list<ObstacleImposedVelocity*>::iterator chargement;
    for (chargement=m_AllImposedVelocitiesOnObstacles.begin();
  	chargement!=m_AllImposedVelocitiesOnObstacles.end(); )
      if ( (*chargement)->isCompleted( time, dt_obstacle ) )
        chargement = m_AllImposedVelocitiesOnObstacles.erase( chargement );
      else chargement++;

    list<ObstacleImposedForce*>::iterator chargement_F;
    for (chargement_F=m_AllImposedForcesOnObstacles.begin();
  	chargement_F!=m_AllImposedForcesOnObstacles.end(); )
    {
      if ( (*chargement_F)->isCompleted( time, dt_obstacle ) )
        chargement_F = m_AllImposedForcesOnObstacles.erase( chargement_F );
      else chargement_F++;
    }
  }

  return ( displacedObstacles );
  }
  catch (const DisplacementError&) {
    throw DisplacementError();
  }
}




// ----------------------------------------------------------------------------
// Computes particles acceleration
void AllComponents::computeParticlesAcceleration( double time )
{
  list<Particle*>::iterator particle;
  for (particle=m_ActiveParticles.begin();
      particle!=m_ActiveParticles.end(); particle++)
    if ( (*particle)->getTag() != 2 )
      (*particle)->computeAcceleration( time );
}





// ----------------------------------------------------------------------------
// Advances particles velocity over dt_particle_vel
void AllComponents::advanceParticlesVelocity( double time, 
    	double const& dt_particle_vel )
{
  list<Particle*>::iterator particle;
  for (particle=m_ActiveParticles.begin();
      particle!=m_ActiveParticles.end(); particle++)
    if ( (*particle)->getTag() != 2 )
      (*particle)->advanceVelocity( time, dt_particle_vel );
}




// ----------------------------------------------------------------------------
// Updates obstacles' velocity without actually moving them
void AllComponents::setKinematicsObstacleWithoutMoving(
	double time, double dt )
{
  bool bbb = Obstacle::getMoveObstacle() ;

  Obstacle::setMoveObstacle( false ) ;

  if ( !m_AllImposedVelocitiesOnObstacles.empty() )
  {
    m_obstacle->resetKinematics();
    m_obstacle->Move( time, dt, false, false );

    list<ObstacleImposedVelocity*>::iterator il;
    for (il=m_AllImposedVelocitiesOnObstacles.begin();
  	il!=m_AllImposedVelocitiesOnObstacles.end(); )
      if ( (*il)->isCompleted( time, dt ) )
        il = m_AllImposedVelocitiesOnObstacles.erase( il );
      else il++;
  }

  Obstacle::setMoveObstacle( bbb ) ;
}





// ----------------------------------------------------------------------------
// Initializes forces exerted on all components and set coordination number
// to 0
void AllComponents::InitializeForces( double time, double dt,
	bool const& withWeight )
{
  list<Particle*>::iterator particle;

  // Particles actives
  for (particle=m_ActiveParticles.begin();
  	particle!=m_ActiveParticles.end(); particle++)
    (*particle)->InitializeForce( withWeight );

  // Obstacles
  m_obstacle->InitializeForce( false );
}




// ----------------------------------------------------------------------------
// Initializes the transformation with crust of all components to not computed
void AllComponents::InitializeRBTransformWithCrustState( double time,
	double dt )
{
  list<Particle*>::iterator particle;

  // Actives particles
  for (particle=m_ActiveParticles.begin();
  	particle!=m_ActiveParticles.end(); particle++)
    (*particle)->initialize_transformWithCrust_to_notComputed();

  // Obstacles
  list<SimpleObstacle*> list_obstacles = m_obstacle->getObstacles();
  list<SimpleObstacle*>::iterator myObs;
  for (myObs=list_obstacles.begin();myObs!=list_obstacles.end();myObs++)
    (*myObs)->getRigidBody()->initialize_transformWithCrust_to_notComputed();
}




// ----------------------------------------------------------------------------
// Compute all particles weight
void AllComponents::computeWeight( double time, double dt )
{
  list<Particle*>::iterator particle;
  vector<Particle*>::iterator ivp;

  // Classes de reference
  for (ivp=m_ReferenceParticles.begin();
  	ivp!=m_ReferenceParticles.end(); ivp++)
    (*ivp)->computeWeight();

  // Particles en attente
  for (particle=m_InactiveParticles.begin();
  	particle!=m_InactiveParticles.end(); particle++)
    (*particle)->computeWeight();

  // Particles actives
  for (particle=m_ActiveParticles.begin();
  	particle!=m_ActiveParticles.end(); particle++)
    (*particle)->computeWeight();
}




// ----------------------------------------------------------------------------
// Returns an obstacle using its name as an input
Obstacle const* AllComponents::getObstacle( string const& name ) const
{
  return ( m_obstacle->getObstacleFromName( name ) );
}




// ----------------------------------------------------------------------------
// Returns the root obstacle
Obstacle* AllComponents::getObstacles()
{
  return ( m_obstacle );
}




// ----------------------------------------------------------------------------
// Returns particle with ID number id
Particle* AllComponents::getParticle( int id )
{
  Particle *particle = NULL;
  list<Particle*>::iterator iter;
  bool found = false;
  for (iter=m_ActiveParticles.begin();iter!=m_ActiveParticles.end()
  	&& !found;iter++)
    if ( (*iter)->getID() == id )
    {
      particle = *iter;
      found = true;
    }

  return ( particle );
}




// ----------------------------------------------------------------------------
// Returns component with ID number id
Component* AllComponents::getComponent( int id )
{
  Component *composant = getParticle( id );
  if ( composant == NULL )
  {
    bool found = false;
    list<SimpleObstacle*> obstacles = m_obstacle->getObstacles();
    list<SimpleObstacle*>::iterator myObs;
    for (myObs=obstacles.begin(); myObs!=obstacles.end() && !found; myObs++)
      if ( (*myObs)->getID() == id )
      {
        composant = *myObs;
        found = true;
      }
  }

  return ( composant );
}




// ----------------------------------------------------------------------------
// Returns a particle from the list of inactive particles
Particle* AllComponents::getParticle( PullMode mode,
	GrainsMPIWrapper const* wrapper )
{
  if ( !m_InactiveParticles.empty() )
  {
    switch ( mode )
    {
      case PM_ORDERED:
        m_wait = m_InactiveParticles.front();
        break;

      case PM_RANDOM:
        if ( m_wait == NULL )
	{
	  double v = double(random()) / double(INT_MAX);
	  int id = int( double(m_InactiveParticles.size()) * v );

	  // Only the master proc randomly selects and then the random
	  // number is broadcasted to all other processes
	  if ( wrapper ) id = wrapper->Broadcast_INT( id );

	  list<Particle*>::iterator p = m_InactiveParticles.begin();
	  for (int i=0; i<id && p!=m_InactiveParticles.end(); i++, p++) {}
	  m_wait = *p;
        }
        break;
    }
  }
  else m_wait = NULL;

  return ( m_wait );
}




// ----------------------------------------------------------------------------
// Returns the list of active particles
list<Particle*>* AllComponents::getActiveParticles()
{
  return ( &m_ActiveParticles );
}




// ----------------------------------------------------------------------------
// Returns a const pointer to the list of active particles
list<Particle*> const* AllComponents::getActiveParticles() const
{
  return ( &m_ActiveParticles );
}




// ----------------------------------------------------------------------------
// Returns a pointer to the list of inactive particles
list<Particle*>* AllComponents::getInactiveParticles()
{
  return ( &m_InactiveParticles );
}




// ----------------------------------------------------------------------------
// Returns a const pointer to the list of inactive particles
list<Particle*> const* AllComponents::getInactiveParticles() const
{
  return ( &m_InactiveParticles );
}




// ----------------------------------------------------------------------------
// Returns a pointer to the list of particles in the halozone
list<Particle*>* AllComponents::getParticlesInHalozone()
{
  return ( &m_ParticlesInHalozone );
}




// ----------------------------------------------------------------------------
// Returns a pointer to the list of particles in the halozone
list<Particle*> const* AllComponents::getParticlesInHalozone() const
{
  return ( &m_ParticlesInHalozone );
}




// ----------------------------------------------------------------------------
// Returns a pointer to the list of clone particles
list<Particle*>* AllComponents::getCloneParticles()
{
  return ( &m_CloneParticles );
}




// ----------------------------------------------------------------------------
// Returns a pointer to the map of serial clone particles */
multimap<int,Particle*>* AllComponents::getPeriodicCloneParticles()
{
  return ( &m_PeriodicCloneParticles );
}




// ----------------------------------------------------------------------------
// Returns a pointer to the map of serial clone particles */
multimap<int,Particle*> const* AllComponents::getPeriodicCloneParticles() const
{
  return ( &m_PeriodicCloneParticles );
}




// ----------------------------------------------------------------------------
// Returns a pointer to the vector of reference particles
vector<Particle*>* AllComponents::getReferenceParticles()
{
  return ( &m_ReferenceParticles );
}




// ----------------------------------------------------------------------------
// Returns a pointer to the vector of reference particles
vector<Particle*> const* AllComponents::getReferenceParticles() const
{
  return ( &m_ReferenceParticles );
}




// ----------------------------------------------------------------------------
// Returns the list of obstacles to send to the fluid
list<Obstacle*> AllComponents::getObstaclesToFluid() const
{
  return ( m_obstacle->getObstaclesToFluid() );
}




// ----------------------------------------------------------------------------
// Returns the maximum particle circumscribed radius
double AllComponents::getCircumscribedRadiusMax()
{
  double radiusMax = 0.0;
  double radius;

  vector<Particle*>::iterator particle;
  for (particle=m_ReferenceParticles.begin();
  	particle!=m_ReferenceParticles.end(); particle++)
  {
    radius = (*particle)->getCircumscribedRadius();
    radiusMax = max( radiusMax, radius );
  }

  return ( radiusMax );
}




// ----------------------------------------------------------------------------
// Returns the minimum particle circumscribed radius
double AllComponents::getCircumscribedRadiusMin()
{
  double radiusMin = 1.e10;
  double radius;

  vector<Particle*>::iterator particle;
  for (particle=m_ReferenceParticles.begin();
  	particle!=m_ReferenceParticles.end(); particle++)
  {
    radius = (*particle)->getCircumscribedRadius();
    radiusMin = radiusMin > radius ? radius : radiusMin;
  }

  return ( radiusMin );
}




// ----------------------------------------------------------------------------
// Returns the maximum particle crust thickness
double AllComponents::getCrustThicknessMax()
{
  double crustThicknessMax = 0.0;
  double crustThickness;

  vector<Particle*>::iterator particle;
  for (particle=m_ReferenceParticles.begin();
  	particle!=m_ReferenceParticles.end(); particle++)
  {
    crustThickness = (*particle)->getCrustThickness();
    crustThicknessMax = crustThicknessMax > crustThickness ?
    	crustThicknessMax : crustThickness;
  }

  return ( crustThicknessMax );
}




// ----------------------------------------------------------------------------
// Returns the minimum particle crust thickness
double AllComponents::getCrustThicknessMin()
{
  double crustThicknessMin = 1.e10;
  double crustThickness;

  vector<Particle*>::iterator particle;
  for (particle=m_ReferenceParticles.begin();
  	particle!=m_ReferenceParticles.end(); particle++)
  {
    crustThickness = (*particle)->getCrustThickness();
    crustThicknessMin = crustThicknessMin > crustThickness ?
    	crustThickness : crustThicknessMin;
  }

  return ( crustThicknessMin );
}




// ----------------------------------------------------------------------------
// Returns the cumulative volume of all particles, both active and inactive
double AllComponents::getVolume() const
{
  double volume = 0.;
  list<Particle*>::const_iterator particle;

  for (particle=m_ActiveParticles.begin();
  	particle!=m_ActiveParticles.end(); particle++)
    volume += (*particle)->getVolume();

  for (particle=m_InactiveParticles.begin();
  	particle!=m_InactiveParticles.end(); particle++)
    volume += (*particle)->getVolume();

  return ( volume );
}




// ----------------------------------------------------------------------------
// Returns the cumulative volume of all actives particles
double AllComponents::getVolumeIn() const
{
  double volume = 0.;
  list<Particle*>::const_iterator particle;

  for (particle=m_ActiveParticles.begin();
  	particle!=m_ActiveParticles.end(); particle++)
    volume += (*particle)->getVolume();

  return ( volume );
}




// ----------------------------------------------------------------------------
// Returns the cumulative volume of all inactives particles
double AllComponents::getVolumeOut() const
{
  double volume = 0.;
  list<Particle*>::const_iterator particle;

  for (particle=m_InactiveParticles.begin();
  	particle!=m_InactiveParticles.end(); particle++)
    volume += (*particle)->getVolume();

  return ( volume );
}




// ----------------------------------------------------------------------------
// Links all active particles to the application
void AllComponents::Link( AppCollision& app )
{
  list<Particle*>::iterator particle;
  for (particle=m_ActiveParticles.begin();
  	particle!=m_ActiveParticles.end(); particle++)
    app.Link( *particle );
}




// ----------------------------------------------------------------------------
// Set all active particles velocity to 0 if reset == "Reset"
void AllComponents::resetKinematics( string const& reset )
{
  if ( reset == "Reset" )
  {
    list<Particle*>::iterator particle;
    for (particle=m_ActiveParticles.begin();
    	particle!=m_ActiveParticles.end(); particle++)
      (*particle)->resetKinematics();
  }
}




// ----------------------------------------------------------------------------
// Transfer the inactive particle waiting to be inserted to the list
// of active particles
void AllComponents::ShiftParticleOutIn( bool const& parallel )
{
  m_wait->setActivity( COMPUTE );
  removeParticleFromList( m_InactiveParticles, m_wait );
  m_ActiveParticles.push_back( m_wait );
  if ( parallel )
  {
    if ( m_wait->getTag() == 1 ) m_ParticlesInHalozone.push_back( m_wait );
    else if ( m_wait->getTag() == 2 ) m_CloneParticles.push_back( m_wait );
  }
  m_wait = NULL;
}




// ----------------------------------------------------------------------------
// Removes the inactive particle waiting to be inserted from the
// list of inactive particles and destroys it
void AllComponents::DeleteAndDestroyWait()
{
  removeParticleFromList( m_InactiveParticles, m_wait );
  delete m_wait;
  m_wait = NULL;  
}




// ----------------------------------------------------------------------------
// Delete the fist instance of a pointer to a particle in a list of
// pointers to particle
bool removeParticleFromList( list<Particle*>& pointerslist, Particle* value )
{
  list<Particle*>::iterator particle;
  bool found =false;
  for (particle=pointerslist.begin();particle!=pointerslist.end() && !found; )
    if ( *particle == value )
    {
      particle = pointerslist.erase( particle );
      found = true;
    }
    else particle++;

  return ( found );
}




// ----------------------------------------------------------------------------
// Delete the fist instance of a pointer to a particle in a set of
// pointers to particle
bool removeParticleFromSet( set<Particle*> &pointersSet, Particle* value )
{
  set<Particle*>::iterator particle;
  bool found =false;
  for (particle=pointersSet.begin();particle!=pointersSet.end() && !found; )
    if ( *particle == value )
    {
      pointersSet.erase( particle );
      found = true;
    }
    else particle++;

  return ( found );
}




// ----------------------------------------------------------------------------
// Delete the fist instance of a pointer to a simple obstacle in a list
// of pointers to simple obstacle
bool removeObstacleFromList( list<SimpleObstacle*>& pointerslist,
	SimpleObstacle* value )
{
  list<SimpleObstacle*>::iterator obs;
  bool found =false;
  for (obs=pointerslist.begin();obs!=pointerslist.end() && !found; )
    if (*obs == value)
    {
      obs = pointerslist.erase( obs );
      found = true;
    }
    else obs++;

  return ( found );
}




// ---------------------------------------------------------------------------
// Reloads components from an input stream
void AllComponents::read( istream& fileSave, string const& filename )
{
  string buffer, readingMode ;
  int nbreParticles_, nbParticleTypes_, ParticleTag, ParticleGeomType;
  Particle *particle;

  // Number of particle types
  fileSave >> buffer >> nbParticleTypes_;

  // Reading the reference particles
  m_ReferenceParticles.reserve( nbParticleTypes_ );
  for (int i=0; i<nbParticleTypes_; i++)
  {
    // Read the buffer "<Particle>" or "<CompositeParticle>"
    fileSave >> buffer;

    // Construct an empty particle
    if ( buffer == "<Particle>" )
      particle = new Particle( false );
    else
    {
      fileSave >> buffer >> buffer;
      if ( buffer == "SpheroCylinder" )
        particle = new SpheroCylinder( false );
      else        
        particle = new CompositeParticle( false );
    }

    // Read from stream
    particle->read( fileSave );

    // Add to the vector of reference particles
    m_ReferenceParticles.push_back( particle );

    // Read the buffer "</Particle>" or "</CompositeParticle>"
    fileSave >> buffer;
  }

  // Number of particles in to read
  fileSave >> readingMode >> nbreParticles_;
  if ( readingMode == "Hybrid" ) GrainsExec::m_writingModeHybrid = true ;

  // Reload of particles using the reference particles
  ifstream FILEbin;
  if ( GrainsExec::m_writingModeHybrid )
  {
    string binary_filename = filename + ".bin";
    FILEbin.open( binary_filename.c_str(), ios::in | ios::binary );
  }

  for (int i=0; i<nbreParticles_; i++)
  {
    // Read the geometric type of the particle
    if ( GrainsExec::m_writingModeHybrid )
      FILEbin.read( reinterpret_cast<char*>( &ParticleGeomType ), sizeof(int) );
    else
      fileSave >> ParticleGeomType;

    // Particle construction
    if ( m_ReferenceParticles[ParticleGeomType]->isCompositeParticle() )
    {    
      if ( m_ReferenceParticles[ParticleGeomType]
      		->getSpecificCompositeShapeName() == "SpheroCylinder" )
        particle = new SpheroCylinder( false );
      else   
        particle = new CompositeParticle( false );
    }
    else
      particle = new Particle( false );

    // Set the geometric type of the particle
    particle->setGeometricType( ParticleGeomType );

    // Read the particle features
    if ( GrainsExec::m_writingModeHybrid )
      particle->read2014_binary( FILEbin, &m_ReferenceParticles );
    else
      particle->read2014( fileSave, &m_ReferenceParticles );

    // Add to lists
    switch ( particle->getActivity() )
    {
      case COMPUTE:
        m_ActiveParticles.push_back( particle );
        ParticleTag = particle->getTag();

	// MPI mode
	if ( GrainsExec::m_MPI )
	  switch ( ParticleTag )
          {
            case 1:
              m_ParticlesInHalozone.push_back( particle );
              break;

            case 2:
              m_CloneParticles.push_back( particle );
              break;

	    default:
              break;
          }
	// Serial mode: periodicity
	else
	  if ( ParticleTag == 2 )
	    m_PeriodicCloneParticles.insert(
	    	pair<int,Particle*>( particle->getID(), particle ) );
        break;

      default:
        m_InactiveParticles.push_back( particle );
        break;
    }
  }
  if ( GrainsExec::m_writingModeHybrid ) FILEbin.close();

  // Reload the obstacle tree
  string name;
  fileSave >> buffer;
  fileSave >> buffer >> buffer;
  fileSave >> buffer;
  while ( buffer != "</Composite>" )
  {
    ObstacleBuilderFactory::reload( buffer, *m_obstacle, fileSave );
    fileSave >> buffer;
  }
  fileSave >> buffer;

  assert( buffer == "</Obstacle>" );
}




// ---------------------------------------------------------------------------
// Writes components to an output stream
void AllComponents::write( ostream &fileSave, string const& filename ) const
{
  fileSave << endl << "NumberOfParticleTypes\t" <<
  	m_ReferenceParticles.size() << endl;

  vector<Particle*>::const_iterator iv;
  for (iv=m_ReferenceParticles.begin();
  	iv!=m_ReferenceParticles.end();iv++) (*iv)->write( fileSave );

  size_t nbActivesNonPer = m_ActiveParticles.size();
  list<Particle*>::const_iterator particle;

  fileSave << endl << ( GrainsExec::m_writingModeHybrid ? "Hybrid" :
  	"Text" ) << " " << nbActivesNonPer + m_InactiveParticles.size()
	<< endl;

  if ( GrainsExec::m_writingModeHybrid )
  {
    string binary_filename = filename + ".bin";
    ofstream FILEbin( binary_filename.c_str(), ios::out | ios::binary );

    for (particle=m_ActiveParticles.begin();
  	particle!=m_ActiveParticles.end(); particle++)
      (*particle)->write2014_binary( FILEbin );

    for (particle=m_InactiveParticles.begin();
    	particle!=m_InactiveParticles.end(); particle++)
      (*particle)->write2014_binary( FILEbin );

    FILEbin.close();
  }
  else
  {
    for (particle=m_ActiveParticles.begin();
  	particle!=m_ActiveParticles.end(); particle++)
      (*particle)->write2014( fileSave );

    for (particle=m_InactiveParticles.begin();
    	particle!=m_InactiveParticles.end(); particle++)
      (*particle)->write2014( fileSave );
  }

  fileSave << endl << "<Obstacle>" << endl;
  m_obstacle->write( fileSave );
  fileSave << endl << "</Obstacle>" << endl;
}




// ----------------------------------------------------------------------------
// Debugging method
void AllComponents::debug( char* s )
{
  cout << s << '\n'
       << "Particles " << m_ActiveParticles.size()
       		+ m_InactiveParticles.size()  << '\n'
       << "   Actives " << m_ActiveParticles.size() << '\t'
       << "   Wait    " << m_InactiveParticles.size()      << endl;
}




// ----------------------------------------------------------------------------
// Output operator
ostream& operator << ( ostream& f, AllComponents const& EC )
{
  f << "Number of particles on all processes = "
  	<< EC.m_total_nb_particles << endl;
  f << "Number of particles = " <<
  	EC.m_ActiveParticles.size() + EC.m_InactiveParticles.size() << endl;
  f << "Number of active particles = " << EC.m_ActiveParticles.size()
  	<< endl;
  f << "Number of inactive particles = " << EC.m_InactiveParticles.size()
  	<< endl;
  f << "Number of particles in halozone = " <<
  	EC.m_ParticlesInHalozone.size() << endl;
  f << "Number of clone particles = " << EC.m_CloneParticles.size() << endl;
  list<Particle*>::const_iterator il;
  for (il=EC.m_ActiveParticles.begin();il!=EC.m_ActiveParticles.end();il++)
    f << *(*il) << endl;

  return f;
}




// ----------------------------------------------------------------------------
// Writes components for Post-Processing at the start of the simulation
void AllComponents::PostProcessing_start( double time, double dt,
	LinkedCell const* LC, vector<Window> const& insert_windows,
	int rank, int nprocs,
	GrainsMPIWrapper const* wrapper )
{
  list<Particle*>* postProcessingPeriodic = NULL;
  list<PostProcessingWriter*>::iterator pp;
  bool written = false;

  // Message on display
  if ( rank == 0 )
    for (pp=m_postProcessors.begin();pp!=m_postProcessors.end() && !written;
    	pp++)
    {
      cout << "Writing results in post-processing files: START" << endl 
      	<< std::flush;
      written = true;
    }

  // Periodic clones
  if ( GrainsExec::m_periodic )
  {
    postProcessingPeriodic = new list<Particle*>;
    if ( GrainsExec::m_MPI ) 
    {
      list<Particle*>::const_iterator particle;
      for (particle=m_CloneParticles.begin();
  	particle!=m_CloneParticles.end(); particle++)
        if ( !LC->isInDomain( (*particle)->getPosition() ) )
	  postProcessingPeriodic->push_back( *particle );
    }
    else
    {
      multimap<int,Particle*>::iterator imm;
      for (imm=m_PeriodicCloneParticles.begin();
  	imm!=m_PeriodicCloneParticles.end();imm++)
        postProcessingPeriodic->push_back( imm->second );
    }
  }

  // Post processing writers
  for (pp=m_postProcessors.begin();pp!=m_postProcessors.end();pp++)
    (*pp)->PostProcessing_start( time, dt, &m_ActiveParticles,
	&m_InactiveParticles, postProcessingPeriodic,
	&m_ReferenceParticles, m_obstacle, LC, insert_windows );

  // Destruction of local containers
  if ( GrainsExec::m_periodic )
  {
    postProcessingPeriodic->clear();
    delete postProcessingPeriodic;
  }

  // Message on display
  written = false ;
  if ( rank == 0 )
    for (pp=m_postProcessors.begin();pp!=m_postProcessors.end() && !written;
    	pp++)
    {
      cout << "Writing results in post-processing files: COMPLETED" << endl
      	<< std::flush;
      written = true;
    }
}




// ----------------------------------------------------------------------------
// Writes components for Post-Processing over the simulation
void AllComponents::PostProcessing( double time, double dt,
	LinkedCell const* LC, int rank,
	int nprocs, GrainsMPIWrapper const* wrapper )
{
  list<Particle*>* postProcessingPeriodic = NULL;
  list<PostProcessingWriter*>::iterator pp;
  bool written = false;

  // Message on display
  if ( rank == 0 )
    for (pp=m_postProcessors.begin();pp!=m_postProcessors.end() && !written;
    	pp++)
    {
      cout << "Writing results in post-processing files: START" << endl;
      written = true;
    }

  // Periodic clones
  if ( GrainsExec::m_periodic )
  {
    postProcessingPeriodic = new list<Particle*>;
    if ( GrainsExec::m_MPI )
    {
      list<Particle*>::const_iterator particle;
      for (particle=m_CloneParticles.begin();
  	particle!=m_CloneParticles.end(); particle++)
        if ( !LC->isInDomain( (*particle)->getPosition() ) )
	  postProcessingPeriodic->push_back( *particle );
    }    
    else
    {
      multimap<int,Particle*>::iterator imm;
      for (imm=m_PeriodicCloneParticles.begin();
  	imm!=m_PeriodicCloneParticles.end();imm++)
        postProcessingPeriodic->push_back( imm->second );
    }
  }

  // Post processing writers
  for (pp=m_postProcessors.begin();pp!=m_postProcessors.end();pp++)
    (*pp)->PostProcessing( time, dt, &m_ActiveParticles,
	&m_InactiveParticles, postProcessingPeriodic,
	&m_ReferenceParticles, m_obstacle, LC );

  // Destruction of local containers
  if ( GrainsExec::m_periodic )
  {
    postProcessingPeriodic->clear();
    delete postProcessingPeriodic;
  }

  // Message on display
  written = false ;
  if ( rank == 0 )
    for (pp=m_postProcessors.begin();pp!=m_postProcessors.end() && !written;
    	pp++)
    {
      cout << "Writing results in post-processing files: COMPLETED" << endl;
      written = true;
    }
}




// ----------------------------------------------------------------------------
// Finalizes Post-Processing at the end of the simulation
void AllComponents::PostProcessing_end()
{
  // Post processing writers
  list<PostProcessingWriter*>::iterator pp;
  for (pp=m_postProcessors.begin();pp!=m_postProcessors.end();pp++)
    (*pp)->PostProcessing_end();
}




// ----------------------------------------------------------------------------
// Writes components for Post-Processing in case of an error in
// contact or displacement
void AllComponents::PostProcessingErreurComponents( string const& filename,
	list<Component*> const& errcomposants )
{
  // Post processing writers
  list<PostProcessingWriter*>::iterator pp;
  for (pp=m_postProcessors.begin();pp!=m_postProcessors.end();pp++)
    (*pp)->writeErreurComponentsPostProcessing( filename, errcomposants );
}




// ----------------------------------------------------------------------------
// Computes and returns the maximum and mean translational velocity
// of all components i.e. particles and obstacles
void AllComponents::ComputeMaxMeanVelocity( double& vmax, double& vmean,
  	GrainsMPIWrapper const* wrapper ) const
{
  double vit;
  vmax = vmean = 0. ;
  size_t ncomp = 0 ;

  // Particles
  list<Particle*>::const_iterator particle;
  for (particle=m_ActiveParticles.begin();
  	particle!=m_ActiveParticles.end(); particle++)
  {
    vit = Norm( *(*particle)->getTranslationalVelocity() );
    vmax = vit > vmax ? vit : vmax;
    vmean += vit;
    ++ncomp;
  }

  // Obstacles
  list<SimpleObstacle*> obstacles = m_obstacle->getObstacles();
  list<SimpleObstacle*>::const_iterator myObs;
  for (myObs=obstacles.begin(); myObs!=obstacles.end(); myObs++)
  {
    if ( (*myObs)->hasMoved() )
    {
      vit = Norm( *(*myObs)->getTranslationalVelocity() );
      vmax = vit > vmax ? vit : vmax;
      vmean += vit;
      ++ncomp;
    }
  }

  if ( wrapper )
  {
    vmax = wrapper->max_DOUBLE_master( vmax );
    ncomp = wrapper->sum_UNSIGNED_INT_master( ncomp );
    vmean = wrapper->sum_DOUBLE_master( vmean );
  }

  if ( ncomp ) vmean /= double(ncomp) ;
}




// ----------------------------------------------------------------------------
// Computes and writes in a file the minimum, maximum and mean
// translational velocity of all particles
void AllComponents::monitorParticlesVelocity( double time, ofstream& fileOut,
  	int rank, GrainsMPIWrapper const* wrapper ) const
{
  Vector3 vmin( 1.20 ), vmax( -1.e20 ), vmean;
  size_t ncomp = m_ActiveParticles.size() ;
  Vector3 const* vtrans = NULL ;
  double vit = 0. ;

  // Particles
  list<Particle*>::const_iterator particle;
  for (particle=m_ActiveParticles.begin();
  	particle!=m_ActiveParticles.end(); particle++)
  {
    vtrans = (*particle)->getTranslationalVelocity();
    for (int i=0;i<3;++i)
    {
      vit = (*vtrans)[i];
      vmin[i] = vmin[i] < vit ? vmin[i] : vit;
      vmax[i] = vmax[i] > vit ? vmax[i] : vit;
      vmean[i] += vit;
    }
  }

  if ( wrapper )
  {
    ncomp = wrapper->sum_UNSIGNED_INT_master( ncomp );
    for (int i=0;i<3;++i)
    {
      vmin[i] = wrapper->min_DOUBLE_master( vmin[i] );
      vmax[i] = wrapper->max_DOUBLE_master( vmax[i] );
      vmean[i] = wrapper->sum_DOUBLE_master( vmean[i] );
    }
  }

  if ( rank == 0 )
  {
    fileOut << time;
    for (int i=0;i<3;++i)
    {
      vmean[i] /= double(ncomp) ;
      fileOut << " " << vmin[i] << " " << vmax[i] << " " << vmean[i];
    }
    fileOut << endl;
  }
}




// ----------------------------------------------------------------------------
// Updates geographic position of particles in the halozone
void AllComponents::updateGeoPositionParticlesHalozone()
{
  list<Particle*>::iterator particle;
  for (particle=m_ParticlesInHalozone.begin();
  	particle!=m_ParticlesInHalozone.end();particle++)
    (*particle)->updateGeoPosition();
}




// ----------------------------------------------------------------------------
// Adds a post-processing writer
void AllComponents::addPostProcessingWriter( PostProcessingWriter* ppw )
{
  m_postProcessors.push_back(ppw);
}




// ----------------------------------------------------------------------------
// Sets the initial post-processing cycle number
void AllComponents::setInitialCycleNumber( int const& cycle0)
{
  list<PostProcessingWriter*>::iterator pp;
  for (pp=m_postProcessors.begin();pp!=m_postProcessors.end();pp++)
    (*pp)->setInitialCycleNumber( cycle0 );
}




// ----------------------------------------------------------------------------
// Checks that the Paraview post-processing writer exists, and if not
// creates it
void AllComponents::checkParaviewPostProcessing( int const& rank,
  	int const& nprocs,
  	string const& name_,
	string const& root_,
  	const bool& isBinary )
{
  PostProcessingWriter *ParaviewPP = NULL;
  list<PostProcessingWriter*>::iterator pp;
  for (pp=m_postProcessors.begin();pp!=m_postProcessors.end();pp++)
    if ( (*pp)->getPostProcessingWriterType() == "Paraview" )
      ParaviewPP = *pp;

  // If no post-processing writer, all processors write data  
  if ( m_postProcessors.empty() )
    PostProcessingWriter::allocate_PostProcessingWindow( nprocs );

  // If a Paraview post-processing writer already exists, it is destroyed and
  // recreated with the proper parameters
  if ( ParaviewPP )
  {
    for (pp=m_postProcessors.begin();pp!=m_postProcessors.end();pp++)
      if ( (*pp) == ParaviewPP )
      {
        delete *pp;
	pp = m_postProcessors.erase( pp );
      }
      else pp++;
    ParaviewPP = NULL;
  }

  // Creation of the Paraview post-processing writer with the proper parameters
  if ( rank == 0 )
  {
    cout << GrainsExec::m_shift3 << "Creation of the Paraview"
    	" post-processing writer with the proper parameters" << endl;
    cout << GrainsExec::m_shift6 << "Postprocessing" << endl;	
  }  
  ParaviewPP = new ParaviewPostProcessingWriter( rank, nprocs, name_, root_,
  	isBinary );
  m_postProcessors.push_back( ParaviewPP );
}




// ----------------------------------------------------------------------------
// Returns the highest particle ID number
int AllComponents::getMaxParticleIDnumber() const
{
  int numeroMax = - 1;
  list<Particle*>::const_iterator particle;

  for (particle=m_ActiveParticles.begin();
  	particle!=m_ActiveParticles.end(); particle++)
    numeroMax = numeroMax < (*particle)->getID() ?
    	(*particle)->getID() : numeroMax;

  for (particle=m_InactiveParticles.begin();
  	particle!=m_InactiveParticles.end(); particle++)
    numeroMax = numeroMax < (*particle)->getID() ?
    	(*particle)->getID() : numeroMax;

  return ( numeroMax );
}




// ----------------------------------------------------------------------------
// Sets the frequency at which the relationship between obstacles
// and linked cell grid is updated
void AllComponents::setObstaclesLinkedCellUpdateFrequency(
	int const& updateFreq )
{
  list<SimpleObstacle*> obstacles = m_obstacle->getObstacles();
  list<SimpleObstacle*>::iterator myObs;
  for (myObs=obstacles.begin(); myObs!=obstacles.end(); myObs++)
    (*myObs)->setObstacleLinkedCellUpdateFrequency( updateFreq );
}




// ----------------------------------------------------------------------------
// Sets the parameters to output load exerted on obstacles
void AllComponents::setOutputObstaclesLoadParameters( string const& root_,
  	int const& freq_,
	list<string> const& ObsNames )
{
  m_outputTorsorObstacles_dir = root_;
  m_outputTorsorObstacles_frequency = freq_;

  for (list<string>::const_iterator il=ObsNames.begin();il!=ObsNames.end();il++)
  {
    Obstacle* pobs = const_cast<Obstacle*>(m_obstacle->getObstacleFromName(
    	*il ));
    if ( pobs )
      m_outputTorsorObstacles.push_back( pobs );
  }
}




// ----------------------------------------------------------------------------
// Writes load on obstacles in a file
void AllComponents::outputObstaclesLoad( double time, double dt,
	bool enforceOutput, bool increaseCounterOnly, int rank, int nprocs,
	GrainsMPIWrapper const* wrapper )
{
  if ( !increaseCounterOnly )
  {
    if ( m_outputTorsorObstacles_counter == 0 || enforceOutput )
    {
      if ( nprocs > 1 )
        wrapper->sumObstaclesLoad( m_obstacle->getObstacles() );

      if ( rank == 0 )
      {
        Torsor const* torseur = NULL;
        Vector3 const* force = NULL;
        Vector3 const* torque = NULL;
        for (list<Obstacle*>::iterator
    		obstacle=m_outputTorsorObstacles.begin();
		obstacle!=m_outputTorsorObstacles.end();obstacle++)
        {
          ofstream OUT( ( m_outputTorsorObstacles_dir
      	+ "/Loading_" + (*obstacle)->getName() + ".res" ).c_str(), ios::app );
          torseur = (*obstacle)->getTorsor();
	  force = torseur->getForce();
          torque = torseur->getTorque();
	  OUT << time << " " <<
		GrainsExec::doubleToString( ios::scientific, 6, (*force)[X] )
		<< " " <<
		GrainsExec::doubleToString( ios::scientific, 6, (*force)[Y] )
		<< " " <<
		GrainsExec::doubleToString( ios::scientific, 6, (*force)[Z] )
		<< " " <<
		GrainsExec::doubleToString( ios::scientific, 6, (*torque)[X] )
		<< " " <<
		GrainsExec::doubleToString( ios::scientific, 6, (*torque)[Y] )
		<< " " <<
		GrainsExec::doubleToString( ios::scientific, 6, (*torque)[Z] )
		<< " " << endl;
          OUT.close();
        }
      }
    }
  }

  if ( !enforceOutput )
  {
    ++m_outputTorsorObstacles_counter;
    if ( m_outputTorsorObstacles_counter ==
    	m_outputTorsorObstacles_frequency )
    m_outputTorsorObstacles_counter = 0 ;
  }
}




// ----------------------------------------------------------------------------
// Initialises output files to write loads on obstacles
void AllComponents::initialiseOutputObstaclesLoadFiles( int rank,
	bool coupledFluid, double time )
{
  m_outputTorsorObstacles_counter = coupledFluid ;

  if ( rank == 0 )
  {
    if ( GrainsExec::m_ReloadType == "new" )
    {
      string cmd = "bash " + GrainsExec::m_GRAINS_HOME
     	+ "/Tools/ExecScripts/ObstaclesLoadFiles_clear.exec "
	+ m_outputTorsorObstacles_dir;
      GrainsExec::m_return_syscmd = system( cmd.c_str() );
    }
    else
      for (list<Obstacle*>::iterator
    	obstacle=m_outputTorsorObstacles.begin();
	obstacle!=m_outputTorsorObstacles.end();obstacle++)
        GrainsExec::checkTime_outputFile( m_outputTorsorObstacles_dir
      		+ "/Loading_" + (*obstacle)->getName() + ".res",
    		time ) ;
  }

}




// ----------------------------------------------------------------------------
// Sets a random translational and angular velocity to all particles
void AllComponents::setRandomMotion( double const& coefTrans,
	double const& coefRot )
{
  for (list<Particle*>::iterator particle=m_ActiveParticles.begin();
	particle!=m_ActiveParticles.end();particle++)
    (*particle)->setRandomMotion( coefTrans, coefRot );
}




// ----------------------------------------------------------------------------
// Initialize all contact map entries to false in all particles
// and all elementary obstacles
void AllComponents::setAllContactMapToFalse()
{
  for (list<Particle*>::iterator particle=m_ActiveParticles.begin();
	particle!=m_ActiveParticles.end();particle++)
    (*particle)->setContactMapToFalse();

  list<SimpleObstacle*> obstacles = m_obstacle->getObstacles();
  list<SimpleObstacle*>::iterator myObs;
  for( myObs=obstacles.begin(); myObs!=obstacles.end(); myObs++ )
    (*myObs)->setContactMapToFalse();
}




// ----------------------------------------------------------------------------
// Update all contact map entries in all particles
// and all elementary obstacles
void AllComponents::updateAllContactMaps()
{
  for (list<Particle*>::iterator particle=m_ActiveParticles.begin();
	particle!=m_ActiveParticles.end();particle++)
    (*particle)->updateContactMap();

  list<SimpleObstacle*> obstacles = m_obstacle->getObstacles();
  list<SimpleObstacle*>::iterator myObs;
  for( myObs=obstacles.begin(); myObs!=obstacles.end(); myObs++ )
    (*myObs)->updateContactMap();
}




// ----------------------------------------------------------------------------
// Set all contact map entry features to zero in all particles
// and all elementary obstacles */
void AllComponents::setAllContactMapFeaturesToZero()
{
  for (list<Particle*>::iterator particle=m_ActiveParticles.begin();
	particle!=m_ActiveParticles.end();particle++)
    (*particle)->setContactMapFeaturesToZero();

  list<SimpleObstacle*> obstacles = m_obstacle->getObstacles();
  list<SimpleObstacle*>::iterator myObs;
  for( myObs=obstacles.begin(); myObs!=obstacles.end(); myObs++ )
    (*myObs)->setContactMapFeaturesToZero();
}




// ----------------------------------------------------------------------------
// Returns the number of particles in this process
size_t AllComponents::getNumberParticles() const
{
  return ( m_nb_particles );
}




// ----------------------------------------------------------------------------
// Returns the number of particles total number of particles in the
// system (i.e. on all subdomains/processes)
size_t AllComponents::getTotalNumberParticles() const
{
  return ( m_total_nb_particles );
}




// ----------------------------------------------------------------------------
// Returns the number of active particles
size_t AllComponents::getNumberActiveParticles() const
{
  return ( m_nb_active_particles );
}




// ----------------------------------------------------------------------------
// Returns the total number of active particles in the system (i.e. 
// on all subdomains/processes)
size_t AllComponents::getTotalNumberActiveParticles() const
{
  return ( m_total_nb_active_particles );
}




// ----------------------------------------------------------------------------
// Returns the number of active particles in this process with tag 0 ou 1
size_t AllComponents::getNumberActiveParticlesOnProc() const
{
  return ( m_nb_active_particles_on_proc );
}




// ----------------------------------------------------------------------------
// Returns the total number of active particles with tag 0 ou 1 in the system 
// (i.e. on all subdomains/processes)
size_t AllComponents::getNumberActiveParticlesOnAllProc() const
{
  return ( m_total_nb_active_particles_on_all_procs );
}




// ----------------------------------------------------------------------------
// Returns the number of inactive particles
size_t AllComponents::getNumberInactiveParticles() const
{
  return ( m_nb_inactive_particles );
}




// ----------------------------------------------------------------------------
// Returns the total number of particles in the physical system 
// (i.e. on all subdomains/processes), i.e. sum of total number of active 
// particles with tag 0 or 1 and inactive particles
size_t AllComponents::getTotalNumberPhysicalParticles() const
{
  return ( m_total_nb_physical_particles );
}




// ----------------------------------------------------------------------------
// Computes and sets the numbers of particles in the system */
void AllComponents::computeNumberParticles( GrainsMPIWrapper const* wrapper )
{
  m_nb_inactive_particles = m_InactiveParticles.size();  

  m_nb_active_particles = m_ActiveParticles.size();

  if ( wrapper ) m_total_nb_active_particles = wrapper->sum_UNSIGNED_INT( 
  	 m_nb_active_particles );
  else m_total_nb_active_particles = m_nb_active_particles;  

  m_nb_particles = m_nb_active_particles + m_nb_inactive_particles;

  m_total_nb_particles = m_total_nb_active_particles + m_nb_inactive_particles;

  m_nb_active_particles_on_proc = m_nb_active_particles;
  for (list<Particle*>::const_iterator il=m_ActiveParticles.begin();
  	il!=m_ActiveParticles.end();il++)
    if ( (*il)->getTag() == 2 ) m_nb_active_particles_on_proc--; 

  if ( wrapper ) m_total_nb_active_particles_on_all_procs = 
  	wrapper->sum_UNSIGNED_INT( m_nb_active_particles_on_proc );
  else m_total_nb_active_particles_on_all_procs = m_nb_active_particles_on_proc;
  
  m_total_nb_physical_particles = m_total_nb_active_particles_on_all_procs
  	+ m_nb_inactive_particles;
	
  GrainsExec::setTotalNumberPhysicalParticles( 
  	m_total_nb_physical_particles );	
}




// ----------------------------------------------------------------------------
// Updates list of particles in parallel
void AllComponents::updateParticleLists( double time, 
	list<Particle*>* newHaloPart )
{
  newHaloPart->clear();
  
  list<Particle*>::iterator particle;
  int tag = 0, tagnm1 = 0;
  for (particle=m_ActiveParticles.begin();particle!=m_ActiveParticles.end(); 
  	particle++)
  { 
    tag = (*particle)->getTag();
    tagnm1 = (*particle)->getTagNm1();
    
    switch ( tagnm1 )
    {
      case 0:
        // Interior to halozone (0 -> 1)
	if ( tag == 1 ) 
	{
	  m_ParticlesInHalozone.push_back( *particle);
	  newHaloPart->push_back( *particle);
          if ( GrainsExec::m_MPI_verbose )
	  {
	    ostringstream oss;
	    oss << "   t=" << GrainsExec::doubleToString( time, TIMEFORMAT ) <<
		" Interior to Halozone (0 -> 1)               Id = " <<
      		(*particle)->getID() << " " << *(*particle)->getPosition()
		<< endl;
	    GrainsMPIWrapper::addToMPIString( oss.str() );
	  } 
	}
	break;
	
      case 1:
        switch ( tag )
	{
	  // Halozone to interior (1 -> 0)
	  case 0:
	    removeParticleFromList( m_ParticlesInHalozone, *particle );
            if ( GrainsExec::m_MPI_verbose )
	    {
              ostringstream oss;
              oss << "   t=" << GrainsExec::doubleToString( time, TIMEFORMAT )
      		<< " Halozone to Interior (1 -> 0)               Id = " <<
      		(*particle)->getID() << " " << *(*particle)->getPosition()
		<< endl;
              GrainsMPIWrapper::addToMPIString( oss.str() );
	    } 
            break;
	    
	  // Halozone to halozone ( 1 -> 1)
	  case 1:
	    if ( (*particle)->getGeoPosition() 
	    	!= (*particle)->getGeoPositionNm1() )
	    {
	      newHaloPart->push_back( *particle);
              if ( GrainsExec::m_MPI_verbose )
	      {
                ostringstream oss;
                oss << "   t=" << GrainsExec::doubleToString( time, TIMEFORMAT )
      		<< " Halozone to Halozone (1 -> 1)               Id = " <<
      		(*particle)->getID() << " " << *(*particle)->getPosition()
		<< endl;
		oss << "                From " <<
		Cell::getGeoPositionName((*particle)->getGeoPositionNm1())
		<< " to " << Cell::getGeoPositionName(
			(*particle)->getGeoPosition()) << endl;
                GrainsMPIWrapper::addToMPIString( oss.str() );
	      } 
	    }	      	  
	    break;
	    
	  // Halozone to Clone (1 -> 2)
	  case 2:
            removeParticleFromList( m_ParticlesInHalozone, *particle );
	    m_CloneParticles.push_back( *particle);
	    if ( GrainsExec::m_MPI_verbose )
            {
              ostringstream oss;
              oss << "   t=" << GrainsExec::doubleToString( time, TIMEFORMAT )
      		<< " Halozone to Clone (1 -> 2)                  Id = " <<
      		(*particle)->getID() << " " << *(*particle)->getPosition()
		<< endl;
              GrainsMPIWrapper::addToMPIString( oss.str() );
            }
	    break;
	}
	break;      
    
      case 2:
        // Clone to halozone (2 -> 1)
	if ( tag == 1 ) 
	{
          removeParticleFromList( m_CloneParticles, *particle );
	  m_ParticlesInHalozone.push_back( *particle);          
          if ( GrainsExec::m_MPI_verbose )
	  {
            ostringstream oss;
            oss << "   t=" << GrainsExec::doubleToString( time, TIMEFORMAT ) <<
      		" Clone to Halozone (2 -> 1)                  Id = " <<
      		(*particle)->getID() << " " << *(*particle)->getPosition()
		<< endl;
            GrainsMPIWrapper::addToMPIString( oss.str() );
	  } 
	}	
	break;
    }
  }   
}
