#include "GrainsMPIWrapper.hh"
#include "GrainsExec.hh"
#include "AllComponents.hh"
#include "AllInsertionWindows.hh"
#include "App.hh"
#include "AppCollision.hh"
#include "RigidBodyWithCrust.hh"
#include "CompositeObstacle.hh"
#include "ObstacleBuilderFactory.hh"
#include "PostProcessingWriter.hh"
#include "ParaviewPostProcessingWriter.hh"
#include "GrainsBuilderFactory.hh"
#include "TrilobeCylinder.hh"
#include "QuadrilobeCylinder.hh"
#include <math.h>
#include <stdlib.h>



// ----------------------------------------------------------------------------
// Default constructor
AllComponents::AllComponents()
  : m_nb_physical_particles_to_insert( 0 )
  , m_wait( NULL )
  , m_nb_particles( 0 )
  , m_total_nb_particles( 0 )
  , m_nb_active_particles( 0 )
  , m_total_nb_active_particles( 0 )
  , m_nb_active_particles_on_proc( 0 )
  , m_total_nb_active_particles_on_all_procs( 0 )
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

  for (particle=m_RemovedParticles.begin(); 
  	particle!=m_RemovedParticles.end();particle++)  delete *particle;
	
  if ( m_wait ) delete m_wait;

  m_ActiveParticles.clear();
  m_RemovedParticles.clear();
  m_PeriodicCloneParticles.clear();
  m_ParticlesInBufferzone.clear();
  m_CloneParticles.clear();

  vector<Particle*>::iterator ivp;
  for (ivp=m_ReferenceParticles.begin();ivp!=m_ReferenceParticles.end(); ivp++)
    delete *ivp;
  m_ReferenceParticles.clear();

  delete m_obstacle;

  list<PostProcessingWriter*>::iterator pp;
  for (pp=m_postProcessors.begin();pp!=m_postProcessors.end();pp++)
    delete *pp;
  m_postProcessors.clear();
 
  // Note: the loadings are destroyed by the destructors of the classes
  // ObstacleKinematicsVelocity and ObstacleKinematicsForce
  // Hence we free the lists here but do not destroy the pointed objects
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
        m_RemovedParticles.push_back(*particle);
        particle = m_ActiveParticles.erase( particle );
        break;

      default:
        break;
    }
  }
}




// ----------------------------------------------------------------------------
// Adds a reference particle
void AllComponents::AddReferenceParticle( Particle* particle, size_t const &n )
{
  bool exist = false;
  size_t ntypes = m_ReferenceParticles.size(), i, i0;

  // Check whether such a type of particle already exists
  for (i=0;i<ntypes && !exist;++i)
    if ( particle->equalType( m_ReferenceParticles[i] ) ) 
    {
      exist = true;
      i0 = i;
    }  
  
  // If it does not exist, we add it to the vector of reference particles
  // If it exists, we add the number of new particles to this type
  if ( !exist )
  {
    m_ReferenceParticles.reserve( m_ReferenceParticles.size() + 1 );
    m_ReferenceParticles.push_back( particle );
    m_NbRemainingParticlesToInsert.reserve( 
  	m_NbRemainingParticlesToInsert.size() + 1 );
    m_NbRemainingParticlesToInsert.push_back( n );
  }
  else
  {
    m_NbRemainingParticlesToInsert[i0] += n;
    delete particle;
    particle = m_ReferenceParticles[i0];
  }
  m_nb_physical_particles_to_insert += n;  
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
// Associates the imposed force to the obstacle
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
	double const& dt_obstacle,
	LinkedCell const* LC )
{
  double subinterval = 0.;
  bool anyactive = false;
  static bool anyactive_previousdt = false;
  
  try 
  {
    // Particles motion
    list<Particle*>::iterator particle;
    for (particle=m_ActiveParticles.begin();
      particle!=m_ActiveParticles.end(); particle++)
      if ( (*particle)->getTag() != 2 )
        (*particle)->Move( time, dt_particle_vel, dt_particle_disp );

    // Obstacles motion
    list<SimpleObstacle*> displacedObstacles;
    if ( !m_AllImposedVelocitiesOnObstacles.empty()
  	|| !m_AllImposedForcesOnObstacles.empty() )
    {            
      anyactive = false;
    
      // Update stress dependent imposed velocity and check if any imposed 
      // velocity/force is active over this time interval
      list<ObstacleImposedVelocity*>::iterator il;
      for (il=m_AllImposedVelocitiesOnObstacles.begin();
  	il!=m_AllImposedVelocitiesOnObstacles.end();il++)
        if ( (*il)->isActif( time - dt_obstacle, time, dt_obstacle, 
      		subinterval ) )
        {
	  (*il)->updateImposedVelocity( LC );
	  anyactive = true;
        }
      list<ObstacleImposedForce*>::iterator il_F;
      for (il_F=m_AllImposedForcesOnObstacles.begin();
  	il_F!=m_AllImposedForcesOnObstacles.end() && !anyactive;il_F++) 
        if ( (*il_F)->isActif( time - dt_obstacle, time, dt_obstacle, 
      		subinterval ) )	
	  anyactive = true;        
    	           
      // Move obstacles
      if ( anyactive || anyactive_previousdt ) m_obstacle->resetKinematics();   
      anyactive_previousdt = anyactive;
      if ( anyactive ) displacedObstacles = 
    	m_obstacle->Move( time, dt_obstacle, false, false );

      // Update imposed velocity kinematics
      for (il=m_AllImposedVelocitiesOnObstacles.begin();
  	il!=m_AllImposedVelocitiesOnObstacles.end(); )
        if ( (*il)->isCompleted( time, dt_obstacle ) )
          il = m_AllImposedVelocitiesOnObstacles.erase( il );
        else il++;

      // Update imposed force kinematics    
      for (il_F=m_AllImposedForcesOnObstacles.begin();
  	il_F!=m_AllImposedForcesOnObstacles.end(); )
      {
        if ( (*il_F)->isCompleted( time, dt_obstacle ) )
          il_F = m_AllImposedForcesOnObstacles.erase( il_F );
        else il_F++;
      }
    }

    return ( displacedObstacles );
  }
  catch ( MotionError const& ) { throw; }
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
void AllComponents::Initialize_RBWCTransformState( double time, double dt )
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

  // Reference types
  for (ivp=m_ReferenceParticles.begin();
  	ivp!=m_ReferenceParticles.end(); ivp++)
    (*ivp)->computeWeight();

  // Active particles
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
// Returns a pointer to the particle to be inserted
Particle* AllComponents::getParticleToInsert( PullMode mode )
{
  if ( m_nb_physical_particles_to_insert )
  {
    size_t type = 0;
    if ( !m_wait )
    {
      switch ( mode )
      {
        case PM_ORDERED:                
          while ( !m_NbRemainingParticlesToInsert[type] ) type++;
          break;

        case PM_RANDOM:
          if ( m_wait == NULL )
	  {
	    // Find the set of particle classes that still have particles to
	    // insert
	    size_t nrem = 0, j = 0, i;
	    for (i=0;i<m_NbRemainingParticlesToInsert.size();++i)
	      if ( m_NbRemainingParticlesToInsert[i] ) nrem++;
	    vector<size_t> setrem( nrem, 0 );
	    for (i=0;i<m_NbRemainingParticlesToInsert.size();++i)
	      if ( m_NbRemainingParticlesToInsert[i] ) 
	      { 
	        setrem[j] = i;
		++j;
              }
	      
	    // Randomly pick a class
	    double v = double(random()) / double(INT_MAX);   	    
	    j = size_t( double(setrem.size()) * v );
	    type = setrem[j];
          }
          break;
      }
      m_wait = m_ReferenceParticles[type]->createCloneCopy( true );           
    }
  }
  else m_wait = NULL;

  return ( m_wait );
}




// ----------------------------------------------------------------------------
// Returns a pointer to the particle to be inserted
Particle* AllComponents::getParticleToInsert( int const& geomtype )
{
  if ( !m_wait )
    m_wait = m_ReferenceParticles[geomtype]->createCloneCopy( true );

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
// Returns a pointer to the list of removed particles
list<Particle*>* AllComponents::getRemovedParticles()
{
  return ( &m_RemovedParticles );
}




// ----------------------------------------------------------------------------
// Returns a pointer to the list of particles in the buffer zone
list<Particle*>* AllComponents::getParticlesInBufferzone()
{
  return ( &m_ParticlesInBufferzone );
}




// ----------------------------------------------------------------------------
// Returns a pointer to the list of particles in the buffer zone
list<Particle*> const* AllComponents::getParticlesInBufferzone() const
{
  return ( &m_ParticlesInBufferzone );
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
// Returns the cumulative volume of all particles, both active and to insert
double AllComponents::getVolume() const
{
  double volume = 0.;
  list<Particle*>::const_iterator particle;
  size_t i, n = m_NbRemainingParticlesToInsert.size();

  for (particle=m_ActiveParticles.cbegin();
  	particle!=m_ActiveParticles.cend(); particle++)
    volume += (*particle)->getVolume();

  for (i=0;i<n;++i)
    volume += double(m_NbRemainingParticlesToInsert[i])
    	* m_ReferenceParticles[i]->getVolume();

  return ( volume );
}




// ----------------------------------------------------------------------------
// Returns the cumulative volume of all actives particles
double AllComponents::getVolumeIn() const
{
  double volume = 0.;
  list<Particle*>::const_iterator particle;

  for (particle=m_ActiveParticles.cbegin();
  	particle!=m_ActiveParticles.cend(); particle++)
    if ( (*particle)->getTag() != 2 )	
      volume += (*particle)->getVolume();

  return ( volume );
}




// ----------------------------------------------------------------------------
// Returns the cumulative volume of all particles to insert
double AllComponents::getVolumeOut() const
{
  double volume = 0.;
  size_t i, n = m_NbRemainingParticlesToInsert.size();

  for (i=0;i<n;++i)
    volume += double(m_NbRemainingParticlesToInsert[i])
    	* m_ReferenceParticles[i]->getVolume();

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
// Sets all active particles velocity to 0 if reset == "Reset"
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
// Adds the waiting particle to the list of active particles
void AllComponents::WaitToActive( bool const& parallel )
{
  m_wait->setActivity( COMPUTE );
  m_ActiveParticles.push_back( m_wait );
  m_nb_physical_particles_to_insert--;
  m_NbRemainingParticlesToInsert[m_wait->getGeometricType()]--;
  if ( parallel )
  {
    if ( m_wait->getTag() == 1 ) m_ParticlesInBufferzone.push_back( m_wait );
    else if ( m_wait->getTag() == 2 ) m_CloneParticles.push_back( m_wait );
  }
  m_wait->setPositionDir_nm1( Z ); 
  m_wait = NULL;
}




// ----------------------------------------------------------------------------
// Destroys the waiting particle
void AllComponents::DeleteAndDestroyWait()
{
  m_NbRemainingParticlesToInsert[m_wait->getGeometricType()]--;
  m_nb_physical_particles_to_insert--;    
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
// Note: back-compatibility for old pre-2024 restart files that features 
// particles to be inserted (previously called inactive particles) is not 
// handled.
void AllComponents::read_pre2024( istream& fileSave, string const& filename,
	GrainsMPIWrapper const* wrapper )
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
      if ( buffer == "TrilobeCylinder" )
        particle = new TrilobeCylinder( false );  
      else if ( buffer == "QuadrilobeCylinder" )
        particle = new QuadrilobeCylinder( false );        	
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
  if ( readingMode == "Hybrid" ) GrainsExec::m_readingModeHybrid = true ;

  // Reload of particles using the reference particles
  ifstream FILEbin;
  if ( GrainsExec::m_readingModeHybrid )
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
      		->getSpecificCompositeShapeName() == "TrilobeCylinder" )
        particle = new TrilobeCylinder( false );
      else if ( m_ReferenceParticles[ParticleGeomType]
      		->getSpecificCompositeShapeName() == "QuadrilobeCylinder" )
        particle = new QuadrilobeCylinder( false );	
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
              m_ParticlesInBufferzone.push_back( particle );
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
        // Note: back-compatibility for old pre-2024 restart files that
	// features particles to be inserted (previously called inactive
	// particles) is not handled.
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
  
  // Set the maximum particle ID number and minimum obstacle ID number
  setParticleMaxIDObstacleMinID( wrapper );
}




// ---------------------------------------------------------------------------
// Reloads reference particles and obstacles from an input stream
size_t AllComponents::read( istream& fileSave, vector<Point3>* known_positions, 
	int const& rank, int const& nprocs )
{
  string buffer, readingMode ;
  size_t i, nbreParticles_, nbParticleTypes_, nfiles;
  Particle *particle;

  // Number of particle types
  fileSave >> buffer >> nbParticleTypes_;

  // Read the reference particles
  m_ReferenceParticles.reserve( nbParticleTypes_ );
  for (i=0; i<nbParticleTypes_; i++)
  {
    // Read the buffer "<Particle>" or "<CompositeParticle>"
    fileSave >> buffer;

    // Construct an empty particle
    if ( buffer == "<Particle>" )
      particle = new Particle( false );
    else
    {
      fileSave >> buffer >> buffer;
      if ( buffer == "TrilobeCylinder" )
        particle = new TrilobeCylinder( false ); 
      else if ( buffer == "QuadrilobeCylinder" )
        particle = new QuadrilobeCylinder( false ); 
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
  
  // Read the number of remaining particles to insert per class
  m_NbRemainingParticlesToInsert.reserve( nbParticleTypes_ );
  for (i=0; i<nbParticleTypes_; i++) 
    m_NbRemainingParticlesToInsert.push_back( 0 );
  fileSave >> buffer;
  m_nb_physical_particles_to_insert = 0;
  for (i=0; i<nbParticleTypes_; i++) 
  {      
    fileSave >> m_NbRemainingParticlesToInsert[i];
    m_nb_physical_particles_to_insert += m_NbRemainingParticlesToInsert[i];
  }
    
  // Read the number of remaining known insertion positions
  // Only the master process stores the positions
  size_t npos = 0;
  fileSave >> buffer >> npos;
  if ( npos )
  {
    Point3 P;
    if ( rank == 0 ) known_positions->reserve( npos );
    for (i=0; i<npos; i++) 
    {
      fileSave >> P[X] >> P[Y] >> P[Z];
      if ( rank == 0 ) known_positions->push_back( P );
    }
  }      
  
  // Reload the obstacle tree
  // The highest level composite obstacle is a default composite obstacle that 
  // is created as the root of the obstacle tree, its torsor is physically 
  // meaningless and is not written, therefore not read here
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


  // Read the number of files
  fileSave >> buffer >> nfiles;
  assert( nfiles == 1 || nfiles == size_t(nprocs) );  
  if ( nprocs > 1 && nfiles == 1 )
    GrainsExec::m_ReadMPIInASingleFile = true;

  // Reload of particles using the reference particles
  if ( nprocs > 1 && !GrainsExec::m_ReadMPIInASingleFile )
  {
    vector<size_t> work( nprocs, 0 );
    for (i=0;i<size_t(nprocs);++i) 
      fileSave >> readingMode >> work[i];
    nbreParticles_ = work[rank];
  }
  else
    fileSave >> readingMode >> nbreParticles_;
    
  if ( readingMode == "Hybrid" ) GrainsExec::m_readingModeHybrid = true ;
  
  return ( nbreParticles_ );
}




// ---------------------------------------------------------------------------
// Reloads particles from an input stream 
void AllComponents::read_particles( string const& filename, size_t const& npart,
    	LinkedCell const* LC, int const& rank,
  	int const& nprocs, GrainsMPIWrapper const* wrapper )
{
  ifstream FILEin;
  int ParticleTag, ParticleGeomType, ParticleID;
  unsigned int ParticleActivity;
  Particle *particle;
  bool construct = false;
  string tline, buffer;
  multimap<int,int> PartToCell;
  pair < multimap<int,int>::iterator, 
  	multimap<int,int>::iterator > crange;
  multimap<int,int>::iterator imm;	  
  Point3 gc;
  size_t nb;
    
  // Particle file
  string pfile = filename;
  size_t pos = filename.find_first_of(".");
  pfile.erase( pos );
  pfile += "_Particles";  
  pfile = GrainsExec::fullResultFileName( pfile, nprocs > 1 
  	&& !GrainsExec::m_ReadMPIInASingleFile );  
  if ( GrainsExec::m_readingModeHybrid ) 
  {
    pfile += ".bin";
    FILEin.open( pfile.c_str(), ios::in | ios::binary );
  } 
  else FILEin.open( pfile.c_str(), ios::in ); 


  // Reload of particles using the reference particles
  for (size_t i=0; i<npart; i++)
  {
    istringstream iss;
    construct = false;  

    // Read data corresponding to this particle
    // Binary file
    if ( GrainsExec::m_readingModeHybrid )
    {
      FILEin.read( reinterpret_cast<char*>( &nb ), sizeof(size_t) );
      char* particleData = new char[nb];
      FILEin.read( particleData, nb );
      string spData( particleData, nb );
      iss.str( spData );
      delete [] particleData;
    }
    // Text file
    else
    {   
      getline( FILEin, tline );
      iss.str( tline ); 
    }
    
    // In parallel, with one restart file per sub-domain, all particles are 
    // constructed by default
    if ( !GrainsExec::m_ReadMPIInASingleFile ) construct = true;
    // With 1 MPI restart file for all sub-domains, we need to check which
    // particles geometrically belong to this sub-domain 
    else
    {        
      // Extract its ID number and position
      // Binary file
      if ( GrainsExec::m_readingModeHybrid )
      {
	iss.read( reinterpret_cast<char*>( &ParticleGeomType ), sizeof(int) );
	iss.read( reinterpret_cast<char*>( &ParticleID ), sizeof(int) );
	iss.read( reinterpret_cast<char*>( &ParticleTag ), sizeof(int) );	
	iss.read( reinterpret_cast<char*>( &ParticleActivity ), 
		sizeof(unsigned int) );
	iss.read( reinterpret_cast<char*>( &gc[X] ), sizeof(double) );
	iss.read( reinterpret_cast<char*>( &gc[Y] ), sizeof(double) );
	iss.read( reinterpret_cast<char*>( &gc[Z] ), sizeof(double) );
      }
      // Text file
      else
        iss >> buffer >> ParticleID >> buffer >> ParticleActivity >> gc[X] >> 
		gc[Y] >> gc[Z];
      
      // Check whether this particle is in the LinkedCell of this sub-domain
      // If construct, reset the input stream position to the beginning       
      if ( LC->isInLinkedCell( gc ) ) 
      {
        construct = true; 
        iss.seekg( 0, iss.beg );
      }
    }
    
    if ( construct )
    {
      // Read the geometric type of the particle
      if ( GrainsExec::m_readingModeHybrid )
        iss.read( reinterpret_cast<char*>( &ParticleGeomType ), 
		sizeof(int) );
      else
        iss >> ParticleGeomType;

      // Particle construction
      if ( m_ReferenceParticles[ParticleGeomType]->isCompositeParticle() )
      {    
        if ( m_ReferenceParticles[ParticleGeomType]
      		->getSpecificCompositeShapeName() == "TrilobeCylinder" )
          particle = new TrilobeCylinder( false );
        else if ( m_ReferenceParticles[ParticleGeomType]
      		->getSpecificCompositeShapeName() == "QuadrilobeCylinder" )
          particle = new QuadrilobeCylinder( false );	  
        else   
          particle = new CompositeParticle( false );
      }
      else
        particle = new Particle( false );

      // Set the geometric type of the particle
      particle->setGeometricType( ParticleGeomType );

      // Read the particle features
      if ( GrainsExec::m_readingModeHybrid )
        particle->read2014_binary( iss, &m_ReferenceParticles );
      else
        particle->read2014( iss, &m_ReferenceParticles );

      // Add to lists
      // We need to reset ParticleTag in case the linked cell changed
      // from the previous simulation
      m_ActiveParticles.push_back( particle );
      ParticleTag = LC->getCell( *(particle->getPosition()) )->getTag(); 
      particle->setTag( ParticleTag );
 
      // MPI mode
      if ( GrainsExec::m_MPI )
	switch ( ParticleTag )
        {
          case 1:
            m_ParticlesInBufferzone.push_back( particle );
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
    }
  }
  FILEin.close();
  // Set the maximum particle ID number and minimum obstacle ID number
  setParticleMaxIDObstacleMinID( wrapper ); 
} 




// ---------------------------------------------------------------------------
// Writes components to an output stream. Only active particles are written 
// to the stream
void AllComponents::write( ostream &fileSave, string const& filename, 
	vector<Point3> const* known_positions, CloneInReload cir, 
	LinkedCell const* LC, int const& rank,
	int const& nprocs, GrainsMPIWrapper const* wrapper ) const
{
  list<Particle*>::const_iterator particle;
  size_t nb = 0;
  
  if ( rank == 0 )
  {
    // Reference particles and obstacles
    fileSave << endl << "NumberOfParticleTypes " <<
  	m_ReferenceParticles.size() << endl;

    vector<Particle*>::const_iterator iv;
    for (iv=m_ReferenceParticles.begin();
  	iv!=m_ReferenceParticles.end();iv++) (*iv)->write( fileSave );

    fileSave << endl;
    fileSave << "RemainingNbParticlesToInsertPerType";
    for (size_t i=0;i<m_NbRemainingParticlesToInsert.size();++i)
      fileSave << " " << m_NbRemainingParticlesToInsert[i];
    fileSave << endl << endl;
    
    fileSave << "RemainingKnownInsertionPositions " << known_positions->size()
    	<< endl;
    if ( known_positions->size() )
    {
      vector<Point3>::const_iterator il;
      for (il=known_positions->cbegin();il!=known_positions->cend();il++)
        fileSave << GrainsExec::doubleToString(ios::scientific,FORMAT16DIGITS,
  		(*il)[X]) << " " << 
		GrainsExec::doubleToString(ios::scientific,FORMAT16DIGITS,
  		(*il)[Y]) << " " << 
		GrainsExec::doubleToString(ios::scientific,FORMAT16DIGITS,
  		(*il)[Z]) << endl;
    }     
              
    fileSave << endl << "<Obstacle>" << endl;
    m_obstacle->write( fileSave );
    fileSave << endl << "</Obstacle>" << endl << endl;
  }
  
  // Number of particles to write
  vector<bool> write( m_nb_active_particles, true );
  size_t m_write = 0;
  int i;
  for (particle=m_ActiveParticles.begin(),i=0;
	particle!=m_ActiveParticles.end(); particle++,++i)
  {	
    switch(cir)
    {
      case CIR_NONE: if ( (*particle)->getTag() > 1 ) write[i] = false; break;
      case CIR_NOPERIODIC: 
	if ( !LC->isInDomain( (*particle)->getPosition() ) ) write[i] = false; 
	break;
      case CIR_ALL: write[i] = true; break;
    }
    if ( write[i] ) ++m_write;
  }    

  size_t* nbpart_proc = NULL;
  if ( nprocs > 1 )
    nbpart_proc =  wrapper->Gather_UNSIGNED_INT_master( m_write );
  else
  {
    nbpart_proc = new size_t[1];
    nbpart_proc[0] = m_write;
  }
  
  if ( rank == 0 )
  {      
    fileSave << "NbFiles " << nprocs << endl;    
    for (i=0;i<nprocs;++i)  
      fileSave << ( GrainsExec::m_writingModeHybrid ? "Hybrid" :
  	"Text" ) << " " << nbpart_proc[i] << endl;
  }
  
  if ( nbpart_proc ) delete [] nbpart_proc;

	
  // Particle file
  string pfile = filename;
  size_t pos = filename.find_first_of(".");
  pfile.erase( pos );
  pfile += "_Particles";
  pfile = GrainsExec::fullResultFileName( pfile, nprocs > 1 );

  if ( GrainsExec::m_writingModeHybrid )
  {
    string binary_filename = pfile + ".bin";
    ofstream FILEbin( binary_filename.c_str(), ios::out | ios::binary );

    // Note: we first write the number of bytes corresponding to each particle
    // to be able to later read exactly the amount of bytes corresponding
    // to that particle in AllComponents::read_particles
    for (particle=m_ActiveParticles.begin(),i=0;
  	particle!=m_ActiveParticles.end(); particle++,++i)      
      if ( write[i] )
      {
        nb = (*particle)->get_numberOfBytes();
        FILEbin.write( reinterpret_cast<char*>( &nb ), sizeof(size_t) );      
        (*particle)->write2014_binary( FILEbin );
      }

    FILEbin.close();
  }
  else
  {
    ofstream FILEtext( pfile.c_str(), ios::out );    

    for (particle=m_ActiveParticles.begin(),i=0;
  	particle!=m_ActiveParticles.end(); particle++,++i)
      if ( write[i] ) (*particle)->write2014( FILEtext );
      
    FILEtext.close();    
  }
}




// ---------------------------------------------------------------------------
// Writes components to a single file MPI File in parallel
void AllComponents::write_singleMPIFile( ostream &fileSave, 
	string const& filename, vector<Point3> const* known_positions, 
	CloneInReload cir, LinkedCell const* LC, 
	GrainsMPIWrapper const* wrapper, bool periodic ) const
{  
  int rank = wrapper->get_rank(), tag, nb_particles_to_write,
  	nprocs = wrapper->get_total_number_of_active_processes();  
  MPI_Comm MPI_COMM_activeProc = wrapper->get_active_procs_comm();
  MPI_File file;
  MPI_Status status; 
  list<Particle*> to_write;
  list<Particle*>::const_iterator particle;
  Point3 const* gc;
  multimap< int, Point3 >* doNotWrite = NULL; 
  multimap< int, Point3 >::iterator imm;      

  // Particles to write to the reload file: 
  // * active particles interior and buffer
  // * actives particles periodic clones if clone writing mode is ALL
  // Some periodic clones may be duplicated and we first need to determine on 
  // each process the list of periodic clones that must not be written by this 
  // process as another process is in charge
  if ( periodic && cir == CIR_ALL ) 
  {
    doNotWrite = wrapper->doNotWritePeriodicClones(
  	m_ActiveParticles, LC );

    for (particle=m_ActiveParticles.begin();
  	particle!=m_ActiveParticles.end(); particle++)
    {
      gc = (*particle)->getPosition(); 
      tag = (*particle)->getTag();
      if ( tag < 2 ) to_write.push_back( *particle );
      else if ( !LC->isInDomain( gc ) )
      {
        if ( doNotWrite->count( (*particle)->getID() ) )
        {
          bool found = false;
          pair < multimap<int,Point3>::iterator, 
  		multimap<int,Point3>::iterator > crange = 
		doNotWrite->equal_range( (*particle)->getID() );
	  for (imm=crange.first; imm!=crange.second && !found; )
	    if ( gc->DistanceTo( imm->second ) 
	  	< 2. * (*particle)->getCrustThickness() )
	      found = true;
	    else imm++;
	  if ( !found ) to_write.push_back( *particle );
	  else doNotWrite->erase( imm );
        }
        else to_write.push_back( *particle );	
      }      
    }
    
    delete doNotWrite;      
  }		
  else
    for (particle=m_ActiveParticles.begin();
  	particle!=m_ActiveParticles.end(); particle++)
      if ( (*particle)->getTag() < 2 ) to_write.push_back( *particle );  
  
  nb_particles_to_write = wrapper->sum_INT( int(to_write.size()) );

  if ( rank == 0 )
  {
    // Reference particles and obstacles
    fileSave << endl << "NumberOfParticleTypes\t" <<
  	m_ReferenceParticles.size() << endl;

    vector<Particle*>::const_iterator iv;
    for (iv=m_ReferenceParticles.begin();
  	iv!=m_ReferenceParticles.end();iv++) (*iv)->write( fileSave );

    fileSave << endl;
    fileSave << "RemainingNbParticlesToInsertPerType";
    for (size_t i=0;i<m_NbRemainingParticlesToInsert.size();++i)
      fileSave << " " << m_NbRemainingParticlesToInsert[i];
    fileSave << endl << endl;
    
    fileSave << "RemainingKnownInsertionPositions " << known_positions->size()
    	<< endl;
    if ( known_positions->size() )
    {
      vector<Point3>::const_iterator il;
      for (il=known_positions->cbegin();il!=known_positions->cend();il++)
        fileSave << GrainsExec::doubleToString(ios::scientific,FORMAT16DIGITS,
  		(*il)[X]) << " " << 
		GrainsExec::doubleToString(ios::scientific,FORMAT16DIGITS,
  		(*il)[Y]) << " " << 
		GrainsExec::doubleToString(ios::scientific,FORMAT16DIGITS,
  		(*il)[Z]) << endl;
    }

    fileSave << endl << "<Obstacle>" << endl;
    m_obstacle->write( fileSave );
    fileSave << endl << "</Obstacle>" << endl << endl;
    fileSave << "NbFiles 1" << endl; 
    fileSave << ( GrainsExec::m_writingModeHybrid ? "Hybrid" :
  	"Text" ) << " " << nb_particles_to_write << endl;     
  }


  // Particle MPI file  
  string pfile = filename;
  size_t pos = filename.find_first_of(".");
  pfile.erase( pos );
  pfile += "_Particles";
  pfile = GrainsExec::fullResultFileName( pfile, false );
  if ( GrainsExec::m_writingModeHybrid ) pfile += ".bin"; 
   
  MPI_File_open( MPI_COMM_activeProc, pfile.c_str(), 
    	MPI_MODE_CREATE | MPI_MODE_WRONLY,
	MPI_INFO_NULL, &file ); 
 
  // Particle stream
  ostringstream oss_particles;
  vector<int> mpifile_offsets( nprocs, 0 ); 
  size_t nb = 0;
  if ( GrainsExec::m_writingModeHybrid ) 
    for (particle=to_write.begin();particle!=to_write.end(); particle++)
    { 
      // Note: we first write the number of bytes corresponding to the particle
      // to be able to later read exactly the amount of bytes corresponding
      // to that particle in AllComponents::read_particles
      nb = (*particle)->get_numberOfBytes();
      oss_particles.write( reinterpret_cast<char*>( &nb ), sizeof(size_t) );
      (*particle)->write2014_binary( oss_particles );
    } 
  else
    for (particle=to_write.begin();particle!=to_write.end(); particle++) 
      (*particle)->write2014( oss_particles );       

  int out_length = int(oss_particles.str().size());
  int* out_length_per_proc = wrapper->AllGather_INT( out_length );
  
  mpifile_offsets[0] = 0;
  for (int i=1;i<nprocs;i++)
    mpifile_offsets[i] = mpifile_offsets[i-1] 
    	+ out_length_per_proc[i-1] * sizeof(char) ;

  MPI_File_write_at_all( file, mpifile_offsets[rank], 
  	oss_particles.str().c_str(), out_length, MPI_CHAR, &status );
  
  MPI_File_close( &file );
  delete out_length_per_proc; 

}




// ----------------------------------------------------------------------------
// Debugging method
void AllComponents::debug( char* s )
{
  cout << s << '\n'
       << "Particles " << m_ActiveParticles.size()
       		+ m_nb_physical_particles_to_insert  << '\n'
       << "   Actives " << m_ActiveParticles.size() << '\t'
       << "   Wait    " << m_nb_physical_particles_to_insert << endl;
}




// ----------------------------------------------------------------------------
// Output operator
ostream& operator << ( ostream& f, AllComponents const& EC )
{
  f << "Number of particles on all processes = "
  	<< EC.m_total_nb_particles << endl;
  f << "Number of particles = " <<
  	EC.m_ActiveParticles.size() + EC.m_nb_physical_particles_to_insert 
	<< endl;
  f << "Number of active particles = " << EC.m_ActiveParticles.size()
  	<< endl;
  f << "Number of particles to insert = " << 
  	EC.m_nb_physical_particles_to_insert << endl;
  f << "Number of particles in buffer zone = " <<
  	EC.m_ParticlesInBufferzone.size() << endl;
  f << "Number of clone particles = " << EC.m_CloneParticles.size() << endl;
  list<Particle*>::const_iterator il;
  for (il=EC.m_ActiveParticles.begin();il!=EC.m_ActiveParticles.end();il++)
    f << *(*il) << endl;

  return f;
}




// ----------------------------------------------------------------------------
// Writes components for Post-Processing at the start of the simulation
void AllComponents::PostProcessing_start( double time, double dt,
	LinkedCell const* LC, AllInsertionWindows const& insert_windows,
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
	&m_RemovedParticles, postProcessingPeriodic,
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
	LinkedCell const* LC, AllInsertionWindows const& insert_windows,
	int rank, int nprocs, GrainsMPIWrapper const* wrapper )
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
	&m_RemovedParticles, postProcessingPeriodic,
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
// contact or motion
void AllComponents::PostProcessingErreurComponents( string const& filename,
	list<Component*> const& errcomposants )
{
  // Post processing writers
  list<PostProcessingWriter*>::iterator pp;
  for (pp=m_postProcessors.begin();pp!=m_postProcessors.end();pp++)
    (*pp)->writeErreurComponents_Paraview( filename, errcomposants );
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
    if ( pobs ) m_outputTorsorObstacles.push_back( pobs );
  }
}




// ----------------------------------------------------------------------------
// Computes load on obstacles
void AllComponents::computeObstaclesLoad( double time, double dt,
      	GrainsMPIWrapper const* wrapper )
{
  // Load on simple obstacles
  if ( wrapper ) wrapper->sumObstaclesLoad( m_obstacle->getObstacles() );
    
  // Load on composite obstacles 
  m_obstacle->getTorsor(); 
}



// ----------------------------------------------------------------------------
// Writes load on obstacles in a file
void AllComponents::outputObstaclesLoad( double time, double dt,
	bool enforceOutput, bool increaseCounterOnly, int rank )
{
  if ( !increaseCounterOnly )
  {
    if ( m_outputTorsorObstacles_counter == 0 || enforceOutput )
    {
      if ( rank == 0 )
      {
        Vector3 const* force = NULL;
        Vector3 const* torque = NULL;
        for (list<Obstacle*>::iterator
    		obstacle=m_outputTorsorObstacles.begin();
		obstacle!=m_outputTorsorObstacles.end();obstacle++)
        {
          ofstream OUT( ( m_outputTorsorObstacles_dir
      	+ "/Loading_" + (*obstacle)->getName() + ".res" ).c_str(), ios::app );
	  force = (*obstacle)->getForce();
          torque = (*obstacle)->getTorque();
	  OUT << GrainsExec::doubleToString( ios::scientific, 6, time ) 
	  	<< " " <<
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
// Set all contact map cumulative features to zero in all particles
// and all elementary obstacles */
void AllComponents::setAllContactMapCumulativeFeaturesToZero()
{
  for (list<Particle*>::iterator particle=m_ActiveParticles.begin();
	particle!=m_ActiveParticles.end();particle++)
    (*particle)->setContactMapCumulativeFeaturesToZero();

  list<SimpleObstacle*> obstacles = m_obstacle->getObstacles();
  list<SimpleObstacle*>::iterator myObs;
  for( myObs=obstacles.begin(); myObs!=obstacles.end(); myObs++ )
    (*myObs)->setContactMapCumulativeFeaturesToZero();
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
// Returns the total number of particles in the physical system 
// (i.e. on all subdomains/processes), i.e. sum of total number of active 
// particles with tag 0 or 1 and particles to be inserted
size_t AllComponents::getTotalNumberPhysicalParticles() const
{
  return ( m_total_nb_physical_particles );
}





// ----------------------------------------------------------------------------
// Returns the number of particles to insert in the physical system
size_t AllComponents::getNumberPhysicalParticlesToInsert() const
{
  return ( m_nb_physical_particles_to_insert );
}




// ----------------------------------------------------------------------------
// Computes and sets the numbers of particles in the system */
void AllComponents::computeNumberParticles( GrainsMPIWrapper const* wrapper )
{
  m_nb_active_particles = m_ActiveParticles.size();

  if ( wrapper ) 
    m_total_nb_active_particles = wrapper->sum_UNSIGNED_INT( 
  	 m_nb_active_particles );
  else 
    m_total_nb_active_particles = m_nb_active_particles; 

  m_nb_particles = m_nb_active_particles + m_nb_physical_particles_to_insert;

  m_total_nb_particles = m_total_nb_active_particles 
  	+ m_nb_physical_particles_to_insert;

  m_nb_active_particles_on_proc = m_nb_active_particles;
  for (list<Particle*>::const_iterator il=m_ActiveParticles.begin();
  	il!=m_ActiveParticles.end();il++)
    if ( (*il)->getTag() == 2 ) m_nb_active_particles_on_proc--; 

  if ( wrapper ) m_total_nb_active_particles_on_all_procs = 
  	wrapper->sum_UNSIGNED_INT( m_nb_active_particles_on_proc );
  else m_total_nb_active_particles_on_all_procs = m_nb_active_particles_on_proc;
  
  m_total_nb_physical_particles = m_total_nb_active_particles_on_all_procs
  	+ m_nb_physical_particles_to_insert;
	
  GrainsExec::setTotalNumberPhysicalParticles( 
  	m_total_nb_physical_particles );	
}




// ----------------------------------------------------------------------------
// Updates list of particles in parallel
void AllComponents::updateParticleLists( double time, 
	list<Particle*>* newBufPart, GrainsMPIWrapper const* wrapper )
{
  newBufPart->clear();
  
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
        // Interior to buffer (0 -> 1)
	if ( tag == 1 ) 
	{
	  m_ParticlesInBufferzone.push_back( *particle);
	  newBufPart->push_back( *particle );
          if ( GrainsExec::m_MPI_verbose )
	  {
	    ostringstream oss;
	    oss << "   t=" << GrainsExec::doubleToString( time, 
	    		FORMAT10DIGITS ) <<
		" Interior to Buffer (0 -> 1) Id = " <<
      		(*particle)->getID() << " " << *(*particle)->getPosition()
		<< endl;
	    GrainsMPIWrapper::addToMPIString( oss.str() );
	  } 
	}
	break;
	
      case 1:
        switch ( tag )
	{
	  // Buffer to interior (1 -> 0)
	  case 0:
	    removeParticleFromList( m_ParticlesInBufferzone, *particle );
            if ( GrainsExec::m_MPI_verbose )
	    {
              ostringstream oss;
              oss << "   t=" << GrainsExec::doubleToString( time, 
	      		FORMAT10DIGITS )
      		<< " Buffer to Interior (1 -> 0) Id = " <<
      		(*particle)->getID() << " " << *(*particle)->getPosition()
		<< endl;
              GrainsMPIWrapper::addToMPIString( oss.str() );
	    } 
            break;
	    
	  // Buffer to buffer ( 1 -> 1)
	  case 1:
	    if ( (*particle)->getGeoPosition() 
	    	!= (*particle)->getGeoPositionNm1() )
	    {
	      newBufPart->push_back( *particle );
              if ( GrainsExec::m_MPI_verbose )
	      {
                ostringstream oss;
                oss << "   t=" << GrainsExec::doubleToString( time, 
			FORMAT10DIGITS )
      		<< " Buffer to Buffer (1 -> 1)   Id = " <<
      		(*particle)->getID() << " " << *(*particle)->getPosition()
		<< endl;
		oss << "                From " <<
		Cell::getGeoPositionName((*particle)->getGeoPositionNm1())
		<< " to " << Cell::getGeoPositionName(
			(*particle)->getGeoPosition()) << endl;
                GrainsMPIWrapper::addToMPIString( oss.str() );
	      } 
	    }
	    else if ( GrainsExec::m_partialPer_is_active )
	    { 
	      if ( wrapper->isGeoPositionPeriodic( 
	      	(*particle)->getCell()->getGeoPosition() ) )
	        if ( GrainsExec::partialPeriodicityCompTest( 
			(*particle)->getPosition() ) &&
		!GrainsExec::partialPeriodicityCompTest(
			(*particle)->getPositionDir_nm1() ) )
	        {
	          newBufPart->push_back( *particle );
                  if ( GrainsExec::m_MPI_verbose )
	          {
                    ostringstream oss;
                    oss << "   t=" << GrainsExec::doubleToString( time, 
			FORMAT10DIGITS )
      		  	<< " Buffer to Buffer (1 -> 1)   Id = " <<
      			(*particle)->getID() << " " << 
			*(*particle)->getPosition() << endl;
		    oss << "                Periodic transition zone" << endl;
                    GrainsMPIWrapper::addToMPIString( oss.str() );
	          } 
	        }
	    }	    
	    	      	  
	    break;
	    
	  // Buffer to clone (1 -> 2)
	  case 2:
            removeParticleFromList( m_ParticlesInBufferzone, *particle );
	    m_CloneParticles.push_back( *particle );
	    if ( GrainsExec::m_MPI_verbose )
            {
              ostringstream oss;
              oss << "   t=" << GrainsExec::doubleToString( time, 
	      		FORMAT10DIGITS )
      		<< " Buffer to Clone (1 -> 2)    Id = " <<
      		(*particle)->getID() << " " << *(*particle)->getPosition()
		<< endl;
              GrainsMPIWrapper::addToMPIString( oss.str() );
            }
	    break;
	}
	break;      
    
      case 2:
        // Clone to buffer (2 -> 1)
	if ( tag == 1 ) 
	{
          removeParticleFromList( m_CloneParticles, *particle );
	  m_ParticlesInBufferzone.push_back( *particle );
	  newBufPart->push_back( *particle );          
          if ( GrainsExec::m_MPI_verbose )
	  {
            ostringstream oss;
            oss << "   t=" << GrainsExec::doubleToString( time, 
	    		FORMAT10DIGITS ) <<
      		" Clone to Buffer (2 -> 1)    Id = " <<
      		(*particle)->getID() << " " << *(*particle)->getPosition()
		<< endl;
            GrainsMPIWrapper::addToMPIString( oss.str() );
	  } 
	}	
	break;
    }
  }   
}




// ----------------------------------------------------------------------------
// Computes and sets the maximum particle ID number and minimum obstacle ID 
// number
void AllComponents::setParticleMaxIDObstacleMinID( 
	GrainsMPIWrapper const* wrapper )
{
  list<Particle*>::iterator particle;
  int maxID = 0;
  for (particle=m_ActiveParticles.begin();particle!=m_ActiveParticles.end(); 
  	particle++)
    maxID = max( maxID, (*particle)->getID() );
  
  if ( wrapper ) maxID = wrapper->max_INT( maxID );
  Particle::setMaxIDnumber( maxID );

  m_obstacle->setMinIDnumber();
}




// ----------------------------------------------------------------------------
// Sets time integration scheme in all particles using the macro variable 
// GrainsExec::m_TIScheme
void AllComponents::setTimeIntegrationScheme()
{
  list<Particle*>::iterator particle;
  vector<Particle*>::iterator ivp;

  // Reference types
  for (ivp=m_ReferenceParticles.begin();
  	ivp!=m_ReferenceParticles.end(); ivp++)
    (*ivp)->setTimeIntegrationScheme();

  // Active particles
  for (particle=m_ActiveParticles.begin();
  	particle!=m_ActiveParticles.end(); particle++)
    (*particle)->setTimeIntegrationScheme(); 
}




// ----------------------------------------------------------------------------
// Sets the center of mass coordinate in a given direction at the previous time 
// of all particles to the current value (generally called at the start of a 
// time step)
void AllComponents::setPositionDir_nm1( Direction dir )
{
  for (list<Particle*>::iterator particle=m_ActiveParticles.begin();
  	particle!=m_ActiveParticles.end();particle++) 
    (*particle)->setPositionDir_nm1( dir ); 
}




// ----------------------------------------------------------------------------
// Sets all particles' momentum equation coupling factor
void AllComponents::setParticlesCouplingFactor()
{
  for (list<Particle*>::iterator particle=m_ActiveParticles.begin();
  	particle!=m_ActiveParticles.end();particle++) 
    (*particle)->setCouplingFactor(); 
}
