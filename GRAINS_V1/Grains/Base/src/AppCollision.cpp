#include "GrainsMPIWrapper.hh"
#include "GrainsExec.hh"
#include "AppCollision.hh"
#include "Component.hh"
#include "ContactBuilderFactory.hh"
#include "PointContact.hh"
#include "SimpleObstacle.hh"


size_t AppCollision::m_allforces_blocksize = 128;


// ----------------------------------------------------------------------------
// Default constructor
AppCollision::AppCollision() 
  : App()
  , m_obstacles( NULL )
  , m_overlap_max( 0. )
  , m_overlap_mean( 0. )
  , m_time_overlapMax( 0. )
  , m_nbIterGJK_mean( 0. )
  , m_nbParticles_mean( 0. )
  , m_allforces_index( 0 )
{
  struct PointForcePostProcessing pfpp;
  pfpp.comp0 = NULL;
  pfpp.comp1 = NULL;  
  m_allforces.reserve( m_allforces_blocksize );
  for (size_t i=0;i<m_allforces_blocksize;++i) m_allforces.push_back( pfpp );
}




// ----------------------------------------------------------------------------
// Destructeur
AppCollision::~AppCollision()
{
  m_allforces.clear(); 
}




// ----------------------------------------------------------------------------
// Returns whether a particle is in contact with another component
// using the method Component::isContact
bool AppCollision::isContact( Particle const* particle ) const
{
  return ( false );
}




// ----------------------------------------------------------------------------
// Returns whether a particle is in contact with another component
// using the method Component::isContactWithCrust
bool AppCollision::isContactWithCrust( Particle const* particle ) const
{
  return ( false );
}




// ----------------------------------------------------------------------------
// Returns whether a particle is in contact with another component
// using the method Component::isClose
bool AppCollision::isClose( Particle const* particle ) const
{
  return ( false );
}




// ----------------------------------------------------------------------------
// Returns whether a particle is in contact with another component
// using the method Component::isCloseWithCrust
bool AppCollision::isCloseWithCrust( Particle const* particle ) const
{
  return ( false );
}




// ----------------------------------------------------------------------------
// Links the parent obstacle with the contact detection algorithm
void AppCollision::Link( Obstacle* obstacle )
{
  if ( !m_obstacles )
  {
    m_obstacles = obstacle;
    m_allSimpleObstacles = obstacle->getObstacles();
  }
}




// ----------------------------------------------------------------------------
// Resets the list of simple obstacles */
void AppCollision::resetListSimpleObstacles()
{
  m_allSimpleObstacles = m_obstacles->getObstacles();
}
 



// ----------------------------------------------------------------------------
// Removes an obstacle from the contact detection algorithm
void AppCollision::remove( SimpleObstacle* obs )
{
  list<SimpleObstacle*>::iterator il;
  for (il=m_allSimpleObstacles.begin();il!=m_allSimpleObstacles.end();)
    if ( *il == obs ) il = m_allSimpleObstacles.erase( il );
    else il++;   
}




// ----------------------------------------------------------------------------
// Computes the average number of particles in the simulation, the
// average is performed over the number of times this method is called
void AppCollision::computeMeanNbParticles( size_t const& nbPart )
{
  static double global_counter_nbParticles_mean = 0.;
  m_nbParticles_mean = ( m_nbParticles_mean * global_counter_nbParticles_mean
  	+ double(nbPart) ) / ( global_counter_nbParticles_mean + 1. );
  global_counter_nbParticles_mean += 1.;			
}




// ----------------------------------------------------------------------------
// Adds a contact point to the contact statistics
void AppCollision::addToContactsFeatures( double time,
	PointContact const& contactPoint )
{
  static double global_counter_overlap_mean = 0.; 
  double overlap = contactPoint.getOverlapDistance();
  if ( overlap < 0. )
  {
    overlap = fabs( overlap ) ;
    
    // Overlap moyen
    m_overlap_mean = ( m_overlap_mean * global_counter_overlap_mean + overlap )
    	/ ( global_counter_overlap_mean + 1. );
    global_counter_overlap_mean += 1.;
    
    // Overlap max
    if ( overlap > m_overlap_max )
    {
      m_overlap_max = overlap;
      m_time_overlapMax = time;
    }	 
  }
  
  static double global_counter_nbIterGJK_mean = 0.;
  int nbGJK = contactPoint.getNbIterGJK();
  if ( nbGJK )
  {
    m_nbIterGJK_mean = ( m_nbIterGJK_mean * global_counter_nbIterGJK_mean 
    	+ double(nbGJK) ) / ( global_counter_nbIterGJK_mean + 1. );
    global_counter_nbIterGJK_mean += 1.; 
  }     
}




// ----------------------------------------------------------------------------
// Returns maximum overlap
double AppCollision::getOverlapMean()
{
  return ( m_overlap_mean );
} 




// ----------------------------------------------------------------------------
// Returns average overlap
double AppCollision::getOverlapMax()
{
  return ( m_overlap_max );
} 




// ----------------------------------------------------------------------------
// Returns time of maximum overlap
double AppCollision::getTimeOverlapMax()
{
  return ( m_time_overlapMax );
}




// ----------------------------------------------------------------------------
// Returns average number of iterations of GJK for convergence
double AppCollision::getNbIterGJKMean()
{
  return ( m_nbIterGJK_mean );
}




// ----------------------------------------------------------------------------
// Returns average number of particles
double AppCollision::getNbParticlesPerProcMean()
{
  return ( m_nbParticles_mean );
}




// ----------------------------------------------------------------------------
// Sets the contact statistics
void AppCollision::setContactsFeatures( double const& overlap_max_,
	  double const& overlap_mean_,
	  double const& time_overlapMax_,
	  double const& nbIterGJK_ )
{
  m_overlap_max = overlap_max_;
  m_overlap_mean = overlap_mean_;
  m_time_overlapMax = time_overlapMax_;
  m_nbIterGJK_mean = nbIterGJK_;
}




// ----------------------------------------------------------------------------
// Returns a pointer to a simple obstacle based on its ID number
SimpleObstacle* AppCollision::getSimpleObstacle( int const& num ) const
{
  list<SimpleObstacle*>::const_iterator il;
  SimpleObstacle* pp = NULL;
  bool found = false;
  for (il=m_allSimpleObstacles.begin();il!=m_allSimpleObstacles.end() && !found;il++)
    if ( (*il)->getID() == num )
    {
      pp = *il;
      found = true;
    }
  
  return ( pp );
}




// ----------------------------------------------------------------------------
// Resets postprocessing force index to 0
void AppCollision::resetPPForceIndex()
{
  m_allforces_index = 0;
}




// ----------------------------------------------------------------------------
// Adds a postprocessing force
void AppCollision::addPPForce( Point3 const& pc, Vector3 const& force,
	Component* comp0_, Component* comp1_)
{
  // Re-allocate 
  if ( m_allforces_index >= m_allforces.size() )
  {
    struct PointForcePostProcessing pfpp;
    pfpp.comp0 = NULL;
    pfpp.comp1 = NULL;
    size_t prevsize = m_allforces.size();  
    m_allforces.reserve( prevsize + m_allforces_blocksize );
    for (size_t i=0;i<m_allforces_blocksize;++i) m_allforces.push_back( pfpp );
  }
  
  // Sets contact force features
  m_allforces[m_allforces_index].geometricPointOfContact = pc;
  m_allforces[m_allforces_index].contactForce = force;
  m_allforces[m_allforces_index].comp0 = comp0_;
  m_allforces[m_allforces_index].comp1 = comp1_;
  ++m_allforces_index;
}



    
// ----------------------------------------------------------------------------
// Returns the number of postprocessing forces
size_t AppCollision::getNbPPForces() const
{
  return ( m_allforces_index );
}




// ----------------------------------------------------------------------------
// Returns a pointer to the vector of postprocessing forces
vector<struct PointForcePostProcessing> const* AppCollision::getPPForces() const
{
  return ( &m_allforces );
}
