#include "GrainsMPIWrapper.hh"
#include "AppCollision.hh"
#include "GrainsExec.hh"
#include "Component.hh"
#include "ContactBuilderFactory.hh"
#include "PointContact.hh"
#include "SimpleObstacle.hh"
#include "Matrix.hh"
#include "GrainsBuilderFactory.hh"


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
  , m_outputForceStats( false )
  , m_outputForceStats_dir( "none" )
  , m_outputForceStats_counter( 0 )
  , m_outputForceStats_frequency( 0 )
{
  struct PointForcePostProcessing pfpp; 
  m_allforces.reserve( m_allforces_blocksize );
  for (size_t i=0;i<m_allforces_blocksize;++i) m_allforces.push_back( pfpp );
  m_stressTensor.setValue( 0., 0., 0., 0., 0., 0., 0., 0., 0. );
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
// Links the root obstacle with the contact detection algorithm at
// the start of the simulation
void AppCollision::Link( Obstacle* root_obstacle )
{
  if ( !m_obstacles )
  {
    m_obstacles = root_obstacle;
    m_allSimpleObstacles = m_obstacles->getObstacles();
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
  for (il=m_allSimpleObstacles.begin();il!=m_allSimpleObstacles.end() && !found;
  	il++)
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
  // Test if we store this force to avoif duplicate force due to multi-procs 
  // or periodic BC
  if ( comp0_->storePPForce( comp1_ ) )
  {  
    // Re-allocate 
    if ( m_allforces_index >= m_allforces.size() )
    {
      struct PointForcePostProcessing pfpp;
      size_t prevsize = m_allforces.size();  
      m_allforces.reserve( prevsize + m_allforces_blocksize );
      for (size_t i=0;i<m_allforces_blocksize;++i) 
        m_allforces.push_back( pfpp );
    }
  
    // Sets contact force features
    m_allforces[m_allforces_index].geometricPointOfContact = pc;
    m_allforces[m_allforces_index].contactForceComp0 = force;
    m_allforces[m_allforces_index].PPptComp0 = comp0_->isObstacle() ? 
  	pc : *(comp0_->getPosition());
    m_allforces[m_allforces_index].PPptComp1 = comp1_->isObstacle() ? 
  	pc : *(comp1_->getPosition());
    ++m_allforces_index;
  }
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




// ----------------------------------------------------------------------------
// Computes the stress tensor in the whole domain */
void AppCollision::computeStressTensor( GrainsMPIWrapper const* wrapper )
{
  size_t m;
  int i, j;
  Vector3 rc0, rc1;
  
  // Reset the tensor to 0
  m_stressTensor.setValue( 0., 0., 0., 0., 0., 0., 0., 0., 0. );
  
  // Loop over all contact forces
  for (m=0;m<m_allforces_index;++m)
  {
    rc0 = m_allforces[m].geometricPointOfContact - m_allforces[m].PPptComp0;
    rc1 = m_allforces[m].geometricPointOfContact - m_allforces[m].PPptComp1;
    for (i=0;i<3;++i)
      for (j=0;j<3;++j)
      {
	m_stressTensor[i][j] += rc0[i] * m_allforces[m].contactForceComp0[j];
	m_stressTensor[i][j] -= rc1[i] * m_allforces[m].contactForceComp0[j];	
      }
  }

  // If in MPI, sum contributions from each subdomain
  if ( wrapper ) m_stressTensor = wrapper->sum_Matrix( m_stressTensor );

  // Divide by the domain volume
  double volume = m_domain_global_size[X] * m_domain_global_size[Y]
  	* ( GrainsBuilderFactory::getContext() == DIM_2 ? 1. : 
		m_domain_global_size[Z] );
  m_stressTensor /= volume;  
}




// ----------------------------------------------------------------------------
// Sets the parameters to output force statistics
void AppCollision::setForceStatsParameters( string const& root_,
  	size_t const& freq_ )
{
  m_outputForceStats = true;
  m_outputForceStats_dir = root_;
  m_outputForceStats_frequency = freq_;
}




// ----------------------------------------------------------------------------
// Returns whether to output force statistics at this time
bool AppCollision::outputForceStatsAtThisTime( bool enforceOutput, 
    	bool increaseCounterOnly )
{
  bool output = false;
  if ( m_outputForceStats )
  {
    if ( ( m_outputForceStats_counter == 0 || enforceOutput )
  	&& !increaseCounterOnly ) output = true;

    if ( !enforceOutput )
    {
      ++m_outputForceStats_counter;
      if ( m_outputForceStats_counter == m_outputForceStats_frequency )
      m_outputForceStats_counter = 0 ;
    }
  }
  
  return ( output );
}




// ----------------------------------------------------------------------------
// Writes load on obstacles in a file
void AppCollision::outputForceStats( double time, double dt, int rank,
    	GrainsMPIWrapper const* wrapper )
{
  // Compute macro stress tensor 
  computeStressTensor( wrapper );

  // Output to a file
  if ( rank == 0 )
  {
    ofstream OUT( ( m_outputForceStats_dir + "/ForceStats.res" ).c_str(), 
  	ios::app );
    OUT << GrainsExec::doubleToString( ios::scientific, 6, time ) 
	<< " " <<
	GrainsExec::doubleToString( ios::scientific, 6, m_stressTensor[X][X] )
	<< " " <<
	GrainsExec::doubleToString( ios::scientific, 6, m_stressTensor[X][Y] )
	<< " " <<
	GrainsExec::doubleToString( ios::scientific, 6, m_stressTensor[X][Z] )
	<< " " <<
	GrainsExec::doubleToString( ios::scientific, 6, m_stressTensor[Y][Y] )
	<< " " <<
	GrainsExec::doubleToString( ios::scientific, 6, m_stressTensor[Y][Z] )
	<< " " <<
	GrainsExec::doubleToString( ios::scientific, 6, m_stressTensor[Z][Z] )
	<< " " <<
	GrainsExec::doubleToString( ios::scientific, 6, 
		- m_stressTensor.trace() / 3. )	<< " " << endl;
    OUT.close();
  }
}




// ----------------------------------------------------------------------------
// Initialises output files to write force statistics
void AppCollision::initialiseForceStatsFiles( int rank,
	bool coupledFluid, double time )
{
  m_outputForceStats_counter = coupledFluid ;

  if ( rank == 0 )
  {
    if ( GrainsExec::m_ReloadType == "new" )
    {
      string cmd = "bash " + GrainsExec::m_GRAINS_HOME
     	+ "/Tools/ExecScripts/ForceStatsFiles_clear.exec "
	+ m_outputForceStats_dir;
      GrainsExec::m_return_syscmd = system( cmd.c_str() );
    }
    else
       GrainsExec::checkTime_outputFile( m_outputForceStats_dir
      		+ "/ForceStats.res", time ) ;
  }

}




// ----------------------------------------------------------------------------
// Returns the macroscopic stress tensor in the whole domain */
Matrix const* AppCollision::getStressTensor() const
{
  return ( &m_stressTensor );
}




// ----------------------------------------------------------------------------
// Returns a component of the macroscopic stress tensor in the whole domain 
double AppCollision::getStressTensorComponent( int k, int l ) const
{
  return ( m_stressTensor[k][l] );
}
