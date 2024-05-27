#include "GrainsMPIWrapper.hh"
#include "App.hh"
#include "ContactBuilderFactory.hh"
#include "GrainsExec.hh"
#include "Cell.hh"
#include "RawDataPostProcessingWriter.hh"
#include <assert.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>


string *GrainsMPIWrapper::m_MPILogString = new string;
vector< vector<int> > GrainsMPIWrapper::m_particleBufferzoneToNeighboringProcs;
vector<int> GrainsMPIWrapper::m_GeoLocReciprocity;


// ----------------------------------------------------------------------------
// Constructor with domain decomposition and periodicity as input parameters
GrainsMPIWrapper::GrainsMPIWrapper( int NX, int NY, int NZ,
  	int PERX, int PERY, int PERZ, string const& oshift ) 
  : m_coords( NULL )
  , m_dim( NULL )
  , m_period( NULL )
  , m_isMPIperiodic( false )   
  , m_rank( 0 ) 
  , m_rank_master( 0 )
  , m_nprocs( 0 ) 
  , m_nprocs_world( 0 )
  , m_is_active( false )
  , m_neighbors( NULL )
  , m_commgrainsMPI_3D( NULL )
  , m_tag_INT( 1 )
  , m_tag_DOUBLE( 2 )
  , m_tag_CHAR( 3 )
{ 
  // We set all communicators such that the rank in each communicator is the 
  // same and we test that this is true. To do so, if the code runs with a 
  // subset of m_nprocs processes only, it runs with the processes numbered
  // 0 to m_nprocs-1, and the process 0 is always considered as the master
  // process

  // Total number of processes and data in the world communicator
  MPI_Comm_size( MPI_COMM_WORLD, &m_nprocs_world );
  MPI_Comm_rank( MPI_COMM_WORLD, &m_rank );  
  if ( NX * NY * NZ > m_nprocs_world )
  {
    if ( m_rank == m_rank_master )
      cout << endl << "!!! WARNING !!!" << endl << 
      	"Domain decomposition does not match total number of processus" 
	<< endl << endl;
    int error_code = 0;
    MPI_Abort( MPI_COMM_WORLD, error_code );
  }
  else
  {
    m_nprocs = NX * NY * NZ;
    if ( m_rank == m_rank_master )
    {
      cout << oshift << "Total number of processes = " << m_nprocs_world 
      	<< endl;
      cout << oshift << "Number of active processes = " << m_nprocs << endl;
      cout << oshift << "Number of sleeping processes = " << 
      	m_nprocs_world - m_nprocs << endl;
    }
  }  

  // Active process group
  MPI_Group world_group;
  MPI_Comm_group( MPI_COMM_WORLD, &world_group ); 
  int *activ_proc = new int[m_nprocs];
  for (int i=0;i<m_nprocs;++i) activ_proc[i] = i;
  MPI_Group_incl( world_group, m_nprocs, activ_proc, &m_MPI_GROUP_activProc );
  if ( m_nprocs < m_nprocs_world )
    MPI_Comm_create( MPI_COMM_WORLD, m_MPI_GROUP_activProc, 
    	&m_MPI_COMM_activeProc );
  else
    MPI_Comm_dup( MPI_COMM_WORLD, &m_MPI_COMM_activeProc );  
  if ( m_rank < m_nprocs ) m_is_active = true;   
  else m_is_active = false;  
  MPI_Group_free( &world_group );

  // Among the active processes
  if ( m_is_active )
  {     
    // Number of processes per direction
    m_dim = new int[3];
    m_dim[0] = NX, m_dim[1] = NY, m_dim[2] = NZ;
  
    // Periodicity of the domain
    m_period = new int[3];
    m_period[0] = PERX, m_period[1] = PERY, m_period[2] = PERZ;  
    if ( m_period[0] || m_period[1] || m_period[2] ) m_isMPIperiodic = true;
    
    // Creates vectors of periodicity of the domain
    m_MPIperiodes.reserve( 27 );
    for (int i=0;i<27;++i) m_MPIperiodes.push_back( Vector3Null ); 

    // Sets the relationship between the GeoPosition in a buffer
    // zone from which data are sent and the GeoPosition of the 
    // neighboring processes that receive the data
    setParticleBufferzoneToNeighboringProcs();

    // Re-number processes: no 
    int reorganisation = 0; 
      
    // Creates the MPI Cartesian topology with periodicity if any
    m_commgrainsMPI_3D = new MPI_Comm;
    MPI_Cart_create( m_MPI_COMM_activeProc, 3, m_dim, m_period, reorganisation, 
   	m_commgrainsMPI_3D );
     
    // Rank and MPI cartesian coordinates
    // We assert that the rank of active processors in all communicators is 
    // the same
    int rank_active, rank_cart; 
    MPI_Comm_rank( MPI_COMM_WORLD, &rank_active );  
    MPI_Comm_rank( *m_commgrainsMPI_3D, &rank_cart );
    assert( rank_active == m_rank );
    assert( rank_cart == m_rank );        
    m_coords = new int[3];
    MPI_Cart_coords( *m_commgrainsMPI_3D, m_rank, 3, m_coords );

    // Process neighbors in the MPI cartesian topology
    m_neighbors = new MPINeighbors( *m_commgrainsMPI_3D, m_coords, m_dim, 
    	m_period ); 
 
    // Timer
    SCT_insert_app( "BuffersCopy" );
    SCT_insert_app( "MPIComm" );
    SCT_insert_app( "UpdateCreateClones" );          
  }
    
}




// ----------------------------------------------------------------------------
// Destructor
GrainsMPIWrapper::~GrainsMPIWrapper()
{
  MPI_Group_free( &m_MPI_GROUP_activProc );

  if ( m_is_active ) 
  {
    MPI_Comm_free( &m_MPI_COMM_activeProc ); 
    delete [] m_coords;
    delete [] m_dim;
    delete [] m_period;
    MPI_Comm_free( m_commgrainsMPI_3D );
    delete m_commgrainsMPI_3D;
    delete m_neighbors;
    if ( m_MPILogString ) delete m_MPILogString;
    vector< vector<int> >::iterator ivv;
    for (ivv=m_particleBufferzoneToNeighboringProcs.begin();
  	ivv!=m_particleBufferzoneToNeighboringProcs.end();ivv++)
      ivv->clear();
    m_particleBufferzoneToNeighboringProcs.clear();    
    m_MPIperiodes.clear();  
  }
}




// ----------------------------------------------------------------------------
// Sets the relationship between the GeoPosition in a buffer
// zone from which data are sent and the GeoPosition of the neighboring
// processes that receive the data
void GrainsMPIWrapper::setParticleBufferzoneToNeighboringProcs()
{
  m_particleBufferzoneToNeighboringProcs.reserve(26);
  vector<int> emptyVECINT;
  for (int i=0;i<26;++i) 
    m_particleBufferzoneToNeighboringProcs.push_back(emptyVECINT);
  vector<int> *work;
  
  // NORTH => NORTH
  work = new vector<int>(1,0);
  (*work)[0] = GEOPOS_NORTH;
  m_particleBufferzoneToNeighboringProcs[GEOPOS_NORTH] = *work;
  work->clear();
  delete work;
  
  // NORTH_EAST => NORTH, EAST, NORTH_EAST
  work = new vector<int>(3,0);
  (*work)[0] = GEOPOS_NORTH;
  (*work)[1] = GEOPOS_EAST;  
  (*work)[2] = GEOPOS_NORTH_EAST;
  m_particleBufferzoneToNeighboringProcs[GEOPOS_NORTH_EAST] = *work;    
  work->clear();
  delete work;
  
  // NORTH_WEST => NORTH, WEST, NORTH_WEST
  work = new vector<int>(3,0);
  (*work)[0] = GEOPOS_NORTH;
  (*work)[1] = GEOPOS_WEST;  
  (*work)[2] = GEOPOS_NORTH_WEST;
  m_particleBufferzoneToNeighboringProcs[GEOPOS_NORTH_WEST] = *work;    
  work->clear();
  delete work;    

  // NORTH_FRONT => NORTH, TOP, NORTH_FRONT
  work = new vector<int>(3,0);
  (*work)[0] = GEOPOS_NORTH;
  (*work)[1] = GEOPOS_FRONT;  
  (*work)[2] = GEOPOS_NORTH_FRONT;
  m_particleBufferzoneToNeighboringProcs[GEOPOS_NORTH_FRONT] = *work;    
  work->clear();
  delete work;  

  // NORTH_BEHIND => NORTH, BEHIND, NORTH_BEHIND
  work = new vector<int>(3,0);
  (*work)[0] = GEOPOS_NORTH;
  (*work)[1] = GEOPOS_BEHIND;  
  (*work)[2] = GEOPOS_NORTH_BEHIND;
  m_particleBufferzoneToNeighboringProcs[GEOPOS_NORTH_BEHIND] = *work;    
  work->clear();
  delete work; 

  // NORTH_EAST_FRONT => NORTH, EAST, TOP, EAST_FRONT, NORTH_EAST, NORTH_FRONT,
  // NORTH_EAST_FRONT
  work = new vector<int>(7,0);
  (*work)[0] = GEOPOS_NORTH;
  (*work)[1] = GEOPOS_EAST;  
  (*work)[2] = GEOPOS_FRONT;
  (*work)[3] = GEOPOS_EAST_FRONT;  
  (*work)[4] = GEOPOS_NORTH_EAST;  
  (*work)[5] = GEOPOS_NORTH_FRONT;  
  (*work)[6] = GEOPOS_NORTH_EAST_FRONT;  
  m_particleBufferzoneToNeighboringProcs[GEOPOS_NORTH_EAST_FRONT] = *work;    
  work->clear();
  delete work; 

  // NORTH_EAST_BEHIND => NORTH, EAST, BEHIND, EAST_BEHIND, NORTH_EAST, 
  // NORTH_BEHIND, NORTH_EAST_BEHIND
  work = new vector<int>(7,0);
  (*work)[0] = GEOPOS_NORTH;
  (*work)[1] = GEOPOS_EAST;  
  (*work)[2] = GEOPOS_BEHIND;
  (*work)[3] = GEOPOS_EAST_BEHIND;  
  (*work)[4] = GEOPOS_NORTH_EAST;  
  (*work)[5] = GEOPOS_NORTH_BEHIND;  
  (*work)[6] = GEOPOS_NORTH_EAST_BEHIND;  
  m_particleBufferzoneToNeighboringProcs[GEOPOS_NORTH_EAST_BEHIND] = *work;    
  work->clear();
  delete work; 

  // NORTH_WEST_FRONT => NORTH, WEST, TOP, WEST_FRONT, NORTH_WEST, NORTH_FRONT,
  // NORTH_WEST_FRONT
  work = new vector<int>(7,0);
  (*work)[0] = GEOPOS_NORTH;
  (*work)[1] = GEOPOS_WEST;  
  (*work)[2] = GEOPOS_FRONT;
  (*work)[3] = GEOPOS_WEST_FRONT;  
  (*work)[4] = GEOPOS_NORTH_WEST;  
  (*work)[5] = GEOPOS_NORTH_FRONT;  
  (*work)[6] = GEOPOS_NORTH_WEST_FRONT;  
  m_particleBufferzoneToNeighboringProcs[GEOPOS_NORTH_WEST_FRONT] = *work;    
  work->clear();
  delete work; 

  // NORTH_WEST_BEHIND => NORTH, WEST, BEHIND, WEST_BEHIND, NORTH_WEST, 
  // NORTH_BEHIND, NORTH_WEST_BEHIND
  work = new vector<int>(7,0);
  (*work)[0] = GEOPOS_NORTH;
  (*work)[1] = GEOPOS_WEST;  
  (*work)[2] = GEOPOS_BEHIND;
  (*work)[3] = GEOPOS_WEST_BEHIND;  
  (*work)[4] = GEOPOS_NORTH_WEST;  
  (*work)[5] = GEOPOS_NORTH_BEHIND;  
  (*work)[6] = GEOPOS_NORTH_WEST_BEHIND;  
  m_particleBufferzoneToNeighboringProcs[GEOPOS_NORTH_WEST_BEHIND] = *work;    
  work->clear();
  delete work; 

  // SOUTH => SOUTH
  work = new vector<int>(1,0);
  (*work)[0] = GEOPOS_SOUTH;
  m_particleBufferzoneToNeighboringProcs[GEOPOS_SOUTH] = *work;
  work->clear();
  delete work;
  
  // SOUTH_EAST => SOUTH, EAST, SOUTH_EAST
  work = new vector<int>(3,0);
  (*work)[0] = GEOPOS_SOUTH;
  (*work)[1] = GEOPOS_EAST;  
  (*work)[2] = GEOPOS_SOUTH_EAST;
  m_particleBufferzoneToNeighboringProcs[GEOPOS_SOUTH_EAST] = *work;    
  work->clear();
  delete work;
  
  // SOUTH_WEST => SOUTH, WEST, SOUTH_WEST
  work = new vector<int>(3,0);
  (*work)[0] = GEOPOS_SOUTH;
  (*work)[1] = GEOPOS_WEST;  
  (*work)[2] = GEOPOS_SOUTH_WEST;
  m_particleBufferzoneToNeighboringProcs[GEOPOS_SOUTH_WEST] = *work;    
  work->clear();
  delete work;    

  // SOUTH_FRONT => SOUTH, TOP, SOUTH_FRONT
  work = new vector<int>(3,0);
  (*work)[0] = GEOPOS_SOUTH;
  (*work)[1] = GEOPOS_FRONT;  
  (*work)[2] = GEOPOS_SOUTH_FRONT;
  m_particleBufferzoneToNeighboringProcs[GEOPOS_SOUTH_FRONT] = *work;    
  work->clear();
  delete work;  

  // SOUTH_BEHIND => SOUTH, BEHIND, SOUTH_BEHIND
  work = new vector<int>(3,0);
  (*work)[0] = GEOPOS_SOUTH;
  (*work)[1] = GEOPOS_BEHIND;  
  (*work)[2] = GEOPOS_SOUTH_BEHIND;
  m_particleBufferzoneToNeighboringProcs[GEOPOS_SOUTH_BEHIND] = *work;    
  work->clear();
  delete work; 

  // SOUTH_EAST_FRONT => SOUTH, EAST, TOP, EAST_FRONT, SOUTH_EAST, SOUTH_FRONT,
  // SOUTH_EAST_FRONT
  work = new vector<int>(7,0);
  (*work)[0] = GEOPOS_SOUTH;
  (*work)[1] = GEOPOS_EAST;  
  (*work)[2] = GEOPOS_FRONT;
  (*work)[3] = GEOPOS_EAST_FRONT;  
  (*work)[4] = GEOPOS_SOUTH_EAST;  
  (*work)[5] = GEOPOS_SOUTH_FRONT;  
  (*work)[6] = GEOPOS_SOUTH_EAST_FRONT;  
  m_particleBufferzoneToNeighboringProcs[GEOPOS_SOUTH_EAST_FRONT] = *work;    
  work->clear();
  delete work; 

  // SOUTH_EAST_BEHIND => SOUTH, EAST, BEHIND, EAST_BEHIND, SOUTH_EAST, 
  // SOUTH_BEHIND, SOUTH_EAST_BEHIND
  work = new vector<int>(7,0);
  (*work)[0] = GEOPOS_SOUTH;
  (*work)[1] = GEOPOS_EAST;  
  (*work)[2] = GEOPOS_BEHIND;
  (*work)[3] = GEOPOS_EAST_BEHIND;  
  (*work)[4] = GEOPOS_SOUTH_EAST;  
  (*work)[5] = GEOPOS_SOUTH_BEHIND;  
  (*work)[6] = GEOPOS_SOUTH_EAST_BEHIND;  
  m_particleBufferzoneToNeighboringProcs[GEOPOS_SOUTH_EAST_BEHIND] = *work;    
  work->clear();
  delete work; 

  // SOUTH_WEST_FRONT => SOUTH, WEST, TOP, WEST_FRONT, SOUTH_WEST, SOUTH_FRONT,
  // SOUTH_WEST_FRONT
  work = new vector<int>(7,0);
  (*work)[0] = GEOPOS_SOUTH;
  (*work)[1] = GEOPOS_WEST;  
  (*work)[2] = GEOPOS_FRONT;
  (*work)[3] = GEOPOS_WEST_FRONT;  
  (*work)[4] = GEOPOS_SOUTH_WEST;  
  (*work)[5] = GEOPOS_SOUTH_FRONT;  
  (*work)[6] = GEOPOS_SOUTH_WEST_FRONT;  
  m_particleBufferzoneToNeighboringProcs[GEOPOS_SOUTH_WEST_FRONT] = *work;    
  work->clear();
  delete work; 

  // SOUTH_WEST_BEHIND => SOUTH, WEST, BEHIND, WEST_BEHIND, SOUTH_WEST, 
  // SOUTH_BEHIND, SOUTH_WEST_BEHIND
  work = new vector<int>(7,0);
  (*work)[0] = GEOPOS_SOUTH;
  (*work)[1] = GEOPOS_WEST;  
  (*work)[2] = GEOPOS_BEHIND;
  (*work)[3] = GEOPOS_WEST_BEHIND;  
  (*work)[4] = GEOPOS_SOUTH_WEST;  
  (*work)[5] = GEOPOS_SOUTH_BEHIND;  
  (*work)[6] = GEOPOS_SOUTH_WEST_BEHIND;  
  m_particleBufferzoneToNeighboringProcs[GEOPOS_SOUTH_WEST_BEHIND] = *work;    
  work->clear();
  delete work; 
  
  // EAST => EAST
  work = new vector<int>(1,0);
  (*work)[0] = GEOPOS_EAST;
  m_particleBufferzoneToNeighboringProcs[GEOPOS_EAST] = *work;
  work->clear();
  delete work;  

  // WEST => WEST
  work = new vector<int>(1,0);
  (*work)[0] = GEOPOS_WEST;
  m_particleBufferzoneToNeighboringProcs[GEOPOS_WEST] = *work;
  work->clear();
  delete work;  

  // EAST_FRONT => EAST, TOP, EAST_FRONT
  work = new vector<int>(3,0);
  (*work)[0] = GEOPOS_EAST;
  (*work)[1] = GEOPOS_FRONT;  
  (*work)[2] = GEOPOS_EAST_FRONT;
  m_particleBufferzoneToNeighboringProcs[GEOPOS_EAST_FRONT] = *work;    
  work->clear();
  delete work;  

  // EAST_BEHIND => EAST, BEHIND, EAST_BEHIND
  work = new vector<int>(3,0);
  (*work)[0] = GEOPOS_EAST;
  (*work)[1] = GEOPOS_BEHIND;  
  (*work)[2] = GEOPOS_EAST_BEHIND;
  m_particleBufferzoneToNeighboringProcs[GEOPOS_EAST_BEHIND] = *work;    
  work->clear();
  delete work; 

  // WEST_FRONT => WEST, TOP, WEST_FRONT
  work = new vector<int>(3,0);
  (*work)[0] = GEOPOS_WEST;
  (*work)[1] = GEOPOS_FRONT;  
  (*work)[2] = GEOPOS_WEST_FRONT;
  m_particleBufferzoneToNeighboringProcs[GEOPOS_WEST_FRONT] = *work;    
  work->clear();
  delete work;  

  // WEST_BEHIND => WEST, BEHIND, WEST_BEHIND
  work = new vector<int>(3,0);
  (*work)[0] = GEOPOS_WEST;
  (*work)[1] = GEOPOS_BEHIND;  
  (*work)[2] = GEOPOS_WEST_BEHIND;
  m_particleBufferzoneToNeighboringProcs[GEOPOS_WEST_BEHIND] = *work;    
  work->clear();
  delete work; 
  
  // TOP => TOP
  work = new vector<int>(1,0);
  (*work)[0] = GEOPOS_FRONT;
  m_particleBufferzoneToNeighboringProcs[GEOPOS_FRONT] = *work;
  work->clear();
  delete work;  

  // BEHIND => BEHIND
  work = new vector<int>(1,0);
  (*work)[0] = GEOPOS_BEHIND;
  m_particleBufferzoneToNeighboringProcs[GEOPOS_BEHIND] = *work;
  work->clear();
  delete work;  
  
  
  // Reciprocal Geolocalisation
  m_GeoLocReciprocity.reserve(26);
  for (int i=0;i<26;++i) m_GeoLocReciprocity.push_back( 0 );
  m_GeoLocReciprocity[GEOPOS_NORTH] = GEOPOS_SOUTH ;  
  m_GeoLocReciprocity[GEOPOS_NORTH_EAST] = GEOPOS_SOUTH_WEST ;    
  m_GeoLocReciprocity[GEOPOS_NORTH_WEST] = GEOPOS_SOUTH_EAST ;    
  m_GeoLocReciprocity[GEOPOS_NORTH_FRONT] = GEOPOS_SOUTH_BEHIND ;    
  m_GeoLocReciprocity[GEOPOS_NORTH_BEHIND] = GEOPOS_SOUTH_FRONT ;    
  m_GeoLocReciprocity[GEOPOS_NORTH_EAST_FRONT] = GEOPOS_SOUTH_WEST_BEHIND ;    
  m_GeoLocReciprocity[GEOPOS_NORTH_EAST_BEHIND] = GEOPOS_SOUTH_WEST_FRONT ;  
  m_GeoLocReciprocity[GEOPOS_NORTH_WEST_FRONT] = GEOPOS_SOUTH_EAST_BEHIND ;    
  m_GeoLocReciprocity[GEOPOS_NORTH_WEST_BEHIND] = GEOPOS_SOUTH_EAST_FRONT ;    
  m_GeoLocReciprocity[GEOPOS_SOUTH] = GEOPOS_NORTH ;    
  m_GeoLocReciprocity[GEOPOS_SOUTH_EAST] = GEOPOS_NORTH_WEST ;    
  m_GeoLocReciprocity[GEOPOS_SOUTH_WEST] = GEOPOS_NORTH_EAST ;    
  m_GeoLocReciprocity[GEOPOS_SOUTH_FRONT] = GEOPOS_NORTH_BEHIND ;  
  m_GeoLocReciprocity[GEOPOS_SOUTH_BEHIND] = GEOPOS_NORTH_FRONT ;    
  m_GeoLocReciprocity[GEOPOS_SOUTH_EAST_FRONT] = GEOPOS_NORTH_WEST_BEHIND ;    
  m_GeoLocReciprocity[GEOPOS_SOUTH_EAST_BEHIND] = GEOPOS_NORTH_WEST_FRONT ;    
  m_GeoLocReciprocity[GEOPOS_SOUTH_WEST_FRONT] = GEOPOS_NORTH_EAST_BEHIND ;    
  m_GeoLocReciprocity[GEOPOS_SOUTH_WEST_BEHIND] = GEOPOS_NORTH_EAST_FRONT ;    
  m_GeoLocReciprocity[GEOPOS_EAST] = GEOPOS_WEST ;  
  m_GeoLocReciprocity[GEOPOS_WEST] = GEOPOS_EAST ;    
  m_GeoLocReciprocity[GEOPOS_EAST_FRONT] = GEOPOS_WEST_BEHIND ;    
  m_GeoLocReciprocity[GEOPOS_EAST_BEHIND] = GEOPOS_WEST_FRONT ;    
  m_GeoLocReciprocity[GEOPOS_WEST_FRONT] = GEOPOS_EAST_BEHIND ;    
  m_GeoLocReciprocity[GEOPOS_WEST_BEHIND] = GEOPOS_EAST_FRONT ;    
  m_GeoLocReciprocity[GEOPOS_FRONT] = GEOPOS_BEHIND ;
  m_GeoLocReciprocity[GEOPOS_BEHIND] = GEOPOS_FRONT ;    
} 




// ----------------------------------------------------------------------------
// Returns the GeoPosition as a function of the relative position in the 
// MPI cartesian topology
GeoPosition GrainsMPIWrapper::getGeoPosition(int i,int j,int k)
{
  GeoPosition geoLoc = GEOPOS_NONE;
  switch( i )
  {
    case -1:
      switch( j )
      {
        case -1: 
	  switch( k )
	  {
	    case -1:
	      geoLoc = GEOPOS_SOUTH_WEST_BEHIND;
	      break;
	    case 0:
	      geoLoc = GEOPOS_SOUTH_WEST;	    
	      break;
	    case 1:
	      geoLoc = GEOPOS_SOUTH_WEST_FRONT;	    
	      break;
	  }
	  break;
	  
	case 0:
	  switch( k )
	  {
	    case -1:
	      geoLoc = GEOPOS_WEST_BEHIND;	    
	      break;
	    case 0:
	      geoLoc = GEOPOS_WEST;	    
	      break;
	    case 1:
	      geoLoc = GEOPOS_WEST_FRONT;	    
	      break;
	  }	
	  break;
	  
	case 1:
	  switch( k )
	  {
	    case -1:
	      geoLoc = GEOPOS_NORTH_WEST_BEHIND;	    
	      break;
	    case 0:
	      geoLoc = GEOPOS_NORTH_WEST;	    
	      break;
	    case 1:
	      geoLoc = GEOPOS_NORTH_WEST_FRONT;	    
	      break;
	  }	
	  break;
      }
      break;
        
    case 0:
      switch( j )
      {
        case -1: 
	  switch( k )
	  {
	    case -1:
	      geoLoc = GEOPOS_SOUTH_BEHIND;	    
	      break;
	    case 0:
	      geoLoc = GEOPOS_SOUTH;	    
	      break;
	    case 1:
	      geoLoc = GEOPOS_SOUTH_FRONT;	    
	      break;
	  }
	  break;
	  
	case 0:
	  switch( k )
	  {
	    case -1:
	      geoLoc = GEOPOS_BEHIND;	    
	      break;
	    case 0:
	      geoLoc = GEOPOS_NONE;	    
	      break;
	    case 1:
	      geoLoc = GEOPOS_FRONT;	    
	      break;
	  }	
	  break;
	  
	case 1:
	  switch( k )
	  {
	    case -1:
	      geoLoc = GEOPOS_NORTH_BEHIND;	    
	      break;
	    case 0:
	      geoLoc = GEOPOS_NORTH;	    
	      break;
	    case 1:
	      geoLoc = GEOPOS_NORTH_FRONT;	    
	      break;
	  }	
	  break;
      }    
      break;
          
    case 1:
      switch( j )
      {
        case -1: 
	  switch( k )
	  {
	    case -1:
	      geoLoc = GEOPOS_SOUTH_EAST_BEHIND;	    
	      break;
	    case 0:
	      geoLoc = GEOPOS_SOUTH_EAST;	    
	      break;
	    case 1:
	      geoLoc = GEOPOS_SOUTH_EAST_FRONT;	    
	      break;
	  }
	  break;
	  
	case 0:
	  switch( k )
	  {
	    case -1:
	      geoLoc = GEOPOS_EAST_BEHIND;	    
	      break;
	    case 0:
	      geoLoc = GEOPOS_EAST;	    
	      break;
	    case 1:
	      geoLoc = GEOPOS_EAST_FRONT;	    
	      break;
	  }	
	  break;
	  
	case 1:
	  switch( k )
	  {
	    case -1:
	      geoLoc = GEOPOS_NORTH_EAST_BEHIND;	    
	      break;
	    case 0:
	      geoLoc = GEOPOS_NORTH_EAST;	    
	      break;
	    case 1:
	      geoLoc = GEOPOS_NORTH_EAST_FRONT;	    
	      break;
	  }	
	  break;
      }    
      break;    
  }
  
  return ( geoLoc );
}   




// ----------------------------------------------------------------------------
// Sets periodicity vectors
void GrainsMPIWrapper::setMPIperiodicVectors( const double& lx, 
	const double& ly, const double& lz )
{
  // West
  if ( m_neighbors->rank( -1, 0, 0 ) != -1 )
  {
    if ( m_period[0] && m_coords[0] == 0 )
    {
      m_MPIperiodes[GEOPOS_WEST][X] = lx ;
      m_MPIperiodes[GEOPOS_WEST][Y] = 0. ;
      m_MPIperiodes[GEOPOS_WEST][Z] = 0. ;
    }  
  }

  // East
  if ( m_neighbors->rank( 1, 0, 0 ) != -1 )
  {
    if ( m_period[0] && m_coords[0] == m_dim[0] - 1 )
    {
      m_MPIperiodes[GEOPOS_EAST][X] = -lx ;
      m_MPIperiodes[GEOPOS_EAST][Y] = 0. ;
      m_MPIperiodes[GEOPOS_EAST][Z] = 0. ;
    }  
  }

  // South
  if ( m_neighbors->rank( 0, -1, 0 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == 0 )
    {
      m_MPIperiodes[GEOPOS_SOUTH][X] = 0. ;
      m_MPIperiodes[GEOPOS_SOUTH][Y] = ly ;
      m_MPIperiodes[GEOPOS_SOUTH][Z] = 0. ;
    }  
  }

  // North
  if ( m_neighbors->rank( 0, 1, 0 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == m_dim[1] - 1 )
    {
      m_MPIperiodes[GEOPOS_NORTH][X] = 0. ;
      m_MPIperiodes[GEOPOS_NORTH][Y] = -ly ;
      m_MPIperiodes[GEOPOS_NORTH][Z] = 0. ;
    }  
  }

  // Bottom
  if ( m_neighbors->rank( 0, 0, -1 ) != -1 )
  {
    if ( m_period[2] && m_coords[2] == 0 )
    {
      m_MPIperiodes[GEOPOS_BEHIND][X] = 0. ;
      m_MPIperiodes[GEOPOS_BEHIND][Y] = 0. ;
      m_MPIperiodes[GEOPOS_BEHIND][Z] = lz ;
    }  
  }

  // Top
  if ( m_neighbors->rank( 0, 0, 1 ) != -1 )
  {
    if ( m_period[2] && m_coords[2] == m_dim[2] - 1 )
    {
      m_MPIperiodes[GEOPOS_FRONT][X] = 0. ;
      m_MPIperiodes[GEOPOS_FRONT][Y] = 0. ;
      m_MPIperiodes[GEOPOS_FRONT][Z] = -lz ;
    }  
  }


  // South West
  if ( m_neighbors->rank( -1, -1, 0 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == 0 )
    {
      m_MPIperiodes[GEOPOS_SOUTH_WEST][X] = 0. ;
      m_MPIperiodes[GEOPOS_SOUTH_WEST][Y] = ly ;
      m_MPIperiodes[GEOPOS_SOUTH_WEST][Z] = 0. ;
    } 

    if ( m_period[0] && m_coords[0] == 0 )
      m_MPIperiodes[GEOPOS_SOUTH_WEST][X] += lx ;
  }
    
  // South East
  if ( m_neighbors->rank( 1, -1, 0 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == 0 )
    {
      m_MPIperiodes[GEOPOS_SOUTH_EAST][X] = 0. ;
      m_MPIperiodes[GEOPOS_SOUTH_EAST][Y] = ly ;
      m_MPIperiodes[GEOPOS_SOUTH_EAST][Z] = 0. ;
    } 

    if ( m_period[0] && m_coords[0] == m_dim[0] - 1 )
      m_MPIperiodes[GEOPOS_SOUTH_EAST][X] += -lx ;
  }  

  // South Bottom
  if ( m_neighbors->rank( 0, -1, -1 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == 0 )
    {
      m_MPIperiodes[GEOPOS_SOUTH_BEHIND][X] = 0. ;
      m_MPIperiodes[GEOPOS_SOUTH_BEHIND][Y] = ly ;
      m_MPIperiodes[GEOPOS_SOUTH_BEHIND][Z] = 0. ;
    } 

    if ( m_period[2] && m_coords[2] == 0 )
      m_MPIperiodes[GEOPOS_SOUTH_BEHIND][Z] += lz ;
  }
    
  // South Top
  if ( m_neighbors->rank( 0, -1, 1 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == 0 )
    {
      m_MPIperiodes[GEOPOS_SOUTH_FRONT][X] = 0. ;
      m_MPIperiodes[GEOPOS_SOUTH_FRONT][Y] = ly ;
      m_MPIperiodes[GEOPOS_SOUTH_FRONT][Z] = 0. ;
    } 

    if ( m_period[2] && m_coords[2] == m_dim[2] - 1 )
      m_MPIperiodes[GEOPOS_SOUTH_FRONT][Z] += -lz ;
  } 

  // North West
  if ( m_neighbors->rank( -1, 1, 0 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == m_dim[1] - 1 )
    {
      m_MPIperiodes[GEOPOS_NORTH_WEST][X] = 0. ;
      m_MPIperiodes[GEOPOS_NORTH_WEST][Y] = -ly ;
      m_MPIperiodes[GEOPOS_NORTH_WEST][Z] = 0. ;
    }
    
    if ( m_period[0] && m_coords[0] == 0 )
      m_MPIperiodes[GEOPOS_NORTH_WEST][X] = +lx ;      
  }
    
  // North East
  if ( m_neighbors->rank( 1, 1, 0 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == m_dim[1] - 1 )
    {
      m_MPIperiodes[GEOPOS_NORTH_EAST][X] = 0. ;
      m_MPIperiodes[GEOPOS_NORTH_EAST][Y] = -ly ;
      m_MPIperiodes[GEOPOS_NORTH_EAST][Z] = 0. ;
    }
      
    if ( m_period[0] && m_coords[0] == m_dim[0] - 1 )
      m_MPIperiodes[GEOPOS_NORTH_EAST][X] += -lx ;      
  }
  
  // North Bottom
  if ( m_neighbors->rank( 0, 1, -1 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == m_dim[1] - 1 )
    {
      m_MPIperiodes[GEOPOS_NORTH_BEHIND][X] = 0. ;
      m_MPIperiodes[GEOPOS_NORTH_BEHIND][Y] = -ly ;
      m_MPIperiodes[GEOPOS_NORTH_BEHIND][Z] = 0. ;
    } 

    if ( m_period[2] && m_coords[2] == 0 )
      m_MPIperiodes[GEOPOS_NORTH_BEHIND][Z] += lz ;
  }
    
  // North Top
  if ( m_neighbors->rank( 0, 1, 1 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == m_dim[1] - 1 )
    {
      m_MPIperiodes[GEOPOS_NORTH_FRONT][X] = 0. ;
      m_MPIperiodes[GEOPOS_NORTH_FRONT][Y] = -ly ;
      m_MPIperiodes[GEOPOS_NORTH_FRONT][Z] = 0. ;
    } 

    if ( m_period[2] && m_coords[2] == m_dim[2] - 1 )
      m_MPIperiodes[GEOPOS_NORTH_FRONT][Z] += -lz ;
  }    

  // West Bottom
  if ( m_neighbors->rank( -1, 0, -1 ) != -1 )
  {
    if ( m_period[0] && m_coords[0] == 0 )
    {
      m_MPIperiodes[GEOPOS_WEST_BEHIND][X] = lx ;
      m_MPIperiodes[GEOPOS_WEST_BEHIND][Y] = 0. ;
      m_MPIperiodes[GEOPOS_WEST_BEHIND][Z] = 0. ;
    } 
    
    if ( m_period[2] && m_coords[2] == 0 )
      m_MPIperiodes[GEOPOS_WEST_BEHIND][Z] += lz ;     
  }
  
  // West Top
  if ( m_neighbors->rank( -1, 0, 1 ) != -1 )
  {
    if ( m_period[0] && m_coords[0] == 0 )
    {
      m_MPIperiodes[GEOPOS_WEST_FRONT][X] = lx ;
      m_MPIperiodes[GEOPOS_WEST_FRONT][Y] = 0. ;
      m_MPIperiodes[GEOPOS_WEST_FRONT][Z] = 0. ;
    } 
    
    if ( m_period[2] && m_coords[2] == m_dim[2] - 1 )
      m_MPIperiodes[GEOPOS_WEST_FRONT][Z] += -lz ;     
  }  

  // East Bottom
  if ( m_neighbors->rank( 1, 0, -1 ) != -1 )
  {
    if ( m_period[0] && m_coords[0] == m_dim[0] - 1 )
    {
      m_MPIperiodes[GEOPOS_EAST_BEHIND][X] = -lx ;
      m_MPIperiodes[GEOPOS_EAST_BEHIND][Y] = 0. ;
      m_MPIperiodes[GEOPOS_EAST_BEHIND][Z] = 0. ;
    } 
    
    if ( m_period[2] && m_coords[2] == 0 )
      m_MPIperiodes[GEOPOS_EAST_BEHIND][Z] += lz ;     
  }
  
  // East Top
  if ( m_neighbors->rank( 1, 0, 1 ) != -1 )
  {
    if ( m_period[0] && m_coords[0] == m_dim[0] - 1 )
    {
      m_MPIperiodes[GEOPOS_EAST_FRONT][X] = -lx ;
      m_MPIperiodes[GEOPOS_EAST_FRONT][Y] = 0. ;
      m_MPIperiodes[GEOPOS_EAST_FRONT][Z] = 0. ;
    } 
    
    if ( m_period[2] && m_coords[2] == m_dim[2] - 1 )
      m_MPIperiodes[GEOPOS_EAST_FRONT][Z] += -lz ;     
  }


  // South West Bottom
  if ( m_neighbors->rank( -1, -1, -1 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == 0 )
    {
      m_MPIperiodes[GEOPOS_SOUTH_WEST_BEHIND][X] = 0. ;
      m_MPIperiodes[GEOPOS_SOUTH_WEST_BEHIND][Y] = ly ;
      m_MPIperiodes[GEOPOS_SOUTH_WEST_BEHIND][Z] = 0. ;
    } 

    if ( m_period[0] && m_coords[0] == 0 )
      m_MPIperiodes[GEOPOS_SOUTH_WEST_BEHIND][X] += lx ;
      
    if ( m_period[2] && m_coords[2] == 0 )
      m_MPIperiodes[GEOPOS_SOUTH_WEST_BEHIND][Z] += lz ;       
  }  

  // South West Top
  if ( m_neighbors->rank( -1, -1, 1 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == 0 )
    {
      m_MPIperiodes[GEOPOS_SOUTH_WEST_FRONT][X] = 0. ;
      m_MPIperiodes[GEOPOS_SOUTH_WEST_FRONT][Y] = ly ;
      m_MPIperiodes[GEOPOS_SOUTH_WEST_FRONT][Z] = 0. ;
    } 

    if ( m_period[0] && m_coords[0] == 0 )
      m_MPIperiodes[GEOPOS_SOUTH_WEST_FRONT][X] += lx ;
      
    if ( m_period[2] && m_coords[2] == m_dim[2] - 1 )
      m_MPIperiodes[GEOPOS_SOUTH_WEST_FRONT][Z] += -lz ;        
  }

  // North West Bottom
  if ( m_neighbors->rank( -1, 1, -1 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == m_dim[1] - 1 )
    {
      m_MPIperiodes[GEOPOS_NORTH_WEST_BEHIND][X] = 0. ;
      m_MPIperiodes[GEOPOS_NORTH_WEST_BEHIND][Y] = -ly ;
      m_MPIperiodes[GEOPOS_NORTH_WEST_BEHIND][Z] = 0. ;
    } 

    if ( m_period[0] && m_coords[0] == 0 )
      m_MPIperiodes[GEOPOS_NORTH_WEST_BEHIND][X] += lx ;
      
    if ( m_period[2] && m_coords[2] == 0 )
      m_MPIperiodes[GEOPOS_NORTH_WEST_BEHIND][Z] += lz ;       
  }  

  // North West Top
  if ( m_neighbors->rank( -1, 1, 1 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == m_dim[1] - 1 )
    {
      m_MPIperiodes[GEOPOS_NORTH_WEST_FRONT][X] = 0. ;
      m_MPIperiodes[GEOPOS_NORTH_WEST_FRONT][Y] = -ly ;
      m_MPIperiodes[GEOPOS_NORTH_WEST_FRONT][Z] = 0. ;
    } 

    if ( m_period[0] && m_coords[0] == 0 )
      m_MPIperiodes[GEOPOS_NORTH_WEST_FRONT][X] += lx ;
      
    if ( m_period[2] && m_coords[2] == m_dim[2] - 1 )
      m_MPIperiodes[GEOPOS_NORTH_WEST_FRONT][Z] += -lz ;        
  }

  // South East Bottom
  if ( m_neighbors->rank( 1, -1, -1 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == 0 )
    {
      m_MPIperiodes[GEOPOS_SOUTH_EAST_BEHIND][X] = 0. ;
      m_MPIperiodes[GEOPOS_SOUTH_EAST_BEHIND][Y] = ly ;
      m_MPIperiodes[GEOPOS_SOUTH_EAST_BEHIND][Z] = 0. ;
    } 

    if ( m_period[0] && m_coords[0] == m_dim[0] - 1 )
      m_MPIperiodes[GEOPOS_SOUTH_EAST_BEHIND][X] += -lx ;
      
    if ( m_period[2] && m_coords[2] == 0 )
      m_MPIperiodes[GEOPOS_SOUTH_EAST_BEHIND][Z] += lz ;       
  }  

  // South East Top
  if ( m_neighbors->rank( 1, -1, 1 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == 0 )
    {
      m_MPIperiodes[GEOPOS_SOUTH_EAST_FRONT][X] = 0. ;
      m_MPIperiodes[GEOPOS_SOUTH_EAST_FRONT][Y] = ly ;
      m_MPIperiodes[GEOPOS_SOUTH_EAST_FRONT][Z] = 0. ;
    } 

    if ( m_period[0] && m_coords[0] == m_dim[0] - 1 )
      m_MPIperiodes[GEOPOS_SOUTH_EAST_FRONT][X] += -lx ;
      
    if ( m_period[2] && m_coords[2] == m_dim[2] - 1 )
      m_MPIperiodes[GEOPOS_SOUTH_EAST_FRONT][Z] += -lz ;        
  }

  // North East Bottom
  if ( m_neighbors->rank( 1, 1, -1 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == m_dim[1] - 1 )
    {
      m_MPIperiodes[GEOPOS_NORTH_EAST_BEHIND][X] = 0. ;
      m_MPIperiodes[GEOPOS_NORTH_EAST_BEHIND][Y] = -ly ;
      m_MPIperiodes[GEOPOS_NORTH_EAST_BEHIND][Z] = 0. ;
    } 

    if ( m_period[0] && m_coords[0] == m_dim[0] - 1 )
      m_MPIperiodes[GEOPOS_NORTH_EAST_BEHIND][X] += -lx ;
      
    if ( m_period[2] && m_coords[2] == 0 )
      m_MPIperiodes[GEOPOS_NORTH_EAST_BEHIND][Z] += lz ;       
  }  

  // North East Top
  if ( m_neighbors->rank( 1, 1, 1 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == m_dim[1] - 1 )
    {
      m_MPIperiodes[GEOPOS_NORTH_EAST_FRONT][X] = 0. ;
      m_MPIperiodes[GEOPOS_NORTH_EAST_FRONT][Y] = -ly ;
      m_MPIperiodes[GEOPOS_NORTH_EAST_FRONT][Z] = 0. ;
    } 

    if ( m_period[0] && m_coords[0] == m_dim[0] - 1 )
      m_MPIperiodes[GEOPOS_NORTH_EAST_FRONT][X] += -lx ;
      
    if ( m_period[2] && m_coords[2] == m_dim[2] - 1 )
      m_MPIperiodes[GEOPOS_NORTH_EAST_FRONT][Z] += -lz ;        
  }
}	 




// ----------------------------------------------------------------------------
// Returns the number of processes in one direction
int GrainsMPIWrapper::get_nb_procs_direction( int i ) const
{
  return ( m_dim[i] );
}




// ----------------------------------------------------------------------------
// Returns the number of processes in all directions
int const* GrainsMPIWrapper::get_nb_procs_direction() const
{
  return ( m_dim );
}




// ----------------------------------------------------------------------------
// Returns the total number of processes in the MPI_COMM_WORLD communicator
int GrainsMPIWrapper::get_total_number_of_processes() const
{
  return ( m_nprocs_world );
}




// ----------------------------------------------------------------------------
// Returns the total number of active processes in the MPI_COMM_activProc 
// communicator
int GrainsMPIWrapper::get_total_number_of_active_processes() const
{
  return ( m_nprocs );
}




// ----------------------------------------------------------------------------
// Returns the process rank in all communicators
int GrainsMPIWrapper::get_rank() const
{
  return ( m_rank );
}




// ----------------------------------------------------------------------------
// Returns whether the process is active
bool GrainsMPIWrapper::isActive() const
{
  return ( m_is_active );
} 




// ----------------------------------------------------------------------------
// Returns the MPI cartesian coordinates of the process
int const* GrainsMPIWrapper::get_MPI_coordinates() const
{
  return ( m_coords );
}




// ----------------------------------------------------------------------------
// Returns the MPINeighbors of the process
MPINeighbors const* GrainsMPIWrapper::get_MPI_neighbors() const
{
  return ( m_neighbors );
}




// ----------------------------------------------------------------------------
// Returns the periodicity of the domain that is the same as the
// MPI periodicity of the MPI cartesian topology
int const* GrainsMPIWrapper::get_MPI_periodicity() const
{
  return ( m_period );
}




// ----------------------------------------------------------------------------
// Writes the MPI wrapper features in a stream
void GrainsMPIWrapper::display( ostream& f, string const& oshift ) const
{
  if ( m_rank == m_rank_master )
  {
    f << oshift << "Domain decomposition = ";
    for (int j=0;j<3;++j) f << "N[" << j << "]=" << m_dim[j] << " ";
    f << endl;
    f << oshift << "MPI periods = ";
    for (int j=0;j<3;++j) f << "P[" << j << "]=" << m_period[j] << " ";
    f << endl;    
  }

  ostringstream out;
  out << oshift << GrainsExec::m_shift3 << "Process PID = " << getpid() << endl;
  out << oshift << GrainsExec::m_shift3 << "Position in the MPI cartesian"
      	" topology = ";
  for (int j=0;j<3;++j) out << m_coords[j] << " ";
  out << endl;
  out << oshift << GrainsExec::m_shift3 << "Neighboring processes in the MPI "
      	"cartesian topology" << endl;
  out << oshift << GrainsExec::m_shift3 
      	<< "---------------------------------------------------";
  for (int i=-1;i<2;i++)
    for (int j=-1;j<2;j++)
      for (int k=-1;k<2;k++)
        if ( m_neighbors->rank( i, j, k ) != -1 )
        {  
          out << endl << oshift << GrainsExec::m_shift3 << "Neighbor (" << i 
	  	<< "," << j << "," << k << ") GEOLOC = " << 
		Cell::getGeoPositionName( getGeoPosition( i, j, k ) ) << endl;
	  int const* coords_ = m_neighbors->coordinates( i, j, k );
	  out << oshift << GrainsExec::m_shift3 
	      	<< "Position in MPI topology = " << coords_[0] << " "
		<< coords_[1] << " " << coords_[2] << endl;
	  out << oshift << GrainsExec::m_shift3 << "Rank in MPI topology = " 
	      	<< m_neighbors->rank( i, j, k ) << endl;
	  out << oshift << GrainsExec::m_shift3 << "MPI periodic vector = " 
	      	<< m_MPIperiodes[ getGeoPosition( i, j, k ) ][X] << " " << 
	      	m_MPIperiodes[ getGeoPosition( i, j, k ) ][Y] << " " << 
	      	m_MPIperiodes[ getGeoPosition( i, j, k ) ][Z];	
        }

  writeStringPerProcess( f, out.str(), true, oshift + GrainsExec::m_shift3 );
}




// ----------------------------------------------------------------------------
// Gathers all particle data on the master process for post-processing purposes
vector< vector<double> >* GrainsMPIWrapper::
    GatherParticleData_PostProcessing( list<Particle*> const& particles,
	size_t const& nb_total_particles ) const
{
  vector< vector<double> >* data_Global = NULL;
  int NB_DOUBLE_PART = 11, recvsize = 0;
  list<Particle*>::const_iterator il;
  int i = 0, nb_part_loc = int( particles.size() ), id = 0;  
  double intTodouble = 0.1 ;  
  MPI_Status status;
  MPI_Request idreq;
  
  
  // Allocate the receiving vector on master process
  if ( m_rank == m_rank_master )
  { 
    vector<double> work( NB_DOUBLE_PART - 1, 0. );
    data_Global = new vector< vector<double> >( nb_total_particles, work );
  }
    
  // Exclude clone particles
  for (il=particles.begin();il!=particles.end();il++)
    if ( (*il)->getTag() == 2 ) nb_part_loc--;
    
  // Buffer size depend on the number of particles per core
  double* buffer = new double[ NB_DOUBLE_PART * nb_part_loc ];
  
  for (il=particles.begin(), i=0; il!=particles.end(); il++)
  {
    if( (*il)->getTag() != 2 )
    {
      buffer[i] = double((*il)->getID()) + intTodouble;
      (*il)->copyPosition( buffer, i+1 );
      (*il)->copyTranslationalVelocity( buffer, i+4 );
      (*il)->copyAngularVelocity( buffer, i+7 );    
      buffer[i+10] = double((*il)->getCoordinationNumber()) + intTodouble;
      i += 11; 
    }
  }

  // Process sends buffer to the master process
  MPI_Isend( buffer, NB_DOUBLE_PART * nb_part_loc, MPI_DOUBLE, 
  	m_rank_master, m_tag_DOUBLE, m_MPI_COMM_activeProc, &idreq );           

  // Reception by the master process
  if ( m_rank == m_rank_master )
  {
    // Allocate the receiving vector on master process
    vector<double> work( NB_DOUBLE_PART, 0. );
    data_Global = new vector< vector<double> >( nb_total_particles, work );

    for (int irank=0; irank<m_nprocs; ++irank)
    {
      // Size of the message
      MPI_Probe( irank, m_tag_DOUBLE, m_MPI_COMM_activeProc, &status );  
      MPI_Get_count( &status, MPI_DOUBLE, &recvsize );

      // Reception of the actual message	
      double* recvbuf_DOUBLE = new double[recvsize];
      MPI_Recv( recvbuf_DOUBLE, recvsize, MPI_DOUBLE, 
          irank, m_tag_DOUBLE, m_MPI_COMM_activeProc, &status );	    
      
      // Copy in data_Global
      for (int j=0; j<recvsize; j+=NB_DOUBLE_PART)
      {
        id = int(recvbuf_DOUBLE[j]);
	for (int k=1;k<NB_DOUBLE_PART;k++)
	  (*data_Global)[id][k-1] = recvbuf_DOUBLE[j+k]; 
      }	
      
      delete [] recvbuf_DOUBLE; 
    }
  }

  // Verify that all non-blocking sends completed
  MPI_Wait( &idreq, &status );    

  delete [] buffer ;
  
  return ( data_Global );
}




// ----------------------------------------------------------------------------
// Gathers the class of all particles on the master process 
vector<int>* GrainsMPIWrapper::
    GatherParticlesClass_PostProcessing(
    const list<Particle*> &particles,
    const size_t& nb_total_particles ) const
{
  int i=0, recvsize = 0;
  MPI_Status status;
  MPI_Request idreq;
  vector<int>* class_Global = NULL;    
  list<Particle*>::const_iterator il;
  int nb_part_loc = int( particles.size() );
  
  // Exclude clone particles
  for (il=particles.begin();il!=particles.end();il++)
    if ( (*il)->getTag() == 2 ) nb_part_loc--;

  // Buffer size depend on the number of particles per core
  int* buffer = new int[ 2 * nb_part_loc ];
  
  for (il=particles.begin(), i=0; il!=particles.end(); il++)
  {
    if ( (*il)->getTag() != 2 )
    {
      buffer[i] = (*il)->getID();
      buffer[i+1] = (*il)->getGeometricType();
      i += 2;
    }
  }

  // Process sends buffer to the master process
  MPI_Isend( buffer, 2 * nb_part_loc, MPI_INT, m_rank_master,
      m_tag_INT, m_MPI_COMM_activeProc, &idreq );	

  // Reception by the master process
  if ( m_rank == m_rank_master )
  {
    class_Global = new vector<int>( nb_total_particles, 0 ) ;

    for (int irank=0; irank<m_nprocs; ++irank)
    {
      // Size of the message
      MPI_Probe( irank, m_tag_INT, m_MPI_COMM_activeProc, &status );  
      MPI_Get_count( &status, MPI_INT, &recvsize );

      // Reception of the actual message	
      int* recvbuf_INT = new int[recvsize];
      MPI_Recv( recvbuf_INT, recvsize, MPI_INT, 
          irank, m_tag_INT, m_MPI_COMM_activeProc, &status );	    
      
      // Copy in class_Global
      for (int j=0; j<recvsize; j+=2)
        (*class_Global)[recvbuf_INT[j]] = recvbuf_INT[j+1];
      
      delete [] recvbuf_INT; 
    }
  }

  // Verify that all non-blocking sends completed
  MPI_Wait( &idreq, &status );    

  delete [] buffer ;

  return ( class_Global );
  
}  




// ----------------------------------------------------------------------------
// Creates and updates clone particles using a Send-Recv strategy
// with neighboring processes in the MPI cartesian topology
void GrainsMPIWrapper::UpdateOrCreateClones_SendRecvLocal_GeoLoc( double time,
	list<Particle*>* particles,
  	list<Particle*> const* particlesBufferzone,
  	list<Particle*>* particlesClones,
	vector<Particle*> const* referenceParticles,
	LinkedCell* LC, bool update,
	bool update_velocity_only )
{
  list<Particle*>::const_iterator il;
  int i, j, recvsize = 0, geoLoc, ireq = 0,
  	contact_map_size = 0, NB_DOUBLE_PER_CONTACT = 13, NB_DOUBLE_PART = 25;
  MPI_Status status;
  MPI_Request sreq = 0;
  list<int> const* neighborsRank = m_neighbors->rankMPINeighborsOnly();
  list<int>::const_iterator irn;
  list<GeoPosition> const* neighborsGeoloc = 
  	m_neighbors->geolocMPINeighborsOnly();
  list<GeoPosition>::const_iterator ign;  
  vector<MPI_Request> idreq( neighborsRank->size(), sreq );    
  double intTodouble = 0.1 ;
  
  SCT_set_start( "BuffersCopy" );

  // Fill the multimap to access clones by their ID number
  // -----------------------------------------------------
  m_AccessToClones.clear();
  for (il=particlesClones->begin();il!=particlesClones->end();il++)
    m_AccessToClones.insert( pair<int,Particle*>( (*il)->getID(), *il ) );     
          
  // Copy particles in buffer zone into local buffers
  // ------------------------------------------------
  vector<int> nbBufGeoLoc(26,0);
  vector<int>::iterator iv;
  for (il=particlesBufferzone->begin();il!=particlesBufferzone->end();il++)
  {
    geoLoc = (*il)->getGeoPosition();
    contact_map_size = (*il)->getContactMapSize();
    for (iv=m_particleBufferzoneToNeighboringProcs[geoLoc].begin();
    	iv!=m_particleBufferzoneToNeighboringProcs[geoLoc].end();iv++)
      nbBufGeoLoc[*iv] += NB_DOUBLE_PART + 1 
      	+ contact_map_size * NB_DOUBLE_PER_CONTACT;
  }             

  // Buffer of doubles: kinematics and configuration as follows 
  // [ID number, class, rank of sending process, translational
  //  velocity, rotation quaternion, angular velocity, transform]
  // ------------------------------------------------------------
  vector<int> index( 26, 0 );
  double *pDOUBLE = NULL;  
  vector<double*> features( 26, pDOUBLE );
  for (i=0;i<26;i++) features[i] = new double[ nbBufGeoLoc[i] ];
  double ParticleID = 0., ParticleClass = 0.;

  for (il=particlesBufferzone->begin(),i=0;il!=particlesBufferzone->end();il++)
  {
    geoLoc = (*il)->getGeoPosition();
    ParticleID = double((*il)->getID()) + intTodouble ;
    ParticleClass = double((*il)->getGeometricType()) + intTodouble ;
    contact_map_size = (*il)->getContactMapSize();
    for (iv=m_particleBufferzoneToNeighboringProcs[geoLoc].begin();
    	iv!=m_particleBufferzoneToNeighboringProcs[geoLoc].end();iv++)
    {
      j = index[*iv]; 
      features[*iv][j] = ParticleID;             
      features[*iv][j+1] = ParticleClass;
      features[*iv][j+2] = double(m_rank) + intTodouble ;
      (*il)->copyTranslationalVelocity( features[*iv], j+3 );
      (*il)->copyQuaternionRotation( features[*iv], j+6 );    
      (*il)->copyAngularVelocity( features[*iv], j+10 );
      (*il)->copyTransform( features[*iv], j+13, m_MPIperiodes[*iv] );
      (*il)->copyContactMap( features[*iv], j+NB_DOUBLE_PART );     
      index[*iv] += NB_DOUBLE_PART + 1 
      	+ contact_map_size * NB_DOUBLE_PER_CONTACT;
    }
  }           
  SCT_get_elapsed_time("BuffersCopy");
  
  // Communication
  // -------------
  bool first_update = true;

  // Send data to neighboring processes based on their geolocalisation
  SCT_set_start( "MPIComm" );
  for (ireq=0,irn=neighborsRank->begin(),ign=neighborsGeoloc->begin();
  	irn!=neighborsRank->end();irn++,ign++,++ireq)
    MPI_Isend( features[*ign], nbBufGeoLoc[*ign], MPI_DOUBLE, 
	*irn, m_tag_DOUBLE + m_GeoLocReciprocity[*ign], 
	m_MPI_COMM_activeProc, &idreq[ireq] );
  SCT_get_elapsed_time( "MPIComm" );

  // Receive data sent by neighboring processes and processing of these data
  // Reception par le processus des messages envoyés par ses voisins
  for (irn=neighborsRank->begin(),ign=neighborsGeoloc->begin();
  	irn!=neighborsRank->end();irn++,ign++)  
  {
    SCT_set_start( "MPIComm" );
	    
    // Reception
    // ---------
    // Size of the message -> number of particles
    MPI_Probe( *irn, m_tag_DOUBLE + *ign, m_MPI_COMM_activeProc, &status );  
    MPI_Get_count( &status, MPI_DOUBLE, &recvsize );

    // Reception of the actual message	    
    double *recvbuf_DOUBLE = new double[recvsize];
    MPI_Recv( recvbuf_DOUBLE, recvsize, MPI_DOUBLE, 
	*irn, m_tag_DOUBLE + *ign, m_MPI_COMM_activeProc, &status );	    

    SCT_add_elapsed_time( "MPIComm" );
    SCT_set_start( "UpdateCreateClones" );
 
    // Creation or update of clone particles
    // -------------------------------------
    if ( update )    
      UpdateClones( time, recvsize, recvbuf_DOUBLE,
		NB_DOUBLE_PART, NB_DOUBLE_PER_CONTACT, particlesClones,
		particles, particlesBufferzone, referenceParticles, LC,
		update_velocity_only );
    else
      CreateClones( time, recvsize, recvbuf_DOUBLE,
		NB_DOUBLE_PART, NB_DOUBLE_PER_CONTACT, particlesClones,
		particles, particlesBufferzone, referenceParticles, LC );      

    delete [] recvbuf_DOUBLE;

    if ( first_update ) 
    {
      SCT_get_elapsed_time( "UpdateCreateClones" );
      first_update = false;
    }
    else SCT_add_elapsed_time( "UpdateCreateClones" );            
  }         

  // Verify that all non-blocking sends completed
  for (ireq=0;ireq<int(idreq.size());++ireq)
    MPI_Wait( &idreq[ireq], &status );  

  for (vector<double*>::iterator ivpd=features.begin();ivpd!=features.end();
  	ivpd++) delete [] *ivpd;
  features.clear(); 
}




// ----------------------------------------------------------------------------
// Updates clones with the data sent by the neighboring processes 
void GrainsMPIWrapper::UpdateClones(double time,
 	int const &recvsize, double const* recvbuf_DOUBLE,
	const int& NB_DOUBLE_PART,
	int const& NB_DOUBLE_PER_CONTACT,  
  	list<Particle*>* particlesClones,
	list<Particle*>* particles,
  	list<Particle*> const* particlesBufferzone,
	vector<Particle*> const* referenceParticles,
	LinkedCell* LC, bool update_velocity_only )
{
  int j, id, contact_map_size = 0;
  bool found = false;
  double distGC = 0. ; 
  Point3 const* GC = NULL ; 
  multimap<int,Particle*>::iterator imm;
  size_t ncid = 0;
  Particle* pClone = NULL ;
  pair < multimap<int,Particle*>::iterator, 
  	multimap<int,Particle*>::iterator > crange;
  
  for (j=0; j<recvsize; )
  {
    found = false;
    id = int( recvbuf_DOUBLE[j] );
    contact_map_size = int( recvbuf_DOUBLE[j+NB_DOUBLE_PART] );  
      
    // Search whether the clone already exists in this local domain
    ncid = m_AccessToClones.count( id );
    switch( ncid )
    {
      case 0: // no clone with this ID number in this local domain
        break;
      case 1: // 1 clone with this ID number in this local domain
        imm = m_AccessToClones.find( id );
        if ( m_isMPIperiodic )
        {
          GC = imm->second->getPosition();
          distGC = sqrt( 
          	pow( recvbuf_DOUBLE[j+22] - (*GC)[X], 2. ) +
            	pow( recvbuf_DOUBLE[j+23] - (*GC)[Y], 2. ) +
            	pow( recvbuf_DOUBLE[j+24] - (*GC)[Z], 2. ) ) ;
          if ( distGC < 1.1 * imm->second->getCrustThickness() )  
              found = true;	
        }
        else found = true;
        break;
      default: // more than 1 clone with this ID number in this local domain
        // Periodic case with 1 multi-proc clone and 1 periodic clone in the
	// same local domain
        crange = m_AccessToClones.equal_range( id );
        for (imm=crange.first; imm!=crange.second && !found; )
        {
          GC = imm->second->getPosition();
          distGC = sqrt( 
          	pow( recvbuf_DOUBLE[j+22] - (*GC)[X], 2. ) +
          	pow( recvbuf_DOUBLE[j+23] - (*GC)[Y], 2. ) +
          	pow( recvbuf_DOUBLE[j+24] - (*GC)[Z], 2. ) ) ;
          if ( distGC < 1.1 * imm->second->getCrustThickness() )  
            found = true;
          else imm++;
        }
        break;    
    }   
       
    // If found, update the clone features
    if ( found )
    {
      // We get the pointer to the clone and erase it from the map, because
      // each clone has a single master only and can therefore be updated 
      // only once
      pClone = imm->second;
      m_AccessToClones.erase( imm );
      
      if ( !update_velocity_only )
      {
        pClone->setPosition( &recvbuf_DOUBLE[j+13] );
        pClone->setQuaternionRotation( recvbuf_DOUBLE[j+6],
		recvbuf_DOUBLE[j+7],
		recvbuf_DOUBLE[j+8],
		recvbuf_DOUBLE[j+9] );
        if ( contact_map_size )
	{
          std::tuple<int,int,int> key;
          Vector3 kdelta, prev_normal, cumulSpringTorque;
	  int ii = 0;
          for (int m=0; m < contact_map_size; m++)
          {
            ii = j + NB_DOUBLE_PART + 1 + NB_DOUBLE_PER_CONTACT * m;
	    key = std::make_tuple(int(recvbuf_DOUBLE[ii]),
                  int(recvbuf_DOUBLE[ii+1]),
                  int(recvbuf_DOUBLE[ii+2]));
	    kdelta = Vector3( recvbuf_DOUBLE[ii+4],
                  recvbuf_DOUBLE[ii+5],
                  recvbuf_DOUBLE[ii+6]) ;
            prev_normal = Vector3( recvbuf_DOUBLE[ii+7],
                  recvbuf_DOUBLE[ii+8],
                  recvbuf_DOUBLE[ii+9]) ;
            cumulSpringTorque = Vector3( recvbuf_DOUBLE[ii+10],
                  recvbuf_DOUBLE[ii+11],
                  recvbuf_DOUBLE[ii+12]) ;
            pClone->copyContactInMap( key, bool(recvbuf_DOUBLE[ii+3]), kdelta,
                  prev_normal, cumulSpringTorque) ;
          }	
	}
      }      
      pClone->setTranslationalVelocity( recvbuf_DOUBLE[j+3],
	  	recvbuf_DOUBLE[j+4],
		recvbuf_DOUBLE[j+5] );
      pClone->setAngularVelocity(recvbuf_DOUBLE[j+10],
	  	recvbuf_DOUBLE[j+11],
		recvbuf_DOUBLE[j+12] );
    }
    else
    {
      cout << "!!! Warning !!! Clone " << id << " not found in proc " 
      	<< m_rank << endl;
    } 
    
    j += NB_DOUBLE_PART + 1 + contact_map_size * NB_DOUBLE_PER_CONTACT; 
  } 
}




// ----------------------------------------------------------------------------
// Creates and updates clones with the data sent by the neighboring processes 
void GrainsMPIWrapper::CreateClones(double time,
 	int const &recvsize, double const* recvbuf_DOUBLE,
	const int& NB_DOUBLE_PART,
	int const& NB_DOUBLE_PER_CONTACT,  
  	list<Particle*>* particlesClones,
	list<Particle*>* particles,
  	list<Particle*> const* particlesBufferzone,
	vector<Particle*> const* referenceParticles,
	LinkedCell* LC )
{
  int j, id, classe, contact_map_size = 0;
  bool found = false;
  double distGC = 0. ; 
  Point3 const* GC = NULL ; 
  multimap<int,Particle*>::iterator imm;
  size_t ncid = 0;
  pair < multimap<int,Particle*>::iterator, 
  	multimap<int,Particle*>::iterator > crange;
	  
  // Note: although we create new clones here, we need to check if they
  // do not already exist. This scenario happens in the 2 following cases:
  // 1) if a master particle moves from a buffer zone cell to another buffer 
  // zone cell and the two buffer zone cells have a different GeoPosition, 
  // then such a master particle is added to the list of new buffer zone 
  // particles but the corresponding clone particle might already exist on the 
  // local process
  // 2) in case of a reloaded simulation, if the linked cell has changed
  // from the previous simulation (the linked cell size has become larger), 
  // some periodic clones may not exist and we need te create them 

  for( j=0; j<recvsize; )
  {
    found = false;
    id = int( recvbuf_DOUBLE[j] ); 
    contact_map_size = int( recvbuf_DOUBLE[j+NB_DOUBLE_PART] );       
      
    // Search whether the clone already exists in this local domain
    ncid = m_AccessToClones.count( id );
    switch( ncid )
    {
      case 0: // no clone with this ID number in this local domain
        break;
      case 1: // 1 clone with this ID number in this local domain
        imm = m_AccessToClones.find( id );
        if ( m_isMPIperiodic )
        {
          GC = imm->second->getPosition();
          distGC = sqrt( 
          	pow( recvbuf_DOUBLE[j+22] - (*GC)[X], 2. ) +
            pow( recvbuf_DOUBLE[j+23] - (*GC)[Y], 2. ) +
            pow( recvbuf_DOUBLE[j+24] - (*GC)[Z], 2. ) ) ;
            if ( distGC < 1.1 * imm->second->getCrustThickness() )  
              found = true;	
        }
        else found = true;
        break;
      default: // more than 1 clone with this ID number in this local domain
        // Periodic case with 1 multi-proc clone and 1 periodic clone in the
	// same local domain
        crange = m_AccessToClones.equal_range( id );
        for (imm=crange.first; imm!=crange.second && !found; )
        {
          GC = imm->second->getPosition();
          distGC = sqrt( 
          pow( recvbuf_DOUBLE[j+22] - (*GC)[X], 2. ) +
          pow( recvbuf_DOUBLE[j+23] - (*GC)[Y], 2. ) +
          pow( recvbuf_DOUBLE[j+24] - (*GC)[Z], 2. ) ) ;
          if ( distGC < 1.1 * imm->second->getCrustThickness() )  
            found = true;
          else imm++;
        }
        break;    
    }   
       
    // If not found, create the clone particle
    if ( !found )
    {
      classe = int( recvbuf_DOUBLE[j+1] );
	
      if ( GrainsExec::m_MPI_verbose )
      {
        ostringstream oss;
        oss << "   t=" << GrainsExec::doubleToString(time,TIMEFORMAT)
		<< " Create Clone                Id = " 
		<< id
		<< " Type = " << classe << " " 
		<< recvbuf_DOUBLE[j+22] << " " 
		<< recvbuf_DOUBLE[j+23] << " " 
		<< recvbuf_DOUBLE[j+24]
		<< endl;
        GrainsMPIWrapper::addToMPIString(oss.str());
      }

      // Creation of the clone particle
      Particle *new_clone = NULL ;
      if ( (*referenceParticles)[classe]->isCompositeParticle() )
        new_clone = new CompositeParticle( id,
              (*referenceParticles)[classe],
              recvbuf_DOUBLE[j+3],
              recvbuf_DOUBLE[j+4],
              recvbuf_DOUBLE[j+5],
              recvbuf_DOUBLE[j+6],
              recvbuf_DOUBLE[j+7],
              recvbuf_DOUBLE[j+8],	
              recvbuf_DOUBLE[j+9],
              recvbuf_DOUBLE[j+10],
              recvbuf_DOUBLE[j+11],
              recvbuf_DOUBLE[j+12],
              &recvbuf_DOUBLE[j+13],
              COMPUTE, 2, 0 );   
      else
        new_clone = new Particle( id,
              (*referenceParticles)[classe],
              recvbuf_DOUBLE[j+3],
              recvbuf_DOUBLE[j+4],
              recvbuf_DOUBLE[j+5],
              recvbuf_DOUBLE[j+6],
              recvbuf_DOUBLE[j+7],
              recvbuf_DOUBLE[j+8],	
              recvbuf_DOUBLE[j+9],
              recvbuf_DOUBLE[j+10],
              recvbuf_DOUBLE[j+11],
              recvbuf_DOUBLE[j+12],
              &recvbuf_DOUBLE[j+13],
              COMPUTE, 2, 0 );

      if ( contact_map_size )
      {
        std::tuple<int,int,int> key;
        Vector3 kdelta, prev_normal, cumulSpringTorque;
	int ii = 0;
        for (int m=0; m < contact_map_size; m++)
        {
          ii = j + NB_DOUBLE_PART + 1 + NB_DOUBLE_PER_CONTACT * m;
	  key = std::make_tuple(int(recvbuf_DOUBLE[ii]),
                  int(recvbuf_DOUBLE[ii+1]),
                  int(recvbuf_DOUBLE[ii+2]));
	  kdelta = Vector3( recvbuf_DOUBLE[ii+4],
                  recvbuf_DOUBLE[ii+5],
                  recvbuf_DOUBLE[ii+6]) ;
          prev_normal = Vector3( recvbuf_DOUBLE[ii+7],
                  recvbuf_DOUBLE[ii+8],
                  recvbuf_DOUBLE[ii+9]) ;
          cumulSpringTorque = Vector3( recvbuf_DOUBLE[ii+10],
                  recvbuf_DOUBLE[ii+11],
                  recvbuf_DOUBLE[ii+12]) ;
          new_clone->copyContactInMap( key, bool(recvbuf_DOUBLE[ii+3]), kdelta,
                  prev_normal, cumulSpringTorque) ;
        }	
      }

      // Add to the LinkedCell
      LC->Link( new_clone );
	
      // Add to active particle and clone particle lists
      particlesClones->push_back( new_clone );
      particles->push_back( new_clone ); 
    }
    else m_AccessToClones.erase( imm ); 
    
    j += NB_DOUBLE_PART + 1 + contact_map_size * NB_DOUBLE_PER_CONTACT; 
  } 
}




// ----------------------------------------------------------------------------
// Broadcasts a double from the master to all processes within the 
// MPI_COMM_activProc communicator
double GrainsMPIWrapper::Broadcast_DOUBLE( double const& d, int source ) const
{
  double collective_d = d;
  
  MPI_Bcast( &collective_d, 1, MPI_DOUBLE, source, m_MPI_COMM_activeProc );
  
  return ( collective_d ); 
}




// ----------------------------------------------------------------------------
// Broadcasts an integer from the master to all processes within the 
// MPI_COMM_activProc communicator
int GrainsMPIWrapper::Broadcast_INT( int const& i, int source ) const
{
  int collective_i = i;
  
  MPI_Bcast( &collective_i, 1, MPI_INT, source, m_MPI_COMM_activeProc );
  
  return ( collective_i ); 
}




// ----------------------------------------------------------------------------
// Broadcasts an integer from the master to all processes within the 
// MPI_COMM_activProc communicator
size_t GrainsMPIWrapper::Broadcast_UNSIGNED_INT( size_t const& i ) const
{
  size_t collective_i = i;
  
  MPI_Bcast( &collective_i, 1, MPI_UNSIGNED_LONG, 0, m_MPI_COMM_activeProc );
  
  return collective_i; 
}




// ----------------------------------------------------------------------------
// Sums a double from all processes on all processes within the 
// MPI_COMM_activProc communicator
double GrainsMPIWrapper::sum_DOUBLE( double const& x ) const
{
  double sum = 0;
  
  MPI_Allreduce( &x, &sum, 1, MPI_DOUBLE, MPI_SUM, m_MPI_COMM_activeProc );
  
  return ( sum );
}




// ----------------------------------------------------------------------------
// Sums a double from all processes on the master process within the 
// MPI_COMM_activProc communicator
double GrainsMPIWrapper::sum_DOUBLE_master( double const& x ) const
{
  double sum = 0;
  
  MPI_Reduce( &x, &sum, 1, MPI_DOUBLE, MPI_SUM, m_rank_master,
  	m_MPI_COMM_activeProc );
  
  return ( sum );
}




// ----------------------------------------------------------------------------
// Sums an integer from all processes on all processes within the 
// MPI_COMM_activProc communicator
int GrainsMPIWrapper::sum_INT( int const& i ) const
{
  int sum = 0;
  
  MPI_Allreduce( &i, &sum, 1, MPI_INT, MPI_SUM, m_MPI_COMM_activeProc );
  
  return ( sum );
}




// ----------------------------------------------------------------------------
// Sums an integer from all processes on all processes within the 
// MPI_COMM_activProc communicator
size_t GrainsMPIWrapper::sum_UNSIGNED_INT( size_t const& i ) const
{
  size_t sum = 0;
  
  MPI_Allreduce( &i, &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 
  	m_MPI_COMM_activeProc );
  
  return ( sum );
}




// ----------------------------------------------------------------------------
// Sums an integer from all processes on the master process within the 
// MPI_COMM_activProc communicator
int GrainsMPIWrapper::sum_INT_master( int const& i ) const
{
  int sum = 0;
  
  MPI_Reduce( &i, &sum, 1, MPI_INT, MPI_SUM, m_rank_master, 
  	m_MPI_COMM_activeProc );
  
  return ( sum );
}




// ----------------------------------------------------------------------------
// Sums an unsigned integer from all processes on the master process 
// within the MPI_COMM_activProc communicator
size_t GrainsMPIWrapper::sum_UNSIGNED_INT_master( size_t const& i ) const
{
  int sum = 0;
  
  MPI_Reduce( &i, &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, m_rank_master, 
  	m_MPI_COMM_activeProc );
  
  return ( sum );
}




// ----------------------------------------------------------------------------
// Performs a "logical and" operation on input boolean value from 
// all processes on all processes
bool GrainsMPIWrapper::logical_and( bool const& input ) const
{
  int land=0;
  
  MPI_Allreduce( &input, &land, 1, MPI_UNSIGNED_SHORT, MPI_LAND,
  	m_MPI_COMM_activeProc );
  
  return ( land );
}




// ----------------------------------------------------------------------------
// Minimum of an integer from all processes on all processes within the 
// MPI_COMM_activProc communicator
int GrainsMPIWrapper::min_INT( int const& i ) const
{
  int collective_i = 0;

  MPI_Allreduce( &i, &collective_i, 1, MPI_INT, MPI_MIN, 
  	m_MPI_COMM_activeProc ) ;

  return ( collective_i );
}  




// ----------------------------------------------------------------------------
// Max d'un entier sur tous les proc
int GrainsMPIWrapper::max_INT( int const& i ) const
{
  int collective_i = 0;

  MPI_Allreduce( &i, &collective_i, 1, MPI_INT, MPI_MAX, 
  	m_MPI_COMM_activeProc ) ;

  return ( collective_i );
}  




// ----------------------------------------------------------------------------
// Min d'un entier non signe sur tous les proc
size_t GrainsMPIWrapper::min_UNSIGNED_INT( size_t const& i ) const
{
  size_t collective_i = 0;

  MPI_Allreduce( &i, &collective_i, 1, MPI_UNSIGNED_LONG, MPI_MIN, 
  	m_MPI_COMM_activeProc ) ;

  return ( collective_i );
}  




// ----------------------------------------------------------------------------
// Maximum of an unsigned integer from all processes on all processes within 
// the MPI_COMM_activProc communicator
size_t GrainsMPIWrapper::max_UNSIGNED_INT( size_t const& i ) const
{
  size_t collective_i = 0;

  MPI_Allreduce( &i, &collective_i, 1, MPI_UNSIGNED_LONG, MPI_MAX, 
  	m_MPI_COMM_activeProc ) ;

  return ( collective_i );
}  




// ----------------------------------------------------------------------------
// Maximum of a double from all processes on all processes within the 
// MPI_COMM_activProc communicator
double GrainsMPIWrapper::max_DOUBLE( double const& x ) const 
{
  double max = 0.;

  MPI_Allreduce( &x, &max, 1, MPI_DOUBLE, MPI_MAX, m_MPI_COMM_activeProc ) ;

  return ( max );
} 



 
// ----------------------------------------------------------------------------
// Maximum of a double from all processes on the master process within the 
// MPI_COMM_activProc communicator
double GrainsMPIWrapper::max_DOUBLE_master( double const& x ) const 
{
  double max = 0.;

  MPI_Reduce( &x, &max, 1, MPI_DOUBLE, MPI_MAX, m_rank_master,
  	m_MPI_COMM_activeProc ) ;

  return ( max );
} 




// ----------------------------------------------------------------------------
// Minimum of a double from all processes on all processes within the 
// MPI_COMM_activProc communicator
double GrainsMPIWrapper::min_DOUBLE( double const& x ) const 
{
  double min = 0.;

  MPI_Allreduce( &x, &min, 1, MPI_DOUBLE, MPI_MIN, m_MPI_COMM_activeProc ) ;

  return ( min );
} 




// ----------------------------------------------------------------------------
// Minimum of a double from all processes on the master process within the 
// MPI_COMM_activProc communicator
double GrainsMPIWrapper::min_DOUBLE_master( double const& x ) const 
{
  double min = 0.;

  MPI_Reduce( &x, &min, 1, MPI_DOUBLE, MPI_MIN, m_rank_master,
  	m_MPI_COMM_activeProc ) ;

  return ( min );
} 




// ----------------------------------------------------------------------------
// AllGather of an unsigned integer from all processes on all processes within 
// the MPI_COMM_activProc communicator
size_t* GrainsMPIWrapper::AllGather_UNSIGNED_INT( size_t const& i ) const
{
  size_t* recv = new size_t[m_nprocs];
  
  MPI_Allgather( &i, 1, MPI_UNSIGNED_LONG, recv, 1, MPI_UNSIGNED_LONG, 
  	m_MPI_COMM_activeProc ); 

  return ( recv );
} 




// ----------------------------------------------------------------------------
// AllGather of an unsigned integer from all processes on all processes within 
// the MPI_COMM_activProc communicator
int* GrainsMPIWrapper::AllGather_INT( int const& i ) const
{
  int* recv = new int[m_nprocs];
  
  MPI_Allgather( &i, 1, MPI_INT, recv, 1, MPI_INT, m_MPI_COMM_activeProc ); 

  return ( recv );
} 




// ----------------------------------------------------------------------------
// Gather of an integer from all processes on the master process 
// within the MPI_COMM_activProc communicator
int* GrainsMPIWrapper::Gather_INT_master( int const& i ) const
{
  int* recv = NULL;
  if ( m_rank == 0 ) recv = new int[m_nprocs];
  
  MPI_Gather( &i, 1, MPI_INT, recv, 1, MPI_INT, 0, m_MPI_COMM_activeProc ); 

  return ( recv );
} 




// ----------------------------------------------------------------------------
// Gather of an unsigned integer from all processes on the master process 
// within the MPI_COMM_activProc communicator
size_t* GrainsMPIWrapper::Gather_UNSIGNED_INT_master( size_t const& i ) const
{
  size_t* recv = NULL;
  if ( m_rank == 0 ) recv = new size_t[m_nprocs];
  
  MPI_Gather( &i, 1, MPI_UNSIGNED_LONG, recv, 1, MPI_UNSIGNED_LONG, 0, 
  	m_MPI_COMM_activeProc ); 

  return ( recv );
} 




// ----------------------------------------------------------------------------
// Broadcasts a 3D point from the master to all processes within 
// the MPI_COMM_activProc communicator
Point3 GrainsMPIWrapper::Broadcast_Point3( Point3 const& pt ) const
{
  double *coordinates = new double[3]; 
  coordinates[0] = pt[X];
  coordinates[1] = pt[Y];  
  coordinates[2] = pt[Z];  

  MPI_Bcast( coordinates, 3, MPI_DOUBLE, 0, m_MPI_COMM_activeProc );  
  
  Point3 cpt( coordinates[0], coordinates[1], coordinates[2] );
  delete [] coordinates;
  
  return ( cpt );
}




// ----------------------------------------------------------------------------
// Broadcasts a 3D vector from the master to all processes within the 
// MPI_COMM_activProc communicator
Vector3 GrainsMPIWrapper::Broadcast_Vector3( Vector3 const& v ) const
{
  double *coordinates = new double[3]; 
  coordinates[0] = v[X];
  coordinates[1] = v[Y];  
  coordinates[2] = v[Z];  

  MPI_Bcast( coordinates, 3, MPI_DOUBLE, 0, m_MPI_COMM_activeProc );  
  
  Vector3 cv( coordinates[0], coordinates[1], coordinates[2] );
  delete [] coordinates;
  
  return ( cv );
}




// ----------------------------------------------------------------------------
// Broadcasts a 3D vector from the master to all processes within the 
// MPI_COMM_activProc communicator
Matrix GrainsMPIWrapper::Broadcast_Matrix( Matrix const& mat ) const
{
  double *mat_coef = new double[9];
  Mat3 const& mmat = mat.getValue(); 
  mat_coef[0] = mmat[X][X];
  mat_coef[1] = mmat[X][Y];  
  mat_coef[2] = mmat[X][Z];  
  mat_coef[3] = mmat[Y][X];
  mat_coef[4] = mmat[Y][Y];  
  mat_coef[5] = mmat[Y][Z];   
  mat_coef[6] = mmat[Z][X];
  mat_coef[7] = mmat[Z][Y];  
  mat_coef[8] = mmat[Z][Z];     

  MPI_Bcast( mat_coef, 9, MPI_DOUBLE, 0, m_MPI_COMM_activeProc );  
  
  Matrix bmat( mat_coef[0], mat_coef[1], mat_coef[2],
  	mat_coef[3], mat_coef[4], mat_coef[5],
	mat_coef[6], mat_coef[7], mat_coef[8]);
  delete [] mat_coef;
  
  return ( bmat );
}




// ----------------------------------------------------------------------------
// Outputs timer summary
void GrainsMPIWrapper::timerSummary( ostream &f ) const
{
  double cputime = SCT_get_total_elapsed_time( "BuffersCopy" )
     	+ SCT_get_total_elapsed_time( "MPIComm" )
     	+ SCT_get_total_elapsed_time( "UpdateCreateClones" );
  f << "MPI wrapper timer" << endl;
  SCT_get_summary( f, cputime );
}




// ----------------------------------------------------------------------------
// Adds a string to the MPI log string
void GrainsMPIWrapper::addToMPIString( string const& add )
{
  if ( m_MPILogString ) *m_MPILogString += add;
} 




// ----------------------------------------------------------------------------
// Outputs the MPI log string per process and reinitialize it to empty
void GrainsMPIWrapper::writeAndFlushMPIString( ostream &f )
{  
  writeStringPerProcess( f, *m_MPILogString, true );

  delete m_MPILogString;
  m_MPILogString = new string;
}




// ----------------------------------------------------------------------------
// MPI_Barrier pour les processus actifs uniquement
void GrainsMPIWrapper::MPI_Barrier_ActivProc() const
{
  MPI_Barrier( m_MPI_COMM_activeProc );

}




// ----------------------------------------------------------------------------
// Shares contact features among all active processes
void GrainsMPIWrapper::ContactsFeatures( double& overlap_max,
	double& overlap_mean,
	double& time_overlapMax,
	double& nbIterGJK_mean ) const
{
  double *recvbuf_overlap = new double[m_nprocs];
  double *recvbuf_time = new double[m_nprocs];   

  MPI_Allgather( &overlap_max, 1, MPI_DOUBLE, recvbuf_overlap,
  	1, MPI_DOUBLE, m_MPI_COMM_activeProc ); 
  MPI_Allgather( &time_overlapMax, 1, MPI_DOUBLE, recvbuf_time,
  	1, MPI_DOUBLE, m_MPI_COMM_activeProc );	 
	
  double ovmax=0.;
  for (int i=0;i<m_nprocs;++i)
    if ( recvbuf_overlap[i] > ovmax )
    {
      overlap_max = recvbuf_overlap[i];
      time_overlapMax = recvbuf_time[i];
      ovmax = overlap_max;
    }

  overlap_mean = sum_DOUBLE( overlap_mean );
  overlap_mean /= m_nprocs;
  
  nbIterGJK_mean = sum_DOUBLE( nbIterGJK_mean );
  nbIterGJK_mean /= m_nprocs;

  delete [] recvbuf_overlap;
  delete [] recvbuf_time;   
}




// ----------------------------------------------------------------------------
// Sums force & torque exerted on obstacles on the master process
void GrainsMPIWrapper::sumObstaclesLoad( 
	list<SimpleObstacle*> const& allMyObs ) const
{
  list<SimpleObstacle*>::const_iterator obstacle ;
  int nobs = int(allMyObs.size()), i = 0 ;
  Vector3 const* force = NULL; 
  Vector3 const* torque = NULL; 
  Vector3 collective_force, collective_torque; 
  double* forcetorque = new double[6*nobs];
  double* forcetorque_collective = new double[6*nobs];   

  // Copy in a local buffer except on the master process where buffer is set 
  // to 0
  if ( m_rank == m_rank_master )
    for (i=0;i<6*nobs;++i) forcetorque[i] = 0.;
  else
    for (obstacle=allMyObs.begin(); obstacle!=allMyObs.end(); obstacle++,i+=6)
    {
      force = (*obstacle)->getForce();
      torque = (*obstacle)->getTorque();
      forcetorque[i] = (*force)[X];
      forcetorque[i+1] = (*force)[Y];    
      forcetorque[i+2] = (*force)[Z];    
      forcetorque[i+3] = (*torque)[X];    
      forcetorque[i+4] = (*torque)[Y];   
      forcetorque[i+5] = (*torque)[Z];        
    } 
  
  // Sum all contributions from other processes on the master
  MPI_Reduce( forcetorque, forcetorque_collective, 6*nobs, MPI_DOUBLE, 
  	MPI_SUM, m_rank_master, m_MPI_COMM_activeProc ); 
	
  // Add all contributions from other processes on the master
  if ( m_rank == m_rank_master )
  {
    i = 0;
    for (obstacle=allMyObs.begin(); obstacle!=allMyObs.end(); obstacle++,i+=6)
    {
      collective_force[X] = forcetorque_collective[i] ;
      collective_force[Y] = forcetorque_collective[i+1] ;    
      collective_force[Z] = forcetorque_collective[i+2] ;    
      collective_torque[X] = forcetorque_collective[i+3] ;
      collective_torque[Y] = forcetorque_collective[i+4] ;    
      collective_torque[Z] = forcetorque_collective[i+5] ; 
      (*obstacle)->addBodyForce( collective_force );
      (*obstacle)->addTorque( collective_torque );           
    }
  }  	     
  
  delete [] forcetorque;
  delete [] forcetorque_collective;  
}




// ----------------------------------------------------------------------------
// Writes the memory consumption per process in a stream
void GrainsMPIWrapper::display_used_memory( ostream& f ) const
{
  if ( m_rank == m_rank_master ) f << "Memory used by Grains3D" << endl;
  
  ostringstream out;
  GrainsExec::display_memory( out, GrainsExec::used_memory() ); 
  writeStringPerProcess( f, out.str(), false, GrainsExec::m_shift3 );   
}




// ----------------------------------------------------------------------------
// Distributes the number of particles in each class and on each
// process in the case of the block structured insertion
void GrainsMPIWrapper::distributeParticlesClassProc( 
  	list< pair<Particle*,size_t> > const& newPart,
	list< pair<Particle*,size_t> >&newPartProc,
	size_t const& npartproc,
	size_t const& ntotalinsert ) const
{
  list< pair<Particle*,size_t> >::const_iterator ipart;
  list< pair<Particle*,size_t> >::iterator ipartProc;
  size_t nclasseproc = 0, npartproc_ = npartproc, i, 
  	nbClasses = newPart.size() ;
  MPI_Status status;
  

  // Initialisation of lists and arrays
  newPartProc = newPart;
  size_t* tabNbPartRestantesParClasse = new size_t[nbClasses];
  for (i=0,ipart=newPart.begin();i<nbClasses;++i,ipart++) 
    tabNbPartRestantesParClasse[i] = ipart->second;
  
  
  // Distribution is performed sequentially from 0 to the last
  // process
  if ( m_rank == m_rank_master )
  {
    // All types except the last type
    ipart = newPart.begin();
    ipartProc = newPartProc.begin();  
    for (i=0,ipart=newPart.begin(),ipartProc=newPartProc.begin();
    	i<nbClasses-1;++i,ipart++,ipartProc++)
    {
      nclasseproc = size_t( npartproc * ipart->second / ntotalinsert ) ;
      // If ratio is integer, we keep it, otherwise we add 1
      if ( fabs( double( npartproc * ipart->second / ntotalinsert )
      	- double( size_t( npartproc * ipart->second / ntotalinsert ) ) ) 
	> 1.e-14 )
	++nclasseproc;
      if ( nclasseproc > tabNbPartRestantesParClasse[i] ) 
        nclasseproc = tabNbPartRestantesParClasse[i];
      if ( nclasseproc > npartproc_ ) nclasseproc = npartproc_; 
      ipartProc->second = nclasseproc;
      npartproc_ -= nclasseproc;
      tabNbPartRestantesParClasse[i] -= nclasseproc;
    }
  
    // Correction for the last type to satisfy that the sum of the numbers of
    // particles per type is equal to the total number of particles to insert on
    // this process
    nclasseproc = npartproc_;
    ipartProc->second = nclasseproc;
    tabNbPartRestantesParClasse[nbClasses-1] -= nclasseproc;
  
    MPI_Send( tabNbPartRestantesParClasse, int(nbClasses), MPI_UNSIGNED_LONG, 
    	1, m_rank, m_MPI_COMM_activeProc );
  }
  else
  {
    if ( m_rank != m_nprocs - 1 )
    {
      MPI_Recv( tabNbPartRestantesParClasse, int(nbClasses), MPI_UNSIGNED_LONG, 
      	m_rank - 1, m_rank - 1, m_MPI_COMM_activeProc, &status );

      // All types except the last type
      for (i=0,ipart=newPart.begin(),ipartProc=newPartProc.begin();
    	i<nbClasses-1;++i,ipart++,ipartProc++)      
      {
        nclasseproc = size_t( npartproc * ipart->second / ntotalinsert ) ;
        // If ratio is integer, we keep it, otherwise we add 1
        if ( fabs( double( npartproc * ipart->second / ntotalinsert )
      	- double( size_t( npartproc * ipart->second / ntotalinsert ) ) ) 
	> 1.e-14 )
	  ++nclasseproc;        
	if ( nclasseproc > tabNbPartRestantesParClasse[i] ) 
          nclasseproc = tabNbPartRestantesParClasse[i];
        if ( nclasseproc > npartproc_ ) nclasseproc = npartproc_; 
        ipartProc->second = nclasseproc;
        npartproc_ -= nclasseproc;
        tabNbPartRestantesParClasse[i] -= nclasseproc;
      }
  
      // Correction for the last class to satisfy that the sum of the numbers of
      // particles per class is equal to the total number of particles to insert
      // on this process
      nclasseproc = npartproc_;
      ipartProc->second = nclasseproc;
      tabNbPartRestantesParClasse[nbClasses-1] -= nclasseproc;
  
      MPI_Send( tabNbPartRestantesParClasse, int(nbClasses), MPI_UNSIGNED_LONG, 
      	m_rank + 1, m_rank, m_MPI_COMM_activeProc );
    }
    else
    {
      MPI_Recv( tabNbPartRestantesParClasse, int(nbClasses), MPI_UNSIGNED_LONG, 
      	m_rank - 1, m_rank - 1, m_MPI_COMM_activeProc, &status );      
	
      for (i=0,ipartProc=newPartProc.begin();i<nbClasses;++i,ipartProc++)
        ipartProc->second = tabNbPartRestantesParClasse[i];
    }
  }    	
  
  // Output and verify distribution per class and per process
  if ( m_rank == m_rank_master )
    cout << endl << "Distribution " << endl << 
    	"Per process and per class" << endl;
  for (int m=0;m<m_nprocs;++m)
  {
    if ( m == m_rank && m_is_active )
    {      
      cout << "Proc " << m_rank << endl;
      size_t ntot_ = 0 ;
      for (i=0,ipart=newPart.begin(),ipartProc=newPartProc.begin();
    	i<nbClasses;++i,ipart++,ipartProc++)
      {
        ntot_ += ipartProc->second;
	cout << "   Class " << i << " = " << ipartProc->second << endl;
      }
      cout << "   Total: " << ntot_ << " = " << npartproc << endl;

    }
    MPI_Barrier( m_MPI_COMM_activeProc );
  }
  MPI_Barrier( m_MPI_COMM_activeProc );

  if ( m_rank == m_rank_master )
    cout << "Per class on all processes" << endl;
  for (i=0,ipart=newPart.begin(),ipartProc=newPartProc.begin();
    	i<nbClasses;++i,ipart++,ipartProc++)
  {
    size_t ntotc =  sum_UNSIGNED_INT_master( ipartProc->second ) ;
    if ( m_rank == m_rank_master )
      cout << "   Class " << i << " " << ntotc << " = " << ipart->second 
      	<< endl;
  }
  if ( m_rank == m_rank_master ) cout << endl;
} 




// ----------------------------------------------------------------------------
// Writes a string per process in a process-id ordered manner
void GrainsMPIWrapper::writeStringPerProcess( ostream& f, string const& out, 
    	bool creturn, string const& shift ) const
{
  string allout = shift + "Process 0";
  bool write = false;
  int dim = 0, recvsize = 0 ;
  MPI_Status status;
  
  if ( m_rank == m_rank_master ) 
  {
    if ( creturn ) allout += "\n";
    else allout += ": ";
    allout += out + "\n";
    if ( !out.empty() ) write = true;   
    for ( int irank=1;irank<m_nprocs;irank++)
    {
      // Size of the message
      MPI_Probe( irank, m_tag_CHAR, m_MPI_COMM_activeProc, &status );  
      MPI_Get_count( &status, MPI_CHAR, &recvsize );

      // Reception of the actual message	
      char *recvbuf_char = new char[recvsize];
      MPI_Recv( recvbuf_char, recvsize, MPI_CHAR, 
          irank, m_tag_CHAR, m_MPI_COMM_activeProc, &status );

      if ( recvsize != 1 )
      { 
        write = true; 
	string recstring = string( recvbuf_char );
        allout += shift + "Process " + GrainsExec::intToString( irank );
        if ( creturn ) allout += "\n";
        else allout += ": ";
        allout += recstring;
        if ( irank != m_nprocs - 1 ) allout += "\n";
      }
      delete recvbuf_char;       	
    }
  }
  else
  {
    // Send the string to the master proc
    dim = int(out.size()+1);
    MPI_Send( out.c_str(), dim, MPI_CHAR, 0, m_tag_CHAR, 
    	m_MPI_COMM_activeProc );
  }

  // The master proc writes to the stream
  if ( m_rank == m_rank_master && write ) f << allout << endl; 
}




// ----------------------------------------------------------------------------
// Sends an array of integers
void GrainsMPIWrapper::send( int const* tab, int const& dim, int const& to ) 
	const
{
  MPI_Send( tab, dim, MPI_INT, to, m_tag_INT, m_MPI_COMM_activeProc );
}




// ----------------------------------------------------------------------------
// Receives an array of integers
void GrainsMPIWrapper::receive( int* &recvbuf, int &dim, int const& from ) 
	const
{
  MPI_Status status;

  // Size of the message
  MPI_Probe( from, m_tag_INT, m_MPI_COMM_activeProc, &status );  
  MPI_Get_count( &status, MPI_INT, &dim );

  // Reception of the actual message	
  recvbuf = new int[dim];
  MPI_Recv( recvbuf, dim, MPI_INT, from, m_tag_INT, 
  	m_MPI_COMM_activeProc, &status );
}




// ----------------------------------------------------------------------------
// Sends an array of doubles
void GrainsMPIWrapper::send( double const* tab, int const& dim, int const& to ) 
	const
{
  MPI_Send( tab, dim, MPI_DOUBLE, to, m_tag_DOUBLE, m_MPI_COMM_activeProc );
}




// ----------------------------------------------------------------------------
// Receives an array of doubles
void GrainsMPIWrapper::receive( double* &recvbuf, int &dim, int const& from ) 
	const
{
  MPI_Status status;

  // Size of the message
  MPI_Probe( from, m_tag_DOUBLE, m_MPI_COMM_activeProc, &status );  
  MPI_Get_count( &status, MPI_DOUBLE, &dim );

  // Reception of the actual message	
  recvbuf = new double[dim];
  MPI_Recv( recvbuf, dim, MPI_DOUBLE, from, m_tag_DOUBLE, 
  	m_MPI_COMM_activeProc, &status );
}




// ----------------------------------------------------------------------------
// Returns the active process communicator
MPI_Comm GrainsMPIWrapper::get_active_procs_comm() const
{
  return ( m_MPI_COMM_activeProc ) ;
} 
