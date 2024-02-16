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
vector< vector<int> > GrainsMPIWrapper::m_particleHalozoneToNeighboringProcs;
vector<int> GrainsMPIWrapper::m_GeoLocReciprocity;


// ----------------------------------------------------------------------------
// Constructor with domain decomposition and periodicity as input parameters
GrainsMPIWrapper::GrainsMPIWrapper( int NX, int NY, int NZ,
  	int PERX, int PERY, int PERZ, string const& oshift ) : 
   m_coords( NULL ),
   m_dim( NULL ),
   m_period( NULL ),
   m_isMPIperiodic( false ),   
   m_rank( 0 ), 
   m_rank_world( 0 ),
   m_rank_masterWorld( 0 ),
   m_nprocs( 0 ), 
   m_nprocs_world( 0 ),
   m_is_activ( false ),
   m_neighbors( NULL ),
   m_commgrainsMPI_3D( NULL ),
   m_nprocs_localComm( 0 ),
   m_rank_localComm( 0 ), 
   m_master_localComm( NULL )
{ 
  // Total number of processes and data in the world communicator
  MPI_Comm_size( MPI_COMM_WORLD, &m_nprocs_world );
  MPI_Comm_rank( MPI_COMM_WORLD, &m_rank_world );  
  if ( NX * NY * NZ > m_nprocs_world )
  {
    if ( m_rank_world == 0 )
      cout << endl << "!!! WARNING !!!" << endl << 
      	"Domain decomposition does not match total number of processus" 
	<< endl << endl;
    int error_code = 0;
    MPI_Abort( MPI_COMM_WORLD, error_code );
  }
  else
  {
    m_nprocs = NX * NY * NZ;
    if ( m_rank_world == 0 )
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
    	&m_MPI_COMM_activProc );
  else
    MPI_Comm_dup( MPI_COMM_WORLD, &m_MPI_COMM_activProc );  
  if ( m_rank_world < m_nprocs ) m_is_activ = true;   
  else m_is_activ = false;  
  MPI_Group_free( &world_group );

  // Among the active processes
  if ( m_is_activ )
  {     
    // Rank of the process in the group of active processes
    MPI_Comm_rank( m_MPI_COMM_activProc, &m_rank ); 

    // Rank in the m_MPI_COMM_activProc communicator of the process that has
    // rank 0 (master) in the MPI_COMM_WORLD communicator
    int loc_m_rank_masterWorld = -1;
    if ( m_rank_world == 0 ) loc_m_rank_masterWorld = m_rank;
    MPI_Allreduce( &loc_m_rank_masterWorld, &m_rank_masterWorld, 1, MPI_INT, 
    	MPI_MAX, m_MPI_COMM_activProc ) ;   

    // Number of processes per direction
    m_dim = new int[3];
    m_dim[0] = NX, m_dim[1] = NY, m_dim[2] = NZ;
  
    // Periodicity of the domain
    m_period = new int[3];
    m_period[0] = PERX, m_period[1] = PERY, m_period[2] = PERZ;  
    if ( m_period[0] || m_period[1] || m_period[2] ) m_isMPIperiodic = true;
    
    // Creates vectors of periodicity of the domain
    m_MPIperiodes.reserve( 27 );
    for (int i=0;i<27;++i) m_MPIperiodes.push_back( Vector3Nul ); 

    // Sets the relationship between the GeoPosition in a halo
    // zone from which data are sent and the GeoPosition of the 
    // neighboring processes that receive the data
    setParticleHalozoneToNeighboringProcs();

    // Re-number processes: no 
    int reorganisation = 0; 
      
    // Creates the MPI Cartesian topology with periodicity if any
    m_commgrainsMPI_3D = new MPI_Comm;
    MPI_Cart_create( m_MPI_COMM_activProc, 3, m_dim, m_period, reorganisation, 
   	m_commgrainsMPI_3D );
     
    // Rank and MPI cartesian coordinates
    // We assume that the process rank is the same in m_MPI_COMM_activProc
    // and in m_commgrainsMPI_3D. This may need to be verified with an assert.
    MPI_Comm_rank( *m_commgrainsMPI_3D, &m_rank );
    if ( m_nprocs == m_nprocs_world ) assert( m_rank == m_rank_world );
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
  else m_rank = -1;
    
}




// ----------------------------------------------------------------------------
// Destructor
GrainsMPIWrapper::~GrainsMPIWrapper()
{
  MPI_Group_free( &m_MPI_GROUP_activProc );
  if ( m_is_activ ) 
  {
    MPI_Comm_free( &m_MPI_COMM_activProc ); 
    delete [] m_coords;
    delete [] m_dim;
    delete [] m_period;
    MPI_Comm_free( m_commgrainsMPI_3D );
    delete m_commgrainsMPI_3D;
    delete m_neighbors;
    vector<MPI_Group*>::iterator ivg;
    for (ivg=m_groupMPINeighbors.begin();ivg!=m_groupMPINeighbors.end();ivg++) 
    {
      MPI_Group_free( *ivg );    
      delete *ivg;
    }
    m_groupMPINeighbors.clear();
    if (!m_isInCommMPINeighbors.empty())
      for (int j=0;j<m_nprocs;j++)
      {
        if ( m_isInCommMPINeighbors[j] ) MPI_Comm_free( m_commMPINeighbors[j] );
        delete m_commMPINeighbors[j]; 
      }
    m_commMPINeighbors.clear();
    m_isInCommMPINeighbors.clear();  
    if ( m_master_localComm ) delete [] m_master_localComm; 
    if ( m_MPILogString ) delete m_MPILogString;
    vector< vector<int> >::iterator ivv;
    for (ivv=m_particleHalozoneToNeighboringProcs.begin();
  	ivv!=m_particleHalozoneToNeighboringProcs.end();ivv++)
      ivv->clear();
    m_particleHalozoneToNeighboringProcs.clear();
    
    m_MPIperiodes.clear();  
  }
}




// ----------------------------------------------------------------------------
// Sets the local MPI communicators involving the MPI cartesian neighbors
void GrainsMPIWrapper::setCommLocal()
{
  int i,j;

  // Size of messages is proportional to the number of neighbors
  int nbv = m_neighbors->nbMPINeighbors();

  // Size of messages
  int* recvcounts = new int[m_nprocs];
  for (i=0; i<m_nprocs; ++i) recvcounts[i]=0;
  MPI_Allgather( &nbv, 1, MPI_INT, recvcounts, 1, MPI_INT, 
  	m_MPI_COMM_activProc );

  // Partition and size of the reception buffer
  int* displs = new int[m_nprocs];
  displs[0] = 0; 
  for (i=1; i<m_nprocs; ++i) displs[i] = displs[i-1]+recvcounts[i-1];  
  int recvsize = displs[m_nprocs-1] + recvcounts[m_nprocs-1];

    
  // Communicate MPI cartesian neighbor ranks
  int* rankMPINeighbors = m_neighbors->rankMPINeighbors();
  int* recvbuf_rang = new int[recvsize];  
  MPI_Allgatherv( rankMPINeighbors, nbv, MPI_INT, recvbuf_rang, 
  	recvcounts, displs, MPI_INT, m_MPI_COMM_activProc );


  // Create local MPI communicators involving the MPI cartesian neighbors
  m_groupMPINeighbors.reserve( m_nprocs );  
  MPI_Group* empty_MPI_Group = NULL;
  for (j=0;j<m_nprocs;j++) m_groupMPINeighbors.push_back(empty_MPI_Group);
  m_commMPINeighbors.reserve( m_nprocs );  
  MPI_Comm* empty_MPI_Comm = NULL;
  for (j=0;j<m_nprocs;j++) m_commMPINeighbors.push_back(empty_MPI_Comm);
  m_isInCommMPINeighbors.reserve( m_nprocs ); 
  for (j=0;j<m_nprocs;j++) m_isInCommMPINeighbors.push_back(false); 
  MPI_Group activ_group;
  MPI_Comm_group( m_MPI_COMM_activProc, &activ_group );   
  for (j=0;j<m_nprocs;j++)
  {
    for (i=displs[j];i<displs[j]+recvcounts[j];++i)
      if ( recvbuf_rang[i] == m_rank ) m_isInCommMPINeighbors[j] = true;
    m_groupMPINeighbors[j] = new MPI_Group;    
    m_commMPINeighbors[j] = new MPI_Comm;
    MPI_Group_incl( activ_group, recvcounts[j], &recvbuf_rang[displs[j]], 
    	m_groupMPINeighbors[j] );
    MPI_Comm_create( m_MPI_COMM_activProc, *m_groupMPINeighbors[j], 
    	m_commMPINeighbors[j] );
  }   
  MPI_Group_free(&activ_group);
  

  // Rank in and size of the local MPI communicator
  MPI_Comm_rank( *m_commMPINeighbors[m_rank], &m_rank_localComm );
  MPI_Comm_size( *m_commMPINeighbors[m_rank], &m_nprocs_localComm );
  m_master_localComm = new int[m_nprocs];
  MPI_Allgather( &m_rank_localComm, 1, MPI_INT, m_master_localComm, 1, MPI_INT,
  	m_MPI_COMM_activProc ); 

		         
  delete [] recvcounts;
  delete [] displs; 
  delete [] rankMPINeighbors;
  delete [] recvbuf_rang;


  // Set the GeoPosition of the master process in the local
  // communicators involving neighbors only to which this process belongs to
  setMasterGeoLocInLocalComm();
  

  // Output on screen
  if ( 0 == m_rank ) cout << "Local MPI communicators" << endl; 
  for (i=0;i<m_nprocs;++i)
  {
    if ( i == m_rank )
    {      
      cout << "Process = " << m_rank << endl;
      cout << "   Rank in local comm = " << m_rank_localComm << endl;
      cout << "   Coordinates in local comm = " << m_coords[0] << " " <<
      	m_coords[1] << " " << m_coords[2] << endl;      
      cout << "   Nb procs in local comm = " << m_nprocs_localComm << endl;
      cout << "   Rank of master in local comms: ";      
      for (j=0;j<m_nprocs;++j) cout << m_master_localComm[j] << " ";
      cout << endl; 
      cout << "   Geolocalisation of master in local comms: ";      
      for (j=0;j<m_nprocs;++j) 
        if ( m_masterGeoPos[j] != -1 ) cout << 
		Cell::getGeoPositionName(m_masterGeoPos[j]) << " "; 
	else cout << m_masterGeoPos[j] << " ";     
      cout << endl; 
      cout << "   Rank in Activ Comm of master in World Comm = " << 
      	m_rank_masterWorld << endl;
    }
    MPI_Barrier( m_MPI_COMM_activProc );
  } 
  if ( 0 == m_rank ) cout << endl;
} 




// ----------------------------------------------------------------------------
// Definition de la geolocalisation du master dans les communicateurs
// lies aux voisins auquel appartient le proc
void GrainsMPIWrapper::setMasterGeoLocInLocalComm()
{
  int i,j,k,globalRank,ii;
  
  m_masterGeoPos.reserve(m_nprocs);
  for (i=0;i<m_nprocs;++i) m_masterGeoPos.push_back(-1); 
  
  // Geolocalisation du master pour ses voisins sur le proc
  // ------------------------------------------------------
  int *geoPosMasterOnProc = new int[m_nprocs];
  for (i=0;i<m_nprocs;++i) geoPosMasterOnProc[i] = -1;
  for (i=-1;i<2;++i)
    for (j=-1;j<2;++j)    
      for (k=-1;k<2;++k)
        if ( i || j || k )
	{
	  globalRank = m_neighbors->rank( i, j, k );
	  if ( globalRank != -1 ) 
	    geoPosMasterOnProc[globalRank] = 
	    	getGeoPosition( -i, -j, -k );	  
	}  

  // Communication des geolocalisations des masters
  // ----------------------------------------------
  int *recvcounts = new int[m_nprocs];
  for (i=0; i<m_nprocs; ++i) recvcounts[i] = m_nprocs;
  int *displs = new int[m_nprocs];
  displs[0] = 0; 
  for (i=1; i<m_nprocs; ++i) displs[i] = displs[i-1] + recvcounts[i-1];  
  int recvsize = displs[m_nprocs-1] + recvcounts[m_nprocs-1];

  int *allGeoPosMasterOnProc = new int[recvsize];  
  MPI_Allgatherv( geoPosMasterOnProc, m_nprocs, MPI_INT, allGeoPosMasterOnProc, 
  	recvcounts, displs, MPI_INT, m_MPI_COMM_activProc ); 
	
  // Stockage dans masterGeoLoc
  // --------------------------
  for (i=0;i<m_nprocs*m_nprocs;++i)
  {
    ii = int(i/m_nprocs);
    j = i - ii * m_nprocs;
    if ( j == m_rank )
      if ( allGeoPosMasterOnProc[i] != -1 ) 
        m_masterGeoPos[ii] = allGeoPosMasterOnProc[i];  
  } 	   
}




// ----------------------------------------------------------------------------
// Definition de particleHalozoneToNeighboringProcs
void GrainsMPIWrapper::setParticleHalozoneToNeighboringProcs()
{
  m_particleHalozoneToNeighboringProcs.reserve(26);
  vector<int> emptyVECINT;
  for (int i=0;i<26;++i) 
    m_particleHalozoneToNeighboringProcs.push_back(emptyVECINT);
  vector<int> *work;
  
  // NORTH => NORTH
  work = new vector<int>(1,0);
  (*work)[0] = GEOPOS_NORTH;
  m_particleHalozoneToNeighboringProcs[GEOPOS_NORTH] = *work;
  work->clear();
  delete work;
  
  // NORTH_EAST => NORTH, EAST, NORTH_EAST
  work = new vector<int>(3,0);
  (*work)[0] = GEOPOS_NORTH;
  (*work)[1] = GEOPOS_EAST;  
  (*work)[2] = GEOPOS_NORTH_EAST;
  m_particleHalozoneToNeighboringProcs[GEOPOS_NORTH_EAST] = *work;    
  work->clear();
  delete work;
  
  // NORTH_WEST => NORTH, WEST, NORTH_WEST
  work = new vector<int>(3,0);
  (*work)[0] = GEOPOS_NORTH;
  (*work)[1] = GEOPOS_WEST;  
  (*work)[2] = GEOPOS_NORTH_WEST;
  m_particleHalozoneToNeighboringProcs[GEOPOS_NORTH_WEST] = *work;    
  work->clear();
  delete work;    

  // NORTH_FRONT => NORTH, TOP, NORTH_FRONT
  work = new vector<int>(3,0);
  (*work)[0] = GEOPOS_NORTH;
  (*work)[1] = GEOPOS_FRONT;  
  (*work)[2] = GEOPOS_NORTH_FRONT;
  m_particleHalozoneToNeighboringProcs[GEOPOS_NORTH_FRONT] = *work;    
  work->clear();
  delete work;  

  // NORTH_BEHIND => NORTH, BEHIND, NORTH_BEHIND
  work = new vector<int>(3,0);
  (*work)[0] = GEOPOS_NORTH;
  (*work)[1] = GEOPOS_BEHIND;  
  (*work)[2] = GEOPOS_NORTH_BEHIND;
  m_particleHalozoneToNeighboringProcs[GEOPOS_NORTH_BEHIND] = *work;    
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
  m_particleHalozoneToNeighboringProcs[GEOPOS_NORTH_EAST_FRONT] = *work;    
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
  m_particleHalozoneToNeighboringProcs[GEOPOS_NORTH_EAST_BEHIND] = *work;    
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
  m_particleHalozoneToNeighboringProcs[GEOPOS_NORTH_WEST_FRONT] = *work;    
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
  m_particleHalozoneToNeighboringProcs[GEOPOS_NORTH_WEST_BEHIND] = *work;    
  work->clear();
  delete work; 

  // SOUTH => SOUTH
  work = new vector<int>(1,0);
  (*work)[0] = GEOPOS_SOUTH;
  m_particleHalozoneToNeighboringProcs[GEOPOS_SOUTH] = *work;
  work->clear();
  delete work;
  
  // SOUTH_EAST => SOUTH, EAST, SOUTH_EAST
  work = new vector<int>(3,0);
  (*work)[0] = GEOPOS_SOUTH;
  (*work)[1] = GEOPOS_EAST;  
  (*work)[2] = GEOPOS_SOUTH_EAST;
  m_particleHalozoneToNeighboringProcs[GEOPOS_SOUTH_EAST] = *work;    
  work->clear();
  delete work;
  
  // SOUTH_WEST => SOUTH, WEST, SOUTH_WEST
  work = new vector<int>(3,0);
  (*work)[0] = GEOPOS_SOUTH;
  (*work)[1] = GEOPOS_WEST;  
  (*work)[2] = GEOPOS_SOUTH_WEST;
  m_particleHalozoneToNeighboringProcs[GEOPOS_SOUTH_WEST] = *work;    
  work->clear();
  delete work;    

  // SOUTH_FRONT => SOUTH, TOP, SOUTH_FRONT
  work = new vector<int>(3,0);
  (*work)[0] = GEOPOS_SOUTH;
  (*work)[1] = GEOPOS_FRONT;  
  (*work)[2] = GEOPOS_SOUTH_FRONT;
  m_particleHalozoneToNeighboringProcs[GEOPOS_SOUTH_FRONT] = *work;    
  work->clear();
  delete work;  

  // SOUTH_BEHIND => SOUTH, BEHIND, SOUTH_BEHIND
  work = new vector<int>(3,0);
  (*work)[0] = GEOPOS_SOUTH;
  (*work)[1] = GEOPOS_BEHIND;  
  (*work)[2] = GEOPOS_SOUTH_BEHIND;
  m_particleHalozoneToNeighboringProcs[GEOPOS_SOUTH_BEHIND] = *work;    
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
  m_particleHalozoneToNeighboringProcs[GEOPOS_SOUTH_EAST_FRONT] = *work;    
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
  m_particleHalozoneToNeighboringProcs[GEOPOS_SOUTH_EAST_BEHIND] = *work;    
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
  m_particleHalozoneToNeighboringProcs[GEOPOS_SOUTH_WEST_FRONT] = *work;    
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
  m_particleHalozoneToNeighboringProcs[GEOPOS_SOUTH_WEST_BEHIND] = *work;    
  work->clear();
  delete work; 
  
  // EAST => EAST
  work = new vector<int>(1,0);
  (*work)[0] = GEOPOS_EAST;
  m_particleHalozoneToNeighboringProcs[GEOPOS_EAST] = *work;
  work->clear();
  delete work;  

  // WEST => WEST
  work = new vector<int>(1,0);
  (*work)[0] = GEOPOS_WEST;
  m_particleHalozoneToNeighboringProcs[GEOPOS_WEST] = *work;
  work->clear();
  delete work;  

  // EAST_FRONT => EAST, TOP, EAST_FRONT
  work = new vector<int>(3,0);
  (*work)[0] = GEOPOS_EAST;
  (*work)[1] = GEOPOS_FRONT;  
  (*work)[2] = GEOPOS_EAST_FRONT;
  m_particleHalozoneToNeighboringProcs[GEOPOS_EAST_FRONT] = *work;    
  work->clear();
  delete work;  

  // EAST_BEHIND => EAST, BEHIND, EAST_BEHIND
  work = new vector<int>(3,0);
  (*work)[0] = GEOPOS_EAST;
  (*work)[1] = GEOPOS_BEHIND;  
  (*work)[2] = GEOPOS_EAST_BEHIND;
  m_particleHalozoneToNeighboringProcs[GEOPOS_EAST_BEHIND] = *work;    
  work->clear();
  delete work; 

  // WEST_FRONT => WEST, TOP, WEST_FRONT
  work = new vector<int>(3,0);
  (*work)[0] = GEOPOS_WEST;
  (*work)[1] = GEOPOS_FRONT;  
  (*work)[2] = GEOPOS_WEST_FRONT;
  m_particleHalozoneToNeighboringProcs[GEOPOS_WEST_FRONT] = *work;    
  work->clear();
  delete work;  

  // WEST_BEHIND => WEST, BEHIND, WEST_BEHIND
  work = new vector<int>(3,0);
  (*work)[0] = GEOPOS_WEST;
  (*work)[1] = GEOPOS_BEHIND;  
  (*work)[2] = GEOPOS_WEST_BEHIND;
  m_particleHalozoneToNeighboringProcs[GEOPOS_WEST_BEHIND] = *work;    
  work->clear();
  delete work; 
  
  // TOP => TOP
  work = new vector<int>(1,0);
  (*work)[0] = GEOPOS_FRONT;
  m_particleHalozoneToNeighboringProcs[GEOPOS_FRONT] = *work;
  work->clear();
  delete work;  

  // BEHIND => BEHIND
  work = new vector<int>(1,0);
  (*work)[0] = GEOPOS_BEHIND;
  m_particleHalozoneToNeighboringProcs[GEOPOS_BEHIND] = *work;
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
// Returns the process rank in the MPI_COMM_WORLD communicator
int GrainsMPIWrapper::get_rank_world() const
{
  return ( m_rank_world );
}




// ----------------------------------------------------------------------------
// Returns the process rank in the MPI_COMM_activProc communicator
int GrainsMPIWrapper::get_rank_active() const
{
  return ( m_rank );
}




// ----------------------------------------------------------------------------
// Returns whether the process is active
bool GrainsMPIWrapper::isActive() const
{
  return ( m_is_activ );
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
  if ( m_rank == 0 )
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
// Gathers all particles on the master process for post-processing purposes
vector<Particle*>* GrainsMPIWrapper::GatherParticles_PostProcessing(
  	list<Particle*> const& particles,
	list<Particle*> const& pwait,
	vector<Particle*> const& referenceParticles,
	size_t const& nb_total_particles ) const
{
  vector<Particle*>* allparticles = NULL;
  
  // TO DO
  
  return ( allparticles );
}   




// ----------------------------------------------------------------------------
// Gathers all particle velocity-position data on the master process for 
// post-processing purposes
vector< vector<double> >* GrainsMPIWrapper::
    GatherPositionVelocity_PostProcessing(
        list<Particle*> const& particles,
	size_t const& nb_total_particles ) const
{
  vector< vector<double> >* cinematique_Global = NULL;

  // TO DO
  
  return ( cinematique_Global );
}




// ----------------------------------------------------------------------------
// Gathers the class of all particles on the master process 
vector< vector<double> >* GrainsMPIWrapper::
    GatherParticlesClass_PostProcessing(
    const list<Particle*> &particles,
    const size_t& nb_total_particles ) const
{
  double intTodouble = 0.1 ;
  int i=0, tag_DOUBLE = 1, recvsize = 0, ID_part=0;
  MPI_Status status;
  MPI_Request idreq;

  vector< vector<double> >* class_Global = NULL;
    
  list<Particle*>::const_iterator il;
  int nb_part_loc = int(particles.size());
  
  // We do not care about particles in  halozone
  for (il=particles.begin();il!=particles.end();il++)
    if ((*il)->getTag()==2) nb_part_loc--;

  // Buffer size depend on the number of particles per core
  double *buffer = new double[2*nb_part_loc];
  
  for (il=particles.begin(), i=0; il!=particles.end(); il++, i+=2)
  {
    if( (*il)->getTag()==2 ) i-=2;
    else
    {
      buffer[i] = (*il)->getID() + intTodouble ;
      buffer[i+1] = (*il)->getGeometricType();
    }
  }

  MPI_Isend( buffer, 2*nb_part_loc, MPI_DOUBLE, m_rank_masterWorld,
      tag_DOUBLE, m_MPI_COMM_activProc, &idreq );	

  // Reception by the master process
  // -------------------------------
  if( m_rank == m_rank_masterWorld )
  {
    vector<double> work( nb_total_particles, 0. ) ; 
    class_Global = new vector< vector<double> >( 1, work ) ;

    for (int irank=0; irank<m_nprocs; ++irank)
    {
      // Size of the message
      MPI_Probe( irank, tag_DOUBLE, m_MPI_COMM_activProc, &status );  
      MPI_Get_count( &status, MPI_DOUBLE, &recvsize );

      // Reception of the actual message	
      double *recvbuf_DOUBLE = new double[recvsize];
      MPI_Recv( recvbuf_DOUBLE, recvsize, MPI_DOUBLE, 
          irank, tag_DOUBLE, m_MPI_COMM_activProc, &status );	    
      
      // Copy in class_Global
      for(int j=0; j<recvsize; j+=2)
      {
        ID_part = int(recvbuf_DOUBLE[j]) ;
        (*class_Global)[0][ID_part] = recvbuf_DOUBLE[j+1];
      }
      
      delete [] recvbuf_DOUBLE; 
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
  	list<Particle*> const* particlesHalozone,
  	list<Particle*>* particlesClones,
	vector<Particle*> const* referenceParticles,
	LinkedCell* LC, bool update )
{
  list<Particle*>::const_iterator il;
  int i, j, tag_DOUBLE = 1, recvsize = 0, geoLoc, ireq = 0;
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
  
        
  // Copy particles in halozone data into local buffers
  // --------------------------------------------------
  vector<int> nbHzGeoLoc(26,0);
  vector<int>::iterator iv;
  for (il=particlesHalozone->begin();il!=particlesHalozone->end();il++)
  {
    geoLoc = (*il)->getGeoPosition();
    for (iv=m_particleHalozoneToNeighboringProcs[geoLoc].begin();
    	iv!=m_particleHalozoneToNeighboringProcs[geoLoc].end();iv++)
      nbHzGeoLoc[*iv]++;
  }             

  // Buffer of doubles: kinematics and configuration as follows 
  // [ID number, class, rank of sending process, translational
  //  velocity, rotation quaternion, angular velocity, transform]
  // ------------------------------------------------------------
  int NB_DOUBLE_PART = 25;

  vector<int> index( 26, 0 );
  double *pDOUBLE = NULL;  
  vector<double*> features( 26, pDOUBLE );
  for (i=0;i<26;i++) features[i] = new double[ NB_DOUBLE_PART * nbHzGeoLoc[i] ];
  double ParticleID = 0., ParticleClass = 0.;

  for (il=particlesHalozone->begin(),i=0;il!=particlesHalozone->end();il++)
  {
    geoLoc = (*il)->getGeoPosition();
    ParticleID = (*il)->getID() + intTodouble ;
    ParticleClass = (*il)->getGeometricType() + intTodouble ;
    for (iv=m_particleHalozoneToNeighboringProcs[geoLoc].begin();
    	iv!=m_particleHalozoneToNeighboringProcs[geoLoc].end();iv++)
    {
      j = index[*iv]; 
      features[*iv][j] = ParticleID;             
      features[*iv][j+1] = ParticleClass;
      features[*iv][j+2] = m_rank + intTodouble ;
      (*il)->copyTranslationalVelocity( features[*iv], j+3 );
      (*il)->copyQuaternionRotation( features[*iv], j+6 );    
      (*il)->copyAngularVelocity( features[*iv], j+10 );
      (*il)->copyTransform( features[*iv], j+13, m_MPIperiodes[*iv] );     
      index[*iv] += NB_DOUBLE_PART;
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
    MPI_Isend( features[*ign], nbHzGeoLoc[*ign] * NB_DOUBLE_PART, MPI_DOUBLE, 
	*irn, tag_DOUBLE + m_GeoLocReciprocity[*ign], 
	m_MPI_COMM_activProc, &idreq[ireq] );
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
    MPI_Probe( *irn, tag_DOUBLE + *ign, m_MPI_COMM_activProc, &status );  
    MPI_Get_count( &status, MPI_DOUBLE, &recvsize );
    recvsize /= NB_DOUBLE_PART;  

    // Reception of the actual message	    
    double *recvbuf_DOUBLE = new double[recvsize * NB_DOUBLE_PART];
    MPI_Recv( recvbuf_DOUBLE, recvsize * NB_DOUBLE_PART, MPI_DOUBLE, 
	*irn, tag_DOUBLE + *ign, m_MPI_COMM_activProc, &status );	    

    SCT_add_elapsed_time( "MPIComm" );
    SCT_set_start( "UpdateCreateClones" );
 
    // Creation or update of clone particles
    // -------------------------------------
    if ( update )    
      UpdateClones( time, recvsize, recvbuf_DOUBLE,
		NB_DOUBLE_PART, particlesClones,
		particles, particlesHalozone, referenceParticles, LC );
    else
      CreateClones( time, recvsize, recvbuf_DOUBLE,
		NB_DOUBLE_PART, particlesClones,
		particles, particlesHalozone, referenceParticles, LC );      

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
// Broadcasts a double from the master to all processes within the 
// MPI_COMM_activProc communicator
double GrainsMPIWrapper::Broadcast_DOUBLE( double const& d ) const
{
  double collective_d= d;
  
  MPI_Bcast( &collective_d, 1, MPI_DOUBLE, 0, m_MPI_COMM_activProc );
  
  return ( collective_d ); 
}




// ----------------------------------------------------------------------------
// Broadcasts an integer from the master to all processes within the 
// MPI_COMM_activProc communicator
int GrainsMPIWrapper::Broadcast_INT( int const& i ) const
{
  int collective_i = i;
  
  MPI_Bcast( &collective_i, 1, MPI_INT, 0, m_MPI_COMM_activProc );
  
  return ( collective_i ); 
}




// ----------------------------------------------------------------------------
// Broadcasts an integer from the master to all processes within the 
// MPI_COMM_activProc communicator
size_t GrainsMPIWrapper::Broadcast_UNSIGNED_INT( size_t const& i ) const
{
  size_t collective_i = i;
  
  MPI_Bcast( &collective_i, 1, MPI_UNSIGNED_LONG, 0, m_MPI_COMM_activProc );
  
  return collective_i; 
}




// ----------------------------------------------------------------------------
// Sums a double from all processes on all processes within the 
// MPI_COMM_activProc communicator
double GrainsMPIWrapper::sum_DOUBLE( double const& x ) const
{
  double sum = 0;
  
  MPI_Allreduce( &x, &sum, 1, MPI_DOUBLE, MPI_SUM, m_MPI_COMM_activProc );
  
  return ( sum );
}




// ----------------------------------------------------------------------------
// Sums a double from all processes on the master process within the 
// MPI_COMM_activProc communicator
double GrainsMPIWrapper::sum_DOUBLE_master( double const& x ) const
{
  double sum = 0;
  
  MPI_Reduce( &x, &sum, 1, MPI_DOUBLE, MPI_SUM, m_rank_masterWorld,
  	m_MPI_COMM_activProc );
  
  return ( sum );
}




// ----------------------------------------------------------------------------
// Sums an integer from all processes on all processes within the 
// MPI_COMM_activProc communicator
int GrainsMPIWrapper::sum_INT( int const& i ) const
{
  int sum = 0;
  
  MPI_Allreduce( &i, &sum, 1, MPI_INT, MPI_SUM, m_MPI_COMM_activProc );
  
  return ( sum );
}




// ----------------------------------------------------------------------------
// Sums an integer from all processes on all processes within the 
// MPI_COMM_activProc communicator
size_t GrainsMPIWrapper::sum_UNSIGNED_INT( size_t const& i ) const
{
  size_t sum = 0;
  
  MPI_Allreduce( &i, &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 
  	m_MPI_COMM_activProc );
  
  return ( sum );
}




// ----------------------------------------------------------------------------
// Sums an integer from all processes on the master process within the 
// MPI_COMM_activProc communicator
int GrainsMPIWrapper::sum_INT_master( int const& i ) const
{
  int sum = 0;
  
  MPI_Reduce( &i, &sum, 1, MPI_INT, MPI_SUM, m_rank_masterWorld, 
  	m_MPI_COMM_activProc );
  
  return ( sum );
}




// ----------------------------------------------------------------------------
// Sums an unsigned integer from all processes on the master process 
// within the MPI_COMM_activProc communicator
size_t GrainsMPIWrapper::sum_UNSIGNED_INT_master( size_t const& i ) const
{
  int sum = 0;
  
  MPI_Reduce( &i, &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, m_rank_masterWorld, 
  	m_MPI_COMM_activProc );
  
  return ( sum );
}




// ----------------------------------------------------------------------------
// Performs a "logical and" operation on input boolean value from 
// all processes on all processes
bool GrainsMPIWrapper::logical_and( bool const& input ) const
{
  int land=0;
  
  MPI_Allreduce( &input, &land, 1, MPI_UNSIGNED_SHORT, MPI_LAND,
  	m_MPI_COMM_activProc );
  
  return ( land );
}




// ----------------------------------------------------------------------------
// Minimum of an integer from all processes on all processes within the 
// MPI_COMM_activProc communicator
int GrainsMPIWrapper::min_INT( int const& i ) const
{
  int collective_i = 0;

  MPI_Allreduce( &i, &collective_i, 1, MPI_INT, MPI_MIN, 
  	m_MPI_COMM_activProc ) ;

  return ( collective_i );
}  




// ----------------------------------------------------------------------------
// Max d'un entier sur tous les proc
int GrainsMPIWrapper::max_INT( int const& i ) const
{
  int collective_i = 0;

  MPI_Allreduce( &i, &collective_i, 1, MPI_INT, MPI_MAX, 
  	m_MPI_COMM_activProc ) ;

  return ( collective_i );
}  




// ----------------------------------------------------------------------------
// Min d'un entier non signe sur tous les proc
size_t GrainsMPIWrapper::min_UNSIGNED_INT( size_t const& i ) const
{
  size_t collective_i = 0;

  MPI_Allreduce( &i, &collective_i, 1, MPI_UNSIGNED_LONG, MPI_MIN, 
  	m_MPI_COMM_activProc ) ;

  return ( collective_i );
}  




// ----------------------------------------------------------------------------
// Maximum of an unsigned integer from all processes on all processes within 
// the MPI_COMM_activProc communicator
size_t GrainsMPIWrapper::max_UNSIGNED_INT( size_t const& i ) const
{
  size_t collective_i = 0;

  MPI_Allreduce( &i, &collective_i, 1, MPI_UNSIGNED_LONG, MPI_MAX, 
  	m_MPI_COMM_activProc ) ;

  return ( collective_i );
}  




// ----------------------------------------------------------------------------
// Maximum of a double from all processes on all processes within the 
// MPI_COMM_activProc communicator
double GrainsMPIWrapper::max_DOUBLE( double const& x ) const 
{
  double max = 0.;

  MPI_Allreduce( &x, &max, 1, MPI_DOUBLE, MPI_MAX, m_MPI_COMM_activProc ) ;

  return ( max );
} 



 
// ----------------------------------------------------------------------------
// Maximum of a double from all processes on the master process within the 
// MPI_COMM_activProc communicator
double GrainsMPIWrapper::max_DOUBLE_master( double const& x ) const 
{
  double max = 0.;

  MPI_Reduce( &x, &max, 1, MPI_DOUBLE, MPI_MAX, m_rank_masterWorld,
  	m_MPI_COMM_activProc ) ;

  return ( max );
} 




// ----------------------------------------------------------------------------
// Minimum of a double from all processes on all processes within the 
// MPI_COMM_activProc communicator
double GrainsMPIWrapper::min_DOUBLE( double const& x ) const 
{
  double min = 0.;

  MPI_Allreduce( &x, &min, 1, MPI_DOUBLE, MPI_MIN, m_MPI_COMM_activProc ) ;

  return ( min );
} 




// ----------------------------------------------------------------------------
// Minimum of a double from all processes on the master process within the 
// MPI_COMM_activProc communicator
double GrainsMPIWrapper::min_DOUBLE_master( double const& x ) const 
{
  double min = 0.;

  MPI_Reduce( &x, &min, 1, MPI_DOUBLE, MPI_MIN, m_rank_masterWorld,
  	m_MPI_COMM_activProc ) ;

  return ( min );
} 




// ----------------------------------------------------------------------------
// AllGather of an unsigned integer from all processes on all processes within 
// the MPI_COMM_activProc communicator
size_t* GrainsMPIWrapper::AllGather_UNSIGNED_INT( size_t const& i ) const
{
  size_t* recv = new size_t[m_nprocs];
  
  MPI_Allgather( &i, 1, MPI_UNSIGNED_LONG, recv, 1, MPI_UNSIGNED_LONG, 
  	m_MPI_COMM_activProc ); 

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

  MPI_Bcast( coordinates, 3, MPI_DOUBLE, 0, m_MPI_COMM_activProc );  
  
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

  MPI_Bcast( coordinates, 3, MPI_DOUBLE, 0, m_MPI_COMM_activProc );  
  
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

  MPI_Bcast( mat_coef, 9, MPI_DOUBLE, 0, m_MPI_COMM_activProc );  
  
  Matrix bmat( mat_coef[0], mat_coef[1], mat_coef[2],
  	mat_coef[3], mat_coef[4], mat_coef[5],
	mat_coef[6], mat_coef[7], mat_coef[8]);
  delete [] mat_coef;
  
  return ( bmat );
}




// ----------------------------------------------------------------------------
// Test avec communicateur local: validation
void GrainsMPIWrapper::testCommLocal_Gather_INT() const
{
  int nb_hz = - m_rank, i, j, ii;
  int *recvbuf_INT = new int[m_nprocs_localComm];
  for (i=0; i<m_nprocs_localComm; ++i) recvbuf_INT[i]=0;   
  for (ii=0;ii<m_nprocs;++ii)  
    if ( m_isInCommMPINeighbors[ii] )
      MPI_Gather( &nb_hz, 1, MPI_INT, recvbuf_INT, 1, MPI_INT,
  	m_master_localComm[ii], *m_commMPINeighbors[ii] ); 
    

  for (i=0;i<m_nprocs;++i)
  {
    if ( i == m_rank )
    {      
      cout << "Test comm local : "; 
      for (j=0;j<m_nprocs_localComm;++j) cout << recvbuf_INT[j] << " ";
      cout << endl;        
    }
    MPI_Barrier( m_MPI_COMM_activProc );
  } 
  
  delete [] recvbuf_INT;
}




// ----------------------------------------------------------------------------
// AllGather of a 1D array of n integers within the MPI_COMM_activProc 
// communicator
void GrainsMPIWrapper::test_AllGatherv_INT( int const& n ) const
{
  int i;

  // Local array
  int *local_vec = new int[n];
  for (i=0;i<n;++i) local_vec[i] = 100000 * m_rank + i;

  // Partition and size of the reception buffer
  int *recvcounts_INT = new int[m_nprocs];
  for (i=0; i<m_nprocs; ++i) recvcounts_INT[i] = n;
  int *displs_INT = new int[m_nprocs];
  displs_INT[0] = 0; 
  for (i=1; i<m_nprocs; ++i) 
    displs_INT[i] = displs_INT[i-1] + recvcounts_INT[i-1];  
  int recvsize_INT = displs_INT[m_nprocs-1] + recvcounts_INT[m_nprocs-1];
  int *recvbuf_INT = new int[recvsize_INT];
  
  // MPI communication 
  MPI_Allgatherv( local_vec, n, MPI_INT, recvbuf_INT, 
  	recvcounts_INT, displs_INT, MPI_INT, m_MPI_COMM_activProc );
	
  delete [] local_vec;
  delete [] recvcounts_INT;
  delete [] displs_INT;
  delete [] recvbuf_INT; 
}




// ----------------------------------------------------------------------------
// AllGather of a 1D array of n integers within the m_commMPINeighbors local 
// communicator
void GrainsMPIWrapper::testCommLocal_AllGatherv_INT( const int &n ) const
{
  int i,ii;

  // Local array
  int *local_vec = new int[n];
  for (i=0;i<n;++i) local_vec[i] = 100000 * m_rank + i;

  // Partition and size of the reception buffer
  int *recvcounts_INT = new int[m_nprocs_localComm];
  for (i=0; i<m_nprocs_localComm; ++i) recvcounts_INT[i] = n;
  int *displs_INT = new int[m_nprocs_localComm];
  displs_INT[0] = 0; 
  for (i=1; i<m_nprocs_localComm; ++i) 
    displs_INT[i] = displs_INT[i-1] + recvcounts_INT[i-1];  
  int recvsize_INT = displs_INT[m_nprocs_localComm-1] 
  	+ recvcounts_INT[m_nprocs_localComm-1];
  int *recvbuf_INT = new int[recvsize_INT];
  
  // MPI communication 
  for (ii=0;ii<m_nprocs;++ii)  
    if ( m_isInCommMPINeighbors[ii] )  
      MPI_Gatherv( local_vec, n, MPI_INT, recvbuf_INT, 
  	recvcounts_INT, displs_INT, MPI_INT, m_master_localComm[ii],
	*m_commMPINeighbors[ii] );
	
  delete [] local_vec;
  delete [] recvcounts_INT;
  delete [] displs_INT;
  delete [] recvbuf_INT; 
}




// ----------------------------------------------------------------------------
// AllGather of a 1D array of n doubles within the MPI_COMM_activProc 
// communicator
void GrainsMPIWrapper::test_AllGatherv_DOUBLE( int const& n ) const
{
  int i;

  // Local array
  double *local_vec = new double[n];
  for (i=0;i<n;++i) local_vec[i] = 20000. * double(m_rank) + 1.1 * double(i);

  // Partition and size of the reception buffer
  int *recvcounts_DOUBLE = new int[m_nprocs];
  for (i=0; i<m_nprocs; ++i) recvcounts_DOUBLE[i] = n;
  int *displs_DOUBLE = new int[m_nprocs];
  displs_DOUBLE[0] = 0; 
  for (i=1; i<m_nprocs; ++i) 
    displs_DOUBLE[i] = displs_DOUBLE[i-1] + recvcounts_DOUBLE[i-1];  
  int recvsize_DOUBLE = displs_DOUBLE[m_nprocs-1] 
  	+ recvcounts_DOUBLE[m_nprocs-1];
  double *recvbuf_DOUBLE = new double[recvsize_DOUBLE];
  
  // MPI communication 
  MPI_Allgatherv( local_vec, n, MPI_DOUBLE, recvbuf_DOUBLE, 
  	recvcounts_DOUBLE, displs_DOUBLE, MPI_DOUBLE, m_MPI_COMM_activProc );
	
  delete [] local_vec;
  delete [] recvcounts_DOUBLE;
  delete [] displs_DOUBLE;
  delete [] recvbuf_DOUBLE; 
}




// ----------------------------------------------------------------------------
// AllGather of a 1D array of n doubles within the m_commMPINeighbors local 
// communicator
void GrainsMPIWrapper::testCommLocal_AllGatherv_DOUBLE( int const& n ) const
{
  int i,ii;

  // Local array
  double *local_vec = new double[n];
  for (i=0;i<n;++i) local_vec[i] = 20000. * double(m_rank) + 1.1 * double(i);

  // Partition and size of the reception buffer
  int *recvcounts_DOUBLE = new int[m_nprocs_localComm];
  for (i=0; i<m_nprocs_localComm; ++i) recvcounts_DOUBLE[i] = n;
  int *displs_DOUBLE = new int[m_nprocs_localComm];
  displs_DOUBLE[0] = 0; 
  for (i=1; i<m_nprocs_localComm; ++i) 
    displs_DOUBLE[i] = displs_DOUBLE[i-1] + recvcounts_DOUBLE[i-1];  
  int recvsize_DOUBLE = displs_DOUBLE[m_nprocs_localComm-1] 
  	+ recvcounts_DOUBLE[m_nprocs_localComm-1];
  double *recvbuf_DOUBLE = new double[recvsize_DOUBLE];
  
  // MPI communication 
  for (ii=0;ii<m_nprocs;++ii)  
    if ( m_isInCommMPINeighbors[ii] )  
      MPI_Gatherv( local_vec, n, MPI_DOUBLE, recvbuf_DOUBLE, 
  	recvcounts_DOUBLE, displs_DOUBLE, MPI_DOUBLE, m_master_localComm[ii],
	*m_commMPINeighbors[ii] );
	
  delete [] local_vec;
  delete [] recvcounts_DOUBLE;
  delete [] displs_DOUBLE;
  delete [] recvbuf_DOUBLE; 
}




// ----------------------------------------------------------------------------
// Send-Recv of a 1D array of n doubles within the m_commMPINeighbors local 
// communicator
void GrainsMPIWrapper::testCommLocal_SendRecv_DOUBLE( int const& n ) const
{
  int i,ii;
  MPI_Status status;
  int nloc = n + m_rank;

  // Local array
  double *local_vec = new double[nloc];
  for (i=0;i<nloc;++i) local_vec[i] = 20000. * double(m_rank) + 1.1 * double(i);
  
  // MPI communication 
  for (ii=0;ii<m_nprocs;++ii)  
    if ( m_isInCommMPINeighbors[ii] )
    {
      if ( ii != m_rank )
        MPI_Ssend( local_vec, nloc, MPI_DOUBLE, m_master_localComm[ii], 0, 
		*m_commMPINeighbors[ii] );	
      else
        for (i=0;i<m_nprocs_localComm;++i)
          if ( i != m_rank_localComm )
          {
            MPI_Probe( i, 0, *m_commMPINeighbors[ii], &status );
	    int sizemess = 0;
	    MPI_Get_count( &status, MPI_DOUBLE, &sizemess );
	    double *recvbuf_DOUBLE = new double[sizemess];
	    MPI_Recv( recvbuf_DOUBLE, sizemess, MPI_DOUBLE, i, 0, 
		*m_commMPINeighbors[ii], &status );	    
            delete [] recvbuf_DOUBLE;	      
          }               
    }

  delete [] local_vec;
}




// ----------------------------------------------------------------------------
// Outputs timer summary
void GrainsMPIWrapper::timerSummary() const
{
  double cputime=SCT_get_total_elapsed_time("BuffersCopy")
     	+SCT_get_total_elapsed_time("MPIComm")
     	+SCT_get_total_elapsed_time("UpdateCreateClones");
  cout << "MPI wrapper timer" << endl;
  SCT_get_summary(cout,cputime);
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
//   for (int i=0;i<get_total_number_of_active_processes();++i)
//   {
//     if ( i == m_rank )
//       if ( !m_MPILogString->empty() )
//       {
//         f << "Processor = " << i << endl << std::flush;
// 	f << *m_MPILogString << std::flush;
//       }
//     MPI_Barrier( m_MPI_COMM_activProc );    
//   }
  
  writeStringPerProcess( f, *m_MPILogString, true );

  delete m_MPILogString;
  m_MPILogString = new string;
}




// ----------------------------------------------------------------------------
// MPI_Barrier pour les processus actifs uniquement
void GrainsMPIWrapper::MPI_Barrier_ActivProc() const
{
  MPI_Barrier( m_MPI_COMM_activProc );

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
  	1, MPI_DOUBLE, m_MPI_COMM_activProc ); 
  MPI_Allgather( &time_overlapMax, 1, MPI_DOUBLE, recvbuf_time,
  	1, MPI_DOUBLE, m_MPI_COMM_activProc );	 
	
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
// Creates and updates clones with the data sent by the neighboring processes 
void GrainsMPIWrapper::UpdateOrCreateClones(double time,
 	int const &recvsize, double const* recvbuf_DOUBLE,
	const int& NB_DOUBLE_PART,  
  	list<Particle*>* particlesClones,
	list<Particle*>* particles,
  	list<Particle*> const* particlesHalozone,
	vector<Particle*> const* referenceParticles,
	LinkedCell* LC )
{
  int j, id, classe;
  bool found = false;
  double distGC = 0. ; 
  Point3 const* GC = NULL ; 
  multimap<int,Particle*>::iterator imm;
  size_t ncid = 0;
  Particle* pClone = NULL ;
  pair < multimap<int,Particle*>::iterator, 
  	multimap<int,Particle*>::iterator > crange;
  
  for( j=0; j<recvsize; ++j )
  {
    found = false;
    id = int( recvbuf_DOUBLE[NB_DOUBLE_PART*j] );  
      
    // Searh whether the clone already exists in this local domain
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
          	pow( recvbuf_DOUBLE[NB_DOUBLE_PART*j+22] - (*GC)[X], 2. ) +
            pow( recvbuf_DOUBLE[NB_DOUBLE_PART*j+23] - (*GC)[Y], 2. ) +
            pow( recvbuf_DOUBLE[NB_DOUBLE_PART*j+24] - (*GC)[Z], 2. ) ) ;
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
          pow( recvbuf_DOUBLE[NB_DOUBLE_PART*j+22] - (*GC)[X], 2. ) +
          pow( recvbuf_DOUBLE[NB_DOUBLE_PART*j+23] - (*GC)[Y], 2. ) +
          pow( recvbuf_DOUBLE[NB_DOUBLE_PART*j+24] - (*GC)[Z], 2. ) ) ;
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
      
      pClone->setPosition( &recvbuf_DOUBLE[NB_DOUBLE_PART*j+13] );
      Vector3 trans( recvbuf_DOUBLE[NB_DOUBLE_PART*j+3],
	  	recvbuf_DOUBLE[NB_DOUBLE_PART*j+4],
		recvbuf_DOUBLE[NB_DOUBLE_PART*j+5] );
      pClone->setTranslationalVelocity( trans );
      pClone->setQuaternionRotation( recvbuf_DOUBLE[NB_DOUBLE_PART*j+6],
		recvbuf_DOUBLE[NB_DOUBLE_PART*j+7],
		recvbuf_DOUBLE[NB_DOUBLE_PART*j+8],
		recvbuf_DOUBLE[NB_DOUBLE_PART*j+9] );
      Vector3 rot(recvbuf_DOUBLE[NB_DOUBLE_PART*j+10],
	  	recvbuf_DOUBLE[NB_DOUBLE_PART*j+11],
		recvbuf_DOUBLE[NB_DOUBLE_PART*j+12] );
      pClone->setAngularVelocity( rot );
    }
    else // is not found
    {
      if ( LC->isInLinkedCell( recvbuf_DOUBLE[NB_DOUBLE_PART*j+22],
                              recvbuf_DOUBLE[NB_DOUBLE_PART*j+23],
                              recvbuf_DOUBLE[NB_DOUBLE_PART*j+24] ) &&
          ( int( recvbuf_DOUBLE[NB_DOUBLE_PART*j+2] ) != m_rank ||
            m_isMPIperiodic ) )
      {
        classe = int( recvbuf_DOUBLE[NB_DOUBLE_PART*j+1] );
	
        if ( GrainsExec::m_MPI_verbose )
        {
          ostringstream oss;
          oss << "   t=" << GrainsExec::doubleToString(time,TIMEFORMAT)
              << " Create Clone                                Id = " 
              << id
              << " Classe = " << classe << " " 
              << recvbuf_DOUBLE[NB_DOUBLE_PART*j+25] << " " 
              << recvbuf_DOUBLE[NB_DOUBLE_PART*j+26] << " " 
              << recvbuf_DOUBLE[NB_DOUBLE_PART*j+27]
              << endl;
          GrainsMPIWrapper::addToMPIString(oss.str());
        }

        // Creation du clone
        Particle *new_clone = NULL ;
        if ( (*referenceParticles)[classe]->isCompositeParticle() )
          new_clone = new CompositeParticle( id,
              (*referenceParticles)[classe],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+3],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+4],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+5],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+6],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+7],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+8],	
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+9],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+10],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+11],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+12],
              &recvbuf_DOUBLE[NB_DOUBLE_PART*j+13],
              COMPUTE, 2, 0 );  
        else
          new_clone = new Particle( id,
              (*referenceParticles)[classe],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+3],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+4],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+5],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+6],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+7],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+8],	
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+9],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+10],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+11],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+12],
              &recvbuf_DOUBLE[NB_DOUBLE_PART*j+13],
              COMPUTE, 2, 0 );

        // Add to the LinkedCell
        LC->Link( new_clone );
	
        // Add to active particle and clone particle lists
        particlesClones->push_back( new_clone );
        particles->push_back( new_clone ); 
	
        // Add to the clone map
        m_AccessToClones.insert( pair<int,Particle*>( id, new_clone ) );
      }    
    }    
  } 
}




// ----------------------------------------------------------------------------
// Updates clones with the data sent by the neighboring processes 
void GrainsMPIWrapper::UpdateClones(double time,
 	int const &recvsize, double const* recvbuf_DOUBLE,
	const int& NB_DOUBLE_PART,  
  	list<Particle*>* particlesClones,
	list<Particle*>* particles,
  	list<Particle*> const* particlesHalozone,
	vector<Particle*> const* referenceParticles,
	LinkedCell* LC )
{
  int j, id;
  bool found = false;
  double distGC = 0. ; 
  Point3 const* GC = NULL ; 
  multimap<int,Particle*>::iterator imm;
  size_t ncid = 0;
  Particle* pClone = NULL ;
  pair < multimap<int,Particle*>::iterator, 
  	multimap<int,Particle*>::iterator > crange;
  
  for( j=0; j<recvsize; ++j )
  {
    found = false;
    id = int( recvbuf_DOUBLE[NB_DOUBLE_PART*j] );  
      
    // Searh whether the clone already exists in this local domain
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
          	pow( recvbuf_DOUBLE[NB_DOUBLE_PART*j+22] - (*GC)[X], 2. ) +
            	pow( recvbuf_DOUBLE[NB_DOUBLE_PART*j+23] - (*GC)[Y], 2. ) +
            	pow( recvbuf_DOUBLE[NB_DOUBLE_PART*j+24] - (*GC)[Z], 2. ) ) ;
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
          	pow( recvbuf_DOUBLE[NB_DOUBLE_PART*j+22] - (*GC)[X], 2. ) +
          	pow( recvbuf_DOUBLE[NB_DOUBLE_PART*j+23] - (*GC)[Y], 2. ) +
          	pow( recvbuf_DOUBLE[NB_DOUBLE_PART*j+24] - (*GC)[Z], 2. ) ) ;
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
      
      pClone->setPosition( &recvbuf_DOUBLE[NB_DOUBLE_PART*j+13] );
      Vector3 trans( recvbuf_DOUBLE[NB_DOUBLE_PART*j+3],
	  	recvbuf_DOUBLE[NB_DOUBLE_PART*j+4],
		recvbuf_DOUBLE[NB_DOUBLE_PART*j+5] );
      pClone->setTranslationalVelocity( trans );
      pClone->setQuaternionRotation( recvbuf_DOUBLE[NB_DOUBLE_PART*j+6],
		recvbuf_DOUBLE[NB_DOUBLE_PART*j+7],
		recvbuf_DOUBLE[NB_DOUBLE_PART*j+8],
		recvbuf_DOUBLE[NB_DOUBLE_PART*j+9] );
      Vector3 rot(recvbuf_DOUBLE[NB_DOUBLE_PART*j+10],
	  	recvbuf_DOUBLE[NB_DOUBLE_PART*j+11],
		recvbuf_DOUBLE[NB_DOUBLE_PART*j+12] );
      pClone->setAngularVelocity( rot );
    }
    else
    {
      cout << "!!! Warning !!! Clone " << id << " not found in proc " 
      	<< m_rank << endl;
    }  
  } 
}




// ----------------------------------------------------------------------------
// Creates and updates clones with the data sent by the neighboring processes 
void GrainsMPIWrapper::CreateClones(double time,
 	int const &recvsize, double const* recvbuf_DOUBLE,
	const int& NB_DOUBLE_PART,  
  	list<Particle*>* particlesClones,
	list<Particle*>* particles,
  	list<Particle*> const* particlesHalozone,
	vector<Particle*> const* referenceParticles,
	LinkedCell* LC )
{
  int j, id, classe;
  bool found = false;
  double distGC = 0. ; 
  Point3 const* GC = NULL ; 
  multimap<int,Particle*>::iterator imm;
  size_t ncid = 0;
  pair < multimap<int,Particle*>::iterator, 
  	multimap<int,Particle*>::iterator > crange;
	  
  // Note: although we create new clones here, we need to check if they
  // do not already exist. This scenario happens if a master particle moves from
  // a halozone cell to another halozone cell and the two halozone cells have
  // a different GeoPosition, then such a master particle is added to the list 
  // of new halozone particle but might already exist on the local process

  for( j=0; j<recvsize; ++j )
  {
    found = false;
    id = int( recvbuf_DOUBLE[NB_DOUBLE_PART*j] );  
      
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
          	pow( recvbuf_DOUBLE[NB_DOUBLE_PART*j+22] - (*GC)[X], 2. ) +
            pow( recvbuf_DOUBLE[NB_DOUBLE_PART*j+23] - (*GC)[Y], 2. ) +
            pow( recvbuf_DOUBLE[NB_DOUBLE_PART*j+24] - (*GC)[Z], 2. ) ) ;
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
          pow( recvbuf_DOUBLE[NB_DOUBLE_PART*j+22] - (*GC)[X], 2. ) +
          pow( recvbuf_DOUBLE[NB_DOUBLE_PART*j+23] - (*GC)[Y], 2. ) +
          pow( recvbuf_DOUBLE[NB_DOUBLE_PART*j+24] - (*GC)[Z], 2. ) ) ;
          if ( distGC < 1.1 * imm->second->getCrustThickness() )  
            found = true;
          else imm++;
        }
        break;    
    }   
       
    // If not found, create the clone particle
    if ( !found )
    {
      classe = int( recvbuf_DOUBLE[NB_DOUBLE_PART*j+1] );
	
      if ( GrainsExec::m_MPI_verbose )
      {
        ostringstream oss;
        oss << "   t=" << GrainsExec::doubleToString(time,TIMEFORMAT)
		<< " Create Clone                                Id = " 
		<< id
		<< " Classe = " << classe << " " 
		<< recvbuf_DOUBLE[NB_DOUBLE_PART*j+22] << " " 
		<< recvbuf_DOUBLE[NB_DOUBLE_PART*j+23] << " " 
		<< recvbuf_DOUBLE[NB_DOUBLE_PART*j+24]
		<< endl;
        GrainsMPIWrapper::addToMPIString(oss.str());
      }

      // Creation of the clone particle
      Particle *new_clone = NULL ;
      if ( (*referenceParticles)[classe]->isCompositeParticle() )
        new_clone = new CompositeParticle( id,
              (*referenceParticles)[classe],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+3],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+4],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+5],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+6],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+7],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+8],	
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+9],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+10],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+11],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+12],
              &recvbuf_DOUBLE[NB_DOUBLE_PART*j+13],
              COMPUTE, 2, 0 );   
      else
        new_clone = new Particle( id,
              (*referenceParticles)[classe],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+3],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+4],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+5],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+6],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+7],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+8],	
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+9],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+10],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+11],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+12],
              &recvbuf_DOUBLE[NB_DOUBLE_PART*j+13],
              COMPUTE, 2, 0 );

      // Add to the LinkedCell
      LC->Link( new_clone );
	
      // Add to active particle and clone particle lists
      particlesClones->push_back( new_clone );
      particles->push_back( new_clone ); 
    }
    else m_AccessToClones.erase( imm ); 
  } 
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
  if ( m_rank == m_rank_masterWorld )
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
  	MPI_SUM, m_rank_masterWorld, m_MPI_COMM_activProc ); 
	
  // Add all contributions from other processes on the master
  if ( m_rank == m_rank_masterWorld )
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
  if ( m_rank == 0 ) f << "Memory used by Grains3D" << endl;
  
  ostringstream out;
  GrainsExec::display_memory( out, GrainsExec::used_memory() ); 
  writeStringPerProcess( f, out.str(), false, GrainsExec::m_shift3 );   
}




// ----------------------------------------------------------------------------
// Distributes the number of particles in each class and on each
// process in the case of the block structured insertion
void GrainsMPIWrapper::distributeParticlesClassProc( 
  	list< pair<Particle*,int> > const& newPart,
	list< pair<Particle*,int> >&newPartProc,
	size_t const& npartproc,
	size_t const& ntotalinsert ) const
{
  list< pair<Particle*,int> >::const_iterator ipart;
  list< pair<Particle*,int> >::iterator ipartProc;
  int nclasseproc = 0, npartproc_ = int(npartproc), i, 
  	nbClasses = int(newPart.size()) ;
  MPI_Status status;
  

  // Initialisation of lists and arrays
  newPartProc = newPart;
  int* tabNbPartRestantesParClasse = new int[nbClasses];
  for (i=0,ipart=newPart.begin();i<nbClasses;++i,ipart++) 
    tabNbPartRestantesParClasse[i] = ipart->second;
  
  
  // Distribution is performed sequentially from 0 to the last
  // process
  if ( m_rank == 0 )
  {
    // All types except the last type
    ipart = newPart.begin();
    ipartProc = newPartProc.begin();  
    for (i=0,ipart=newPart.begin(),ipartProc=newPartProc.begin();
    	i<nbClasses-1;++i,ipart++,ipartProc++)
    {
      nclasseproc = int( npartproc * ipart->second / ntotalinsert ) ;
      // If ratio is integer, we keep it, otherwise we add 1
      if ( fabs( double( npartproc * ipart->second / ntotalinsert )
      	- double( int( npartproc * ipart->second / ntotalinsert ) ) ) > 1.e-14 )
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
  
    MPI_Send( tabNbPartRestantesParClasse, nbClasses, MPI_INT, 1, m_rank, 
  	m_MPI_COMM_activProc );
  }
  else
  {
    if ( m_rank != m_nprocs - 1 )
    {
      MPI_Recv( tabNbPartRestantesParClasse, nbClasses, MPI_INT, m_rank - 1,
    	m_rank - 1, m_MPI_COMM_activProc, &status );

    // All types except the last type
      for (i=0,ipart=newPart.begin(),ipartProc=newPartProc.begin();
    	i<nbClasses-1;++i,ipart++,ipartProc++)      
      {
        nclasseproc = int( npartproc * ipart->second / ntotalinsert ) ;
        // If ratio is integer, we keep it, otherwise we add 1
        if ( fabs( double( npartproc * ipart->second / ntotalinsert )
      	- double( int( npartproc * ipart->second / ntotalinsert ) ) ) > 1.e-14 )
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
  
      MPI_Send( tabNbPartRestantesParClasse, nbClasses, MPI_INT, m_rank + 1, 
  		m_rank, m_MPI_COMM_activProc );
    }
    else
    {
      MPI_Recv( tabNbPartRestantesParClasse, nbClasses, MPI_INT, m_rank - 1,
    	m_rank - 1, m_MPI_COMM_activProc, &status );      
	
      for (i=0,ipartProc=newPartProc.begin();i<nbClasses;++i,ipartProc++)
        ipartProc->second = tabNbPartRestantesParClasse[i];
    }
  }    	
  
  // Output and verify distribution per class and per process
  if ( m_rank == 0 )
    cout << endl << "Distribution " << endl << 
    	"Per process and per class" << endl;
  for (int m=0;m<m_nprocs;++m)
  {
    if ( m == m_rank && m_is_activ )
    {      
      cout << "Proc " << m_rank << endl;
      int ntot_ = 0 ;
      for (i=0,ipart=newPart.begin(),ipartProc=newPartProc.begin();
    	i<nbClasses;++i,ipart++,ipartProc++)
      {
        ntot_ += ipartProc->second;
	cout << "   Class " << i << " = " << ipartProc->second << endl;
      }
      cout << "   Total: " << ntot_ << " = " << npartproc << endl;

    }
    MPI_Barrier( m_MPI_COMM_activProc );
  }
  MPI_Barrier( m_MPI_COMM_activProc );

  if ( m_rank == 0 )
    cout << "Per class on all processes" << endl;
  for (i=0,ipart=newPart.begin(),ipartProc=newPartProc.begin();
    	i<nbClasses;++i,ipart++,ipartProc++)
  {
    int ntotc =  sum_INT_master( ipartProc->second ) ;
    if ( m_rank == 0 )
      cout << "   Class " << i << " " << ntotc << " = " << ipart->second 
      	<< endl;
  }
  if ( m_rank == 0 ) cout << endl;
} 




// ----------------------------------------------------------------------------
// Gathers all periodic clone particles on the master 
// process for post-processing purposes
list<Particle*>* GrainsMPIWrapper::GatherPeriodicClones_PostProcessing(
  	list<Particle*> const& periodicCloneParticles,
	vector<Particle*> const& referenceParticles ) const
{
  list<Particle*>* allClonesPeriodiques = NULL;
//   list<Particle*>::const_iterator il;
//   int i,j;
// 
//   // Taille des messages proportionnel au nombre de particles dans 
//   // ParticlesReference
//   // ---------------------------------------------------------------
//   int nb_part = int(particlesClonesPeriodiques.size());
//   for (il=particlesClonesPeriodiques.begin();
//   	il!=particlesClonesPeriodiques.end();il++)
//     if ((*il)->getTag()==2) nb_part--;
// 
//   // Taille des messages à passer
//   int *recvcounts = new int[m_nprocs];
//   for (i=0; i<m_nprocs; ++i) recvcounts[i]=0;
//   MPI_Gather( &nb_part, 1, MPI_INT, recvcounts, 1, MPI_INT, 
//   	m_rank_masterWorld, m_MPI_COMM_activProc ); 
// 
//   // Partition et taille du buffer de réception
//   int *displs = new int[m_nprocs];
//   displs[0] = 0; 
//   for (i=1; i<m_nprocs; ++i) displs[i] = displs[i-1]+recvcounts[i-1];  
//   int recvsize = displs[m_nprocs-1]+recvcounts[m_nprocs-1];
// 
// 
// 
//   // Communication des entiers: 
//   // Ordre par particle: 
//   // [numero de particle, classe]
//   // -----------------------------
//   int NB_INT_PART = 2;
//   int *numClass = new int[NB_INT_PART*nb_part];
//   for (il=particlesClonesPeriodiques.begin(),i=0;
//   	il!=particlesClonesPeriodiques.end();il++)
//   {
//     if ( (*il)->getTag() == 0 || (*il)->getTag() == 1 )
//     {    
//       numClass[i] = (*il)->getID();
//       numClass[i+1] = (*il)->getGeometricType();
//       i+=NB_INT_PART;
//     }
//   }            
// 
//   // Partition et taille du buffer de réception
//   int *recvcounts_INT = new int[m_nprocs];
//   for (i=0; i<m_nprocs; ++i) recvcounts_INT[i] = NB_INT_PART * recvcounts[i];
//   int *displs_INT = new int[m_nprocs];
//   displs_INT[0] = 0; 
//   for (i=1; i<m_nprocs; ++i) 
//     displs_INT[i] = displs_INT[i-1] + recvcounts_INT[i-1];  
//   int recvsize_INT = displs_INT[m_nprocs-1] + recvcounts_INT[m_nprocs-1];
//   int *recvbuf_INT = new int[recvsize_INT];
// 
//   // Communication des entiers: 2 integer par particle
//   MPI_Gatherv( numClass, NB_INT_PART * nb_part, MPI_INT, recvbuf_INT, 
//   	recvcounts_INT, displs_INT, MPI_INT, 
// 	m_rank_masterWorld, m_MPI_COMM_activProc );
// 
// 
// 
//   // Communication des doubles: cinématique & configuration
//   // Ordre par particle: 
//   // [position, vitesse translation, quaternion rotation, vitesse rotation]
//   // ----------------------------------------------------------------------
//   int NB_DOUBLE_PART = 26;  
//   double *features = new double[NB_DOUBLE_PART*nb_part];
//   for (il=particlesClonesPeriodiques.begin(),i=0;
//   	il!=particlesClonesPeriodiques.end();il++)
//   {
//     if ( (*il)->getTag() == 0 || (*il)->getTag() == 1 )
//     {
//       (*il)->copyTranslationalVelocity( features, i );
//       (*il)->copyQuaternionRotation( features, i+3 );    
//       (*il)->copyAngularVelocity( features, i+7 );
//       (*il)->copyTransform( features, i+10 ); 
//       i += NB_DOUBLE_PART;
//     }                  
//   }  
// 
//   // Partition et taille du buffer de réception
//   int *recvcounts_DOUBLE = new int[m_nprocs];
//   for (i=0; i<m_nprocs; ++i) 
//     recvcounts_DOUBLE[i] = NB_DOUBLE_PART * recvcounts[i];
//   int *displs_DOUBLE = new int[m_nprocs];
//   displs_DOUBLE[0] = 0; 
//   for (i=1; i<m_nprocs; ++i) 
//     displs_DOUBLE[i] = displs_DOUBLE[i-1] + recvcounts_DOUBLE[i-1];  
//   int recvsize_DOUBLE = displs_DOUBLE[m_nprocs-1] 
//   	+ recvcounts_DOUBLE[m_nprocs-1];
//   double *recvbuf_DOUBLE = new double[recvsize_DOUBLE];
//       
//   // Communication de la cinématique & configuration: 26 double par 
//   // particle
//   MPI_Gatherv( features, NB_DOUBLE_PART * nb_part, MPI_DOUBLE, recvbuf_DOUBLE, 
//   	recvcounts_DOUBLE, displs_DOUBLE, MPI_DOUBLE, 
// 	m_rank_masterWorld, m_MPI_COMM_activProc ); 
// 
// 
//   // Creation des Particles pour Post-processing
//   // --------------------------------------------
//   if ( m_rank == m_rank_masterWorld )
//   {
//     allClonesPeriodiques = new list<Particle*>;
// 
//     // Clones periodiques sur tous les proc
//     for (j=0;j<recvsize;++j)
//     {
//       // Creation de la particle
//       Particle *part_post=new Particle( recvbuf_INT[NB_INT_PART*j],
//       		referenceParticles[recvbuf_INT[NB_INT_PART*j+1]],
//       		recvbuf_DOUBLE[NB_DOUBLE_PART*j],
//       		recvbuf_DOUBLE[NB_DOUBLE_PART*j+1],
//     		recvbuf_DOUBLE[NB_DOUBLE_PART*j+2],
//       		recvbuf_DOUBLE[NB_DOUBLE_PART*j+3],
//       		recvbuf_DOUBLE[NB_DOUBLE_PART*j+4],
//     		recvbuf_DOUBLE[NB_DOUBLE_PART*j+5],		
//       		recvbuf_DOUBLE[NB_DOUBLE_PART*j+6],
//       		recvbuf_DOUBLE[NB_DOUBLE_PART*j+7],
//     		recvbuf_DOUBLE[NB_DOUBLE_PART*j+8],
//     		recvbuf_DOUBLE[NB_DOUBLE_PART*j+9],		
//       		&recvbuf_DOUBLE[NB_DOUBLE_PART*j+10],		
// 		COMPUTE,
// 		0 ); 		    
//       allClonesPeriodiques->push_back(part_post);    
//     }
//   } 
// 
//   delete [] numClass;
//   delete [] features;    
//   delete [] recvcounts;
//   delete [] recvcounts_INT;
//   delete [] recvcounts_DOUBLE;    
//   delete [] displs;
//   delete [] displs_INT;
//   delete [] displs_DOUBLE;      
//   delete [] recvbuf_INT;
//   delete [] recvbuf_DOUBLE;

  return ( allClonesPeriodiques );
} 




// ----------------------------------------------------------------------------
// Writes a string per process in a process-id ordered manner
void GrainsMPIWrapper::writeStringPerProcess( ostream& f, string const& out, 
    	bool creturn, string const& shift ) const
{
  string allout = shift + "Process 0";
  int dim = 0, recvsize = 0 ;
  int tag_CHAR = 111;
  MPI_Status status;
  
  if ( m_rank == 0 ) 
  {
    if ( creturn ) allout += "\n";
    else allout += ": ";
    allout += out + "\n"; 
    for ( int irank=1;irank<m_nprocs;irank++)
    {
      // Size of the message
      MPI_Probe( irank, tag_CHAR, m_MPI_COMM_activProc, &status );  
      MPI_Get_count( &status, MPI_CHAR, &recvsize );

      // Reception of the actual message	
      char *recvbuf_char = new char[recvsize];
      MPI_Recv( recvbuf_char, recvsize, MPI_CHAR, 
          irank, tag_CHAR, m_MPI_COMM_activProc, &status );

      if ( recvsize != 1 )
      { 
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
    MPI_Send( out.c_str(), dim, MPI_CHAR, 0, tag_CHAR, m_MPI_COMM_activProc );
  }

  // The master proc writes to the stream
  if ( m_rank == 0 ) f << allout << endl; 
}
