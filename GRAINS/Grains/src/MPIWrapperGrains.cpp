#include "MPIWrapperGrains.hh"
#include "App.H"
#include "Contact_BuilderFactory.hh"
#include "ObstaclePeriodique.hh"
#include "ParticulePeriodique.hh"
#include "Grains_Exec.hh"
#include "Cellule.H"
#include "Text_PostProcessingWriter.hh"
#include <assert.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>


string *MPIWrapperGrains::m_MPILogString = new string;
vector< vector<int> > MPIWrapperGrains::m_particuleHalozoneToNeighboringProcs;
vector<int> MPIWrapperGrains::m_GeoLocReciprocity;


//-----------------------------------------------------------------------------
// Constructeur par defaut
MPIWrapperGrains::MPIWrapperGrains() :
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
   m_voisins( NULL ),
   m_commgrainsMPI_3D( NULL ),
   m_nprocs_localComm( 0 ),
   m_rank_localComm( 0 ),
   m_master_localComm( NULL ),
   m_MPI_GROUP_Periodic( NULL ),
   m_MPI_COMM_Periodic( NULL ),
   m_isInCommPeriodic( false )
{}




//-----------------------------------------------------------------------------
// Constructeur avec argument la decomposition de domaine
MPIWrapperGrains::MPIWrapperGrains( int NX, int NY, int NZ,
  	int PERX, int PERY, int PERZ ) :
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
   m_voisins( NULL ),
   m_commgrainsMPI_3D( NULL ),
   m_nprocs_localComm( 0 ),
   m_rank_localComm( 0 ),
   m_master_localComm( NULL ),
   m_MPI_GROUP_Periodic( NULL ),
   m_MPI_COMM_Periodic( NULL ),
   m_isInCommPeriodic( false )
{
  // Nombre total de processus
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
      cout << "Total number of processus = " << m_nprocs_world << endl;
      cout << "Number of active processus = " << m_nprocs << endl;
      cout << "Number of sleeping processus = " << m_nprocs_world - m_nprocs
      	<< endl;
    }
  }

  // Groupe de processus actifs
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

  if ( m_is_activ )
  {
    // Rang du processus dans le groupe de processus actifs
    MPI_Comm_rank( m_MPI_COMM_activProc, &m_rank );

    // rang dans le communicateur m_MPI_COMM_activProc du processus ayant
    // le rang 0 (master) dans le communicateur MPI_COMM_WORLD
    int loc_m_rank_masterWorld = -1;
    if ( m_rank_world == 0 ) loc_m_rank_masterWorld = m_rank;
    MPI_Allreduce( &loc_m_rank_masterWorld, &m_rank_masterWorld, 1, MPI_INT,
    	MPI_MAX, m_MPI_COMM_activProc ) ;

    // Nombre de processus par direction
    m_dim = new int[3];
    m_dim[0] = NX, m_dim[1] = NY, m_dim[2] = NZ;

    // Periodicite du domaine global
    m_period = new int[3];
    m_period[0] = PERX, m_period[1] = PERY, m_period[2] = PERZ;
    if ( m_period[0] || m_period[1] || m_period[2] ) m_isMPIperiodic = true;

    // Creation des periodes liees a la topologie MPI dans le cas ou la
    // periodicite geometrique du domaine est traitee par MPI
    m_MPIperiodes.reserve( 27 );
    for (int i=0;i<27;++i) m_MPIperiodes.push_back( VecteurNul );

    // Correspondance entre la geolocalisation d'une particule dans la zone de
    // halo et les procs � qui les infos d'une particule doivent etre envoyees
    setParticuleHalozoneToNeighboringProcs();

    // R�-num�rotation des processus: non => rank == m_rank_world
    int reorganisation = 0;

    // Cr�ation d'une topologie cart�sienne
    m_commgrainsMPI_3D = new MPI_Comm;
    MPI_Cart_create( m_MPI_COMM_activProc, 3, m_dim, m_period, reorganisation,
   	m_commgrainsMPI_3D );

    // Acc�s au rang d'un proc donn� et � ses coordonn�es
    MPI_Comm_rank( *m_commgrainsMPI_3D, &m_rank );
    if ( m_nprocs == m_nprocs_world ) assert( m_rank == m_rank_world );
    m_coords = new int[3];
    MPI_Cart_coords( *m_commgrainsMPI_3D, m_rank, 3, m_coords );

    // Definition des voisins d'un processus dans la topologie cartesienne
    m_voisins = new Voisins( *m_commgrainsMPI_3D, m_coords, m_dim, m_period );

    // Timer
    SCT_insert_app( "Copie_Buffers" );
    SCT_insert_app( "MPIComm" );
    SCT_insert_app( "UpdateCreateClones" );
  }
  else m_rank = -1;

}




//-----------------------------------------------------------------------------
// Destructeur
MPIWrapperGrains::~MPIWrapperGrains()
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
    delete m_voisins;
    vector<MPI_Group*>::iterator ivg;
    for (ivg=m_groupVoisins.begin();ivg!=m_groupVoisins.end();ivg++)
    {
      MPI_Group_free( *ivg );
      delete *ivg;
    }
    m_groupVoisins.clear();
    if (!m_isInCommVoisins.empty())
      for (int j=0;j<m_nprocs;j++)
      {
        if ( m_isInCommVoisins[j] ) MPI_Comm_free( m_commVoisins[j] );
        delete m_commVoisins[j];
      }
    m_commVoisins.clear();
    m_isInCommVoisins.clear();
    if ( m_master_localComm ) delete [] m_master_localComm;
    if ( m_MPILogString ) delete m_MPILogString;
    vector< vector<int> >::iterator ivv;
    for (ivv=m_particuleHalozoneToNeighboringProcs.begin();
  	ivv!=m_particuleHalozoneToNeighboringProcs.end();ivv++)
      ivv->clear();
    m_particuleHalozoneToNeighboringProcs.clear();

    if ( m_MPI_GROUP_Periodic )
    {
      MPI_Group_free( m_MPI_GROUP_Periodic );
      delete m_MPI_GROUP_Periodic;
    }
    if ( m_MPI_COMM_Periodic )
    {
      if ( m_isInCommPeriodic ) MPI_Comm_free( m_MPI_COMM_Periodic );
      delete m_MPI_COMM_Periodic;
    }

    m_MPIperiodes.clear();
  }
}




// ----------------------------------------------------------------------------
// Definition des communicateurs locaux lies aux voisins
void MPIWrapperGrains::setCommLocal()
{
  int i,j;

  // Taille des messages proportionnel au nombre de voisins
  // ------------------------------------------------------
  int nbv = m_voisins->nbVoisins();

  // Taille des messages � passer
  int *recvcounts = new int[m_nprocs];
  for (i=0; i<m_nprocs; ++i) recvcounts[i]=0;
  MPI_Allgather( &nbv, 1, MPI_INT, recvcounts, 1, MPI_INT,
  	m_MPI_COMM_activProc );

  // Partition et taille du buffer de r�ception
  int *displs = new int[m_nprocs];
  displs[0] = 0;
  for (i=1; i<m_nprocs; ++i) displs[i] = displs[i-1]+recvcounts[i-1];
  int recvsize = displs[m_nprocs-1] + recvcounts[m_nprocs-1];


  // Communication des rangs des voisins
  // -----------------------------------
  int *rangVoisins = m_voisins->rangVoisins();
  int *recvbuf_rang = new int[recvsize];
  MPI_Allgatherv( rangVoisins, nbv, MPI_INT, recvbuf_rang,
  	recvcounts, displs, MPI_INT, m_MPI_COMM_activProc );


  // Cr�ation des communicateurs lies aux voisins
  // --------------------------------------------
  m_groupVoisins.reserve( m_nprocs );
  MPI_Group* empty_MPI_Group = NULL;
  for (j=0;j<m_nprocs;j++) m_groupVoisins.push_back(empty_MPI_Group);
  m_commVoisins.reserve( m_nprocs );
  MPI_Comm* empty_MPI_Comm = NULL;
  for (j=0;j<m_nprocs;j++) m_commVoisins.push_back(empty_MPI_Comm);
  m_isInCommVoisins.reserve( m_nprocs );
  for (j=0;j<m_nprocs;j++) m_isInCommVoisins.push_back(false);
  MPI_Group activ_group;
  MPI_Comm_group( m_MPI_COMM_activProc, &activ_group );
  for (j=0;j<m_nprocs;j++)
  {
    for (i=displs[j];i<displs[j]+recvcounts[j];++i)
      if ( recvbuf_rang[i] == m_rank ) m_isInCommVoisins[j] = true;
    m_groupVoisins[j] = new MPI_Group;
    m_commVoisins[j] = new MPI_Comm;
    MPI_Group_incl( activ_group, recvcounts[j], &recvbuf_rang[displs[j]],
    	m_groupVoisins[j] );
    MPI_Comm_create( m_MPI_COMM_activProc, *m_groupVoisins[j],
    	m_commVoisins[j] );
  }
  MPI_Group_free(&activ_group);


  // Rang & taille du communicateur local
  // ------------------------------------
  MPI_Comm_rank( *m_commVoisins[m_rank], &m_rank_localComm );
  MPI_Comm_size( *m_commVoisins[m_rank], &m_nprocs_localComm );
  m_master_localComm = new int[m_nprocs];
  MPI_Allgather( &m_rank_localComm, 1, MPI_INT, m_master_localComm, 1, MPI_INT,
  	m_MPI_COMM_activProc );


  delete [] recvcounts;
  delete [] displs;
  delete [] rangVoisins;
  delete [] recvbuf_rang;


  // Definition de la geolocalisation du master dans les communicateurs
  // lies aux voisins auquel appartient le proc
  // ------------------------------------------------------------------
  setMasterGeoLocInLocalComm();


  // Sortie ecran
  // ------------
  if ( 0 == m_rank ) cout << "Definition des communicateurs locaux" << endl;
  for (i=0;i<m_nprocs;++i)
  {
    if ( i == m_rank )
    {
      cout << "Processeur = " << m_rank << endl;
      cout << "   Rang dans Comm local = " << m_rank_localComm << endl;
      cout << "   Coordinates dans Comm local = " << m_coords[0] << " " <<
      	m_coords[1] << " " << m_coords[2] << endl;
      cout << "   Nb procs dans Comm local = " << m_nprocs_localComm << endl;
      cout << "   Rang du master dans Comms locaux : ";
      for (j=0;j<m_nprocs;++j) cout << m_master_localComm[j] << " ";
      cout << endl;
      cout << "   Geolocalisation du master dans Comms locaux : ";
      for (j=0;j<m_nprocs;++j)
        if ( m_masterGeoLoc[j] != -1 ) cout <<
		Cellule::getMPIGeoLocalisationName(m_masterGeoLoc[j]) << " ";
	else cout << m_masterGeoLoc[j] << " ";
      cout << endl;
      cout << "   Rang dans Comm Activ du master de Comm World = " <<
      	m_rank_masterWorld << endl;
    }
    MPI_Barrier( m_MPI_COMM_activProc );
  }
  if ( 0 == m_rank ) cout << endl;
}




//-----------------------------------------------------------------------------
// Definition de la geolocalisation du master dans les communicateurs
// lies aux voisins auquel appartient le proc
void MPIWrapperGrains::setMasterGeoLocInLocalComm()
{
  int i,j,k,globalRank,ii;

  m_masterGeoLoc.reserve(m_nprocs);
  for (i=0;i<m_nprocs;++i) m_masterGeoLoc.push_back(-1);

  // Geolocalisation du master pour ses voisins sur le proc
  // ------------------------------------------------------
  int *geoLocMasterOnProc = new int[m_nprocs];
  for (i=0;i<m_nprocs;++i) geoLocMasterOnProc[i] = -1;
  for (i=-1;i<2;++i)
    for (j=-1;j<2;++j)
      for (k=-1;k<2;++k)
        if ( i || j || k )
	{
	  globalRank = m_voisins->rank( i, j, k );
	  if ( globalRank != -1 )
	    geoLocMasterOnProc[globalRank] =
	    	getMPIGeoLocalisation( -i, -j, -k );
	}

  // Communication des geolocalisations des masters
  // ----------------------------------------------
  int *recvcounts = new int[m_nprocs];
  for (i=0; i<m_nprocs; ++i) recvcounts[i] = m_nprocs;
  int *displs = new int[m_nprocs];
  displs[0] = 0;
  for (i=1; i<m_nprocs; ++i) displs[i] = displs[i-1] + recvcounts[i-1];
  int recvsize = displs[m_nprocs-1] + recvcounts[m_nprocs-1];

  int *allGeoLocMasterOnProc = new int[recvsize];
  MPI_Allgatherv( geoLocMasterOnProc, m_nprocs, MPI_INT, allGeoLocMasterOnProc,
  	recvcounts, displs, MPI_INT, m_MPI_COMM_activProc );

  // Stockage dans masterGeoLoc
  // --------------------------
  for (i=0;i<m_nprocs*m_nprocs;++i)
  {
    ii = int(i/m_nprocs);
    j = i - ii * m_nprocs;
    if ( j == m_rank )
      if ( allGeoLocMasterOnProc[i] != -1 )
        m_masterGeoLoc[ii] = allGeoLocMasterOnProc[i];
  }
}




//-----------------------------------------------------------------------------
// Definition de particuleHalozoneToNeighboringProcs
void MPIWrapperGrains::setParticuleHalozoneToNeighboringProcs()
{
  m_particuleHalozoneToNeighboringProcs.reserve(26);
  vector<int> emptyVECINT;
  for (int i=0;i<26;++i)
    m_particuleHalozoneToNeighboringProcs.push_back(emptyVECINT);
  vector<int> *work;

  // NORTH => NORTH
  work = new vector<int>(1,0);
  (*work)[0] = MPIGEO_NORTH;
  m_particuleHalozoneToNeighboringProcs[MPIGEO_NORTH] = *work;
  work->clear();
  delete work;

  // NORTH_EAST => NORTH, EAST, NORTH_EAST
  work = new vector<int>(3,0);
  (*work)[0] = MPIGEO_NORTH;
  (*work)[1] = MPIGEO_EAST;
  (*work)[2] = MPIGEO_NORTH_EAST;
  m_particuleHalozoneToNeighboringProcs[MPIGEO_NORTH_EAST] = *work;
  work->clear();
  delete work;

  // NORTH_WEST => NORTH, WEST, NORTH_WEST
  work = new vector<int>(3,0);
  (*work)[0] = MPIGEO_NORTH;
  (*work)[1] = MPIGEO_WEST;
  (*work)[2] = MPIGEO_NORTH_WEST;
  m_particuleHalozoneToNeighboringProcs[MPIGEO_NORTH_WEST] = *work;
  work->clear();
  delete work;

  // NORTH_TOP => NORTH, TOP, NORTH_TOP
  work = new vector<int>(3,0);
  (*work)[0] = MPIGEO_NORTH;
  (*work)[1] = MPIGEO_TOP;
  (*work)[2] = MPIGEO_NORTH_TOP;
  m_particuleHalozoneToNeighboringProcs[MPIGEO_NORTH_TOP] = *work;
  work->clear();
  delete work;

  // NORTH_BOTTOM => NORTH, BOTTOM, NORTH_BOTTOM
  work = new vector<int>(3,0);
  (*work)[0] = MPIGEO_NORTH;
  (*work)[1] = MPIGEO_BOTTOM;
  (*work)[2] = MPIGEO_NORTH_BOTTOM;
  m_particuleHalozoneToNeighboringProcs[MPIGEO_NORTH_BOTTOM] = *work;
  work->clear();
  delete work;

  // NORTH_EAST_TOP => NORTH, EAST, TOP, EAST_TOP, NORTH_EAST, NORTH_TOP,
  // NORTH_EAST_TOP
  work = new vector<int>(7,0);
  (*work)[0] = MPIGEO_NORTH;
  (*work)[1] = MPIGEO_EAST;
  (*work)[2] = MPIGEO_TOP;
  (*work)[3] = MPIGEO_EAST_TOP;
  (*work)[4] = MPIGEO_NORTH_EAST;
  (*work)[5] = MPIGEO_NORTH_TOP;
  (*work)[6] = MPIGEO_NORTH_EAST_TOP;
  m_particuleHalozoneToNeighboringProcs[MPIGEO_NORTH_EAST_TOP] = *work;
  work->clear();
  delete work;

  // NORTH_EAST_BOTTOM => NORTH, EAST, BOTTOM, EAST_BOTTOM, NORTH_EAST,
  // NORTH_BOTTOM, NORTH_EAST_BOTTOM
  work = new vector<int>(7,0);
  (*work)[0] = MPIGEO_NORTH;
  (*work)[1] = MPIGEO_EAST;
  (*work)[2] = MPIGEO_BOTTOM;
  (*work)[3] = MPIGEO_EAST_BOTTOM;
  (*work)[4] = MPIGEO_NORTH_EAST;
  (*work)[5] = MPIGEO_NORTH_BOTTOM;
  (*work)[6] = MPIGEO_NORTH_EAST_BOTTOM;
  m_particuleHalozoneToNeighboringProcs[MPIGEO_NORTH_EAST_BOTTOM] = *work;
  work->clear();
  delete work;

  // NORTH_WEST_TOP => NORTH, WEST, TOP, WEST_TOP, NORTH_WEST, NORTH_TOP,
  // NORTH_WEST_TOP
  work = new vector<int>(7,0);
  (*work)[0] = MPIGEO_NORTH;
  (*work)[1] = MPIGEO_WEST;
  (*work)[2] = MPIGEO_TOP;
  (*work)[3] = MPIGEO_WEST_TOP;
  (*work)[4] = MPIGEO_NORTH_WEST;
  (*work)[5] = MPIGEO_NORTH_TOP;
  (*work)[6] = MPIGEO_NORTH_WEST_TOP;
  m_particuleHalozoneToNeighboringProcs[MPIGEO_NORTH_WEST_TOP] = *work;
  work->clear();
  delete work;

  // NORTH_WEST_BOTTOM => NORTH, WEST, BOTTOM, WEST_BOTTOM, NORTH_WEST,
  // NORTH_BOTTOM, NORTH_WEST_BOTTOM
  work = new vector<int>(7,0);
  (*work)[0] = MPIGEO_NORTH;
  (*work)[1] = MPIGEO_WEST;
  (*work)[2] = MPIGEO_BOTTOM;
  (*work)[3] = MPIGEO_WEST_BOTTOM;
  (*work)[4] = MPIGEO_NORTH_WEST;
  (*work)[5] = MPIGEO_NORTH_BOTTOM;
  (*work)[6] = MPIGEO_NORTH_WEST_BOTTOM;
  m_particuleHalozoneToNeighboringProcs[MPIGEO_NORTH_WEST_BOTTOM] = *work;
  work->clear();
  delete work;

  // SOUTH => SOUTH
  work = new vector<int>(1,0);
  (*work)[0] = MPIGEO_SOUTH;
  m_particuleHalozoneToNeighboringProcs[MPIGEO_SOUTH] = *work;
  work->clear();
  delete work;

  // SOUTH_EAST => SOUTH, EAST, SOUTH_EAST
  work = new vector<int>(3,0);
  (*work)[0] = MPIGEO_SOUTH;
  (*work)[1] = MPIGEO_EAST;
  (*work)[2] = MPIGEO_SOUTH_EAST;
  m_particuleHalozoneToNeighboringProcs[MPIGEO_SOUTH_EAST] = *work;
  work->clear();
  delete work;

  // SOUTH_WEST => SOUTH, WEST, SOUTH_WEST
  work = new vector<int>(3,0);
  (*work)[0] = MPIGEO_SOUTH;
  (*work)[1] = MPIGEO_WEST;
  (*work)[2] = MPIGEO_SOUTH_WEST;
  m_particuleHalozoneToNeighboringProcs[MPIGEO_SOUTH_WEST] = *work;
  work->clear();
  delete work;

  // SOUTH_TOP => SOUTH, TOP, SOUTH_TOP
  work = new vector<int>(3,0);
  (*work)[0] = MPIGEO_SOUTH;
  (*work)[1] = MPIGEO_TOP;
  (*work)[2] = MPIGEO_SOUTH_TOP;
  m_particuleHalozoneToNeighboringProcs[MPIGEO_SOUTH_TOP] = *work;
  work->clear();
  delete work;

  // SOUTH_BOTTOM => SOUTH, BOTTOM, SOUTH_BOTTOM
  work = new vector<int>(3,0);
  (*work)[0] = MPIGEO_SOUTH;
  (*work)[1] = MPIGEO_BOTTOM;
  (*work)[2] = MPIGEO_SOUTH_BOTTOM;
  m_particuleHalozoneToNeighboringProcs[MPIGEO_SOUTH_BOTTOM] = *work;
  work->clear();
  delete work;

  // SOUTH_EAST_TOP => SOUTH, EAST, TOP, EAST_TOP, SOUTH_EAST, SOUTH_TOP,
  // SOUTH_EAST_TOP
  work = new vector<int>(7,0);
  (*work)[0] = MPIGEO_SOUTH;
  (*work)[1] = MPIGEO_EAST;
  (*work)[2] = MPIGEO_TOP;
  (*work)[3] = MPIGEO_EAST_TOP;
  (*work)[4] = MPIGEO_SOUTH_EAST;
  (*work)[5] = MPIGEO_SOUTH_TOP;
  (*work)[6] = MPIGEO_SOUTH_EAST_TOP;
  m_particuleHalozoneToNeighboringProcs[MPIGEO_SOUTH_EAST_TOP] = *work;
  work->clear();
  delete work;

  // SOUTH_EAST_BOTTOM => SOUTH, EAST, BOTTOM, EAST_BOTTOM, SOUTH_EAST,
  // SOUTH_BOTTOM, SOUTH_EAST_BOTTOM
  work = new vector<int>(7,0);
  (*work)[0] = MPIGEO_SOUTH;
  (*work)[1] = MPIGEO_EAST;
  (*work)[2] = MPIGEO_BOTTOM;
  (*work)[3] = MPIGEO_EAST_BOTTOM;
  (*work)[4] = MPIGEO_SOUTH_EAST;
  (*work)[5] = MPIGEO_SOUTH_BOTTOM;
  (*work)[6] = MPIGEO_SOUTH_EAST_BOTTOM;
  m_particuleHalozoneToNeighboringProcs[MPIGEO_SOUTH_EAST_BOTTOM] = *work;
  work->clear();
  delete work;

  // SOUTH_WEST_TOP => SOUTH, WEST, TOP, WEST_TOP, SOUTH_WEST, SOUTH_TOP,
  // SOUTH_WEST_TOP
  work = new vector<int>(7,0);
  (*work)[0] = MPIGEO_SOUTH;
  (*work)[1] = MPIGEO_WEST;
  (*work)[2] = MPIGEO_TOP;
  (*work)[3] = MPIGEO_WEST_TOP;
  (*work)[4] = MPIGEO_SOUTH_WEST;
  (*work)[5] = MPIGEO_SOUTH_TOP;
  (*work)[6] = MPIGEO_SOUTH_WEST_TOP;
  m_particuleHalozoneToNeighboringProcs[MPIGEO_SOUTH_WEST_TOP] = *work;
  work->clear();
  delete work;

  // SOUTH_WEST_BOTTOM => SOUTH, WEST, BOTTOM, WEST_BOTTOM, SOUTH_WEST,
  // SOUTH_BOTTOM, SOUTH_WEST_BOTTOM
  work = new vector<int>(7,0);
  (*work)[0] = MPIGEO_SOUTH;
  (*work)[1] = MPIGEO_WEST;
  (*work)[2] = MPIGEO_BOTTOM;
  (*work)[3] = MPIGEO_WEST_BOTTOM;
  (*work)[4] = MPIGEO_SOUTH_WEST;
  (*work)[5] = MPIGEO_SOUTH_BOTTOM;
  (*work)[6] = MPIGEO_SOUTH_WEST_BOTTOM;
  m_particuleHalozoneToNeighboringProcs[MPIGEO_SOUTH_WEST_BOTTOM] = *work;
  work->clear();
  delete work;

  // EAST => EAST
  work = new vector<int>(1,0);
  (*work)[0] = MPIGEO_EAST;
  m_particuleHalozoneToNeighboringProcs[MPIGEO_EAST] = *work;
  work->clear();
  delete work;

  // WEST => WEST
  work = new vector<int>(1,0);
  (*work)[0] = MPIGEO_WEST;
  m_particuleHalozoneToNeighboringProcs[MPIGEO_WEST] = *work;
  work->clear();
  delete work;

  // EAST_TOP => EAST, TOP, EAST_TOP
  work = new vector<int>(3,0);
  (*work)[0] = MPIGEO_EAST;
  (*work)[1] = MPIGEO_TOP;
  (*work)[2] = MPIGEO_EAST_TOP;
  m_particuleHalozoneToNeighboringProcs[MPIGEO_EAST_TOP] = *work;
  work->clear();
  delete work;

  // EAST_BOTTOM => EAST, BOTTOM, EAST_BOTTOM
  work = new vector<int>(3,0);
  (*work)[0] = MPIGEO_EAST;
  (*work)[1] = MPIGEO_BOTTOM;
  (*work)[2] = MPIGEO_EAST_BOTTOM;
  m_particuleHalozoneToNeighboringProcs[MPIGEO_EAST_BOTTOM] = *work;
  work->clear();
  delete work;

  // WEST_TOP => WEST, TOP, WEST_TOP
  work = new vector<int>(3,0);
  (*work)[0] = MPIGEO_WEST;
  (*work)[1] = MPIGEO_TOP;
  (*work)[2] = MPIGEO_WEST_TOP;
  m_particuleHalozoneToNeighboringProcs[MPIGEO_WEST_TOP] = *work;
  work->clear();
  delete work;

  // WEST_BOTTOM => WEST, BOTTOM, WEST_BOTTOM
  work = new vector<int>(3,0);
  (*work)[0] = MPIGEO_WEST;
  (*work)[1] = MPIGEO_BOTTOM;
  (*work)[2] = MPIGEO_WEST_BOTTOM;
  m_particuleHalozoneToNeighboringProcs[MPIGEO_WEST_BOTTOM] = *work;
  work->clear();
  delete work;

  // TOP => TOP
  work = new vector<int>(1,0);
  (*work)[0] = MPIGEO_TOP;
  m_particuleHalozoneToNeighboringProcs[MPIGEO_TOP] = *work;
  work->clear();
  delete work;

  // BOTTOM => BOTTOM
  work = new vector<int>(1,0);
  (*work)[0] = MPIGEO_BOTTOM;
  m_particuleHalozoneToNeighboringProcs[MPIGEO_BOTTOM] = *work;
  work->clear();
  delete work;


  // Geolocalisation reciproque
  m_GeoLocReciprocity.reserve(26);
  for (int i=0;i<26;++i) m_GeoLocReciprocity.push_back( 0 );
  m_GeoLocReciprocity[MPIGEO_NORTH] = MPIGEO_SOUTH ;
  m_GeoLocReciprocity[MPIGEO_NORTH_EAST] = MPIGEO_SOUTH_WEST ;
  m_GeoLocReciprocity[MPIGEO_NORTH_WEST] = MPIGEO_SOUTH_EAST ;
  m_GeoLocReciprocity[MPIGEO_NORTH_TOP] = MPIGEO_SOUTH_BOTTOM ;
  m_GeoLocReciprocity[MPIGEO_NORTH_BOTTOM] = MPIGEO_SOUTH_TOP ;
  m_GeoLocReciprocity[MPIGEO_NORTH_EAST_TOP] = MPIGEO_SOUTH_WEST_BOTTOM ;
  m_GeoLocReciprocity[MPIGEO_NORTH_EAST_BOTTOM] = MPIGEO_SOUTH_WEST_TOP ;
  m_GeoLocReciprocity[MPIGEO_NORTH_WEST_TOP] = MPIGEO_SOUTH_EAST_BOTTOM ;
  m_GeoLocReciprocity[MPIGEO_NORTH_WEST_BOTTOM] = MPIGEO_SOUTH_EAST_TOP ;
  m_GeoLocReciprocity[MPIGEO_SOUTH] = MPIGEO_NORTH ;
  m_GeoLocReciprocity[MPIGEO_SOUTH_EAST] = MPIGEO_NORTH_WEST ;
  m_GeoLocReciprocity[MPIGEO_SOUTH_WEST] = MPIGEO_NORTH_EAST ;
  m_GeoLocReciprocity[MPIGEO_SOUTH_TOP] = MPIGEO_NORTH_BOTTOM ;
  m_GeoLocReciprocity[MPIGEO_SOUTH_BOTTOM] = MPIGEO_NORTH_TOP ;
  m_GeoLocReciprocity[MPIGEO_SOUTH_EAST_TOP] = MPIGEO_NORTH_WEST_BOTTOM ;
  m_GeoLocReciprocity[MPIGEO_SOUTH_EAST_BOTTOM] = MPIGEO_NORTH_WEST_TOP ;
  m_GeoLocReciprocity[MPIGEO_SOUTH_WEST_TOP] = MPIGEO_NORTH_EAST_BOTTOM ;
  m_GeoLocReciprocity[MPIGEO_SOUTH_WEST_BOTTOM] = MPIGEO_NORTH_EAST_TOP ;
  m_GeoLocReciprocity[MPIGEO_EAST] = MPIGEO_WEST ;
  m_GeoLocReciprocity[MPIGEO_WEST] = MPIGEO_EAST ;
  m_GeoLocReciprocity[MPIGEO_EAST_TOP] = MPIGEO_WEST_BOTTOM ;
  m_GeoLocReciprocity[MPIGEO_EAST_BOTTOM] = MPIGEO_WEST_TOP ;
  m_GeoLocReciprocity[MPIGEO_WEST_TOP] = MPIGEO_EAST_BOTTOM ;
  m_GeoLocReciprocity[MPIGEO_WEST_BOTTOM] = MPIGEO_EAST_TOP ;
  m_GeoLocReciprocity[MPIGEO_TOP] = MPIGEO_BOTTOM ;
  m_GeoLocReciprocity[MPIGEO_BOTTOM] = MPIGEO_TOP ;
}




//-----------------------------------------------------------------------------
// Renvoie la MPIGeoLocalisation en fonction de la position relative
MPIGeoLocalisation MPIWrapperGrains::getMPIGeoLocalisation(int i,int j,int k)
{
  MPIGeoLocalisation geoLoc = MPIGEO_NONE;
  switch( i )
  {
    case -1:
      switch( j )
      {
        case -1:
	  switch( k )
	  {
	    case -1:
	      geoLoc = MPIGEO_SOUTH_WEST_BOTTOM;
	      break;
	    case 0:
	      geoLoc = MPIGEO_SOUTH_WEST;
	      break;
	    case 1:
	      geoLoc = MPIGEO_SOUTH_WEST_TOP;
	      break;
	  }
	  break;

	case 0:
	  switch( k )
	  {
	    case -1:
	      geoLoc = MPIGEO_WEST_BOTTOM;
	      break;
	    case 0:
	      geoLoc = MPIGEO_WEST;
	      break;
	    case 1:
	      geoLoc = MPIGEO_WEST_TOP;
	      break;
	  }
	  break;

	case 1:
	  switch( k )
	  {
	    case -1:
	      geoLoc = MPIGEO_NORTH_WEST_BOTTOM;
	      break;
	    case 0:
	      geoLoc = MPIGEO_NORTH_WEST;
	      break;
	    case 1:
	      geoLoc = MPIGEO_NORTH_WEST_TOP;
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
	      geoLoc = MPIGEO_SOUTH_BOTTOM;
	      break;
	    case 0:
	      geoLoc = MPIGEO_SOUTH;
	      break;
	    case 1:
	      geoLoc = MPIGEO_SOUTH_TOP;
	      break;
	  }
	  break;

	case 0:
	  switch( k )
	  {
	    case -1:
	      geoLoc = MPIGEO_BOTTOM;
	      break;
	    case 0:
	      geoLoc = MPIGEO_NONE;
	      break;
	    case 1:
	      geoLoc = MPIGEO_TOP;
	      break;
	  }
	  break;

	case 1:
	  switch( k )
	  {
	    case -1:
	      geoLoc = MPIGEO_NORTH_BOTTOM;
	      break;
	    case 0:
	      geoLoc = MPIGEO_NORTH;
	      break;
	    case 1:
	      geoLoc = MPIGEO_NORTH_TOP;
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
	      geoLoc = MPIGEO_SOUTH_EAST_BOTTOM;
	      break;
	    case 0:
	      geoLoc = MPIGEO_SOUTH_EAST;
	      break;
	    case 1:
	      geoLoc = MPIGEO_SOUTH_EAST_TOP;
	      break;
	  }
	  break;

	case 0:
	  switch( k )
	  {
	    case -1:
	      geoLoc = MPIGEO_EAST_BOTTOM;
	      break;
	    case 0:
	      geoLoc = MPIGEO_EAST;
	      break;
	    case 1:
	      geoLoc = MPIGEO_EAST_TOP;
	      break;
	  }
	  break;

	case 1:
	  switch( k )
	  {
	    case -1:
	      geoLoc = MPIGEO_NORTH_EAST_BOTTOM;
	      break;
	    case 0:
	      geoLoc = MPIGEO_NORTH_EAST;
	      break;
	    case 1:
	      geoLoc = MPIGEO_NORTH_EAST_TOP;
	      break;
	  }
	  break;
      }
      break;
  }

  return geoLoc;
}




//-----------------------------------------------------------------------------
// Definition des vecteurs de periodicite MPI dans le communicateur
// standard commgrainsMPI_3D
void MPIWrapperGrains::setMPIperiodicVectors( const Scalar& lx,
	const Scalar& ly, const Scalar& lz )
{
  // West
  if ( m_voisins->rank( -1, 0, 0 ) != -1 )
  {
    if ( m_period[0] && m_coords[0] == 0 )
    {
      m_MPIperiodes[MPIGEO_WEST][X] = lx ;
      m_MPIperiodes[MPIGEO_WEST][Y] = 0. ;
      m_MPIperiodes[MPIGEO_WEST][Z] = 0. ;
    }
  }

  // East
  if ( m_voisins->rank( 1, 0, 0 ) != -1 )
  {
    if ( m_period[0] && m_coords[0] == m_dim[0] - 1 )
    {
      m_MPIperiodes[MPIGEO_EAST][X] = -lx ;
      m_MPIperiodes[MPIGEO_EAST][Y] = 0. ;
      m_MPIperiodes[MPIGEO_EAST][Z] = 0. ;
    }
  }

  // South
  if ( m_voisins->rank( 0, -1, 0 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == 0 )
    {
      m_MPIperiodes[MPIGEO_SOUTH][X] = 0. ;
      m_MPIperiodes[MPIGEO_SOUTH][Y] = ly ;
      m_MPIperiodes[MPIGEO_SOUTH][Z] = 0. ;
    }
  }

  // North
  if ( m_voisins->rank( 0, 1, 0 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == m_dim[1] - 1 )
    {
      m_MPIperiodes[MPIGEO_NORTH][X] = 0. ;
      m_MPIperiodes[MPIGEO_NORTH][Y] = -ly ;
      m_MPIperiodes[MPIGEO_NORTH][Z] = 0. ;
    }
  }

  // Bottom
  if ( m_voisins->rank( 0, 0, -1 ) != -1 )
  {
    if ( m_period[2] && m_coords[2] == 0 )
    {
      m_MPIperiodes[MPIGEO_BOTTOM][X] = 0. ;
      m_MPIperiodes[MPIGEO_BOTTOM][Y] = 0. ;
      m_MPIperiodes[MPIGEO_BOTTOM][Z] = lz ;
    }
  }

  // Top
  if ( m_voisins->rank( 0, 0, 1 ) != -1 )
  {
    if ( m_period[2] && m_coords[2] == m_dim[2] - 1 )
    {
      m_MPIperiodes[MPIGEO_TOP][X] = 0. ;
      m_MPIperiodes[MPIGEO_TOP][Y] = 0. ;
      m_MPIperiodes[MPIGEO_TOP][Z] = -lz ;
    }
  }


  // South West
  if ( m_voisins->rank( -1, -1, 0 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == 0 )
    {
      m_MPIperiodes[MPIGEO_SOUTH_WEST][X] = 0. ;
      m_MPIperiodes[MPIGEO_SOUTH_WEST][Y] = ly ;
      m_MPIperiodes[MPIGEO_SOUTH_WEST][Z] = 0. ;
    }

    if ( m_period[0] && m_coords[0] == 0 )
      m_MPIperiodes[MPIGEO_SOUTH_WEST][X] += lx ;
  }

  // South East
  if ( m_voisins->rank( 1, -1, 0 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == 0 )
    {
      m_MPIperiodes[MPIGEO_SOUTH_EAST][X] = 0. ;
      m_MPIperiodes[MPIGEO_SOUTH_EAST][Y] = ly ;
      m_MPIperiodes[MPIGEO_SOUTH_EAST][Z] = 0. ;
    }

    if ( m_period[0] && m_coords[0] == m_dim[0] - 1 )
      m_MPIperiodes[MPIGEO_SOUTH_EAST][X] += -lx ;
  }

  // South Bottom
  if ( m_voisins->rank( 0, -1, -1 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == 0 )
    {
      m_MPIperiodes[MPIGEO_SOUTH_BOTTOM][X] = 0. ;
      m_MPIperiodes[MPIGEO_SOUTH_BOTTOM][Y] = ly ;
      m_MPIperiodes[MPIGEO_SOUTH_BOTTOM][Z] = 0. ;
    }

    if ( m_period[2] && m_coords[2] == 0 )
      m_MPIperiodes[MPIGEO_SOUTH_BOTTOM][Z] += lz ;
  }

  // South Top
  if ( m_voisins->rank( 0, -1, 1 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == 0 )
    {
      m_MPIperiodes[MPIGEO_SOUTH_TOP][X] = 0. ;
      m_MPIperiodes[MPIGEO_SOUTH_TOP][Y] = ly ;
      m_MPIperiodes[MPIGEO_SOUTH_TOP][Z] = 0. ;
    }

    if ( m_period[2] && m_coords[2] == m_dim[2] - 1 )
      m_MPIperiodes[MPIGEO_SOUTH_TOP][Z] += -lz ;
  }

  // North West
  if ( m_voisins->rank( -1, 1, 0 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == m_dim[1] - 1 )
    {
      m_MPIperiodes[MPIGEO_NORTH_WEST][X] = 0. ;
      m_MPIperiodes[MPIGEO_NORTH_WEST][Y] = -ly ;
      m_MPIperiodes[MPIGEO_NORTH_WEST][Z] = 0. ;
    }

    if ( m_period[0] && m_coords[0] == 0 )
      m_MPIperiodes[MPIGEO_NORTH_WEST][X] = +lx ;
  }

  // North East
  if ( m_voisins->rank( 1, 1, 0 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == m_dim[1] - 1 )
    {
      m_MPIperiodes[MPIGEO_NORTH_EAST][X] = 0. ;
      m_MPIperiodes[MPIGEO_NORTH_EAST][Y] = -ly ;
      m_MPIperiodes[MPIGEO_NORTH_EAST][Z] = 0. ;
    }

    if ( m_period[0] && m_coords[0] == m_dim[0] - 1 )
      m_MPIperiodes[MPIGEO_NORTH_EAST][X] += -lx ;
  }

  // North Bottom
  if ( m_voisins->rank( 0, 1, -1 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == m_dim[1] - 1 )
    {
      m_MPIperiodes[MPIGEO_NORTH_BOTTOM][X] = 0. ;
      m_MPIperiodes[MPIGEO_NORTH_BOTTOM][Y] = -ly ;
      m_MPIperiodes[MPIGEO_NORTH_BOTTOM][Z] = 0. ;
    }

    if ( m_period[2] && m_coords[2] == 0 )
      m_MPIperiodes[MPIGEO_NORTH_BOTTOM][Z] += lz ;
  }

  // North Top
  if ( m_voisins->rank( 0, 1, 1 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == m_dim[1] - 1 )
    {
      m_MPIperiodes[MPIGEO_NORTH_TOP][X] = 0. ;
      m_MPIperiodes[MPIGEO_NORTH_TOP][Y] = -ly ;
      m_MPIperiodes[MPIGEO_NORTH_TOP][Z] = 0. ;
    }

    if ( m_period[2] && m_coords[2] == m_dim[2] - 1 )
      m_MPIperiodes[MPIGEO_NORTH_TOP][Z] += -lz ;
  }

  // West Bottom
  if ( m_voisins->rank( -1, 0, -1 ) != -1 )
  {
    if ( m_period[0] && m_coords[0] == 0 )
    {
      m_MPIperiodes[MPIGEO_WEST_BOTTOM][X] = lx ;
      m_MPIperiodes[MPIGEO_WEST_BOTTOM][Y] = 0. ;
      m_MPIperiodes[MPIGEO_WEST_BOTTOM][Z] = 0. ;
    }

    if ( m_period[2] && m_coords[2] == 0 )
      m_MPIperiodes[MPIGEO_WEST_BOTTOM][Z] += lz ;
  }

  // West Top
  if ( m_voisins->rank( -1, 0, 1 ) != -1 )
  {
    if ( m_period[0] && m_coords[0] == 0 )
    {
      m_MPIperiodes[MPIGEO_WEST_TOP][X] = lx ;
      m_MPIperiodes[MPIGEO_WEST_TOP][Y] = 0. ;
      m_MPIperiodes[MPIGEO_WEST_TOP][Z] = 0. ;
    }

    if ( m_period[2] && m_coords[2] == m_dim[2] - 1 )
      m_MPIperiodes[MPIGEO_WEST_TOP][Z] += -lz ;
  }

  // East Bottom
  if ( m_voisins->rank( 1, 0, -1 ) != -1 )
  {
    if ( m_period[0] && m_coords[0] == m_dim[0] - 1 )
    {
      m_MPIperiodes[MPIGEO_EAST_BOTTOM][X] = -lx ;
      m_MPIperiodes[MPIGEO_EAST_BOTTOM][Y] = 0. ;
      m_MPIperiodes[MPIGEO_EAST_BOTTOM][Z] = 0. ;
    }

    if ( m_period[2] && m_coords[2] == 0 )
      m_MPIperiodes[MPIGEO_EAST_BOTTOM][Z] += lz ;
  }

  // East Top
  if ( m_voisins->rank( 1, 0, 1 ) != -1 )
  {
    if ( m_period[0] && m_coords[0] == m_dim[0] - 1 )
    {
      m_MPIperiodes[MPIGEO_EAST_TOP][X] = -lx ;
      m_MPIperiodes[MPIGEO_EAST_TOP][Y] = 0. ;
      m_MPIperiodes[MPIGEO_EAST_TOP][Z] = 0. ;
    }

    if ( m_period[2] && m_coords[2] == m_dim[2] - 1 )
      m_MPIperiodes[MPIGEO_EAST_TOP][Z] += -lz ;
  }


  // South West Bottom
  if ( m_voisins->rank( -1, -1, -1 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == 0 )
    {
      m_MPIperiodes[MPIGEO_SOUTH_WEST_BOTTOM][X] = 0. ;
      m_MPIperiodes[MPIGEO_SOUTH_WEST_BOTTOM][Y] = ly ;
      m_MPIperiodes[MPIGEO_SOUTH_WEST_BOTTOM][Z] = 0. ;
    }

    if ( m_period[0] && m_coords[0] == 0 )
      m_MPIperiodes[MPIGEO_SOUTH_WEST_BOTTOM][X] += lx ;

    if ( m_period[2] && m_coords[2] == 0 )
      m_MPIperiodes[MPIGEO_SOUTH_WEST_BOTTOM][Z] += lz ;
  }

  // South West Top
  if ( m_voisins->rank( -1, -1, 1 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == 0 )
    {
      m_MPIperiodes[MPIGEO_SOUTH_WEST_TOP][X] = 0. ;
      m_MPIperiodes[MPIGEO_SOUTH_WEST_TOP][Y] = ly ;
      m_MPIperiodes[MPIGEO_SOUTH_WEST_TOP][Z] = 0. ;
    }

    if ( m_period[0] && m_coords[0] == 0 )
      m_MPIperiodes[MPIGEO_SOUTH_WEST_TOP][X] += lx ;

    if ( m_period[2] && m_coords[2] == m_dim[2] - 1 )
      m_MPIperiodes[MPIGEO_SOUTH_WEST_TOP][Z] += -lz ;
  }

  // North West Bottom
  if ( m_voisins->rank( -1, 1, -1 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == m_dim[1] - 1 )
    {
      m_MPIperiodes[MPIGEO_NORTH_WEST_BOTTOM][X] = 0. ;
      m_MPIperiodes[MPIGEO_NORTH_WEST_BOTTOM][Y] = -ly ;
      m_MPIperiodes[MPIGEO_NORTH_WEST_BOTTOM][Z] = 0. ;
    }

    if ( m_period[0] && m_coords[0] == 0 )
      m_MPIperiodes[MPIGEO_NORTH_WEST_BOTTOM][X] += lx ;

    if ( m_period[2] && m_coords[2] == 0 )
      m_MPIperiodes[MPIGEO_NORTH_WEST_BOTTOM][Z] += lz ;
  }

  // North West Top
  if ( m_voisins->rank( -1, 1, 1 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == m_dim[1] - 1 )
    {
      m_MPIperiodes[MPIGEO_NORTH_WEST_TOP][X] = 0. ;
      m_MPIperiodes[MPIGEO_NORTH_WEST_TOP][Y] = -ly ;
      m_MPIperiodes[MPIGEO_NORTH_WEST_TOP][Z] = 0. ;
    }

    if ( m_period[0] && m_coords[0] == 0 )
      m_MPIperiodes[MPIGEO_NORTH_WEST_TOP][X] += lx ;

    if ( m_period[2] && m_coords[2] == m_dim[2] - 1 )
      m_MPIperiodes[MPIGEO_NORTH_WEST_TOP][Z] += -lz ;
  }

  // South East Bottom
  if ( m_voisins->rank( 1, -1, -1 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == 0 )
    {
      m_MPIperiodes[MPIGEO_SOUTH_EAST_BOTTOM][X] = 0. ;
      m_MPIperiodes[MPIGEO_SOUTH_EAST_BOTTOM][Y] = ly ;
      m_MPIperiodes[MPIGEO_SOUTH_EAST_BOTTOM][Z] = 0. ;
    }

    if ( m_period[0] && m_coords[0] == m_dim[0] - 1 )
      m_MPIperiodes[MPIGEO_SOUTH_EAST_BOTTOM][X] += -lx ;

    if ( m_period[2] && m_coords[2] == 0 )
      m_MPIperiodes[MPIGEO_SOUTH_EAST_BOTTOM][Z] += lz ;
  }

  // South East Top
  if ( m_voisins->rank( 1, -1, 1 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == 0 )
    {
      m_MPIperiodes[MPIGEO_SOUTH_EAST_TOP][X] = 0. ;
      m_MPIperiodes[MPIGEO_SOUTH_EAST_TOP][Y] = ly ;
      m_MPIperiodes[MPIGEO_SOUTH_EAST_TOP][Z] = 0. ;
    }

    if ( m_period[0] && m_coords[0] == m_dim[0] - 1 )
      m_MPIperiodes[MPIGEO_SOUTH_EAST_TOP][X] += -lx ;

    if ( m_period[2] && m_coords[2] == m_dim[2] - 1 )
      m_MPIperiodes[MPIGEO_SOUTH_EAST_TOP][Z] += -lz ;
  }

  // North East Bottom
  if ( m_voisins->rank( 1, 1, -1 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == m_dim[1] - 1 )
    {
      m_MPIperiodes[MPIGEO_NORTH_EAST_BOTTOM][X] = 0. ;
      m_MPIperiodes[MPIGEO_NORTH_EAST_BOTTOM][Y] = -ly ;
      m_MPIperiodes[MPIGEO_NORTH_EAST_BOTTOM][Z] = 0. ;
    }

    if ( m_period[0] && m_coords[0] == m_dim[0] - 1 )
      m_MPIperiodes[MPIGEO_NORTH_EAST_BOTTOM][X] += -lx ;

    if ( m_period[2] && m_coords[2] == 0 )
      m_MPIperiodes[MPIGEO_NORTH_EAST_BOTTOM][Z] += lz ;
  }

  // North East Top
  if ( m_voisins->rank( 1, 1, 1 ) != -1 )
  {
    if ( m_period[1] && m_coords[1] == m_dim[1] - 1 )
    {
      m_MPIperiodes[MPIGEO_NORTH_EAST_TOP][X] = 0. ;
      m_MPIperiodes[MPIGEO_NORTH_EAST_TOP][Y] = -ly ;
      m_MPIperiodes[MPIGEO_NORTH_EAST_TOP][Z] = 0. ;
    }

    if ( m_period[0] && m_coords[0] == m_dim[0] - 1 )
      m_MPIperiodes[MPIGEO_NORTH_EAST_TOP][X] += -lx ;

    if ( m_period[2] && m_coords[2] == m_dim[2] - 1 )
      m_MPIperiodes[MPIGEO_NORTH_EAST_TOP][Z] += -lz ;
  }
}




//-----------------------------------------------------------------------------
// Nombre de processus dans la direction i
int MPIWrapperGrains::nb_procs_direction(int i) const
{
  return m_dim[i];
}




//-----------------------------------------------------------------------------
// Nombre de processus dans la direction i
int const* MPIWrapperGrains::nb_procs_direction() const
{
  return m_dim;
}




//-----------------------------------------------------------------------------
// Nombre total de processus dans le communicateur MPI_COMM_WORLD
int MPIWrapperGrains::nombre_total_procs() const
{
  return m_nprocs;
}




//-----------------------------------------------------------------------------
// Nombre total de processus dans le communicateur MPI_COMM_WORLD
int MPIWrapperGrains::nombreTotalProcs()
{
  int nprocs_ = 0;
  MPI_Comm_size( MPI_COMM_WORLD, &nprocs_ );
  return nprocs_;
}




//-----------------------------------------------------------------------------
// Nombre de processus dans le communicateur m_MPI_COMM_activProc
int MPIWrapperGrains::nombre_total_procs_ACTIV() const
{
  int nprocs_ = 0;
  MPI_Comm_size( m_MPI_COMM_activProc, &nprocs_ );
  return nprocs_;
}




//-----------------------------------------------------------------------------
// Rang de processus dans le communicateur MPI_COMM_WORLD
int MPIWrapperGrains::rank_WORLD() const
{
  return m_rank_world;
}




//-----------------------------------------------------------------------------
// Rang de processus dans le communicateur m_MPI_COMM_activProc
int MPIWrapperGrains::rank_ACTIV() const
{
  return m_rank;
}




//-----------------------------------------------------------------------------
// Le processus est il actif ? */
bool MPIWrapperGrains::isActiv() const
{
  return m_is_activ;
}




//-----------------------------------------------------------------------------
// Rang de processus dans le communicateur MPI_COMM_WORLD
int MPIWrapperGrains::rankOf_WORLD()
{
  int rankproc = 0;
  MPI_Comm_rank( MPI_COMM_WORLD, &rankproc );
  return rankproc;
}




//-----------------------------------------------------------------------------
// Coordonnees du processus
int const* MPIWrapperGrains::MPI_coordonnees() const
{
  return m_coords;
}




//-----------------------------------------------------------------------------
// Voisins du processus
Voisins const* MPIWrapperGrains::MPI_voisins() const
{
  return m_voisins;
}




//-----------------------------------------------------------------------------
// Voisins du processus
int const* MPIWrapperGrains::MPI_periodicite() const
{
  return m_period;
}




//-----------------------------------------------------------------------------
// Affichage des caracteristiques
void MPIWrapperGrains::display( ostream &f ) const
{
  if ( m_rank == 0 )
  {
    f << endl << "Domain decomposition = ";
    for (int j=0;j<3;++j) cout << "N[" << j << "]=" << m_dim[j] << " ";
    f << endl;
    f << "MPI periods = ";
    for (int j=0;j<3;++j) cout << "P[" << j << "]=" << m_period[j] << " ";
    f << endl;
    App::affiche_attributs_statiques(f);
    f << endl;
  }

  for (int m=0;m<m_nprocs;++m)
  {
    if ( m == m_rank && m_is_activ )
    {
      f << "Processeur = " << m_rank << " PID = " << getpid() << endl;
      f << "Position dans la topologie = ";
      for (int j=0;j<3;++j) cout << m_coords[j] << " ";
      f << endl;
      f << "Processeurs voisins dans la topologie" << endl;
      f << "-------------------------------------" << endl;
      for (int i=-1;i<2;i++)
        for (int j=-1;j<2;j++)
          for (int k=-1;k<2;k++)
	    if ( m_voisins->rank( i, j, k ) != -1 )
            {
              f << "Neighbor (" << i << "," << j << "," << k << ") GEOLOC = "
	      	<< Cellule::getMPIGeoLocalisationName(
			getMPIGeoLocalisation( i, j, k ) ) << endl;
	      int const* coords_ = m_voisins->coordinates( i, j, k );
	      f << "Position in MPI topology = " << coords_[0] << " " <<
		coords_[1] << " " << coords_[2] << endl;
	      f << "Rank in MPI topology = " << m_voisins->rank( i, j, k )
	      	<< endl;
	      f << "MPI periodic vector = " <<
	      	m_MPIperiodes[ getMPIGeoLocalisation( i, j, k ) ][X] << " " <<
	      	m_MPIperiodes[ getMPIGeoLocalisation( i, j, k ) ][Y] << " " <<
	      	m_MPIperiodes[ getMPIGeoLocalisation( i, j, k ) ][Z] << endl;
            }
      f << endl;
    }
    MPI_Barrier( m_MPI_COMM_activProc );
  }

  if ( m_rank == 0 ) f << endl;
}




// ----------------------------------------------------------------------------
// Collecte sur le processeur master de l'ensemble des particules
// sur les diff�rents processeurs pour post-processing
// WARNING : NOT USED ANYMORE, NOW WE TRANSFERT PROPORTIES, NOT PARTICLES
vector<Particule*>* MPIWrapperGrains::GatherParticules_PostProcessing(
  	const list<Particule*> &particules,
	const list<Particule*> &pwait,
	vector<Particule*> const& ParticuleClassesReference,
	const size_t& nb_total_particules) const
{
  vector<Particule*>* allparticules=NULL;
  list<Particule*>::const_iterator il;
  int i,j;

  // Taille des messages proportionnel au nombre de particules dans
  // ParticulesReference
  // ---------------------------------------------------------------
  int nb_part = int(particules.size());
  for (il=particules.begin();il!=particules.end();il++)
    if ( (*il)->getTag() == 2 ) nb_part--;

  // Taille des messages � passer
  int *recvcounts = new int[m_nprocs];
  for (i=0; i<m_nprocs; ++i) recvcounts[i] = 0;
  MPI_Gather( &nb_part, 1, MPI_INT, recvcounts, 1, MPI_INT,
  	m_rank_masterWorld, m_MPI_COMM_activProc );

  // Partition et taille du buffer de r�ception
  int *displs = new int[m_nprocs];
  displs[0] = 0;
  for (i=1; i<m_nprocs; ++i) displs[i] = displs[i-1] + recvcounts[i-1];
  int recvsize = displs[m_nprocs-1] + recvcounts[m_nprocs-1];



  // Communication des entiers:
  // Ordre par particule:
  // [numero de particule, classe]
  // -----------------------------
  int NB_INT_PART = 3;
  int *numClass = new int[NB_INT_PART*nb_part];
  for (il=particules.begin(),i=0;il!=particules.end();
  	il++)
  {
    if ( (*il)->getTag() == 0 || (*il)->getTag() == 1 )
    {
      numClass[i] = (*il)->getID();
      numClass[i+1] = (*il)->getParticuleClasse();
      numClass[i+2] = (*il)->getCoordinationNumber();
      i += NB_INT_PART;
    }
  }

  // Partition et taille du buffer de r�ception
  int *recvcounts_INT = new int[m_nprocs];
  for (i=0; i<m_nprocs; ++i) recvcounts_INT[i] = NB_INT_PART * recvcounts[i];
  int *displs_INT = new int[m_nprocs];
  displs_INT[0] = 0;
  for (i=1; i<m_nprocs; ++i)
    displs_INT[i] = displs_INT[i-1] + recvcounts_INT[i-1];
  int recvsize_INT = displs_INT[m_nprocs-1] + recvcounts_INT[m_nprocs-1];
  int *recvbuf_INT = new int[recvsize_INT];

  // Communication des entiers: 3 integer par particule
  MPI_Gatherv( numClass, NB_INT_PART * nb_part, MPI_INT, recvbuf_INT,
  	recvcounts_INT, displs_INT, MPI_INT,
	m_rank_masterWorld, m_MPI_COMM_activProc );



  // Communication des doubles: cin�matique & configuration
  // Ordre par particule:
  // [position, vitesse translation, quaternion rotation, vitesse rotation]
  // ----------------------------------------------------------------------
  int NB_DOUBLE_PART = 26;
  double *features = new double[NB_DOUBLE_PART*nb_part];
  for (il=particules.begin(),i=0;il!=particules.end();
  	il++)
  {
    if ( (*il)->getTag() == 0 || (*il)->getTag() == 1 )
    {
      (*il)->copyVitesseTranslation( features, i );
      (*il)->copyQuaternionRotation( features, i+3 );
      (*il)->copyVitesseRotation( features, i+7 );
      (*il)->copyTransform( features, i+10 );
      i += NB_DOUBLE_PART;
    }
  }

  // Partition et taille du buffer de r�ception
  int *recvcounts_DOUBLE = new int[m_nprocs];
  for (i=0; i<m_nprocs; ++i)
    recvcounts_DOUBLE[i] = NB_DOUBLE_PART * recvcounts[i];
  int *displs_DOUBLE = new int[m_nprocs];
  displs_DOUBLE[0] = 0;
  for (i=1; i<m_nprocs; ++i)
    displs_DOUBLE[i] = displs_DOUBLE[i-1] + recvcounts_DOUBLE[i-1];
  int recvsize_DOUBLE = displs_DOUBLE[m_nprocs-1]
  	+ recvcounts_DOUBLE[m_nprocs-1];
  double *recvbuf_DOUBLE = new double[recvsize_DOUBLE];

  // Communication de la cin�matique & configuration: 26 double par
  // particule
  MPI_Gatherv( features, NB_DOUBLE_PART * nb_part, MPI_DOUBLE, recvbuf_DOUBLE,
  	recvcounts_DOUBLE, displs_DOUBLE, MPI_DOUBLE,
	m_rank_masterWorld, m_MPI_COMM_activProc );


  // Creation des Particules pour Post-processing
  // --------------------------------------------
  if ( m_rank == m_rank_masterWorld )
  {
    Particule* particule = NULL;
    allparticules = new vector<Particule*>( nb_total_particules, particule );

    // Particules actives sur tous les proc
    for (j=0;j<recvsize;++j)
    {
      // Creation de la particule
      if ( (ParticuleClassesReference)[recvbuf_INT[NB_INT_PART*j+1]]
	->isCompParticule() )
      {
	Particule *part_post = new CompParticule( recvbuf_INT[NB_INT_PART*j],
      		ParticuleClassesReference[recvbuf_INT[NB_INT_PART*j+1]],
      		recvbuf_DOUBLE[NB_DOUBLE_PART*j],
      		recvbuf_DOUBLE[NB_DOUBLE_PART*j+1],
    		recvbuf_DOUBLE[NB_DOUBLE_PART*j+2],
      		recvbuf_DOUBLE[NB_DOUBLE_PART*j+3],
      		recvbuf_DOUBLE[NB_DOUBLE_PART*j+4],
    		recvbuf_DOUBLE[NB_DOUBLE_PART*j+5],
      		recvbuf_DOUBLE[NB_DOUBLE_PART*j+6],
      		recvbuf_DOUBLE[NB_DOUBLE_PART*j+7],
    		recvbuf_DOUBLE[NB_DOUBLE_PART*j+8],
    		recvbuf_DOUBLE[NB_DOUBLE_PART*j+9],
      		&recvbuf_DOUBLE[NB_DOUBLE_PART*j+10],
		COMPUTE,
		0,
		recvbuf_INT[NB_INT_PART*j+2] );
	(*allparticules)[recvbuf_INT[NB_INT_PART*j]] = part_post;
      }
      else
      {
	Particule *part_post = new Particule( recvbuf_INT[NB_INT_PART*j],
      		ParticuleClassesReference[recvbuf_INT[NB_INT_PART*j+1]],
      		recvbuf_DOUBLE[NB_DOUBLE_PART*j],
      		recvbuf_DOUBLE[NB_DOUBLE_PART*j+1],
    		recvbuf_DOUBLE[NB_DOUBLE_PART*j+2],
      		recvbuf_DOUBLE[NB_DOUBLE_PART*j+3],
      		recvbuf_DOUBLE[NB_DOUBLE_PART*j+4],
    		recvbuf_DOUBLE[NB_DOUBLE_PART*j+5],
      		recvbuf_DOUBLE[NB_DOUBLE_PART*j+6],
      		recvbuf_DOUBLE[NB_DOUBLE_PART*j+7],
    		recvbuf_DOUBLE[NB_DOUBLE_PART*j+8],
    		recvbuf_DOUBLE[NB_DOUBLE_PART*j+9],
      		&recvbuf_DOUBLE[NB_DOUBLE_PART*j+10],
		COMPUTE,
		0,
		recvbuf_INT[NB_INT_PART*j+2] );
	(*allparticules)[recvbuf_INT[NB_INT_PART*j]] = part_post;
      }
    }

    // Particules en attente
    if ( !pwait.empty() )
    {
      double *vt = new double[3];
      double *vrot = new double[3];
      double *qrot = new double[4];
      double *transform = new double[16];
      for (il=pwait.begin();il!=pwait.end();il++)
      {
        (*il)->copyVitesseTranslation( vt, 0 );
        (*il)->copyQuaternionRotation( qrot, 0 );
        (*il)->copyVitesseRotation( vrot, 0 );
        (*il)->copyTransform( transform, 0 );

        // Creation de la particule
	if (ParticuleClassesReference[(*il)->getParticuleClasse()]
	  ->isCompParticule())
	{
	  Particule *part_post = new CompParticule( (*il)->getID(),
      		ParticuleClassesReference[(*il)->getParticuleClasse()],
      		vt[0],vt[1],vt[2],
      		qrot[0],qrot[1],qrot[2],qrot[3],
		vrot[0],vrot[1],vrot[2],
      		transform,
		WAIT,
		0 );
	  (*allparticules)[(*il)->getID()] = part_post;
	}
	else
	{
	  Particule *part_post = new Particule( (*il)->getID(),
      		ParticuleClassesReference[(*il)->getParticuleClasse()],
      		vt[0],vt[1],vt[2],
      		qrot[0],qrot[1],qrot[2],qrot[3],
		vrot[0],vrot[1],vrot[2],
      		transform,
		WAIT,
		0 );
	  (*allparticules)[(*il)->getID()] = part_post;
	}
      }
      delete [] vt;
      delete [] vrot;
      delete [] qrot;
      delete [] transform;
    }
  }

  // Debug
//   if ( m_rank == m_rank_masterWorld )
//     for (i=0;i<int((*allparticules).size());++i)
//       if ( (*allparticules)[i] == NULL ) cout << "Part " << i << endl;
  // End debug

  delete [] numClass;
  delete [] features;
  delete [] recvcounts;
  delete [] recvcounts_INT;
  delete [] recvcounts_DOUBLE;
  delete [] displs;
  delete [] displs_INT;
  delete [] displs_DOUBLE;
  delete [] recvbuf_INT;
  delete [] recvbuf_DOUBLE;

  return allparticules;
}




// ----------------------------------------------------------------------------
// Collecte sur le processeur master de l'ensemble des positions, vitesses de
// translation et vitesses de rotation des particules
// sur les diff�rents processeurs pour post-processing
// Thus, we do not need to create then delete all particles on the master proc
// M.BERNARD - 2013
vector< vector<double> >* MPIWrapperGrains::
    GatherPositionVitesse_PostProcessing(
        const list<Particule*> &particules,
        const size_t& nb_total_particules ) const
{
  bool b_totalForceOP = Text_PostProcessingWriter::b_totalForce;
  bool b_cumulatedContactForceOP = Text_PostProcessingWriter::
      b_cumulatedContactForce;
  bool b_instantaneousContactForceOP = Text_PostProcessingWriter::
      b_instantaneousContactForce;
  bool b_cumulatedLubriForceOP = Text_PostProcessingWriter::
      b_cumulatedLubriForce;
  bool b_hydroForceOP = Text_PostProcessingWriter::b_hydroForce;
  bool b_slipVelocityOP = Text_PostProcessingWriter::b_slipVelocity;
  bool b_temperatureOP = Text_PostProcessingWriter::b_temperature;
  bool b_stressTensor = Grains_Exec::m_stressTensor; // Macroscale
  bool b_particleStressTensor =
    Text_PostProcessingWriter::b_particleStressTensor; // Microscale
  int nb_infos=0, nTF=0, nCCF=0, nICF=0, nHF=0, nSV=0, nT=0, nCLF=0,
	 nSTM=0;

  Vecteur const* vitesseT = NULL;
  Vecteur const* vitesseR = NULL;
  Vecteur const* force = NULL;
  Vecteur const* contactF = NULL;
  Vecteur const* contactF_inst = NULL;
  Vecteur const* lubriF = NULL;
  Vecteur const* demcfdHydroForce = NULL;
  Vecteur const* demcfdSlipVelocity = NULL;
  double const* temperature = NULL;
  // Internal moments or individual stress tensor
  vector<Scalar> const* InternalFeatures = NULL;
  double const* heatflux = NULL;
  double const* nusselt = NULL;
  double const* fluidTemperature = NULL;

  Point const* position = NULL;
  double intTodouble = 0.1 ;
  int i=0, tag_DOUBLE = 1, recvsize = 0, ID_part=0;
  MPI_Status status;
  MPI_Request idreq;

  vector< vector<double> >* cinematique_Global = NULL;

  list<Particule*>::const_iterator il;
  int nb_part_loc = int(particules.size());
  // On comptabilise pas les particules de la halo zone
  for (il=particules.begin();il!=particules.end();il++)
    if ((*il)->getTag()==2) nb_part_loc--;

  // infos : ID + 3xPosition + 3xVTrans + 3xVRot + nbContact
  nb_infos = 11;
  if( b_totalForceOP ) nTF = 3;
  if( b_cumulatedContactForceOP ) nCCF = 3;
  if( b_instantaneousContactForceOP ) nICF = 3;
  if( b_hydroForceOP ) nHF = 3;
  if( b_slipVelocityOP ) nSV = 3;
  if( b_temperatureOP ) nT = 4;
  if( b_cumulatedLubriForceOP ) nCLF = 3;
  if( b_stressTensor || b_particleStressTensor ) nSTM = 9;
  nb_infos += nTF+nCCF+nICF+nHF+nSV+nT+nCLF+nSTM;

  // Taille des messages proportionnelle au nombre de particules
  // sur chaque proc
  double *buffer = new double[nb_infos*nb_part_loc];

  for (il=particules.begin(), i=0; il!=particules.end(); il++, i+=nb_infos)
  {
    if( (*il)->getTag() == 2 ) i-=nb_infos;
    else
    {
      position = (*il)->getPosition();
      vitesseT = (*il)->getVitesseTranslation();
      vitesseR = (*il)->getVitesseRotation();

      buffer[i] = (*il)->getID() + intTodouble ;
      buffer[i+1] = (*position)[0];
      buffer[i+2] = (*position)[1];
      buffer[i+3] = (*position)[2];
      buffer[i+4] = (*vitesseT)[0];
      buffer[i+5] = (*vitesseT)[1];
      buffer[i+6] = (*vitesseT)[2];
      buffer[i+7] = (*vitesseR)[0];
      buffer[i+8] = (*vitesseR)[1];
      buffer[i+9] = (*vitesseR)[2];
      buffer[i+10] = (*il)->getCoordinationNumber();
      if( b_totalForceOP )
      {
        force = (*il)->getForce();
        buffer[i+11] = (*force)[0];
        buffer[i+12] = (*force)[1];
        buffer[i+13] = (*force)[2];
      }
      if( b_cumulatedContactForceOP )
      {
        contactF = (*il)->getForceContactPP();
        buffer[i+11+nTF] = (*contactF)[0];
        buffer[i+12+nTF] = (*contactF)[1];
        buffer[i+13+nTF] = (*contactF)[2];
      }
      if( b_instantaneousContactForceOP )
      {
        contactF_inst = (*il)->getForceContactPP_instantaneous();
        buffer[i+11+nTF+nCCF] = (*contactF_inst)[0];
        buffer[i+12+nTF+nCCF] = (*contactF_inst)[1];
        buffer[i+13+nTF+nCCF] = (*contactF_inst)[2];
      }
      if( b_hydroForceOP )
      {
        demcfdHydroForce = (*il)->getParticleHydroForce();
        buffer[i+11+nTF+nCCF+nICF] = (*demcfdHydroForce)[0];
        buffer[i+12+nTF+nCCF+nICF] = (*demcfdHydroForce)[1];
        buffer[i+13+nTF+nCCF+nICF] = (*demcfdHydroForce)[2];
      }
      if( b_slipVelocityOP )
      {
        demcfdSlipVelocity = (*il)->getParticleSlipVel();
        buffer[i+11+nTF+nCCF+nICF+nHF] = (*demcfdSlipVelocity)[0];
        buffer[i+12+nTF+nCCF+nICF+nHF] = (*demcfdSlipVelocity)[1];
        buffer[i+13+nTF+nCCF+nICF+nHF] = (*demcfdSlipVelocity)[2];
      }
      if( b_temperatureOP )
      {
        temperature = (*il)->get_solidTemperature();
	heatflux = (*il)->get_fluidSolidHeatFlux();
	nusselt = (*il)->get_solidNusselt();
	fluidTemperature = (*il)->get_DEMCFD_fluidTemperature();
//        cout << " TEMPORARY will copy tempS in buffer "
//             << *temperature << endl;
        buffer[i+11+nTF+nCCF+nICF+nHF+nSV] = *temperature;
	buffer[i+12+nTF+nCCF+nICF+nHF+nSV] = *heatflux;
	buffer[i+13+nTF+nCCF+nICF+nHF+nSV] = *nusselt;
	buffer[i+14+nTF+nCCF+nICF+nHF+nSV] = *fluidTemperature;
      }
      if( b_cumulatedLubriForceOP )
      {
        lubriF = (*il)->getForceLubriPP();
        buffer[i+11+nTF+nCCF+nICF+nHF+nSV+nT] = (*lubriF)[0];
        buffer[i+12+nTF+nCCF+nICF+nHF+nSV+nT] = (*lubriF)[1];
        buffer[i+13+nTF+nCCF+nICF+nHF+nSV+nT] = (*lubriF)[2];
      }
      if( b_stressTensor || b_particleStressTensor )
      {
	if ( b_stressTensor ) InternalFeatures = (*il)->getInternalMoment();
	else InternalFeatures = (*il)->getStressTensor();
        buffer[i+11+nTF+nCCF+nICF+nHF+nSV+nT+nCLF] = (*InternalFeatures)[0];
        buffer[i+12+nTF+nCCF+nICF+nHF+nSV+nT+nCLF] = (*InternalFeatures)[1];
        buffer[i+13+nTF+nCCF+nICF+nHF+nSV+nT+nCLF] = (*InternalFeatures)[2];
        buffer[i+14+nTF+nCCF+nICF+nHF+nSV+nT+nCLF] = (*InternalFeatures)[3];
        buffer[i+15+nTF+nCCF+nICF+nHF+nSV+nT+nCLF] = (*InternalFeatures)[4];
        buffer[i+16+nTF+nCCF+nICF+nHF+nSV+nT+nCLF] = (*InternalFeatures)[5];
        buffer[i+17+nTF+nCCF+nICF+nHF+nSV+nT+nCLF] = (*InternalFeatures)[6];
        buffer[i+18+nTF+nCCF+nICF+nHF+nSV+nT+nCLF] = (*InternalFeatures)[7];
        buffer[i+19+nTF+nCCF+nICF+nHF+nSV+nT+nCLF] = (*InternalFeatures)[8];
      }
    }
  }
  MPI_Isend( buffer, nb_infos*nb_part_loc, MPI_DOUBLE, m_rank_masterWorld,
      tag_DOUBLE, m_MPI_COMM_activProc, &idreq );


  // RECEPTION par le Master
  // -----------------------
  if( m_rank == m_rank_masterWorld )
  {
    vector<double> work( nb_total_particules, 0. ) ;
    cinematique_Global = new vector< vector<double> >( nb_infos-1, work ) ;

    for (int irank=0; irank<m_nprocs; ++irank)
    {
      // Evaluation de la taille du message
      MPI_Probe( irank, tag_DOUBLE, m_MPI_COMM_activProc, &status );
      MPI_Get_count( &status, MPI_DOUBLE, &recvsize );

      // Reception du message de doubles
      double *recvbuf_DOUBLE = new double[recvsize];
      MPI_Recv( recvbuf_DOUBLE, recvsize, MPI_DOUBLE,
	irank, tag_DOUBLE, m_MPI_COMM_activProc, &status );

      // Copie dans cinematique_Global
      for(int j=0; j<recvsize; j+=nb_infos)
      {
        ID_part = int(recvbuf_DOUBLE[j]) ;
        (*cinematique_Global)[0][ID_part] = recvbuf_DOUBLE[j+1];
        (*cinematique_Global)[1][ID_part] = recvbuf_DOUBLE[j+2];
        (*cinematique_Global)[2][ID_part] = recvbuf_DOUBLE[j+3];
        (*cinematique_Global)[3][ID_part] = recvbuf_DOUBLE[j+4];
        (*cinematique_Global)[4][ID_part] = recvbuf_DOUBLE[j+5];
        (*cinematique_Global)[5][ID_part] = recvbuf_DOUBLE[j+6];
        (*cinematique_Global)[6][ID_part] = recvbuf_DOUBLE[j+7];
        (*cinematique_Global)[7][ID_part] = recvbuf_DOUBLE[j+8];
        (*cinematique_Global)[8][ID_part] = recvbuf_DOUBLE[j+9];
        (*cinematique_Global)[9][ID_part] = recvbuf_DOUBLE[j+10];
        if( b_totalForceOP )
        {
          (*cinematique_Global)[10][ID_part] = recvbuf_DOUBLE[j+11];
          (*cinematique_Global)[11][ID_part] = recvbuf_DOUBLE[j+12];
          (*cinematique_Global)[12][ID_part] = recvbuf_DOUBLE[j+13];
        }
        if( b_cumulatedContactForceOP )
        {
          (*cinematique_Global)[10+nTF][ID_part] = recvbuf_DOUBLE[j+11+nTF];
          (*cinematique_Global)[11+nTF][ID_part] = recvbuf_DOUBLE[j+12+nTF];
          (*cinematique_Global)[12+nTF][ID_part] = recvbuf_DOUBLE[j+13+nTF];
        }
        if( b_instantaneousContactForceOP )
        {
          (*cinematique_Global)[10+nTF+nCCF][ID_part] = recvbuf_DOUBLE[j+11+nTF+nCCF];
          (*cinematique_Global)[11+nTF+nCCF][ID_part] = recvbuf_DOUBLE[j+12+nTF+nCCF];
          (*cinematique_Global)[12+nTF+nCCF][ID_part] = recvbuf_DOUBLE[j+13+nTF+nCCF];
        }
        if( b_hydroForceOP )
        {
          (*cinematique_Global)[10+nTF+nICF+nCCF][ID_part] =
              recvbuf_DOUBLE[j+11+nTF+nICF+nCCF];
          (*cinematique_Global)[11+nTF+nICF+nCCF][ID_part] =
              recvbuf_DOUBLE[j+12+nTF+nICF+nCCF];
          (*cinematique_Global)[12+nTF+nICF+nCCF][ID_part] =
              recvbuf_DOUBLE[j+13+nTF+nICF+nCCF];
        }
        if( b_slipVelocityOP )
        {
          (*cinematique_Global)[10+nTF+nCCF+nICF+nHF][ID_part] =
              recvbuf_DOUBLE[j+11+nTF+nCCF+nICF+nHF];
          (*cinematique_Global)[11+nTF+nCCF+nICF+nHF][ID_part] =
              recvbuf_DOUBLE[j+12+nTF+nCCF+nICF+nHF];
          (*cinematique_Global)[12+nTF+nCCF+nICF+nHF][ID_part] =
              recvbuf_DOUBLE[j+13+nTF+nCCF+nICF+nHF];
        }
        if( b_temperatureOP )
        {
          (*cinematique_Global)[10+nTF+nCCF+nICF+nHF+nSV][ID_part] =
              recvbuf_DOUBLE[j+11+nTF+nCCF+nICF+nHF+nSV];
	  (*cinematique_Global)[11+nTF+nCCF+nHF+nSV][ID_part] =
              recvbuf_DOUBLE[j+12+nTF+nCCF+nHF+nSV];
	  (*cinematique_Global)[12+nTF+nCCF+nHF+nSV][ID_part] =
              recvbuf_DOUBLE[j+13+nTF+nCCF+nHF+nSV];
	  (*cinematique_Global)[13+nTF+nCCF+nHF+nSV][ID_part] =
              recvbuf_DOUBLE[j+14+nTF+nCCF+nHF+nSV];
        }
        if( b_cumulatedLubriForceOP )
        {
          (*cinematique_Global)[10+nTF+nCCF+nICF+nHF+nSV+nT][ID_part] =
	      recvbuf_DOUBLE[j+11+nTF+nCCF+nICF+nHF+nSV+nT];
          (*cinematique_Global)[11+nTF+nCCF+nICF+nHF+nSV+nT][ID_part] =
	      recvbuf_DOUBLE[j+12+nTF+nCCF+nICF+nHF+nSV+nT];
          (*cinematique_Global)[12+nTF+nCCF+nICF+nHF+nSV+nT][ID_part] =
	      recvbuf_DOUBLE[j+13+nTF+nCCF+nICF+nHF+nSV+nT];
        }
	if( b_stressTensor || b_particleStressTensor )
	{
	  // Particle internal moments divided by the volume (function of
	  // the time). It is a 3x3 matrix.
          (*cinematique_Global)[10+nTF+nCCF+nICF+nHF+nSV+nT+nCLF][ID_part] =
	      recvbuf_DOUBLE[j+11+nTF+nCCF+nICF+nHF+nSV+nT+nCLF];
          (*cinematique_Global)[11+nTF+nCCF+nICF+nHF+nSV+nT+nCLF][ID_part] =
	      recvbuf_DOUBLE[j+12+nTF+nCCF+nICF+nHF+nSV+nT+nCLF];
          (*cinematique_Global)[12+nTF+nCCF+nICF+nHF+nSV+nT+nCLF][ID_part] =
	      recvbuf_DOUBLE[j+13+nTF+nCCF+nICF+nHF+nSV+nT+nCLF];
          (*cinematique_Global)[13+nTF+nCCF+nICF+nHF+nSV+nT+nCLF][ID_part] =
	      recvbuf_DOUBLE[j+14+nTF+nCCF+nICF+nHF+nSV+nT+nCLF];
          (*cinematique_Global)[14+nTF+nCCF+nICF+nHF+nSV+nT+nCLF][ID_part] =
	      recvbuf_DOUBLE[j+15+nTF+nCCF+nICF+nHF+nSV+nT+nCLF];
          (*cinematique_Global)[15+nTF+nCCF+nICF+nHF+nSV+nT+nCLF][ID_part] =
	      recvbuf_DOUBLE[j+16+nTF+nCCF+nICF+nHF+nSV+nT+nCLF];
          (*cinematique_Global)[16+nTF+nCCF+nICF+nHF+nSV+nT+nCLF][ID_part] =
	      recvbuf_DOUBLE[j+17+nTF+nCCF+nICF+nHF+nSV+nT+nCLF];
          (*cinematique_Global)[17+nTF+nCCF+nICF+nHF+nSV+nT+nCLF][ID_part] =
	      recvbuf_DOUBLE[j+18+nTF+nCCF+nICF+nHF+nSV+nT+nCLF];
          (*cinematique_Global)[18+nTF+nCCF+nICF+nHF+nSV+nT+nCLF][ID_part] =
	      recvbuf_DOUBLE[j+19+nTF+nCCF+nICF+nHF+nSV+nT+nCLF];
	}
      }
      delete [] recvbuf_DOUBLE;
    }
  }

  // Verifie que les envois non bloquants sont termin�s
  MPI_Wait( &idreq, &status );

  delete [] buffer;

  return cinematique_Global;

}




// ----------------------------------------------------------------------------
// Collecte sur le processeur master de la classe de toutes les particules
vector< vector<double> >* MPIWrapperGrains::
    GatherParticlesClass_PostProcessing(
    const list<Particule*> &particules,
    const size_t& nb_total_particules ) const
{
  double intTodouble = 0.1 ;
  int i=0, tag_DOUBLE = 1, recvsize = 0, ID_part=0;
  MPI_Status status;
  MPI_Request idreq;

  vector< vector<double> >* class_Global = NULL;

  list<Particule*>::const_iterator il;
  int nb_part_loc = int(particules.size());

  // We do not care about particles in  halozone
  for (il=particules.begin();il!=particules.end();il++)
    if ((*il)->getTag()==2) nb_part_loc--;

  // Buffer size depend on the number of particles per core
  double *buffer = new double[2*nb_part_loc];

  for (il=particules.begin(), i=0; il!=particules.end(); il++, i+=2)
  {
    if( (*il)->getTag()==2 ) i-=2;
    else
    {
      buffer[i] = (*il)->getID() + intTodouble ;
      buffer[i+1] = (*il)->getParticuleClasse();
    }
  }

  MPI_Isend( buffer, 2*nb_part_loc, MPI_DOUBLE, m_rank_masterWorld,
      tag_DOUBLE, m_MPI_COMM_activProc, &idreq );

  // RECEPTION par le Master
  // -----------------------
  if( m_rank == m_rank_masterWorld )
  {
    vector<double> work( nb_total_particules, 0. ) ;
    class_Global = new vector< vector<double> >( 1, work ) ;

    for (int irank=0; irank<m_nprocs; ++irank)
    {
      // Evaluation de la taille du message
      MPI_Probe( irank, tag_DOUBLE, m_MPI_COMM_activProc, &status );
      MPI_Get_count( &status, MPI_DOUBLE, &recvsize );

      // Reception du message de doubles
      double *recvbuf_DOUBLE = new double[recvsize];
      MPI_Recv( recvbuf_DOUBLE, recvsize, MPI_DOUBLE,
          irank, tag_DOUBLE, m_MPI_COMM_activProc, &status );

      // Copie dans class_Global
      for(int j=0; j<recvsize; j+=2)
      {
        ID_part = int(recvbuf_DOUBLE[j]) ;
        (*class_Global)[0][ID_part] = recvbuf_DOUBLE[j+1];
      }

      delete [] recvbuf_DOUBLE;
    }
  }

  // Verifie que les envois non bloquants sont termin�s
  MPI_Wait( &idreq, &status );

  delete [] buffer ;

  return class_Global;

}




// ----------------------------------------------------------------------------
// Cr�ation des nouveaux clones
void MPIWrapperGrains::UpdateOrCreateClones_AllGatherGlobal(Scalar time,
    list<Particule*>* particulesClones,
    list<Particule*>* particules,
    list<Particule*> const* particulesHalozone,
    vector<Particule*> const* ParticuleClassesReference,
    LinkedCell* LC)
{
  list<Particule*>::const_iterator il;
  int i;
  bool AdamsBashforth = Grains_Exec::m_TIScheme == "SecondOrderAdamsBashforth";
  bool b_hydroForce = Grains_Exec::m_withHydroForce ;
  bool b_liftForce = Grains_Exec::m_withLiftForce ;
  double intTodouble = 0.1 ;

  SCT_set_start("MPIComm");

  // Remplissage de la multimap pour acces aux clones par numero
  // -----------------------------------------------------------
  AccessToClones.clear();
  for (il=particulesClones->begin();il!=particulesClones->end();il++)
    AccessToClones.insert( pair<int,Particule*>( (*il)->getID(), *il ) );

  // Taille des messages proportionnel au nombre de particules dans
  // particulesHalozone
  // --------------------------------------------------------------
  int nb_hz = int(particulesHalozone->size());

  // Taille des messages � passer
  int *recvcounts = new int[m_nprocs];
  for (i=0; i<m_nprocs; ++i) recvcounts[i]=0;
  MPI_Allgather( &nb_hz, 1, MPI_INT, recvcounts, 1, MPI_INT,
  	m_MPI_COMM_activProc );

  // Partition et taille du buffer de r�ception
  int *displs = new int[m_nprocs];
  displs[0] = 0;
  for (i=1; i<m_nprocs; ++i) displs[i] = displs[i-1] + recvcounts[i-1];
  int recvsize = displs[m_nprocs-1] + recvcounts[m_nprocs-1];

  SCT_get_elapsed_time( "MPIComm" );
  SCT_set_start( "Copie_Buffers" );

  // Buffer de doubles: cin�matique & configuration
  // Ordre par particule: [numero de particule, classe, rang exp�diteur,
  //	position, vitesse translation, quaternion rotation,
  // 	vitesse rotation]
  // -------------------------------------------------------------------------
//  int NB_DOUBLE_PART = AdamsBashforth ? 41 : 29;
  int NB_DOUBLE_PART = 29;
  if( AdamsBashforth ) NB_DOUBLE_PART += 12;
  if( b_hydroForce ) NB_DOUBLE_PART += 7;
  if( b_liftForce ) NB_DOUBLE_PART += 3;
  double *features = new double[NB_DOUBLE_PART*nb_hz];
  for (il=particulesHalozone->begin(),i=0;il!=particulesHalozone->end();
      il++,i+=NB_DOUBLE_PART)
  {
    features[i] = (*il)->getID() + intTodouble;
    features[i+1] = (*il)->getParticuleClasse() + intTodouble;
    features[i+2] = m_rank + intTodouble;
    (*il)->copyVitesseTranslation( features, i+3 );
    (*il)->copyQuaternionRotation( features, i+6 );
    (*il)->copyVitesseRotation( features, i+10 );
    (*il)->copyTransform( features, i+13 );
    if( AdamsBashforth )
    {
      (*il)->copyCinematiqueNm2( features, i+29 );
      if ( b_hydroForce )
        (*il)->copyFluidInformations( features, i+41 );
    }
    else if( b_hydroForce )
       (*il)->copyFluidInformations( features, i+29 );
  }

  SCT_get_elapsed_time( "Copie_Buffers" );
  SCT_set_start( "MPIComm" );

  // Partition et taille du buffer de r�ception
  int *recvcounts_DOUBLE = new int[m_nprocs];
  for (i=0; i<m_nprocs; ++i)
    recvcounts_DOUBLE[i] = NB_DOUBLE_PART * recvcounts[i];
  int *displs_DOUBLE = new int[m_nprocs];
  displs_DOUBLE[0] = 0;
  for (i=1; i<m_nprocs; ++i)
    displs_DOUBLE[i] = displs_DOUBLE[i-1] + recvcounts_DOUBLE[i-1];
  int recvsize_DOUBLE = displs_DOUBLE[m_nprocs-1]
  	+ recvcounts_DOUBLE[m_nprocs-1];
  double *recvbuf_DOUBLE = new double[recvsize_DOUBLE];

  // Communication de la cin�matique & configuration: 26 double par
  // particule
  MPI_Allgatherv( features, NB_DOUBLE_PART*nb_hz, MPI_DOUBLE, recvbuf_DOUBLE,
  	recvcounts_DOUBLE, displs_DOUBLE, MPI_DOUBLE, m_MPI_COMM_activProc );

  SCT_add_elapsed_time("MPIComm");
  SCT_set_start("UpdateCreateClones");

  // Creation ou maj des clones
  // --------------------------
  UpdateOrCreateClones( time, recvsize, recvbuf_DOUBLE,
	NB_DOUBLE_PART, particulesClones,
	particules, particulesHalozone, ParticuleClassesReference, LC );

  delete [] features;
  delete [] recvcounts;
  delete [] recvcounts_DOUBLE;
  delete [] displs;
  delete [] displs_DOUBLE;
  delete [] recvbuf_DOUBLE;

  SCT_get_elapsed_time( "UpdateCreateClones" );
}





// ----------------------------------------------------------------------------
// Cr�ation des nouveaux clones
void MPIWrapperGrains::UpdateOrCreateClones_AllGatherLocal(Scalar time,
	list<Particule*>* particulesClones,
	list<Particule*>* particules,
  	list<Particule*> const* particulesHalozone,
	vector<Particule*> const* ParticuleClassesReference,
	LinkedCell* LC)
{
  list<Particule*>::const_iterator il;
  int i, ii;
  bool AdamsBashforth = Grains_Exec::m_TIScheme == "SecondOrderAdamsBashforth";
  bool b_hydroForce = Grains_Exec::m_withHydroForce ;
  bool b_liftForce = Grains_Exec::m_withLiftForce ;
  double intTodouble = 0.1 ;

  SCT_set_start("MPIComm");

  // Remplissage de la multimap pour acces aux clones par numero
  // -----------------------------------------------------------
  AccessToClones.clear();
  for (il=particulesClones->begin();il!=particulesClones->end();il++)
    AccessToClones.insert( pair<int,Particule*>( (*il)->getID(), *il ) );

  // Taille des messages proportionnel au nombre de particules dans
  // particulesHalozone
  // --------------------------------------------------------------
  int nb_hz = int(particulesHalozone->size());

  // Taille des messages � passer
  int *recvcounts = new int[m_nprocs_localComm];
  for (i=0; i<m_nprocs_localComm; ++i) recvcounts[i]=0;
  for (ii=0;ii<m_nprocs;++ii)
    if ( m_isInCommVoisins[ii] )
      MPI_Gather( &nb_hz, 1, MPI_INT, recvcounts, 1, MPI_INT,
      	m_master_localComm[ii], *m_commVoisins[ii] );

  // Partition et taille du buffer de r�ception
  int *displs = new int[m_nprocs_localComm];
  displs[0] = 0;
  for (i=1; i<m_nprocs_localComm; ++i)
    displs[i] = displs[i-1] + recvcounts[i-1];
  int recvsize = displs[m_nprocs_localComm-1]
  	+ recvcounts[m_nprocs_localComm-1];

  SCT_get_elapsed_time( "MPIComm" );
  SCT_set_start( "Copie_Buffers" );

  // Buffer de doubles: cin�matique & configuration
  // Ordre par particule: [numero de particule, classe, rang exp�diteur,
  //	position, vitesse translation, quaternion rotation,
  // 	vitesse rotation]
  // -------------------------------------------------------------------------
  int NB_DOUBLE_PART = 29;
  if( AdamsBashforth ) NB_DOUBLE_PART += 12;
  if( b_hydroForce ) NB_DOUBLE_PART += 7;
  if( b_liftForce ) NB_DOUBLE_PART += 3;
  double *features = new double[NB_DOUBLE_PART*nb_hz];
  for (il=particulesHalozone->begin(),i=0; il!=particulesHalozone->end();
  	il++,i+=NB_DOUBLE_PART)
  {
    features[i] = (*il)->getID() + intTodouble;
    features[i+1] = (*il)->getParticuleClasse() + intTodouble;
    features[i+2] = m_rank + intTodouble;
    (*il)->copyVitesseTranslation( features, i+3 );
    (*il)->copyQuaternionRotation( features, i+6 );
    (*il)->copyVitesseRotation( features, i+10 );
    (*il)->copyTransform( features, i+13 );
    if( AdamsBashforth )
    {
      (*il)->copyCinematiqueNm2( features, i+29 );
      if( b_hydroForce )
        (*il)->copyFluidInformations( features, i+41 );
    }
    else if( b_hydroForce )
       (*il)->copyFluidInformations( features, i+29 );
  }

  SCT_get_elapsed_time( "Copie_Buffers" );
  SCT_set_start( "MPIComm" );

  // Partition et taille du buffer de r�ception
  int *recvcounts_DOUBLE = new int[m_nprocs_localComm];
  for (i=0; i<m_nprocs_localComm; ++i)
    recvcounts_DOUBLE[i] = NB_DOUBLE_PART * recvcounts[i];
  int *displs_DOUBLE = new int[m_nprocs_localComm];
  displs_DOUBLE[0] = 0;
  for (i=1; i<m_nprocs_localComm; ++i)
    displs_DOUBLE[i] = displs_DOUBLE[i-1] + recvcounts_DOUBLE[i-1];
  int recvsize_DOUBLE = displs_DOUBLE[m_nprocs_localComm-1]
  	+ recvcounts_DOUBLE[m_nprocs_localComm-1];
  double *recvbuf_DOUBLE = new double[recvsize_DOUBLE];

  // Communication de la cin�matique & configuration: 26 double par
  // particule
  for (ii=0;ii<m_nprocs;++ii)
    if( m_isInCommVoisins[ii] )
      MPI_Gatherv( features, NB_DOUBLE_PART * nb_hz, MPI_DOUBLE, recvbuf_DOUBLE,
  	recvcounts_DOUBLE, displs_DOUBLE, MPI_DOUBLE, m_master_localComm[ii],
	*m_commVoisins[ii] );

  SCT_add_elapsed_time( "MPIComm" );
  SCT_set_start( "UpdateCreateClones" );

  // Creation ou maj des clones
  // --------------------------
  UpdateOrCreateClones( time, recvsize, recvbuf_DOUBLE,
	NB_DOUBLE_PART, particulesClones,
	particules, particulesHalozone, ParticuleClassesReference, LC );

  delete [] features;
  delete [] recvcounts;
  delete [] recvcounts_DOUBLE;
  delete [] displs;
  delete [] displs_DOUBLE;
  delete [] recvbuf_DOUBLE;

  SCT_get_elapsed_time( "UpdateCreateClones" );
}




// ----------------------------------------------------------------------------
// Cr�ation des nouveaux clones
void MPIWrapperGrains::UpdateOrCreateClones_SendRecvLocal(Scalar time,
	list<Particule*>* particulesClones,
	list<Particule*>* particules,
  	list<Particule*> const* particulesHalozone,
	vector<Particule*> const* ParticuleClassesReference,
	LinkedCell* LC)
{
  list<Particule*>::const_iterator il;
  int i, tag_DOUBLE = 1, recvsize = 0, ireq = 0;
  MPI_Status status;
  MPI_Request sreq = 0;
  list<int> const* neighborsRank = m_voisins->rangVoisinsSeuls();
  list<int>::const_iterator irn;
  vector<MPI_Request> idreq( neighborsRank->size(), sreq );
  bool AdamsBashforth = Grains_Exec::m_TIScheme == "SecondOrderAdamsBashforth";
  bool b_hydroForce = Grains_Exec::m_withHydroForce ;
  bool b_liftForce = Grains_Exec::m_withLiftForce ;
  double intTodouble = 0.1 ;

  SCT_set_start("Copie_Buffers");

  // Remplissage de la multimap pour acces aux clones par numero
  // -----------------------------------------------------------
  AccessToClones.clear();
  for (il=particulesClones->begin();il!=particulesClones->end();il++)
    AccessToClones.insert( pair<int,Particule*>( (*il)->getID(), *il ) );

  // Copie des infos de particulesHalozone dans des buffers locaux
  // -------------------------------------------------------------
  int nb_hz = int(particulesHalozone->size());

  // Buffer de doubles: cin�matique & configuration
  // Ordre par particule: [numero de particule, classe, rang exp�diteur,
  //	position, vitesse translation, quaternion rotation,
  // 	vitesse rotation]
  // -------------------------------------------------------------------------
  int NB_DOUBLE_PART = 29;
  if( AdamsBashforth ) NB_DOUBLE_PART += 12;
  if( b_hydroForce ) NB_DOUBLE_PART += 7; // 4
  if( b_liftForce ) NB_DOUBLE_PART += 3;
  double *features = new double[NB_DOUBLE_PART*nb_hz];
  for (il=particulesHalozone->begin(),i=0;il!=particulesHalozone->end();
  	il++,i+=NB_DOUBLE_PART)
  {
    features[i] = (*il)->getID() + intTodouble;
    features[i+1] = (*il)->getParticuleClasse() + intTodouble;
    features[i+2] = m_rank + intTodouble;
    (*il)->copyVitesseTranslation( features, i+3 );
    (*il)->copyQuaternionRotation( features, i+6 );
    (*il)->copyVitesseRotation( features, i+10 );
    (*il)->copyTransform( features, i+13 );
    if ( AdamsBashforth )
    {
      (*il)->copyCinematiqueNm2( features, i+29 );
      if( b_hydroForce )
        (*il)->copyFluidInformations( features, i+41 );
    }
    else if( b_hydroForce )
       (*il)->copyFluidInformations( features, i+29 );
  }

  SCT_get_elapsed_time( "Copie_Buffers" );

  // Communication
  // -------------
  bool first_update = true;

  // Envoi par le processus � tous ses voisins
  SCT_set_start( "MPIComm" );
  for (ireq=0,irn=neighborsRank->begin();irn!=neighborsRank->end();irn++,++ireq)
    MPI_Isend( features, nb_hz * NB_DOUBLE_PART, MPI_DOUBLE,
	*irn, tag_DOUBLE, m_MPI_COMM_activProc, &idreq[ireq] );
  SCT_get_elapsed_time( "MPIComm" );

  // Reception par le processus des messages envoy�s par ses voisins
  for (irn=neighborsRank->begin();irn!=neighborsRank->end();irn++)
  {
    SCT_set_start( "MPIComm" );

    // Reception
    // ---------
    // Taille du message � recevoir -> nb de particules
    MPI_Probe( *irn, tag_DOUBLE, m_MPI_COMM_activProc, &status );
    MPI_Get_count( &status, MPI_DOUBLE, &recvsize );
    recvsize /= NB_DOUBLE_PART;

    // Reception du message de doubles
    double *recvbuf_DOUBLE = new double[recvsize * NB_DOUBLE_PART];
    MPI_Recv( recvbuf_DOUBLE, recvsize * NB_DOUBLE_PART, MPI_DOUBLE,
	*irn, tag_DOUBLE, m_MPI_COMM_activProc, &status );

    SCT_add_elapsed_time( "MPIComm" );
    SCT_set_start( "UpdateCreateClones" );

    // Creation ou maj des clones
    // --------------------------
    UpdateOrCreateClones( time, recvsize, recvbuf_DOUBLE,
		NB_DOUBLE_PART, particulesClones,
		particules, particulesHalozone, ParticuleClassesReference, LC );

    delete [] recvbuf_DOUBLE;

    if ( first_update )
    {
      SCT_get_elapsed_time( "UpdateCreateClones" );
      first_update = false;
    }
    else SCT_add_elapsed_time( "UpdateCreateClones" );
  }

  // Verifie que toutes les send non bloquants ont termine
  for (ireq=0;ireq<int(idreq.size());++ireq)
    MPI_Wait( &idreq[ireq], &status );

  delete [] features;
}




// ----------------------------------------------------------------------------
// Cr�ation des nouveaux clones
void MPIWrapperGrains::UpdateOrCreateClones_SendRecvLocal_GeoLoc(Scalar time,
	list<Particule*>* particulesClones,
	list<Particule*>* particules,
  	list<Particule*> const* particulesHalozone,
	vector<Particule*> const* ParticuleClassesReference,
	LinkedCell* LC)
{
  list<Particule*>::const_iterator il;
  int i, j, tag_DOUBLE = 1, recvsize = 0, geoLoc, ireq = 0;
  MPI_Status status;
  MPI_Request sreq = 0;
  list<int> const* neighborsRank = m_voisins->rangVoisinsSeuls();
  list<int>::const_iterator irn;
  list<MPIGeoLocalisation> const* neighborsGeoloc =
  	m_voisins->geolocVoisinsSeuls();
  list<MPIGeoLocalisation>::const_iterator ign;
  vector<MPI_Request> idreq( neighborsRank->size(), sreq );
  bool AdamsBashforth = Grains_Exec::m_TIScheme == "SecondOrderAdamsBashforth";
  bool b_hydroForce = Grains_Exec::m_withHydroForce ;
  bool b_liftForce = Grains_Exec::m_withLiftForce ;
  bool b_cohesiveForce = Grains_Exec::m_withCohesion ;
  bool b_fluidTemperature = Grains_Exec::m_withFluidTemperature ;
  bool b_solidTemperature = Grains_Exec::m_withSolidTemperature ;
  bool b_stochDrag = Grains_Exec::m_withStochasticDrag ;
  bool b_stochNu = Grains_Exec::m_withStochasticNusselt ;
  double intTodouble = 0.1 ;

  SCT_set_start( "Copie_Buffers" );

  // Remplissage de la multimap pour acces aux clones par numero
  // -----------------------------------------------------------
  AccessToClones.clear();
  for (il=particulesClones->begin();il!=particulesClones->end();il++)
    AccessToClones.insert( pair<int,Particule*>( (*il)->getID(), *il ) );


  // Copie des infos de particulesHalozone dans des buffers locaux
  // -------------------------------------------------------------
  vector<int> nbHzGeoLoc(26,0);
  vector<int>::iterator iv;
  for (il=particulesHalozone->begin();il!=particulesHalozone->end();il++)
  {
    geoLoc = (*il)->getGeoLocalisation();
    for (iv=m_particuleHalozoneToNeighboringProcs[geoLoc].begin();
    	iv!=m_particuleHalozoneToNeighboringProcs[geoLoc].end();iv++)
      nbHzGeoLoc[*iv]++;
  }

  // Buffer de doubles: cin�matique & configuration
  // Ordre par particule: [numero de particule, classe, rang exp�diteur,
  //	position, vitesse translation, quaternion rotation,
  // 	vitesse rotation]
  // -------------------------------------------------------------------------
  int NB_DOUBLE_PART=0, nAB=0, nHF=0, nLF=0, nCF=0, nT=0, nST = 0, nSTN=0;
  if( AdamsBashforth ) nAB = 12;
  if( b_hydroForce ) nHF = 7; // Eps + (Ux,Uy,Uz) + grad(Px,Py,Pz)
  if( b_liftForce ) nLF = 3; // (OMx, OMy, OMz)
  if( b_cohesiveForce ) nCF = 50; // 10*5
  if( b_solidTemperature && !b_fluidTemperature ) nT = 2;
  else if( b_fluidTemperature ) nT = 3;
  if( b_stochDrag ) nST =3;
  if( b_stochNu ) nSTN = 3;
  NB_DOUBLE_PART = 29 + nAB + nHF + nLF + nCF + nT + nST + nSTN;

  vector<int> index( 26, 0 );
  double *pDOUBLE = NULL;
  vector<double*> features( 26, pDOUBLE );
  for (i=0;i<26;i++) features[i] = new double[ NB_DOUBLE_PART * nbHzGeoLoc[i] ];
  double ParticuleID=0, ParticuleClasse=0;

  for (il=particulesHalozone->begin(),i=0;il!=particulesHalozone->end();il++)
  {
    geoLoc = (*il)->getGeoLocalisation();
    ParticuleID = (*il)->getID() + intTodouble ;
    ParticuleClasse = (*il)->getParticuleClasse() + intTodouble ;
    for (iv=m_particuleHalozoneToNeighboringProcs[geoLoc].begin();
    	iv!=m_particuleHalozoneToNeighboringProcs[geoLoc].end();iv++)
    {
      j = index[*iv];
      features[*iv][j] = ParticuleID;
      features[*iv][j+1] = ParticuleClasse;
      features[*iv][j+2] = m_rank + intTodouble ;
      (*il)->copyVitesseTranslation( features[*iv], j+3 );
      (*il)->copyQuaternionRotation( features[*iv], j+6 );
      (*il)->copyVitesseRotation( features[*iv], j+10 );
      (*il)->copyTransform( features[*iv], j+13, m_MPIperiodes[*iv] );

      if( AdamsBashforth )
        (*il)->copyCinematiqueNm2( features[*iv], j+29 );
      if( b_hydroForce )
        (*il)->copyFluidInformations( features[*iv], j+29+nAB );
      if( b_liftForce )
        (*il)->copyFluidVorticity( features[*iv], j+29+nAB+nHF );
      if( b_cohesiveForce )
      {
        (*il)->copy_VectFmaxDist( features[*iv], j+29+nAB+nHF+nLF );
        (*il)->copy_VectIdParticle( features[*iv], j+29+10+nAB+nHF+nLF );
        (*il)->copy_VectKnElast( features[*iv], j+29+20+nAB+nHF+nLF );
        (*il)->copy_VectKtElast( features[*iv], j+29+30+nAB+nHF+nLF );
        (*il)->copy_VectInitialOverlap( features[*iv], j+29+40+nAB+nHF+nLF );
      }
      if( b_solidTemperature && !b_fluidTemperature )
      {
        (*il)->copy_solidTemperature( features[*iv], j+29+nAB+nHF+nLF+nCF );
	(*il)->copy_solidNusselt( features[*iv], j+30+nAB+nHF+nLF+nCF );
//        cout << "TEMPORARY : MPI WG sending proc "<< m_rank
//             << " ID " << ParticuleID
//             << " copying  " << features[*iv][j+29+nAB+nHF+nLF+nCF] <<endl;
      }
      else if( b_fluidTemperature )
      {
        (*il)->copy_solidTemperature( features[*iv], j+29+nAB+nHF+nLF+nCF );
	(*il)->copy_solidNusselt( features[*iv], j+30+nAB+nHF+nLF+nCF );
        (*il)->copy_fluidTemperature( features[*iv], j+31+nAB+nHF+nLF+nCF );
      }
      if( b_stochDrag )
	(*il)->copy_rnd( features[*iv], j+29+nAB+nHF+nLF+nCF+nT);
      if( b_stochNu )
	(*il)->copy_rnd_Nu( features[*iv], j+29+nAB+nHF+nLF+nCF+nT+nST);

      index[*iv] += NB_DOUBLE_PART;
    }
  }
  SCT_get_elapsed_time("Copie_Buffers");

  // Communication
  // -------------
  bool first_update = true;

  // Envoi par le processus � tous ses voisins les infos liees a leur
  // geolocalisation
  SCT_set_start( "MPIComm" );
  for (ireq=0,irn=neighborsRank->begin(),ign=neighborsGeoloc->begin();
  	irn!=neighborsRank->end();irn++,ign++,++ireq)
    {
      MPI_Isend( features[*ign], nbHzGeoLoc[*ign] * NB_DOUBLE_PART, MPI_DOUBLE,
	*irn, tag_DOUBLE + m_GeoLocReciprocity[*ign],
	m_MPI_COMM_activProc, &idreq[ireq] );
      // if (nbHzGeoLoc[*ign]!=0)
      // {
      //   // Routine to store MPI messages in file
      //   ofstream MPI_log ;
      //   char filepath[]="Grains/Init/MPI_log_?_to_?.txt";
      //   filepath[20] = m_rank + '0';
      //   filepath[25] = *irn + '0';
      //   MPI_log.open(filepath, std::ios_base::app) ;
      //   for(int k=0;k<nbHzGeoLoc[*ign]*NB_DOUBLE_PART;k++)
      //   {
      //     MPI_log << std::fixed << std::setprecision(15) << features[*ign][k] ;
      //     MPI_log << ",";
      //   }
      //   MPI_log << "\n";
      // }
    }
  SCT_get_elapsed_time( "MPIComm" );

  // Reception par le processus des messages envoy�s par ses voisins
  for (irn=neighborsRank->begin(),ign=neighborsGeoloc->begin();
  	irn!=neighborsRank->end();irn++,ign++)
  {
    SCT_set_start( "MPIComm" );

    // Reception
    // ---------
    // Taille du message � recevoir -> nb de particules
    MPI_Probe( *irn, tag_DOUBLE + *ign, m_MPI_COMM_activProc, &status );
    MPI_Get_count( &status, MPI_DOUBLE, &recvsize );
    recvsize /= NB_DOUBLE_PART;

    // Reception du message de doubles
    double *recvbuf_DOUBLE = new double[recvsize * NB_DOUBLE_PART];
    MPI_Recv( recvbuf_DOUBLE, recvsize * NB_DOUBLE_PART, MPI_DOUBLE,
	*irn, tag_DOUBLE + *ign, m_MPI_COMM_activProc, &status );
    // if (recvsize != 0)
    // {
    //   // Routine to store MPI messages in file
    //   ofstream MPI_log ;
    //   char filepath[]="Grains/Init/MPI_log_?_from_?.txt";
    //   filepath[20] = m_rank + '0';
    //   filepath[27] = *irn + '0';
    //   MPI_log.open(filepath, std::ios_base::app) ;
    //   for(int k=0;k<recvsize*NB_DOUBLE_PART;k++)
    //   {
    //     MPI_log << std::fixed << std::setprecision(15) << recvbuf_DOUBLE[k] ;
    //     MPI_log << ",";
    //   }
    //   MPI_log << "\n";
    // }

    SCT_add_elapsed_time( "MPIComm" );
    SCT_set_start( "UpdateCreateClones" );

    // Creation ou maj des clones
    // --------------------------
    UpdateOrCreateClones( time, recvsize, recvbuf_DOUBLE,
		NB_DOUBLE_PART, particulesClones,
		particules, particulesHalozone, ParticuleClassesReference, LC );

    delete [] recvbuf_DOUBLE;

    if ( first_update )
    {
      SCT_get_elapsed_time( "UpdateCreateClones" );
      first_update = false;
    }
    else SCT_add_elapsed_time( "UpdateCreateClones" );
  }

  // Verifie que toutes les send non bloquants ont termine
  for (ireq=0;ireq<int(idreq.size());++ireq)
    MPI_Wait( &idreq[ireq], &status );

  for (vector<double*>::iterator ivpd=features.begin();ivpd!=features.end();
  	ivpd++) delete [] *ivpd;
  features.clear();
}




// ----------------------------------------------------------------------------
// Broadcast un double � partir du master
double MPIWrapperGrains::Broadcast_DOUBLE( const double &i ) const
{
  double collective_i = i;

  MPI_Bcast( &collective_i, 1, MPI_DOUBLE, 0, m_MPI_COMM_activProc );

  return collective_i;
}




// ----------------------------------------------------------------------------
// Broadcast un entier � partir du master
int MPIWrapperGrains::Broadcast_INT( const int &i ) const
{
  int collective_i = i;

  MPI_Bcast( &collective_i, 1, MPI_INT, 0, m_MPI_COMM_activProc );

  return collective_i;
}




// ----------------------------------------------------------------------------
// Broadcast un entier non signe � partir du master
size_t MPIWrapperGrains::Broadcast_UNSIGNED_INT( const size_t &i ) const
{
  size_t collective_i = i;

  MPI_Bcast( &collective_i, 1, MPI_UNSIGNED_LONG, 0, m_MPI_COMM_activProc );

  return collective_i;
}




// ----------------------------------------------------------------------------
// Somme d'un double sur tous les proc
double MPIWrapperGrains::sum_DOUBLE( double x ) const
{
  double sum = 0;

  MPI_Allreduce( &x, &sum, 1, MPI_DOUBLE, MPI_SUM, m_MPI_COMM_activProc );

  return sum;
}




// ----------------------------------------------------------------------------
// Somme d'un double de tous les proc sur le master
double MPIWrapperGrains::sum_DOUBLE_master( double x ) const
{
  double sum = 0;

  MPI_Reduce( &x, &sum, 1, MPI_DOUBLE, MPI_SUM, m_rank_masterWorld,
  	m_MPI_COMM_activProc );

  return sum;
}




// ----------------------------------------------------------------------------
// Somme un entier sur tous les proc
int MPIWrapperGrains::sum_INT( int i ) const
{
  int sum = 0;

  MPI_Allreduce( &i, &sum, 1, MPI_INT, MPI_SUM, m_MPI_COMM_activProc );

  return sum;
}




// ----------------------------------------------------------------------------
// Somme un entier non signe sur tous les proc
size_t MPIWrapperGrains::sum_UNSIGNED_INT( size_t i ) const
{
  size_t sum = 0;

  MPI_Allreduce( &i, &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM,
  	m_MPI_COMM_activProc );

  return sum;
}




// ----------------------------------------------------------------------------
// Somme un entier de tous les proc sur le master
int MPIWrapperGrains::sum_INT_master( int i ) const
{
  int sum = 0;

  MPI_Reduce( &i, &sum, 1, MPI_INT, MPI_SUM, m_rank_masterWorld,
  	m_MPI_COMM_activProc );

  return sum;
}




// ----------------------------------------------------------------------------
// Somme un entier de tous les proc sur le master
size_t MPIWrapperGrains::sum_UNSIGNED_INT_master( size_t i ) const
{
  size_t sum = 0;

  MPI_Reduce( &i, &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, m_rank_masterWorld,
  	m_MPI_COMM_activProc );

  return sum;
}




// ----------------------------------------------------------------------------
// Logical AND
bool MPIWrapperGrains::logical_and( bool input ) const
{
  int land=0;

  MPI_Allreduce( &input, &land, 1, MPI_UNSIGNED_SHORT, MPI_LAND,
  	m_MPI_COMM_activProc );

  return land;
}




// ----------------------------------------------------------------------------
// Min d'un entier sur tous les proc
int MPIWrapperGrains::min_INT( int i ) const
{
  int collective_i = 0;

  MPI_Allreduce( &i, &collective_i, 1, MPI_INT, MPI_MIN,
  	m_MPI_COMM_activProc ) ;

  return collective_i;
}




// ----------------------------------------------------------------------------
// Max d'un entier sur tous les proc
int MPIWrapperGrains::max_INT( int i ) const
{
  int collective_i = 0;

  MPI_Allreduce( &i, &collective_i, 1, MPI_INT, MPI_MAX, m_MPI_COMM_activProc ) ;

  return collective_i;
}




// ----------------------------------------------------------------------------
// Min d'un entier non signe sur tous les proc
size_t MPIWrapperGrains::min_UNSIGNED_INT( size_t i ) const
{
  size_t collective_i = 0;

  MPI_Allreduce( &i, &collective_i, 1, MPI_UNSIGNED_LONG, MPI_MIN,
  	m_MPI_COMM_activProc ) ;

  return collective_i;
}




// ----------------------------------------------------------------------------
// Max d'un entier non signe sur tous les proc
size_t MPIWrapperGrains::max_UNSIGNED_INT( size_t i ) const
{
  size_t collective_i = 0;

  MPI_Allreduce( &i, &collective_i, 1, MPI_UNSIGNED_LONG, MPI_MAX,
  	m_MPI_COMM_activProc ) ;

  return collective_i;
}




// ----------------------------------------------------------------------------
// Max d'un double sur tous les proc
double MPIWrapperGrains::max_DOUBLE( double x ) const
{
  double collective_x = 0.;

  MPI_Allreduce( &x, &collective_x, 1, MPI_DOUBLE, MPI_MAX,
  	m_MPI_COMM_activProc ) ;

  return collective_x;
}




// ----------------------------------------------------------------------------
// Max d'un double de tous les proc sur le master
double MPIWrapperGrains::max_DOUBLE_master( double x ) const
{
  double max = 0.;

  MPI_Reduce( &x, &max, 1, MPI_DOUBLE, MPI_MAX, m_rank_masterWorld,
  	m_MPI_COMM_activProc ) ;

  return max;
}




// ----------------------------------------------------------------------------
// Min d'un double sur tous les proc
double MPIWrapperGrains::min_DOUBLE( double x ) const
{
  double collective_x = 0.;

  MPI_Allreduce( &x, &collective_x, 1, MPI_DOUBLE, MPI_MIN,
  	m_MPI_COMM_activProc ) ;

  return collective_x;
}




// ----------------------------------------------------------------------------
// Min d'un double de tous les proc sur le master
double MPIWrapperGrains::min_DOUBLE_master( double x ) const
{
  double min = 0.;

  MPI_Reduce( &x, &min, 1, MPI_DOUBLE, MPI_MIN, m_rank_masterWorld,
  	m_MPI_COMM_activProc ) ;

  return min;
}




// ----------------------------------------------------------------------------
// AllGather d'un entier non signe
size_t* MPIWrapperGrains::AllGather_UNSIGNED_INT( size_t i ) const
{
  size_t* recv = new size_t[m_nprocs];

  MPI_Allgather( &i, 1, MPI_UNSIGNED_LONG, recv, 1, MPI_UNSIGNED_LONG,
  	m_MPI_COMM_activProc );

  return recv;
}




// ----------------------------------------------------------------------------
// Broadcast un point � partir du master
Point MPIWrapperGrains::Broadcast_Point( const Point &pt ) const
{
  double *coordinates = new double[3];
  coordinates[0] = pt[X];
  coordinates[1] = pt[Y];
  coordinates[2] = pt[Z];

  MPI_Bcast( coordinates, 3, MPI_DOUBLE, 0, m_MPI_COMM_activProc );

  Point cpt( coordinates[0],coordinates[1],coordinates[2] );
  delete [] coordinates;

  return cpt;
}




// ----------------------------------------------------------------------------
// Broadcast un vecteur � partir du master
Vecteur MPIWrapperGrains::Broadcast_Vecteur( const Vecteur &v ) const
{
  double *coordinates = new double[3];
  coordinates[0] = v[X];
  coordinates[1] = v[Y];
  coordinates[2] = v[Z];

  MPI_Bcast( coordinates, 3, MPI_DOUBLE, 0, m_MPI_COMM_activProc );

  Vecteur cv( coordinates[0],coordinates[1],coordinates[2] );
  delete [] coordinates;

  return cv;
}




// ----------------------------------------------------------------------------
// Broadcast une matrice � partir du master
Matrix MPIWrapperGrains::Broadcast_Matrix( const Matrix &mat ) const
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

  return bmat;
}




// ----------------------------------------------------------------------------
// Test avec communicateur local: validation
void MPIWrapperGrains::testCommLocal() const
{
  int nb_hz = -m_rank,i,j,ii;
  int *recvbuf_INT = new int[m_nprocs_localComm];
  for (i=0; i<m_nprocs_localComm; ++i) recvbuf_INT[i]=0;
  for (ii=0;ii<m_nprocs;++ii)
    if ( m_isInCommVoisins[ii] )
      MPI_Gather( &nb_hz, 1, MPI_INT, recvbuf_INT, 1, MPI_INT,
  	m_master_localComm[ii], *m_commVoisins[ii] );


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
// Test: AllGather d'un vecteur de n entiers
void MPIWrapperGrains::test_AllGatherv_INT( const int &n ) const
{
  int i;

  // Tableau local � envoyer
  int *local_vec = new int[n];
  for (i=0;i<n;++i) local_vec[i] = 100000 * m_rank + i;

  // Partition et taille du buffer de r�ception
  int *recvcounts_INT = new int[m_nprocs];
  for (i=0; i<m_nprocs; ++i) recvcounts_INT[i] = n;
  int *displs_INT = new int[m_nprocs];
  displs_INT[0] = 0;
  for (i=1; i<m_nprocs; ++i)
    displs_INT[i] = displs_INT[i-1] + recvcounts_INT[i-1];
  int recvsize_INT = displs_INT[m_nprocs-1] + recvcounts_INT[m_nprocs-1];
  int *recvbuf_INT = new int[recvsize_INT];

  // Communication des entiers
  MPI_Allgatherv( local_vec, n, MPI_INT, recvbuf_INT,
  	recvcounts_INT, displs_INT, MPI_INT, m_MPI_COMM_activProc );

  delete [] local_vec;
  delete [] recvcounts_INT;
  delete [] displs_INT;
  delete [] recvbuf_INT;
}




// ----------------------------------------------------------------------------
// Test avec communicateur local: AllGather d'un vecteur de n entiers
void MPIWrapperGrains::testCommLocal_AllGatherv_INT( const int &n ) const
{
  int i,ii;

  // Tableau local � envoyer
  int *local_vec = new int[n];
  for (i=0;i<n;++i) local_vec[i] = 100000 * m_rank + i;

  // Partition et taille du buffer de r�ception
  int *recvcounts_INT = new int[m_nprocs_localComm];
  for (i=0; i<m_nprocs_localComm; ++i) recvcounts_INT[i] = n;
  int *displs_INT = new int[m_nprocs_localComm];
  displs_INT[0] = 0;
  for (i=1; i<m_nprocs_localComm; ++i)
    displs_INT[i] = displs_INT[i-1] + recvcounts_INT[i-1];
  int recvsize_INT = displs_INT[m_nprocs_localComm-1]
  	+ recvcounts_INT[m_nprocs_localComm-1];
  int *recvbuf_INT = new int[recvsize_INT];

  // Communication des entiers
  for (ii=0;ii<m_nprocs;++ii)
    if ( m_isInCommVoisins[ii] )
      MPI_Gatherv( local_vec, n, MPI_INT, recvbuf_INT,
  	recvcounts_INT, displs_INT, MPI_INT, m_master_localComm[ii],
	*m_commVoisins[ii] );

  delete [] local_vec;
  delete [] recvcounts_INT;
  delete [] displs_INT;
  delete [] recvbuf_INT;
}




// ----------------------------------------------------------------------------
// Test: AllGather d'un vecteur de n doubles
void MPIWrapperGrains::test_AllGatherv_DOUBLE( const int &n ) const
{
  int i;

  // Tableau local � envoyer
  double *local_vec = new double[n];
  for (i=0;i<n;++i) local_vec[i] = 20000. * double(m_rank) + 1.1 * double(i);

  // Partition et taille du buffer de r�ception
  int *recvcounts_DOUBLE = new int[m_nprocs];
  for (i=0; i<m_nprocs; ++i) recvcounts_DOUBLE[i] = n;
  int *displs_DOUBLE = new int[m_nprocs];
  displs_DOUBLE[0] = 0;
  for (i=1; i<m_nprocs; ++i)
    displs_DOUBLE[i] = displs_DOUBLE[i-1] + recvcounts_DOUBLE[i-1];
  int recvsize_DOUBLE = displs_DOUBLE[m_nprocs-1]
  	+ recvcounts_DOUBLE[m_nprocs-1];
  double *recvbuf_DOUBLE = new double[recvsize_DOUBLE];

  // Communication des doubles
  MPI_Allgatherv( local_vec, n, MPI_DOUBLE, recvbuf_DOUBLE,
  	recvcounts_DOUBLE, displs_DOUBLE, MPI_DOUBLE, m_MPI_COMM_activProc );

  delete [] local_vec;
  delete [] recvcounts_DOUBLE;
  delete [] displs_DOUBLE;
  delete [] recvbuf_DOUBLE;
}




// ----------------------------------------------------------------------------
// Test avec communicateur local: AllGather d'un vecteur de n doubles
void MPIWrapperGrains::testCommLocal_AllGatherv_DOUBLE( const int &n ) const
{
  int i,ii;

  // Tableau local � envoyer
  double *local_vec = new double[n];
  for (i=0;i<n;++i) local_vec[i] = 20000. * double(m_rank) + 1.1 * double(i);

  // Partition et taille du buffer de r�ception
  int *recvcounts_DOUBLE = new int[m_nprocs_localComm];
  for (i=0; i<m_nprocs_localComm; ++i) recvcounts_DOUBLE[i] = n;
  int *displs_DOUBLE = new int[m_nprocs_localComm];
  displs_DOUBLE[0] = 0;
  for (i=1; i<m_nprocs_localComm; ++i)
    displs_DOUBLE[i] = displs_DOUBLE[i-1] + recvcounts_DOUBLE[i-1];
  int recvsize_DOUBLE = displs_DOUBLE[m_nprocs_localComm-1]
  	+ recvcounts_DOUBLE[m_nprocs_localComm-1];
  double *recvbuf_DOUBLE = new double[recvsize_DOUBLE];

  // Communication des doubles
  for (ii=0;ii<m_nprocs;++ii)
    if ( m_isInCommVoisins[ii] )
      MPI_Gatherv( local_vec, n, MPI_DOUBLE, recvbuf_DOUBLE,
  	recvcounts_DOUBLE, displs_DOUBLE, MPI_DOUBLE, m_master_localComm[ii],
	*m_commVoisins[ii] );

  delete [] local_vec;
  delete [] recvcounts_DOUBLE;
  delete [] displs_DOUBLE;
  delete [] recvbuf_DOUBLE;
}




// ----------------------------------------------------------------------------
// Test avec communicateur local: Send-Recv d'un vecteur de n doubles
void MPIWrapperGrains::testCommLocal_SendRecv_DOUBLE( const int &n ) const
{
  int i,ii;
  MPI_Status status;
  int nloc = n + m_rank;
//  MPI_Request sreq,rreq;

  // Tableau local � envoyer
  double *local_vec = new double[nloc];
  for (i=0;i<nloc;++i) local_vec[i] = 20000. * double(m_rank) + 1.1 * double(i);

  // Communication des doubles
  for (ii=0;ii<m_nprocs;++ii)
    if ( m_isInCommVoisins[ii] )
    {
      if ( ii != m_rank )
        MPI_Ssend( local_vec, nloc, MPI_DOUBLE, m_master_localComm[ii], 0,
		*m_commVoisins[ii] );
//         MPI_Isend( local_vec, n, MPI_DOUBLE, master_localComm[ii], 0,
// 		*commVoisins[ii], &sreq );
      else
        for (i=0;i<m_nprocs_localComm;++i)
          if ( i != m_rank_localComm )
          {
            MPI_Probe( i, 0, *m_commVoisins[ii], &status );
	    int sizemess = 0;
	    MPI_Get_count( &status, MPI_DOUBLE, &sizemess );
	    double *recvbuf_DOUBLE = new double[sizemess];
	    MPI_Recv( recvbuf_DOUBLE, sizemess, MPI_DOUBLE, i, 0,
		*m_commVoisins[ii], &status );
// 	    MPI_Irecv( recvbuf_DOUBLE, n, MPI_DOUBLE, i, 0,
// 		*commVoisins[ii], &rreq );
// 	    MPI_Wait( &rreq, &status );
            delete [] recvbuf_DOUBLE;
          }
    }

  delete [] local_vec;
}




// ----------------------------------------------------------------------------
// Bilan timer
void MPIWrapperGrains::bilanTimer() const
{
  double cputime=SCT_get_total_elapsed_time("Copie_Buffers")
     	+SCT_get_total_elapsed_time("MPIComm")
     	+SCT_get_total_elapsed_time("UpdateCreateClones");
  cout << "Bilan de la partie MPIComm" << endl;
  SCT_get_summary(cout,cputime);
}




// ----------------------------------------------------------------------------
// Complete la chaine de debug
void MPIWrapperGrains::addToMPIString( const string &add )
{
  if ( m_MPILogString ) *m_MPILogString += add;
}




// ----------------------------------------------------------------------------
// Affiche la chaine de log des comm MPI par proc et la reinitialise
void MPIWrapperGrains::writeAndFlushMPIString( ostream &f )
{
  for (int i=0;i<nombre_total_procs_ACTIV();++i)
  {
    if ( i == m_rank )
      if ( !m_MPILogString->empty() )
        f << "Processor = " << i << endl << *m_MPILogString;
    MPI_Barrier( m_MPI_COMM_activProc );
  }

  delete m_MPILogString;
  m_MPILogString = new string;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// MPI_Barrier pour les processus actifs uniquement
void MPIWrapperGrains::MPI_Barrier_ActivProc() const
{
  MPI_Barrier( m_MPI_COMM_activProc );

}




//-----------------------------------------------------------------------------
// Max d'un MaxOverlap sur tous les proc
void MPIWrapperGrains::ContactsFeatures( Scalar& overlap_max,
	Scalar& overlap_mean,
	Scalar& time_overlapMax,
	Scalar& nbIterGJK_mean ) const
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




//-----------------------------------------------------------------------------
// Cr�ation & mise � jour des clones sur la base des infos
// communiquees par les autres proc
void MPIWrapperGrains::UpdateOrCreateClones(Scalar time,
 	int const &recvsize, double const* recvbuf_DOUBLE,
	const int& NB_DOUBLE_PART,
  	list<Particule*>* particulesClones,
	list<Particule*>* particules,
  	list<Particule*> const* particulesHalozone,
	vector<Particule*> const* ParticuleClassesReference,
	LinkedCell* LC )
{
  int j, id, classe;
  bool found = false;
  bool AdamsBashforth = Grains_Exec::m_TIScheme == "SecondOrderAdamsBashforth";
  bool b_hydroForce = Grains_Exec::m_withHydroForce ;
  bool b_liftForce = Grains_Exec::m_withLiftForce ;
  bool b_cohesiveForce = Grains_Exec::m_withCohesion;
  bool b_fluidTemperature = Grains_Exec::m_withFluidTemperature ;
  bool b_solidTemperature = Grains_Exec::m_withSolidTemperature ;
  bool b_stochDrag = Grains_Exec::m_withStochasticDrag ;
  bool b_stochNu = Grains_Exec::m_withStochasticNusselt ;
  Scalar distGC = 0. ;
  Point const* GC = NULL ;
  multimap<int,Particule*>::iterator imm;
  size_t ncid = 0;
  Particule* pClone = NULL ;
  pair < multimap<int,Particule*>::iterator,
  	multimap<int,Particule*>::iterator > crange;

  int nAB=0, nHF=0, nLF=0, nCF=0, nT=0, nST=0;
  if ( AdamsBashforth ) nAB = 12;
  if ( b_hydroForce ) nHF = 7; // Eps + (Ux,Uy,Uz) + grad(Px,Py,Pz)
  if ( b_liftForce ) nLF = 3; // (OMx, OMy, OMz)
  if ( b_cohesiveForce ) nCF = 50; // 10*5
  if ( b_solidTemperature && !b_fluidTemperature ) nT = 2;
  else if ( b_fluidTemperature ) nT = 3;
  if ( b_stochDrag ) nST = 3;

  for( j=0; j<recvsize; ++j )
  {
    found = false;
    id = int( recvbuf_DOUBLE[NB_DOUBLE_PART*j] );

    // Recherche si le clone existe deja sur ce processeur
    ncid = AccessToClones.count( id );
    switch( ncid )
    {
      case 0: // pas de clone de numero id sur ce processeur
        break;
      case 1: // 1 clone de numero id sur ce processeur
        imm = AccessToClones.find( id );
        if ( m_isMPIperiodic )
        {
          GC = imm->second->getPosition();
          distGC = sqrt(
          	pow( recvbuf_DOUBLE[NB_DOUBLE_PART*j+25] - (*GC)[X], 2. ) +
            pow( recvbuf_DOUBLE[NB_DOUBLE_PART*j+26] - (*GC)[Y], 2. ) +
            pow( recvbuf_DOUBLE[NB_DOUBLE_PART*j+27] - (*GC)[Z], 2. ) ) ;
            if ( distGC < 1.1 * imm->second->getRayonInteraction() )
              found = true;
        }
        else found = true;
        break;
      default: // plus de 1 clone de numero id sur ce processeur
        // Cas multiperiodique avec 1 clone multi-processeur et 1 clone
        // periodique de meme numero sur le meme processeur
        crange = AccessToClones.equal_range( id );
        for (imm=crange.first; imm!=crange.second && !found; )
        {
          GC = imm->second->getPosition();
          distGC = sqrt(
          pow( recvbuf_DOUBLE[NB_DOUBLE_PART*j+25] - (*GC)[X], 2. ) +
          pow( recvbuf_DOUBLE[NB_DOUBLE_PART*j+26] - (*GC)[Y], 2. ) +
          pow( recvbuf_DOUBLE[NB_DOUBLE_PART*j+27] - (*GC)[Z], 2. ) ) ;
          if ( distGC < 1.1 * imm->second->getRayonInteraction() )
            found = true;
          else imm++;
        }
        break;
    }

    // Si oui => mise � jour
    if( found )
    {
      // Recuperation du pointeur sur le clone et effacement de la map
      // car un clone n'a qu'un seul et unique maitre et donc ne peut etre
      // mis a jour qu'une seule et unique fois
      pClone = imm->second;
      AccessToClones.erase( imm );

      pClone->setPosition( &recvbuf_DOUBLE[NB_DOUBLE_PART*j+13] );
      Vecteur trans( recvbuf_DOUBLE[NB_DOUBLE_PART*j+3],
	  	recvbuf_DOUBLE[NB_DOUBLE_PART*j+4],
		recvbuf_DOUBLE[NB_DOUBLE_PART*j+5] );
      pClone->setVitesseTranslation( trans );
      pClone->setQuaternionRotation( recvbuf_DOUBLE[NB_DOUBLE_PART*j+6],
		recvbuf_DOUBLE[NB_DOUBLE_PART*j+7],
		recvbuf_DOUBLE[NB_DOUBLE_PART*j+8],
		recvbuf_DOUBLE[NB_DOUBLE_PART*j+9] );
      Vecteur rot(recvbuf_DOUBLE[NB_DOUBLE_PART*j+10],
	  	recvbuf_DOUBLE[NB_DOUBLE_PART*j+11],
		recvbuf_DOUBLE[NB_DOUBLE_PART*j+12] );
      pClone->setVitesseRotation( rot );


      if( AdamsBashforth )
        pClone->setCinematiqueNm2( &recvbuf_DOUBLE[NB_DOUBLE_PART*j+29] );
      if( b_hydroForce )
      {
        pClone->set_DEMCFD_volumeFraction(
            recvbuf_DOUBLE[NB_DOUBLE_PART*j+29+nAB] );
        pClone->setVitesseTr_fluide(
            recvbuf_DOUBLE[NB_DOUBLE_PART*j+30+nAB],
            recvbuf_DOUBLE[NB_DOUBLE_PART*j+31+nAB],
            recvbuf_DOUBLE[NB_DOUBLE_PART*j+32+nAB] );
        pClone->setGradientPression_fluide(
            recvbuf_DOUBLE[NB_DOUBLE_PART*j+33+nAB],
            recvbuf_DOUBLE[NB_DOUBLE_PART*j+34+nAB],
            recvbuf_DOUBLE[NB_DOUBLE_PART*j+35+nAB] );
      }
      if( b_liftForce )
        pClone->setVorticity_fluide(
            recvbuf_DOUBLE[NB_DOUBLE_PART*j+29+nAB+nHF],
            recvbuf_DOUBLE[NB_DOUBLE_PART*j+30+nAB+nHF],
            recvbuf_DOUBLE[NB_DOUBLE_PART*j+31+nAB+nHF] );
      if( b_cohesiveForce )
      {
        for( int k=0; k<10; k++ )
        {
          pClone->set_VectFmaxDist(
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+29+k+nAB+nHF+nLF],k);
          pClone->set_VectIdParticle(
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+39+k+nAB+nHF+nLF],k);
          pClone->set_VectKnElast(
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+49+k+nAB+nHF+nLF],k);
          pClone->set_VectKtElast(
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+59+k+nAB+nHF+nLF],k);
          pClone->set_VectInitialOverlap(
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+69+k+nAB+nHF+nLF],k);
        }
      }
      if( b_solidTemperature && !b_fluidTemperature )
      {
//        cout << " TEMPORARY : MPI WG updating proc "<< m_rank
//             << " id " << id
//             << " update with " << recvbuf_DOUBLE[NB_DOUBLE_PART*j+29+nAB+nHF+nLF+nCF]
//             << endl;

        pClone->set_solidTemperature(
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+29+nAB+nHF+nLF+nCF]);
	pClone->set_solidNusselt(
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+30+nAB+nHF+nLF+nCF]);
      }
      else if( b_fluidTemperature )
      {
        pClone->set_solidTemperature(
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+29+nAB+nHF+nLF+nCF]);
	pClone->set_solidNusselt(
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+30+nAB+nHF+nLF+nCF]);
        pClone->set_DEMCFD_fluidTemperature(
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+31+nAB+nHF+nLF+nCF]);
      }
      if (b_stochDrag)
	pClone->set_rnd(recvbuf_DOUBLE[NB_DOUBLE_PART*j+29+nAB+nHF+nLF+nCF+nT],
			recvbuf_DOUBLE[NB_DOUBLE_PART*j+30+nAB+nHF+nLF+nCF+nT],
			recvbuf_DOUBLE[NB_DOUBLE_PART*j+31+nAB+nHF+nLF+nCF+nT]);
      if (b_stochNu)
	pClone->set_rnd_Nu(recvbuf_DOUBLE[NB_DOUBLE_PART*j+29+nAB+nHF+nLF+nCF+nT+nST],
                           recvbuf_DOUBLE[NB_DOUBLE_PART*j+30+nAB+nHF+nLF+nCF+nT+nST],
                           recvbuf_DOUBLE[NB_DOUBLE_PART*j+31+nAB+nHF+nLF+nCF+nT+nST]);
    }
    else // is not found
    {
      if( LC->isInLinkedCell( recvbuf_DOUBLE[NB_DOUBLE_PART*j+25],
                              recvbuf_DOUBLE[NB_DOUBLE_PART*j+26],
                              recvbuf_DOUBLE[NB_DOUBLE_PART*j+27] ) &&
          ( int( recvbuf_DOUBLE[NB_DOUBLE_PART*j+2] ) != m_rank ||
            m_isMPIperiodic ) )
      {
        classe = int( recvbuf_DOUBLE[NB_DOUBLE_PART*j+1] );

        if( Grains_Exec::m_MPI_verbose )
        {
          ostringstream oss;
          oss << "   t=" << Grains_Exec::doubleToString(time,TIMEFORMAT)
              << " Create Clone                                Id = "
              << id
              << " Classe = " << classe << " "
              << recvbuf_DOUBLE[NB_DOUBLE_PART*j+25] << " "
              << recvbuf_DOUBLE[NB_DOUBLE_PART*j+26] << " "
              << recvbuf_DOUBLE[NB_DOUBLE_PART*j+27]
              << endl;
          MPIWrapperGrains::addToMPIString(oss.str());
        }

        // Creation du clone
        Particule *new_clone = NULL ;
        if( (*ParticuleClassesReference)[classe]->isCompParticule() )
          new_clone = new CompParticule( id,
              (*ParticuleClassesReference)[classe],
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
              COMPUTE,
              2, 0 );
        else
          new_clone = new Particule( id,
              (*ParticuleClassesReference)[classe],
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
              COMPUTE,
              2 );
//        // TO BE MODIFIED, GET A POINTER ONCE, THEN PERFORM TEST ON IT !!!
//        list<App*> allApp = Grains_Exec::get_listApp();
//        list<App*>::iterator app;
//        for (app=allApp.begin(); app!=allApp.end(); app++)
//          if( (*app)->isName("TraineeHydro") && !b_hydroForce )
//            new_clone->allocateDEMCFD_FluidInfos();
        if( b_hydroForce )
          new_clone->allocateDEMCFD_FluidInfos();

        if( AdamsBashforth )
          new_clone->setCinematiqueNm2( &recvbuf_DOUBLE[NB_DOUBLE_PART*j+29] );
        if( b_hydroForce )
        {
          new_clone->set_DEMCFD_volumeFraction(
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+29+nAB] );
          new_clone->setVitesseTr_fluide(
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+30+nAB],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+31+nAB],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+32+nAB] );
          new_clone->setGradientPression_fluide(
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+33+nAB],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+34+nAB],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+35+nAB] );
        }
        if( b_liftForce )
          new_clone->setVorticity_fluide(
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+29+nAB+nHF],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+30+nAB+nHF],
              recvbuf_DOUBLE[NB_DOUBLE_PART*j+31+nAB+nHF] );
        if( b_cohesiveForce )
        {
          for( int k=0;k<10;k++ )
          {
            new_clone->set_VectFmaxDist(
                recvbuf_DOUBLE[NB_DOUBLE_PART*j+29+k+nAB+nHF+nLF],k);
            new_clone->set_VectIdParticle(
                recvbuf_DOUBLE[NB_DOUBLE_PART*j+39+k+nAB+nHF+nLF],k);
            new_clone->set_VectKnElast(
                recvbuf_DOUBLE[NB_DOUBLE_PART*j+49+k+nAB+nHF+nLF],k);
            new_clone->set_VectKtElast(
                recvbuf_DOUBLE[NB_DOUBLE_PART*j+59+k+nAB+nHF+nLF],k);
            new_clone->set_VectInitialOverlap(
                recvbuf_DOUBLE[NB_DOUBLE_PART*j+69+k+nAB+nHF+nLF],k);
          }
        }
        if( b_solidTemperature && !b_fluidTemperature )
        {
//          cout << "TEMPORARY : MPI WG creating proc " << m_rank
//               << " id " << id
//               << " initialize  "<< recvbuf_DOUBLE[NB_DOUBLE_PART*j+29+nAB+nHF+nLF+nCF] << endl;
          new_clone->set_solidTemperature(
                recvbuf_DOUBLE[NB_DOUBLE_PART*j+29+nAB+nHF+nLF+nCF]);
	  new_clone->set_solidNusselt(
                recvbuf_DOUBLE[NB_DOUBLE_PART*j+30+nAB+nHF+nLF+nCF]);
        }
        else if( b_fluidTemperature )
        {
          new_clone->set_solidTemperature(
                recvbuf_DOUBLE[NB_DOUBLE_PART*j+29+nAB+nHF+nLF+nCF]);
	  new_clone->set_solidNusselt(
                recvbuf_DOUBLE[NB_DOUBLE_PART*j+30+nAB+nHF+nLF+nCF]);
          new_clone->set_DEMCFD_fluidTemperature(
                recvbuf_DOUBLE[NB_DOUBLE_PART*j+31+nAB+nHF+nLF+nCF]);
        }
        if (b_stochDrag)
        {
	 new_clone->set_rnd(recvbuf_DOUBLE[NB_DOUBLE_PART*j+29+nAB+nHF+nLF+nCF+nT],
			recvbuf_DOUBLE[NB_DOUBLE_PART*j+30+nAB+nHF+nLF+nCF+nT],
			recvbuf_DOUBLE[NB_DOUBLE_PART*j+31+nAB+nHF+nLF+nCF+nT]);
        }
         if (b_stochNu)
        {
	 new_clone->set_rnd_Nu(recvbuf_DOUBLE[NB_DOUBLE_PART*j+29+nAB+nHF+nLF+nCF+nT+nST],
                               recvbuf_DOUBLE[NB_DOUBLE_PART*j+30+nAB+nHF+nLF+nCF+nT+nST],
                               recvbuf_DOUBLE[NB_DOUBLE_PART*j+31+nAB+nHF+nLF+nCF+nT+nST]);
        }
       // Ajout dans le LinkedCell
        LC->Link( new_clone );

        // Ajout dans les differentes listes de EnsComposants
        particulesClones->push_back(new_clone);
        particules->push_back(new_clone);

        // Ajout dans la map des clones
        AccessToClones.insert( pair<int,Particule*>( id, new_clone ) );
      }
    }
  }

  /* Mise a jour des donnees des particules elementaires */
  for(list<Particule*>::iterator il=particulesClones->begin();
  	il!=particulesClones->end();il++)
    if ( (*il)->isCompParticule() )
      (*il)->setElementPosition();
}



//-----------------------------------------------------------------------------
// Somme des efforst sur les obstacles sur le master
void MPIWrapperGrains::sumObstaclesLoad( list<MonObstacle*> const& allMyObs )
	const
{
  list<MonObstacle*>::const_iterator obstacle ;
  int nobs = int(allMyObs.size()), i = 0 ;
  Vecteur const* force = NULL;
  Vecteur const* torque = NULL;
  Vecteur collective_force, collective_torque;
  double* forcetorque = new double[6*nobs];
  double* forcetorque_collective = new double[6*nobs];

  // Copie dans un buffer local
  // sauf dans le master o� on initialise � 0
  if ( m_rank == m_rank_masterWorld )
    for (i=0;i<6*nobs;++i) forcetorque[i] = 0.;
  else
    for (obstacle=allMyObs.begin(); obstacle!=allMyObs.end(); obstacle++,i+=6)
    {
      force = (*obstacle)->getForce();
      torque = (*obstacle)->getMoment();
      forcetorque[i] = (*force)[X];
      forcetorque[i+1] = (*force)[Y];
      forcetorque[i+2] = (*force)[Z];
      forcetorque[i+3] = (*torque)[X];
      forcetorque[i+4] = (*torque)[Y];
      forcetorque[i+5] = (*torque)[Z];
    }

  // Communication du vecteur force&torque
  MPI_Reduce( forcetorque, forcetorque_collective, 6*nobs, MPI_DOUBLE,
  	MPI_SUM, m_rank_masterWorld, m_MPI_COMM_activProc );

  // Ajout sur le master des contributions des autres processeurs
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
      (*obstacle)->addMoment( collective_torque );
    }
  }

  delete [] forcetorque;
  delete [] forcetorque_collective;
}




//-----------------------------------------------------------------------------
// Ecriture de la memoire utilisee par la simulation sur chaque processeur
void MPIWrapperGrains::display_used_memory( ostream &f ) const
{
  if ( m_rank == 0 )
    f << "Memoire utilisee par Grains3D" << endl;

  for (int m=0;m<m_nprocs;++m)
  {
    if ( m == m_rank && m_is_activ )
    {
      f << "   Processeur = " << m_rank << " = ";
      Grains_Exec::display_memory( f, Grains_Exec::used_memory() );
      f << endl;
    }
    MPI_Barrier( m_MPI_COMM_activProc );
  }
}




//-----------------------------------------------------------------------------
// Repartition des nombres de particules par proc et par classe
// dans l'insertion par bloc structure
void MPIWrapperGrains::distributeParticulesClassProc(
  	list< pair<Particule*,int> > const& newPart,
	list< pair<Particule*,int> >&newPartProc,
	size_t const& npartproc,
	size_t const& ntotalinsert ) const
{
  list< pair<Particule*,int> >::const_iterator ipart;
  list< pair<Particule*,int> >::iterator ipartProc;
  int nclasseproc = 0, npartproc_ = int(npartproc), i,
  	nbClasses = int(newPart.size()) ;
  MPI_Status status;


  // Initialisation des listes & tableaux
  newPartProc = newPart;
  int* tabNbPartRestantesParClasse = new int[nbClasses];
  for (i=0,ipart=newPart.begin();i<nbClasses;++i,ipart++)
    tabNbPartRestantesParClasse[i] = ipart->second;


  // Schema saute-mouton: les actions sont faites un proc apres l'autre dans
  // l'autre croissant en debutant au 0
  if ( m_rank == 0 )
  {
    // Pour toutes les classe sauf la derniere
    ipart = newPart.begin();
    ipartProc = newPartProc.begin();
    for (i=0,ipart=newPart.begin(),ipartProc=newPartProc.begin();
    	i<nbClasses-1;++i,ipart++,ipartProc++)
    {
      nclasseproc = int( npartproc * ipart->second / ntotalinsert ) ;
      // Si le ratio est entier, on le garde tel quel
      // sinon on ajoute 1
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

    // Correction pour la derniere pour assurer que la somme de toutes les
    // classes sur ce proc est egale au nb de particules a inserer sur ce proc
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

      // Pour toutes les classe sauf la derniere
      for (i=0,ipart=newPart.begin(),ipartProc=newPartProc.begin();
    	i<nbClasses-1;++i,ipart++,ipartProc++)
      {
        nclasseproc = int( npartproc * ipart->second / ntotalinsert ) ;
        // Si le ratio est entier, on le garde tel quel
        // sinon on ajoute 1
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

      // Correction pour la derniere pour assurer que la somme de toutes les
      // classes sur ce proc est egale au nb de particules a inserer sur ce proc
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

  // Verification de la repartition par proc et par classe
  if ( m_rank == 0 )
    cout << endl << "Verification de la repartition" << endl <<
    	"Par proc et par classe" << endl;
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
	cout << "   Classe " << i << " = " << ipartProc->second << endl;
      }
      cout << "   Total: " << ntot_ << " = " << npartproc << endl;

    }
    MPI_Barrier( m_MPI_COMM_activProc );
  }
  MPI_Barrier( m_MPI_COMM_activProc );

  if ( m_rank == 0 )
    cout << "Par classe sur tous les procs" << endl;
  for (i=0,ipart=newPart.begin(),ipartProc=newPartProc.begin();
    	i<nbClasses;++i,ipart++,ipartProc++)
  {
    int ntotc =  sum_INT_master( ipartProc->second ) ;
    if ( m_rank == 0 )
      cout << "   Classe " << i << " " << ntotc << " = " << ipart->second
      	<< endl;
  }
  if ( m_rank == 0 ) cout << endl;
}
