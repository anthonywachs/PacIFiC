#include "GrainsMPITest.hh"
#include "Grains_Exec.hh"
#include "Contact_BuilderFactory.hh"
#include "Obstacle_BuilderFactory.H"
#include "ObstacleChargement.H"
#include "Voisins.hh"


//-----------------------------------------------------------------------------
// Constructeur par defaut
GrainsMPITest::GrainsMPITest() :
    GrainsMPI()
{}




//-----------------------------------------------------------------------------
// Destructeur
GrainsMPITest::~GrainsMPITest()
{}




// ----------------------------------------------------------------------------
// Gestion de la simulation
void GrainsMPITest::Simulation(bool predict, bool isPredictorCorrector,
	bool explicit_added_mass)
{
  if ( m_processorIsActiv )
  {
    m_temps = m_tdeb;
    list<App*>::iterator app;
    size_t ninsertedPart = 0;
    bool b_lastTime_save = false; 
    bool b_error_occured = false;
    bool b_perfTiming = m_rank == 0 || m_allProcTiming; 
    list<int> ClonestoDestroy, ClonestoParticules, NewPartRefPerHalozone,
  	PartRefPerOutDomainHalozone, InNotRefPerHalozone;
    string appName;
    
    if ( b_perfTiming ) 
    { 
      CT_set_start();
      SCT_insert_app( "ParticulesInsertion" );
      SCT_set_start( "ParticulesInsertion" );
      SCT_insert_app( "CalculerForces" );  
      SCT_insert_app( "Deplacer" );
      SCT_insert_app( "Actualiser" );               
      SCT_insert_app( "LinkUpdate" );
      SCT_insert_app( "InitialisationSortieResultats" );
      SCT_insert_app( "SortieResultats" );                   
    }

    // Creation, insertion et link des particules
    InsertCreateNewParticules();

    // Nombre de particules inserees et nombre total de particules dans le
    // système sur tous les procs
    ninsertedPart = m_wrapper->sum_UNSIGNED_INT( 
    	m_composants.nbreParticulesActivesOnProc() );
    m_composants.setNbreParticulesOnAllProc( ninsertedPart
  	+ m_composants.nbreParticulesWait() );
    if ( m_rank == 0 )
      cout << "Nb total de particules sur tous les procs = " << 
    	m_composants.nbreParticulesOnAllProc() << endl
	<< "Nb total de particules inserees a l'initialisation = "
	<< ninsertedPart << endl << endl;

    if ( b_perfTiming )
    { 
      SCT_get_elapsed_time( "ParticulesInsertion" );
      SCT_set_start( "InitialisationSortieResultats" );
    }

    // Initialisation des particules periodiques
    if ( Grains_Exec::m_periodique == true ) periodicParticles( b_perfTiming );

    // Par defaut l'etat initial est conserve pour Post-Processing
    m_composants.PostProcessing_start( m_temps, m_dt, m_sec, m_fenetres, 
    	m_rank, m_nprocs, m_wrapper );

    // Pour coherence entre temps initial et temps de sauvegarde
    list<Scalar>::iterator tempsSave = m_save.begin();
    while ( *tempsSave - m_temps < 0.01 * m_dt ) tempsSave++; 
    
    // Affichage temps initial
    if ( m_rank == 0 ) cout << "Temps = " << m_tdeb << "s" << endl;
    if ( b_perfTiming ) SCT_get_elapsed_time( "InitialisationSortieResultats" );

  
    // Algorithme de simulation
    // ------------------------
    while ( m_tfin - m_temps > 0.01 * m_dt ) 
    {
      m_temps += m_dt;
      b_lastTime_save = false; 


      // Initialisation de l'indicateur de calcul
      // de la transformation avec epaiseur de croute a faux 
      if ( b_perfTiming ) SCT_set_start( "CalculerForces" );
      m_composants.InitializeVdWState( m_temps, m_dt );
      if ( b_perfTiming ) SCT_get_elapsed_time( "CalculerForces" );

	
      // Insertion des particules en attente
      if ( b_perfTiming ) SCT_set_start( "ParticulesInsertion" );
      if( m_methode_insertion_particules == INLINE )
      insertParticule( m_mode_insertion_particules );
      m_sec->computeMeanNbParticules( 
      	m_composants.getParticulesActives()->size() );
      if ( b_perfTiming ) SCT_get_elapsed_time( "ParticulesInsertion" );	


      // Calcul des forces 
      if ( b_perfTiming ) SCT_set_start( "CalculerForces" );      
      // Initialisation des torseurs de force 
      m_composants.InitializeForces( m_temps, m_dt, true );
      
      // Calcul des forces exterieures : contact - fluide...      
      try {
        for (app=m_allApp.begin(); app!=m_allApp.end(); app++)
          (*app)->CalculerForces( m_temps, m_dt, 
	  	m_composants.getParticulesActives() );
	if ( b_perfTiming ) SCT_add_elapsed_time( "CalculerForces" );
      }
      catch (ErreurContact & chocCroute) 
      {
	cout << "Rank " << m_rank << " has caught exception" << endl << endl;
        chocCroute.Message(cout);
	m_composants.PostProcessingErreurComposants( "ErreurContact",
		chocCroute.getComposants() );
      } 


      // Deplacement des composants
      if ( b_perfTiming ) SCT_set_start( "Deplacer" );
      try {	
	m_composants.Deplacer( m_temps, m_dt );
      }
      catch (ErreurDeplacement &errDeplacement) 
      {
	cout << "Rank " << m_rank << " has caught exception" << endl << endl;
	errDeplacement.Message( cout );
	m_composants.PostProcessingErreurComposants( "ErreurDeplacement",
		errDeplacement.getComposant() );	
      }	 	
      if ( b_perfTiming ) SCT_get_elapsed_time( "Deplacer" );

      if ( b_perfTiming ) SCT_set_start( "Actualiser" );
      m_composants.Actualiser();
      if ( b_perfTiming ) SCT_get_elapsed_time( "Actualiser" );


      // Actualisation des particules & obstacles dans les cellules
      if ( b_perfTiming ) SCT_set_start( "LinkUpdate" );
      try {
        m_sec->LinkUpdate( m_temps, m_dt, 
		m_composants.getParticulesActives() );
      }
      catch (ErreurSimulation &errSimulation) 
      {
	cout << "Rank " << m_rank << " has caught exception" << endl << endl;
        errSimulation.Message( cout );
      }
      if ( b_perfTiming ) SCT_get_elapsed_time( "LinkUpdate" );
	
           
      // Sauvegarde des resultats
      if ( *tempsSave - m_temps < 0.01 * m_dt && tempsSave != m_save.end() ) 
      {
	// Affichage ecran du temps
	if ( b_perfTiming ) SCT_set_start( "SortieResultats" );
	if ( m_rank == 0 )
	{
          cout << "Temps = " << m_temps << "s" << endl;
	  if ( m_composants.nbreParticulesWait() )
	    cout << "Nombre de particules en attente = "
    		<< m_composants.nbreParticulesWait() << endl;
	}

        // Reload
	saveReload( m_temps );
	
	// Ecriture des résultats
	m_composants.PostProcessing( m_temps, m_dt, m_sec, m_rank, m_nprocs, 
		m_wrapper );

        // Affichage ecran des operations de communication
        if ( Grains_Exec::m_MPI_verbose ) 
	  m_wrapper->writeAndFlushMPIString( cout );
	
	// Prochain temps de sortie résultats
	while ( *tempsSave - m_temps < 0.01 * m_dt 
		&& tempsSave != m_save.end() ) tempsSave++;
	  
	b_lastTime_save = true; 
	if ( b_perfTiming ) SCT_get_elapsed_time( "SortieResultats" );
      }          
    } 


    if ( b_perfTiming ) SCT_set_start( "SortieResultats" );
    if ( m_rank == 0 ) 
    {
      cout << endl; 
      if ( m_composants.nbreParticulesWait() )
        cout << "Nombre de particules en attente = "
    	<< m_composants.nbreParticulesWait() << endl; 
    }
    
    if ( !b_lastTime_save )
    {      
      // Post processing
      if (m_rank == 0)  cout << "Ecriture des resultats pour post-processing"
    	<< " au dernier temps par defaut" << endl;
      m_composants.PostProcessing( m_temps, m_dt, m_sec, m_rank, m_nprocs, 
      	m_wrapper );
      m_composants.PostProcessing_end();	

      // Affichage ecran des operations de communication
      if ( Grains_Exec::m_MPI_verbose ) 
        m_wrapper->writeAndFlushMPIString( cout );

      // Sauvegarde du fichier de fin pour Reload
      if ( !b_error_occured )
      {
        if ( m_rank == 0 ) cout << "Sauvegarde du fichier de reload" << endl;
        saveReload( m_temps );
      } 
    }

    // Caracteristiques globales des contacts pendant la simulation
    Scalar overlap_max_ = AppSec::getOverlapMax();
    Scalar overlap_mean_ = AppSec::getOverlapMean();
    Scalar time_overlapMax_ = AppSec::getTimeOverlapMax(); 
    Scalar nbIterGJK_mean_ = AppSec::getNbIterGJKMean();     
    m_wrapper->ContactsFeatures( overlap_max_, overlap_mean_,
    	time_overlapMax_, nbIterGJK_mean_ );
    AppSec::setContactsFeatures( overlap_max_, overlap_mean_,
    	time_overlapMax_, nbIterGJK_mean_ );
    if ( m_rank == 0 ) 
    {
      cout << "Rayon d'interaction minimal = " << 
    	m_composants.getRayonInteractionMin() << endl;
      cout << "Overlap moyen = " << AppSec::getOverlapMean() << endl;
      cout << "Overlap maximal = " << AppSec::getOverlapMax() << endl;
      cout << "Temps de l'overlap maximal = " << 
	AppSec::getTimeOverlapMax() << endl;     
      cout << "Nb d'iterations moyen de GJK = " << 
    	AppSec::getNbIterGJKMean() << endl; 		
    }    
  
    if ( b_perfTiming ) SCT_get_elapsed_time("SortieResultats");

    for (int i=0;i<m_wrapper->nombre_total_procs_ACTIV();++i)
    {    
      if ( m_rank == i )
        if ( b_perfTiming )
        {
          cout << endl;
	  if ( m_allProcTiming && m_wrapper->nombre_total_procs_ACTIV() > 1 ) 
	    cout << "Processor " << m_rank << endl;
          m_wrapper->bilanTimer();
          double cputime = CT_get_elapsed_time();
          cout << endl << "Full problem" << endl;
          write_elapsed_time_smhd(cout,cputime,"Computation time"); 
	  cout << "Mean number of particles on this sub-domain = " << 
	  	AppSec::getNbParticulesPerProcMean() << endl; 	      
          SCT_get_summary( cout, cputime );
        }
      m_wrapper->MPI_Barrier_ActivProc(); 
    }   
  }
}
