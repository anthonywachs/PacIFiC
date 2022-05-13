#include "GrainsMPI.H"
#include "Grains_Exec.hh"
#include "Contact_BuilderFactory.hh"
#include "Obstacle_BuilderFactory.H"
#include "ObstacleChargement.H"
#include "Voisins.hh"
#include "PostProcessingWriter.hh"
#include "CohContact.H"


//-----------------------------------------------------------------------------
// Constructeur par defaut
GrainsMPI::GrainsMPI() :
    Grains(),
  m_wrapper( NULL ),
  m_MPIstrategie( "AllgatherGlobal" )
{}




//-----------------------------------------------------------------------------
// Destructeur
GrainsMPI::~GrainsMPI()
{
  delete m_wrapper;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Construction du probleme.
void GrainsMPI::Construction( DOMElement* rootElement )
{
  Grains::Construction( rootElement );
  if ( Grains_Exec::m_periodique == true )
    m_wrapper->setCommPeriodic( m_sec->intersectObstaclePeriodique() );
}




// ----------------------------------------------------------------------------
// Gestion de la simulation
void GrainsMPI::Simulation( bool predict, bool isPredictorCorrector,
	bool explicit_added_mass )
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
    int synchro_error_counter = 0, synchro_error_freq = 50;
    Scalar vmax = 0., vmean = 0. ;
    vector<Point> obsPos;

    if ( b_perfTiming )
    {
      CT_set_start();
      SCT_insert_app( "ParticulesInsertion" );
      SCT_set_start( "ParticulesInsertion" );
      SCT_insert_app( "CalculerForces" );
      SCT_insert_app( "ParticulesPeriodiques" );
      SCT_insert_app( "Deplacer" );
      SCT_insert_app( "Actualiser" );
      SCT_insert_app( "DestroyOutOfDomainClones" );
      SCT_insert_app( "updateInteriorsTag" );
      SCT_insert_app( "updateClonesHalozonesTag" );
      SCT_insert_app( "LinkUpdate" );
      SCT_insert_app( "CreateClones" );
      SCT_insert_app( "InitialisationSortieResultats" );
      SCT_insert_app( "SortieResultats" );
    }

    // Creation, insertion et link des particules
    InsertCreateNewParticules();

    // Nombre de particules inserees et nombre total de particules dans le
    // systeme sur tous les procs
    ninsertedPart = m_wrapper->sum_UNSIGNED_INT(
    	m_composants.nbreParticulesActivesOnProc() );
    m_composants.setNbreParticulesOnAllProc( ninsertedPart
  	+ m_composants.nbreParticulesWait() );
    if ( m_rank == 0 )
      cout << "Nb total de particules sur tous les procs = "
           << m_composants.nbreParticulesOnAllProc() << endl
           << "Nb total de particules inserees a l'initialisation = "
           << ninsertedPart << endl << endl;

//      if( (*app)->isName("TraineeHydro") && !Grains_Exec::m_withdemcfd )
//    for( app=m_allApp.begin(); app!=m_allApp.end(); app++ )
//      if( Grains_Exec::m_withHydroForce )
//      {
//        list<Particule*>* allParticles = m_composants.getParticulesActives();
//        list<Particule*>::iterator il;
//        for( il=allParticles->begin(); il!=allParticles->end(); il++ )
//          (*il)->allocateDEMCFD_FluidInfos();
//      }
    if( Grains_Exec::m_withHydroForce )
      allocateDEMCFD_FluidInfos();

    if( Grains_Exec::m_withFluidTemperature ||
        Grains_Exec::m_withSolidTemperature )
      app_FluidTemperature->InitializeTemperature( 0., m_dt,
          m_composants.getParticulesActives(), 0. );

    if ( b_perfTiming )
    {
      SCT_get_elapsed_time( "ParticulesInsertion" );
      SCT_set_start( "InitialisationSortieResultats" );
    }

    // Initialisation des particules periodiques
    if ( Grains_Exec::m_periodique == true ) periodicParticles( b_perfTiming );

    // Initialisation de la cinematique des obstacles
    m_composants.setCinematiqueObstacleSansDeplacement( m_temps, m_dt );

    // Initialise post-processing des contraintes
    ofstream fStrS, fVoidRatio;
    if ( Grains_Exec::m_stressTensor )
    {
      m_VolumeIn = 0.;
      //m_composants.setStressTensorDomain( m_constrainteFenetres );
      m_composants.ComputeVoidRatio( m_VolumeIn, obsPos, m_wrapper );
      if ( m_rank == 0 )
      {
	fVoidRatio.open( (m_fileSave + "_VoidRatio.dat").c_str(), ios::out );
	fVoidRatio << m_temps << " "
	<< Grains_Exec::doubleToString( ios::scientific, 6, m_VolumeIn );
	for ( size_t i=0; i<obsPos.size(); ++i )
	  for ( int j=0; j<3; ++j)
	    fVoidRatio << " "
	    << Grains_Exec::doubleToString( ios::scientific, 6, obsPos[i][j] );
	fVoidRatio << endl;
      }
    }

    // Par defaut l'etat initial est conserve pour Post-Processing
    m_composants.PostProcessing_start( m_temps, m_dt, m_sec, m_fenetres,
    	m_rank, m_nprocs, m_wrapper );

    m_composants.ComputeMaxMeanVelocity( vmax, vmean, m_wrapper );
    ofstream fVitMax;
    if ( m_rank == 0 )
    {
      fVitMax.open( (m_fileSave + "_VitesseMaxMean.dat").c_str(), ios::out );
      cout << "Vitesse des composants : max = " << vmax
           << " moyenne = " << vmean << endl;
      fVitMax << m_temps
              << "\t"
              << Grains_Exec::doubleToString( ios::scientific, 6, vmax )
              << "\t"
              << Grains_Exec::doubleToString( ios::scientific, 6, vmean )
              << endl;
    }
    display_used_memory();

    // Post-processing des efforts sur les obstacles
    m_composants.initialiseOutputObstaclesLoadFiles( m_rank, false, m_temps );
    m_composants.outputObstaclesLoad( m_temps, m_dt, false,
      	Grains_Exec::m_ReloadType == "same", m_rank, m_nprocs,
    	m_wrapper );

    // Pour coherence entre temps initial et temps de sauvegarde
    list<Scalar>::iterator tempsSave = m_save.begin();
    while ( *tempsSave - m_temps < 0.01 * m_dt ) tempsSave++;

    // Affichage temps initial
    if ( m_rank == 0 ) cout << "Temps = " << m_tdeb << "s" << endl;
    if ( b_perfTiming ) SCT_get_elapsed_time( "InitialisationSortieResultats" );

    // If cohesive contact, check for particle initially glued
    if (Grains_Exec::m_withCohesion)
    {
      m_sec->InitializeCohesiveForces( 0., m_dt,
          m_composants.getParticulesActives() );
    }

    // Algorithme de simulation
    // ------------------------
    while ( m_tfin - m_temps > 0.01 * m_dt )
    {
      m_temps += m_dt;
      b_lastTime_save = false;
      ++synchro_error_counter;

      // Initialisation de l'indicateur de calcul
      // de la transformation avec epaiseur de croute a faux
      if ( b_perfTiming ) SCT_set_start( "CalculerForces" );
      if ( !b_error_occured ) m_composants.InitializeVdWState( m_temps, m_dt );
      if ( b_perfTiming ) SCT_get_elapsed_time( "CalculerForces" );


      // Insertion des particules en attente
      if( b_perfTiming ) SCT_set_start( "ParticulesInsertion" );
      if( m_methode_insertion_particules == INLINE && !b_error_occured )
        insertParticule( m_mode_insertion_particules );
      if( !b_error_occured ) m_sec->computeMeanNbParticules(
          m_composants.getParticulesActives()->size() );
      if( b_perfTiming ) SCT_get_elapsed_time( "ParticulesInsertion" );

      // Initialisation des maps de contact
      m_composants.setAllContactMapToFalse();

      // Calcul des forces
      if ( b_perfTiming ) SCT_set_start( "CalculerForces" );
      // Initialisation des torseurs de force
      if ( !b_error_occured )
        m_composants.InitializeForces( m_temps, m_dt, true );

      // Calcul des forces exterieures : contact - fluide...
      try {
        if ( !b_error_occured )
          for (app=m_allApp.begin(); app!=m_allApp.end(); app++)
            (*app)->CalculerForces( m_temps, m_dt,
                m_composants.getParticulesActives() );
        if ( b_perfTiming ) SCT_add_elapsed_time( "CalculerForces" );

        // Calcul des forces de contact sur les clones periodiques
        if ( b_perfTiming ) SCT_set_start( "ParticulesPeriodiques" );
        if ( Grains_Exec::m_periodique == true && !b_error_occured )
          m_sec->CalculerForcesClonesPeriodiques( m_temps, m_dt,
        m_composants.getParticulesClonesPeriodiques() );
        if ( b_perfTiming ) SCT_get_elapsed_time( "ParticulesPeriodiques" );
      }
      catch (ErreurContact & chocCroute)
      {
        cout << "Rank " << m_rank << " has caught exception" << endl << endl;
        chocCroute.Message(cout);
        m_composants.PostProcessingErreurComposants( "ErreurContact",
            chocCroute.getComposants() );
        b_error_occured = true;
      }

      // Calcul de la temperature en fonction des flux de chaleur
//      if( app_FluidTemperature ) // change name to app_Temperature or use bool
      if( Grains_Exec::m_withFluidTemperature ||
          Grains_Exec::m_withSolidTemperature )
      {
          m_composants.ComputeTemperature( m_temps, m_dt );
      }

      // Mise a jour des maps de contact
      m_composants.updateAllContactMaps();

      // Deplacement des composants
      if ( b_perfTiming ) SCT_set_start( "Deplacer" );
      if ( Grains_Exec::m_periodique == true && !b_error_occured )
        m_wrapper->commParticulesClonesPer_AllGatherGlobal_Forces( m_temps,
  		m_composants.getParticulesReferencesPeriodiques(),
  		m_composants.getParticulesClonesPeriodiques() );

      if ( !b_error_occured )
        try {
          m_composants.Deplacer( m_temps, m_dt );
        }
        catch (ErreurDeplacement &errDeplacement)
        {
          cout << "Rank " << m_rank << " has caught exception" << endl << endl;
          errDeplacement.Message(cout);
          m_composants.PostProcessingErreurComposants("ErreurDeplacement",
            errDeplacement.getComposant());
          b_error_occured = true;
        }

      if ( b_perfTiming ) SCT_get_elapsed_time( "Deplacer" );

      if ( b_perfTiming ) SCT_set_start( "Actualiser" );
      if ( !b_error_occured ) m_composants.Actualiser();
      if ( b_perfTiming ) SCT_get_elapsed_time( "Actualiser" );

      // Mise � jour du tag des particules int�rieures
      if ( b_perfTiming ) SCT_set_start( "updateInteriorsTag" );

      if ( !b_error_occured ) m_sec->updateInteriorTag( m_temps,
            m_composants.getParticulesActives(),
            m_composants.getParticulesHalozone(),
            m_wrapper );

      if ( b_perfTiming ) SCT_get_elapsed_time( "updateInteriorsTag" );

      // Cr�ation ou Mise � jour des clones
      if ( b_perfTiming ) SCT_set_start( "CreateClones" );
      if ( m_MPIstrategie == "AllgatherGlobal" )
        m_wrapper->UpdateOrCreateClones_AllGatherGlobal( m_temps,
            m_composants.getParticulesClones(),
            m_composants.getParticulesActives(),
            m_composants.getParticulesHalozone(),
            m_composants.getParticuleClassesReference(),
            m_sec );
      else if ( m_MPIstrategie == "AllgatherLocal" )
        m_wrapper->UpdateOrCreateClones_AllGatherLocal( m_temps,
            m_composants.getParticulesClones(),
            m_composants.getParticulesActives(),
            m_composants.getParticulesHalozone(),
            m_composants.getParticuleClassesReference(),
            m_sec );
      else if ( m_MPIstrategie == "SRLocal" )
        m_wrapper->UpdateOrCreateClones_SendRecvLocal( m_temps,
            m_composants.getParticulesClones(),
            m_composants.getParticulesActives(),
            m_composants.getParticulesHalozone(),
            m_composants.getParticuleClassesReference(),
            m_sec );
      else if ( m_MPIstrategie == "SRLocalCommOpt" )
        m_wrapper->UpdateOrCreateClones_SendRecvLocal_GeoLoc( m_temps,
            m_composants.getParticulesClones(),
            m_composants.getParticulesActives(),
            m_composants.getParticulesHalozone(),
            m_composants.getParticuleClassesReference(),
            m_sec );
      if ( b_perfTiming ) SCT_get_elapsed_time( "CreateClones" );

      // Suppresion des clones hors du Linked Cell (car passes sur un autre
      // processeur)
      if ( b_perfTiming ) SCT_set_start( "DestroyOutOfDomainClones" );
      if ( !b_error_occured ) m_sec->DestroyOutOfDomainClones( m_temps,
		m_composants.getParticulesClones(),
		m_composants.getParticulesReferencesPeriodiques(),
		m_composants.getParticulesActives() );
      if ( b_perfTiming ) SCT_get_elapsed_time( "DestroyOutOfDomainClones" );


      // Mise � jour du tag des particules clones et halozone
      if ( b_perfTiming ) SCT_set_start( "updateClonesHalozonesTag" );
      if ( !b_error_occured ) m_sec->updateHalozoneCloneTag( m_temps,
		m_composants.getParticulesHalozone(),
		m_composants.getParticulesClones(),
		m_wrapper );
      if ( b_perfTiming ) SCT_get_elapsed_time( "updateClonesHalozonesTag" );


      // Actualisation des particules & obstacles dans les cellules
      if ( b_perfTiming ) SCT_set_start( "LinkUpdate" );
      if ( !b_error_occured )
        try {
          m_sec->LinkUpdate( m_temps, m_dt,
	  	m_composants.getParticulesActives() );
        }
        catch (ErreurSimulation &errSimulation)
        {
	  cout << "Rank " << m_rank << " has caught exception" << endl << endl;
          errSimulation.Message(cout);
	  b_error_occured = true;
        }
      if ( !b_error_occured )
        m_composants.updateGeoLocalisationParticulesHalozone();
      if ( b_perfTiming ) SCT_get_elapsed_time( "LinkUpdate" );


      // Synchronization des erreurs potentielles
      if ( synchro_error_counter == synchro_error_freq )
      {
        if ( b_perfTiming ) SCT_set_start( "CreateClones" );
        if ( m_wrapper->max_INT( Grains_Exec::m_exception_Contact +
      		Grains_Exec::m_exception_Deplacement +
            Grains_Exec::m_exception_Simulation ) )
          { b_error_occured = true; break; }
        if ( b_perfTiming ) SCT_add_elapsed_time( "CreateClones" );
          synchro_error_counter = 0;
      }


      // Traitement des particules periodiques
      if ( Grains_Exec::m_periodique == true && !b_error_occured )
        periodicParticles( b_perfTiming );


      // Post-processing des efforts sur les obstacles
      m_composants.outputObstaclesLoad( m_temps, m_dt, false, false,
      	m_rank, m_nprocs, m_wrapper );

      // Sauvegarde des resultats
      if ( *tempsSave - m_temps < 0.01 * m_dt && tempsSave != m_save.end() )
      {
	// Affichage ecran du temps
	if ( b_perfTiming ) SCT_set_start( "SortieResultats" );
	m_composants.ComputeMaxMeanVelocity( vmax, vmean, m_wrapper );

	// Post-processing des contraintes
	if ( Grains_Exec::m_stressTensor )
	{
	  vector<Scalar> stressTensor(9,0.);
	  //m_composants.setStressTensorDomain( m_constrainteFenetres );
	  //m_composants.computeStressTensor( stressTensor, m_wrapper );
	  m_composants.ComputeVoidRatio( m_VolumeIn, obsPos, m_wrapper );
	  if ( m_rank == 0 )
	  {
	    fVoidRatio << m_temps << "\t" <<
	      Grains_Exec::doubleToString( ios::scientific, 6, m_VolumeIn );
	    for ( size_t i=0; i<obsPos.size(); ++i )
	      for ( int j=0; j<3; ++j)
	        fVoidRatio << " "
		<< Grains_Exec::doubleToString( ios::scientific, 6, obsPos[i][j] );
	    fVoidRatio << endl;
	  }
	}

        if ( m_rank == 0 )
        {
          fVitMax << Grains_Exec::doubleToString( ios::scientific, 6, m_temps )
    		<< "\t" << Grains_Exec::doubleToString( ios::scientific, 6,
		vmax ) << "\t" << Grains_Exec::doubleToString( ios::scientific,
		6, vmean ) << endl;
          cout << "Temps = " << m_temps << "s" << endl;
	  cout << "Vitesse des composants : max = " << vmax << " moyenne = " <<
    		vmean << endl;
	  if ( m_composants.nbreParticulesWait() )
	    cout << "Nombre de particules en attente = "
    		<< m_composants.nbreParticulesWait() << endl;
	}

        // Utilisation memoire
        display_used_memory();

	// Synchronization des erreurs potentielles avant sauvegarde resultats
        if ( b_perfTiming ) SCT_set_start( "CreateClones" );
	if ( m_wrapper->max_INT( Grains_Exec::m_exception_Contact +
      		Grains_Exec::m_exception_Deplacement +
		Grains_Exec::m_exception_Simulation ) )
	  b_error_occured = true;
	if ( b_perfTiming ) SCT_add_elapsed_time( "CreateClones" );

        // Reload
	if ( !b_error_occured ) saveReload( m_temps );

	// Ecriture des r�sultats
	m_composants.PostProcessing( m_temps, m_dt,
		b_error_occured ? NULL : m_sec, m_rank, m_nprocs, m_wrapper );

        // Affichage ecran des operations de communication
        if ( Grains_Exec::m_MPI_verbose )
	  m_wrapper->writeAndFlushMPIString( cout );

	// Prochain temps de sortie r�sultats
	while ( *tempsSave - m_temps < 0.01 * m_dt
		&& tempsSave != m_save.end() ) tempsSave++;

	b_lastTime_save = true;
	if ( b_perfTiming ) SCT_get_elapsed_time( "SortieResultats" );
      } // end if tempsSave
    } // end if t < EndTime


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
      // Vitesses max et moyenne
      m_composants.ComputeMaxMeanVelocity( vmax, vmean, m_wrapper );
      if ( m_rank == 0 )
      {
        fVitMax << Grains_Exec::doubleToString( ios::scientific, 6, m_temps )
    	<< "\t" << Grains_Exec::doubleToString( ios::scientific, 6, vmax )
	<< "\t" << Grains_Exec::doubleToString( ios::scientific, 6, vmean )
	<< endl;
	cout << "Vitesse des composants : max = " << vmax << " moyenne = " <<
    		vmean << endl;
      }

      // Post processing
      if ( m_rank == 0 )  cout << "Ecriture des resultats pour post-processing"
    	<< " au dernier temps par defaut" << endl;
      m_composants.PostProcessing( m_temps, m_dt,
      	b_error_occured ? NULL : m_sec, m_rank, m_nprocs, m_wrapper );
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
    if ( m_rank == 0 ) fVitMax.close();
    if ( Grains_Exec::m_stressTensor )
      if ( m_rank == 0 )
      {
	fStrS.close();
	fVoidRatio.close();
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

    if ( b_perfTiming ) SCT_get_elapsed_time( "SortieResultats" );

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
          write_elapsed_time_smhd(cout,cputime,"Computation time" );
	  cout << "Mean number of particles on this sub-domain = " <<
	  	AppSec::getNbParticulesPerProcMean() << endl;
          SCT_get_summary(cout,cputime);
        }
      m_wrapper->MPI_Barrier_ActivProc();
    }
  }
}




// ----------------------------------------------------------------------------
// Cree, insert et link les nouvelles particules
void GrainsMPI::InsertCreateNewParticules()
{
  // IMPORTANT: quel que soit le systeme etudie, si celui ci contient
  // N particules et M obstacles, les particules sont numerotees de 0 � N-1
  // et les obstacles de N � N+M-1
  // Dans le cas de reload avec insertion supplementaire, cela necessite de
  // renumeroter les obstacles

  int numPartMax = numeroMaxParticules();
  numPartMax = m_wrapper->max_INT( numPartMax );
  if (numPartMax) ++numPartMax;
  list< pair<Particule*,int> >::iterator ipart;

  // Construction des nouvelles particules
  Composant::setNbComposantsCrees( numPartMax );

  // Si position par defaut, on positionne les particules
  if ( m_position != "" )
  {
    // Positions dans un bloc structure
    if ( m_position == "STRUCTURED" )
      setPositionParticulesBloc( m_mode_insertion_particules );
    // Positions dans un fichier
    else setPositionParticulesFichier( m_mode_insertion_particules );

    m_wrapper->MPI_Barrier_ActivProc();
    if ( m_rank == 0 )
      cout << "Positions affectees aux particules" << endl;
  }
  else
  {
    for (ipart=m_newParticules.begin();ipart!=m_newParticules.end();ipart++)
    {
      int nbre = ipart->second;
      for (int ii=0; ii<nbre; ii++)
      {
	  Particule* particule = ipart->first->createCloneCopy();
	  m_composants.Ajouter( particule );
      }
    }
  }

  // Renumerotation des obstacles
  int numInitObstacles = numPartMax;
  for (ipart=m_newParticules.begin();ipart!=m_newParticules.end();ipart++)
    numInitObstacles += ipart->second;
  list<MonObstacle*> lisObstaclesPrimaires =
    	m_composants.getObstacles()->getObstacles();
  list<MonObstacle*>::iterator iobs;
  int ObstacleID = numInitObstacles;
  for (iobs=lisObstaclesPrimaires.begin();iobs!=lisObstaclesPrimaires.end();
    	iobs++)
  {
    (*iobs)->setID( ObstacleID );
    ++ObstacleID;
  }
  Composant::setNbComposantsCrees( ObstacleID );

  // Affectation des composants aux forces exterieures
  m_composants.Link( *m_sec );

  // Insertion des particules dont la position est donnee
  size_t nbPW = 0 ;
  bool b_insertion_BEFORE = true;
  if ( m_methode_insertion_particules == BEFORE )
  {
    nbPW = m_composants.nbreParticulesWait() ;
    for (size_t i=0;i<nbPW && b_insertion_BEFORE;++i)
      b_insertion_BEFORE = insertParticule( m_mode_insertion_particules );
    if ( !b_insertion_BEFORE )
    {
      cout << "Processeur " << m_rank << endl;
      cout << "Prob d'insertion par methode positions" << endl;
      cout << "Nb de particules restantes a inserer = " <<
      	m_composants.nbreParticulesWait() << endl;
    }

    // il faut a present creer les clones correspondants qui ne sont pas pris
    // en charge par les methodes setPositionParticulesBloc et
    // setPositionParticulesFichier
    if ( m_position != "" )
    {
      if ( m_MPIstrategie == "AllgatherGlobal" )
        m_wrapper->UpdateOrCreateClones_AllGatherGlobal( m_temps,
      		m_composants.getParticulesClones(),
		m_composants.getParticulesActives(),
  		m_composants.getParticulesHalozone(),
		m_composants.getParticuleClassesReference(),
		m_sec );
      else if ( m_MPIstrategie == "AllgatherLocal" )
        m_wrapper->UpdateOrCreateClones_AllGatherLocal( m_temps,
      		m_composants.getParticulesClones(),
		m_composants.getParticulesActives(),
  		m_composants.getParticulesHalozone(),
		m_composants.getParticuleClassesReference(),
		m_sec );
      else if ( m_MPIstrategie == "SRLocal" )
        m_wrapper->UpdateOrCreateClones_SendRecvLocal( m_temps,
      		m_composants.getParticulesClones(),
		m_composants.getParticulesActives(),
  		m_composants.getParticulesHalozone(),
		m_composants.getParticuleClassesReference(),
		m_sec );
      else if ( m_MPIstrategie == "SRLocalCommOpt" )
        m_wrapper->UpdateOrCreateClones_SendRecvLocal_GeoLoc( m_temps,
      		m_composants.getParticulesClones(),
		m_composants.getParticulesActives(),
  		m_composants.getParticulesHalozone(),
		m_composants.getParticuleClassesReference(),
		m_sec );
    }
  }

  // Si probleme d'insertion par methode BEFORE, arret du calcul
  if ( !b_insertion_BEFORE ) grainsAbort();

  double volIN = m_composants.getVolumeIn(),
  	volOUT = m_composants.getVolumeOut();
  volIN = m_wrapper->sum_DOUBLE( volIN );
  volOUT = m_wrapper->sum_DOUBLE( volOUT );

  if ( m_rank == 0 )
    cout << endl << "Volume des particules IN  : " << volIN  << '\n'
      	<< "                      OUT : " << volOUT << '\n'
      	<< endl;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Insertion d'une particule dans les algorithmes
bool GrainsMPI::insertParticule( const PullMode& mode )
{
  Particule *particule = NULL;
  list<App*>::iterator app;
  static int insert_counter = 0;
  bool collective_contact = true, inLinkCell = false;
  Point position;
  Transform trot;
  Matrix mrot;
  Vecteur vtrans, vrot;

  if ( insert_counter == 0 )
  {
    // En parall�le, si on insert des particules avec des positions predefinies
    // (fichier ou Bloc), les particules sont crees et inserees "a la volee"
    // Donc chaque processeur a une liste differente de particules en attente
    // Ainsi, le tirage aleatoire est realise dans la methode
    // setPositionParticulesFichier, donc a ce niveau on force a parcourir la
    // liste des particules en attente de maniere sequentielle (on parcourt
    // toute la liste de toute facon)
    //
    // Si l'insertion est de type tir-fenetre, l'operation est similaire au
    // sequentiel, tous les processeurs ont la meme liste de particules en
    // attente et le tirage aleatoire est realise dans composants.getParticule
    // si le mode est RANDOM
    if ( m_position != "" )
      particule = m_composants.getParticule( ORDER, m_wrapper );
    else particule = m_composants.getParticule( mode, m_wrapper );

    if ( particule )
    {
      bool contact = false;

      // Cas d'une insertion tir-fenetre
      if ( m_position == "" )
      {
        // Afin que le tirage al�atoire soit le m�me sur tous les procs
        // seul le master tire au hasard et envoie les infos aux autres procs

        // Initialisation de la vitesse de la particule
        if ( m_rank == 0 ) computeInitVit( vtrans, vrot );
        vtrans = m_wrapper->Broadcast_Vecteur( vtrans );
        vrot = m_wrapper->Broadcast_Vecteur( vrot );
        particule->setVitesseTranslation( vtrans );
        particule->setVitesseRotation( vrot );

        // Initialisation de la position du centre de gravite de la particule
        if ( m_rank == 0 ) position = getPoint();
        position = m_wrapper->Broadcast_Point( position );
        particule->setPosition( position );

        // Initialisation de la position angulaire de la particule
        if ( m_configAleatoire )
        {
          if ( m_rank == 0 ) mrot = getRandomRotation();
          mrot = m_wrapper->Broadcast_Matrix( mrot );
          trot.setBasis( mrot );
          particule->composePosition( trot );
        }
      }
      // Cas d'une insertion par positions predefinies
      else
      {
        // Initialisation de la vitesse de la particule
        computeInitVit( vtrans, vrot );
        particule->setVitesseTranslation( vtrans );
        particule->setVitesseRotation( vrot );

        // Initialisation de la position angulaire de la particule
        if ( m_configAleatoire )
        {
          trot.setBasis( getRandomRotation() );
          particule->composePosition( trot );
        }
      }

      particule->initializeVdWtransform_to_notComputed();

      if ( m_sec->isInLinkedCell( *particule->getPosition() ) )
      {
        inLinkCell = true;
        // Recherche de contact avec les autres particules & obstacles
        if ( !m_force_insertion )
          contact = m_sec->isContactVdW( particule );
        else contact = false;
      }

      // Bilan sur tous les procs
      if ( !m_force_insertion )
        collective_contact = m_wrapper->max_INT( contact );
      else collective_contact = 0;

      // Si pas de contact => insertion sur le proc sur lequel la particule est
      // situee et �limination sur les autres
      // sinon pas d'insertion et une nouvelle position (aleatoire) sera
      // tiree au prochain appel
      if ( !collective_contact )
      {
        if ( inLinkCell )
        {
          m_sec->Link( particule );
          m_composants.ShiftParticuleOutIn();

          // Verifie l'allocation de la structure des infos du fluide pour
          // le calcul de la force de trainee avec le fluide au repos
          // (fluide non resolu)
          if ( Grains_Exec::m_withHydroForce )
            particule->allocateDEMCFD_FluidInfos();
        }
        else
          m_composants.DeleteAndDestroyWait();
      }
    }
  }

  ++insert_counter;
  if ( insert_counter == m_insertion_frequency ) insert_counter = 0;

  return( inLinkCell && ( !collective_contact || m_force_insertion ) );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Positionne les particules en attente pour la simulation a partir
// d'un fichier de positions
void GrainsMPI::setPositionParticulesFichier( const PullMode& mode )
{
  // Verification de l'existence du fichier
  ifstream filePos;
  size_t notok = 0;
  if ( m_rank == 0 )
  {
    filePos.open( m_position.c_str() );
    if ( !filePos.is_open() ) notok = 1;
  }
  notok = m_wrapper->sum_UNSIGNED_INT( notok );
  if ( notok )
  {
    if ( m_rank == 0 ) cout << "ERR : Fichier absent " << m_position
    	<< endl << endl;
    m_wrapper->MPI_Barrier_ActivProc();
    grainsAbort();
  }

  Particule* newPart = NULL, *particuleClasse = NULL;
  int id;
  Point position;
  list< pair<Particule*,int> > newParticules_ = m_newParticules;
  list< pair<Particule*,int> >::iterator ipart;

  // Verifie que le nb de positions est egal au nombre de particules a inserer
  size_t nwait = 0, npos = 0;
  if ( m_rank == 0 )
  {
    for (ipart=m_newParticules.begin();ipart!=m_newParticules.end();ipart++)
      nwait += ipart->second;
    while ( filePos >> position ) ++npos;
    if ( nwait != npos ) notok = 1;
    filePos.close();
  }
  notok = m_wrapper->sum_UNSIGNED_INT( notok );
  if ( notok )
  {
    if ( m_rank == 0 )
      cout << "ERR: nombre de particules a inserer est "
           << "different du nombre de positions dans le fichier " << m_position
           << " : " << nwait << " != " << npos << endl << endl;
    m_wrapper->MPI_Barrier_ActivProc();
    grainsAbort();
  }


  // Initialisation de la numerotation des particules
  int numPartMax = numeroMaxParticules();
  numPartMax = m_wrapper->max_INT( numPartMax );
  if ( numPartMax ) ++numPartMax;
  id = numPartMax;


  // Lit la position des particules et cree une particule quand la position
  // est dans le domaine de calcul (sans les zones de halo du LinkedCell)
  filePos.open( m_position.c_str() );
  while ( filePos >> position )
  {
    // Selection d'une classe de particules
    particuleClasse = getParticuleClasseForCreation( mode, newParticules_,
    	false );

    // Si la particule est sur ce processeur, on la creee pour l'inserer
    // ulterieurement
    if ( App::isInLocalDomain( &position ) )
    {
      newPart = particuleClasse->createCloneCopy();
      // Previous implementation: not compatible with the compParticule
      // newPart = new Particule( *particuleClasse );
      newPart->setPosition( position );
      newPart->setID( id );
      m_composants.Ajouter( newPart );
    }
    id++;
  }
  filePos.close();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Positionne les particules en attente pour la simulation a partir
// d'un bloc structure de positions
void GrainsMPI::setPositionParticulesBloc( const PullMode& mode )
{
  Particule* newPart = NULL, *particuleClasse = NULL;
  int id, k, l, m, kmin = 0, kmax = int(m_blocInsert->NX) - 1,
  	lmin = 0, lmax = int(m_blocInsert->NY) - 1,
	mmin = 0, mmax = int(m_blocInsert->NZ) - 1,
	npartproc;
  Point position;
  list< pair<Particule*,int> > newParticules_ = m_newParticules;
  bool found = false, no_overlap = false ;

  // Test du nb de particules a inserer
  size_t ntotalinsert = 0;
  list< pair<Particule*,int> >::iterator ipart,ipart2;
  for (ipart=m_newParticules.begin();ipart!=m_newParticules.end();ipart++)
    ntotalinsert += ipart->second;

  if ( ntotalinsert != m_blocInsert->NX * m_blocInsert->NY * m_blocInsert->NZ )
  {
    if ( m_rank == 0 )
      cout << "ERR: nombre de particules a inserer est different du"
    	<< " nombre de positions du bloc structure: " <<
	ntotalinsert << " != " <<
	m_blocInsert->NX * m_blocInsert->NY * m_blocInsert->NZ << endl;
    m_wrapper->MPI_Barrier_ActivProc();
    grainsAbort();
  }


  // Taille de maille du bloc structure entre 2 particles
  double deltax = ( m_blocInsert->bloc.ptB[X] - m_blocInsert->bloc.ptA[X] )
  	/ double(m_blocInsert->NX) ;
  double deltay = ( m_blocInsert->bloc.ptB[Y] - m_blocInsert->bloc.ptA[Y] )
  	/ double(m_blocInsert->NY) ;
  double deltaz = ( m_blocInsert->bloc.ptB[Z] - m_blocInsert->bloc.ptA[Z] )
  	/ double(m_blocInsert->NZ) ;


  // Coordonnees du domaine local
  Point origineLocale;
  App::getOrigineLocale( origineLocale[X], origineLocale[Y],
  	origineLocale[Z] );
  Vecteur dimensionsLocales;
  App::getDimensionsLocales( dimensionsLocales[X], dimensionsLocales[Y],
  	dimensionsLocales[Z] );
  Point MaxLocal = origineLocale + dimensionsLocales;


  // Test qu'aucune coordonnee a inserer n'est exactement egale aux limites
  // du domaine locale
  // si oui translation de 1e-14
  vector<size_t> coorMatchLocLim( 3, 0 );
  for (k=0;k<int(m_blocInsert->NX) && !coorMatchLocLim[0];++k)
  {
    position[X] = m_blocInsert->bloc.ptA[X] + ( double(k) + 0.5 ) * deltax;
    if ( fabs( position[X] - origineLocale[X] ) < 1.e-14
    	|| fabs( position[X] - MaxLocal[X] ) < 1.e-14 )
      coorMatchLocLim[0] = 1;
  }

  for (l=0;l<int(m_blocInsert->NY) && !coorMatchLocLim[1];++l)
  {
    position[Y] = m_blocInsert->bloc.ptA[Y] + ( double(l) + 0.5 ) * deltay;
    if ( fabs( position[Y] - origineLocale[Y] ) < 1.e-14
    	|| fabs( position[Y] - MaxLocal[Y] ) < 1.e-14 )
      coorMatchLocLim[1] = 1;
  }

  for (m=0;m<int(m_blocInsert->NZ) && !coorMatchLocLim[2];++m)
  {
    position[Z] = m_blocInsert->bloc.ptA[Z] + ( double(m) + 0.5 ) * deltaz;
    if ( fabs( position[Z] - origineLocale[Z] ) < 1.e-14
    	|| fabs( position[Z] - MaxLocal[Z] ) < 1.e-14 )
      coorMatchLocLim[2] = 1;
  }

  coorMatchLocLim[0] = m_wrapper->max_UNSIGNED_INT( coorMatchLocLim[0] );
  if ( coorMatchLocLim[0] ) m_blocInsert->bloc.ptA[X] += 1.e-14;
  coorMatchLocLim[1] = m_wrapper->max_UNSIGNED_INT( coorMatchLocLim[1] );
  if ( coorMatchLocLim[1] ) m_blocInsert->bloc.ptA[Y] += 1.e-14;
  coorMatchLocLim[2] = m_wrapper->max_UNSIGNED_INT( coorMatchLocLim[2] );
  if ( coorMatchLocLim[2] ) m_blocInsert->bloc.ptA[Z] += 1.e-14;

  if ( ( coorMatchLocLim[0] || coorMatchLocLim[1] || coorMatchLocLim[2] )
  	&& m_rank == 0 )
  {
    cout << endl << "Warning: Position en bloc structure: certaines coordonnees"
    	<< " sont exactement egales aux limites du domaine dans les directions:"
	<< endl;
    for (size_t i=0;i<3;++i)
      if ( coorMatchLocLim[i] )
        cout << "   * " << ( i == 0 ? "X" : i == 1 ? "Y" : "Z" ) <<
      	" translation automatique de 1.e-14" << endl;
  }
  m_wrapper->MPI_Barrier_ActivProc();

  // Nb de particules a inserer sur le proc
  // Recherche des indices min et max des particules dont le centre de gravite
  // est sur le proc
  k = 0;
  found = false;
  while( !found && k < int(m_blocInsert->NX) )
  {
    position[X] = m_blocInsert->bloc.ptA[X] + ( double(k) + 0.5 ) * deltax;
    if ( position[X] > origineLocale[X] ) {kmin = k; found = true;}
    ++k;
  }
  if ( found )
  {
    k = kmax;
    found = false;
    while( !found && k >= 0 )
    {
      position[X] = m_blocInsert->bloc.ptA[X] + ( double(k) + 0.5 ) * deltax;
      if ( position[X] < MaxLocal[X] ) {kmax = k; found = true;}
      --k;
    }
    if ( !found ) no_overlap = true ;
  }
  else no_overlap = true ;

  if ( !no_overlap )
  {
    l = 0;
    found = false;
    while( !found && l < int(m_blocInsert->NY) )
    {
      position[Y] = m_blocInsert->bloc.ptA[Y] + ( double(l) + 0.5 ) * deltay;
      if ( position[Y] > origineLocale[Y] ) {lmin = l; found = true;}
      ++l;
    }
    if ( found )
    {
      l = lmax;
      found = false;
      while( !found && l >= 0 )
      {
        position[Y] = m_blocInsert->bloc.ptA[Y] + ( double(l) + 0.5 ) * deltay;
        if ( position[Y] < MaxLocal[Y] ) {lmax = l; found = true;}
        --l;
      }
      if ( !found ) no_overlap = true ;
    }
    else no_overlap = true ;
  }

  if ( !no_overlap )
  {
    m = 0;
    found = false;
    while( !found && m < int(m_blocInsert->NZ) )
    {
      position[Z] = m_blocInsert->bloc.ptA[Z] + ( double(m) + 0.5 ) * deltaz;
      if ( position[Z] > origineLocale[Z] ) {mmin = m; found = true;}
      ++m;
    }
    if ( found )
    {
      m = mmax;
      found = false;
      while( !found && m >= 0 )
      {
        position[Z] = m_blocInsert->bloc.ptA[Z] + ( double(m) + 0.5 ) * deltaz;
        if ( position[Z] < MaxLocal[Z] ) {mmax = m; found = true;}
        --m;
      }
      if ( !found ) no_overlap = true ;
    }
    else no_overlap = true ;
  }

  if ( no_overlap )
  {
    npartproc = kmin = lmin = mmin = 0;
    kmax = lmax = mmax = -1;
  }
  else npartproc = ( kmax - kmin + 1 ) * ( lmax - lmin + 1 )
  	* ( mmax - mmin + 1 );


  // Nombre de particules de chaque classe a inserer sur le proc
  // Ne marche que si le nombre de particles par classe est superieur
  // au nombre de proc
  if ( m_newParticules.size() == 1 )
    newParticules_.front().second = int(npartproc);
  else
    m_wrapper->distributeParticulesClassProc( m_newParticules,
    	newParticules_, npartproc, ntotalinsert );


  // Initialisation de la numerotation des particules
  size_t* vnbpartproc = m_wrapper->AllGather_UNSIGNED_INT( npartproc );
  int numPartMax = numeroMaxParticules();
  numPartMax = m_wrapper->max_INT( numPartMax );
  if ( numPartMax ) ++numPartMax;
  id = numPartMax;
  for (k=0;k<m_rank;++k) id += int(vnbpartproc[k]);
  delete [] vnbpartproc;


  // Affectation des positions
  for (k=kmin;k<=kmax;++k)
    for (l=lmin;l<=lmax;++l)
      for (m=mmin;m<=mmax;++m)
      {
	// Position
        position[X] = m_blocInsert->bloc.ptA[X] + ( double(k) + 0.5 ) * deltax;
        position[Y] = m_blocInsert->bloc.ptA[Y] + ( double(l) + 0.5 ) * deltay;
        position[Z] = m_blocInsert->bloc.ptA[Z] + ( double(m) + 0.5 ) * deltaz;

        // Selection d'une classe de particules
        particuleClasse = getParticuleClasseForCreation( mode, newParticules_,
		true );

        // Creation de la particule pour insertion ulterieure
        newPart = particuleClasse->createCloneCopy();
        // Previous implementation: not compatible with the compParticule
        // newPart = new Particule( *particuleClasse );
        newPart->setPosition( position );
        newPart->setID( id );
        m_composants.Ajouter( newPart );
	id++;
     }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Lecture des donnees pour les simulations MPI
void GrainsMPI::readDomainDecomposition( DOMNode* root,
  	const Scalar& lx, const Scalar& ly, const Scalar& lz )
{
  // Decomposition de domaine
  int nx, ny, nz = 1;
  DOMNode* decomp = ReaderXML::getNode( root, "DomainDecomposition" );
  nx = ReaderXML::getNodeAttr_Int( decomp, "NX" );
  ny = ReaderXML::getNodeAttr_Int( decomp, "NY" );
  if ( m_dimension == 3 ) nz = ReaderXML::getNodeAttr_Int( decomp, "NZ" );

  int perx = 0, pery = 0, perz = 0;
  DOMNode* mpiperiode = ReaderXML::getNode( root, "MPIperiodes" );
  if ( mpiperiode )
  {
    perx = ReaderXML::getNodeAttr_Int( mpiperiode, "X" );
    if ( perx != 1 ) perx = 0;
    pery = ReaderXML::getNodeAttr_Int( mpiperiode, "Y" );
    if ( pery != 1 ) pery = 0;
    if ( m_dimension == 3 )
      perz = ReaderXML::getNodeAttr_Int( mpiperiode, "Z" );
    if ( perz != 1 ) perz = 0;
  }
  if ( perx || pery || perz ) Grains_Exec::m_MPIperiodique = true;

  // Strategie de communication
  DOMNode* mpiNode = ReaderXML::getNode( root, "MPI" );
  m_MPIstrategie = ReaderXML::getNodeAttr_String( mpiNode, "Strategie" );
  Grains_Exec::m_MPI_verbose =
  	ReaderXML::getNodeAttr_Int( mpiNode, "VerbosityLevel" );
  if ( Grains_Exec::m_MPIperiodique ) m_MPIstrategie = "SRLocalCommOpt";

  // Construction du wrapper MPI
  m_wrapper = new MPIWrapperGrains( nx, ny, nz, perx, pery, perz );
  Grains_Exec::setComm( m_wrapper );
  Grains_Exec::m_MPI = true;
  m_processorIsActiv = m_wrapper->isActiv();
  m_rank = m_wrapper->rank_ACTIV();
  m_nprocs = m_wrapper->nombre_total_procs_ACTIV();
  if ( m_processorIsActiv )
  {
    // Communicateur local
    if ( m_MPIstrategie == "AllgatherLocal" ) m_wrapper->setCommLocal();

    // Grandeurs geometriques locales au processus (domaine local)
    App::setDlocale( lx / m_wrapper->nb_procs_direction(0),
  	ly / m_wrapper->nb_procs_direction(1),
	lz / m_wrapper->nb_procs_direction(2) );
    App::setOriginelocale( m_wrapper->nb_procs_direction(),
  	m_wrapper->MPI_coordonnees() );

    // MPI periodes
    if ( Grains_Exec::m_MPIperiodique )
      m_wrapper->setMPIperiodicVectors( lx, ly, lz );

    // Display le wrapper & la decomposition de domaine
    m_wrapper->display( cout );
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie le nom complet du fichier de resultat .result
string GrainsMPI::fullResultFileName( const string &rootname ) const
{
  string fullname = rootname;
  ostringstream oss;
  oss << "_" << m_rank;
  fullname += oss.str()+".result";

  return fullname;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Lecture des donnees pour les simulations periodiques
void GrainsMPI::readPeriodic( DOMElement* rootElement )
{
  Grains::readPeriodic( rootElement );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Definition du LinkedCell
void GrainsMPI::defineLinkedCell( Scalar const& rayon )
{
  m_sec->define( 2. * rayon, m_wrapper->nb_procs_direction(),
  	m_wrapper->MPI_coordonnees(), m_wrapper->MPI_voisins(),
	m_wrapper->MPI_periodicite() );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Arret du code en cas de probleme
void GrainsMPI::grainsAbort() const
{
  int error_code = 0;
  MPI_Abort( MPI_COMM_WORLD, error_code );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Numero maximum de particules
// A.WACHS - Mars.2010 - Creation
int GrainsMPI::numeroMaxParticules() const
{
  int numMax = m_composants.numeroMaxParticules();
  int collective_numMax = m_wrapper->max_INT( numMax );
  return collective_numMax;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Tire une classe de particules parmi les nouvelles particules a inserer
// A.WACHS - Mars.2010 - Creation
Particule* GrainsMPI::getParticuleClasseForCreation( const PullMode& mode,
  	list< pair<Particule*,int> >& ParticuleClassesForCreation,
	bool const& random_local )
{
  Particule* particuleClasse = NULL;
  list< pair<Particule*,int> >::iterator il;
  list<int> classesDispo;
  bool found = false;
  int i,i0;
  double v;

  switch (mode)
  {
    case ORDER:
      for (il=ParticuleClassesForCreation.begin();
      	il!=ParticuleClassesForCreation.end() && !found;il++)
	if ( il->second != 0 )
	{
	  found = true;
	  particuleClasse = il->first;
	  --il->second;
	}

      if ( !found )
      {
        cout << "Plus de particules disponibles dans aucune classe pour "
		<< "creation" << endl;
	grainsAbort();
      }
      break;

    case RANDOM:
      // On recherche d'abord les classes de particules pour lesquelles
      // il reste des particules a inserer
      i=0;
      for (il=ParticuleClassesForCreation.begin();
      	il!=ParticuleClassesForCreation.end();il++,i++)
	if ( il->second != 0 ) classesDispo.push_back(i);

      // En random local, chaque proc fait son propre tirage aleatoire
      if ( random_local )
      {
        v = double(random()) / double(INT_MAX);
        i0 = int(double(classesDispo.size()) * v);
      }
      // En random global, afin que le tirage al�atoire soit le m�me sur tous
      // les procs seul le master tire au hasard et envoie la classe aux autres
      // procs
      else
      {
        if ( m_rank == 0 )
        {
          v = double(random()) / double(INT_MAX);
          i0 = int(double(classesDispo.size()) * v);
        }
        i0 = m_wrapper->Broadcast_INT( i0 );
      }

      il=ParticuleClassesForCreation.begin();
      for (i=0;i<i0;i++) il++;
      particuleClasse = il->first;
      --il->second;
      break;

    case NONE:
      break;
  }

  return particuleClasse;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Traitement des particules periodiques
void GrainsMPI::periodicParticles( bool b_perfTiming )
{
  if ( b_perfTiming )
    SCT_set_start( "ParticulesPeriodiques" );

  // Ajout des nouvelles references periodiques
  m_sec->addNewPeriodicReference_MPI( m_temps,
		m_composants.getParticulesReferencesPeriodiques() );

  // Cr�ation ou Mise � jour des clones periodiques
  m_wrapper->commParticulesRefPer_AllGatherGlobal_UpdateOrCreate( m_temps,
		m_composants.getParticulesReferencesPeriodiques(),
		m_composants.getParticulesClonesPeriodiques(),
		m_composants.getParticuleClassesReference(),
		m_sec );

  // Creation/Destruction des clones periodiques
  list<int> ClonestoDestroy,ClonestoParticules,PartRefPerHalozone,
  	PartRefPerOutDomainHalozone,InNotRefPerHalozone;
  m_sec->LinkUpdateParticulesPeriodiques_MPI_Step1( m_temps,
		ClonestoDestroy,
  		ClonestoParticules,
		PartRefPerHalozone,
  		PartRefPerOutDomainHalozone,
  		InNotRefPerHalozone,
  	  	m_composants.getParticulesClonesPeriodiques(),
		m_composants.getParticulesReferencesPeriodiques(),
	  	m_composants.getParticulesActives(),
  		m_composants.getParticulesHalozone() );

  // Comm des listes de changement de statut des clones periodiques
  m_wrapper->commPer_AllGatherGlobal_listINT( m_temps,
	  	PartRefPerHalozone );
  m_wrapper->commPer_AllGatherGlobal_listINT( m_temps,
		PartRefPerOutDomainHalozone );
  m_wrapper->commPer_AllGatherGlobal_listINT( m_temps, InNotRefPerHalozone );
  m_wrapper->commPer_AllGatherGlobal_listINT( m_temps,
		ClonestoDestroy );
  m_wrapper->commPer_AllGatherGlobal_listINT( m_temps,
		ClonestoParticules );

  // Changement de statut
  m_composants.statutClonesReferencesPeriodiques_MPI_Step2( m_temps,
  		PartRefPerHalozone, PartRefPerOutDomainHalozone,
		InNotRefPerHalozone, m_sec );
  m_composants.statutClonesPeriodiques_MPI_Step3( m_temps,
  		ClonestoDestroy, ClonestoParticules, m_sec );

  if ( b_perfTiming )
    SCT_add_elapsed_time( "ParticulesPeriodiques" );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecriture de la memoire utilisee par la simulation
void GrainsMPI::display_used_memory() const
{
  m_wrapper->display_used_memory( cout );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Synchronize the PPWindow boolean relative to each sub-domain
void GrainsMPI::synchronize_PPWindow()
{
  vector<bool> b;
  b = PostProcessingWriter::get_PostProcessingWindow();

  for( int i=0; i<m_nprocs; i++ )
  {
    b[i] = m_wrapper->logical_and( b[i] );
    PostProcessingWriter::set_PostProcessingWindow( i, b[i] );
  }

  if( m_rank == 0 )
  {
    cout << "\nProcessors which will write down Paraview particles are : ";
    for( int i=0; i<m_nprocs; i++ )
      if( b[i] > 0.5 )
        cout <<" "<< i <<endl;
  }
}
