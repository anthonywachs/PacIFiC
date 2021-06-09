#include "MPIWrapperGrains.hh"
#include "Grains_Exec.hh"
#include "GrainsCoupledWithFluid.hh"
#include "InterfaceFluide_BuilderFactory.hh"
#include "Contact_BuilderFactory.hh"
#include "LinkedCell.H"
#include "AddedMass.H"
//#include "AppFluide_Drag.H"
// #include "AppFluide_Temperature.H"
#include "AppFluide_LubricationCorrection.H"
#include "Obstacle_BuilderFactory.H"
#include "ObstacleChargement.H"
#include "PostProcessingWriter.hh"
#include "Text_PostProcessingWriter.hh"
#include "PostProcessingWriter_BuilderFactory.hh"
#include "Sphere.H"
#include "Disque.H"


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
AppFluide_LubricationCorrection* GrainsCoupledWithFluid::m_lubricationForce;


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructeur par defaut
GrainsCoupledWithFluid::GrainsCoupledWithFluid( Scalar rhoFluide,
	double grid_size ) :
  Grains(),
  m_InterfaceFluide( NULL ),
  m_explicitAddedMass( NULL ),
  m_bd_completed( false ),
  m_etatSysteme( NULL ),
  m_forceReloadSame( false ),
//  app_HydroForce( NULL ),
//  app_FluidTemperature( NULL ),
  b_fixed_particles( false )
{
  m_gridsize = grid_size;
  Particule::setFluideMasseVolumique(rhoFluide);
  m_InterfaceFluide = InterfaceFluide_BuilderFactory::create();
  Disque::SetvisuNodeNb(40);
  Sphere::SetvisuNodeNbPerQar(5);
  Particule::setMassCorrection( true );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Destructeur
GrainsCoupledWithFluid::~GrainsCoupledWithFluid()
{
  if ( !m_bd_completed ) BeforeDestructor();
  delete m_InterfaceFluide;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Construction du probleme.
void GrainsCoupledWithFluid::Construction( DOMElement* rootElement )
{
  assert(rootElement != NULL);
  DOMNode* root = ReaderXML::getNode( rootElement, "Construction" );

  string restart;

  // Recipient
  DOMNode* recipient = ReaderXML::getNode( root, "Recipient" );
  double lx = ReaderXML::getNodeAttr_Double( recipient, "LX" );
  double ly = ReaderXML::getNodeAttr_Double( recipient, "LY" );
  double lz = ReaderXML::getNodeAttr_Double( recipient, "LZ" );
  DOMNode* origine_recipient = ReaderXML::getNode( root, "Origine" );
  double ox = 0., oy = 0., oz = 0. ;
  if ( origine_recipient )
  {
    ox = ReaderXML::getNodeAttr_Double( origine_recipient, "OX" );
    oy = ReaderXML::getNodeAttr_Double( origine_recipient, "OY" );
    oz = ReaderXML::getNodeAttr_Double( origine_recipient, "OZ" );
  }
  App::setD( lx, ly, lz, ox, oy, oz );

  // Decomposition de domaine
  readDomainDecomposition( root, lx - ox, ly - oy, lz - oz );

  // Schema d'integration en temps
  DOMNode* timeIntegrator = ReaderXML::getNode( root, "TimeIntegration" );

  if (timeIntegrator)
    Grains_Exec::m_TIScheme = ReaderXML::getNodeAttr_String( timeIntegrator,
    	"Type" );

  if ( m_processorIsActiv )
  {
    // Definition des "applications" accessibles
    // L'algorithme sec est present par defaut.
    m_sec = new LinkedCell();
    m_sec->setName( "LinkedCell" );
    m_allApp.push_back( m_sec );

    // Reload ?
    DOMNode* reload = ReaderXML::getNode( root, "Reload" );
    if ( reload )
    {
      // Mode de restart
      string reload_type;
      if ( m_forceReloadSame ) reload_type = "same" ;
      else reload_type = ReaderXML::getNodeAttr_String( reload, "Type" );

      if ( reload_type == "new" || reload_type == "same" )
        Grains_Exec::m_ReloadType = reload_type ;

      // Lecture du nom de fichier de simulation precedent pour reload
      // Si le mode est "same", le fichier de reload est le m�me que
      // le fichier de sortie
      if ( Grains_Exec::m_ReloadType == "new" )
        restart  = ReaderXML::getNodeAttr_String( reload, "Fichier" );
      else
      {
        DOMNode* rootSimu = ReaderXML::getNode( rootElement, "Simulation" );
	DOMNode* fileRestartOutput = ReaderXML::getNode( rootSimu, "Fichier" );
        restart = ReaderXML::getNodeValue_String( fileRestartOutput );
	restart = Grains_Exec::restartFileName_AorB( restart, "_RFTable.txt" );
	Grains_Exec::m_reloadFile_suffix =
		restart.substr( restart.size()-1, 1 );
      }
      restart = fullResultFileName( restart );
      cout << "  Restart du fichier " << restart << endl;

      // Extrait le repertoire de reload a partir du fichier principal de
      // restart
      Grains_Exec::m_ReloadDirectory = Grains_Exec::extractRoot( restart );

      string   cle;
      ifstream simulLoad( restart.c_str() );
      simulLoad >> cle >> m_temps;
			m_new_reload_format = false ;
      if ( cle == "NEW_RELOAD_FORMAT" )
      {
        m_new_reload_format = true ;
        simulLoad >> cle >> m_temps;
      }
      else if (cle == "NEW_RELOAD_FORMAT_AND_CONTACT_HISTORY")
      {
        // Contact history is the new standard - just as reload 2014 - but these
        // conditions ensure previous cases are still compatible.
        m_new_reload_format = true ;
        m_history_storage = true ;
        simulLoad >> cle >> m_temps;
      }
      Contact_BuilderFactory::reload( simulLoad );
      m_composants.read( simulLoad, restart, m_new_reload_format,
			 										m_history_storage);
      Contact_BuilderFactory::set_materialsForObstaclesOnly_reload(
      		m_composants.getParticuleClassesReference() );
      simulLoad >> cle;

      string check_matA, check_matB;
      bool contactLaws_ok = Contact_BuilderFactory::checkContactLawsExist(
    	check_matA, check_matB );
      if ( !contactLaws_ok )
      {
        if ( m_rank == 0 )
          cout << "Pas de loi de contact disponible pour les materiaux : "
	   << check_matA << " & " << check_matB << endl;
         grainsAbort();
      }

      string reset = ReaderXML::getNodeAttr_String( reload, "Vitesse" );
      m_composants.ResetCinematique( reset );
    }

    // Construction du probleme et affectation des composants
    m_rayon = m_composants.getRayonMax();
    DOMNode* rootForces = ReaderXML::getNode( rootElement, "Forces" );
    DOMNode* nLubrication = ReaderXML::getNode( rootForces, "LubricationForce");
    // When the particle is described with 16 grid cells dx, the critical
    // distance dc to activate lubrication correction is 1 grid cells,
    // above N=16, dc/dx increases slightly and linearly with N.
    // The size of Linkedcell is increased wrt dc.
    // ATTENTION: this is implemented differently in
    // GrainsCoupledWithFluidMPI since dc = 0.5R in meso-scale model
    // We assume that Grains MPI is always used for meso-scale and
    // Grains Serial is always used for micro-scale model
    if ( nLubrication )
    {
      if ( m_rank == 0 ) cout << "   LinkedCell size increased "<< 100. *
      m_gridsize * ( 1. +  ( 2. * m_rayon / m_gridsize - 16. ) * 0.1 / 8. )
      /(2. * m_rayon) <<" percent due to the lubrication correction " << endl;
      m_rayon = m_rayon + 0.5 * m_gridsize *
      	( 1. +  ( 2. * m_rayon / m_gridsize - 16. ) * 0.1 / 8. );
      Grains_Exec::m_withlubrication = true;
      double eps_cut;
      if ( ReaderXML::hasNodeAttr_Double( nLubrication,"eps_cut" ) )
        eps_cut = ReaderXML::getNodeAttr_Double( nLubrication,"eps_cut" );
      else
        eps_cut = 0.02;
      GrainsCoupledWithFluid::m_lubricationForce =
       new AppFluide_LubricationCorrection( m_gridsize, true, eps_cut );
      m_allApp.push_back( GrainsCoupledWithFluid::m_lubricationForce );
    }
    else GrainsCoupledWithFluid::m_lubricationForce = NULL;

    // Cells can be widened by the user
    Scalar LC_coef = 1.;
    DOMNode* cellsize = ReaderXML::getNode( root, "CellSize" );
    if ( cellsize )
      LC_coef = ReaderXML::getNodeAttr_Double( cellsize, "Factor" );
    if ( LC_coef < 1. ) LC_coef = 1.;
    defineLinkedCell( LC_coef * m_rayon );

    m_sec->Link( m_composants.getObstacles() );

    if ( m_rank == 0 )
    {
      cout << "Traitement des contacts particule-obstacle "
    	<< "dans le LinkedCell" << endl;
      cout << endl << "Schema d'integration en temps = " <<
      	Grains_Exec::m_TIScheme << endl << endl;
    }

    // Nb de particules sur tous les proc
    m_composants.setNbreParticulesOnAllProc( m_composants.nbreParticules() );

    // Postprocess contact energy dissipation rate
    DOMNode* rootSimu = ReaderXML::getNode( rootElement, "Simulation" );
    DOMNode* nContDiss = ReaderXML::getNode( rootSimu,
    "ContactDissipationRate" );
    if (nContDiss) Grains_Exec::m_ContactDissipation = true;
  }

  // Nb of particles per class
  int nbPC = int(m_composants.getParticuleClassesReference()->size());
  m_composants.prepare_Polydisperse_vectors( nbPC );
  m_composants.set_NbParticlesPerClass();
  m_composants.compute_SauterMeanDiameter();
  m_composants.compute_ParticleClassesConcentration();
}



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Construction de la simulation.
void GrainsCoupledWithFluid::Chargement( DOMElement* rootElement )
{
  if ( m_processorIsActiv )
  {
    assert(rootElement != NULL);
    DOMNode* root = ReaderXML::getNode( rootElement, "Simulation" );

    // Increment
    DOMNode* increment = ReaderXML::getNode( root, "IncreTemps" );
    m_dt = ReaderXML::getNodeAttr_Double( increment, "dt" );

    // Simul Temps
    DOMNode* simulTemps = ReaderXML::getNode( root, "SimulTemps" );
    m_simulTime = ReaderXML::getNodeAttr_Double( simulTemps, "t" );
    if (app_HydroForce)
      app_HydroForce->set_simultime( m_simulTime );
    if (app_FluidTemperature)
      app_FluidTemperature->set_simultime( m_simulTime );

		// Basilisk
    DOMNode* basilisk = ReaderXML::getNode( root, "Basilisk" );
    m_dtgmin = ReaderXML::getNodeAttr_Double( basilisk, "dtgmin" );
    m_dtgmax = ReaderXML::getNodeAttr_Double( basilisk, "dtgmax" );
    m_N = ReaderXML::getNodeAttr_Int( basilisk, "N" );

    // Fichier
    DOMNode* file = ReaderXML::getNode( root, "Fichier" );
    m_fileSave = ReaderXML::getNodeValue_String( file );
    if ( Grains_Exec::m_ReloadType == "new" ) clearResultXmlFiles();
    Grains_Exec::m_SaveDirectory = Grains_Exec::extractRoot( m_fileSave );
    DOMNode* writingMode = ReaderXML::getNode( root, "ModeEcriture" );
    if ( writingMode )
    {
      string wmode = ReaderXML::getNodeValue_String( writingMode );
      if ( wmode == "Hybride" ) Grains_Exec::m_writingModeHybrid = true ;
    }

    if ( m_rank == 0 )
    {
      cout << "Sauvegarde des fichiers de restart" << endl;
      cout << "   Nom du fichier = " << m_fileSave << endl;
      cout << "   Repertoire = " << Grains_Exec::m_SaveDirectory << endl;
      cout << endl;
    }

    // Chargements
    DOMNode* nChargements = ReaderXML::getNode( root, "Chargements" );
    if ( nChargements )
    {
      DOMNode* nDeplaceObs = ReaderXML::getNode( nChargements,
      		"DeplacementGeometrique" );
      if ( nDeplaceObs )
      {
        string deplobsval = ReaderXML::getNodeAttr_String( nDeplaceObs,
      		"Value" );
	if ( deplobsval == "false" ) Obstacle::setDeplaceObstacle( false ) ;
      }
      DOMNodeList* allChargements = ReaderXML::getNodes( nChargements );
      for (XMLSize_t i=0; i<allChargements->getLength(); i++)
      {
        DOMNode* nChargement = allChargements->item( i );
        ObstacleChargement* chargement = new ObstacleChargement( nChargement,
		m_dt, m_rank );
        m_composants.Associer( *chargement );
      }
    }

    // Post-processing writers
    DOMNode* nPostProcessors = ReaderXML::getNode( root,
    	"PostProcessingWriters" );
    if ( nPostProcessors )
    {
      DOMNodeList* allPPW = ReaderXML::getNodes( nPostProcessors );
      for (XMLSize_t i=0; i<allPPW->getLength(); i++)
      {
        DOMNode* nPPW = allPPW->item( i );
        PostProcessingWriter* ppw = PostProcessingWriter_BuilderFactory::create(
      		nPPW, m_rank, m_nprocs );
        if ( ppw ) m_composants.addPostProcessingWriter( ppw );
      }
      // Tell app fluid drag to record fluid/solid slip velocity
      if( app_HydroForce )
        app_HydroForce->set_slipveloutput( Text_PostProcessingWriter::
            b_slipVelocity );

      PostProcessingWriter::allocate_PostProcessingWindow( m_nprocs );

      // Eventual Window where Post-processing Writers are active
      DOMNode* nRestrictedParaviewWindow = ReaderXML::getNode( root,
      	"RestrictedParaviewWindow" );
      if ( nRestrictedParaviewWindow )
      {
        DOMNodeList* nWindowPoints = ReaderXML::getNodes(
		nRestrictedParaviewWindow );

 	DOMNode* pointA = nWindowPoints->item( 0 );
 	DOMNode* pointB = nWindowPoints->item( 1 );

 	Fenetre PPWindow;
	PPWindow.ftype = FENETRE_BOX;
	PPWindow.radius = PPWindow.radius_int = PPWindow.hauteur = 0. ;
	PPWindow.axisdir = W ;
 	PPWindow.ptA[X] = ReaderXML::getNodeAttr_Double( pointA, "X" );
	PPWindow.ptA[Y] = ReaderXML::getNodeAttr_Double( pointA, "Y" );
	PPWindow.ptA[Z] = ReaderXML::getNodeAttr_Double( pointA, "Z" );
	PPWindow.ptB[X] = ReaderXML::getNodeAttr_Double( pointB, "X" );
	PPWindow.ptB[Y] = ReaderXML::getNodeAttr_Double( pointB, "Y" );
	PPWindow.ptB[Z] = ReaderXML::getNodeAttr_Double( pointB, "Z" );

	double Ox, Oy, Oz, lx, ly, lz;
	bool b_X=false, b_Y=false, b_Z=false, b_PPWindow=false;
	App::getOrigineLocale( Ox, Oy, Oz );
	App::getDimensionsLocales( lx, ly, lz );

	if( (PPWindow.ptA[X]>=Ox && PPWindow.ptA[X]<Ox+lx)||
	    (PPWindow.ptB[X]>=Ox && PPWindow.ptB[X]<Ox+lx)||
	    (Ox>PPWindow.ptA[X] && Ox<PPWindow.ptB[X])||
	    (Ox>PPWindow.ptB[X] && Ox<PPWindow.ptA[X]) )
	  b_X = true;
	if( (PPWindow.ptA[Y]>=Oy && PPWindow.ptA[Y]<Oy+ly)||
	    (PPWindow.ptB[Y]>=Oy && PPWindow.ptB[Y]<Oy+ly)||
	    (Oy>PPWindow.ptA[Y] && Oy<PPWindow.ptB[Y])||
	    (Oy>PPWindow.ptB[Y] && Oy<PPWindow.ptA[Y]) )
	  b_Y = true;
	if( (PPWindow.ptA[Z]>=Oz && PPWindow.ptA[Z]<Oz+lz)||
	    (PPWindow.ptB[Z]>=Oz && PPWindow.ptB[Z]<Oz+lz)||
	    (Oz>PPWindow.ptA[Z] && Oz<PPWindow.ptB[Z])||
	    (Oz>PPWindow.ptB[Z] && Oz<PPWindow.ptA[Z]) )
	  b_Z = true;

	if( b_X && b_Y && b_Z )
	  b_PPWindow = true;

	PostProcessingWriter::set_PostProcessingWindow( m_rank, b_PPWindow );

	if ( m_nprocs > 1 )
	  synchronize_PPWindow();
      }
    }

    // Frequence de mise a jour du lien entre obstacles et LinkedCell
    DOMNode* nodeUpdateFreq = ReaderXML::getNode( root, "LinkUpdate" );
    int updateFreq = 1;
    if ( nodeUpdateFreq )
      updateFreq = ReaderXML::getNodeAttr_Int( nodeUpdateFreq, "frequence" );
    m_composants.setObstaclesLinkedCellUpdateFrequency( updateFreq );

    // Deplacement geometrique des obstacles
    DOMNode* nodeDeplaceObs = ReaderXML::getNode( root,
      		"DeplacementGeometrique" );
    if ( nodeDeplaceObs )
    {
      string deplobsval = ReaderXML::getNodeAttr_String( nodeDeplaceObs,
      		"value" );
      if ( deplobsval == "false" ) Obstacle::setDeplaceObstacle( false ) ;
    }

    // Post-processing des efforts sur les obstacles
    DOMNode* ppObstacles = ReaderXML::getNode( root, "EffortsObstacles" );
    if ( ppObstacles )
    {
      int outputFreq = ReaderXML::getNodeAttr_Int( ppObstacles, "Frequence" );
      string ppObsroot = ReaderXML::getNodeAttr_String( ppObstacles, "Root" );
      list<string> allppObsName;
      DOMNodeList* allppObs = ReaderXML::getNodes( ppObstacles );
      for (XMLSize_t i=0; i<allppObs->getLength(); i++)
      {
        DOMNode* nppObs = allppObs->item( i );
        allppObsName.push_back(
		ReaderXML::getNodeAttr_String( nppObs, "Name" ) );
      }
      m_composants.setOutputObstaclesLoadParameters( ppObsroot,
  	outputFreq, allppObsName );
    }
  }
}





// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Set the value Delta_t for one time (for coupling with Basilisk which handles adaptive time-steps)
void GrainsCoupledWithFluid::set_timeStep( double const dtfluid )
{

  if ( m_processorIsActiv ) {
    double dteff;
    dteff = max(min (m_dtgmax, dtfluid/m_N), m_dtgmin);
    m_N = int(dtfluid/dteff);

    m_dt = dtfluid/m_N;
    m_simulTime = dtfluid;

    // cout << "Can: granular time-step: dtg " << m_dt << endl;
    // cout << "Can: number of granular time-steps " << m_N << endl;
    // cout << "Can: fluid time-step: dtf within grains " << dtfluid << endl;

  }

}
void GrainsCoupledWithFluid::set_ExplicitAddedMasstimeStep( double dtfluid )
{

  if ( m_processorIsActiv )
    if ( m_explicitAddedMass )
      dynamic_cast<AddedMass*>(m_explicitAddedMass)->setsimulTime( dtfluid );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Construction des forces actives
void GrainsCoupledWithFluid::Forces( DOMElement* rootElement )
{
  if ( m_processorIsActiv )
  {
    assert(rootElement != NULL);
    DOMNode* root = ReaderXML::getNode( rootElement, "Forces" );

    if ( m_rank == 0 ) cout << "Simulation avec :\n";

    // Gravite ?
    DOMNode* nGravite = ReaderXML::getNode( root, "Gravite" );
    if (nGravite) {
      Grains_Exec::m_vgravite[X] = ReaderXML::getNodeAttr_Double(
      	nGravite, "GX" );
      Grains_Exec::m_vgravite[Y] = ReaderXML::getNodeAttr_Double(
      	nGravite, "GY" );
      Grains_Exec::m_vgravite[Z] = ReaderXML::getNodeAttr_Double(
      	nGravite, "GZ" );
      if ( m_rank == 0 ) cout << "  Gravite\n";
    }
    else
    {
      if ( m_rank == 0 ) cout << "Gravite obligatoire !!";
      grainsAbort();
    }

    // Calcul du poids des particules
    m_composants.computeWeight( 0., 0. );

    // Hydro Force ?
    DOMNode* nHydro = ReaderXML::getNode( root, "DragForce" );
    if ( nHydro )
    {
      string isLift, isPressureGradient, is_added_mass_demcfd,isStochasticDrag;
      if ( ReaderXML::hasNodeAttr_String( nHydro, "WithLift" ) )
        isLift = ReaderXML::getNodeAttr_String( nHydro, "WithLift" );
      if( isLift == "yes" )
      {
        Grains_Exec::m_withLiftForce = true;
        if ( m_rank == 0 ) cout << "  Lift force and hydroTorque\n";
      }
      if ( ReaderXML::hasNodeAttr_String( nHydro, "WithAddedMass" ) )
        is_added_mass_demcfd = ReaderXML::getNodeAttr_String( nHydro, "WithAddedMass" );

      if ( is_added_mass_demcfd =="yes")
      {
        Grains_Exec::m_addedmass_demcfd = true;
        AddAddedMassInfos2Particles_demcfd();
        if ( m_rank == 0 ) cout << "  Added mass\n";
      }

      if ( ReaderXML::hasNodeAttr_String( nHydro, "WithPressureGradient" ) )
        isPressureGradient = ReaderXML::getNodeAttr_String( nHydro,
                "WithPressureGradient" );
      if( isPressureGradient == "yes" )
      {
        Grains_Exec::m_withPressureGradient = true;
        if ( m_rank == 0 ) cout << "  Pressure Gradient in hydro force\n";
      }
      if (ReaderXML::hasNodeAttr_String( nHydro, "WithStochasticDrag") )
	isStochasticDrag = ReaderXML::getNodeAttr_String( nHydro, "WithStochasticDrag");
      if ( isStochasticDrag == "yes" )
      {
        Grains_Exec::m_withStochasticDrag = true;
        if ( m_rank == 0 ) cout << "  Stochastic Drag force\n";
      }

      app_HydroForce = new AppFluide_Drag( nHydro );
      app_HydroForce->setName( "hydroForce" );
      m_allApp.push_back( app_HydroForce );
      Grains_Exec::m_withHydroForce = true;
      if ( m_rank == 0 ) cout << "  Trainee hydrodynamique" << endl;
    }

    // Fluid-Solid heat transfert ?
    DOMNode* nTemperature = ReaderXML::getNode( root, "Temperature" );
    if( nTemperature )
    {
      app_FluidTemperature = new AppFluide_Temperature( nTemperature );
      app_FluidTemperature->setName( "Temperature" );
      m_allApp.push_back( app_FluidTemperature );
      if ( m_rank == 0 ) cout << "  Temperature" << endl;
//      Grains_Exec::m_withFluidTemperature = true;
    }

    // Lubrication correction
    if ( m_lubricationForce )
      if ( m_rank == 0 ) cout << "  Lubrication force\n";

    // Affectation des particules aux forces de contact
    m_composants.Link( *m_sec );

    if ( m_rank == 0 ) cout << endl;
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Sauvegarde de l'etat de simulation */
void GrainsCoupledWithFluid::Save( const string &ext ) const
{
  // Sauvegarde du fichier de caracteristiques pour le fluide
  string fileOUT( "Res/particles_features" );
  fileOUT += ext + ".dat";
  WriteParticulesInFluid( fileOUT );
}




/* Forces initialization step
-----------------------------*/
void GrainsCoupledWithFluid::InitializeForcesT0(
    const double &fluidThermalConductivity_ )
{
  if( Grains_Exec::m_withFluidTemperature )
    app_FluidTemperature->InitializeTemperature( 0., m_dt,
        m_composants.getParticulesActives(), fluidThermalConductivity_ );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Gestion de la simulation
void GrainsCoupledWithFluid::Simulation( bool predict,
	bool isPredictorCorrector, bool explicit_added_mass )
{
  if ( m_processorIsActiv )
  {
    // Initialisation de la cinematique des obstacles
    // Fait au debut du 1er appel par le fluide uniquement, d'ou l'utilisation
    // d'un compteur statique
    static size_t init_counter = 0 ;
    if ( !init_counter )
      m_composants.setCinematiqueObstacleSansDeplacement( m_temps, m_dt );
    ++init_counter;

    Scalar time = 0.;

    // Predicteur/Correcteur
    if ( predict )
    {
      Grains::setMode(true);
      if ( isPredictorCorrector ) saveState();
    }
    else Grains::setMode(false);

    // Boucle sur un pas de temps fluide
    while( m_simulTime - time > 0.01 * m_dt )
    {
      try {
        time  += m_dt;
        m_temps += m_dt;

	if( m_composants.IsShrinking() ) m_composants.ShrinkingRate( m_temps );


        // Initialisation de l'indicateur de calcul
        // de la transformation avec epaiseur de croute a faux
        m_composants.InitializeVdWState( m_temps, m_dt );

        // Initialisation des torseurs de force
        m_composants.InitializeForces( m_temps, m_dt, predict );

        // Creation/Destruction des clones periodiques
        m_sec->LinkUpdateParticulesPeriodiques( m_temps,
            m_composants.getParticulesActives(),
            m_composants.getParticulesClonesPeriodiques(),
            m_composants.getParticuleClassesReference() );

        // Calcul des forces de contact
        m_sec->CalculerForces( m_temps, m_dt,
            m_composants.getParticulesActives() );

        // Calcul des forces de masse ajout�e
        if ( predict && m_explicitAddedMass )
          m_explicitAddedMass->CalculerForces( m_temps, m_dt,
              m_composants.getParticulesActives() );

        // Caclul des forces de trainee hydro en DEM-CFD
        if( app_HydroForce )
          app_HydroForce->CalculerForces( m_temps, m_dt,
              m_composants.getParticulesActives() );

        // Caclul de la temperature en DEM-CFD
        if( app_FluidTemperature )
          app_FluidTemperature->CalculerForces( m_temps, m_dt,
              m_composants.getParticulesActives() );

        // Calcul des contacts de contact sur les clones periodiques
        m_sec->CalculerForcesClonesPeriodiques( m_temps, m_dt,
                m_composants.getParticulesClonesPeriodiques() );
        m_composants.AddForcesFromPeriodicClonesToParticules( m_temps, m_dt );

        if( !b_fixed_particles )
        {
          // Deplacement et actualisation des composants
          m_composants.Deplacer( m_temps, m_dt );
          m_composants.Actualiser();
        }

        // Caclul de la temperature en DEM-CFD
        if( app_FluidTemperature )
       {
          m_composants.ComputeTemperature( m_temps, m_dt );
        }

        // Mise � jour des clones periodiques
        m_composants.updateClonesPeriodiques( m_sec );

        // Actualisation des particules & obstacles dans les cellules
        m_sec->LinkUpdate( m_temps, m_dt,
            m_composants.getParticulesActives() );
      }
      catch (ErreurContact &choc)
      {
        // Fin de simulation sur choc
        cerr << '\n';
        m_composants.PostProcessingErreurComposants( "ErreurContact",
            choc.getComposants() );
        choc.Message(cerr);
        break;
      }
      catch (ErreurDeplacement &errDeplacement)
      {
        // Fin de simulation sur deplacement trop grand
        cerr << '\n';
        errDeplacement.Message(cerr);
        m_composants.PostProcessingErreurComposants( "ErreurDeplacement",
            errDeplacement.getComposant() );
        break;
      }
      catch (ErreurSimulation &erreur)
      {
        // Fin de simulation sur erreur
        cerr << '\n';
        erreur.Message(cerr);
        break;
      }
    }

    cout << "Grains Time = " <<
    	Grains_Exec::doubleToString( ios::scientific, 8, m_temps ) << endl;
    cout << "Grains dtfluid = " <<
    	Grains_Exec::doubleToString( ios::scientific, 8, time ) << endl;
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Modification de la vitesse des particules
void GrainsCoupledWithFluid::UpdateParticulesVelocities(
	const bool &b_set_velocity_nm1_and_diff)
{
  if ( m_processorIsActiv )
  {
    // Mise � jour de la vitesse des particules par le fluide
    m_InterfaceFluide->UpdateParticulesVelocities(
    	*m_composants.getParticulesActives(),
  	m_dt, b_set_velocity_nm1_and_diff );

    // Copie vers les clones periodiques
    // On passe NULL pour LinkedCell car a ce stade les clones multi-periodiques
    // n'existent pas
    m_composants.updateClonesPeriodiques( NULL );
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void GrainsCoupledWithFluid::WriteParticulesInFluid(const string &filename)
	const
{
  if ( m_processorIsActiv )
    m_InterfaceFluide->WriteParticulesInFluid(
    	*m_composants.getParticulesActives(), filename );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Modification de la vitesse des particules
void GrainsCoupledWithFluid::UpdateParticulesVelocities(
	const vector<vector<double> > &velocities,
	const bool &b_set_velocity_nm1_and_diff )
{
  if ( m_processorIsActiv )
  {
    // Mise � jour de la vitesse des particules par le fluide
    m_InterfaceFluide->UpdateParticulesVelocities(
    	*m_composants.getParticulesActives(),
  	m_dt, velocities, b_set_velocity_nm1_and_diff );

    // Copie vers les clones periodiques
    // On passe NULL pour LinkedCell car a ce stade on souhaite simplement
    // mettre a jour les vitesse & position
    m_composants.updateClonesPeriodiques( NULL );
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void GrainsCoupledWithFluid::WriteParticulesInFluid( istringstream &is ) const
{
  if ( m_processorIsActiv )
    m_InterfaceFluide->WriteParticulesInFluid(
    	*m_composants.getParticulesActives(),
    	m_composants.getObstaclesToFluid(), is );
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void GrainsCoupledWithFluid::WriteParticulesInDSFluid( istringstream &is ) const
{
  if ( m_processorIsActiv )
    m_InterfaceFluide->WriteParticulesInDSFluid(
    	*m_composants.getParticulesActives(), is );
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void GrainsCoupledWithFluid::ReadParticulesFromDSFluid( istringstream &is ) const
{
  if ( m_processorIsActiv )
    m_InterfaceFluide->ReadParticulesFromDSFluid(
    	*m_composants.getParticulesActives(), is );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void GrainsCoupledWithFluid::WritePVGCInFluid(const string &filename) const
{
  if ( m_processorIsActiv )
    m_InterfaceFluide->WritePVGCInFluid( *m_composants.getParticulesActives(),
    	filename );
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void GrainsCoupledWithFluid::WriteParticulesInFluid( BasiliskDataStructure * is ) const
{
  if ( m_processorIsActiv )
    m_InterfaceFluide->WriteParticulesInFluid(
    	*m_composants.getParticulesActives(),
    	m_composants.getObstaclesToFluid(), is );
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Modification de la vitesse des particules de Grains par Basilisk
void GrainsCoupledWithFluid::UpdateParticulesVelocities(
	BasiliskDataStructure * is,
	const bool &b_set_velocity_nm1_and_diff )
{
  if ( m_processorIsActiv )
  {
    // Mise a jour de la vitesse des particules par le fluide
    m_InterfaceFluide->UpdateParticulesVelocities(
    	*m_composants.getParticulesActives(),
  	m_dt, is, b_set_velocity_nm1_and_diff );

    // Copie vers les clones periodiques
    // On passe NULL pour LinkedCell car a ce stade on souhaite simplement
    // mettre a jour les vitesse & position
    m_composants.updateClonesPeriodiques( NULL );
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void GrainsCoupledWithFluid::WritePVGCInFluid( istringstream &is ) const
{
  if ( m_processorIsActiv )
    m_InterfaceFluide->WritePVGCInFluid( *m_composants.getParticulesActives(),
    	is );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void GrainsCoupledWithFluid::AddExplicitAddedMass(
	bool const& restart, string const& PelDirRes )
{
  if ( m_processorIsActiv )
  {
    m_explicitAddedMass =
        new AddedMass( Particule::getFluideMasseVolumique(), m_simulTime );
    list<Particule*>* particules = m_composants.getParticulesActives();
    for (list<Particule*>::iterator particule=particules->begin();
        particule!=particules->end();particule++)
      (*particule)->createAddedMassInfos();
    m_allApp.push_back( m_explicitAddedMass );
    Particule::setExplicitMassCorrection( true );
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void GrainsCoupledWithFluid::InitializeExplicitAddedMassRestart(
	bool const& restart, string const& dirRes,
	string const& fluidsolver )
{
  if ( m_processorIsActiv )
  {
    if ( restart )
    {
      if ( fluidsolver == "PeliGRIFF" )
        m_composants.setVelocityAndVelocityDifferencePreviousTimeRestart(
      		dirRes );

      else
      {
        double previousdtfluid = m_composants.
		setVelocityAndVelocityDifferencePreviousTimeRestart_Basilisk(
      		dirRes );
	dynamic_cast<AddedMass*>(m_explicitAddedMass)->setsimulTime(
		previousdtfluid );
      }
    }
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Sauvegarde de l'etat de simulation */
void GrainsCoupledWithFluid::saveState()
{
  m_etatSysteme = new SystemState( m_temps, m_composants );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Restauration de l'etat de simulation
void GrainsCoupledWithFluid::restaureState()
{
  m_etatSysteme->RestaureState( m_temps, m_composants, m_sec, m_rank );
  delete m_etatSysteme;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void GrainsCoupledWithFluid::InitializeCorrectorStep(
	const vector<vector<double> > &velocities )
{
  if ( m_processorIsActiv )
  {
    // Recupere l'etat precedent
    restaureState();

    // Mise a jour des vitesses venant du fluide
    UpdateParticulesVelocities( velocities, false );
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Insertion d'une particule dans les algorithmes
bool GrainsCoupledWithFluid::insertParticule(const PullMode& mode)
{
  cout << "!!! Warning: Inserting particules is not permitted in "
  	<< "GrainsCoupledWithFluid and GrainsCoupledWithFluidMPI !!!" << endl;
  cout << "Only reloading a previous simulation file is implemented" << endl;
  grainsAbort();

  return false;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Positionne les particules en attente pour la simulation a partir
// d'un fichier de positions
void GrainsCoupledWithFluid::setPositionParticulesFichier(const PullMode& mode)
{
  cout << "!!! Warning: Setting particules position is not permitted in "
  	<< "GrainsCoupledWithFluid and GrainsCoupledWithFluidMPI !!!" << endl;
  cout << "Only reloading a previous simulation file is implemented" << endl;
  grainsAbort();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Positionne les particules en attente pour la simulation a partir
// d'un bloc structure de positions
void GrainsCoupledWithFluid::setPositionParticulesBloc(const PullMode& mode)
{
  cout << "!!! Warning: Setting particules position is not permitted in "
  	<< "GrainsCoupledWithFluid and GrainsCoupledWithFluidMPI !!!" << endl;
  cout << "Only reloading a previous simulation file is implemented" << endl;
  grainsAbort();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Sauvegarde par defaut de l'etat initial pour post-processing
void GrainsCoupledWithFluid::InitialPostProcessing()
{
  // Par defaut l'etat initial est sauvegarde
  if ( m_processorIsActiv )
    m_composants.PostProcessing_start( m_temps, m_dt, m_sec, m_fenetres );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Sauvegarde pour post-processing et restart
void GrainsCoupledWithFluid::doPostProcessing()
{
  if ( m_processorIsActiv )
  {
    // Post processing
    m_composants.PostProcessing( m_temps, m_dt, m_sec );

    // Sauvegarde du fichier de fin pour Reload
    saveReload( m_temps );
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Sauvegarde pour post-processing et restart
void GrainsCoupledWithFluid::doPostProcessingEvo( const double &time )
{
  if ( m_processorIsActiv )
  {
    // Post processing
    m_composants.PostProcessing( time, m_dt, m_sec );

    // Sauvegarde du fichier de fin pour Reload
    saveReload( m_temps );
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Operations a effectuer avant appel du destructeur
void GrainsCoupledWithFluid::BeforeDestructor()
{
  if ( m_processorIsActiv )
  {
    // Avant destruction, sauvegarde de l'etat final
    m_composants.PostProcessing_end();

    // Sauvegarde du fichier de fin pour Reload
    // !!! le temps pour le reload n'a pas vraiment d'importance !!!
//    saveReload( m_temps );
  }
  m_bd_completed = true;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Temps de simulation (pas de temps fluide)
double GrainsCoupledWithFluid::getSimulTime() const
{
  return m_simulTime;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoi le vecteur gravite
vector<double> GrainsCoupledWithFluid::getVecteurGravite() const
{
  vector<double> vg(3,0.);
  for (int i=0;i<3;++i) vg[i] = Grains_Exec::m_vgravite[i];

  return vg;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Send the number of particles on all proc
size_t GrainsCoupledWithFluid::getNumberOfParticleOnAllProc() const
{
  return Grains_Exec::nbreParticulesOnAllProc();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoi la dimension de l'espace physique
int GrainsCoupledWithFluid::getDimension() const
{
  return m_dimension;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoi le booleen pour savoir si la force de portance demcfd est active
bool GrainsCoupledWithFluid::getIsLiftActive() const
{
  return Grains_Exec::m_withLiftForce;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoi le booleen pour savoir si latemperature demcfd est active
bool GrainsCoupledWithFluid::getIsTemperature() const
{
  return Grains_Exec::m_withFluidTemperature;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Verifie que le post processing Paraview est actif, sinon le cree
void GrainsCoupledWithFluid::checkParaviewPostProcessing( const string &name_,
	const string &root_,
  	const bool &isBinary )
{
  if ( m_processorIsActiv )
    m_composants.checkParaviewPostProcessing( m_rank, m_nprocs, name_,
    	root_, isBinary );
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Verifie que le post processing Paraview est actif pour Basilisk, sinon le cree
// Can Selcuk 2018
void GrainsCoupledWithFluid::checkParaviewPostProcessing( const char  * name_,
	const char * root_,
  	const bool &isBinary ) {

  string str(name_);
  string str_root(root_);

  if ( m_processorIsActiv ) {
    m_composants.checkParaviewPostProcessing( m_rank, m_nprocs, str,
					      str_root, isBinary );
  }
}






// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Verifie que le post processing Matlab est actif, sinon le cree
void GrainsCoupledWithFluid::checkMatlabPostProcessing( const string &name_,
	const string &root_,
  	const bool &isBinary )
{
  if ( m_processorIsActiv )
    m_composants.checkMatlabPostProcessing( m_rank, m_nprocs, name_,
  	root_, isBinary );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Definit le numero initial de cycle d'ecriture des fichiers de post
// processing
void GrainsCoupledWithFluid::setInitialCycleNumber( const int& cycle0 )
{
  if ( m_processorIsActiv ) m_composants.setInitialCycleNumber( cycle0 );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Definit la valeur de la translation de post processing dans
// Paraview pour le cas projection-translation
void GrainsCoupledWithFluid::setParaviewPostProcessingTranslationVector(
	double const& tvx, double const& tvy, double const& tvz )
{
  if ( m_processorIsActiv )
  {
    if ( !Grains_Exec::m_translationParaviewPostProcessing )
      Grains_Exec::m_translationParaviewPostProcessing = new Vecteur();
    (*Grains_Exec::m_translationParaviewPostProcessing)[X] = tvx;
    (*Grains_Exec::m_translationParaviewPostProcessing)[Y] = tvy;
    (*Grains_Exec::m_translationParaviewPostProcessing)[Z] = tvz;
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Retourne la liste des particules
list<Particule*>* GrainsCoupledWithFluid::getListOfParticles()
{
  return m_composants.getParticulesActives();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Allocate DEM-CFD_FluidInfos structure
void GrainsCoupledWithFluid::allocateDEMCFD_FluidInfos()
{
  Particule::setMassCorrection( false );
  Grains::allocateDEMCFD_FluidInfos();
  Grains_Exec::m_withdemcfd = true;
//  Particule::setMassCorrection( false );
//  list<Particule*>* allParticles = m_composants.getParticulesActives();
//  list<Particule*>::iterator il;
//  for(il=allParticles->begin(); il!=allParticles->end(); il++)
//    (*il)->allocateDEMCFD_FluidInfos();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Definit la viscosite du fluide environnant
void GrainsCoupledWithFluid::setFluidViscosity( double mu )
{
  Particule::setFluidViscosity( mu );
  if ( GrainsCoupledWithFluid::m_lubricationForce ) GrainsCoupledWithFluid::
  	m_lubricationForce->set_viscosity( mu );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Set boolean "b_fixed_particles" to true in grains3D.
//    Then, we still go through GRAINS3D, we still compute forces, but
//    we dont move  particles
void GrainsCoupledWithFluid::set_fixed_particles( )
{
  b_fixed_particles = true;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Set boolean "b_fixed_particles_temp" to true in grains3D.
//    Then, we still go through GRAINS3D, we still compute forces, but
//    we dont move  particles
void GrainsCoupledWithFluid::set_fixed_particles_temp( )
{
  b_fixed_particles_temp = true;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Get processor's origin coordinates
void GrainsCoupledWithFluid::getProcessorOrigin( double& x_solid,
		double& y_solid, double& z_solid )
{
  App::getOrigineLocale( x_solid, y_solid, z_solid );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Initialise les clones periodiques
void GrainsCoupledWithFluid::initializeClonesPeriodiques()
{
  if ( m_processorIsActiv )
  {
    // Creation/Destruction des clones periodiques
    m_sec->LinkUpdateParticulesPeriodiques( m_temps,
		m_composants.getParticulesActives(),
  		m_composants.getParticulesClonesPeriodiques(),
		m_composants.getParticuleClassesReference() );
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Calcule I.w et I.dw pour chaque particule pour evaluer le moment
// hydrodynamique exerce par le fluide
void GrainsCoupledWithFluid::compute_Iw_Idw(
	vector< vector<double> > const* w_dw,
  	vector< vector<double> >* Iw_Idw )
{
  if ( m_processorIsActiv )
  {
    list<Particule*> const* pactives = m_composants.getParticulesActives();
    list<Particule*>::const_iterator particule;
    Vecteur w, res ;
    int pID = 0 ;

    for (particule=pactives->begin(); particule!=pactives->end();particule++)
    {
      pID = (*particule)->getID() ;

      // Calcul de I.w
      w[X] = (*w_dw)[pID][X];
      w[Y] = (*w_dw)[pID][Y];
      w[Z] = (*w_dw)[pID][Z];

      res = (*particule)->getCinematique()->calculIdwExplicite(
	w, (*particule)->getInertie() );

      (*Iw_Idw)[pID][0] = res[X];
      (*Iw_Idw)[pID][1] = res[Y];
      (*Iw_Idw)[pID][2] = res[Z];

      // Calcul de I.dw
      w[X] = (*w_dw)[pID][3];
      w[Y] = (*w_dw)[pID][4];
      w[Z] = (*w_dw)[pID][5];

      res = (*particule)->getCinematique()->calculIdwExplicite(
	w, (*particule)->getInertie() );

      (*Iw_Idw)[pID][3] = res[X];
      (*Iw_Idw)[pID][4] = res[Y];
      (*Iw_Idw)[pID][5] = res[Z];
    }
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecrit les efforts sur les obstacles dans des fichiers
void GrainsCoupledWithFluid::outputObstaclesLoad( double temps )
{
  if ( m_processorIsActiv )
    m_composants.outputObstaclesLoad( temps, m_dt, true );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Initialise les fichiers de sortie des efforts sur les obstacles
void GrainsCoupledWithFluid::initialiseOutputObstaclesLoadFiles( double temps )
{
  if ( m_processorIsActiv )
    m_composants.initialiseOutputObstaclesLoadFiles( m_rank, true, temps );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Return pointer to the lubrication class
AppFluide_LubricationCorrection* GrainsCoupledWithFluid::LubricationCorrection()
{
return GrainsCoupledWithFluid::m_lubricationForce;
}





// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void GrainsCoupledWithFluid::AddAddedMassInfos2Particles_demcfd()
{
  if ( m_processorIsActiv )
  {
    list<Particule*>* particules = m_composants.getParticulesActives();
    for (list<Particule*>::iterator particule=particules->begin();
        particule!=particules->end();particule++)
      (*particule)->createAddedMassInfos();
  }
}
