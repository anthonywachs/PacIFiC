#include "Grains.H"
#include "Contact_BuilderFactory.hh"
#include "LinkedCell.H"
#include "AppFluide_Drag.H"
#include "Obstacle_BuilderFactory.H"
#include "ObstacleChargement.H"
#include "Grains_BuilderFactory.H"
#include "PostProcessingWriter.hh"
#include "PostProcessingWriter_BuilderFactory.hh"
#include "Text_PostProcessingWriter.hh"
#include "ObstaclePeriodique.hh"
#include "CompParticule.hh"
#include "stdlib.h"


// Initialisation des attributs static
bool Grains::m_predictor_mode = true;


//-----------------------------------------------------------------------------
// Constructeur par defaut
Grains::Grains() :
  ComputingTime("Solver"),
  m_sec( NULL ),
  app_HydroForce( NULL ),
  app_FluidTemperature( NULL ),
  app_SolidTemperature( NULL ),
  m_dimension( 3 ),
  m_allProcTiming( true ),
  m_new_reload_format( false ),
  m_history_storage( false ),
//  m_dragForce_with_fluidAtRest( false ),
  m_mode_insertion_particules( NONE ),
  m_methode_insertion_particules( NOINSERT ),
  m_methode_initvit_particules( ZERO ),
  m_blocInsert( NULL ),
  m_configAleatoire( false ),
  m_insertion_frequency( 1 ),
  m_force_insertion( false ),
  m_RandomMotionCoefTrans( 0. ),
  m_RandomMotionCoefRot( 0. ),
  m_rank( 0 ),
  m_nprocs( 1 ),
  m_processorIsActiv( true )
{
  if ( Grains_BuilderFactory::getContext() == DIM_2 ) m_dimension = 2;
}




//-----------------------------------------------------------------------------
// Destructeur
Grains::~Grains()
{
  if ( m_blocInsert ) delete m_blocInsert;
  list<App*>::iterator app;
  for (app=m_allApp.begin(); app!=m_allApp.end(); app++) delete *app;
  m_newParticules.clear();
  Contact_BuilderFactory::eraseAllContactLaws();
  Grains_Exec::GarbageCollector();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Construction du probleme.
void Grains::Construction( DOMElement* rootElement )
{
  assert(rootElement != NULL);
  DOMNode* root = ReaderXML::getNode( rootElement, "Construction" );

  bool brestart = false, bnewpart = false, bnewobst = false, bnewper = false;
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
  if ( timeIntegrator )
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
      brestart = true;

      // Mode de restart
      string reload_type = ReaderXML::getNodeAttr_String( reload, "Type" );
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

      // Extrait le repertoire de reload a partir du fichier principal de
      // restart
      Grains_Exec::m_ReloadDirectory = Grains_Exec::extractRoot( restart );
      // Lecture du fichier de reload
      // WRN : gestion des contacts pas simple & utilise un max de pointeurs
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
                          m_history_storage );
      Contact_BuilderFactory::set_materialsForObstaclesOnly_reload(
          m_composants.getParticuleClassesReference() );
      simulLoad >> cle;

      // Remise a zero ou non de la vitesse des particules
      string reset = ReaderXML::getNodeAttr_String( reload, "Vitesse" );
      m_composants.ResetCinematique( reset );
    }

    // Particules
    DOMNode* particules = ReaderXML::getNode( root, "Particules" );
    if ( particules )
    {
      bnewpart = true;

      // /!\ At this stage, nbPC=0 because ParticuleClassesReference is empty
      int nbPC = int(m_composants.getParticuleClassesReference()->size());

      DOMNodeList* allParticules = ReaderXML::getNodes( rootElement,
      	"Particule" );
      m_composants.prepare_Polydisperse_vectors(
          int(allParticules->getLength()) );

      for (XMLSize_t i=0; i<allParticules->getLength(); i++)
      {
        DOMNode* nParticule = allParticules->item( i );
        // /!\ We read the number of particles in insert.xml, but it can be
        // larger than the number of positions in Grains/Init/init_position.data
        int nbre = ReaderXML::getNodeAttr_Int( nParticule, "Nombre" );

        // Remarque: les particules de r�f�rence ont un num�ro g�n�rique -1
        // d'o� "false" dans le constructeur pour auto_numbering = false
        Particule* particuleRef = new Particule( nParticule, false,
            nbPC+int(i) );
        m_composants.AjouterClasseParticules( particuleRef );
        pair<Particule*,int> ppp( particuleRef, nbre );
        m_newParticules.push_back( ppp );

        m_composants.set_NbParticlesPerClass( int(i), nbre );
      }
      m_composants.compute_SauterMeanDiameter();
      m_composants.compute_ParticleClassesConcentration();
    }

    // Particules Composites
    DOMNode* compParticules = ReaderXML::getNode( root, "CompParticules" );

    if ( compParticules )
    {
      bnewpart = true;

      int nbPC  = int(m_composants.getParticuleClassesReference()->size());

      DOMNodeList* allCompParticules = ReaderXML::getNodes( rootElement,
          "CompParticule");

      // nombre de particules elemtaires
      size_t nbreCompPart = allCompParticules->getLength();

      for (XMLSize_t i=0; i<nbreCompPart; i++)
      {
        DOMNode* nCompParticule = allCompParticules->item( i );
        int nbre = ReaderXML::getNodeAttr_Int( nCompParticule, "Nombre" );
        bool nName = ReaderXML::hasNodeAttr_String( nCompParticule, "Name" );
        string name;

        // Remarque: les particules de r�f�rence ont un num�ro g�n�rique -1
        // d'o� "false" dans le constructeur pour auto_numbering = false
        if ( nName )
        {
          name = ReaderXML::getNodeAttr_String(nCompParticule, "Name");
          Particule* particuleRef = new CompParticule( nCompParticule,
              false, nbPC+int(i), name );
          m_composants.AjouterClasseParticules( particuleRef );
          pair<Particule*,int> ppp( particuleRef, nbre );
          m_newParticules.push_back( ppp );
        }
        else
        {
          Particule* particuleRef = new CompParticule( nCompParticule,
	      false, nbPC+int(i) );

          m_composants.AjouterClasseParticules( particuleRef );
          pair<Particule*,int> ppp( particuleRef, nbre );
          m_newParticules.push_back( ppp );
        }
      }
    }

    // Obstacles
    DOMNode* obstacles = ReaderXML::getNode( root, "Obstacles" );
    if ( obstacles )
    {
      bnewobst = true;

      DOMNodeList* allCompObstacles = ReaderXML::getNodes( obstacles );
      for (XMLSize_t i=0; i<allCompObstacles->getLength(); i++)
      {
        DOMNode* nCompObs = allCompObstacles->item( i );
        Obstacle *obstacle = Obstacle_BuilderFactory::create( nCompObs );
        m_composants.Ajouter( obstacle );
      }
    }

    // Cas de simulation Periodique ?
    DOMNode* periodes = ReaderXML::getNode( root, "Periodes" );
    if ( periodes )
    {
      bnewper = false;

      if ( Grains_Exec::m_MPIperiodique )
      {
        if ( m_rank == 0 )
          cout << "Si la periodicite est geree par le pattern MPI, il n'est pas"
	   " possible de definir en plus des obstacles periodiques"
	   << " ( bloc xml <Periodes> </Periodes> )" << endl << endl;
        grainsAbort();
      }
      readPeriodic( rootElement );
    }

    // Contact
    DOMNode* contact = ReaderXML::getNode( root, "Contacts" );
    if ( contact )
      Contact_BuilderFactory::define( contact );
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
  }

  // Validation & affichage informatif
  if ( !brestart && !bnewpart && !bnewobst && !bnewper )
  {
    if ( m_rank == 0 )
      cout << "ERR : Mise en donnees incomplete sur <Contruction>" << endl;
    grainsAbort();
  }

  if ( m_rank == 0 ) cout << "Etude avec :" << endl;
  if ( brestart )
  {
    if ( m_rank == 0 )
    {
      cout << "   Type de restart = " << Grains_Exec::m_ReloadType << endl;
      cout << "   Repertoire de restart = " << Grains_Exec::m_ReloadDirectory
           << endl;
      cout << "   Format du fichier de restart = " << ( m_new_reload_format ?
    	"2014" : "ancien" ) << endl;
      cout << "   Storage of the particle's contact history = "
      << ( m_history_storage ? "yes" : "no" ) << endl;
    }
    cout << "   Restart du fichier " << restart << endl;
  }
  if ( bnewpart )
    if ( m_rank == 0 ) cout << "  Ajout de particules" << endl;
  if ( bnewobst )
    if ( m_rank == 0 ) cout << "  Ajout de parois" << endl;

  if ( m_processorIsActiv )
  {
    // Construction du probleme et affectation des composants
    Scalar rayon = m_composants.getRayonMax();
    if ( rayon < 1e-8 ) grainsAbort();
    else if ( m_rank == 0 )
      cout << endl << "Rayon max des particules = " << rayon << endl;

    // If coehsive contact, particles in relation can be further than dp
    // Also cells can be widened by the user
    Scalar LC_coef = 1.;
    DOMNode* cellsize = ReaderXML::getNode( root, "CellSize" );
    if ( cellsize )
      LC_coef = ReaderXML::getNodeAttr_Double( cellsize, "Factor" );
    if ( LC_coef < 1. ) LC_coef = 1.;
    if( Grains_Exec::m_withCohesion && LC_coef < 1.2 )
      LC_coef *= 1.2;
    defineLinkedCell( LC_coef * rayon );
    m_sec->Link( m_composants.getObstacles() );
    if ( m_rank == 0 )
    {
      cout << "Traitement des contacts particule-obstacle "
    	<< "dans le LinkedCell" << endl;
      cout << endl << "Schema d'integration en temps = " <<
      	Grains_Exec::m_TIScheme << endl << endl;
    }

    // Construction de la liste des periodes
    list<MonObstacle*> allobstacles = m_composants.getObstacles()
    	->getObstacles();
    list<Vecteur const*> allPeriodes ;
    for (list<MonObstacle*>::iterator il=allobstacles.begin();
    	il!=allobstacles.end();il++)
      if ( (*il)->getObstaclePeriodic() )
        allPeriodes.push_back( (*il)->getObstaclePeriodic()->getPeriode() );
    if ( !allPeriodes.empty() )
    {
      if ( m_rank == 0 ) cout << "Periodes" << endl;
      m_periodicVectors.reserve( allPeriodes.size() );
      for (list<Vecteur const*>::iterator iv=allPeriodes.begin();
      	iv!=allPeriodes.end();iv++)
      {
        m_periodicVectors.push_back( *iv );
        cout << "   " << *(*iv);
      }
      if ( m_rank == 0 ) cout << endl;
    }
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Construction des forces actives
void Grains::Forces( DOMElement* rootElement )
{
  double rhoF=0., muF=0.;
  string isLift;

  if ( m_processorIsActiv )
  {
    assert(rootElement != NULL);
    DOMNode* root = ReaderXML::getNode( rootElement, "Forces" );

    if( m_rank == 0 ) cout << "Simulation avec :" << endl;

    // Gravite ?
    DOMNode* nGravite = ReaderXML::getNode( root, "Gravite" );
    if( nGravite )
    {
      Grains_Exec::m_vgravite[X] = ReaderXML::getNodeAttr_Double(
      	nGravite, "GX" );
      Grains_Exec::m_vgravite[Y] = ReaderXML::getNodeAttr_Double(
      	nGravite, "GY" );
      Grains_Exec::m_vgravite[Z] = ReaderXML::getNodeAttr_Double(
      	nGravite, "GZ" );
      if( m_rank == 0 ) cout << "  Gravite" << endl;
    }
    else
    {
      if( m_rank == 0 ) cout << "Gravite obligatoire !!";
      grainsAbort();
    }

    // Masse volumique du fluide
    // si 0 => force = poids des particules
    // si != 0 => Archimede
    DOMNode* nFluidProperties = ReaderXML::getNode( root, "FluidProperties" );
    if( nFluidProperties )
    {
      rhoF = ReaderXML::getNodeAttr_Double( nFluidProperties, "rhoF" );

      if( ReaderXML::hasNodeAttr_String( nFluidProperties, "muF" ) )
        muF = ReaderXML::getNodeAttr_Double( nFluidProperties, "muF" );

//      if( ReaderXML::hasNodeAttr_String( nFluidProperties, "TempF" ) )
//        TempF = ReaderXML::getNodeAttr_Double( nFluidProperties, "TempF" );
      if( m_rank == 0 ) cout << "  Archimede effect" << endl;
    }

    // We keep this paragraph for backward compatibility
    DOMNode* nMasseVolFluide = ReaderXML::getNode( root, "MasseVolFluide" );
    if( nMasseVolFluide )
      rhoF = ReaderXML::getNodeAttr_Double( nMasseVolFluide, "Value" );

    Particule::setFluideMasseVolumique( rhoF ) ;
    Particule::setFluidViscosity( muF ) ;
//    Particule::setFluidInitialTemperature( TempF ) ;

    // Calcul du poids des particules
    m_composants.computeWeight( 0., 0. );

    // Drag Force ?
    DOMNode* nDrag = ReaderXML::getNode( root, "DragForce" );
    if( nDrag )
    {
      if( muF<1.e-12 && m_rank==0 )
      {
        cout << " FATAL ERROR : Drag module requested with muF<1.e-12 " << endl;
        exit(0);
      }
      if ( ReaderXML::hasNodeAttr_String( nDrag, "WithLift" ) )
        isLift = ReaderXML::getNodeAttr_String( nDrag, "WithLift" );
      if( isLift == "yes" )
      {
        Grains_Exec::m_withLiftForce = true;
        if ( m_rank == 0 ) cout << "  Lift force" << endl;
      }

      app_HydroForce = new AppFluide_Drag( nDrag );
      app_HydroForce->setName( "TraineeHydro" );
      m_allApp.push_back( app_HydroForce );
      Grains_Exec::set_listApp( m_allApp );
      Grains_Exec::m_withHydroForce = true;
      if ( m_rank == 0 ) cout << "  Trainee hydrodynamique" << endl;
    }

    // Heat transfert ?
    DOMNode* nTemperature = ReaderXML::getNode( root, "Temperature" );
    if( nTemperature )
    {
      app_FluidTemperature = new AppFluide_Temperature( nTemperature );
      app_FluidTemperature->setName( "FluidTemperature" );
      m_allApp.push_back( app_FluidTemperature );
      if ( m_rank == 0 ) cout << "  Heat Transfer" << endl;
    }

    if ( m_rank == 0 ) cout << endl;
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Construction de la simulation.
void Grains::Chargement( DOMElement* rootElement )
{
  if ( m_processorIsActiv )
  {
    assert(rootElement != NULL);
    DOMNode* root = ReaderXML::getNode( rootElement, "Simulation" );

    // Temps
    DOMNode* temps = ReaderXML::getNode( root, "Temps" );
    m_tdeb = ReaderXML::getNodeAttr_Double( temps, "Debut" );
    m_tfin = ReaderXML::getNodeAttr_Double( temps, "Fin" );
    if ( Grains_Exec::m_ReloadType == "same" ) m_tdeb = m_temps;

    // Pas de temps
    DOMNode* increment = ReaderXML::getNode( root, "IncreTemps" );
    m_dt = ReaderXML::getNodeAttr_Double( increment, "dt" );

    // Sauvegarde
    DOMNode* tsave = ReaderXML::getNode( root, "TempsSave" );
    double debutSave = ReaderXML::getNodeAttr_Double( tsave, "Debut" );
    double finSave   = ReaderXML::getNodeAttr_Double( tsave, "Fin" );
    double dtSave    = ReaderXML::getNodeAttr_Double( tsave, "Dt" );
    if ( dtSave < m_dt ) dtSave = m_dt;
    for (Scalar t=debutSave; t-finSave < 0.01 * m_dt; t+=dtSave)
      m_save.push_back(t);

    // Fichier de reload et mode d'ecriture
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
      cout << "   Mode d'ecriture = " << ( Grains_Exec::m_writingModeHybrid ?
      	"Hybride" : "Texte" ) << endl;
      cout << endl;
    }

    // Insertion & Position | Fenetres ?
    DOMNode* insertion = ReaderXML::getNode( root, "Insertion" );
    if ( insertion )
    {
      string type = ReaderXML::getNodeAttr_String( insertion, "Mode" );
      if ( m_rank == 0 ) cout << "Insertion :" << endl;
      if ( type == "Sequentiel" )
      {
        m_mode_insertion_particules = ORDER;
        if ( m_rank == 0 ) cout << "   * mode sequentiel" << endl;
      }
      else
      {
        m_mode_insertion_particules = RANDOM;
        if ( m_rank == 0 ) cout << "   * mode aleatoire" << endl;
      }

      string methode = ReaderXML::getNodeAttr_String( insertion, "Methode" );
      if ( methode == "Preinstall" )
      {
        m_methode_insertion_particules = BEFORE;
        if ( m_rank == 0 ) cout << "   * Avant le debut" << endl;
      }
      else
      {
        m_methode_insertion_particules = INLINE;
        if ( m_rank == 0 ) cout << "   * En dynamique pendant la simulation"
		<< endl;
      }

      string config = ReaderXML::getNodeAttr_String( insertion,
      		"Configuration" );
      if ( config == "Aleatoire" )
      {
        m_configAleatoire = true;
        if ( m_rank == 0 ) cout << "   * configuration aleatoire" << endl;
      }
      else
      {
        m_configAleatoire = false;
        if ( m_rank == 0 ) cout << "   * configuration fixe" << endl;
      }

      m_insertion_frequency = ReaderXML::getNodeAttr_Int( insertion,
      	"Frequence" );
      if ( m_rank == 0 )
        cout << "   * frequence = " << m_insertion_frequency << endl;

      string random_mode = ReaderXML::getNodeAttr_String( insertion,
      	"Aleatoire" );
      if ( random_mode == "Total" ) srand( (unsigned int)(time(NULL)) );
      else random_mode = "Reproductible";
      if ( m_rank == 0 )
        cout << "   * mode aleatoire = " << random_mode << endl;

      // Position
      DOMNode* nPosition = ReaderXML::getNode( insertion, "Position" );
      if ( nPosition )
      {
        if ( m_methode_insertion_particules == BEFORE )
          m_force_insertion = true;
        m_position = ReaderXML::getNodeValue_String( nPosition );
        if ( m_position == "STRUCTURED" )
        {
          m_blocInsert = new struct BlocInsertion;
          m_blocInsert->bloc.ftype = FENETRE_BOX;
          m_blocInsert->bloc.radius = m_blocInsert->bloc.hauteur = 0. ;
          m_blocInsert->bloc.axisdir = W ;

          DOMNode* nFenetre = ReaderXML::getNode( insertion, "Fenetre" );
              DOMNodeList* points = ReaderXML::getNodes( nFenetre );
          DOMNode* pointA = points->item( 0 );
          DOMNode* pointB = points->item( 1 );
          m_blocInsert->bloc.ptA[X] =
            ReaderXML::getNodeAttr_Double( pointA, "X" );
          m_blocInsert->bloc.ptA[Y] =
            ReaderXML::getNodeAttr_Double( pointA, "Y" );
          m_blocInsert->bloc.ptA[Z] =
            ReaderXML::getNodeAttr_Double( pointA, "Z" );
          m_blocInsert->bloc.ptB[X] =
            ReaderXML::getNodeAttr_Double( pointB, "X" );
          m_blocInsert->bloc.ptB[Y] =
            ReaderXML::getNodeAttr_Double( pointB, "Y" );
          m_blocInsert->bloc.ptB[Z] =
            ReaderXML::getNodeAttr_Double( pointB, "Z" );

          DOMNode* nStruct = ReaderXML::getNode( insertion, "Nombre" );
          m_blocInsert->NX = ReaderXML::getNodeAttr_Int( nStruct, "NX" );
          m_blocInsert->NY = ReaderXML::getNodeAttr_Int( nStruct, "NY" );
          m_blocInsert->NZ = ReaderXML::getNodeAttr_Int( nStruct, "NZ" );

          if ( m_rank == 0 )
          {
            cout << endl;
            cout << "Position en bloc structure" << endl;
            cout << "   * Point min = " << m_blocInsert->bloc.ptA[X]
                << " " << m_blocInsert->bloc.ptA[Y] << " " <<
            m_blocInsert->bloc.ptA[Z] << endl;
            cout << "   * Point max = " << m_blocInsert->bloc.ptB[X]
                << " " << m_blocInsert->bloc.ptB[Y] << " " <<
            m_blocInsert->bloc.ptB[Z] << endl;
            cout << "   * Bloc = " << m_blocInsert->NX
                << " x " << m_blocInsert->NY << " x " <<
            m_blocInsert->NZ << endl;
            cout << endl;
          }
        }
      }


      // Mode d'initialisation de la vitesse
      if ( m_rank == 0 )
        cout << endl << "Initialisation de vitesse des particules : " << endl;

      DOMNode* nodeInitVit = ReaderXML::getNode( insertion, "InitVitesse" );
      if ( nodeInitVit )
      {
	string sInitVitmode =
	    ReaderXML::getNodeAttr_String( nodeInitVit, "Mode" );
	if ( sInitVitmode == "Constant" )
	{
	  m_methode_initvit_particules = CONSTANT;

          DOMNode* nVitTransInit = ReaderXML::getNode( nodeInitVit,
	  	"TranslationVelocity" );
          if ( nVitTransInit )
          {
            m_InitVtrans[X] = ReaderXML::getNodeAttr_Double( nVitTransInit,
	      	"VX" );
            m_InitVtrans[Y] = ReaderXML::getNodeAttr_Double( nVitTransInit,
	      	"VY" );
            m_InitVtrans[Z] = ReaderXML::getNodeAttr_Double( nVitTransInit,
	      	"VZ" );
          }

          DOMNode* nVitRotInit = ReaderXML::getNode( nodeInitVit,
	    	"RotationVelocity" );
          if ( nVitRotInit )
          {
            m_InitVrot[X] = ReaderXML::getNodeAttr_Double( nVitRotInit,
	      	"RX" );
            m_InitVrot[Y] = ReaderXML::getNodeAttr_Double( nVitRotInit,
	      	"RY" );
            m_InitVrot[Z] = ReaderXML::getNodeAttr_Double( nVitRotInit,
	      	"RZ" );
          }
	}
	else if ( sInitVitmode == "RandomMotion" )
	{
	  m_methode_initvit_particules = RANDOM_MOTION;

	  DOMNode* nodeRandomTrans = ReaderXML::getNode( nodeInitVit,
	    	"Translation" );
	  if ( nodeRandomTrans )
	    m_RandomMotionCoefTrans =
	 	ReaderXML::getNodeAttr_Double( nodeRandomTrans, "Amplitude" );

	  DOMNode* nodeRandomRot = ReaderXML::getNode( nodeInitVit,
	    	"Rotation" );
	  if ( nodeRandomRot )
	    m_RandomMotionCoefRot =
	 	ReaderXML::getNodeAttr_Double( nodeRandomRot, "Amplitude" ) ;
	}
	else m_methode_initvit_particules = ZERO;
      }

      if ( m_rank == 0 )
      {
	switch( m_methode_initvit_particules )
        {
          case CONSTANT :
            cout << "   * vitesse initiale constante en translation = ( " <<
	   	 m_InitVtrans[X] << ", " << m_InitVtrans[Y] << ", "
		 << m_InitVtrans[Z] << " )" << endl;
            cout << "   * vitesse initiale constante en rotation = ( " <<
	   	 m_InitVrot[X] << ", " << m_InitVrot[Y] << ", "
		 << m_InitVrot[Z] << " )" << endl;
            break;

          case RANDOM_MOTION :
            cout << "   * vitesse initiale aleatoire en translation = " <<
	   	m_RandomMotionCoefTrans << endl;
            cout << "   * vitesse initiale aleatoire en rotation = " <<
	   	m_RandomMotionCoefRot << endl;
            break;

          default :
            cout << "   * vitesse nulle en translation et rotation" << endl;
            break;
        }
	cout << endl;
      }


      // Fenetres pour insertion aleatoire inline
      DOMNode* nFenetres = ReaderXML::getNode( insertion, "Fenetres" );
      if ( nFenetres )
      {
	DOMNodeList* allFenetres = ReaderXML::getNodes( nFenetres );
        for (XMLSize_t i=0; i<allFenetres->getLength(); i++)
	{
	  DOMNode* nFenetre = allFenetres->item( i );
	  string fenetre_type = "Box";
	  if ( ReaderXML::hasNodeAttr_String( nFenetre, "Type" ) )
	    fenetre_type = ReaderXML::getNodeAttr_String( nFenetre, "Type" );

	  Fenetre fenetre;
	  if ( fenetre_type == "Cylindre" ) fenetre.ftype = FENETRE_CYLINDER;
	  else if ( fenetre_type == "Annulaire" )
	    fenetre.ftype = FENETRE_ANNULUS;
	  else if ( fenetre_type == "Ligne" )
	    fenetre.ftype = FENETRE_LINE;
	  else fenetre.ftype = FENETRE_BOX;

	  DOMNodeList* points = NULL;
	  DOMNode* pointA = NULL;
	  DOMNode* pointB = NULL;
	  DOMNode* cylCarac = NULL;
	  string axisdir_str = "W";

	  switch( fenetre.ftype )
	  {
	    case FENETRE_BOX:
 	      points = ReaderXML::getNodes( nFenetre );
	      pointA = points->item( 0 );
	      pointB = points->item( 1 );

	      fenetre.radius = fenetre.radius_int = fenetre.hauteur = 0. ;
	      fenetre.axisdir = W ;
	      fenetre.ptA[X] = ReaderXML::getNodeAttr_Double( pointA, "X" );
	      fenetre.ptA[Y] = ReaderXML::getNodeAttr_Double( pointA, "Y" );
              fenetre.ptA[Z] = ReaderXML::getNodeAttr_Double( pointA, "Z" );
	      fenetre.ptB[X] = ReaderXML::getNodeAttr_Double( pointB, "X" );
	      fenetre.ptB[Y] = ReaderXML::getNodeAttr_Double( pointB, "Y" );
	      fenetre.ptB[Z] = ReaderXML::getNodeAttr_Double( pointB, "Z" );
	    break;

	    case FENETRE_CYLINDER:
	      pointA = ReaderXML::getNode( nFenetre, "Point" );
	      fenetre.ptA[X] = ReaderXML::getNodeAttr_Double( pointA, "X" );
	      fenetre.ptA[Y] = ReaderXML::getNodeAttr_Double( pointA, "Y" );
              fenetre.ptA[Z] = ReaderXML::getNodeAttr_Double( pointA, "Z" );
	      fenetre.ptB[X] = fenetre.ptB[Y] = fenetre.ptB[Z] = 0.;
	      cylCarac = ReaderXML::getNode( nFenetre, "Dimensions" );
	      fenetre.radius = ReaderXML::getNodeAttr_Double( cylCarac,
			"Radius" );
	      fenetre.radius_int = 0.;
	      fenetre.hauteur = ReaderXML::getNodeAttr_Double( cylCarac,
			"Hauteur" );
	      axisdir_str = ReaderXML::getNodeAttr_String( cylCarac,
	          "Direction" );
	      if ( axisdir_str == "X" ) fenetre.axisdir = X;
	      else if ( axisdir_str == "Y" ) fenetre.axisdir = Y;
	      else if ( axisdir_str == "Z" ) fenetre.axisdir = Z;
	      else
	      {
                if ( m_rank == 0 )
                  cout << "Direction erronee dans la defintion de la fenetre"
                       << " cylindrique; valeurs admissibles: X, Y ou Z"
                       << endl;
                grainsAbort();
              }
              break;

            case FENETRE_ANNULUS:
              pointA = ReaderXML::getNode( nFenetre, "Point" );
              fenetre.ptA[X] = ReaderXML::getNodeAttr_Double( pointA, "X" );
              fenetre.ptA[Y] = ReaderXML::getNodeAttr_Double( pointA, "Y" );
              fenetre.ptA[Z] = ReaderXML::getNodeAttr_Double( pointA, "Z" );
              fenetre.ptB[X] = fenetre.ptB[Y] = fenetre.ptB[Z] = 0.;
              cylCarac = ReaderXML::getNode( nFenetre, "Dimensions" );
              fenetre.radius = ReaderXML::getNodeAttr_Double( cylCarac,
                  "RadiusExt" );
              fenetre.radius_int = ReaderXML::getNodeAttr_Double( cylCarac,
                  "RadiusInt" );
              fenetre.hauteur = ReaderXML::getNodeAttr_Double( cylCarac,
                  "Hauteur" );
              axisdir_str = ReaderXML::getNodeAttr_String( cylCarac,
                  "Direction" );
              if ( axisdir_str == "X" ) fenetre.axisdir = X;
              else if ( axisdir_str == "Y" ) fenetre.axisdir = Y;
              else if ( axisdir_str == "Z" ) fenetre.axisdir = Z;
              else
              {
                if ( m_rank == 0 )
                  cout << "Direction erronee dans la defintion de la fenetre"
                       << " cylindrique; valeurs admissibles: X, Y ou Z"
                       << endl;
                grainsAbort();
              }
              break;

	    case FENETRE_LINE:
 	      points = ReaderXML::getNodes( nFenetre );
	      pointA = points->item( 0 );
	      pointB = points->item( 1 );

	      fenetre.radius = fenetre.radius_int = fenetre.hauteur = 0. ;
	      fenetre.axisdir = W ;
	      fenetre.ptA[X] = ReaderXML::getNodeAttr_Double( pointA, "X" );
	      fenetre.ptA[Y] = ReaderXML::getNodeAttr_Double( pointA, "Y" );
	      fenetre.ptA[Z] = ReaderXML::getNodeAttr_Double( pointA, "Z" );
	      fenetre.ptB[X] = ReaderXML::getNodeAttr_Double( pointB, "X" );
	      fenetre.ptB[Y] = ReaderXML::getNodeAttr_Double( pointB, "Y" );
	      fenetre.ptB[Z] = ReaderXML::getNodeAttr_Double( pointB, "Z" );
	    break;

	    default:
	      if ( m_rank == 0 ) cout << "Type de fenetre inconnue" << endl;
	      grainsAbort();
	      break;
	  }

	  m_fenetres.insert( m_fenetres.begin(), fenetre );
        }
      }
    }


    // Chargements
    DOMNode* nChargements = ReaderXML::getNode( root, "Chargements" );
    if ( nChargements )
    {
      DOMNodeList* allChargements = ReaderXML::getNodes( nChargements );
      for (XMLSize_t i=0; i<allChargements->getLength(); i++)
      {
        DOMNode* nChargement = allChargements->item( i );

	// Chargements en Force
	if ( ReaderXML::getNodeAttr_String( nChargement, "Mode") == "Force" )
	{
	  ObstacleChargement_F* chargement = new ObstacleChargement_F(
	      nChargement, m_dt, m_rank );
	  m_composants.Associer( *chargement );
	}
	else
	{
	  ObstacleChargement* chargement = new ObstacleChargement( nChargement,
		m_dt, m_rank );
	  m_composants.Associer( *chargement );
	}

	if ( ReaderXML::hasNodeAttr_String( nChargement, "GeomTranslation" ) )
	{
	  string translated = ReaderXML::getNodeAttr_String( nChargement,
	    "GeomTranslation" );
	  if ( translated == "False" ) Obstacle::setDeplaceObstacle( false );
	}
      }
    }

    // Post-processing des contraintes
    DOMNode* ppStress = ReaderXML::getNode( root, "Contraintes" );
    if ( ppStress )
    {
      Grains_Exec::m_stressTensor = true;
      DOMNode* nFenetres = ReaderXML::getNode( ppStress, "Fenetres" );
      // Take the domain size if no window
      if ( !nFenetres ) m_constrainteFenetres.clear();
      else
      {
        DOMNodeList* allFenetres = ReaderXML::getNodes( nFenetres );
        for (XMLSize_t i=0; i<allFenetres->getLength(); i++)
        {
          DOMNode* nFenetre = allFenetres->item( i );
          DOMNodeList* points = ReaderXML::getNodes( nFenetre );
          DOMNode* pointA = points->item( 0 );
          DOMNode* pointB = points->item( 1 );
          Fenetre fenetre;
          fenetre.ptA[X] = ReaderXML::getNodeAttr_Double( pointA, "X" );
          fenetre.ptA[Y] = ReaderXML::getNodeAttr_Double( pointA, "Y" );
          fenetre.ptA[Z] = ReaderXML::getNodeAttr_Double( pointA, "Z" );
          fenetre.ptB[X] = ReaderXML::getNodeAttr_Double( pointB, "X" );
          fenetre.ptB[Y] = ReaderXML::getNodeAttr_Double( pointB, "Y" );
          fenetre.ptB[Z] = ReaderXML::getNodeAttr_Double( pointB, "Z" );
          m_constrainteFenetres.insert( m_constrainteFenetres.begin(), fenetre );
        }
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

      // Eventual window where Paraview writes down granular outputs
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

	if( m_nprocs > 1 )
	  synchronize_PPWindow();
      }

    }


    // Frequence de mise a jour du lien entre obstacles et LinkedCell
    DOMNode* nodeUpdateFreq = ReaderXML::getNode( root, "LinkUpdate" );
    int updateFreq = 1;
    if ( nodeUpdateFreq )
      updateFreq = ReaderXML::getNodeAttr_Int( nodeUpdateFreq, "frequence" );
    m_composants.setObstaclesLinkedCellUpdateFrequency( updateFreq );

//    // Deplacement geometrique des obstacles
//    DOMNode* nodeDeplaceObs = ReaderXML::getNode( root,
//      		"DeplacementGeometrique" );
//    if ( nodeDeplaceObs )
//    {
//      string deplobsval = ReaderXML::getNodeAttr_String( nodeDeplaceObs,
//      		"value" );
//      if ( deplobsval == "false" ) Obstacle::setDeplaceObstacle( false ) ;
//    }

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




// ----------------------------------------------------------------------------
// Gestion de la simulation
void Grains::Simulation( bool predict, bool isPredictorCorrector,
	bool explicit_added_mass )
{
  if ( m_processorIsActiv )
  {
    m_temps = m_tdeb;
    list<App*>::iterator app;
    size_t npwait_nm1 = 0;
    bool b_lastTime_save = false;
    bool b_error_occured = false ;
    Scalar vmax = 0., vmean = 0. ;

    // Timers
    CT_set_start();
    SCT_insert_app( "ParticulesInsertion" );
    SCT_set_start( "ParticulesInsertion" );
    SCT_insert_app( "CalculerForces" );
    SCT_insert_app( "ParticulesPeriodiques" );
    SCT_insert_app( "Deplacer" );
    SCT_insert_app( "Actualiser" );
    SCT_insert_app( "LinkUpdate" );
    SCT_insert_app( "InitialisationSortieResultats" );
    SCT_insert_app( "SortieResultats" );

    // Creation, insertion et link des particules
    InsertCreateNewParticules();
    SCT_get_elapsed_time( "ParticulesInsertion" );

    // vector<Cellule*>::const_iterator mycell;
    // for(mycell=m_sec->getAllCellules()->begin();
    //     mycell!=m_sec->getAllCellules()->end(); mycell++)
    //   {
    //     if ((*mycell)->nombreObstacles())
    //     {
    //       cout << (*mycell)->nombreObstacles() << " obstacle in the current cell" << endl;
    //       cout << "obs id3 = " << ((*mycell)->getObstacle()->front())->getID() << endl;
    //     }
    //   }

    // Verifie l'allocation de la structure des infos du fluide pour
    // le calcul de la force de trainee avec le fluide au repos (non resolu)
    // dans les particules deja actives (cas de reload)
    if( Grains_Exec::m_withHydroForce )
      allocateDEMCFD_FluidInfos();

    if( Grains_Exec::m_withFluidTemperature ||
        Grains_Exec::m_withSolidTemperature )
    {
      app_FluidTemperature->InitializeTemperature( 0., m_dt,
        m_composants.getParticulesActives(), 0. );
    }
    SCT_set_start( "InitialisationSortieResultats" );
    // Nombre de particules inserees et nombre total de particules dans le
    // syst�me sur tous les procs
    m_composants.setNbreParticulesOnAllProc( m_composants.nbreParticules() );
    npwait_nm1 = m_composants.nbreParticulesWait();

    // Initialisation de la cinematique des obstacles
    m_composants.setCinematiqueObstacleSansDeplacement( m_temps, m_dt );

    // Dans le cas d'un mouvement initial aleatoire
    m_composants.setRandomMotion( m_RandomMotionCoefTrans,
	m_RandomMotionCoefRot );
    // Par defaut l'etat initial est conserve pour Post-Processing
    m_composants.PostProcessing_start( m_temps, m_dt, m_sec, m_fenetres );
    ofstream fVitMax( (m_fileSave + "_VitesseMaxMean.dat").c_str(), ios::out );
    m_composants.ComputeMaxMeanVelocity( vmax, vmean );
    cout << "Vitesse des composants : max = " << vmax << " moyenne = " <<
    	vmean << endl;
    fVitMax << Grains_Exec::doubleToString( ios::scientific, 6, m_temps )
    	<< "\t" << Grains_Exec::doubleToString( ios::scientific, 6, vmax )
	<< "\t" << Grains_Exec::doubleToString( ios::scientific, 6, vmean )
	<< endl;
    display_used_memory();

    // Post-processing des efforts sur les obstacles
    m_composants.initialiseOutputObstaclesLoadFiles( m_rank, false, m_temps );
    m_composants.outputObstaclesLoad( m_temps, m_dt, false,
      Grains_Exec::m_ReloadType == "same" );

    list<Scalar>::iterator tempsSave = m_save.begin();
    while ( *tempsSave - m_temps < 0.01 * m_dt && tempsSave != m_save.end() )
      tempsSave++;
    SCT_get_elapsed_time( "InitialisationSortieResultats" );
    // If cohesive contact, check for particle initially glued
    if( Grains_Exec::m_withCohesion )
    {
      m_sec->InitializeCohesiveForces( 0., m_dt,
          m_composants.getParticulesActives() );
    }

    // Algorithme de simulation
    // ~~~~~~~~~~~~~~~~~~~~~~~~
    cout << "Time \t TO \tend \tParticules \tIn \tOut" << endl;
    while ( m_tfin - m_temps > 0.01 * m_dt )
    {
      try {
        m_temps += m_dt;
        b_lastTime_save = false;

        // If shrinking choice in Insert.xml is set to 1 it enters this if
        if( m_composants.IsShrinking() ) m_composants.ShrinkingRate( m_temps );
        // Initialisation de l'indicateur de calcul
        // de la transformation avec epaiseur de croute a faux
        SCT_set_start( "CalculerForces" );
        m_composants.InitializeVdWState( m_temps, m_dt );
        SCT_get_elapsed_time( "CalculerForces" );

        // Insertion des particules en attente
        SCT_set_start( "ParticulesInsertion" );
        if (npwait_nm1 != m_composants.nbreParticulesWait())
          cout << "\r                                              "
               << "                 " << flush;
        ostringstream oss;
        oss.width(10);
        oss << left << m_temps;
        cout << '\r' << oss.str() << "  \t" << m_tfin << "\t\t\t"
           << m_composants.nbreParticulesActivesOnProc() << '\t'
           << m_composants.nbreParticulesWait()    << flush;

        npwait_nm1 = m_composants.nbreParticulesWait();
        if ( m_methode_insertion_particules == INLINE )
          insertParticule( m_mode_insertion_particules );
        SCT_get_elapsed_time( "ParticulesInsertion" );

        // Creation/Destruction des clones periodiques
        SCT_set_start( "ParticulesPeriodiques" );
        m_sec->LinkUpdateParticulesPeriodiques( m_temps,
        m_composants.getParticulesActives(),
        m_composants.getParticulesClonesPeriodiques(),
        m_composants.getParticuleClassesReference() );
        SCT_get_elapsed_time( "ParticulesPeriodiques" );

        // Initialisation des maps de contact
	      m_composants.setAllContactMapToFalse();

        // Calcul des forces
        SCT_set_start( "CalculerForces" );
        // Initialisation des torseurs de force
        m_composants.InitializeForces( m_temps, m_dt, true );

        // Calcul des forces exterieures sur les particules actives
        // Forces de volume et de contact
        for (app=m_allApp.begin(); app!=m_allApp.end(); app++)
          (*app)->CalculerForces( m_temps, m_dt,
        m_composants.getParticulesActives() );
        SCT_add_elapsed_time( "CalculerForces" );


        // Calcul des forces de contact sur les clones periodiques
        SCT_set_start( "ParticulesPeriodiques" );
        m_sec->CalculerForcesClonesPeriodiques( m_temps, m_dt,
        m_composants.getParticulesClonesPeriodiques() );
        m_composants.AddForcesFromPeriodicClonesToParticules( m_temps, m_dt );
        SCT_add_elapsed_time( "ParticulesPeriodiques" );


        // Caclul de la temperature en fonction des flux de chaleur
        if( app_FluidTemperature )
          m_composants.ComputeTemperature( m_temps, m_dt );

        // Mise a jour des maps de contact
	      m_composants.updateAllContactMaps();

        // Deplacement des particules en fonction des forces exercees
        SCT_set_start( "Deplacer" );
        m_composants.Deplacer( m_temps, m_dt );
        SCT_get_elapsed_time( "Deplacer" );
        SCT_set_start( "Actualiser" );
        m_composants.Actualiser();
        SCT_get_elapsed_time( "Actualiser" );


        // Mise � jour des clones periodiques
        SCT_set_start( "ParticulesPeriodiques" );
        m_composants.updateClonesPeriodiques( m_sec );
        SCT_add_elapsed_time( "ParticulesPeriodiques" );


        // Actualisation des particules & obstacles dans les cellules
        SCT_set_start( "LinkUpdate" );
        m_sec->LinkUpdate( m_temps, m_dt,
        m_composants.getParticulesActives() );
        SCT_get_elapsed_time( "LinkUpdate" );


        // Post-processing des efforts sur les obstacles
        m_composants.outputObstaclesLoad( m_temps, m_dt );

        // Sauvegarde du temps de simulation eventuel
        if ( *tempsSave - m_temps < 0.01 * m_dt && tempsSave != m_save.end() )
        {
          SCT_set_start( "SortieResultats" );
          m_composants.ComputeMaxMeanVelocity( vmax, vmean );
          cout << endl << "Vitesse des composants : max = " << vmax
               << " moyenne = " << vmean << endl;
          fVitMax << Grains_Exec::doubleToString( ios::scientific, 6, m_temps )
               << "\t" << Grains_Exec::doubleToString( ios::scientific, 6,
               vmax ) << "\t" << Grains_Exec::doubleToString( ios::scientific,
               6, vmean ) << endl;

          display_used_memory();
          saveReload( m_temps );
          m_composants.PostProcessing( m_temps, m_dt, m_sec );
          while ( *tempsSave - m_temps < 0.01 * m_dt
              && tempsSave != m_save.end() ) tempsSave++;
          b_lastTime_save = true;
          SCT_get_elapsed_time( "SortieResultats" );
        }
      }
      catch (ErreurContact &chocCroute)
      {
        // Fin de simulation sur choc
        cout << endl;
        m_composants.PostProcessingErreurComposants( "ErreurContact",
            chocCroute.getComposants() );
        chocCroute.Message( cout );
        b_error_occured = true;
        break;
      }
      catch (ErreurDeplacement &errDeplacement)
      {
        // Fin de simulation sur deplacement trop grand
        cout << endl;
        m_composants.PostProcessingErreurComposants( "ErreurDeplacement",
            errDeplacement.getComposant() );
        errDeplacement.Message(cout);
        b_error_occured = true;
        break;
      }
      catch (ErreurSimulation &errSimulation)
      {
        // Fin de simulation sur erreur
        cout << endl;
        errSimulation.Message(cout);
        b_error_occured = true;
        break;
      }
    }

    ostringstream oss;
    oss.width(10);
    oss << left << m_temps;
    cout << "\r                                              "
         << "                 " << flush;
    cout << '\r' << oss.str() << "  \t" << m_tfin << "\t\t\t"
         << m_composants.nbreParticulesActivesOnProc() << '\t'
         << m_composants.nbreParticulesWait() << endl;
    SCT_set_start( "SortieResultats" );


    if ( !b_lastTime_save )
    {
      // Max vitesse
      m_composants.ComputeMaxMeanVelocity( vmax, vmean );
      cout << "Vitesse des composants : max = " << vmax << " moyenne = "
           << vmean << endl;
      fVitMax << Grains_Exec::doubleToString( ios::scientific, 6, m_temps )
              << "\t" << Grains_Exec::doubleToString( ios::scientific, 6, vmax )
              << "\t" << Grains_Exec::doubleToString( ios::scientific, 6, vmean )
              << endl;

      // Post processing
      cout << endl << "Ecriture des resultats pour post-processing"
           << " au dernier temps par defaut" << endl;
      m_composants.PostProcessing( m_temps, m_dt, m_sec );
      m_composants.PostProcessing_end();

      // Sauvegarde du fichier de fin pour Reload
      if ( !b_error_occured )
      {
        cout << "Sauvegarde du fichier de reload" << endl;
        saveReload( m_temps );
      }
    }
    fVitMax.close();

    // Caracteristiques globales des contacts pendant la simulation
    cout << endl << "Rayon d'interaction minimal = " <<
    	m_composants.getRayonInteractionMin() << endl;
    cout << "Overlap moyen = " << AppSec::getOverlapMean() << endl;
    cout << "Overlap maximal = " << AppSec::getOverlapMax() << endl;
    cout << "Temps de l'overlap maximal = " <<
	AppSec::getTimeOverlapMax() << endl;
    cout << "Nb d'iterations moyen de GJK = " <<
    	AppSec::getNbIterGJKMean() << endl;
    SCT_get_elapsed_time("SortieResultats");


    // Bilan timing
    double cputime = CT_get_elapsed_time();
    cout << endl << "Full problem" << endl;
    write_elapsed_time_smhd(cout,cputime,"Computation time");
    SCT_get_summary( cout, cputime );
  }
}




// ----------------------------------------------------------------------------
// Cree, insert et link les nouvelles particules
void Grains::InsertCreateNewParticules()
{
  // IMPORTANT: quel que soit le systeme etudie, si celui ci contient
  // N particules et M obstacles, les particules sont numerotees de 0 � N-1
  // et les obstacles de N � N+M-1
  // Dans le cas de reload avec insertion supplementaire, cela necessite de
  // renumeroter les obstacles

  int first_obstacle_id;
  int numPartMax = numeroMaxParticules();
  if (numPartMax == -1)
  {
    first_obstacle_id = 0;
    numPartMax = 0;
  }
  else first_obstacle_id = numPartMax+1;
  list< pair<Particule*,int> >::iterator ipart;

  // Construction des nouvelles particules
  Composant::setNbComposantsCrees( numPartMax );
  for (ipart=m_newParticules.begin();ipart!=m_newParticules.end();ipart++)
  {
    int nbre = ipart->second;
    for (int ii=0; ii<nbre; ii++)
    {
      Particule* particule = ipart->first->createCloneCopy();
      m_composants.Ajouter( particule );
    }
  }

  // Renumerotation des obstacles
  int numInitObstacles = first_obstacle_id;
  for (ipart=m_newParticules.begin();ipart!=m_newParticules.end();ipart++)
    numInitObstacles += ipart->second;
  list<MonObstacle*> lisObstaclesPrimaires =
    	m_composants.getObstacles()->getObstacles();
  list<MonObstacle*>::iterator iobs;
  int ObstacleID = numInitObstacles;
  for (iobs=lisObstaclesPrimaires.begin();iobs!=lisObstaclesPrimaires.end();
    	iobs++)
  {
    // m_composants.updateContactMapId( (*iobs)->getID(), ObstacleID );
    (*iobs)->setID( ObstacleID );
    ++ObstacleID;
  }

  // Affectation des composants aux forces de contact
  m_composants.Link( *m_sec );

  // Si position par defaut, on positionne les particules
  if ( m_position != "" )
  {
    // Positions dans un bloc structure
    if ( m_position == "STRUCTURED" )
      setPositionParticulesBloc( m_mode_insertion_particules );
    // Positions dans un fichier
    else setPositionParticulesFichier( m_mode_insertion_particules );
  }

  if ( m_methode_insertion_particules == BEFORE )
  {
    cout << "Particules \tIn \tOut" << endl;
    while ( m_composants.nbreParticulesWait() != 0 )
    {
        insertParticule( m_mode_insertion_particules );
        cout << "\r                                              " << flush;
        cout << "\r\t\t"
	   << m_composants.nbreParticulesActivesOnProc() << '\t'
	   << m_composants.nbreParticulesWait()    << flush;
    }
    cout << endl;
  }
  cout << endl;

  if ( m_rank == 0 )
  {
    cout << "Volume des particules IN  : " << m_composants.getVolumeIn()  << endl
         << "                      OUT : " << m_composants.getVolumeOut() << endl
         << endl;
  }
}




// ----------------------------------------------------------------------------
// Donne un point au hasard dans le plan d'insertion.
Point Grains::getPoint() const
{
  Point P;
  int nFenetre = 0;
  int nbreFenetres = int( m_fenetres.size() );

  // Tirage aleatoire d'une fenetre
  if ( nbreFenetres != 1 )
  {
    Scalar n = Scalar(random()) / Scalar(INT_MAX);
    nFenetre = int( n * nbreFenetres );
    if ( nFenetre == nbreFenetres ) nFenetre--;
  }

  // Tirage aleatoire d'une position dans la fenetre selectionnee
  double r = 0., theta = 0., axiscoor = 0.;
  switch( m_fenetres[nFenetre].ftype )
  {
    case FENETRE_BOX:
      P[X] = m_fenetres[nFenetre].ptA[X]
      	+ ( Scalar(random()) / Scalar(INT_MAX) )
	* ( m_fenetres[nFenetre].ptB[X] - m_fenetres[nFenetre].ptA[X] );
      P[Y] = m_fenetres[nFenetre].ptA[Y]
      	+ ( Scalar(random()) / Scalar(INT_MAX) )
	* ( m_fenetres[nFenetre].ptB[Y] - m_fenetres[nFenetre].ptA[Y] );
      P[Z] = m_fenetres[nFenetre].ptA[Z]
      	+ ( Scalar(random()) / Scalar(INT_MAX) )
	* ( m_fenetres[nFenetre].ptB[Z] - m_fenetres[nFenetre].ptA[Z] );
      break;

    case FENETRE_CYLINDER:
      r = ( Scalar(random()) / Scalar(INT_MAX) ) * m_fenetres[nFenetre].radius;
      theta = ( Scalar(random()) / Scalar(INT_MAX) ) * 2. * PI;
      axiscoor = ( Scalar(random()) / Scalar(INT_MAX) )
      	* m_fenetres[nFenetre].hauteur;
      switch ( m_fenetres[nFenetre].axisdir )
      {
        case X:
	  P[X] = m_fenetres[nFenetre].ptA[X] + axiscoor;
	  P[Y] = m_fenetres[nFenetre].ptA[Y] + r * cos( theta );
	  P[Z] = m_fenetres[nFenetre].ptA[Z] + r * sin( theta );
	  break;

        case Y:
	  P[X] = m_fenetres[nFenetre].ptA[X] + r * sin( theta );
	  P[Y] = m_fenetres[nFenetre].ptA[Y] + axiscoor;
	  P[Z] = m_fenetres[nFenetre].ptA[Z] + r * cos( theta );
	  break;

        default:
	  P[X] = m_fenetres[nFenetre].ptA[X] + r * cos( theta );
	  P[Y] = m_fenetres[nFenetre].ptA[Y] + r * sin( theta );
	  P[Z] = m_fenetres[nFenetre].ptA[Z] + axiscoor;
	  break;
      }
      break;

    case FENETRE_ANNULUS:
      r = m_fenetres[nFenetre].radius_int
      	+ ( Scalar(random()) / Scalar(INT_MAX) )
      	* ( m_fenetres[nFenetre].radius - m_fenetres[nFenetre].radius_int );
      theta = ( Scalar(random()) / Scalar(INT_MAX) ) * 2. * PI;
      axiscoor = ( Scalar(random()) / Scalar(INT_MAX) )
      	* m_fenetres[nFenetre].hauteur;
      switch ( m_fenetres[nFenetre].axisdir )
      {
        case X:
	  P[X] = m_fenetres[nFenetre].ptA[X] + axiscoor;
	  P[Y] = m_fenetres[nFenetre].ptA[Y] + r * cos( theta );
	  P[Z] = m_fenetres[nFenetre].ptA[Z] + r * sin( theta );
	  break;

        case Y:
	  P[X] = m_fenetres[nFenetre].ptA[X] + r * sin( theta );
	  P[Y] = m_fenetres[nFenetre].ptA[Y] + axiscoor;
	  P[Z] = m_fenetres[nFenetre].ptA[Z] + r * cos( theta );
	  break;

        default:
	  P[X] = m_fenetres[nFenetre].ptA[X] + r * cos( theta );
	  P[Y] = m_fenetres[nFenetre].ptA[Y] + r * sin( theta );
	  P[Z] = m_fenetres[nFenetre].ptA[Z] + axiscoor;
	  break;
      }
      break;

    case FENETRE_LINE:
      P = m_fenetres[nFenetre].ptA + ( Scalar(random()) / Scalar(INT_MAX) )
      	* ( m_fenetres[nFenetre].ptB - m_fenetres[nFenetre].ptA );
      break;

    default:
      break;
  }

  return P;
}




// ----------------------------------------------------------------------------
// Selection aleatoire d'une transformation de type rotation
Matrix Grains::getRandomRotation() const
{
  Matrix rotation;

  Scalar angleZ = 2. * PI * (Scalar)rand() / RAND_MAX;
  Matrix rZ( cos(angleZ), -sin(angleZ), 0.,
  	sin(angleZ), cos(angleZ), 0.,
	0., 0., 1. );

  if ( m_dimension == 3 )
  {
    Scalar angleX = 2. * PI * (Scalar)rand() / RAND_MAX;
    Matrix rX( 1., 0., 0.,
   	0., cos(angleX), -sin(angleX),
	0., sin(angleX), cos(angleX) );

    Scalar angleY = 2. * PI * (Scalar)rand() / RAND_MAX;
    Matrix rY( cos(angleY), 0., sin(angleY),
   	0., 1., 0.,
	-sin(angleY), 0., cos(angleY) );
    Matrix tmp = rY * rZ;
    rotation = rX * tmp;
  }
  else rotation = rZ;

  return rotation;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Insertion d'une particule dans les algorithmes
bool Grains::insertParticule( const PullMode& mode )
{
  static int insert_counter = 0;
  bool contact = true;
  Vecteur vtrans, vrot ;

  if ( insert_counter == 0 )
  {
    Particule *particule = m_composants.getParticule( mode );
    if ( particule )
    {
      // Initialisation de la vitesse de la particule
      computeInitVit( vtrans, vrot );
      particule->setVitesseTranslation( vtrans );
      particule->setVitesseRotation( vrot );

      // Initialisation de la position du centre de gravite de la particule
      if ( m_position == "" )
      {
        Point position = getPoint();
        particule->setPosition( position );
      }

      // Initialisation de la position angulaire de la particule
      if ( m_configAleatoire )
      {
        Transform trot;
        trot.setBasis( getRandomRotation() );
        particule->composePosition( trot );
      }

      particule->initializeVdWtransform_to_notComputed();
      contact = m_sec->isContactVdW( particule );

      // Si plusieurs periodes, on suppose que le nb est <= 3
      // (multi-periodicite 3D max) mais dans Grains-sequentiel, on cree
      // les obstacles periodiques et leur clone, donc on a au max 6 vecteurs
      // periodiques qui vont toujours par paire v,-v
      if ( !m_periodicVectors.empty() && !contact )
      {
        Point refpos = *(particule->getPosition());
        Vecteur pervec ;
        for (size_t i=0;i<m_periodicVectors.size() && !contact;++i)
        {
          particule->setPosition( refpos + *(m_periodicVectors[i]) ) ;
          particule->initializeVdWtransform_to_notComputed();
          contact = m_sec->isContactVdW( particule );
        }

        if ( m_periodicVectors.size() >= 3 && !contact )
        {
          for (size_t ii=0;ii<3 && !contact;++ii)
          {
            pervec = ii < 2 ? *(m_periodicVectors[ii]) : VecteurNul;
            for (size_t jj=0;jj<3 && !contact;++jj)
            {
              pervec += jj < 2 ? *(m_periodicVectors[jj+2]) : VecteurNul ;
              particule->setPosition( refpos + pervec ) ;
              particule->initializeVdWtransform_to_notComputed();
              contact = m_sec->isContactVdW( particule );
            }

            if ( m_periodicVectors.size() == 6 && !contact )
            {
              for (size_t kk=0;kk<3 && !contact;++kk)
              {
                pervec += kk < 2 ? *(m_periodicVectors[kk+4]) : VecteurNul ;
                particule->setPosition( refpos + pervec ) ;
                particule->initializeVdWtransform_to_notComputed();
                contact = m_sec->isContactVdW( particule );
              }
            }
          }
        }
        particule->setPosition( refpos ) ;
        particule->initializeVdWtransform_to_notComputed();
      }

      if ( !contact || m_force_insertion )
      {
        // Lie la particule dans le LinkCell
        m_sec->Link( particule );

        // Switch la particule de la liste d'attente a la liste active
        m_composants.ShiftParticuleOutIn();

        // Verifie l'allocation de la structure des infos du fluide pour
        // le calcul de la force de trainee avec le fluide au repos
        // (fluide non resolu)
        if ( Grains_Exec::m_withHydroForce )
          particule->allocateDEMCFD_FluidInfos();
      }
    }
  }

  ++insert_counter;
  if ( insert_counter == m_insertion_frequency ) insert_counter = 0;

  return( !contact || m_force_insertion );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie la vitesse initiale d'une particule a partir du mode d'initialisation
void Grains::computeInitVit( Vecteur& vtrans, Vecteur& vrot ) const
{
  switch( m_methode_initvit_particules )
  {
    case ZERO :
      vtrans.reset();
      vrot.reset();
      break;

    case CONSTANT :
      vtrans = m_InitVtrans;
      vrot = m_InitVrot;
      break;

    case RANDOM_MOTION :
      if ( m_RandomMotionCoefTrans )
      {
        vtrans[X] = m_RandomMotionCoefTrans * (
	    	2. * (double(random()) / double(INT_MAX)) - 1. ) ;
        vtrans[Y] = m_RandomMotionCoefTrans * (
	    	2. * (double(random()) / double(INT_MAX)) - 1. ) ;
        if ( m_dimension == 3 )
          vtrans[Z] = m_RandomMotionCoefTrans * (
	      	2. * (double(random()) / double(INT_MAX)) - 1. ) ;
	else vtrans[Z] = 0.;
      }
      else vtrans.reset();

      if ( m_RandomMotionCoefRot )
      {
	vrot[X] = m_RandomMotionCoefRot * (
	    	2. * (double(random()) / double(INT_MAX)) - 1. ) ;
	vrot[Y] = m_RandomMotionCoefRot * (
	    	2. * (double(random()) / double(INT_MAX)) - 1. ) ;
	if ( m_dimension == 3 )
	  vrot[Z] = m_RandomMotionCoefRot * (
	      	2. * (double(random()) / double(INT_MAX)) - 1. ) ;
	else vrot[Z] = 0.;
      }
      else vrot.reset();
      break;

    default :
      vtrans.reset();
      vrot.reset();
      break;
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Positionne les particules en attente pour la simulation a partir
// d'un fichier de positions
void Grains::setPositionParticulesFichier( const PullMode& mode )
{
  ifstream filePos( m_position.c_str() );

  // Verification de l'existence du fichier
  if ( !filePos.is_open() )
  {
    cout << "ERR : Fichier absent " << m_position << endl;
    grainsAbort();
  }

  list<Particule*> copie_ParticulesWait;
  list<Particule*>* particules = m_composants.getParticulesWait();
  list<Particule*>::iterator particule = particules->begin();
  Particule* selected=NULL;
  int id,i;
  double v;
  Point position;
  size_t nwait = particules->size(), npos = 0;

  // Verifie que le nb de positions est egal au nombre de particules a inserer
  while ( filePos >> position ) ++npos;
  if ( nwait != npos )
  {
    cout << "ERR: nombre de particules a inserer est different du"
    	<< " nombre de positions dans le fichier " << m_position << " : "
	<< nwait << " != " << npos << endl;
    grainsAbort();
  }
  filePos.close();

  // Lit et affecte la position des particules
  filePos.open( m_position.c_str() );
  if ( mode == RANDOM ) copie_ParticulesWait = *particules;
  while ( filePos >> position )
  {
    switch (mode)
    {
      case ORDER:
        (*particule)->setPosition( position );
        particule++;
        break;

      case RANDOM:
        v = double(random()) / double(INT_MAX);
        id = int( double(copie_ParticulesWait.size()) * v );
        particule = copie_ParticulesWait.begin();
        for (i=0; i<id && particule!=copie_ParticulesWait.end();
          i++, particule++) {}
        selected = *particule;
        particule = copie_ParticulesWait.erase( particule );
        selected->setPosition( position );
        break;

      case NONE:
        break;
    }
  }
  filePos.close();
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Positionne les particules en attente pour la simulation a partir
// d'un bloc structure de positions
void Grains::setPositionParticulesBloc( const PullMode& mode )
{
  list<Particule*> copie_ParticulesWait;
  list<Particule*>* particules = m_composants.getParticulesWait();
  list<Particule*>::iterator particule = particules->begin();
  Particule* selected=NULL;
  int id, i;
  size_t k, l, m, nwait = particules->size();
  double v;
  Point position;
  if ( mode == RANDOM ) copie_ParticulesWait = *particules;

  double deltax = ( m_blocInsert->bloc.ptB[X] - m_blocInsert->bloc.ptA[X] )
  	/ double(m_blocInsert->NX) ;
  double deltay = ( m_blocInsert->bloc.ptB[Y] - m_blocInsert->bloc.ptA[Y] )
  	/ double(m_blocInsert->NY) ;
  double deltaz = ( m_blocInsert->bloc.ptB[Z] - m_blocInsert->bloc.ptA[Z] )
  	/ double(m_blocInsert->NZ) ;

  // Verifie que le nb de positions est egal au nombre de particules a inserer
  if ( nwait != m_blocInsert->NX * m_blocInsert->NY * m_blocInsert->NZ )
  {
    cout << "ERR: nombre de particules a inserer est different du"
    	<< " nombre de positions du bloc structure: " << nwait << " != " <<
	m_blocInsert->NX * m_blocInsert->NY * m_blocInsert->NZ << endl;
    grainsAbort();
  }

  // Affecte la position des particules
  for (k=0;k<m_blocInsert->NX;++k)
    for (l=0;l<m_blocInsert->NY;++l)
      for (m=0;m<m_blocInsert->NZ;++m)
      {
        position[X] = m_blocInsert->bloc.ptA[X] + ( double(k) + 0.5 ) * deltax;
        position[Y] = m_blocInsert->bloc.ptA[Y] + ( double(l) + 0.5 ) * deltay;
        position[Z] = m_blocInsert->bloc.ptA[Z] + ( double(m) + 0.5 ) * deltaz;
	switch (mode)
        {
          case ORDER:
            (*particule)->setPosition( position );
            particule++;
            break;

          case RANDOM:
            v = double(random()) / double(INT_MAX);
            id = int( double(copie_ParticulesWait.size()) * v );
            particule = copie_ParticulesWait.begin();
            for (i=0; i<id && particule!=copie_ParticulesWait.end();
          	i++, particule++) {}
            selected = *particule;
            particule = copie_ParticulesWait.erase( particule );
            selected->setPosition( position );
            break;

          case NONE:
            break;
       }
     }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Schema de resolution: predicteur ou correcteur
// A.WACHS - Janvier 2009 - Creation
void Grains::setMode( const bool &predictor )
{
  m_predictor_mode = predictor;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Schema de resolution: predicteur ou correcteur
// A.WACHS - Janvier 2009 - Creation
bool Grains::isModePredictor()
{
  return m_predictor_mode;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Lecture des donnees pour les simulations MPI
void Grains::readDomainDecomposition( DOMNode* root,
  	const Scalar& lx, const Scalar& ly, const Scalar& lz)
{
  int nprocs[3] = { 1, 1, 1 };
  int coords[3] = { 0, 0, 0 };
  App::setDlocale( lx, ly, lz );
  App::setOriginelocale( nprocs, coords );
  m_processorIsActiv = true;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Renvoie le nom complet du fichier de resultat .resul
string Grains::fullResultFileName( const string &rootname ) const
{
  string fullname = rootname;
  fullname += ".result";

  return fullname;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Lecture des donnees pour les simulations periodiques
void Grains::readPeriodic(DOMElement* rootElement)
{
  DOMNodeList* allPeriodes = ReaderXML::getNodes( rootElement, "Periode" );
  for (XMLSize_t i=0; i<allPeriodes->getLength(); i++)
  {
    DOMNode* nPeriode = allPeriodes->item( i );
    Obstacle *obstacle = Obstacle_BuilderFactory::create( nPeriode );
    m_composants.Ajouter( obstacle );
  }
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Definition du LinkedCell
void Grains::defineLinkedCell( Scalar const& rayon )
{
  m_sec->define( 2. * rayon );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Arret du code en cas de probleme
void Grains::grainsAbort() const
{
  exit(1);
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Sauvegarde des fichiers de reload
void Grains::saveReload( Scalar const& temps )
{
  static size_t reload_counter =
  	Grains_Exec::m_reloadFile_suffix == "A" ? 1 : 0 ;
  string reload_suffix = reload_counter ? "B" : "A" ;
  string reload = m_fileSave + reload_suffix ;
  reload = fullResultFileName( reload );

  // Save current reload file
  ofstream result( reload.c_str() );
  // The contact history is stored regardless of the contact law considered (if
  // the contact law is without history, it only adds one int per particle)
  result << "NEW_RELOAD_FORMAT_AND_CONTACT_HISTORY 2014" << endl;
  result << "#Temps\t" << temps << endl;
  Contact_BuilderFactory::save( result, m_fileSave, m_rank );
  m_composants.write( result, reload.c_str() );
  result << "#EndTemps" << endl;
  result.close();

  // Check that all reload files are in the same directory
  // Add one linen to RFTable.txt
  if ( m_rank == 0 )
  {
    Grains_Exec::checkAllFilesForReload();
    ofstream RFT_out( ( m_fileSave + "_RFTable.txt" ).c_str(), ios::app );
    RFT_out << temps << " " << m_fileSave + reload_suffix << endl;
    RFT_out.close();
  }

  ++reload_counter;
  if ( reload_counter == 2 ) reload_counter = 0;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Numero maximum de particules
// A.WACHS - Mars.2010 - Creation
int Grains::numeroMaxParticules() const
{
  return ( m_composants.numeroMaxParticules() );
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Efface les fichiers .result et .xml */
void Grains::clearResultXmlFiles() const
{
  if ( m_rank == 0 )
  {
    string cmd = "bash " + Grains_Exec::m_GRAINS_HOME
     	+ "/ExecScripts/init_clear.exec " + m_fileSave;
    system( cmd.c_str() );
  }

}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Ecriture de la memoire utilisee par la simulation
void Grains::display_used_memory() const
{
  cout << "Memoire utilisee par Grains3D = ";
  Grains_Exec::display_memory( cout, Grains_Exec::used_memory() );
  cout << endl;
}




// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Allocate DEM-CFD_FluidInfos structure
void Grains::allocateDEMCFD_FluidInfos()
{
//   Grains_Exec::m_withdemcfd = true;
//   Particule::setMassCorrection( false );
  list<Particule*>* allParticles = m_composants.getParticulesActives();
  list<Particule*>::iterator il;
  for( il=allParticles->begin(); il!=allParticles->end(); il++ )
    (*il)->allocateDEMCFD_FluidInfos();
}
