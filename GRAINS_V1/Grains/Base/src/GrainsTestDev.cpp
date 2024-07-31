#include "MPINeighbors.hh"
#include "GrainsTestDev.hh"
#include "GrainsBuilderFactory.hh"
#include "ObstacleBuilderFactory.hh"
#include "ContactBuilderFactory.hh"
#include "PostProcessingWriter.hh"
#include "ParaviewPostProcessingWriter.hh"
#include "STLObstacle.hh"
#include "ConvexBuilderFactory.hh"
#include "PostProcessingWriterBuilderFactory.hh"
#include "RigidBody.hh"
#include "Sphere.hh"
#include "Rectangle.hh"
#include <stdlib.h>
#include <time.h>
#include <random>
#include <thread>

#include <chrono>
#include "Box.hh"
#include "Cylinder.hh"
#include "Superquadric.hh"
#include "PointContact.hh"
#include "Particle.hh"
#include "KinematicsBuilderFactory.hh"


using namespace std;
using namespace std::chrono;


// ----------------------------------------------------------------------------
// Default constructor
GrainsTestDev::GrainsTestDev() 
  : Grains()
{}




// ----------------------------------------------------------------------------
// Destructor
GrainsTestDev::~GrainsTestDev()
{}




// ----------------------------------------------------------------------------
// Writes an initial message in the standard output only on the process ranked 0
void GrainsTestDev::initialOutputMessage()
{
  int rankproc = 0;
  MPI_Comm_rank( MPI_COMM_WORLD, &rankproc );
  if ( !rankproc )
    cout << "Grains3D Test - Development - For developers only" << endl;
}



// ----------------------------------------------------------------------------
// Tasks to perform before time-stepping 
void GrainsTestDev::do_before_time_stepping( DOMElement* rootElement )
{
  Construction( rootElement );
  Forces( rootElement );
  AdditionalFeatures( rootElement );             
  
  cout << "Entering DBTS AM" << endl;

  // Initialisation
  double vmax = 0., vmean = 0. ;

  // Timers
  CT_set_start();
  SCT_insert_app( "Initialization" );
  SCT_set_start( "Initialization" );
  cout << endl << "Initialization" << endl;

  // Set time to initial time
  m_time = m_tstart;


  // Particle creation, insertion and link to grid
  InsertCreateNewParticles();

  // Number of particles: inserted and in the system
  m_allcomponents.computeNumberParticles( m_wrapper );
  m_npwait_nm1 = m_allcomponents.getNumberPhysicalParticlesToInsert();

  // Initialisation obstacle kinematics
  m_allcomponents.setKinematicsObstacleWithoutMoving( m_time, m_dt );

  // In case of initial random motion
  if ( m_initvit_mode == IV_RANDOM )
    m_allcomponents.setRandomMotion( m_RandomMotionCoefTrans,
	m_RandomMotionCoefRot );

  // Writing results for postprocessing
  m_allcomponents.PostProcessing_start( m_time, m_dt, m_collision,
  	m_insertion_windows );

  // Track component max and mean velocity
  fVitMax.open( (m_fileSave + "_VelocityMaxMean.dat").c_str(), ios::out );
  m_allcomponents.ComputeMaxMeanVelocity( vmax, vmean );
  cout << "Component velocity : max = " << vmax << " average = " <<
    	vmean << endl;
  fVitMax << GrainsExec::doubleToString( ios::scientific, 6, m_time )
    	<< "\t" << GrainsExec::doubleToString( ios::scientific, 6, vmax )
	<< "\t" << GrainsExec::doubleToString( ios::scientific, 6, vmean )
	<< endl;

  // Display memory used by Grains
  display_used_memory();

  // Postprocessing of force & torque on obstacles
  m_allcomponents.initialiseOutputObstaclesLoadFiles( m_rank, false, m_time );
  m_allcomponents.outputObstaclesLoad( m_time, m_dt, false,
      GrainsExec::m_ReloadType == "same", m_rank );

  // Next time of writing results
  m_timeSave = m_save.begin();
  while ( *m_timeSave - m_time < 0.01 * m_dt && m_timeSave != m_save.end() )
      m_timeSave++;

  cout << "Initialization completed" << endl << endl;
  SCT_get_elapsed_time( "Initialization" );
       
}






// ----------------------------------------------------------------------------
// Tasks to perform before time-stepping 
void GrainsTestDev::do_after_time_stepping()
{
  double vmax = 0., vmean = 0. ;

  // Particles in & out at the end of the simulation
  ostringstream oss;
  oss.width(10);
  oss << left << m_time;
  cout << "\r                                              "
         << "                 " << flush;
  cout << '\r' << oss.str() << "  \t" << m_tend << "\t\t\t"
         << m_allcomponents.getNumberActiveParticlesOnProc() << '\t'
         << m_allcomponents.getNumberPhysicalParticlesToInsert() << endl;

  // Write reload files
  if ( !m_lastTime_save )
  {
    SCT_insert_app( "OutputResults" );
    SCT_set_start( "OutputResults" );

    // Track component max and mean velocity at the end of the simulation
    m_allcomponents.ComputeMaxMeanVelocity( vmax, vmean );
    cout << endl << "Component velocity : max = " << vmax
	<< " average = " << vmean << endl;
    fVitMax << GrainsExec::doubleToString( ios::scientific, 6, m_time )
	<< "\t" << GrainsExec::doubleToString( ios::scientific, 6,
	vmax ) << "\t" << GrainsExec::doubleToString( ios::scientific,
	6, vmean ) << endl;
    fVitMax.close();

    // Reload files are alsways written at the end of the simulation
    if ( !m_error_occured ) saveReload( m_time );

    SCT_get_elapsed_time( "OutputResults" );
  }

  // Final tasks performed by postprocessing writers
  m_allcomponents.PostProcessing_end();

  // Contact features over the simulation
  cout << endl << "Contact features over the simulation" << endl;
  cout << GrainsExec::m_shift3 << "Minimal crust thickness = " <<
    	m_allcomponents.getCrustThicknessMin() << endl;
  cout << GrainsExec::m_shift3 << "Average overlap = " <<
  	m_collision->getOverlapMean() << endl;
  cout << GrainsExec::m_shift3 << "Maximum overlap = " <<
  	m_collision->getOverlapMax() << endl;
  cout << GrainsExec::m_shift3 << "Time of maximum overlap = " <<
	m_collision->getTimeOverlapMax() << endl;
  cout << GrainsExec::m_shift3 << "Average number of iterations of GJK = " <<
    	m_collision->getNbIterGJKMean() << endl;

  // Timer outcome
  double cputime = CT_get_elapsed_time();
  cout << endl << "Full problem" << endl;
  write_elapsed_time_smhd(cout,cputime,"Computation time");
  SCT_get_summary( cout, cputime );
}







// ----------------------------------------------------------------------------
// Construction of the simulation: linked cell, particles & obstacles, domain 
// decomposition
void GrainsTestDev::Construction( DOMElement* rootElement )
{
  assert( rootElement != NULL );
  DOMNode* root = ReaderXML::getNode( rootElement, "Construction" );

  bool bnewpart = false, bnewobst = false, b2024 = false;
  string restart;
  size_t npart;

  // Domain size: origin, max coordinates and periodicity
  DOMNode* domain = ReaderXML::getNode( root, "LinkedCell" );
  double mx = ReaderXML::getNodeAttr_Double( domain, "MX" );
  double my = ReaderXML::getNodeAttr_Double( domain, "MY" );
  double mz = ReaderXML::getNodeAttr_Double( domain, "MZ" );

  DOMNode* domain_origin = ReaderXML::getNode( root, "Origin" );
  double ox = 0., oy = 0., oz = 0. ;
  if ( domain_origin )
  {
    ox = ReaderXML::getNodeAttr_Double( domain_origin, "OX" );
    oy = ReaderXML::getNodeAttr_Double( domain_origin, "OY" );
    oz = ReaderXML::getNodeAttr_Double( domain_origin, "OZ" );
  }
  App::set_dimensions( mx, my, mz, ox, oy, oz );

  m_periodicity.reserve(3);
  for (size_t i=0;i<3;++i) m_periodicity.push_back(false);
  int perx = 0, pery = 0, perz = 0;
  DOMNode* nPeriodicity = ReaderXML::getNode( root, "Periodicity" );
  if ( nPeriodicity )
  {
    perx = ReaderXML::getNodeAttr_Int( nPeriodicity, "PX" );
    if ( perx != 1 ) perx = 0;
    m_periodicity[X] = perx;
    pery = ReaderXML::getNodeAttr_Int( nPeriodicity, "PY" );
    if ( pery != 1 ) pery = 0;
    m_periodicity[Y] = pery;
    if ( m_dimension == 3 )
      perz = ReaderXML::getNodeAttr_Int( nPeriodicity, "PZ" );
    if ( perz != 1 ) perz = 0;
    m_periodicity[Z] = perz;
  }
  if ( perx || pery || perz )
    GrainsExec::m_periodic = m_periodic = true;
  App::set_periodicity( m_periodicity );


  // Domain decomposition
  readDomainDecomposition( root, mx - ox, my - oy, mz - oz );


  // Display domain size
  if ( m_rank == 0 )
  {
    cout << GrainsExec::m_shift3 << "Domain size" << endl;
    App::output_domain_features( cout, GrainsExec::m_shift6 );
  }


  // Construction on active processes
  if ( m_processorIsActive )
  {
    // Output message
    if ( m_rank == 0 ) cout << GrainsExec::m_shift3 << "Construction" << endl;


    // Create the LinkedCell collision detection app
    m_collision = new LinkedCell();
    m_collision->setName( "LinkedCell" );
    m_allApp.push_back( m_collision );


    // Reload
    DOMNode* reload = ReaderXML::getNode( root, "Reload" );
    if ( reload )
    {
      m_restart = true;

      // Restart mode
      string reload_type = ReaderXML::getNodeAttr_String( reload, "Type" );
      if ( reload_type == "new" || reload_type == "same" )
        GrainsExec::m_ReloadType = reload_type ;


      // Reload file name depending on the restart mode
      // If the mode is "same", the restart file is the same as the output file
      // and is determined by searching the RFTable file
      if ( GrainsExec::m_ReloadType == "new" )
        restart  = ReaderXML::getNodeAttr_String( reload, "Filename" );
      else
      {
        DOMNode* rootSimu = ReaderXML::getNode( rootElement, "Simulation" );
        DOMNode* fileRestartOutput = ReaderXML::getNode( rootSimu,
		"RestartFile" );
        restart = ReaderXML::getNodeAttr_String( fileRestartOutput, "Name" );
        restart = GrainsExec::restartFileName_AorB( restart, "_RFTable.txt" );
        GrainsExec::m_reloadFile_suffix =
            restart.substr( restart.size()-1, 1 );
      }
      restart = GrainsExec::fullResultFileName( restart );

      // Extract the reload directory from the reload file
      GrainsExec::m_ReloadDirectory = GrainsExec::extractRoot( restart );

      // Read the reload file and check the restart format
      string cle;
      ifstream simulLoad( restart.c_str() );
      simulLoad >> cle; 
      if ( cle == "__Format2024__" ) 
      { 
        b2024 = true;
        simulLoad >> cle >> m_time;
      }
      else simulLoad >> m_time;         
      ContactBuilderFactory::reload( simulLoad );
      if ( !b2024 )
      {
        m_allcomponents.read_pre2024( simulLoad, restart, m_wrapper );
        ContactBuilderFactory::set_materialsForObstaclesOnly_reload(
          m_allcomponents.getReferenceParticles() );
      }
      else
        npart = m_allcomponents.read( simulLoad, m_insertion_position, 
		m_rank, m_nprocs );      
      simulLoad >> cle;
      simulLoad.close();

      // Whether to reset velocity to 0
      string reset = ReaderXML::getNodeAttr_String( reload, "Velocity" );
      m_allcomponents.resetKinematics( reset );
    }


    // Particles
    DOMNode* particles = ReaderXML::getNode( root, "Particles" );
    if ( particles )
    {
      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 <<
      	"Reading new particle types" << endl;
      bnewpart = true;
      int nbPC = int( m_allcomponents.getReferenceParticles()->size() );

      DOMNodeList* allParticles = ReaderXML::getNodes( rootElement,
      	"Particle" );

      for (XMLSize_t i=0; i<allParticles->getLength(); i++)
      {
        DOMNode* nParticle = allParticles->item( i );
        size_t nb = ReaderXML::getNodeAttr_Int( nParticle, "Number" );

        // Remark: reference particles' ID number is -1, which explains
        // auto_numbering = false in the constructor
        Particle* particleRef = new Particle( nParticle, nbPC+int(i) );
        m_allcomponents.AddReferenceParticle( particleRef, nb );
        pair<Particle*,size_t> ppp( particleRef, nb );
        m_newParticles.push_back( ppp );
      }

      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 <<
      	"Reading new particle types completed" << endl;
    }

    // Composite particles
    DOMNode* nCompositeParticles =
    	ReaderXML::getNode( root, "CompositeParticles" );
    if ( nCompositeParticles )
    {
      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 <<
      	"Reading new composite particle types" << endl;
      bnewpart = true;
      int nbPC  = int(m_allcomponents.getReferenceParticles()->size());

      DOMNodeList* allCompParticles = ReaderXML::getNodes( rootElement,
          "CompositeParticle");

      for (XMLSize_t i=0; i<allCompParticles->getLength(); i++)
      {
        DOMNode* nCompParticle = allCompParticles->item( i );
        size_t nb = ReaderXML::getNodeAttr_Int( nCompParticle, "Number" );

        // Remark: reference particles' ID number is -1, which explains
        // auto_numbering = false in the constructor
        Particle* particleRef = NULL;
	string sshape = "none";
	if ( ReaderXML::hasNodeAttr( nCompParticle, "SpecificShape" )  )
	  sshape = ReaderXML::getNodeAttr_String( nCompParticle, 
	  	"SpecificShape" );
	if ( sshape == "SpheroCylinder" )
	  particleRef = new SpheroCylinder( nCompParticle, nbPC+int(i) );
	else 	
	  particleRef = new CompositeParticle( nCompParticle, nbPC+int(i) );
        m_allcomponents.AddReferenceParticle( particleRef, nb );
        pair<Particle*,size_t> ppp( particleRef, nb );
        m_newParticles.push_back( ppp );
      }

      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 <<
      	"Reading new composite particle types completed" << endl;
    }


    // Obstacles
    DOMNode* obstacles = ReaderXML::getNode( root, "Obstacles" );
    if ( obstacles )
    {
      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 <<
      	"Reading new obstacles" << endl;
      bnewobst = true;

      DOMNodeList* allCompObstacles = ReaderXML::getNodes( obstacles );
      for (XMLSize_t i=0; i<allCompObstacles->getLength(); i++)
      {
        DOMNode* nCompObs = allCompObstacles->item( i );
        Obstacle *obstacle = ObstacleBuilderFactory::create( nCompObs );
        m_allcomponents.AddObstacle( obstacle );
      }

      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 <<
      	"Reading new obstacles completed" << endl;
    }


    // Contact force models
    DOMNode* contact = ReaderXML::getNode( root, "ContactForceModels" );
    if ( contact )
    {
      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 <<
      	"Reading new contact force models" << endl;
      ContactBuilderFactory::define( contact );
      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 <<
      	"Reading new contact force models completed" << endl;
    }
    string check_matA, check_matB;
    bool contactForceModels_ok =
    	ContactBuilderFactory::checkContactForceModelsExist( check_matA,
		check_matB );
    if ( !contactForceModels_ok )
    {
      if ( m_rank == 0 )
        cout << GrainsExec::m_shift6 << "No contact force model defined for "
		"materials : " << check_matA << " & " << check_matB << endl;
      grainsAbort();
    }


    // Check that construction is fine
    if ( !m_restart && !bnewpart && !bnewobst )
    {
      if ( m_rank == 0 )
        cout << "ERR : Error in input file in <Contruction>" << endl;
      grainsAbort();
    }


    // Set up the linked cell
    if ( m_rank == 0 ) cout << GrainsExec::m_shift6 <<
    	"Set up the linked cell grid" << endl;

    // Maximum circumscribed radius of particles
    double maxR = m_allcomponents.getCircumscribedRadiusMax();
    if ( maxR < 1.e-12 ) grainsAbort();
    else if ( m_rank == 0 )
      cout << GrainsExec::m_shift9 << "Maximum circumscribed particle radius = "
      	<< maxR << endl;

    // Scaling coefficient of linked cell size
    double LC_coef = 1.;
    DOMNode* nLC = ReaderXML::getNode( root, "LinkedCell" );
    if ( ReaderXML::hasNodeAttr( nLC, "CellSizeFactor" ) )
      LC_coef = ReaderXML::getNodeAttr_Double( nLC, "CellSizeFactor" );
    if ( LC_coef < 1. ) LC_coef = 1.;
    else if ( m_rank == 0 )
      cout << GrainsExec::m_shift9 << "Cell size factor = " << LC_coef << endl;

    // Define the linked cell grid
    defineLinkedCell( LC_coef * maxR, GrainsExec::m_shift9 );

    // If reload with 2024 format, read the particle reload file
    if ( reload && b2024 )
      m_allcomponents.read_particles( restart, npart, m_collision, m_rank, 
      	m_nprocs, m_wrapper );

    // Link obstacles with the linked cell grid
    m_collision->Link( m_allcomponents.getObstacles() );
  }
}




// ----------------------------------------------------------------------------
// External force definition
void GrainsTestDev::Forces( DOMElement* rootElement )
{
  if ( m_processorIsActive )
  {
    assert( rootElement != NULL );
    DOMNode* root = ReaderXML::getNode( rootElement, "Forces" );

    // Output message
    if ( m_rank == 0 ) cout << GrainsExec::m_shift3 << "Forces" << endl;


    // Read the forces
    if ( root )
    {
      // Gravity
      DOMNode* nGravity = ReaderXML::getNode( root, "Gravity" );
      if( nGravity )
      {
        GrainsExec::m_vgravity[X] = ReaderXML::getNodeAttr_Double(
      		nGravity, "GX" );
        GrainsExec::m_vgravity[Y] = ReaderXML::getNodeAttr_Double(
      		nGravity, "GY" );
        GrainsExec::m_vgravity[Z] = ReaderXML::getNodeAttr_Double(
      		nGravity, "GZ" );
        if ( m_rank == 0 ) cout << GrainsExec::m_shift6 << "Gravity = " <<
		GrainsExec::m_vgravity << endl;
      }
      else
      {
        if ( m_rank == 0 ) cout << GrainsExec::m_shift6 <<
		"Gravity is mandatory !!" << endl;
        grainsAbort();
      }
    }
    else
    {
      if ( m_rank == 0 )
      {
        cout << GrainsExec::m_shift6 << "No force specified"
      		<< endl;
        cout << GrainsExec::m_shift6 << "At least gravity is mandatory !!"
      		<< endl;
        grainsAbort();
      }
    }


    // Computes particle weight
    m_allcomponents.computeWeight( 0., 0. );
  }
}





// ----------------------------------------------------------------------------
// Additional features of the simulation: time features, insertion, 
// post-processing
void GrainsTestDev::AdditionalFeatures( DOMElement* rootElement )
{
  if ( m_processorIsActive )
  {
    assert( rootElement != NULL );
    DOMNode* root = ReaderXML::getNode( rootElement, "Simulation" );


    // Output message
    if ( m_rank == 0 ) cout << GrainsExec::m_shift3 << "Simulation" << endl;


    // Check that Simulation node exists
    if ( !root )
    {
      cout << GrainsExec::m_shift6 << "<Simulation> node is mandatory !!"
      		<< endl;
      grainsAbort();
    }


    // Time interval
    DOMNode* nTimeInterval = ReaderXML::getNode( root, "TimeInterval" );
    if ( nTimeInterval )
    {
      m_tstart = ReaderXML::getNodeAttr_Double( nTimeInterval, "Start" );
      m_tend = ReaderXML::getNodeAttr_Double( nTimeInterval, "End" );
      if ( GrainsExec::m_ReloadType == "same" ) m_tstart = m_time;
      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 << "Time interval = ["
      	<< m_tstart << "," << m_tend << "]" << endl;
    }
    else
    {
      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 <<
		"Time Interval is mandatory !!" << endl;
      grainsAbort();
    }

    // Time step
    DOMNode* nTimeStep = ReaderXML::getNode( root, "TimeStep" );
    if ( nTimeStep )
    {
      m_dt = ReaderXML::getNodeAttr_Double( nTimeStep, "dt" );
      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 <<
      	"Time step magnitude = " << m_dt << endl;
    }
    else
    {
      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 <<
		"Time step magnitude is mandatory !!" << endl;
      grainsAbort();
    }

    // Time integrator
    DOMNode* nTimeIntegration = ReaderXML::getNode( root, "TimeIntegration" );
    if ( nTimeIntegration )
      GrainsExec::m_TIScheme = ReaderXML::getNodeAttr_String( nTimeIntegration,
    		"Type" );
    if ( m_rank == 0 ) cout << GrainsExec::m_shift6 <<
      	"Time integration scheme = " << GrainsExec::m_TIScheme << endl;


    // Restart file and writing mode
    DOMNode* nRestartFile = ReaderXML::getNode( root, "RestartFile" );
    if ( nRestartFile )
    {
      m_fileSave = ReaderXML::getNodeAttr_String( nRestartFile, "Name" );
      if ( GrainsExec::m_ReloadType == "new" ) clearResultXmlFiles();
      GrainsExec::m_SaveDirectory = GrainsExec::extractRoot( m_fileSave );
      string wmode = ReaderXML::getNodeAttr_String( nRestartFile,
      	"WritingMode" );
      if ( wmode == "Hybrid" ) GrainsExec::m_writingModeHybrid = true ;
      if ( m_rank == 0 )
      {
        cout << GrainsExec::m_shift6 << "Restart file" << endl;
	cout << GrainsExec::m_shift9 << "File name = " << m_fileSave << endl;
        cout << GrainsExec::m_shift9 << "Directory = " <<
		GrainsExec::m_SaveDirectory << endl;
        cout << GrainsExec::m_shift9 << "Writing mode = " <<
		( GrainsExec::m_writingModeHybrid ? "Hybrid" : "Text" ) << endl;
      }
    }
    else
    {
      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 <<
		"RestartFile features are mandatory !!" << endl;
      grainsAbort();
    }


    // Output data frequency
    DOMNode* nTimeSave = ReaderXML::getNode( root, "TimeSave" );
    if ( nTimeSave )
    {
      double startSave = ReaderXML::getNodeAttr_Double( nTimeSave, "Start" );
      double endSave = ReaderXML::getNodeAttr_Double( nTimeSave, "End" );
      double dtSave = ReaderXML::getNodeAttr_Double( nTimeSave, "Every" );
      if ( dtSave < m_dt ) dtSave = m_dt;
      for (double t=startSave; t-endSave < 0.01 * m_dt; t+=dtSave)
        m_save.push_back(t);
      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 << "Output data every "
      	<< dtSave << " from " << startSave << " to " << endSave << endl;
    }
    else
    {
      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 <<
		"Output data time features are mandatory !!" << endl;
      grainsAbort();
    }


    // Moving obstacles
    DOMNode* nMovingObstacles = ReaderXML::getNode( root,
    	"MovingObstacles" );
    int ObstacleUpdateFreq = 1;
    bool displaceObstacles = true;
    if ( nMovingObstacles )
    {
      // Linked cell grid update frequency
      if ( ReaderXML::hasNodeAttr( nMovingObstacles, "LinkUpdateEvery" ) )
        ObstacleUpdateFreq = ReaderXML::getNodeAttr_Int( nMovingObstacles,
      		"LinkUpdateEvery" );

      // Whether moving obstacles are geometrically displaced
      if ( ReaderXML::hasNodeAttr( nMovingObstacles, "GeometricallyDisplace" ) )
      {
        string disp = ReaderXML::getNodeAttr_String( nMovingObstacles,
     		"GeometricallyDisplace" );
        if ( disp == "False" )
	{
	  displaceObstacles = false;
	  Obstacle::setMoveObstacle( false );
	}
      }
    }

    if ( m_rank == 0 )
    {
      cout << GrainsExec::m_shift6 << "Moving obstacles (if any)" << endl;
      cout << GrainsExec::m_shift9 <<
	"Moving obstacle - linked cell grid update every " <<
	ObstacleUpdateFreq << " time step" <<
	( ObstacleUpdateFreq > 1 ? "s" : "" ) << endl;
      cout << GrainsExec::m_shift9 <<
      	"Displace moving obstacles geometrically = " <<
      	( displaceObstacles ? "True" : "False" ) << endl;
    }


    // Particle insertion
    DOMNode* nInsertion = ReaderXML::getNode( root, "ParticleInsertion" );
    if ( nInsertion )
    {
      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 << "Particle insertion"
      	<< endl;

      // Insertion mode
      DOMNode* nMode = ReaderXML::getNode( nInsertion, "Mode" );
      if ( nMode )
      {
        string type = ReaderXML::getNodeAttr_String( nMode, "Type" );
	if ( type == "InitialTime" ) m_insertion_mode = IM_INITIALTIME;
	else if ( type == "OverTime" ) m_insertion_mode = IM_OVERTIME;
      }
      if ( m_rank == 0 ) cout << GrainsExec::m_shift9 << "Mode = " <<
      	( m_insertion_mode == IM_INITIALTIME ? "At initial time" :
		"Over time" ) << endl;


      // Insertion order
      DOMNode* nOrder = ReaderXML::getNode( nInsertion, "Order" );
      if ( nOrder )
      {
        string type = ReaderXML::getNodeAttr_String( nOrder, "Type" );
	if ( type == "Random" ) m_insertion_order = PM_RANDOM;
      }
      if ( m_rank == 0 ) cout << GrainsExec::m_shift9 << "Order = " <<
      	( m_insertion_order == PM_ORDERED ? "Ordered as defined in input file" :
		"Random" ) << endl;


      // Initial angular position
      DOMNode* nInitAngPos = ReaderXML::getNode( nInsertion,
      	"InitialAngularPosition" );
      if ( nInitAngPos )
      {
        string type = ReaderXML::getNodeAttr_String( nInitAngPos, "Type" );
	if ( type == "Random" ) m_init_angpos = IAP_RANDOM;
      }
      if ( m_rank == 0 ) cout << GrainsExec::m_shift9 << "Initial angular"
      	" position = " << ( m_init_angpos == IAP_FIXED ? "Fixed as defined in "
		"input file particle class" : "Random" ) << endl;


      // Random generator seed
      DOMNode* nRGS = ReaderXML::getNode( nInsertion, "RandomGeneratorSeed" );
      if ( nRGS )
      {
        string type = ReaderXML::getNodeAttr_String( nRGS, "Type" );
	if ( type == "Random" )
	{
	  m_randomseed = RGS_RANDOM;
	  srand( (unsigned int)( time(NULL)) );
	}
      }
      if ( m_rank == 0 ) cout << GrainsExec::m_shift9 << "Random generator"
      	" seed = " << ( m_randomseed == RGS_DEFAULT ? "Default to 1 "
	"(infinitely reproducible)" : "Initialized with running day/time "
	"(non-reproducible)" ) << endl;


      // Insertion attempt frequency
      DOMNode* nFrequency = ReaderXML::getNode( nInsertion,
      	"Frequency" );
      if ( nFrequency )
        m_insertion_frequency = size_t(
		ReaderXML::getNodeAttr_Int( nFrequency, "TryEvery" ));
      if ( m_rank == 0 ) cout << GrainsExec::m_shift9 << "Insertion "
      	"frequency = " << m_insertion_frequency << endl;


      // Force insertion
      DOMNode* nForceInsertion = ReaderXML::getNode( nInsertion,
      	"ForcedInsertion" );
      if ( nForceInsertion )
      {
        string value = ReaderXML::getNodeAttr_String( nForceInsertion,
		"Value" );
	if ( value == "True" ) m_force_insertion = true;
      }
      if ( m_rank == 0 ) cout << GrainsExec::m_shift9 << "Force insertion = " <<
      	( m_force_insertion  ? "True" : "False" ) << endl;


      // Particle positions via a external file OR a structured array OR a
      // collection of insertion windows, in this order of priority
      // Remark: these 3 modes cannot be combined
      DOMNode* nPosition = ReaderXML::getNode( nInsertion, "ParticlePosition" );
      if ( nPosition )
      {
        if ( m_rank == 0 ) cout << GrainsExec::m_shift6 << "Particle positions"
      		<< endl;

        // Fixed particle positions via an external file
        DOMNode* nFile = ReaderXML::getNode( nPosition, "File" );
        if ( nFile )
        {
		
          m_position = ReaderXML::getNodeAttr_String( nFile, "Name" );
          if ( m_rank == 0 ) cout << GrainsExec::m_shift9 << "External file = "
		<< m_position << endl;
        }
        else
        {
          // Fixed particle positions via a structured array
	  DOMNode* nStruct = ReaderXML::getNode( nPosition, "StructuredArray" );
	  if ( nStruct )
	  {
	    m_InsertionArray = new struct StructArrayInsertion;
            m_InsertionArray->box.ftype = WINDOW_BOX;
            m_InsertionArray->box.radius = m_InsertionArray->box.height = 0. ;
            m_InsertionArray->box.axisdir = NONE ;

            DOMNode* nBox = ReaderXML::getNode( nStruct, "Box" );
            if ( nBox )
            {
              DOMNodeList* points = ReaderXML::getNodes( nBox );
              DOMNode* pointA = points->item( 0 );
              DOMNode* pointB = points->item( 1 );
              m_InsertionArray->box.ptA[X] =
            	ReaderXML::getNodeAttr_Double( pointA, "X" );
              m_InsertionArray->box.ptA[Y] =
            	ReaderXML::getNodeAttr_Double( pointA, "Y" );
              m_InsertionArray->box.ptA[Z] =
            	ReaderXML::getNodeAttr_Double( pointA, "Z" );
              m_InsertionArray->box.ptB[X] =
            	ReaderXML::getNodeAttr_Double( pointB, "X" );
              m_InsertionArray->box.ptB[Y] =
            	ReaderXML::getNodeAttr_Double( pointB, "Y" );
              m_InsertionArray->box.ptB[Z] =
            	ReaderXML::getNodeAttr_Double( pointB, "Z" );
            }
	    else
	    {
              if ( m_rank == 0 ) cout << GrainsExec::m_shift6 <<
		"Node Box is required in <StructuredArray> !!" << endl;
              grainsAbort();
	    }

            DOMNode* nNumber = ReaderXML::getNode( nStruct, "Number" );
	    if ( nNumber )
	    {
              m_InsertionArray->NX = ReaderXML::getNodeAttr_Int( nNumber,
	      	"NX" );
              m_InsertionArray->NY = ReaderXML::getNodeAttr_Int( nNumber,
	      	"NY" );
              m_InsertionArray->NZ = ReaderXML::getNodeAttr_Int( nNumber,
	      	"NZ" );
	    }
	    else
	    {
              if ( m_rank == 0 ) cout << GrainsExec::m_shift6 <<
		"Node Number is required in <StructuredArray> !!" << endl;
              grainsAbort();
	    }

            m_position == "STRUCTURED";
	    if ( m_rank == 0 )
	    {
	      cout << GrainsExec::m_shift9 << "Structured array" << endl;
              cout << GrainsExec::m_shift12 << "Point3 min = " <<
	    	m_InsertionArray->box.ptA[X] << " " <<
		m_InsertionArray->box.ptA[Y] << " " <<
                m_InsertionArray->box.ptA[Z] << endl;
              cout << GrainsExec::m_shift12 << "Point3 max = " <<
	    	m_InsertionArray->box.ptB[X] << " " <<
		m_InsertionArray->box.ptB[Y] << " " <<
                m_InsertionArray->box.ptB[Z] << endl;
              cout << GrainsExec::m_shift12 << "Array = " <<
	    	m_InsertionArray->NX << " x " <<
		m_InsertionArray->NY << " x " <<
            	m_InsertionArray->NZ << endl;
	    }
	  }
	  else
	  {
	    // Random particle positions from a collection of insertion windows
	    DOMNode* nWindows = ReaderXML::getNode( nPosition, "Windows" );
            if ( nWindows )
	    {
	      cout << GrainsExec::m_shift9 << "Insertion windows" << endl;
	      DOMNodeList* allWindows = ReaderXML::getNodes( nWindows );
              for (XMLSize_t i=0; i<allWindows->getLength(); i++)
	      {
	        DOMNode* nWindow = allWindows->item( i );
                Window iwindow;
		readWindow( nWindow, iwindow, GrainsExec::m_shift12 );
	        m_insertion_windows.insert( m_insertion_windows.begin(),
	      		iwindow );
              }
	    }
            else
	    {
              if ( m_insertion_mode != IM_NOINSERT )
              {
                if ( m_rank == 0 ) cout << GrainsExec::m_shift6 <<
			"Insertion positions or windows are mandatory !!"
			<< endl;
                grainsAbort();
              }
	    }
	  }
        }
      }
      else
      {
        if ( m_insertion_mode != IM_NOINSERT )
        {
          if ( m_rank == 0 ) cout << GrainsExec::m_shift6 <<
		"Insertion positions/windows are mandatory !!" << endl;
          grainsAbort();
        }
      }

      // Initialization of particle velocity
      if ( m_rank == 0 )
        cout << GrainsExec::m_shift6 << "Particle initial velocity: ";

      DOMNode* nInitVit = ReaderXML::getNode( nInsertion, "InitialVelocity" );
      if ( nInitVit )
      {
        string sInitVitmode =
	    ReaderXML::getNodeAttr_String( nInitVit, "Mode" );

        if ( sInitVitmode == "Constant" )
        {
	  m_initvit_mode = IV_CONSTANT;

          DOMNode* nVitTransInit = ReaderXML::getNode( nInitVit,
	  	"TranslationalVelocity" );
          if ( nVitTransInit )
          {
            m_InitVtrans[X] = ReaderXML::getNodeAttr_Double( nVitTransInit,
	      	"VX" );
            m_InitVtrans[Y] = ReaderXML::getNodeAttr_Double( nVitTransInit,
	      	"VY" );
            m_InitVtrans[Z] = ReaderXML::getNodeAttr_Double( nVitTransInit,
	      	"VZ" );
          }

          DOMNode* nVitRotInit = ReaderXML::getNode( nInitVit,
	    	"AngularVelocity" );
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
        else if ( sInitVitmode == "Random" )
        {
	  m_initvit_mode = IV_RANDOM;

	  DOMNode* nRandomTrans = ReaderXML::getNode( nInitVit,
		"Translational" );
	  if ( nRandomTrans )
	    m_RandomMotionCoefTrans =
		ReaderXML::getNodeAttr_Double( nRandomTrans, "Amplitude" );

	  DOMNode* nRandomRot = ReaderXML::getNode( nInitVit,
	    	"Angular" );
	  if ( nRandomRot )
	    m_RandomMotionCoefRot =
	 	ReaderXML::getNodeAttr_Double( nRandomRot, "Amplitude" ) ;
        }
        else m_initvit_mode = IV_ZERO;
      }

      if ( m_rank == 0 )
      {
        switch( m_initvit_mode )
        {
          case IV_CONSTANT :
            cout << "constant" << endl;
	    cout << GrainsExec::m_shift9 << "translational = ( " <<
	   	 m_InitVtrans[X] << ", " << m_InitVtrans[Y] << ", "
		 << m_InitVtrans[Z] << " )" << endl;
	    cout << GrainsExec::m_shift9 << "angular = ( " <<
	   	 m_InitVrot[X] << ", " << m_InitVrot[Y] << ", "
		 << m_InitVrot[Z] << " )" << endl;
            break;

          case IV_RANDOM :
            cout << "random" << endl;
	    cout << GrainsExec::m_shift9 << "translational amplitude = " <<
	   	m_RandomMotionCoefTrans << endl;
	    cout << GrainsExec::m_shift9 << "angular amplitude = " <<
	   	m_RandomMotionCoefRot << endl;
            break;

          default :
            cout << "uniformly zero" << endl;
            break;
        }
      }
    }
    else
    {
      m_insertion_mode = IM_NOINSERT;
      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 << "No insertion" << endl;
    }


    // Obstacle loadings
    DOMNode* nObstacleLoadings = ReaderXML::getNode( root, "ObstacleLoadings" );
    size_t error = 0;
    if ( nObstacleLoadings )
    {
      if ( m_rank == 0 )
        cout << GrainsExec::m_shift6 << "Obstacle loadings" << endl;

      DOMNodeList* allOLs = ReaderXML::getNodes( nObstacleLoadings );
      for (XMLSize_t i=0; i<allOLs->getLength(); i++)
      {
        DOMNode* nOL = allOLs->item( i );
	string type = ReaderXML::getNodeAttr_String( nOL, "Type" );
	if ( m_rank == 0 )
          cout << GrainsExec::m_shift9 << "Type = " << type << endl;

	// Chargements en Force
	if ( type == "Force" )
	{
	  ObstacleImposedForce* load = new ObstacleImposedForce(
	      nOL, m_dt, m_rank, error );
	  if ( error != 0 ) grainsAbort();
	  else m_allcomponents.LinkImposedMotion( load );
	}
	else if ( type == "Velocity" )
	{
	  ObstacleImposedVelocity* load = new ObstacleImposedVelocity(
	  	nOL, m_dt, m_rank, error );
	  if ( error != 0 ) grainsAbort();
	  else m_allcomponents.LinkImposedMotion( load );
	}
	else
        {
	  if ( m_rank == 0 ) cout << GrainsExec::m_shift6 <<
		"Unknown obstacle loading type; values: Force or Velocity"
		<< endl;
          grainsAbort();
	}
      }
    }


    // Post-processing writers
    DOMNode* nPostProcessing = ReaderXML::getNode( root,
    	"PostProcessing" );
    if ( nPostProcessing )
    {
      if ( m_rank == 0 )
        cout << GrainsExec::m_shift6 << "Postprocessing" << endl;

      // Post-processing subdomain
      PostProcessingWriter::allocate_PostProcessingWindow( m_nprocs );

      DOMNode* nPostProcessingDomain = ReaderXML::getNode( nPostProcessing,
      	"Domain" );
      if ( nPostProcessingDomain )
      {
        DOMNodeList* nWindowPoints = ReaderXML::getNodes(
		nPostProcessingDomain );

 	DOMNode* pointA = nWindowPoints->item( 0 );
 	DOMNode* pointB = nWindowPoints->item( 1 );

 	Window PPWindow;
	PPWindow.ftype = WINDOW_BOX;
	PPWindow.radius = PPWindow.radius_int = PPWindow.height = 0. ;
	PPWindow.axisdir = NONE ;
 	PPWindow.ptA[X] = ReaderXML::getNodeAttr_Double( pointA, "X" );
	PPWindow.ptA[Y] = ReaderXML::getNodeAttr_Double( pointA, "Y" );
	PPWindow.ptA[Z] = ReaderXML::getNodeAttr_Double( pointA, "Z" );
	PPWindow.ptB[X] = ReaderXML::getNodeAttr_Double( pointB, "X" );
	PPWindow.ptB[Y] = ReaderXML::getNodeAttr_Double( pointB, "Y" );
	PPWindow.ptB[Z] = ReaderXML::getNodeAttr_Double( pointB, "Z" );

	double Ox, Oy, Oz, lx, ly, lz;
	bool b_X = false, b_Y = false, b_Z = false, b_PPWindow = false;
	App::get_local_domain_origin( Ox, Oy, Oz );
	App::get_local_domain_size( lx, ly, lz );

	if ( ( PPWindow.ptA[X] >= Ox && PPWindow.ptA[X] < Ox + lx )
		|| ( PPWindow.ptB[X] >= Ox && PPWindow.ptB[X] < Ox + lx )
		|| ( Ox > PPWindow.ptA[X] && Ox < PPWindow.ptB[X] )
		|| ( Ox > PPWindow.ptB[X] && Ox < PPWindow.ptA[X] ) )
	  b_X = true;
	if ( ( PPWindow.ptA[Y] >= Oy && PPWindow.ptA[Y] < Oy + ly )
		|| ( PPWindow.ptB[Y] >= Oy && PPWindow.ptB[Y] < Oy + ly )
		|| ( Oy > PPWindow.ptA[Y] && Oy < PPWindow.ptB[Y] )
		|| ( Oy > PPWindow.ptB[Y] && Oy < PPWindow.ptA[Y] ) )
	  b_Y = true;
	if ( ( PPWindow.ptA[Z] >= Oz && PPWindow.ptA[Z] < Oz + lz )
		|| ( PPWindow.ptB[Z] >= Oz && PPWindow.ptB[Z] < Oz + lz )
		|| ( Oz > PPWindow.ptA[Z] && Oz < PPWindow.ptB[Z] )
		|| ( Oz > PPWindow.ptB[Z] && Oz < PPWindow.ptA[Z] ) )
	  b_Z = true;

	if ( b_X && b_Y && b_Z ) b_PPWindow = true;

	PostProcessingWriter::set_PostProcessingWindow( m_rank, b_PPWindow );

	synchronize_PPWindow();

	if ( m_rank == 0 )
	{
	  cout << GrainsExec::m_shift9 << "Domain" << endl;
          cout << GrainsExec::m_shift12 << "Point3 A = " <<
		PPWindow.ptA[X] << " " << PPWindow.ptA[Y] << " " <<
		PPWindow.ptA[Z] << endl;
          cout << GrainsExec::m_shift12 << "Point3 B = " <<
		PPWindow.ptB[X] << " " << PPWindow.ptB[Y] << " " <<
		PPWindow.ptB[Z] << endl;
	}
      }
      else
	if ( m_rank == 0 )
	  cout << GrainsExec::m_shift9 << "Domain = linked cell grid"
	  	<< endl;


      // Post-processing writers
      DOMNode* nWriters = ReaderXML::getNode( nPostProcessing, "Writers" );
      if ( nWriters )
      {
        DOMNodeList* allPPW = ReaderXML::getNodes( nWriters );
        for (XMLSize_t i=0; i<allPPW->getLength(); i++)
        {
          DOMNode* nPPW = allPPW->item( i );
          PostProcessingWriter* ppw =
	  	PostProcessingWriterBuilderFactory::create(
      		nPPW, m_rank, m_nprocs );
          if ( ppw ) m_allcomponents.addPostProcessingWriter( ppw );
	  else
          {
	    if ( m_rank == 0 ) cout << GrainsExec::m_shift6 <<
		"Unknown postprocessing writer in node <Writers>"
		<< endl;
            grainsAbort();
	  }
        }
      }
      else
        if ( m_rank == 0 ) cout << GrainsExec::m_shift6
      		<< "No postprocessing writers" << endl;


      // Total Force & torque on obstacles
      DOMNode* nForceTorqueObstacles = ReaderXML::getNode( nPostProcessing,
      	"ForceTorqueObstacles" );
      if ( nForceTorqueObstacles )
      {
        int FToutputFreq = ReaderXML::getNodeAttr_Int( nForceTorqueObstacles,
		"Every" );
        string ppObsdir = ReaderXML::getNodeAttr_String( nForceTorqueObstacles,
		"Directory" );
        list<string> allppObsName;
        DOMNodeList* allppObs = ReaderXML::getNodes( nForceTorqueObstacles );
        for (XMLSize_t i=0; i<allppObs->getLength(); i++)
        {
          DOMNode* nppObs = allppObs->item( i );
          allppObsName.push_back(
		ReaderXML::getNodeAttr_String( nppObs, "ObstacleName" ) );
        }
        m_allcomponents.setOutputObstaclesLoadParameters( ppObsdir,
        	FToutputFreq, allppObsName );

	if ( m_rank == 0 )
	{
	  cout << GrainsExec::m_shift9 << "Force & torque on obstacles" << endl;
          cout << GrainsExec::m_shift12 << "Write values in file every " <<
		FToutputFreq << " time step" <<
		( FToutputFreq > 1 ? "s" : "" ) << endl;
          cout << GrainsExec::m_shift12 << "Output file directory name = "
    		<< ppObsdir << endl;
          cout << GrainsExec::m_shift12 << "Obstacle names" << endl;
	  for (list<string>::const_iterator il=allppObsName.begin();
	  	il!=allppObsName.end();il++)
	    cout << GrainsExec::m_shift15 << *il << endl;
	}
      }
    }
    else
      if ( m_rank == 0 ) cout << GrainsExec::m_shift6
      	<< "No postprocessing" << endl;
  }
}



// ----------------------------------------------------------------------------
// Runs the simulation over the prescribed time interval
void GrainsTestDev::Simulation( double time_interval )
{
// 
//   int rankproc = 0, nprocs = 0;
//   MPI_Comm_rank( MPI_COMM_WORLD, &rankproc );
//   MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
// 
//   // Input file
//   string filename = "Grains/Init/insert.xml";
//   size_t error = 0;
//   size_t pos = filename.find(".xml");
//   if ( pos == string::npos )
//   {
//     cout << "ERROR : input file need the .xml extension" << endl;
//     error = 1;
//   }
// 
//   // Creating STL file and reading it
//  // string str1 = "wall_double_ramp.stl";
//   //string str1 = "wall_unit.stl";
//   string str1 = "long_wall.stl";
//     
//   Obstacle *stlob = new STLObstacle( "testSTL", str1 );  
// 
//   cout << "Creating STLObstacle..." << endl;
// 
// 
//   // ******************************
//   //           TESTING 
//   // ******************************
// 
//   int test1 = 0, test2 = 1;
// 
//   if (test1)
//   {
// 
//   double xp = 0.5,   yp = 0.5,   zp = 0.5;
//   double R = 0.1;
// 
//   // Test 1: projection belongs to triangle
//   
//   Point3 Pa(0., 0.5, 0.);
//   Point3 Pb(-1000., 0.5, -1000.);
//   Point3 Pc(-1., 0., 1.);
//   Point3 Pd(-1., 0.1, 1.);
//   Point3 P1(-1., 0., 1.);
//   Point3 P2(1., 0., 1.);
//   Point3 P3(0., 0., -1.);
//  
//   
//   cout << "IsInter(Pa): " << STLObstacle::intersect(Pa, P1, P2, P3) << endl;
//   cout << "IsInter(Pb): " << STLObstacle::intersect(Pb, P1, P2, P3) << endl;
//   cout << "IsInter(Pc): " << STLObstacle::intersect(Pc, P1, P2, P3) << endl;
//   cout << "IsInter(Pd): " << STLObstacle::intersect(Pd, P1, P2, P3) << endl;
// 
//   }
// 
//   cout << "Testing again..." << endl;
// 
//   // Test 2: triangle area
//  
//   if (test2)
//   {
// 
//   double x1, y1, z1, x2, y2, z2, x3, y3, z3;
// 
//   x1 = 0.0;
//   y1 = 0.0;
//   z1 = 0.0;
// 
//   x2 = 1.0;
//   y2 = 0.0;
//   z2 = 0.0;
// 
//   x3 = 0.0;
//   y3 = 0.0;
//   z3 = 1.0;
// 
//   Vector3 n;
//   n[0] = 0.; n[1] = 0.; n[2] = 0.;
// 
//   size_t vid = 0;
//   STLVertex *v1 = new STLVertex(x1,y1,z1,n,vid); 
//   STLVertex *v2 = new STLVertex(x2,y2,z2,n,vid);  
//   STLVertex *v3 = new STLVertex(x3,y3,z3,n,vid);
// 
//   size_t tid = 0;
//   tuple<STLVertex*,STLVertex*,STLVertex*> vetest;
//   vetest = std::make_tuple(v1, v2, v3);
//   STLTriangle trtest(vetest,n,tid); 
// 
//   cout << "Testing area computation: " << trtest.getSurfaceArea() << endl;
// 
//   }
// 
//   // Test 3: RBWC
// //   
// //   Rectangle *Rect = new Rectangle( 1.0, 1.0 );
// //   Rect->ClosestPoint( RigidBodyWithCrust &Rect );

}
