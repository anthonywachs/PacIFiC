#include "Grains.hh"
#include "ContactBuilderFactory.hh"
#include "LinkedCell.hh"
#include "AppBrownian.hh"
#include "ObstacleBuilderFactory.hh"
#include "ObstacleImposedVelocity.hh"
#include "GrainsBuilderFactory.hh"
#include "PostProcessingWriter.hh"
#include "PostProcessingWriterBuilderFactory.hh"
#include "RawDataPostProcessingWriter.hh"
#include "SpheroCylinder.hh"
#include "Quadrilobe.hh"
#include "stdlib.h"


// Initialisation des attributs static
bool Grains::m_predictor_mode = true;


// ----------------------------------------------------------------------------
// Default constructor
Grains::Grains()
  : ComputingTime("Solver")
  , m_collision( NULL )
  , m_tstart( 0. )
  , m_tend( 0. )
  , m_dt( 0. )
  , m_time( 0. )
  , m_lastTime_save( false )
  , m_error_occured( false )
  , m_fileSave( "undefined" )
  , m_clonesInReloadFile( CIR_NONE )
  , m_dimension( 3 )
  , m_allProcTiming( true )
  , m_restart( false )
  , m_insertion_order( PM_ORDERED )
  , m_insertion_mode( IM_NOINSERT )
  , m_initvit_mode( IV_ZERO )
  , m_init_angpos( IAP_FIXED )
  , m_randomseed( RGS_DEFAULT )
  , m_InsertionArray( NULL )
  , m_insertion_position( NULL )
  , m_insertion_angular_position( NULL )  
  , m_insertion_frequency( 1 )
  , m_force_insertion( false )
  , m_RandomMotionCoefTrans( 0. )
  , m_RandomMotionCoefRot( 0. )
  , m_npwait_nm1( 0 )
  , m_rank( 0 )
  , m_nprocs( 1 )
  , m_processorIsActive( true )
  , m_wrapper( NULL )
{
  if ( GrainsBuilderFactory::getContext() == DIM_2 ) m_dimension = 2;
  m_insertion_position = new list<Point3>;
  m_insertion_angular_position = new list<Matrix>;
}




// ----------------------------------------------------------------------------
// Destructor
Grains::~Grains()
{
  if ( m_InsertionArray ) delete m_InsertionArray;
  list<App*>::iterator app;
  for (app=m_allApp.begin(); app!=m_allApp.end(); app++) delete *app;
  m_newParticles.clear();
  ContactBuilderFactory::eraseAllContactForceModels();
  GrainsExec::GarbageCollector();
  m_insertion_position->clear();
  delete m_insertion_position;
  m_insertion_angular_position->clear();
  delete m_insertion_angular_position;  
}




// ----------------------------------------------------------------------------
// Writes an initial message in the standard output only on the process ranked 0
void Grains::initialOutputMessage()
{
  cout << "Grains3D serial" << endl;
}



// ----------------------------------------------------------------------------
// Tasks to perform before time-stepping
void Grains::do_before_time_stepping( DOMElement* rootElement )
{
  // Read the input file
  Construction( rootElement );
  Forces( rootElement );
  AdditionalFeatures( rootElement );

  if ( m_processorIsActive )
  {
    // Initialisation
    double vmax = 0., vmean = 0. ;

    // Timers
    CT_set_start();
    SCT_insert_app( "Initialization" );
    SCT_set_start( "Initialization" );
    if ( m_rank == 0 ) cout << endl << "Initialization" << endl;

    // Set time to initial time
    m_time = m_tstart;

    // Particle creation, insertion and link to grid
    InsertCreateNewParticles();

    // Number of particles: inserted and in the system
    m_allcomponents.computeNumberParticles( m_wrapper );
    m_npwait_nm1 = m_allcomponents.getNumberPhysicalParticlesToInsert();

    // Initialisation of obstacle kinematics
    m_allcomponents.setKinematicsObstacleWithoutMoving( m_time, m_dt );

    // In case of initial random motion
    if ( m_initvit_mode == IV_RANDOM )
      m_allcomponents.setRandomMotion( m_RandomMotionCoefTrans,
	m_RandomMotionCoefRot );

    // Writing results for postprocessing
    m_allcomponents.PostProcessing_start( m_time, m_dt, m_collision,
  	m_insertion_windows, m_rank, m_nprocs, m_wrapper );

    // Track component max and mean velocity
    m_allcomponents.ComputeMaxMeanVelocity( vmax, vmean, m_wrapper );
    if ( m_rank == 0 ) 
    {
      fVitMax.open( (m_fileSave + "_VelocityMaxMean.dat").c_str(), ios::out );
      cout << "Component velocity : max = " << vmax << " average = " <<
    	vmean << endl;
      fVitMax << GrainsExec::doubleToString( ios::scientific, 6, m_time )
    	<< "\t" << GrainsExec::doubleToString( ios::scientific, 6, vmax )
	<< "\t" << GrainsExec::doubleToString( ios::scientific, 6, vmean )
	<< endl;
    }

    // Display memory used by Grains
    display_used_memory();

    // Postprocessing of force & torque on obstacles
    m_allcomponents.initialiseOutputObstaclesLoadFiles( m_rank, false, m_time );
    m_allcomponents.outputObstaclesLoad( m_time, m_dt, false,
      GrainsExec::m_ReloadType == "same", m_rank );
      
    // Postprocessing of force statistics
    m_collision->initialiseForceStatsFiles( m_rank, false, m_time );      
    bool forcestats = m_collision->outputForceStatsAtThisTime( false, 
	GrainsExec::m_ReloadType == "same" );	
    if ( forcestats )
      m_collision->outputForceStats( m_time, m_dt, m_rank, m_wrapper );

    // Next time of writing results
    m_timeSave = m_save.begin();
    while ( *m_timeSave - m_time < 0.01 * m_dt && m_timeSave != m_save.end() )
      m_timeSave++;
      
    // In case of SecondOrderLeapFrog, initialize acceleration
    if ( GrainsExec::m_TIScheme == "SecondOrderLeapFrog" && 
    	!GrainsExec::m_isReloaded )
    {
      computeParticlesForceAndTorque();
      m_allcomponents.setAllContactMapCumulativeFeaturesToZero();
    }

    if ( m_rank == 0 ) cout << "Initialization completed" << endl << endl;
    SCT_get_elapsed_time( "Initialization" );
  }
}




// ----------------------------------------------------------------------------
// Tasks to perform before time-stepping
void Grains::do_after_time_stepping()
{
  if ( m_processorIsActive )
  {
    double vmax = 0., vmean = 0. ;

    // Particles in & out at the end of the simulation 
    if ( m_rank == 0 )
    {  
      cout << endl << "Number of active particles in the simulation = " 
    	<< m_allcomponents.getNumberActiveParticlesOnAllProc() << endl;
      cout << "Number of particles to insert in the simulation = " 
    	<< m_allcomponents.getNumberPhysicalParticlesToInsert() << endl;
    }	

    // Write reload files
    if ( !m_lastTime_save )
    {
      SCT_set_start( "OutputResults" );

      // Track component max and mean velocity at the end of the simulation
      m_allcomponents.ComputeMaxMeanVelocity( vmax, vmean, m_wrapper );
      if ( m_rank == 0 )
      {
        cout << endl << "Component velocity : max = " << vmax
		<< " average = " << vmean << endl;
        fVitMax << GrainsExec::doubleToString( ios::scientific, 6, m_time )
		<< "\t" << GrainsExec::doubleToString( ios::scientific, 6,
		vmax ) << "\t" << GrainsExec::doubleToString( ios::scientific,
		6, vmean ) << endl;
        fVitMax.close();
      }

      // Reload files are always written at the end of the simulation
      if ( !m_error_occured ) saveReload( m_time );

      SCT_get_elapsed_time( "OutputResults" );
    }

    // Final tasks performed by postprocessing writers
    m_allcomponents.PostProcessing_end();

    // Contact features over the simulation
    double omax = m_collision->getOverlapMax(),
  	omean = m_collision->getOverlapMean(),
	timemax = m_collision->getTimeOverlapMax(),
	ngjk = m_collision->getNbIterGJKMean();
    if ( m_wrapper ) m_wrapper->ContactsFeatures( omax, omean, timemax, ngjk );	

    if ( m_rank == 0 )
    {   
      cout << endl << "Contact features over the simulation" << endl;
      cout << GrainsExec::m_shift3 << "Minimal crust thickness = " <<
    	m_allcomponents.getCrustThicknessMin() << endl;
      cout << GrainsExec::m_shift3 << "Maximum overlap = " << omax << endl;
      cout << GrainsExec::m_shift3 << "Average overlap = " << omean << endl;
      cout << GrainsExec::m_shift3 << "Time of maximum overlap = " <<
	timemax << endl;
      cout << GrainsExec::m_shift3 << "Average number of iterations of GJK = " <<
    	ngjk << endl << endl;
    }

    // Timer outcome
    display_timer_summary();
  }
}




// ----------------------------------------------------------------------------
// Runs the simulation over the prescribed time interval
void Grains::Simulation( double time_interval )
{
  double vmax = 0., vmean = 0. ;
  bool forcestats = false;

  // Timers
  SCT_insert_app( "ParticlesInsertion" );
  SCT_insert_app( "ComputeForces" );
  SCT_insert_app( "Move" );
  SCT_insert_app( "UpdateParticleActivity" );
  SCT_insert_app( "LinkUpdate" );
  SCT_insert_app( "OutputResults" );

  // Simulation: time marching algorithm
  cout << "Time \t TO \tend \tParticles \tIn \tOut" << endl;
  while ( m_tend - m_time > 0.01 * m_dt )
  {
    try
    {
      m_time += m_dt;
      GrainsExec::m_time_counter++;

      // Check whether data are output at this time
      m_lastTime_save = false;
      GrainsExec::m_output_data_at_this_time = false;
      GrainsExec::m_postprocess_forces_at_this_time = false;
      if ( m_timeSave != m_save.end() )
        if ( *m_timeSave - m_time < 0.01 * m_dt )
	{
	  // Set the global data output boolean to true
	  GrainsExec::m_output_data_at_this_time = true;

	  // Next time of writing files
	  m_timeSave++;

	  m_lastTime_save = true;
	}
      forcestats = m_collision->outputForceStatsAtThisTime( false, false );
      if ( GrainsExec::m_output_data_at_this_time || forcestats )
      {
        GrainsExec::m_postprocess_forces_at_this_time = true;
	m_collision->resetPPForceIndex();
      }
	

      // Insertion of particles
      SCT_set_start( "ParticlesInsertion" );      
      m_npwait_nm1 = m_allcomponents.getNumberPhysicalParticlesToInsert();      
      if ( m_insertion_mode == IM_OVERTIME )
        insertParticle( m_insertion_order );      
      m_allcomponents.computeNumberParticles( m_wrapper );
      if ( m_npwait_nm1 
      	!= m_allcomponents.getNumberPhysicalParticlesToInsert() )
          cout << "\r                                              "
               << "                 " << flush;
      ostringstream oss;
      oss.width(10);
      oss << left << m_time;
      cout << '\r' << oss.str() << "  \t" << m_tend << "\t\t\t"
           << m_allcomponents.getNumberActiveParticlesOnProc() << '\t'
           << m_allcomponents.getNumberPhysicalParticlesToInsert() 
	   << flush;
      SCT_get_elapsed_time( "ParticlesInsertion" );


      if ( GrainsExec::m_TIScheme == "SecondOrderLeapFrog" )
      {
        // We implement the kick-drift-kick version of the LeapFrog scheme

	// Move particles and obstacles
	// Update particle velocity over dt/2 and particle position over dt,
	// obstacle velocity and position over dt
	// v_i+1/2 = v_i + a_i * dt / 2
	// x_i+1 = x_i + v_i+1/2 * dt
        moveParticlesAndObstacles( 0.5 * m_dt, m_dt, m_dt );
	
        // Compute particle forces and torque
	// Compute f_i+1 and a_i+1 as a function of (x_i+1,v_i+1/2)
	SCT_set_start( "ComputeForces" );
        computeParticlesForceAndTorque();
        SCT_get_elapsed_time( "ComputeForces" );	
	
	// Update particle velocity over dt/2
	// v_i+1 = v_i+1/2 + a_i+1 * dt / 2 
	m_allcomponents.advanceParticlesVelocity( m_time, 0.5 * m_dt );
      }
      else
      {
        // Compute particle forces and torque
	// Compute f_i and a_i as a function of (x_i,v_i) 
	SCT_set_start( "ComputeForces" );
        computeParticlesForceAndTorque();
        SCT_get_elapsed_time( "ComputeForces" );
      
        // Move particles and obstacles
	// x_i+1 = x_i + g(v_i,v_i-1,a_i,dt)
	// v_i+1 = v_i + g(a_i,a_i-1,dt)
        moveParticlesAndObstacles( m_dt, m_dt, m_dt );
      }


      // Compute and write force & torque exerted on obstacles
      m_allcomponents.computeObstaclesLoad( m_time, m_dt ); 
      m_allcomponents.outputObstaclesLoad( m_time, m_dt, false, false, m_rank );
      
      
      // Compute and write force statistics
      if ( forcestats ) m_collision->outputForceStats( m_time, m_dt, m_rank, 
      	m_wrapper );
        

      // Write postprocessing and reload files
      if ( GrainsExec::m_output_data_at_this_time )
      {
	SCT_set_start( "OutputResults" );

	// Track component max and mean velocity
	m_allcomponents.ComputeMaxMeanVelocity( vmax, vmean );
        cout << endl << "Component velocity : max = " << vmax
               << " average = " << vmean << endl;
        fVitMax << GrainsExec::doubleToString( ios::scientific, 6, m_time )
               << "\t" << GrainsExec::doubleToString( ios::scientific, 6,
               vmax ) << "\t" << GrainsExec::doubleToString( ios::scientific,
               6, vmean ) << endl;

	// Display memory used by Grains
	display_used_memory();

	// Write reload files
	saveReload( m_time );

	// Write postprocessing files
        m_allcomponents.PostProcessing( m_time, m_dt, m_collision );

	SCT_get_elapsed_time( "OutputResults" );
      }
    }
    catch ( ContactError& errContact )
    {
      // Max overlap exceeded
      cout << endl;
      m_allcomponents.PostProcessingErreurComponents( "ContactError",
            errContact.getComponents() );
      errContact.Message( cout );
      m_error_occured = true;
      break;
    }
    catch ( MotionError& errMotion )
    {
      // Particle motion over dt is too large
      cout << endl;
      m_allcomponents.PostProcessingErreurComponents( "MotionError",
            errMotion.getComponent() );
      errMotion.Message(cout);
      m_error_occured = true;
      break;
    }
    catch ( SimulationError& errSimulation )
    {
      // Simulation error
      cout << endl;
      errSimulation.Message(cout);
      m_error_occured = true;
      break;
    }
  }
}




// ----------------------------------------------------------------------------
// Computes particle forces and torques
void Grains::computeParticlesForceAndTorque()
{
  // Initiliaze all component transforms with crust to non computed
  m_allcomponents.InitializeRBTransformWithCrustState( m_time, m_dt );

  // Initialize contact maps (for contact model with memory)
  // TODO: add if-statement that bypasses this step when contact model has
  // no memory.
  m_allcomponents.setAllContactMapToFalse();

  // Compute forces from all applications
  // Initialisation torsors with weight only
  m_allcomponents.InitializeForces( m_time, m_dt, true );

  // Compute forces
  for (list<App*>::iterator app=m_allApp.begin(); app!=m_allApp.end(); app++)
    (*app)->ComputeForces( m_time, m_dt,
	m_allcomponents.getActiveParticles() );

  // Update contact maps (for contact model with memory)
  // TODO: add if-statement that bypasses this step when contact model has
  // no memory.
  m_allcomponents.updateAllContactMaps();
}




// ----------------------------------------------------------------------------
// Moves particles and obstacles	
void Grains::moveParticlesAndObstacles( double const& dt_particle_vel, 
    	double const& dt_particle_disp,
	double const& dt_obstacle )
{
  // Solve Newton's law and move particles
  SCT_set_start( "Move" );
  m_allcomponents.Move( m_time, dt_particle_vel, dt_particle_disp, 
      	dt_obstacle, m_collision );


  // In case of periodicity, update periodic clones and destroy periodic
  // clones that are out of the linked cell grid
  if ( m_periodic )
    m_collision->updateDestroyPeriodicClones(
		m_allcomponents.getActiveParticles(),
		m_allcomponents.getPeriodicCloneParticles() );
  SCT_get_elapsed_time( "Move" );


  // Update particle activity
  SCT_set_start( "UpdateParticleActivity" );
  m_allcomponents.UpdateParticleActivity();
  SCT_get_elapsed_time( "UpdateParticleActivity" );


  // Update the particle & obstacles links with the grid
  SCT_set_start( "LinkUpdate" );
  m_collision->LinkUpdate( m_time, m_dt, m_allcomponents.getActiveParticles() );


  // In case of periodicity, create new periodic clones and destroy periodic
  // clones because the master particle has changed tag/geoposition
  if ( m_periodic )
    m_collision->createDestroyPeriodicClones(
	m_allcomponents.getActiveParticles(),
	m_allcomponents.getPeriodicCloneParticles(),
	m_allcomponents.getReferenceParticles() );
  SCT_get_elapsed_time( "LinkUpdate" );
}




// ----------------------------------------------------------------------------
// Construction of the simulation: linked cell, particles & obstacles, domain
// decomposition
void Grains::Construction( DOMElement* rootElement )
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
      GrainsExec::m_isReloaded = m_restart;

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
      restart = GrainsExec::fullResultFileName( restart, false );
      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 <<
      	"Simulation reloaded from " << restart << endl;      

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
      assert( cle == "#EndTime" );
      simulLoad.close();      

      // Whether to reset velocity to 0
      string reset = ReaderXML::getNodeAttr_String( reload, "Velocity" );
      m_allcomponents.resetKinematics( reset );
    }


    // Collision detection
    DOMNode* collision = 
                      ReaderXML::getNode( root, "CollisionDetectionAlgorithm" );
    if ( !collision )
    {
      if ( m_rank == 0 )
        cout << GrainsExec::m_shift6 <<
              "Default collision detection algorithm using GJK, " <<
              "1E-15 tolerance, " << 
              "without acceleration, " <<
              "and without any bounding volume collision detection!" << endl;
    }
    else 
    {
      DOMNode* collisionAlg = 
                       ReaderXML::getNode( collision, "CollisionDetection" );
      if ( collisionAlg )
      {
        string nCollisionAlg = 
                          ReaderXML::getNodeAttr_String( collisionAlg, "Type" );
        if ( nCollisionAlg != "GJK" )
        {
          if ( m_rank == 0 )
	    cout << GrainsExec::m_shift6 <<
              "Collision detection algorithm is not defined!" << endl;
          grainsAbort();
        }

        double tol = ReaderXML::getNodeAttr_Double( collisionAlg, "Tolerance" );
        if ( tol < 1e-15 )
        {
          if ( m_rank == 0 )
	    cout << GrainsExec::m_shift6 <<
              "Tolerance should be greater than 1E-15!" << endl;
          grainsAbort();
        }
        else
          GrainsExec::m_colDetTolerance = tol;
        
        string acc = 
                  ReaderXML::getNodeAttr_String( collisionAlg, "Acceleration" );
        if ( acc == "ON" )
          GrainsExec::m_colDetAcceleration = true;
        else if ( acc == "OFF" )
          GrainsExec::m_colDetAcceleration = false;
        else
        {
          if ( m_rank == 0 )
	    cout << GrainsExec::m_shift6 <<
              "Acceleration should be ON or OFF!" << endl;
          grainsAbort();
        }

        if ( m_rank == 0 )
	  cout << GrainsExec::m_shift6 <<
              "Collision detection algorithm using " <<
              GrainsExec::m_colDetMethod << ", " <<
              GrainsExec::m_colDetTolerance <<
              " tolerance, and acceleration is " << 
              ( GrainsExec::m_colDetAcceleration ? "on." : "off." ) << endl;
      }
      else
      {
        if ( m_rank == 0 )
	  cout << GrainsExec::m_shift6 <<
              "Default collision detection algorithm using GJK, " <<
              "1E-15 tolerance, and without acceleration!" << endl;
      }

      DOMNode* bVolumeAlg = 
                       ReaderXML::getNode( collision, "BoundingVolume" );
      if ( bVolumeAlg )
      {
        string nBVtype = ReaderXML::getNodeAttr_String( bVolumeAlg, "Type" );
        if ( nBVtype == "OBB" )
          GrainsExec::m_colDetBoundingVolume = 1;
        else if ( nBVtype == "OBC" )
          GrainsExec::m_colDetBoundingVolume = 2;
        else
          GrainsExec::m_colDetBoundingVolume = 0;
      }
      if ( m_rank == 0 && GrainsExec::m_colDetBoundingVolume == 1 )
        cout << GrainsExec::m_shift6 <<
        "Pre-collision Test with oriented bounding boxes." << endl;
      else if ( m_rank == 0 && GrainsExec::m_colDetBoundingVolume == 2 )
        cout << GrainsExec::m_shift6 <<
        "Pre-collision Test with oriented bounding cylinders." << endl;
      else if ( m_rank == 0 )
        cout << GrainsExec::m_shift6 <<
        "Pre-collision Test with bounding volumes is off." << endl;
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
        int nb = ReaderXML::getNodeAttr_Int( nParticle, "Number" );
        
        // Remark: reference particles' ID number is -1, which explains
        // auto_numbering = false in the constructor
        Particle* particleRef = new Particle( nParticle, nbPC+int(i) );
        m_allcomponents.AddReferenceParticle( particleRef, nb  );
        pair<Particle*,size_t> ppp( particleRef, nb );
        m_newParticles.push_back( ppp );
      }

      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 <<
      	"Reading new particle types completed" << endl;
    }
    
    if ( GrainsExec::getMinCrustThickness() <= 0. )
    {
      if ( m_rank == 0 )
        cout << "Some Particles have been defined with a negative or"
		" 0 crust thickness" << endl;
      grainsAbort();            
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
	else if ( sshape == "Quadrilobe" )
	  particleRef = new Quadrilobe( nCompParticle, nbPC+int(i) );
	else 	
	  particleRef = new CompositeParticle( nCompParticle, nbPC+int(i) );
        m_allcomponents.AddReferenceParticle( particleRef, nb );
        pair<Particle*,size_t> ppp( particleRef, nb );
        m_newParticles.push_back( ppp );
      }

      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 <<
      	"Reading new composite particle types completed" << endl;
    }
    
    if ( GrainsExec::getMinCrustThickness() <= 0. )
    {
      if ( m_rank == 0 )
        cout << "Some Composite Particles have been defined with a negative or"
		" 0 crust thickness" << endl;
      grainsAbort();            
    }    
    
    if  ( m_rank == 0 ) cout << GrainsExec::m_shift6 <<
      	"Number of particle types = " << 
	m_allcomponents.getReferenceParticles()->size() << endl;


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
    
    
    // Check that crust thickness is positive for all rigid bodies
    if ( GrainsExec::getMinCrustThickness() <= 0. )
    {
      if ( m_rank == 0 )
        cout << "Some Obstacles have been defined with a negative or"
		" 0 crust thickness" << endl;
      grainsAbort();            
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
        cout << "No contact force model defined for "
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
    if ( maxR < 1.e-12 ) maxR = 1.e16;
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
    if ( m_restart && b2024 )
      m_allcomponents.read_particles( restart, npart, m_collision, m_rank, 
      	m_nprocs, m_wrapper );  

    // Link the root obstacle with the linked cell grid
    m_collision->Link( m_allcomponents.getObstacles() );
  }
}




// ----------------------------------------------------------------------------
// External force definition
void Grains::Forces( DOMElement* rootElement )
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
      if ( nGravity )
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
      
      // Brownian force
      DOMNode* nAppBrownian = ReaderXML::getNode( root, "Brownian" );
      if ( nAppBrownian )
      {
        size_t error = 0;
	App* force = new AppBrownian( nAppBrownian, m_rank, error );
	if ( error ) grainsAbort();
	else m_allApp.push_back( force );
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
void Grains::AdditionalFeatures( DOMElement* rootElement )
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
    if ( GrainsExec::m_TIScheme != "SecondOrderLeapFrog" )
      m_allcomponents.setTimeIntegrationScheme();
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
      if ( ReaderXML::hasNodeAttr( nRestartFile, "SingleFile" ) )
      {
        string wsf = ReaderXML::getNodeAttr_String( nRestartFile,
      		"SingleFile" );
        if ( wsf == "True" && m_nprocs > 1 ) 
	  GrainsExec::m_SaveMPIInASingleFile = true ;
      }
      if ( ReaderXML::hasNodeAttr( nRestartFile, "Clones" ) )
      {
        string wclones = ReaderXML::getNodeAttr_String( nRestartFile,
      		"Clones" );
        if ( wclones == "All" ) m_clonesInReloadFile = CIR_ALL;
	else if ( wclones == "NoPer" ) m_clonesInReloadFile = CIR_NOPERIODIC;
	else m_clonesInReloadFile = CIR_NONE;
      }      		        
      if ( m_rank == 0 )
      {
        cout << GrainsExec::m_shift6 << "Restart file" << endl;
	cout << GrainsExec::m_shift9 << "File name = " << m_fileSave << endl;
        cout << GrainsExec::m_shift9 << "Directory = " <<
		GrainsExec::m_SaveDirectory << endl;
        cout << GrainsExec::m_shift9 << "Writing mode = " <<
		( GrainsExec::m_writingModeHybrid ? "Hybrid" : "Text" ) << endl;
        cout << GrainsExec::m_shift9 << "Single file = " <<
		( GrainsExec::m_SaveMPIInASingleFile ? "True" : "False" ) 
		<< endl;
        cout << GrainsExec::m_shift9 << "Clones = " <<
		( m_clonesInReloadFile == CIR_ALL ? "All" : 
		m_clonesInReloadFile == CIR_NOPERIODIC ? "No periodic clones" :
		"None" ) << endl;			
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
    bool geomMoveObstacles = true;
    if ( nMovingObstacles )
    {
      // Linked cell grid update frequency
      if ( ReaderXML::hasNodeAttr( nMovingObstacles, "LinkUpdateEvery" ) )
        ObstacleUpdateFreq = ReaderXML::getNodeAttr_Int( nMovingObstacles,
      		"LinkUpdateEvery" );

      // Whether moving obstacles are geometrically moving
      if ( ReaderXML::hasNodeAttr( nMovingObstacles, "GeometricallyMove" ) )
      {
        string disp = ReaderXML::getNodeAttr_String( nMovingObstacles,
     		"GeometricallyMove" );
        if ( disp == "False" )
	{
	  geomMoveObstacles = false;
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
      	"Move obstacles geometrically = " <<
      	( geomMoveObstacles ? "True" : "False" ) << endl;
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
	else if ( type == "File" ) 
	{
	  m_init_angpos = IAP_FILE;
	  m_angular_position = 
	    ReaderXML::getNodeAttr_String( nInitAngPos, "Name" );
	}
      }
      if ( m_rank == 0 ) cout << GrainsExec::m_shift9 << "Initial angular"
      	" position = " << ( m_init_angpos == IAP_FIXED ? "Fixed as defined in "
		"input file particle class" : ( 
		m_init_angpos == IAP_RANDOM ? "Random" : "Fixed defined in "
		"file " + m_angular_position ) ) << endl;


      // Random generator seed
      DOMNode* nRGS = ReaderXML::getNode( nInsertion, "RandomGeneratorSeed" );
      int seed = 0;
      if ( nRGS )
      {
        string type = ReaderXML::getNodeAttr_String( nRGS, "Type" );
	if ( type == "UserDefined" )
	{
	  m_randomseed = RGS_UDEF;
	  seed = ReaderXML::getNodeAttr_Int( nRGS, "Value" );
	  srand( (unsigned int)( seed ) );
	}
	else if ( type == "Random" )
	{
	  m_randomseed = RGS_RANDOM;
	  srand( (unsigned int)( time(NULL) ) );
	}
      }
      if ( m_rank == 0 ) cout << GrainsExec::m_shift9 << "Random generator"
      	" seed = " << ( m_randomseed == RGS_DEFAULT ? "Default to 1 "
	"(infinitely reproducible)" : ( m_randomseed == RGS_UDEF ? 
	"Set to "+GrainsExec::intToString( seed )+" (infinitely reproducible)" :
	"Initialized with running day/time (non-reproducible)" ) ) << endl;


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
      	"ForceInsertion" );
      if ( nForceInsertion )
      {
        string value = ReaderXML::getNodeAttr_String( nForceInsertion,
		"Value" );
	if ( value == "True" ) m_force_insertion = true;
      }
      if ( m_rank == 0 ) cout << GrainsExec::m_shift9 << "Force insertion = " <<
      	( m_force_insertion  ? "True" : "False" ) << endl;


      // Particle positions via an external file OR a structured array OR a
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

            m_position = "STRUCTURED";
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
	      if ( m_rank == 0 )
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
                if ( m_rank == 0 )
                  cout << GrainsExec::m_shift6 <<
            "Insertion positions or windows are mandatory !!" << endl;
                grainsAbort();
              }
	    }
	  }
        }
      }
      else
      {
        if ( m_insertion_mode != IM_NOINSERT && !getNbInsertionPositions() )
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

	// Imposed force
	if ( type == "Force" )
	{
	  ObstacleImposedForce* load = new ObstacleImposedForce(
	      nOL, m_dt, m_rank, error );
	  if ( error != 0 ) grainsAbort();
	  else m_allcomponents.LinkImposedMotion( load );
	}
	// Imposed velocity	
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
          cout << GrainsExec::m_shift12 << "Obstacle name";
	  if ( allppObsName.size() == 1 ) cout << " = " << *allppObsName.begin()
	  	<< endl;
	  else
	  {
	    cout << "s" << endl;
	    for (list<string>::const_iterator il=allppObsName.begin();
	  	il!=allppObsName.end();il++)
	      cout << GrainsExec::m_shift15 << *il << endl;
	  }
	}
      }
      
      
      // Force statistics (for now macro stress tensor in the whole domain)
      DOMNode* nForceStats = ReaderXML::getNode( nPostProcessing,
      	"ForceStats" );
      if ( nForceStats )
      {
        size_t FSoutputFreq = size_t( 
		ReaderXML::getNodeAttr_Int( nForceStats, "Every" ) );
        string ppObsdir = ReaderXML::getNodeAttr_String( nForceStats,
		"Directory" );
        m_collision->setForceStatsParameters( ppObsdir, FSoutputFreq );		

	if ( m_rank == 0 )
	{
	  cout << GrainsExec::m_shift9 << "Force stats in the whole domain" 
	  	<< endl;
          cout << GrainsExec::m_shift12 << "Write values in file every " <<
		FSoutputFreq << " time step" <<
		( FSoutputFreq > 1 ? "s" : "" ) << endl;
          cout << GrainsExec::m_shift12 << "Output file directory name = "
    		<< ppObsdir << endl;
	}
      }      
    }
    else
      if ( m_rank == 0 ) cout << GrainsExec::m_shift6
      	<< "No postprocessing" << endl;
  }
}




// ----------------------------------------------------------------------------
// Creates, inserts and links new particles in the simulation
void Grains::InsertCreateNewParticles()
{
  // IMPORTANT: for any system with N particles and M obstacles, particles are
  // numbered 1 to N and obstacles are numbered -1 to -M

  // Link all components with the grid
  m_allcomponents.Link( *m_collision );

  // In case of a restarted simulation, if the linked cell changed from the 
  // previous simulation, we need to check that all periodic clones are there
  if ( m_restart && m_periodic ) checkClonesReload();

  // Set particle positions from file or from a structured array
  size_t error = 0;
  if ( m_position != "" )
  {
    // From a structured array
    if ( m_position == "STRUCTURED" ) error = setPositionParticlesArray();
    // From a file
    else error = setPositionParticlesFromFile();
  }
  if ( error ) grainsAbort();
  
  // Set angular particle positions from file
  if ( m_angular_position != "" )
    error = setAngularPositionParticlesFromFile();
  if ( error ) grainsAbort();  

  // Insertion at initial time
  if ( m_insertion_mode == IM_INITIALTIME )
  {
    cout << "Particles \tIn \tOut" << endl;
    while ( m_allcomponents.getNumberPhysicalParticlesToInsert() )
    {
        insertParticle( m_insertion_order );
        cout << "\r                                              " << flush;
        cout << "\r\t\t"
	   << m_allcomponents.getNumberActiveParticlesOnProc() << '\t'
	   << m_allcomponents.getNumberPhysicalParticlesToInsert() 
	   << flush;
    }
    cout << endl;
  }

  if ( m_rank == 0 )
    cout << "Particle volume IN  : " <<
    	m_allcomponents.getVolumeIn()  << endl
	<< "                OUT : " <<
	m_allcomponents.getVolumeOut() << endl;
}




// ----------------------------------------------------------------------------
// Returns a point randomly selected in one of the insertion windows or in 
// the list of positions. Note: list of positions has priority over windows 
// until it is empty
Point3 Grains::getInsertionPoint()
{
  Point3 P;
  
  // List of positions
  if ( !m_insertion_position->empty() )
  {
    // If order, pick the first position
    m_il_sp = m_insertion_position->begin();    
    
    // If random, pick a random position in the list
    if ( m_insertion_order == PM_RANDOM )  
    {
      double n = double(random()) / double(INT_MAX);
      size_t shift = size_t( n * double(m_insertion_position->size()) );
      std::advance( m_il_sp, shift );      
    }
    P = *m_il_sp; 
  }
  // Insertion windows
  else
  {
    int nWindow = 0;
    int nbreWindows = int( m_insertion_windows.size() );

    // Random selection of an insertion window
    if ( nbreWindows != 1 )
    {
      double n = double(random()) / double(INT_MAX);
      nWindow = int( n * nbreWindows );
      if ( nWindow == nbreWindows ) nWindow--;
    }

    // Random position in the selected insertion window
    double r = 0., theta = 0., axiscoor = 0.;
    switch( m_insertion_windows[nWindow].ftype )
    {
      case WINDOW_BOX:
        P[X] = m_insertion_windows[nWindow].ptA[X]
      		+ ( double(random()) / double(INT_MAX) )
		* ( m_insertion_windows[nWindow].ptB[X]
			- m_insertion_windows[nWindow].ptA[X] );
        P[Y] = m_insertion_windows[nWindow].ptA[Y]
      		+ ( double(random()) / double(INT_MAX) )
		* ( m_insertion_windows[nWindow].ptB[Y]
			- m_insertion_windows[nWindow].ptA[Y] );
      	P[Z] = m_insertion_windows[nWindow].ptA[Z]
      		+ ( double(random()) / double(INT_MAX) )
		* ( m_insertion_windows[nWindow].ptB[Z]
			- m_insertion_windows[nWindow].ptA[Z] );
        break;

      case WINDOW_CYLINDER:
        r = ( double(random()) / double(INT_MAX) )
      		* m_insertion_windows[nWindow].radius;
        theta = ( double(random()) / double(INT_MAX) ) * 2. * PI;
        axiscoor = ( double(random()) / double(INT_MAX) )
      		* m_insertion_windows[nWindow].height;
        switch ( m_insertion_windows[nWindow].axisdir )
        {
          case X:
	    P[X] = m_insertion_windows[nWindow].ptA[X] + axiscoor;
	    P[Y] = m_insertion_windows[nWindow].ptA[Y] + r * cos( theta );
	    P[Z] = m_insertion_windows[nWindow].ptA[Z] + r * sin( theta );
	    break;

          case Y:
	    P[X] = m_insertion_windows[nWindow].ptA[X] + r * sin( theta );
	    P[Y] = m_insertion_windows[nWindow].ptA[Y] + axiscoor;
	    P[Z] = m_insertion_windows[nWindow].ptA[Z] + r * cos( theta );
	    break;

          default:
	    P[X] = m_insertion_windows[nWindow].ptA[X] + r * cos( theta );
	    P[Y] = m_insertion_windows[nWindow].ptA[Y] + r * sin( theta );
	    P[Z] = m_insertion_windows[nWindow].ptA[Z] + axiscoor;
	    break;
        }
        break;

      case WINDOW_ANNULUS:
        r = m_insertion_windows[nWindow].radius_int
      		+ ( double(random()) / double(INT_MAX) )
      		* ( m_insertion_windows[nWindow].radius
			- m_insertion_windows[nWindow].radius_int );
        theta = ( double(random()) / double(INT_MAX) ) * 2. * PI;
        axiscoor = ( double(random()) / double(INT_MAX) )
      		* m_insertion_windows[nWindow].height;
        switch ( m_insertion_windows[nWindow].axisdir )
        {
          case X:
	    P[X] = m_insertion_windows[nWindow].ptA[X] + axiscoor;
	    P[Y] = m_insertion_windows[nWindow].ptA[Y] + r * cos( theta );
	    P[Z] = m_insertion_windows[nWindow].ptA[Z] + r * sin( theta );
	    break;

          case Y:
	    P[X] = m_insertion_windows[nWindow].ptA[X] + r * sin( theta );
	    P[Y] = m_insertion_windows[nWindow].ptA[Y] + axiscoor;
	    P[Z] = m_insertion_windows[nWindow].ptA[Z] + r * cos( theta );
	    break;

          default:
	    P[X] = m_insertion_windows[nWindow].ptA[X] + r * cos( theta );
	    P[Y] = m_insertion_windows[nWindow].ptA[Y] + r * sin( theta );
	    P[Z] = m_insertion_windows[nWindow].ptA[Z] + axiscoor;
	    break;
        }
        break;

      case WINDOW_LINE:
        P = m_insertion_windows[nWindow].ptA
      		+ ( double(random()) / double(INT_MAX) )
      		* ( m_insertion_windows[nWindow].ptB
			- m_insertion_windows[nWindow].ptA );
        break;

      default:
        break;
    }
  }

  return ( P );
}




// ----------------------------------------------------------------------------
// Attempts to insert a particle in the simulation
bool Grains::insertParticle( PullMode const& mode )
{
  static size_t insert_counter = 0;
  bool insert = true;
  Vector3 vtrans, vrot ;
  Quaternion qrot;
  size_t npositions = m_insertion_position->size();
  size_t nangpositions = m_insertion_angular_position->size();  

  if ( insert_counter == 0 )
  {
    Particle *particle = m_allcomponents.getParticleToInsert( mode );
    if ( particle )
    {
      // Initialisation of the centre of mass position of the particle
      particle->setPosition( getInsertionPoint() );

      // Initialisation of the angular position of the particle
      // Rem: we compose to the right by a pure rotation as the particle
      // already has a non-zero position that we do not want to change (and
      // that we would change if we would compose to the left)
      if ( m_init_angpos != IAP_FIXED )
      {
        Transform trot;
	if ( m_init_angpos == IAP_RANDOM ) 
          trot.setBasis( GrainsExec::RandomRotationMatrix( m_dimension ) );
	else // m_init_angpos == IAP_FILE
	{
	  m_il_sap = m_insertion_angular_position->begin();
	  trot.setBasis( *m_il_sap );
	}
        particle->composePositionRightByTransform( trot );
      }

      // Initialisation of the particle velocity
      computeInitVit( vtrans, vrot );
      particle->setTranslationalVelocity( vtrans );
      particle->setAngularVelocity( vrot );

      // Transform and quaternion
      particle->initialize_transformWithCrust_to_notComputed();
      qrot.setQuaternion( particle->getRigidBody()->getTransform()
		->getBasis() );
      particle->setQuaternionRotation( qrot );      

      // If insertion if successful, shift particle from wait to inserted
      insert = m_collision->insertParticleSerial( particle,
      	m_allcomponents.getActiveParticles(),
	m_allcomponents.getPeriodicCloneParticles(),
	m_allcomponents.getReferenceParticles(),
	m_periodic, m_force_insertion );
      if ( insert )
      {
        m_allcomponents.WaitToActive();
	particle->InitializeForce( true );
	if ( npositions ) m_insertion_position->erase( m_il_sp );
	if ( nangpositions ) m_insertion_angular_position->erase( m_il_sap );	
      }
    }
  }

  ++insert_counter;
  if ( insert_counter == m_insertion_frequency ) insert_counter = 0;

  return ( insert || m_force_insertion );
}




// ----------------------------------------------------------------------------
// Computes the initial velocity of a particle
void Grains::computeInitVit( Vector3& vtrans, Vector3& vrot ) const
{
  switch( m_initvit_mode )
  {
    case IV_ZERO :
      vtrans.reset();
      vrot.reset();
      break;

    case IV_CONSTANT :
      vtrans = m_InitVtrans;
      vrot = m_InitVrot;
      break;

    case IV_RANDOM :
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
	vrot[Z] = m_RandomMotionCoefRot * (
	      	2. * (double(random()) / double(INT_MAX)) - 1. ) ;	
	if ( m_dimension == 3 )
	{
	  vrot[X] = m_RandomMotionCoefRot * (
	    	2. * (double(random()) / double(INT_MAX)) - 1. ) ;
	  vrot[Y] = m_RandomMotionCoefRot * (
	    	2. * (double(random()) / double(INT_MAX)) - 1. ) ;
	}
	else vrot[X] = vrot[Y] = 0.;
      }
      else vrot.reset();
      break;

    default :
      vtrans.reset();
      vrot.reset();
      break;
  }
}




// ----------------------------------------------------------------------------
// Sets particle initial positions from a file
size_t Grains::setPositionParticlesFromFile()
{
  size_t nnewpart = 0, npos = 0, error = 0;
  list< pair<Particle*,size_t> >::const_iterator il;
  Point3 position;
  ifstream filePos( m_position.c_str() );

  // Check that the file exists
  if ( !filePos.is_open() )
  {
    cout << "ERR : File " << m_position << " does not exist !" << endl;
    error = 1;
  }

  // Check that the number of positions in the file equals the number of 
  // new particles to insert
  if ( !error )
  {
    for (il=m_newParticles.cbegin();il!=m_newParticles.cend();il++)
      nnewpart += il->second;
    while ( filePos >> position ) ++npos;
    filePos.close();
    if ( nnewpart != npos )
    {
      cout << "ERR: number of new particles to insert is different from the"
    	<< " number of positions in file " << m_position << " : "
	<< nnewpart << " != " << npos << endl;
      filePos.close(); 
      error = 1;
    }
    else
    {
      // Read positions from the file and add them to the list
      filePos.open( m_position.c_str() );
      while ( filePos >> position ) m_insertion_position->push_back( position );
      filePos.close();
    }
  }

  return ( error );
}




// ----------------------------------------------------------------------------
// Sets angular particle initial positions from a file
size_t Grains::setAngularPositionParticlesFromFile()
{
  size_t nnewpart = 0, npos = 0, error = 0;
  list< pair<Particle*,size_t> >::const_iterator il;
  Matrix mat;
  ifstream fileAngPos( m_angular_position.c_str() );

  // Check that the file exists
  if ( !fileAngPos.is_open() )
  {
    cout << "ERR : File " << m_angular_position << " does not exist !" << endl;
    error = 1;
  }

  // Check that the number of positions in the file equals the number of 
  // new particles to insert
  if ( !error )
  {
    for (il=m_newParticles.cbegin();il!=m_newParticles.cend();il++)
      nnewpart += il->second;
    while ( fileAngPos >> mat ) ++npos;
    fileAngPos.close();
    if ( nnewpart != npos )
    {
      cout << "ERR: number of new particles to insert is different from the"
    	<< " number of angular positions in file " << m_angular_position 
	<< " : " << nnewpart << " != " << npos << endl;
      fileAngPos.close(); 
      error = 1;
    }
    else
    {
      // Read angular positions from the file and add them to the list
      fileAngPos.open( m_angular_position.c_str() ); 
      while ( fileAngPos >> mat ) 
        m_insertion_angular_position->push_back( mat );
      fileAngPos.close();
    }
  }
  
  return ( error );
}



// ----------------------------------------------------------------------------
// Sets particle initial position with a structured array
size_t Grains::setPositionParticlesArray()
{
  size_t k, l, m, nnewpart = 0, error = 0;
  list< pair<Particle*,size_t> >::const_iterator il;
  Point3 position;
  double deltax = ( m_InsertionArray->box.ptB[X]
  	- m_InsertionArray->box.ptA[X] ) / double(m_InsertionArray->NX) ;
  double deltay = ( m_InsertionArray->box.ptB[Y]
  	- m_InsertionArray->box.ptA[Y] ) / double(m_InsertionArray->NY) ;
  double deltaz = ( m_InsertionArray->box.ptB[Z]
  	- m_InsertionArray->box.ptA[Z] ) / double(m_InsertionArray->NZ) ;

  // Checks that number of positions equals the number of new particles 
  // to insert
  for (il=m_newParticles.cbegin();il!=m_newParticles.cend();il++)
    nnewpart += il->second;
  if ( nnewpart != m_InsertionArray->NX * m_InsertionArray->NY
  	* m_InsertionArray->NZ )
  {
    cout << "ERR: number of new particles to insert is different from the"
    	<< " number of positions in the structured array: " << nnewpart << 
	" != " << m_InsertionArray->NX * m_InsertionArray->NY 
	* m_InsertionArray->NZ << endl;
    error = 1;
  }
  else
  {
    // Compute positions and add them to the list
    for (k=0;k<m_InsertionArray->NX;++k)
      for (l=0;l<m_InsertionArray->NY;++l)
        for (m=0;m<m_InsertionArray->NZ;++m)
        {
          position[X] = m_InsertionArray->box.ptA[X]
		+ ( double(k) + 0.5 ) * deltax;
          position[Y] = m_InsertionArray->box.ptA[Y]
		+ ( double(l) + 0.5 ) * deltay;
          position[Z] = m_InsertionArray->box.ptA[Z]
		+ ( double(m) + 0.5 ) * deltaz;
          m_insertion_position->push_back( position );
        }
  }
  
  return ( error );  
}




// ----------------------------------------------------------------------------
// Sets the time algorithm to predictor or corrector mode
void Grains::setMode( const bool &predictor )
{
  m_predictor_mode = predictor;
}




// ----------------------------------------------------------------------------
// Returns whether the time algorithm is in predictor mode
bool Grains::isModePredictor()
{
  return ( m_predictor_mode );
}




// ----------------------------------------------------------------------------
// Reads data for MPI simulations and creates and sets the MPI wrapper
void Grains::readDomainDecomposition( DOMNode* root,
  	double const& lx, double const& ly, double const& lz )
{
  // Nothing to do in serial mode
  m_processorIsActive = true;
}




// ----------------------------------------------------------------------------
// Sets the linked cell grid
void Grains::defineLinkedCell( double const& radius, string const& oshift )
{
  size_t error = m_collision->set( 2. * radius, oshift );
  if ( error ) grainsAbort();
}




// ----------------------------------------------------------------------------
// Emergency termination in case of an issue
void Grains::grainsAbort() const
{
  exit(1);
}




// ----------------------------------------------------------------------------
// Writes reload files
void Grains::saveReload( double const& time )
{
  static size_t reload_counter =
  	GrainsExec::m_reloadFile_suffix == "A" ? 1 : 0 ;
  string reload_suffix = reload_counter ? "B" : "A" ;
  string reload = m_fileSave + reload_suffix ;
  reload = GrainsExec::fullResultFileName( reload, false );

  // Save current reload file
  ofstream result;
  if ( m_rank == 0 )
  {
    result.open( reload.c_str(), ios::out );
    result << "__Format2024__" << endl;  
    result << "#Time " << time << endl;
  }
  ContactBuilderFactory::save( result, m_fileSave, m_rank );
  if ( m_nprocs > 1 && GrainsExec::m_SaveMPIInASingleFile )
    m_allcomponents.write_singleMPIFile( result, reload, m_insertion_position,
    	m_clonesInReloadFile, m_collision, m_wrapper, m_periodic );
  else
    m_allcomponents.write( result, reload, m_insertion_position, 
    	m_clonesInReloadFile, m_collision, m_rank, m_nprocs, m_wrapper );
  if ( m_rank == 0 )
  {  
    result << "#EndTime" << endl;
    result.close();

    // Check that all reload files are in the same directory
    // Add one line to RFTable.txt
    GrainsExec::checkAllFilesForReload();
    ofstream RFT_out( ( m_fileSave + "_RFTable.txt" ).c_str(), ios::app );
    RFT_out << time << " " << m_fileSave + reload_suffix << endl;
    RFT_out.close();
  }

  ++reload_counter;
  if ( reload_counter == 2 ) reload_counter = 0;
}




// ----------------------------------------------------------------------------
// Returns the maximum particle ID number
int Grains::getMaxParticleIDnumber() const
{
  return ( m_allcomponents.getMaxParticleIDnumber() );
}




// ----------------------------------------------------------------------------
// Deletes .result and .xml result files
void Grains::clearResultXmlFiles() const
{
  if ( m_rank == 0 )
  {
    string cmd = "bash " + GrainsExec::m_GRAINS_HOME
     	+ "/Tools/ExecScripts/init_clear.exec " + m_fileSave;
    GrainsExec::m_return_syscmd = system( cmd.c_str() );
  }

}




// ----------------------------------------------------------------------------
// Displays the memory used by the simulation
void Grains::display_used_memory() const
{
  cout << "Memory used by Grains3D = ";
  GrainsExec::display_memory( cout, GrainsExec::used_memory() );
  cout << endl;
}




// ----------------------------------------------------------------------------
// Synchronizes the PPWindow boolean relative to each sub-domain
void Grains::synchronize_PPWindow()
{}




// ----------------------------------------------------------------------------
// Reads a window
void Grains::readWindow( DOMNode* nWindow, Window& iwindow,
	string const& oshift )
{
  // Insertion window type
  string iwindow_type = "Box";
  if ( ReaderXML::hasNodeAttr( nWindow, "Type" ) )
    iwindow_type = ReaderXML::getNodeAttr_String( nWindow,
		  	"Type" );

  if ( iwindow_type == "Cylinder" )
    iwindow.ftype = WINDOW_CYLINDER;
  else if ( iwindow_type == "Annulus" )
    iwindow.ftype = WINDOW_ANNULUS;
  else if ( iwindow_type == "Line" )
    iwindow.ftype = WINDOW_LINE;
  else if ( iwindow_type == "Box" )
    iwindow.ftype = WINDOW_BOX;
  else iwindow.ftype = WINDOW_NONE;

  DOMNodeList* points = NULL;
  DOMNode* pointA = NULL;
  DOMNode* pointB = NULL;
  DOMNode* cylGeom = NULL;
  DOMNode* annGeom = NULL;
  string axisdir_str = "W";

  // Read features depending on the type
  switch( iwindow.ftype )
  {
    case WINDOW_BOX:
      points = ReaderXML::getNodes( nWindow );
      pointA = points->item( 0 );
      pointB = points->item( 1 );

      iwindow.radius = iwindow.radius_int = iwindow.height = 0. ;
      iwindow.axisdir = NONE ;
      iwindow.ptA[X] = ReaderXML::getNodeAttr_Double( pointA, "X" );
      iwindow.ptA[Y] = ReaderXML::getNodeAttr_Double( pointA, "Y" );
      iwindow.ptA[Z] = ReaderXML::getNodeAttr_Double( pointA, "Z" );
      iwindow.ptB[X] = ReaderXML::getNodeAttr_Double( pointB, "X" );
      iwindow.ptB[Y] = ReaderXML::getNodeAttr_Double( pointB, "Y" );
      iwindow.ptB[Z] = ReaderXML::getNodeAttr_Double( pointB, "Z" );

      if ( m_rank == 0 )
      {
        cout << oshift << "Type = " << iwindow_type << endl;
	cout << oshift << GrainsExec::m_shift3 << "Point3 min = " <<
		iwindow.ptA[X] << " " << iwindow.ptA[Y] << " " <<
		iwindow.ptA[Z] << endl;
        cout << oshift << GrainsExec::m_shift3 << "Point3 max = " <<
		iwindow.ptB[X] << " " << iwindow.ptB[Y] << " " <<
		iwindow.ptB[Z] << endl;
      }
      break;

    case WINDOW_CYLINDER:
      pointA = ReaderXML::getNode( nWindow, "BottomCentre" );
      iwindow.ptA[X] = ReaderXML::getNodeAttr_Double( pointA, "X" );
      iwindow.ptA[Y] = ReaderXML::getNodeAttr_Double( pointA, "Y" );
      iwindow.ptA[Z] = ReaderXML::getNodeAttr_Double( pointA, "Z" );
      iwindow.ptB[X] = iwindow.ptB[Y] = iwindow.ptB[Z] = 0.;
      cylGeom = ReaderXML::getNode( nWindow, "Cylinder" );
      iwindow.radius = ReaderXML::getNodeAttr_Double( cylGeom, "Radius" );
      iwindow.radius_int = 0.;
      iwindow.height = ReaderXML::getNodeAttr_Double( cylGeom, "Height" );
      axisdir_str = ReaderXML::getNodeAttr_String( cylGeom, "Direction" );
      if ( axisdir_str == "X" ) iwindow.axisdir = X;
      else if ( axisdir_str == "Y" ) iwindow.axisdir = Y;
      else if ( axisdir_str == "Z" ) iwindow.axisdir = Z;
      else
      {
	if ( m_rank == 0 )
          cout << "Wrong axis direction in cylindrical "
		<< "insertion window; values: X, Y or Z" << endl;
        grainsAbort();
      }

      if ( m_rank == 0 )
      {
	cout << oshift << "Type = " << iwindow_type << endl;
        cout << oshift << GrainsExec::m_shift3 << "Bottom centre = " <<
		iwindow.ptA[X] << " " << iwindow.ptA[Y] << " " <<
		iwindow.ptA[Z] << endl;
	cout << oshift << GrainsExec::m_shift3 << "Radius = " <<
	  	iwindow.radius << endl;
        cout << oshift << GrainsExec::m_shift3 << "Height = " <<
		iwindow.height << endl;
        cout << oshift << GrainsExec::m_shift3 << "Direction = " <<
		axisdir_str << endl;
      }
      break;

    case WINDOW_ANNULUS:
      pointA = ReaderXML::getNode( nWindow, "BottomCentre" );
      iwindow.ptA[X] = ReaderXML::getNodeAttr_Double( pointA, "X" );
      iwindow.ptA[Y] = ReaderXML::getNodeAttr_Double( pointA, "Y" );
      iwindow.ptA[Z] = ReaderXML::getNodeAttr_Double( pointA, "Z" );
      iwindow.ptB[X] = iwindow.ptB[Y] = iwindow.ptB[Z] = 0.;
      annGeom = ReaderXML::getNode( nWindow, "Annulus" );
      iwindow.radius = ReaderXML::getNodeAttr_Double( annGeom,
	  	"RadiusExt" );
      iwindow.radius_int = ReaderXML::getNodeAttr_Double( annGeom,
		"RadiusInt" );
      iwindow.height = ReaderXML::getNodeAttr_Double( annGeom,
		"Height" );
      axisdir_str = ReaderXML::getNodeAttr_String( annGeom,
		"Direction" );
      if ( axisdir_str == "X" ) iwindow.axisdir = X;
      else if ( axisdir_str == "Y" ) iwindow.axisdir = Y;
      else if ( axisdir_str == "Z" ) iwindow.axisdir = Z;
      else
      {
        if ( m_rank == 0 )
          cout << "Wrong axis direction in cylindrical "
		<< "insertion window; values: X, Y or Z" << endl;
        grainsAbort();
      }

      if ( m_rank == 0 )
      {
	cout << oshift << "Type = " << iwindow_type << endl;
        cout << oshift << GrainsExec::m_shift3 << "Bottom centre = " <<
		iwindow.ptA[X] << " " << iwindow.ptA[Y] << " " <<
		iwindow.ptA[Z] << endl;
	cout << oshift << GrainsExec::m_shift3 << "External radius = " <<
		iwindow.radius << endl;
        cout << oshift << GrainsExec::m_shift3 << "Internal radius = " <<
		iwindow.radius_int << endl;
        cout << oshift << GrainsExec::m_shift3 << "Height = " <<
		iwindow.height << endl;
        cout << oshift << GrainsExec::m_shift3 << "Direction = " <<
		axisdir_str << endl;
      }
      break;

    case WINDOW_LINE:
      points = ReaderXML::getNodes( nWindow );
      pointA = points->item( 0 );
      pointB = points->item( 1 );

      iwindow.radius = iwindow.radius_int = iwindow.height = 0. ;
      iwindow.axisdir = NONE ;
      iwindow.ptA[X] = ReaderXML::getNodeAttr_Double( pointA, "X" );
      iwindow.ptA[Y] = ReaderXML::getNodeAttr_Double( pointA, "Y" );
      iwindow.ptA[Z] = ReaderXML::getNodeAttr_Double( pointA, "Z" );
      iwindow.ptB[X] = ReaderXML::getNodeAttr_Double( pointB, "X" );
      iwindow.ptB[Y] = ReaderXML::getNodeAttr_Double( pointB, "Y" );
      iwindow.ptB[Z] = ReaderXML::getNodeAttr_Double( pointB, "Z" );

      if ( m_rank == 0 )
      {
        cout << oshift << "Type = " << iwindow_type << endl;
        cout << oshift << GrainsExec::m_shift3 << "Point3 A = " <<
		iwindow.ptA[X] << " " << iwindow.ptA[Y] << " " <<
		iwindow.ptA[Z] << endl;
        cout << oshift << GrainsExec::m_shift3 << "Point3 B = " <<
		iwindow.ptB[X] << " " << iwindow.ptB[Y] << " " <<
		iwindow.ptB[Z] << endl;
      }
      break;

    default:
      if ( m_rank == 0 ) cout << "Unknown insertion window "
		"type" << endl;
      grainsAbort();
      break;
  }
}




// ----------------------------------------------------------------------------
// Outputs timer summary
void Grains::display_timer_summary()
{
  double cputime = CT_get_elapsed_time();
  cout << "Full problem" << endl;
  write_elapsed_time_smhd( cout, cputime, "Computing time" );
  SCT_get_summary( cout, cputime );
}




// ----------------------------------------------------------------------------
// Returns the number of insertion positions */
size_t Grains::getNbInsertionPositions() const
{
  return ( m_insertion_position->size() );
}




// ----------------------------------------------------------------------------
// Checks the (periodic) clones when a simulation is reloaded */
void Grains::checkClonesReload()
{
  m_collision->checkPeriodicClonesReload(
	m_allcomponents.getActiveParticles(),
	m_allcomponents.getPeriodicCloneParticles(),
	m_allcomponents.getReferenceParticles(),
	m_time );
}
