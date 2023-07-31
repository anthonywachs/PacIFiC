#include "GrainsMPIWrapper.hh"
#include "GrainsExec.hh"
#include "GrainsCoupledWithFluid.hh"
#include "ContactBuilderFactory.hh"
#include "PostProcessingWriter.hh"
#include "PostProcessingWriterBuilderFactory.hh"
#include "RawDataPostProcessingWriter.hh"
#include "Sphere.hh"
#include "Disc.hh"


// ----------------------------------------------------------------------------
// Default constructor
GrainsCoupledWithFluid::GrainsCoupledWithFluid( double fluid_density_ ) 
  : Grains()
  , m_fluid_density( fluid_density_ )
  , m_fluidflow_dt( 0. )
  , m_forceReloadSame( false )
  , m_PRSHydroFT( NULL )
{
  Disc::SetvisuNodeNb( 40 );
  Sphere::SetvisuNodeNbPerQar( 5 );
  Particle::setFluidDensity( m_fluid_density ); 
}




// ----------------------------------------------------------------------------
// Destructor
GrainsCoupledWithFluid::~GrainsCoupledWithFluid()
{}




// ----------------------------------------------------------------------------
// Writes an initial message in the standard output only on the process ranked 0
void GrainsCoupledWithFluid::initialOutputMessage()
{
  int rankproc = 0;
  MPI_Comm_rank( MPI_COMM_WORLD, &rankproc );
  if ( !rankproc )
    cout << "GrainsCoupledWithFluid serial" << endl;
}




// ----------------------------------------------------------------------------
// Tasks to perform before time-stepping 
void GrainsCoupledWithFluid::do_before_time_stepping( DOMElement* rootElement )
{
  // Read the input file
  Construction( rootElement );
  Forces( rootElement );
  AdditionalFeatures( rootElement ); 
  
  cout << endl << "Initialization" << endl;

  // Set time to initial time 
  m_time = m_tstart;

  // Link all components with the grid
  m_allcomponents.Link( *m_collision );

  // Number of particles: inserted and in the system
  m_allcomponents.setNumberParticlesOnAllProc( 
  	m_allcomponents.getNumberParticles() );
  m_npwait_nm1 = m_allcomponents.getNumberInactiveParticles(); 
      
  // Allocate hydro force and torque arrays in the AppPRSHydroFT app
  if ( m_PRSHydroFT ) m_PRSHydroFT->allocateHydroFT( 
  	m_allcomponents.getNumberActiveParticlesOnProc() );

  cout << "Initialization completed" << endl << endl;                           
}




// ----------------------------------------------------------------------------
// Tasks to perform before time-stepping 
void GrainsCoupledWithFluid::do_after_time_stepping()
{  
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
}




// ----------------------------------------------------------------------------
// Runs the simulation over the prescribed time interval
void GrainsCoupledWithFluid::Simulation( double time_interval )
{
  list<App*>::iterator app;
  
  // Compute the number of granular time steps
  m_fluidflow_dt = time_interval;
  double dteff = max( min ( m_max_dt, m_fluidflow_dt / double(m_ndt) ), 
  	m_min_dt );
  m_ndt = size_t( m_fluidflow_dt / dteff );
  m_dt = m_fluidflow_dt / double(m_ndt);  
//   cout << m_fluidflow_dt << " " << m_dt << " " << m_ndt << endl;
//   cout << m_time << endl;

  // Simulation: time marching algorithm
  for (size_t m=0;m<m_ndt;++m)
  {
    try 
    {        	
      m_time += m_dt;
	
      // Initiliaze all component transforms with crust to non computed
      m_allcomponents.InitializeRBTransformWithCrustState( m_time, m_dt );


      // Compute volume and contact forces
      // Initialisation torsors with weight only
      m_allcomponents.InitializeForces( m_time, m_dt, true );
      
      // Compute forces from all applications     
      for (app=m_allApp.begin(); app!=m_allApp.end(); app++)
        (*app)->ComputeForces( m_time, m_dt, 
      		m_allcomponents.getActiveParticles() );


      // Solve Newton's law and move particles
      m_allcomponents.Move( m_time, m_dt );
      
      
      // In case of periodicity, update periodic clones and destroy periodic
      // clones that are out of the linked cell grid
      if ( m_periodic )
        m_collision->updateDestroyPeriodicClones( 
		m_allcomponents.getActiveParticles(),
		m_allcomponents.getPeriodicCloneParticles() );
      
      	
      // Update particle activity
      m_allcomponents.UpdateParticleActivity();
 

      // Update the particle & obstacles links with the grid
      m_collision->LinkUpdate( m_time, m_dt, 
        	m_allcomponents.getActiveParticles() );


      // In case of periodicity, create new periodic clones and destroy periodic
      // clones because the master particle has changed tag/geoposition
      if ( m_periodic )
        m_collision->createDestroyPeriodicClones( 
		m_allcomponents.getActiveParticles(),
		m_allcomponents.getPeriodicCloneParticles(),
		m_allcomponents.getReferenceParticles() );
  
        
      // Write force & torque exerted on obstacles
      m_allcomponents.outputObstaclesLoad( m_time, m_dt );
    } 
    catch (ContactError &chocCroute) 
    {
      // Max overlap exceeded
      cout << endl;
      m_allcomponents.PostProcessingErreurComponents( "ContactError",
            chocCroute.getComponents() );
      chocCroute.Message( cout );
      m_error_occured = true;
      break;
    } 
    catch (DisplacementError &errDeplacement) 
    {
      // Particle displacement over dt is too large
      cout << endl;
      m_allcomponents.PostProcessingErreurComponents( "DisplacementError",
            errDeplacement.getComponent() );
      errDeplacement.Message(cout);
      m_error_occured = true;	
      break;
    } 
    catch (SimulationError &errSimulation) 
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
// Construction of the simulation: linked cell, particles & obstacles, domain 
// decomposition 
void GrainsCoupledWithFluid::Construction( DOMElement* rootElement )
{
  assert( rootElement != NULL );
  DOMNode* root = ReaderXML::getNode( rootElement, "Construction" );

  bool brestart = false;
  string restart;

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
      brestart = true;

      // Restart mode
      string reload_type;
      if ( m_forceReloadSame ) reload_type = "same" ;
      else reload_type = ReaderXML::getNodeAttr_String( reload, "Type" );      
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
      restart = fullResultFileName( restart );
      
      // Extract the reload directory from the reload file
      GrainsExec::m_ReloadDirectory = GrainsExec::extractRoot( restart ); 

      // Read the reload file
      string cle;
      ifstream simulLoad( restart.c_str() );
      simulLoad >> cle >> m_time;
      ContactBuilderFactory::reload( simulLoad );
      m_allcomponents.read( simulLoad, restart );
      ContactBuilderFactory::set_materialsForObstaclesOnly_reload(
          m_allcomponents.getReferenceParticles() );
      simulLoad >> cle;

      // Whether to reset velocity to 0
      string reset = ReaderXML::getNodeAttr_String( reload, "Velocity" );
      m_allcomponents.resetKinematics( reset );
    }          
   

    // Check that construction is fine
    if ( !brestart ) 
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
    
    // Link obstacles with the linked cell grid
    m_collision->Link( m_allcomponents.getObstacles() );     
  }              
}




// ----------------------------------------------------------------------------
// External force definition
void GrainsCoupledWithFluid::Forces( DOMElement* rootElement )
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
      
      // PRS hydro forces and torques
      DOMNode* nPRSHydroFT = ReaderXML::getNode( root, "PRSHydro" ); 
      if ( nPRSHydroFT )
      {
        m_PRSHydroFT = new AppPRSHydroFT();
	m_PRSHydroFT->setName( "PRSHydroFT" );
        m_allApp.push_back( m_PRSHydroFT );	
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
void GrainsCoupledWithFluid::AdditionalFeatures( DOMElement* rootElement )
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
    

    // Time step
    DOMNode* nTimeStep = ReaderXML::getNode( root, "TimeStep" );
    if ( nTimeStep )
    {
      m_min_dt = ReaderXML::getNodeAttr_Double( nTimeStep, "mindt" );
      m_max_dt = ReaderXML::getNodeAttr_Double( nTimeStep, "maxdt" );
      m_ndt = size_t( ReaderXML::getNodeAttr_Int( nTimeStep, "ndt" ) );      
      if ( m_rank == 0 ) 
      {
        cout << GrainsExec::m_shift6 << 
      		"Minimum time step magnitude = " << m_min_dt << endl;
        cout << GrainsExec::m_shift6 << 
      		"Maximum time step magnitude = " << m_max_dt << endl;		
        cout << GrainsExec::m_shift6 << 
      		"Number of time steps over a fluid time step = " << m_ndt 
		<< endl;
      }
      m_dt = m_min_dt;		
    }
    else
    {
      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 << 
		"Time step features are mandatory !!" << endl;
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
// Sets the initial physical time and initializes what depends 
// on the initial time
void GrainsCoupledWithFluid::setInitialTime( double const& time0 )
{
  // Initial time
  m_tstart = time0;
  m_time = time0;
  
  // Initialisation obstacle kinematics
  m_allcomponents.setKinematicsObstacleWithoutMoving( m_time, m_dt ); 

  // Postprocessing of force & torque on obstacles 
  m_allcomponents.initialiseOutputObstaclesLoadFiles( m_rank, false, m_time );  
  m_allcomponents.outputObstaclesLoad( m_time, m_dt, false, 
      GrainsExec::m_ReloadType == "same" );
}




// ----------------------------------------------------------------------------
// Sets the initial postprocessing cycle number
void GrainsCoupledWithFluid::setInitialCycleNumber( int const& cycle0 )
{
  if ( m_processorIsActive ) m_allcomponents.setInitialCycleNumber( cycle0 );
}




// ----------------------------------------------------------------------------
// Checks that the Paraview post-processing writer exists, otherwise
// creates it
void GrainsCoupledWithFluid::checkParaviewPostProcessing( string const& name_, 
	string const& root_,
	bool const& isBinary )
{
  if ( m_processorIsActive ) 
    m_allcomponents.checkParaviewPostProcessing( m_rank, m_nprocs, name_,
    	root_, isBinary );
}

void GrainsCoupledWithFluid::checkParaviewPostProcessing( char const* name_, 
	char const* root_,
  	bool const& isBinary )
{
  if ( m_processorIsActive ) 
  {
    string str_name_( name_ );
    string str_root_( root_ );
    m_allcomponents.checkParaviewPostProcessing( m_rank, m_nprocs, str_name_,
    	str_root_, isBinary );	
  }
}




// ----------------------------------------------------------------------------
// Sets the boolean m_forceReloadSame to true. This forces the code 
// to restart a simulation as a continuation of a previous simulation
void GrainsCoupledWithFluid::setReloadSame()
{
  m_forceReloadSame = true;
}




// ----------------------------------------------------------------------------
// Initially writes postprocessing files    
void GrainsCoupledWithFluid::InitialPostProcessing( size_t indent_width )
{
  if ( m_processorIsActive )
    m_allcomponents.PostProcessing_start( m_time, m_dt, m_collision, 
    	m_insertion_windows );
}




// ----------------------------------------------------------------------------
// Writes postprocessing files    
void GrainsCoupledWithFluid::doPostProcessing( size_t indent_width )
{
  if ( m_processorIsActive )
  {
    // Write postprocessing files    
    m_allcomponents.PostProcessing( m_time, m_dt, m_collision );

    // Write reload files
    saveReload( m_time );
  }
}	




// ----------------------------------------------------------------------------
// Writes features of moving rigid bodies in a stream to be used by the 
// fluid flow solver
void GrainsCoupledWithFluid::GrainsToFluid( istringstream &is ) const
{
  if ( m_processorIsActive )
  {
    list<Particle*> const* particles = m_allcomponents.getActiveParticles();
    list<Obstacle*> obstaclesToFluid = m_allcomponents.getObstaclesToFluid();
    ostringstream particles_features;
    list<Particle*>::const_iterator particle, clone;
    int componentIDinFluid = 0, particleID = 0, ncorners = 0;
    size_t nclonesper = 0, nparticles = 0;
    Vector3 const* vtrans = NULL;
    Vector3 const* vrot = NULL;
    Point3 const* centre = NULL;
    vector<double> inertia( 6, 0. );
    string particleType = "P";
    double density = 0., mass = 0., radius = 0.;
    Matrix mr;
    multimap<int,Particle*> const* particlesPeriodicClones = 
    	m_allcomponents.getPeriodicCloneParticles();
    multimap<int,Particle*>::const_iterator imm;
    pair < multimap<int,Particle*>::const_iterator,
  	multimap<int,Particle*>::const_iterator > crange;	
    Particle* periodic_clone = NULL ;
    Vector3 periodicVector;
    
    // TO DO: the particle number used to communicate with the fluid is 
    // different from its actual number in Grains. 
    // This needs to be fixed in the future 
    for (particle=particles->begin();particle!=particles->end();particle++)
      if ( (*particle)->getTag() < 2 ) nparticles++;	          
    
    // Total number of rigid bodies to send to the fluid flow solver
    particles_features << nparticles + obstaclesToFluid.size() << endl;

    // Particle features
    for (particle=particles->begin(), componentIDinFluid=0; 
      	particle!=particles->end();particle++)
    {
      if ( (*particle)->getTag() < 2 ) 
      {     
        if ( (*particle)->getActivity() == COMPUTE )      
        {        
	  vtrans = (*particle)->getTranslationalVelocity();
          vrot = (*particle)->getAngularVelocity();
          centre = (*particle)->getPosition();
          density = (*particle)->getDensity();
          mass = (*particle)->getMass();
          (*particle)->computeInertiaTensorSpaceFixed( inertia );
          radius = (*particle)->getCircumscribedRadius();
          ncorners = (*particle)->getNbCorners();
          particleID = (*particle)->getID();
	  nclonesper = particlesPeriodicClones->count( particleID ); 
          particleType = "P";
          if ( nclonesper ) particleType = "PP";
          mr = (*particle)->getRigidBody()->getTransform()->getBasis();

          particles_features << componentIDinFluid << " " << ncorners << endl;
          particles_features << particleType << " " <<
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			(*vtrans)[X] ) << " " << 
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			(*vtrans)[Y] ) << " " << 			
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			(*vtrans)[Z] ) << " " << 			
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			(*vrot)[X] ) << " " << 
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			(*vrot)[Y] ) << " " << 			
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			(*vrot)[Z] ) << " " << 
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			density ) << " " <<
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			mass ) << " " <<			 
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			inertia[0] ) << " " <<
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			inertia[1] ) << " " <<			
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			inertia[2] ) << " " <<
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			inertia[3] ) << " " <<
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			inertia[4] ) << " " <<			
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			inertia[5] ) << " " <<			
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			mr[X][X] ) << " " <<
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			mr[X][Y] ) << " " <<			
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			mr[X][Z] ) << " " <<	
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			mr[Y][X] ) << " " <<
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			mr[Y][Y] ) << " " <<			
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			mr[Y][Z] ) << " " <<			
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			mr[Z][X] ) << " " <<
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			mr[Z][Y] ) << " " <<			
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			mr[Z][Z] ) << " " <<			
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			(*centre)[X] ) << " " <<
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			(*centre)[Y] ) << " " <<				
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			(*centre)[Z] ) << endl;			

          if ( particleType == "PP" )
          {
            particles_features << nclonesper << endl;
	    crange = particlesPeriodicClones->equal_range( particleID );
	    for (imm=crange.first; imm!=crange.second;imm++)
	    {
	      periodic_clone = imm->second;
              periodicVector = *(periodic_clone->getPosition()) - *centre;
              particles_features << 
	      	GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			periodicVector[X] ) << " " <<
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			periodicVector[Y] ) << " " <<				
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			periodicVector[Z] ) << endl;	      
	    }
          }

          particles_features << radius ;
          (*particle)->writePositionInFluid( particles_features );
        }
        else
        {
          particles_features << componentIDinFluid << " " << "1" << endl;
          particles_features << "P " <<
		"0. 0. 0. 0. 0. 0. 1e8 1. " <<
		"1.  1.  1.  1.  1.  1. " <<
		"1. 0. 0. 0. 1. 0. 0. 0. 1." 
		"0. 0. 0. " << endl;
          particles_features << "0.  1"     << endl;
          particles_features << "0. 0. 0." << endl;
	  particles_features << "0" << endl;
        }
      
        componentIDinFluid++;
      }
    }

    // Obstacles to be sent to the fluid
    list<Obstacle*>::const_iterator obst;
    string obstacleType = "O";
    for (obst=obstaclesToFluid.begin();obst!=obstaclesToFluid.end();obst++,
    	componentIDinFluid++)
    {
      vtrans = (*obst)->getTranslationalVelocity();
      vrot = (*obst)->getAngularVelocity();
      centre = (*obst)->getPosition();
      radius = (*obst)->getCircumscribedRadius();
      ncorners = (*obst)->getRigidBody()->getConvex()->getNbCorners();
      mr = (*obst)->getRigidBody()->getTransform()->getBasis();

      particles_features << componentIDinFluid << " " << ncorners << endl;
      particles_features << obstacleType << " " <<
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			(*vtrans)[X] ) << " " << 
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			(*vtrans)[Y] ) << " " << 			
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			(*vtrans)[Z] ) << " " << 			
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			(*vrot)[X] ) << " " << 
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			(*vrot)[Y] ) << " " << 			
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			(*vrot)[Z] ) << " " << 
		"1000. 0. 0. 0. 0. 0. 0. 0. " <<		
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			mr[X][X] ) << " " <<
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			mr[X][Y] ) << " " <<			
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			mr[X][Z] ) << " " <<	
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			mr[Y][X] ) << " " <<
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			mr[Y][Y] ) << " " <<			
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			mr[Y][Z] ) << " " <<			
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			mr[Z][X] ) << " " <<
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			mr[Z][Y] ) << " " <<			
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			mr[Z][Z] ) << " " <<			
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			(*centre)[X] ) << " " <<
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			(*centre)[Y] ) << " " <<				
		GrainsExec::doubleToString( ios::scientific, POSITIONFORMAT,
			(*centre)[Z] ) << endl;			

      particles_features << radius ;
      (*obst)->writePositionInFluid( particles_features );    
    }

    // Transfer from oss to iss
    is.str( particles_features.rdbuf()->str() );     
  }  
}




// ----------------------------------------------------------------------------
// Updates particles velocity with data from the fluid solver
void GrainsCoupledWithFluid::updateParticlesVelocity( 
  	vector<vector<double> > const& velocity_data_array,
  	bool const& b_set_velocity_nm1_and_diff )
{
  if ( m_processorIsActive )
  {
    list<Particle*>* particles = m_allcomponents.getActiveParticles(); 
    list<Particle*>::iterator particle;
    int id = 0;
    size_t vecSize = 0;
    Vector3 vtrans, vrot;
    
    // TO DO: the particle number used to communicate with the fluid is 
    // different from its actual number in Grains. 
    // This needs to be fixed in the future

    if ( m_dimension == 3 )
    {
      if ( velocity_data_array.size() != particles->size() )
        cout << "WARNING: numbers of particles in Grains and in the fluid "
		<< "solver are different" << endl;

      for (particle=particles->begin(), id=0; particle!=particles->end();
       		particle++, id++)
      {
        if ( (*particle)->getActivity() == COMPUTE
    		&& (*particle)->getTag() < 2 )
        {
          vecSize = (velocity_data_array[id]).size();
          if ( vecSize != 6 )
	    cout << "ERROR: the velocity data array does not have the right "
		<< "size, it is " << vecSize << " while it must be 6" << endl;

          vtrans[X] = velocity_data_array[id][0];
          vtrans[Y] = velocity_data_array[id][1];
          vtrans[Z] = velocity_data_array[id][2];
          vrot[X] = velocity_data_array[id][3];
          vrot[Y] = velocity_data_array[id][4];
          vrot[Z] = velocity_data_array[id][5];

          (*particle)->setTranslationalVelocity( vtrans );
          (*particle)->setAngularVelocity( vrot );

          if ( b_set_velocity_nm1_and_diff )
            (*particle)->setVelocityAndVelocityDifferencePreviousTime();
        }
      }    
    }
    else
    {
      // TO DO
    }     
  }
}




// ----------------------------------------------------------------------------
// Updates particles hydro force and torque with data from the fluid solver
void GrainsCoupledWithFluid::updateParticlesHydroFT( 
  	vector< vector<double> > const* hydroft_data_array )
{
  m_PRSHydroFT->setHydroFT( hydroft_data_array );
} 




// ----------------------------------------------------------------------------
// Sets the boolean Particle::setFluidCorrectedAcceleration. Default
// value is True, i.e., the particle acceleration is corrected by
// the factor ( 1 - fluid_density / particle_density )
void GrainsCoupledWithFluid::setFluidCorrectedAcceleration( bool correct )
{
  Particle::setFluidCorrectedAcceleration( correct );
}




// ----------------------------------------------------------------------------
// Sets the Paraview post-processing translation vector in case of
// projection-translation
void GrainsCoupledWithFluid::setParaviewPostProcessingTranslationVector( 
      	double const& tvx, double const& tvy, double const& tvz )
{
  if ( !GrainsExec::m_translationParaviewPostProcessing )
    GrainsExec::m_translationParaviewPostProcessing = new Vector3();
  (*GrainsExec::m_translationParaviewPostProcessing)[X] = tvx;
  (*GrainsExec::m_translationParaviewPostProcessing)[Y] = tvy;
  (*GrainsExec::m_translationParaviewPostProcessing)[Z] = tvz;
}
