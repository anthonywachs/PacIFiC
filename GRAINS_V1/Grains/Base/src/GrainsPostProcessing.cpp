#include "GrainsMPIWrapper.hh"
#include "GrainsExec.hh"
#include "GrainsPostProcessing.hh"
#include "ContactBuilderFactory.hh"
#include "PostProcessingWriter.hh"
#include "PostProcessingWriterBuilderFactory.hh"
#include "RawDataPostProcessingWriter.hh"
#include "Sphere.hh"
#include "Disc.hh"


// ----------------------------------------------------------------------------
// Constructeur par defaut
GrainsPostProcessing::GrainsPostProcessing() 
  : Grains()
{
  m_global_porosity = NULL;
}




// ----------------------------------------------------------------------------
// Destructor
GrainsPostProcessing::~GrainsPostProcessing()
{}




// ----------------------------------------------------------------------------
// Writes an initial message in the standard output only on the process ranked 0
void GrainsPostProcessing::initialOutputMessage()
{
  int rankproc = 0;
  MPI_Comm_rank( MPI_COMM_WORLD, &rankproc );
  if ( !rankproc )
    cout << "GrainsPostProcessing serial" << endl;
}




// ----------------------------------------------------------------------------
// Tasks to perform before time-stepping 
void GrainsPostProcessing::do_before_time_stepping( DOMElement* rootElement )
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

  // Initialisation obstacle kinematics
  m_allcomponents.setKinematicsObstacleWithoutMoving( m_time, m_dt ); 

  cout << "Initialization completed" << endl << endl;                           
}




// ----------------------------------------------------------------------------
// Tasks to perform before time-stepping 
void GrainsPostProcessing::do_after_time_stepping()
{  
  if ( !m_rank )
    cout << "PostProcessing completed" << endl; 	 
}




// ----------------------------------------------------------------------------
// Runs the simulation over the prescribed time interval
void GrainsPostProcessing::Simulation( double time_interval )
{
  // Global porosity
  if ( m_global_porosity )
  {    
    double volparticles = 0.;
    double volporodomain = ( m_global_porosity->domain.ptB[X] 
    	- m_global_porosity->domain.ptA[X] ) *
	( m_global_porosity->domain.ptB[Y] 
    	- m_global_porosity->domain.ptA[Y] ) *
	( m_global_porosity->domain.ptB[Z] 
    	- m_global_porosity->domain.ptA[Z] );
    double dx = ( m_global_porosity->domain.ptB[X] 
    	- m_global_porosity->domain.ptA[X] ) 
	/ double(m_global_porosity->nintervals[X]);
    double dy = ( m_global_porosity->domain.ptB[Y] 
    	- m_global_porosity->domain.ptA[Y] ) 
	/ double(m_global_porosity->nintervals[Y]);
    double dz = ( m_global_porosity->domain.ptB[Z] 
    	- m_global_porosity->domain.ptA[Z] ) 
	/ double(m_global_porosity->nintervals[Z]);
    double dv = dx * dy * dz;
    Point3 elemVolCenter;
    list<Cell*> cells;
    list<Cell*>::const_iterator ic;
    bool found = false;
    
    for (size_t i=0;i<m_global_porosity->nintervals[X];++i)
      for (size_t j=0;j<m_global_porosity->nintervals[Y];++j)    
        for (size_t k=0;k<m_global_porosity->nintervals[Z];++k)
	{
	  // Coordinates of the center of the elementary volume
	  elemVolCenter[X] = m_global_porosity->domain.ptA[X]
	  	+ ( double(i) + 0.5 ) * dx;
	  elemVolCenter[Y] = m_global_porosity->domain.ptA[Y]
	  	+ ( double(j) + 0.5 ) * dy;	  
	  elemVolCenter[Z] = m_global_porosity->domain.ptA[Z]
	  	+ ( double(k) + 0.5 ) * dz;	  
	  
	  // List of cells containing the cell that elemVolCenter belongs
	  // to and its neighboring cells
	  cells = m_collision->getCellAndCellNeighborhood( elemVolCenter );
	  
	  // Loops over the cells and check whether elemVolCenter belongs
	  // to any particle in any of these cells
	  found = false;
	  for (ic=cells.begin();ic!=cells.end() && !found;ic++)
	    found = (*ic)->isInParticle( elemVolCenter );
	    
          if ( found ) volparticles += dv;	  
	}
         
    double porosity = ( volporodomain - volparticles ) / volporodomain;    
    cout << "Average porosity = " << porosity << endl;
    cout << endl;    	
  }
}




// ----------------------------------------------------------------------------
// Construction of the simulation: linked cell, particles & obstacles, domain 
// decomposition 
void GrainsPostProcessing::Construction( DOMElement* rootElement )
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

      // Restart mode is new, not read in XML file
      GrainsExec::m_ReloadType = "new" ;
      restart  = ReaderXML::getNodeAttr_String( reload, "Filename" );	
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
    
    // Link obstacles with the linked cell grid
    m_collision->Link( m_allcomponents.getObstacles() );     
  }              
}




// ----------------------------------------------------------------------------
// External force definition
void GrainsPostProcessing::Forces( DOMElement* rootElement )
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
void GrainsPostProcessing::AdditionalFeatures( DOMElement* rootElement )
{
  if ( m_processorIsActive )
  {
    assert( rootElement != NULL );
    DOMNode* root = ReaderXML::getNode( rootElement, "Simulation" );    
    

    // Check that Simulation node exists
    if ( !root )
    {
      cout << GrainsExec::m_shift6 << "<Simulation> node is mandatory !!" 
      		<< endl;
      grainsAbort();          
    }


    // Output message
    if ( m_rank == 0 ) cout << GrainsExec::m_shift3 << "PostProcessing" << endl;
    
    
    // Post-processing
    DOMNode* nPostProcessing = ReaderXML::getNode( root, 
    	"PostProcessing" );
    if ( nPostProcessing )
    {
      DOMNode* nGlobalPoro = ReaderXML::getNode( nPostProcessing, 
    	"GlobalPorosity" );
      if ( nGlobalPoro ) 
      {
        cout << GrainsExec::m_shift6 << "Global porosity" 
      		<< endl;
	m_global_porosity = new struct GlobalPorosity;	 

        // Domain features
        DOMNode* nWindow = ReaderXML::getNode( nGlobalPoro, "Window" ); 
        readWindow( nWindow, m_global_porosity->domain, GrainsExec::m_shift9 );
	
	// Number of intervals in each direction for numerical integration
	m_global_porosity->nintervals[0] = size_t(ReaderXML::getNodeAttr_Int( 
		nGlobalPoro, "N0" ));
	m_global_porosity->nintervals[1] = size_t(ReaderXML::getNodeAttr_Int( 
		nGlobalPoro, "N1" ));
	m_global_porosity->nintervals[2] = size_t(ReaderXML::getNodeAttr_Int( 
		nGlobalPoro, "N2" ));
	cout << GrainsExec::m_shift9 << "Discretization = " 
		<< m_global_porosity->nintervals[0] << " x "
		<< m_global_porosity->nintervals[1] << " x "
		<< m_global_porosity->nintervals[2] << endl;			
      }
    }
    else
      if ( m_rank == 0 ) cout << GrainsExec::m_shift6 
      	<< "No postprocessing" << endl;             
  }
}
