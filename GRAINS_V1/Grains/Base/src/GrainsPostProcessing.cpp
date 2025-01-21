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
  
  // In case of a periodic simulation, if the linked cell changed from the 
  // previous simulation or periodic clones were not saved in the restart file,
  // we need to check that all periodic clones are there
  if ( m_periodic ) checkClonesReload();  

  // Number of particles: inserted and in the system
  m_allcomponents.computeNumberParticles( m_wrapper );
  m_npwait_nm1 = m_allcomponents.getNumberPhysicalParticlesToInsert();

  // Initialisation of obstacle kinematics
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
    double volporodomain = 0.;
    double dx = 0., dy = 0., dz = 0.;
    Point3 elemVolCenter;
    list<Cell*> cells;
    list<Cell*>::const_iterator ic;
    bool found = false;
    list<Particle*> const* allpart = m_allcomponents.getActiveParticles();
    list<Particle*>::const_iterator il;      
    size_t nx, ny, nz;
    Point3 intVolptA;	
    Vector3 intVolExtent;
    double dxl, dyl, dzl, dvl; 
        
    if ( m_global_porosity->domain.ftype == WINDOW_BOX )
    {       
      volporodomain = ( m_global_porosity->domain.ptB[X] 
    	- m_global_porosity->domain.ptA[X] ) *
	( m_global_porosity->domain.ptB[Y] 
    	- m_global_porosity->domain.ptA[Y] ) *
	( m_global_porosity->domain.ptB[Z] 
    	- m_global_porosity->domain.ptA[Z] );
      dx = ( m_global_porosity->domain.ptB[X] 
    	- m_global_porosity->domain.ptA[X] ) 
	/ double(m_global_porosity->nintervals[X]);
      dy = ( m_global_porosity->domain.ptB[Y] 
    	- m_global_porosity->domain.ptA[Y] ) 
	/ double(m_global_porosity->nintervals[Y]);
      dz = ( m_global_porosity->domain.ptB[Z] 
    	- m_global_porosity->domain.ptA[Z] ) 
	/ double(m_global_porosity->nintervals[Z]);    
      BBox BBdomain( m_global_porosity->domain.ptA,
      	 m_global_porosity->domain.ptB ), BBintVol;
      
      for (il=allpart->begin();il!=allpart->end();il++)
      {
        BBox BBpart = (*il)->BoundingBox(); 
	
	// The particle is not fully contained, we pixelate the part
	// that belongs to the porosity domain
	if ( !BBdomain.fullyContain( BBpart ) )
	{
	  // We 1st check that the particle bounding box is not fully
	  // outside the porosity domain
	  if ( intersect( BBdomain, BBpart ) ) 
	  {
	    BBintVol.closest( BBdomain, BBpart );
	    intVolptA.setValue( BBintVol.getLower( X ), BBintVol.getLower( Y ),
	    	BBintVol.getLower( Z ) );
	    intVolExtent = BBintVol.getExtent();
	    nx = size_t( 2. * intVolExtent[X] / dx ) + 1;
	    ny = size_t( 2. * intVolExtent[Y] / dy ) + 1;
	    nz = size_t( 2. * intVolExtent[Z] / dz ) + 1;
	    dxl = 2. * intVolExtent[X] / double(nx);
	    dyl = 2. * intVolExtent[Y] / double(ny); 
	    dzl = 2. * intVolExtent[Z] / double(nz);
	    dvl = dxl * dyl * dzl;

            for (size_t i=0;i<nx;++i)
              for (size_t j=0;j<ny;++j)    
                for (size_t k=0;k<nz;++k)
	        {
	          // Coordinates of the center of the elementary volume
	          elemVolCenter[X] = intVolptA[X] + ( double(i) + 0.5 ) * dxl;
	          elemVolCenter[Y] = intVolptA[Y] + ( double(j) + 0.5 ) * dyl;
	          elemVolCenter[Z] = intVolptA[Z] + ( double(k) + 0.5 ) * dzl;
	  
	          // List of cells containing the cell that elemVolCenter 
		  // belongs to and its neighboring cells
	          cells = m_collision->getCellAndCellNeighborhood( 
		  	elemVolCenter );
	  
	          // Loops over the cells and check whether elemVolCenter 
		  // belongs to this particle in any of these cells
	          found = false;
	          for (ic=cells.begin();ic!=cells.end() && !found;ic++)
	            found = (*ic)->isInParticle( elemVolCenter, *il );
	    
                  if ( found ) volparticles += dvl;	  
	        }
	     	
	  }
	}
	// The particle is fully contained in the porosity domain
	// We simply add the volume
	else
	  volparticles += (*il)->getVolume();
      }
    }
    else if ( m_global_porosity->domain.ftype == WINDOW_CYLINDER )
    {
      volporodomain = PI * m_global_porosity->domain.radius
      	* m_global_porosity->domain.radius * m_global_porosity->domain.height;
      dx = ( m_global_porosity->domain.axisdir == X ?
      	m_global_porosity->domain.height 
		: 2. * m_global_porosity->domain.radius )
	/ double(m_global_porosity->nintervals[X]);
      dy = ( m_global_porosity->domain.axisdir == Y ?
      	m_global_porosity->domain.height 
		: 2. * m_global_porosity->domain.radius )
	/ double(m_global_porosity->nintervals[Y]);
      dz = ( m_global_porosity->domain.axisdir == Z ?
      	m_global_porosity->domain.height 
		: 2. * m_global_porosity->domain.radius )
	/ double(m_global_porosity->nintervals[Z]); 
      Point3 center;
      double circumscribed_radius = 0.;
      size_t intersect = 0;
      BBox BBdomain, BBintVol;
      Point3 BBCenter = m_global_porosity->domain.ptA;
      BBCenter[m_global_porosity->domain.axisdir] += 
      	0.5 * m_global_porosity->domain.height;
      Vector3 BBextent( m_global_porosity->domain.axisdir == X ?
      	0.5 * m_global_porosity->domain.height 
		: m_global_porosity->domain.radius,
	m_global_porosity->domain.axisdir == Y ?
      	0.5 * m_global_porosity->domain.height 
		: m_global_porosity->domain.radius,
	m_global_porosity->domain.axisdir == Z ?
      	0.5 * m_global_porosity->domain.height 
		: m_global_porosity->domain.radius );
      BBdomain.setCenter( BBCenter );
      BBdomain.setExtent( BBextent );					
	            
      for (il=allpart->begin();il!=allpart->end();il++)
      {
        center = *((*il)->getPosition());
	circumscribed_radius = (*il)->getCircumscribedRadius();
	intersect = GrainsExec::AACylinderSphereIntersection( center,
    		circumscribed_radius,
		m_global_porosity->domain.ptA,
		m_global_porosity->domain.radius,
		m_global_porosity->domain.height,
		m_global_porosity->domain.axisdir );

        // If the circumscribed sphere of the particle is fully contained
	// in the porosity domain, we simply add the volume
	if ( intersect == 0 )
	  volparticles += (*il)->getVolume();

        // If the circumscribed sphere of the particle intersects the cylinder
	// boundary, we pixelate both the part that belongs to the porosity
	// domain and the cylindrical porosity domain
	else if ( intersect == 2 )
	{
	  BBox BBpart = (*il)->BoundingBox(); 
	  BBintVol.closest( BBdomain, BBpart );
	  intVolptA.setValue( BBintVol.getLower( X ), BBintVol.getLower( Y ),
	    	BBintVol.getLower( Z ) );
	  intVolExtent = BBintVol.getExtent();
	  nx = size_t( 2. * intVolExtent[X] / dx ) + 1;
	  ny = size_t( 2. * intVolExtent[Y] / dy ) + 1;
	  nz = size_t( 2. * intVolExtent[Z] / dz ) + 1;
	  dxl = 2. * intVolExtent[X] / double(nx);
	  dyl = 2. * intVolExtent[Y] / double(ny); 
	  dzl = 2. * intVolExtent[Z] / double(nz);
	  dvl = dxl * dyl * dzl;

          for (size_t i=0;i<nx;++i)
            for (size_t j=0;j<ny;++j)    
              for (size_t k=0;k<nz;++k)
              {
	        // Coordinates of the center of the elementary volume
	        elemVolCenter[X] = intVolptA[X] + ( double(i) + 0.5 ) * dxl;
	        elemVolCenter[Y] = intVolptA[Y] + ( double(j) + 0.5 ) * dyl;
	        elemVolCenter[Z] = intVolptA[Z] + ( double(k) + 0.5 ) * dzl;
		
		// If the cell center belongs to the cylindrical porosity domain
		if ( GrainsExec::isPointInAACylinder( elemVolCenter,
			m_global_porosity->domain.ptA,
			m_global_porosity->domain.radius,
			m_global_porosity->domain.height,
			m_global_porosity->domain.axisdir ) )		 
	        {
	          // List of cells containing the cell that elemVolCenter 
		  // belongs to and its neighboring cells
	          cells = m_collision->getCellAndCellNeighborhood( 
		  	elemVolCenter );
	  
	          // Loops over the cells and check whether elemVolCenter 
		  // belongs to this particle in any of these cells
	          found = false;
	          for (ic=cells.begin();ic!=cells.end() && !found;ic++)
	            found = (*ic)->isInParticle( elemVolCenter, *il );
	    
                  if ( found ) volparticles += dvl;	  
	        }
	      }	
	}	
      }
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

  bool b2024 = false;
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

      // Restart mode is new, not read in XML file
      GrainsExec::m_ReloadType = "new" ;
      restart  = ReaderXML::getNodeAttr_String( reload, "Filename" );	
      restart = GrainsExec::fullResultFileName( restart, false );
      
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
   

    // Check that construction is fine
    if ( !m_restart ) 
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
    if ( b2024 )
      m_allcomponents.read_particles( restart, npart, m_collision, m_rank, 
      	m_nprocs, m_wrapper );
    
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
