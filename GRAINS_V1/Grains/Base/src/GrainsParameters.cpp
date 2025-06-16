#include "MPINeighbors.hh"
#include "GrainsParameters.hh"
#include "ContactBuilderFactory.hh"
#include "ObstacleBuilderFactory.hh"
#include "ObstacleImposedVelocity.hh"
#include "GrainsBuilderFactory.hh"
#include "PostProcessingWriter.hh"
#include "ParaviewPostProcessingWriter.hh"
#include "ContactForceModel.hh"
#include "TrilobeCylinder.hh"
#include "QuadrilobeCylinder.hh"


// ----------------------------------------------------------------------------
// Default constructor
GrainsParameters::GrainsParameters() 
  : Grains()
  , m_colRelVel( 0. ) 
{}




// ----------------------------------------------------------------------------
// Destructor
GrainsParameters::~GrainsParameters()
{}




// ----------------------------------------------------------------------------
// Writes an initial message in the standard output only on the process ranked 0
void GrainsParameters::initialOutputMessage()
{
  int rankproc = 0;
  MPI_Comm_rank( MPI_COMM_WORLD, &rankproc );
  if ( !rankproc )
    cout << "Grains3D Parameters" << endl;
}



// ----------------------------------------------------------------------------
// Tasks to perform before time-stepping 
void GrainsParameters::do_before_time_stepping( DOMElement* rootElement )
{
  Construction( rootElement );
  Forces( rootElement );
  AdditionalFeatures( rootElement );                    
}




// ----------------------------------------------------------------------------
// Tasks to perform before time-stepping 
void GrainsParameters::do_after_time_stepping()
{}




// ----------------------------------------------------------------------------
// Construction of the simulation: linked cell, particles & obstacles, domain 
// decomposition
void GrainsParameters::Construction( DOMElement* rootElement )
{
  assert( rootElement != NULL );
  DOMNode* root = ReaderXML::getNode( rootElement, "Construction" );

  bool bnewpart = false, bnewobst = false;

  // Domain size: origin and max coordinates
  DOMNode* domain = ReaderXML::getNode( root, "LinkedCell" );
  double lx = ReaderXML::getNodeAttr_Double( domain, "MX" );
  double ly = ReaderXML::getNodeAttr_Double( domain, "MY" );
  double lz = ReaderXML::getNodeAttr_Double( domain, "MZ" );
  DOMNode* domain_origin = ReaderXML::getNode( root, "Origin" );
  double ox = 0., oy = 0., oz = 0. ;
  if ( domain_origin )
  {
    ox = ReaderXML::getNodeAttr_Double( domain_origin, "OX" );
    oy = ReaderXML::getNodeAttr_Double( domain_origin, "OY" );
    oz = ReaderXML::getNodeAttr_Double( domain_origin, "OZ" ); 
  }   
  App::set_dimensions( lx, ly, lz, ox, oy, oz );

 
  // Domain decomposition
  readDomainDecomposition( root, lx - ox, ly - oy, lz - oz ); 
  
  
  // Construction on active processes
  if ( m_processorIsActive )
  {  
    // Output message
    if ( m_rank == 0 ) cout << GrainsExec::m_shift3 << "Construction" << endl;
    
    
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
	if ( sshape == "TrilobeCylinder" )
	  particleRef = new TrilobeCylinder( nCompParticle, nbPC+int(i) );
	else if ( sshape == "QuadrilobeCylinder" )
	  particleRef = new QuadrilobeCylinder( nCompParticle, nbPC+int(i) );
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
    if ( !bnewpart && !bnewobst ) 
    {
      if ( m_rank == 0 )
        cout << "ERR : Error in input file in <Contruction>" << endl;
      grainsAbort();
    }    
  }    
}




// ----------------------------------------------------------------------------
// External force definition
void GrainsParameters::Forces( DOMElement* rootElement )
{}




// ----------------------------------------------------------------------------
// Additional features of the simulation: time features, insertion, 
// post-processing
void GrainsParameters::AdditionalFeatures( DOMElement* rootElement )
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

    // Velocity relative
    DOMNode* nodeVrel = ReaderXML::getNode( root, 
    	"CollisionalRelativeVelocity" );
    if ( !nodeVrel )
    {
      cout << GrainsExec::m_shift6 << "<CollisionalRelativeVelocity> node "
      		<< "is mandatory !!" << endl;
      grainsAbort();
    }
    else
      m_colRelVel = ReaderXML::getNodeAttr_Double( nodeVrel, "value" );
  }
}




// ----------------------------------------------------------------------------
// Runs the simulation over the prescribed time interval
void GrainsParameters::Simulation( double time_interval )
{
  if ( m_processorIsActive )
  {
    vector<Particle*> const* particleClasses = 
    	m_allcomponents.getReferenceParticles();
    size_t nClasses = particleClasses->size(), i, j;

    cout << endl << "PARTICLE/PARTICLE CONTACTS" << endl;    
    // Between particles of same class
    for (i=0;i<nClasses;++i)
    {
      Particle* particle = new Particle( *((*particleClasses)[i]), true );    
      cout << "Contact Class " << i << " / Class " << i << endl;
      ContactBuilderFactory::contactForceModel( 
       	(*particleClasses)[i]->getMaterial(),
      	particle->getMaterial() )->computeAndWriteEstimates( 
		(*particleClasses)[i], particle, m_colRelVel, m_dt, cout );    
      delete particle;
    }    

    // Between particles of different class
    for (i=0;i<nClasses;++i)
      for (j=i+1;j<nClasses;++j)
      {
        cout << "Contact Class " << i << " / Class " << j << endl;	
	ContactBuilderFactory::contactForceModel(
		(*particleClasses)[i]->getMaterial(),
      		(*particleClasses)[j]->getMaterial() )
		->computeAndWriteEstimates( (*particleClasses)[i],
			 (*particleClasses)[j], m_colRelVel, m_dt, cout ); 	
      }

    cout << "PARTICLE/OBSTACLE CONTACTS" << endl;    
    // Particles-obstacles
    list<SimpleObstacle*> allObs = m_allcomponents.getObstacles()
    	->getObstacles();
    list<SimpleObstacle*>::iterator il;    
    for (i=0;i<nClasses;++i)
      for (il=allObs.begin();il!=allObs.end();il++)
      {
        cout << "Contact Class " << i << " / " << (*il)->getName() 
		<< endl;
	ContactBuilderFactory::contactForceModel(
	  	(*particleClasses)[i]->getMaterial(), (*il)->getMaterial() )
		->computeAndWriteEstimates( (*particleClasses)[i], *il,
		m_colRelVel, m_dt, cout ); 
      } 
	
    cout << "PARTICLE DISPLACEMENT" << endl;    
    for (i=0;i<nClasses;++i)
    {
      cout << "  Class " << i << " : v0*dt/crust = " << 
      	m_colRelVel * m_dt / (*particleClasses)[i]->getCrustThickness() 
	<< endl;
    }   	     
  }
}
