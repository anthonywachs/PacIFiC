#include "GrainsMPIWrapper.hh"
#include "GrainsExec.hh"
#include "GrainsBuilderFactory.hh"
#include "ParaviewPostProcessingWriter.hh"
#include "Particle.hh"
#include "Component.hh"
#include "Obstacle.hh"
#include "Box.hh"
#include "Cylinder.hh"
#include "Segment.hh"
#include "Cell.hh"
#include "Vector3.hh"
#include <zlib.h>
using namespace solid;

static int sizeof_Float32 = 4 ;
static int sizeof_Int32 = 4 ;


// ----------------------------------------------------------------------------
// Default constructor
ParaviewPostProcessingWriter::ParaviewPostProcessingWriter()
{}




// ----------------------------------------------------------------------------
// Constructor with XML node, rank and number of processes as input parameters
ParaviewPostProcessingWriter::ParaviewPostProcessingWriter( DOMNode* dn,
    int const& rank_, int const& nbranks_, bool const& verbose )
  : PostProcessingWriter( dn, rank_, nbranks_ )
  , m_ParaviewCycleNumber( 0 )
  , m_binary( false )
  , m_postProcessObstacle( true )
  , m_initialCycleNumber_forced( false )
  , m_network( false )
  , m_mpiio_singlefile( false )
  , m_pertype( true )
  , BUFFER( NULL )
  , ALLOCATED( 0 )
  , OFFSET( 0 )
{
  m_ParaviewFilename = ReaderXML::getNodeAttr_String( dn, "RootName" );
  m_ParaviewFilename_dir = ReaderXML::getNodeAttr_String( dn, "Directory" );
  if ( ReaderXML::hasNodeAttr( dn, "InitialCycleNumber" ) )
    m_ParaviewCycleNumber = ReaderXML::getNodeAttr_Int( dn, 
  	"InitialCycleNumber" );
  string sm_binary = ReaderXML::getNodeAttr_String( dn, "WritingMode" );
  if ( sm_binary == "Binary" ) m_binary = true;
  if ( ReaderXML::hasNodeAttr( dn, "ForceNetwork" ) ) 
  {
    string sm_network = ReaderXML::getNodeAttr_String( dn, "ForceNetwork" );
    if ( sm_network == "True" ) m_network = true;
  }
  if ( ReaderXML::hasNodeAttr( dn, "Obstacle" ) )
  { 
    string sm_obstacle = ReaderXML::getNodeAttr_String( dn, "Obstacle" );
    if ( sm_obstacle == "False" ) m_postProcessObstacle = false; 
  }
  if ( ReaderXML::hasNodeAttr( dn, "MPIIO" ) && m_nprocs > 1 )
  { 
    string sm_onefile = ReaderXML::getNodeAttr_String( dn, "MPIIO" );
    if ( sm_onefile == "True" ) m_mpiio_singlefile = true; 
  } 
  if ( ReaderXML::hasNodeAttr( dn, "PerType" ) )
  { 
    string sm_pertype = ReaderXML::getNodeAttr_String( dn, "PerType" );
    if ( sm_pertype == "False" ) m_pertype = false; 
  }    
  
  if ( m_rank == 0 && verbose )
  {
    cout << GrainsExec::m_shift9 << "Type = Paraview" << endl;
    cout << GrainsExec::m_shift12 << "Output file root name = " 
    	<< m_ParaviewFilename << endl;
    cout << GrainsExec::m_shift12 << "Output file directory name = " 
    	<< m_ParaviewFilename_dir << endl;
    cout << GrainsExec::m_shift12 << "Writing mode = " 
    	<< ( m_binary ? "Binary" : "Text" ) << endl;
    cout << GrainsExec::m_shift12 << "Write force network = " 
    	<< ( m_network ? "True" : "False" ) << endl;
    cout << GrainsExec::m_shift12 << "Write obstacles = " 
    	<< ( m_postProcessObstacle ? "True" : "False" ) << endl;
    cout << GrainsExec::m_shift12 << "MPIIO = " 
    	<< ( m_mpiio_singlefile ? "True" : "False" ) << endl;
    cout << GrainsExec::m_shift12 << "Per particle type = " 
    	<< ( m_pertype ? "True" : "False" ) << endl;			
  } 
}




// ----------------------------------------------------------------------------
// Constructor with input parameters 
ParaviewPostProcessingWriter::ParaviewPostProcessingWriter(
    int const& rank_,
    int const& nbranks_,
    const string &name_,
    const string &root_,
    const bool &isBinary,
    bool const& verbose )
  : PostProcessingWriter( rank_, nbranks_ )
  , m_ParaviewFilename_dir( root_ )
  , m_ParaviewFilename( name_ )
  , m_ParaviewCycleNumber( 0 )
  , m_binary( isBinary )
  , m_postProcessObstacle( true ) 
  , m_initialCycleNumber_forced( false )
  , m_network( false )
  , m_mpiio_singlefile( false ) 
  , m_pertype( true )      
  , BUFFER( NULL )
  , ALLOCATED( 0 )
  , OFFSET( 0 )
{
  if ( m_rank == 0 && verbose )
  {
    cout << GrainsExec::m_shift9 << "Type = Paraview" << endl;
    cout << GrainsExec::m_shift12 << "Output file root name = " 
    	<< m_ParaviewFilename << endl;
    cout << GrainsExec::m_shift12 << "Output file directory name = " 
    	<< m_ParaviewFilename_dir << endl;
    cout << GrainsExec::m_shift12 << "Writing mode = " 
    	<< ( m_binary ? "Binary" : "Text" ) << endl;
    cout << GrainsExec::m_shift12 << "Write force network = " 
    	<< ( m_network ? "True" : "False" ) << endl;
    cout << GrainsExec::m_shift12 << "Write obstacles = " 
    	<< ( m_postProcessObstacle ? "True" : "False" ) << endl;
    cout << GrainsExec::m_shift12 << "MPIIO = " 
    	<< ( m_mpiio_singlefile ? "True" : "False" ) << endl;
    cout << GrainsExec::m_shift12 << "Per particle type = " 
    	<< ( m_pertype ? "True" : "False" ) << endl;				
  }
}	




// ----------------------------------------------------------------------------
// Destructor
ParaviewPostProcessingWriter::~ParaviewPostProcessingWriter()
{
  vector<ostringstream*>::iterator iv;
  for (iv=m_Paraview_saveParticles_pvd.begin(); 
  	iv!=m_Paraview_saveParticles_pvd.end();iv++) delete *iv;
  m_Paraview_saveParticles_pvd.clear(); 
}




// ----------------------------------------------------------------------------
// Initializes the post-processing writer
void ParaviewPostProcessingWriter::PostProcessing_start(
    double const& time, 
    double const& dt,
    list<Particle*> const* particles,
    list<Particle*> const* pwait,
    list<Particle*> const* periodic_clones,
    vector<Particle*> const* referenceParticles,
    Obstacle *obstacle,
    LinkedCell const* LC,
    vector<Window> const& insert_windows )
{
  size_t nbParticleTypes = referenceParticles->size() ;
  size_t nbPPTypes = m_pertype ? nbParticleTypes : 1;  

  if ( GrainsExec::m_ReloadType == "new" ) 
  {
    clearResultFiles();
    if ( GrainsExec::m_MPI )
      GrainsExec::getComm()->MPI_Barrier_ActivProc();
    
    if ( m_rank == 0 )
    {
      // Obstacles
      if ( m_postProcessObstacle )
      {
        m_Paraview_saveObstacles_pvd << "<?xml version=\"1.0\"?>" << endl;
        m_Paraview_saveObstacles_pvd << 
       		"<VTKFile type=\"Collection\" version=\"0.1\""
       		<< " byte_order=\"LittleEndian\"";
        if ( m_binary ) 
          m_Paraview_saveObstacles_pvd << 
	  	" compressor=\"vtkZLibDataCompressor\"";
        m_Paraview_saveObstacles_pvd << ">" << endl;
        m_Paraview_saveObstacles_pvd << "<Collection>" << endl;
      }
       
      // Particles
      ostringstream *ossNULL = NULL;
      m_Paraview_saveParticles_pvd.reserve( nbPPTypes );
      for (size_t i=0;i<nbPPTypes;++i)
        m_Paraview_saveParticles_pvd.push_back(ossNULL);
      for (size_t i=0;i<nbPPTypes;++i)
        m_Paraview_saveParticles_pvd[i] = new ostringstream;      
    
      for (size_t i=0;i<nbPPTypes;++i)
      {
        *m_Paraview_saveParticles_pvd[i] << "<?xml version=\"1.0\"?>" << endl;
        *m_Paraview_saveParticles_pvd[i] << 
       		"<VTKFile type=\"Collection\" version=\"0.1\""
       		<< " byte_order=\"LittleEndian\"";
        
	if ( m_binary ) 
          *m_Paraview_saveParticles_pvd[i] 
		<< " compressor=\"vtkZLibDataCompressor\"";
			
        *m_Paraview_saveParticles_pvd[i] << ">" << endl;
        *m_Paraview_saveParticles_pvd[i] << "<Collection>" << endl;
      }
    
      // Periodic clone particles
      if ( GrainsExec::m_periodic )
      { 
        m_Paraview_savePeriodicCloneParticles_pvd << "<?xml version=\"1.0\"?>" 
		<< endl;
        m_Paraview_savePeriodicCloneParticles_pvd << 
       		"<VTKFile type=\"Collection\" version=\"0.1\""
       		<< " byte_order=\"LittleEndian\"";
        
	if ( m_binary ) 
          m_Paraview_savePeriodicCloneParticles_pvd 
		<< " compressor=\"vtkZLibDataCompressor\"";
			
        m_Paraview_savePeriodicCloneParticles_pvd << ">" << endl;	
        m_Paraview_savePeriodicCloneParticles_pvd << "<Collection>" << endl; 
      }             

      if ( LC )
      {
        // Particle translational and angular velocity vectors
        m_Paraview_saveParticleVelocityVectors_pvd << "<?xml version=\"1.0\"?>" 
		<< endl;
        m_Paraview_saveParticleVelocityVectors_pvd << 
       		"<VTKFile type=\"Collection\" version=\"0.1\""
       		<< " byte_order=\"LittleEndian\"";
	      
        if ( m_binary ) 
          m_Paraview_saveParticleVelocityVectors_pvd 
		<< " compressor=\"vtkZLibDataCompressor\"";
			
        m_Paraview_saveParticleVelocityVectors_pvd << ">" << endl; 
        m_Paraview_saveParticleVelocityVectors_pvd << "<Collection>" << endl; 
    
        // Contact force vectors
        m_Paraview_saveContactForceVectors_pvd << "<?xml version=\"1.0\"?>" 
		<< endl;
        m_Paraview_saveContactForceVectors_pvd << 
       		"<VTKFile type=\"Collection\" version=\"0.1\""
       		<< " byte_order=\"LittleEndian\"";
	      
        if ( m_binary ) 
          m_Paraview_saveContactForceVectors_pvd 
		<< " compressor=\"vtkZLibDataCompressor\"";
		
        m_Paraview_saveContactForceVectors_pvd << ">" << endl; 
        m_Paraview_saveContactForceVectors_pvd << "<Collection>" << endl; 
    
	
	// Contact force network
	if ( m_network )
	{
	  m_Paraview_saveContactForceChains_pvd << "<?xml version=\"1.0\"?>" 
	  	<< endl;
	  m_Paraview_saveContactForceChains_pvd << 
       		"<VTKFile type=\"Collection\" version=\"0.1\""
       		<< " byte_order=\"LittleEndian\"";
	      
	  if ( m_binary ) 
	    m_Paraview_saveContactForceChains_pvd 
		<< " compressor=\"vtkZLibDataCompressor\"";
		
	  m_Paraview_saveContactForceChains_pvd<< ">" << endl; 
	  m_Paraview_saveContactForceChains_pvd << "<Collection>" << endl; 
	}
	
        // Linked cell grid 
	if ( m_nprocs > 1 && !m_mpiio_singlefile )
          writePVTU_Paraview( m_ParaviewFilename + "_LinkedCell" , 
		&empty_string_list, &empty_string_list, &empty_string_list );
      }                  
    } 

    // Linked cell grid
    if ( LC )
    { 
      ostringstream ossRK;
      ossRK << m_rank;       
      writeLinkedCellPostProcessing_Paraview( LC,
      	m_ParaviewFilename + "_LinkedCell" + 
	( m_nprocs > 1 && !m_mpiio_singlefile ? "_" + ossRK.str() : "" ) 
	+ ".vtu",  m_nprocs > 1 && !m_mpiio_singlefile );
    } 
           
    one_output( time, dt, particles, periodic_clones, referenceParticles,
  	obstacle, LC );
  }
  else
  {
    // In reload mode "same", the obstacle post-processing file is required
    // therefore we force m_postProcessObstacle = true
    m_postProcessObstacle = true ;

    if ( m_rank == 0 )
    {     
      // Obstacles
      readPVDFile( m_ParaviewFilename_dir + "/" + m_ParaviewFilename
       	+ "_Obstacles.pvd", m_Paraview_saveObstacles_pvd );

      // Particles
      ostringstream *ossNULL = NULL;
      m_Paraview_saveParticles_pvd.reserve( nbPPTypes );
      for (size_t i=0;i<nbPPTypes;++i)
        m_Paraview_saveParticles_pvd.push_back(ossNULL);
      for (size_t i=0;i<nbPPTypes;++i)
        m_Paraview_saveParticles_pvd[i] = new ostringstream;      
    
      if ( nbPPTypes == 1 )
        readPVDFile( m_ParaviewFilename_dir + "/" + m_ParaviewFilename
       	+ "_Particles.pvd", *m_Paraview_saveParticles_pvd[0] );
      else
        for (size_t i=0;i<nbPPTypes;++i)
	{
          ostringstream* ossPC = new ostringstream;
          *ossPC << i;
          readPVDFile( m_ParaviewFilename_dir + "/" + m_ParaviewFilename
		+ "_Particles_Type"+ossPC->str() + ".pvd",
		*m_Paraview_saveParticles_pvd[i] );
	  delete ossPC;	
	}     
 
      if ( GrainsExec::m_periodic )
        readPVDFile( m_ParaviewFilename_dir + "/" + m_ParaviewFilename
       	+ "_PeriodicCloneParticles.pvd", 
	m_Paraview_savePeriodicCloneParticles_pvd );

      if ( LC )
      {
        // Particle translational and angular velocity vectors
	readPVDFile( m_ParaviewFilename_dir + "/" + m_ParaviewFilename
       	+ "_ParticleVelocityVectors.pvd", 
	m_Paraview_saveParticleVelocityVectors_pvd );  

        // Contact force vectors
	readPVDFile( m_ParaviewFilename_dir + "/" + m_ParaviewFilename
       	+ "_ContactForceVectors.pvd", m_Paraview_saveContactForceVectors_pvd );

        // Contact force network
	if ( m_network )
        {
          readPVDFile( m_ParaviewFilename_dir + "/" + m_ParaviewFilename
       		+ "_ContactForceChains.pvd", 
		m_Paraview_saveContactForceChains_pvd );  
        }
      }
    }  

    // Cycle number
    if ( !m_initialCycleNumber_forced )
      m_ParaviewCycleNumber = getPreviousCycleNumber();
    ++m_ParaviewCycleNumber;                
  }
  
  // Windows d'insertion
  if ( m_rank == 0 ) 
  {
    writeInsertionPostProcessing_Paraview( insert_windows, 
    	m_ParaviewFilename + "_InsertionWindows" );
	
    if ( GrainsExec::m_periodic == true && LC )
      writePeriodicBoundaryPostProcessing_Paraview(
  	LC, m_ParaviewFilename + "_PeriodicBoundaries" );
  }	
}




// ----------------------------------------------------------------------------
// Reads a pvd file and transfers its content to an output stream
void ParaviewPostProcessingWriter::readPVDFile( string const& filename, 
	ostringstream& ossflux )
{
  string tline, buffer, part, keyword;
  bool keep = true ;
  int cyleNumber = 0;
    
  ifstream fileIN( filename.c_str(), ios::in );
  getline( fileIN, tline );
  while ( !fileIN.eof() ) 
  { 
    if ( tline != "</Collection>" && tline != "</VTKFile>" )
    {
      keep = true ;
      if ( m_initialCycleNumber_forced )
      {
        istringstream iss( tline );
        iss >> keyword;
        if ( keyword == "<DataSet" )
        {
          iss >> buffer >> buffer >> buffer >> part;
          size_t pos = part.find( "_T" );
          string sub = part.substr( pos );
          sub.erase( sub.begin(), sub.begin() + 2 );
	  pos = sub.rfind( "." );
          sub.erase( sub.begin() + pos, sub.end() );	  
          istringstream issNum( sub );
          issNum >> cyleNumber; 
          if ( cyleNumber > m_ParaviewCycleNumber ) keep = false ;
        }
      }       
    }
    else keep = false ;
    if ( keep ) ossflux << tline << endl;
    getline( fileIN, tline );  
  }
  fileIN.close();  
}




// ---------------------------------------------------------------------------
// Gets the last cycle number in case of reload
int ParaviewPostProcessingWriter::getPreviousCycleNumber() const
{
  int cyleNumber = 0;
  string tline, previous_tline, buffer, part;

  // Read the obstacles pvd file
  ifstream fileIN( (m_ParaviewFilename_dir + "/" + m_ParaviewFilename
       	+ "_Obstacles.pvd" ).c_str(), ios::in );
  while ( tline != "</Collection>" )
  {
    previous_tline = tline;
    getline( fileIN, tline );
  }
  
  // Process the penultimate line to extract the cycle number
  istringstream iss( previous_tline ); 
  iss >> buffer >> buffer >> buffer >> buffer >> part;
  size_t pos = part.find( "_T" );
  string sub = part.substr( pos );
  sub.erase( sub.begin(), sub.begin() + 2 );
  sub.erase( sub.end() - 7, sub.end() ); 
  istringstream issNum( sub );
  issNum >> cyleNumber; 

  return cyleNumber;
}




// ----------------------------------------------------------------------------
// Writes data
void ParaviewPostProcessingWriter::PostProcessing(
    double const& time, 
    double const& dt,
    list<Particle*> const* particles,
    list<Particle*> const* pwait,
    list<Particle*> const* periodic_clones,
    vector<Particle*> const* referenceParticles,
    Obstacle* obstacle,
    LinkedCell const* LC )
{
  one_output( time, dt, particles, periodic_clones, referenceParticles,
  	obstacle, LC );
}




// ----------------------------------------------------------------------------
// Finalizes writing data
void ParaviewPostProcessingWriter::PostProcessing_end()
{}




// ----------------------------------------------------------------------------
// Writes data at one physical time
void ParaviewPostProcessingWriter::one_output(
    double const& time,
    double const& dt,
    list<Particle*> const* particles,
    list<Particle*> const* periodic_clones,
    vector<Particle*> const* referenceParticles,
    Obstacle* obstacle,
    LinkedCell const* LC )
{
  size_t nbParticleTypes = referenceParticles->size() ;
  size_t nbPPTypes = m_pertype ? nbParticleTypes : 1;  
  list<string> Scalars;
  Scalars.push_back("NormU");
  Scalars.push_back("NormOm");
  Scalars.push_back("CoordNumb");
  ostringstream ossCN, ossRK;
  ossCN << m_ParaviewCycleNumber; 
  ossRK << m_rank;   

  // Obstacles
  if ( m_postProcessObstacle )
  {
    updateObstaclesIndicator( time, dt, obstacle );

    if ( m_rank == 0 ) 
    {  
      list<SimpleObstacle*> allObstacles = obstacle->getObstacles();
     
      string obsFilename = m_ParaviewFilename + "_Obstacles_T" + ossCN.str() + 
   	".vtu";
      m_Paraview_saveObstacles_pvd << "<DataSet timestep=\"" << time 
      	<< "\" " << "group=\"\" part=\"0\" file=\"" << obsFilename << "\"/>\n"; 
	              
      ofstream f( ( m_ParaviewFilename_dir + "/" + m_ParaviewFilename
       	+ "_Obstacles.pvd" ).c_str(), ios::out );	 
      f << m_Paraview_saveObstacles_pvd.str();
      f << "</Collection>" << endl;
      f << "</VTKFile>" << endl;
      f.close();      
            
      if ( obstacle ) writeObstaclesPostProcessing_Paraview( allObstacles,
     		obsFilename );   
    }
  }
   
  // Particles
  if ( nbPPTypes == 1 )
  { 
    string partFilename = m_ParaviewFilename + "_Particles_T" + ossCN.str();
    if ( m_rank == 0 ) 
    {     
      *m_Paraview_saveParticles_pvd[0] << "<DataSet timestep=\"" << time 
      	<< "\" " << "group=\"\" part=\"0\" file=\"" << partFilename 
	<< ( m_nprocs > 1 && !m_mpiio_singlefile ? 
	  ".pvtu\"/>" : ".vtu\"/>" ) << endl;             
       
      ofstream g( ( m_ParaviewFilename_dir + "/" + m_ParaviewFilename
       	+ "_Particles.pvd" ).c_str(), ios::out );
      g << m_Paraview_saveParticles_pvd[0]->str();
      g << "</Collection>" << endl;
      g << "</VTKFile>" << endl;
      g.close();

      if ( m_nprocs > 1 && !m_mpiio_singlefile )
      {       
        if ( nbParticleTypes == 1 &&
	(*referenceParticles)[0]->getRigidBody()->getConvex()
       	->isSphere() && !GrainsExec::m_SphereAsPolyParaview )
        {
          list<string> ptVec;
          ptVec.push_back("Orientation");
          writePVTU_Paraview( partFilename, &ptVec, &Scalars, 
	 	&empty_string_list );
        }
        else writePVTU_Paraview( partFilename, &empty_string_list,
	 	&empty_string_list, &Scalars ); 
      }    
    }
    
    // VTU files
    if ( nbParticleTypes == 1 &&
    	(*referenceParticles)[0]->getRigidBody()->getConvex()
	->isSphere() && !GrainsExec::m_SphereAsPolyParaview )
    {
      if ( m_mpiio_singlefile && m_nprocs > 1 )
      {
	if ( m_binary )
	  writeSpheresPostProcessing_Paraview_MPIIO_binary( particles,
     	  	partFilename + ".vtu", false, 
		PostProcessingWriter::m_bPPWindow[m_rank] );
	else
	  writeSpheresPostProcessing_Paraview_MPIIO_text( particles,
     	  	partFilename + ".vtu", false,
		PostProcessingWriter::m_bPPWindow[m_rank] );	    	     
      }
      else		
	writeSpheresPostProcessing_Paraview( particles, partFilename 
	  	+ ( m_nprocs > 1 ? "_" + ossRK.str() : "" ) + ".vtu", false,
		PostProcessingWriter::m_bPPWindow[m_rank] );
    }
    else
    {             
      if ( m_mpiio_singlefile && m_nprocs > 1 )
      {
	if ( m_binary )
	  writeParticlesPostProcessing_Paraview_MPIIO_binary( particles,
     	  	partFilename + ".vtu", false, 
		PostProcessingWriter::m_bPPWindow[m_rank] );
	  else
	    writeParticlesPostProcessing_Paraview_MPIIO_text( particles,
     	  	partFilename + ".vtu", false,
		PostProcessingWriter::m_bPPWindow[m_rank] );	    	     
      }
      else		
	writeParticlesPostProcessing_Paraview( particles, partFilename 
	  + ( m_nprocs > 1 ? "_" + ossRK.str() : "" ) + ".vtu", false,
	  PostProcessingWriter::m_bPPWindow[m_rank] );
    }
  }
  else
  {
    list<Particle*> empty_list;
    vector< list<Particle*> > partPerType( nbParticleTypes, 
     	empty_list );
    list<Particle*>::const_iterator particle;
     
    for (particle=particles->begin();particle!=particles->end();
     	particle++)
      partPerType[(*particle)->getGeometricType()].push_back(*particle);
       
    for (size_t i=0;i<nbParticleTypes;++i)
    {
      ostringstream* ossPC = new ostringstream;
      *ossPC << i;
      string partFilename = m_ParaviewFilename + "_Particles_Type" +
      	ossPC->str() + "_T" + ossCN.str();

      if ( m_rank == 0 ) 
      { 	
        *m_Paraview_saveParticles_pvd[i] << "<DataSet timestep=\"" << time 
       		<< "\" " << "group=\"\" part=\"0\" file=\"" << partFilename 
		<< ( m_nprocs > 1 && !m_mpiio_singlefile ? 
	  	".pvtu\"/>" : ".vtu\"/>" ) << endl;  
       
        ofstream g( ( m_ParaviewFilename_dir + "/" + m_ParaviewFilename
       		+ "_Particles_Type" + ossPC->str() + ".pvd" ).c_str(), 
		ios::out );
        g << m_Paraview_saveParticles_pvd[i]->str();
        g << "</Collection>" << endl;
        g << "</VTKFile>" << endl;
        g.close(); 

        if ( m_nprocs > 1 && !m_mpiio_singlefile )
	{
	  if ( (*referenceParticles)[i]->getRigidBody()->getConvex()
	 	->isSphere() && !GrainsExec::m_SphereAsPolyParaview )
          {
            list<string> ptVec;
            ptVec.push_back("Orientation");
            writePVTU_Paraview( partFilename, &ptVec, &Scalars,
	   	&empty_string_list );
          }
          else writePVTU_Paraview( partFilename, &empty_string_list,
	 	&empty_string_list, &Scalars ); 
	}  
      }

      // VTU files
      if ( (*referenceParticles)[i]->getRigidBody()->getConvex()
          ->isSphere() && !GrainsExec::m_SphereAsPolyParaview )
      {
        if ( m_mpiio_singlefile && m_nprocs > 1 )
        {
	  if ( m_binary )
	    writeSpheresPostProcessing_Paraview_MPIIO_binary( &partPerType[i],
     	  	partFilename + ".vtu", false, 
		PostProcessingWriter::m_bPPWindow[m_rank] );
	  else
	    writeSpheresPostProcessing_Paraview_MPIIO_text( &partPerType[i],
     	  	partFilename + ".vtu", false,
		PostProcessingWriter::m_bPPWindow[m_rank] );	    	     
        }
        else		
	  writeSpheresPostProcessing_Paraview( &partPerType[i], partFilename 
	  	+ ( m_nprocs > 1 ? "_" + ossRK.str() : "" ) + ".vtu", false,
		PostProcessingWriter::m_bPPWindow[m_rank] );
      }
      else
      {             
        if ( m_mpiio_singlefile && m_nprocs > 1 )
        {
	  if ( m_binary )
	    writeParticlesPostProcessing_Paraview_MPIIO_binary( &partPerType[i],
     	  	partFilename + ".vtu", false, 
		PostProcessingWriter::m_bPPWindow[m_rank] );
	  else
	    writeParticlesPostProcessing_Paraview_MPIIO_text( &partPerType[i],
     	  	partFilename + ".vtu", false,
		PostProcessingWriter::m_bPPWindow[m_rank] );	    	     
        }
        else		
	  writeParticlesPostProcessing_Paraview( &partPerType[i], partFilename 
	  + ( m_nprocs > 1 ? "_" + ossRK.str() : "" ) + ".vtu", false,
	  PostProcessingWriter::m_bPPWindow[m_rank] );
      }
      
      delete ossPC;
    }
  }
   
  if ( GrainsExec::m_periodic ) 
  { 
    string partFilename = m_ParaviewFilename + "_PeriodicCloneParticles_T" + 
     	ossCN.str();
	
    if ( m_rank == 0 ) 
    {
      m_Paraview_savePeriodicCloneParticles_pvd << "<DataSet timestep=\"" 
      	<< time << "\" " << "group=\"\" part=\"0\" file=\"" << partFilename
	<< ( m_nprocs > 1 && !m_mpiio_singlefile ? ".pvtu\"/>" : ".vtu\"/>" ) 
	<< endl; 	
       
      ofstream g( ( m_ParaviewFilename_dir + "/" + m_ParaviewFilename
       	+ "_PeriodicCloneParticles.pvd" ).c_str(), ios::out );
      g << m_Paraview_savePeriodicCloneParticles_pvd.str();
      g << "</Collection>" << endl;
      g << "</VTKFile>" << endl;
      g.close(); 

      if ( m_nprocs > 1 && !m_mpiio_singlefile )
      {
	if ( nbParticleTypes == 1 
		&& (*referenceParticles)[0]->getRigidBody()->getConvex()
	 	->isSphere() && !GrainsExec::m_SphereAsPolyParaview )
        {
          list<string> ptVec;
          ptVec.push_back("Orientation");
          writePVTU_Paraview( partFilename, &ptVec, &Scalars,
	   	&empty_string_list );
        }
        else writePVTU_Paraview( partFilename, &empty_string_list,
	 	&empty_string_list, &Scalars ); 
      } 
    }

    // VTU files
    if ( nbParticleTypes == 1 
	&& (*referenceParticles)[0]->getRigidBody()->getConvex()
          ->isSphere() && !GrainsExec::m_SphereAsPolyParaview )
    {
      if ( m_mpiio_singlefile && m_nprocs > 1 )
      {
	if ( m_binary )
	  writeSpheresPostProcessing_Paraview_MPIIO_binary( periodic_clones,
     	  	partFilename + ".vtu", true, 
		PostProcessingWriter::m_bPPWindow[m_rank] );
	else
	  writeSpheresPostProcessing_Paraview_MPIIO_text( periodic_clones,
     	  	partFilename + ".vtu", true,
		PostProcessingWriter::m_bPPWindow[m_rank] );	    	     
      }
      else		
	writeSpheresPostProcessing_Paraview( periodic_clones, partFilename 
	  	+ ( m_nprocs > 1 ? "_" + ossRK.str() : "" ) + ".vtu", true,
		PostProcessingWriter::m_bPPWindow[m_rank] );
    }
    else
    {             
      if ( m_mpiio_singlefile && m_nprocs > 1 )
      {
	if ( m_binary )
	  writeParticlesPostProcessing_Paraview_MPIIO_binary( periodic_clones,
     	  	partFilename + ".vtu", true, 
		PostProcessingWriter::m_bPPWindow[m_rank] );
	  else
	    writeParticlesPostProcessing_Paraview_MPIIO_text( periodic_clones,
     	  	partFilename + ".vtu", true,
		PostProcessingWriter::m_bPPWindow[m_rank] );	    	     
      }
      else		
	writeParticlesPostProcessing_Paraview( periodic_clones, partFilename 
	  + ( m_nprocs > 1 ? "_" + ossRK.str() : "" ) + ".vtu", true,
	  PostProcessingWriter::m_bPPWindow[m_rank] );
    }
  }
   
  if ( LC )
  {
    // Particle translational and angular velocity vectors
    string vectFilename = m_ParaviewFilename + "_ParticleVelocityVectors_T" 
    	+ ossCN.str();
	
    if ( m_rank == 0 ) 
    {   
      m_Paraview_saveParticleVelocityVectors_pvd << "<DataSet timestep=\"" 
      	<< time << "\" " << "group=\"\" part=\"0\" file=\"" << vectFilename
	<< ( m_nprocs > 1 && !m_mpiio_singlefile ? ".pvtu\"/>" : ".vtu\"/>" ) 
	<< endl; 	
       
      ofstream h( ( m_ParaviewFilename_dir + "/" + m_ParaviewFilename
       	+ "_ParticleVelocityVectors.pvd" ).c_str(), ios::out );
      h << m_Paraview_saveParticleVelocityVectors_pvd.str();
      h << "</Collection>" << endl;
      h << "</VTKFile>" << endl;
      h.close();      

      list<string> vecMotion;
      vecMotion.push_back("U");
      vecMotion.push_back("Omega");
      if ( m_nprocs > 1 && !m_mpiio_singlefile )
        writePVTU_Paraview( vectFilename, &vecMotion, &empty_string_list,
       		&empty_string_list ); 
    }  
            
    // VTU files
    if ( m_mpiio_singlefile && m_nprocs > 1 )
    {
      if ( m_binary )
	writeParticleVelocityVectorsPostProcessing_Paraview_MPIIO_binary( 
		particles, vectFilename + ".vtu", false,
		PostProcessingWriter::m_bPPWindow[m_rank] );
	else
	  writeParticleVelocityVectorsPostProcessing_Paraview_MPIIO_text( 
	  	particles, vectFilename + ".vtu", false,
		PostProcessingWriter::m_bPPWindow[m_rank] );	    	     
    }
    else		
      writeParticleVelocityVectorsPostProcessing_Paraview( particles,
   	vectFilename + ( m_nprocs > 1 ? "_" + ossRK.str() : "" ) + ".vtu", 
	PostProcessingWriter::m_bPPWindow[m_rank] );


    // Contact force vectors
    string forceFilename = m_ParaviewFilename + "_ContactForceVectors_T" + 
      	ossCN.str();
	
    if ( m_rank == 0 ) 
    {   
      m_Paraview_saveContactForceVectors_pvd << "<DataSet timestep=\"" << time 
	<< "\" " << "group=\"\" part=\"0\" file=\"" << forceFilename
 	<< ( m_nprocs > 1 && !m_mpiio_singlefile ? ".pvtu\"/>" : ".vtu\"/>" ) 
	<< endl; 	
        
      ofstream h( ( m_ParaviewFilename_dir + "/" + m_ParaviewFilename
        	+ "_ContactForceVectors.pvd" ).c_str(), ios::out );
      h << m_Paraview_saveContactForceVectors_pvd.str();   
      h << "</Collection>" << endl;
      h << "</VTKFile>" << endl;
      h.close();      
 
      list<string> vecForce;
      vecForce.push_back("Force");
      if ( m_nprocs > 1 && !m_mpiio_singlefile )
        writePVTU_Paraview( forceFilename, &vecForce, &empty_string_list,
        	&empty_string_list ); 
    }

    // VTU files
    if ( m_mpiio_singlefile && m_nprocs > 1 )
    {
      if ( m_binary )
	writeContactForceVectorsPostProcessing_Paraview_MPIIO_binary( 
		particles, LC, forceFilename + ".vtu", time,
		PostProcessingWriter::m_bPPWindow[m_rank] );
	else
	writeContactForceVectorsPostProcessing_Paraview_MPIIO_text( 
		particles, LC, forceFilename + ".vtu", time,
		PostProcessingWriter::m_bPPWindow[m_rank] );    	     
    }
    else		
      writeContactForceVectorsPostProcessing_Paraview( particles, LC,
   	forceFilename + ( m_nprocs > 1 ? "_" + ossRK.str() : "" ) + ".vtu", 
	time, PostProcessingWriter::m_bPPWindow[m_rank] );


    if ( m_network )
    {
//       // Force chain
//       string forceChainFilename = m_ParaviewFilename + "_ContactForceChains_T" +
//         	ossCN.str();
// 		
//       if ( m_rank == 0 )
//       {
//         m_Paraview_saveContactForceChains_pvd << "<DataSet timestep=\"" << time 
//           	<< "\" " << "group=\"\" part=\"0\" file=\"" << 
// 		forceChainFilename << ".pvtp\"/>" << endl; 	
//           
//         ofstream h( ( m_ParaviewFilename_dir + "/" + m_ParaviewFilename
//           	+ "_ContactForceChains.pvd" ).c_str(), ios::out );
//         h << m_Paraview_saveContactForceChains_pvd.str();      
//         h << "</Collection>" << endl;
//         h << "</VTKFile>" << endl;
//         h.close();      
//  
//         list<string> vecForce;
//         vecForce.push_back("ForceChain");
//         writePVTP_Paraview( forceChainFilename, &vecForce, &empty_string_list,
//           	&empty_string_list ); 
//       }
//  
//       // Does this processor write data ?
//       if ( PostProcessingWriter::m_bPPWindow[m_rank] )
//         writeContactForceChains_Paraview( particles, LC, forceChainFilename 
// 		+ "_" + ossRK.str() + ".vtp", time );
    }
  }

  m_ParaviewCycleNumber++; 
}




// ----------------------------------------------------------------------------
// Updates obstacles indicator
void ParaviewPostProcessingWriter::updateObstaclesIndicator(
	double const& time,
  	double const& dt,
	Obstacle* obstacle )
{
  list<SimpleObstacle*> allObstacles = obstacle->getObstacles();
  for (list<SimpleObstacle*>::iterator iv=allObstacles.begin();
  	iv!=allObstacles.end();iv++) (*iv)->setIndicator( 0. );
  obstacle->updateIndicator( time, dt );
}




// ----------------------------------------------------------------------------
// Writes obstacles data
void ParaviewPostProcessingWriter::writeObstaclesPostProcessing_Paraview(
	list<SimpleObstacle*> const& allObstacles, string const& obsFilename )
{
  ofstream f( ( m_ParaviewFilename_dir + "/" + obsFilename ).c_str(), 
  	ios::out );
  list<SimpleObstacle*>::const_iterator il;

  f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
    	<< "byte_order=\"LittleEndian\" ";
  if ( m_binary ) f << "compressor=\"vtkZLibDataCompressor\"";
  f << ">" << endl;
  f << "<UnstructuredGrid>" << endl;
  int nbpts=0,nbcells=0,i;
  for (il=allObstacles.begin();il!=allObstacles.end();il++)
  {    
    nbpts += (*il)->numberOfPoints_PARAVIEW();
    nbcells += (*il)->numberOfCells_PARAVIEW();      
  }
  f << "<Piece NumberOfPoints=\"" << nbpts << "\""
    	<< " NumberOfCells=\"" << nbcells << "\">" << endl;
  f << "<Points>" << endl;
  f << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    list<Point3> ppp;
    list<Point3>::iterator ilpp;    
    start_output_binary( sizeof_Float32, 3*nbpts ) ;
    for (il=allObstacles.begin();il!=allObstacles.end();il++)
    {
      ppp = (*il)->get_polygonsPts_PARAVIEW();
      for (ilpp=ppp.begin();ilpp!=ppp.end();ilpp++)
        for (int comp=0;comp<3;++comp)
	  write_double_binary( (*ilpp)[comp] ) ;
    }
    flush_binary( f, "writeObstaclesPostProcessing_Paraview/Points" );      
  }
  else
    for (il=allObstacles.begin();il!=allObstacles.end();il++)
       (*il)->write_polygonsPts_PARAVIEW( f );
  f << "</DataArray>" << endl;
  f << "</Points>" << endl;

  list<int> connectivity,offsets,cellstype;
  list<int>::iterator ii;
  int firstpoint_globalnumber=0,last_offset=0;
  for (il=allObstacles.begin();il!=allObstacles.end();il++)
    (*il)->write_polygonsStr_PARAVIEW( connectivity,
    	offsets, cellstype, firstpoint_globalnumber, last_offset );   
  f << "<Cells>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"connectivity\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(connectivity.size()) ) ;
    for (ii=connectivity.begin();ii!=connectivity.end();ii++)
      write_int_binary( *ii );
    flush_binary( f, "writeObstaclesPostProcessing_Paraview/connectivity" );
  }
  else  
    for (ii=connectivity.begin();ii!=connectivity.end();ii++)
      f << *ii << " ";	
  f << endl;      
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"offsets\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(offsets.size()) ) ;
    for (ii=offsets.begin();ii!=offsets.end();ii++)
      write_int_binary( *ii );
    flush_binary( f, "writeObstaclesPostProcessing_Paraview/offsets" );
  }
  else  
    for (ii=offsets.begin();ii!=offsets.end();ii++)
      f << *ii << " ";	
  f << endl; 
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"types\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(cellstype.size()) ) ;
    for (ii=cellstype.begin();ii!=cellstype.end();ii++)
      write_int_binary( *ii );
    flush_binary( f, "writeObstaclesPostProcessing_Paraview/types" );
  }
  else  
    for (ii=cellstype.begin();ii!=cellstype.end();ii++)
      f << *ii << " ";	
  f << endl;
  f << "</DataArray>" << endl;
  f << "</Cells>" << endl;
  
  f << "<CellData Scalars=\"Indicator\">" << endl;
  f << "<DataArray type=\"Float32\" Name=\"Indicator\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">"; 
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary ) start_output_binary( sizeof_Float32, int(cellstype.size()) );
  for (il=allObstacles.begin();il!=allObstacles.end();il++)
  {
    double indic = (*il)->getIndicator();
    int nc =(*il)->numberOfCells_PARAVIEW(); 
    if ( m_binary ) for (i=0;i<nc;++i) write_double_binary( indic );
    else for (i=0;i<nc;++i) f << indic << " ";         
  }        
  if ( m_binary ) flush_binary( f, 
  	"writeObstaclesPostProcessing_Paraview/Indicator" ); 
  f << endl;
  f << "</DataArray>" << endl;   
  f << "</CellData>" << endl;    
    
  f << "</Piece>" << endl;
  f << "</UnstructuredGrid>" << endl;
  if ( m_binary )
  {
    f << "<AppendedData encoding=\"raw\">" << endl << "    _" ;
    f.write( BUFFER, OFFSET ) ;
    delete [] BUFFER ; BUFFER = 0 ;
    ALLOCATED = 0 ;
    OFFSET = 0 ;
    f << endl << "</AppendedData>" << endl;    
  }  
  f << "</VTKFile>" << endl;	
  f.close();	    
} 	




// ----------------------------------------------------------------------------
// Writes linked-cell grid data
void ParaviewPostProcessingWriter::writeLinkedCellPostProcessing_Paraview(
	LinkedCell const* LC, string const& partFilename,
	bool const& local )
{
  vector<double> 
  	x = local ? LC->local_coordinates( X ) : LC->global_coordinates( X ), 
	y = local ? LC->local_coordinates( Y ) : LC->global_coordinates( Y ),
  	z = local ? LC->local_coordinates( Z ) : LC->global_coordinates( Z );
  size_t nx = x.size(), ny = y.size(), nz = z.size(), i, j, k;
    
  ofstream f( ( m_ParaviewFilename_dir + "/" + partFilename ).c_str(), 
  	ios::out );
  
  f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
    	<< "byte_order=\"LittleEndian\" ";
  if ( m_binary ) f << "compressor=\"vtkZLibDataCompressor\"";
  f << ">" << endl;
  f << "<UnstructuredGrid>" << endl;
  size_t nbcells = nx * ny + nx * nz + ny * nz;
  size_t nbpts = 2 * nbcells; 
  f << "<Piece NumberOfPoints=\"" << nbpts << "\""
    	<< " NumberOfCells=\"" << nbcells << "\">" << endl;
  f << "<Points>" << endl;
  f << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Float32, 3*int(nbpts) ) ;
    // X-Y plane
    for (i=0;i<nx;i++)
      for (j=0;j<ny;j++)
      {
        write_double_binary( x[i] ) ;
        write_double_binary( y[j] ) ;	
        write_double_binary( z[0] ) ;	
        write_double_binary( x[i] ) ;	
        write_double_binary( y[j] ) ;	
        write_double_binary( z[nz-1] ) ;		
      }
      
    // X-Z plane
    for (i=0;i<nx;i++)
      for (k=0;k<nz;k++)
      {
        write_double_binary( x[i] ) ;
        write_double_binary( y[0] ) ;	
        write_double_binary( z[k] ) ;	
        write_double_binary( x[i] ) ;	
        write_double_binary( y[ny-1] ) ;	
        write_double_binary( z[k] ) ;	
      }
      
    // Y-Z plane
    for (j=0;j<ny;j++)
      for (k=0;k<nz;k++)
      {
        write_double_binary( x[0] ) ;
        write_double_binary( y[j] ) ;	
        write_double_binary( z[k] ) ;	
        write_double_binary( x[nx-1] ) ;	
        write_double_binary( y[j] ) ;	
        write_double_binary( z[k] ) ;	
      }
                        
    flush_binary( f, "writeLinkedCellPostProcessing_Paraview/Points" );
  }
  else
  {
    // X-Y plane
    for (i=0;i<nx;i++)
      for (j=0;j<ny;j++)
      {
        f << x[i] << " " << y[j] << " " << z[0] << endl;
        f << x[i] << " " << y[j] << " " << z[nz-1] << endl;	
      }

    // X-Z plane
    for (i=0;i<nx;i++)
      for (k=0;k<nz;k++)
      {
        f << x[i] << " " << y[0] << " " << z[k] << endl;
        f << x[i] << " " << y[ny-1] << " " << z[k] << endl;	
      }    

    // Y-Z plane
    for (j=0;j<ny;j++)
      for (k=0;k<nz;k++)
      {
        f << x[0] << " " << y[j] << " " << z[k] << endl;
        f << x[nx-1] << " " << y[j] << " " << z[k] << endl;	
      }
  }           
  f << "</DataArray>" << endl;
  f << "</Points>" << endl;
  f << "<Cells>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"connectivity\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  int counter = 0;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, 2*int(nbpts) ) ;
    // X-Y plane
    for (i=0;i<nx;i++)
      for (j=0;j<ny;j++,counter+=2)
      {
        write_int_binary( counter );
        write_int_binary( counter + 1 );	
      }
      
    // X-Z plane
    for (i=0;i<nx;i++)
      for (k=0;k<nz;k++,counter+=2)
      {
        write_int_binary( counter );
        write_int_binary( counter + 1 );	
      }
	
    // Y-Z plane
    for (j=0;j<ny;j++)
      for (k=0;k<nz;k++,counter+=2)
      {
        write_int_binary( counter );
        write_int_binary( counter + 1 );	
      }
      
    flush_binary( f, "writeLinkedCellPostProcessing_Paraview/connectivity" );
  }
  else
  {  
    // X-Y plane
    for (i=0;i<nx;i++)
      for (j=0;j<ny;j++,counter+=2)
        f << counter << " " << counter + 1 << " ";
      
    // X-Z plane
    for (i=0;i<nx;i++)
      for (k=0;k<nz;k++,counter+=2)
        f << counter << " " << counter + 1 << " ";
	
    // Y-Z plane
    for (j=0;j<ny;j++)
      for (k=0;k<nz;k++,counter+=2)
        f << counter << " " << counter + 1 << " ";
  }	
  f << endl;     
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"offsets\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(nbpts) ) ;
    for (i=0;i<nbpts;i++) write_int_binary( 2*(int(i)+1) );
    flush_binary( f, "writeLinkedCellPostProcessing_Paraview/offsets" );
  }
  else 
    for (i=0;i<nbpts;i++) f << 2*(i+1) << " ";	
  f << endl; 
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"types\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(nbpts) ) ;
    int type = 3;
    for (i=0;i<nbpts;i++) write_int_binary( type );
    flush_binary( f, "writeLinkedCellPostProcessing_Paraview/types" );
  }
  else  
    for (i=0;i<nbpts;i++) f << "3 ";
  f << endl;
  f << "</DataArray>" << endl;
  f << "</Cells>" << endl;
  f << "</Piece>" << endl;
  f << "</UnstructuredGrid>" << endl;
  if ( m_binary )
  {
    f << "<AppendedData encoding=\"raw\">" << endl << "    _" ;
    f.write( BUFFER, OFFSET ) ;
    delete [] BUFFER ; BUFFER = 0 ;
    ALLOCATED = 0 ;
    OFFSET = 0 ;
    f << endl << "</AppendedData>" << endl;    
  }  
  f << "</VTKFile>" << endl;	
  f.close();
}




// ----------------------------------------------------------------------------
// Writes insertion windows data
void ParaviewPostProcessingWriter::writeInsertionPostProcessing_Paraview(
   	vector<Window> const& insert_windows,
  	string const& partFilename )
{
  vector<Window>::const_iterator iv;
  list<RigidBody*> iwlist;
  list<RigidBody*>::iterator il = iwlist.begin(); 
  Transform gcwindow; 
  Convex* ccw = NULL;
  RigidBody* ffw = NULL;
  Vector3 v_axis;
  Matrix mrot;
  double Le = 0., Lh = 0., Lt = 0., angle = 0., 
  	mean_radius = 0., thickness = 0. ;
  Point3 panelCenter, center;
  size_t nbPanels = 32;
  double polygonAngle = 2. * PI / double(nbPanels);
  
  for (iv=insert_windows.begin();iv!=insert_windows.end();iv++)
  {
    switch( iv->ftype )
    {
      case WINDOW_BOX:
        gcwindow.setOrigin( 0.5 * ( iv->ptA + iv->ptB ) );
        ccw = new Box( 0.5 * ( iv->ptB - iv->ptA ) );
        ffw = new RigidBody( ccw, gcwindow );
        iwlist.push_back( ffw );
        break;
      
      case WINDOW_CYLINDER:
        switch ( iv->axisdir )
        {
          case X: 
            v_axis[X] = iv->height; 
            mrot.setValue( cos( 0.5 * PI ), -sin( 0.5 * PI ), 0.,
                sin( 0.5 * PI ), cos( 0.5 * PI ), 0.,
                0., 0., 1. );
            break;

          case Y: 
            v_axis[Y] = iv->height; 
            break;
  
          default: 
            v_axis[Z] = iv->height; 
            mrot.setValue( 1., 0., 0.,
                0., cos( 0.5 * PI ), -sin( 0.5 * PI ),
                0., sin( 0.5 * PI ), cos( 0.5 * PI ) );
            break;
        }
        gcwindow.setOrigin( iv->ptA + 0.5 * v_axis );
        ccw = new Cylinder( iv->radius, iv->height );
        ffw = new RigidBody( ccw, gcwindow );
        ffw->getTransform()->setBasis( mrot );
        iwlist.push_back( ffw );
        break;	

      case WINDOW_ANNULUS:
        mean_radius = iv->radius_int + 0.5 * ( iv->radius - iv->radius_int ) ;
        thickness = iv->radius - iv->radius_int	;
        center = iv->ptA;
        Le = thickness;       
        Lt = iv->radius * tan( polygonAngle ); 
        Lh = iv->height;
        switch ( iv->axisdir )
        {
          case X: 
            center[X] += 0.5 * iv->height;
            for (size_t iNb=0; iNb!=nbPanels; iNb++)
            {
              angle = 2. * PI * double(iNb) / double(nbPanels);
              panelCenter[X] = center[X];
              panelCenter[Y] = center[Y] + mean_radius * cos(angle);
              panelCenter[Z] = center[Z] + mean_radius * sin(angle);
              gcwindow.setOrigin( panelCenter );
 
              mrot.setValue( 1., 0., 0.,
                  0., cos(angle), -sin(angle),
                  0., sin(angle), cos(angle) );
              gcwindow.setBasis( mrot );
              ccw = new Box( Lh, Le, Lt );
              ffw = new RigidBody( ccw, gcwindow );
              iwlist.push_back( ffw );
            }
            break;

          case Y: 
            center[Y] += 0.5 * iv->height;
            for (size_t iNb=0; iNb!=nbPanels; iNb++)
            {
              angle = 2. * PI * double(iNb) / double(nbPanels);
              panelCenter[X] = center[X] + mean_radius * sin(angle);
              panelCenter[Y] = center[Y];
              panelCenter[Z] = center[Z] + mean_radius * cos(angle);
              gcwindow.setOrigin( panelCenter );

              mrot.setValue( cos(angle), 0., sin(angle),
                  0., 1., 0.,
                  -sin(angle), 0., cos(angle) );
              gcwindow.setBasis( mrot );
              ccw = new Box( Lt, Lh, Le );
              ffw = new RigidBody( ccw, gcwindow );
              iwlist.push_back( ffw );
            }
            break;

          default:
            center[Z] += 0.5 * iv->height;
            for (size_t iNb=0; iNb!=nbPanels; iNb++)
            {
              angle = 2. * PI * double(iNb) / double(nbPanels);
              panelCenter[X] = center[X] + mean_radius * cos(angle);
              panelCenter[Y] = center[Y] + mean_radius * sin(angle);
              panelCenter[Z] = center[Z];
              gcwindow.setOrigin( panelCenter );

              mrot.setValue( cos(angle), -sin(angle), 0.,
                  sin(angle), cos(angle), 0.,
                  0., 0., 1. );
              gcwindow.setBasis( mrot );
              ccw = new Box( Le, Lt, Lh );
              ffw = new RigidBody( ccw, gcwindow );
              iwlist.push_back( ffw );
            }
            break;
        }
        break;

      case WINDOW_LINE:
        gcwindow = Segment::computeTransform( 
            0.5 * ( iv->ptB - iv->ptA ), 0.5 * ( iv->ptA + iv->ptB ) );
        ccw = new Segment( Norm( iv->ptB - iv->ptA ) );
        ffw = new RigidBody( ccw, gcwindow );
        iwlist.push_back( ffw );
        break;

      default:
        break;
    }
  }

  ofstream f( ( m_ParaviewFilename_dir + "/" + partFilename + ".vtu" ).c_str(),
      ios::out );
  
  f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
    << "byte_order=\"LittleEndian\" ";
  if ( m_binary ) f << "compressor=\"vtkZLibDataCompressor\"";
  f << ">" << endl;
  f << "<UnstructuredGrid>" << endl;
  int nbpts = 0, nbcells = 0;
  for (il=iwlist.begin();il!=iwlist.end();il++)
  {
    nbpts += (*il)->getConvex()->numberOfPoints_PARAVIEW();
    nbcells += (*il)->getConvex()->numberOfCells_PARAVIEW();
  }
  f << "<Piece NumberOfPoints=\"" << nbpts << "\""
    << " NumberOfCells=\"" << nbcells << "\">" << endl;
  f << "<Points>" << endl;
  f << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    list<Point3> ppp;
    list<Point3>::iterator ilpp;
    start_output_binary( sizeof_Float32, 3*nbpts ) ;
    for (il=iwlist.begin();il!=iwlist.end();il++)
    {
      ppp = (*il)->get_polygonsPts_PARAVIEW();
      for (ilpp=ppp.begin();ilpp!=ppp.end();ilpp++)
        for (int comp=0;comp<3;++comp)
          write_double_binary( (*ilpp)[comp] ) ;
    }
    flush_binary( f, "writeInsertionPostProcessing_Paraview/Points" );
  }
  else
    for (il=iwlist.begin();il!=iwlist.end();il++)
      (*il)->write_polygonsPts_PARAVIEW( f );
  f << "</DataArray>" << endl;
  f << "</Points>" << endl;

  list<int> connectivity,offsets,cellstype;
  list<int>::iterator ii;
  int firstpoint_globalnumber=0,last_offset=0;
  for (il=iwlist.begin();il!=iwlist.end();il++)
    (*il)->getConvex()->write_polygonsStr_PARAVIEW( connectivity,
        offsets, cellstype, firstpoint_globalnumber, last_offset );
  f << "<Cells>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"connectivity\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(connectivity.size()) ) ;
    for (ii=connectivity.begin();ii!=connectivity.end();ii++)
      write_int_binary( *ii );
    flush_binary( f, "writeInsertionPostProcessing_Paraview/connectivity" );
  }
  else  
    for (ii=connectivity.begin();ii!=connectivity.end();ii++)
      f << *ii << " ";
  f << endl;
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"offsets\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(offsets.size()) ) ;
    for (ii=offsets.begin();ii!=offsets.end();ii++)
      write_int_binary( *ii );
    flush_binary( f, "writeInsertionPostProcessing_Paraview/offsets" );
  }
  else  
    for (ii=offsets.begin();ii!=offsets.end();ii++)
      f << *ii << " ";	
  f << endl; 
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"types\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(cellstype.size()) ) ;
    for (ii=cellstype.begin();ii!=cellstype.end();ii++)
      write_int_binary( *ii );
    flush_binary( f, "writeInsertionPostProcessing_Paraview/types" );
  }
  else
    for (ii=cellstype.begin();ii!=cellstype.end();ii++)
      f << *ii << " ";	
  f << endl;
  f << "</DataArray>" << endl;
  f << "</Cells>" << endl;     
  f << "</Piece>" << endl;
  f << "</UnstructuredGrid>" << endl;
  if ( m_binary )
  {
    f << "<AppendedData encoding=\"raw\">" << endl << "    _" ;
    f.write( BUFFER, OFFSET ) ;
    delete [] BUFFER ; BUFFER = 0 ;
    ALLOCATED = 0 ;
    OFFSET = 0 ;
    f << endl << "</AppendedData>" << endl;    
  }  
  f << "</VTKFile>" << endl;	
  f.close();
  
  for (il=iwlist.begin();il!=iwlist.end();il++) delete *il;
}




// ----------------------------------------------------------------------------
// Writes periodic boundary data
void ParaviewPostProcessingWriter::writePeriodicBoundaryPostProcessing_Paraview(
   	LinkedCell const* LC, string const& partFilename )
{
  list<RigidBody*> iwlist;
  list<RigidBody*>::iterator il = iwlist.begin(); 
  Convex* ccw = NULL;
  RigidBody* ffw = NULL;
  Transform tt; 
  double ox, oy, oz, lx, ly, lz;
  
  App::get_origin( ox, oy, oz );
  App::get_size( lx, ly, lz );
    
  // Periodicity in X
  if ( LC->isPeriodic( 0 ) )
  {
    tt.setOrigin( ox, oy + ly / 2., oz + lz / 2. );
    ccw = new Box( 0., ly, lz );
    ffw = new RigidBody( ccw, tt );
    iwlist.push_back( ffw );
    
    tt.setOrigin( ox + lx, oy + ly / 2., oz + lz / 2. );
    ccw = new Box( 0., ly, lz );
    ffw = new RigidBody( ccw, tt );
    iwlist.push_back( ffw );    
  }
  
  // Periodicity in Y
  if ( LC->isPeriodic( 1 ) )
  {
    tt.setOrigin( ox + lx / 2., oy, oz + lz / 2. );
    ccw = new Box( lx, 0., lz );
    ffw = new RigidBody( ccw, tt );
    iwlist.push_back( ffw );
    
    tt.setOrigin( ox + lx / 2., oy + ly, oz + lz / 2. );
    ccw = new Box( lx, 0., lz );
    ffw = new RigidBody( ccw, tt );
    iwlist.push_back( ffw );   
  } 
  
  // Periodicity in Z  
  if ( GrainsBuilderFactory::getContext() == DIM_3 )
    if ( LC->isPeriodic( 2 ) )
    {
      tt.setOrigin( ox + lx / 2., oy + ly / 2., oz );
      ccw = new Box( lx, ly, 0. );
      ffw = new RigidBody( ccw, tt );
      iwlist.push_back( ffw );
    
      tt.setOrigin( ox + lx / 2., oy + ly / 2., oz + lz );
      ccw = new Box( lx, ly, 0. );
      ffw = new RigidBody( ccw, tt );
      iwlist.push_back( ffw );   
    }       

  ofstream f( ( m_ParaviewFilename_dir + "/" + partFilename + ".vtu" ).c_str(),
      ios::out );
  
  f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
    << "byte_order=\"LittleEndian\" ";
  if ( m_binary ) f << "compressor=\"vtkZLibDataCompressor\"";
  f << ">" << endl;
  f << "<UnstructuredGrid>" << endl;
  int nbpts = 0, nbcells = 0;
  for (il=iwlist.begin();il!=iwlist.end();il++)
  {
    nbpts += (*il)->getConvex()->numberOfPoints_PARAVIEW();
    nbcells += (*il)->getConvex()->numberOfCells_PARAVIEW();
  }
  f << "<Piece NumberOfPoints=\"" << nbpts << "\""
    << " NumberOfCells=\"" << nbcells << "\">" << endl;
  f << "<Points>" << endl;
  f << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    list<Point3> ppp;
    list<Point3>::iterator ilpp;
    start_output_binary( sizeof_Float32, 3*nbpts ) ;
    for (il=iwlist.begin();il!=iwlist.end();il++)
    {
      ppp = (*il)->get_polygonsPts_PARAVIEW();
      for (ilpp=ppp.begin();ilpp!=ppp.end();ilpp++)
        for (int comp=0;comp<3;++comp)
          write_double_binary( (*ilpp)[comp] ) ;
    }
    flush_binary( f, "writeInsertionPostProcessing_Paraview/Points" );
  }
  else
    for (il=iwlist.begin();il!=iwlist.end();il++)
      (*il)->write_polygonsPts_PARAVIEW( f );
  f << "</DataArray>" << endl;
  f << "</Points>" << endl;

  list<int> connectivity,offsets,cellstype;
  list<int>::iterator ii;
  int firstpoint_globalnumber=0,last_offset=0;
  for (il=iwlist.begin();il!=iwlist.end();il++)
    (*il)->getConvex()->write_polygonsStr_PARAVIEW( connectivity,
        offsets, cellstype, firstpoint_globalnumber, last_offset );
  f << "<Cells>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"connectivity\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(connectivity.size()) ) ;
    for (ii=connectivity.begin();ii!=connectivity.end();ii++)
      write_int_binary( *ii );
    flush_binary( f, "writeInsertionPostProcessing_Paraview/connectivity" );
  }
  else  
    for (ii=connectivity.begin();ii!=connectivity.end();ii++)
      f << *ii << " ";
  f << endl;
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"offsets\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(offsets.size()) ) ;
    for (ii=offsets.begin();ii!=offsets.end();ii++)
      write_int_binary( *ii );
    flush_binary( f, "writeInsertionPostProcessing_Paraview/offsets" );
  }
  else  
    for (ii=offsets.begin();ii!=offsets.end();ii++)
      f << *ii << " ";	
  f << endl; 
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"types\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(cellstype.size()) ) ;
    for (ii=cellstype.begin();ii!=cellstype.end();ii++)
      write_int_binary( *ii );
    flush_binary( f, "writeInsertionPostProcessing_Paraview/types" );
  }
  else
    for (ii=cellstype.begin();ii!=cellstype.end();ii++)
      f << *ii << " ";	
  f << endl;
  f << "</DataArray>" << endl;
  f << "</Cells>" << endl;     
  f << "</Piece>" << endl;
  f << "</UnstructuredGrid>" << endl;
  if ( m_binary )
  {
    f << "<AppendedData encoding=\"raw\">" << endl << "    _" ;
    f.write( BUFFER, OFFSET ) ;
    delete [] BUFFER ; BUFFER = 0 ;
    ALLOCATED = 0 ;
    OFFSET = 0 ;
    f << endl << "</AppendedData>" << endl;    
  }  
  f << "</VTKFile>" << endl;	
  f.close();
  
  for (il=iwlist.begin();il!=iwlist.end();il++) delete *il;
}	




// ----------------------------------------------------------------------------
// Writes a pvtp file 
void ParaviewPostProcessingWriter::writePVTP_Paraview( string const& filename,
  	list<string> const* pointVector,
	list<string> const* pointScalar,
	list<string> const* cellScalar )
{
  ofstream f( ( m_ParaviewFilename_dir + "/" + filename + ".pvtp" ).c_str(),
  	ios::out );
  f << "<?xml version=\"1.0\"?>" << endl; 
  f << "<VTKFile type=\"PPolyData\" version=\"0.1\" "
    	<< "byte_order=\"LittleEndian\" ";
  if ( m_binary ) f << "compressor=\"vtkZLibDataCompressor\"";
  f << ">" << endl;
  f << "<PPolyData GhostLevel=\"0\">" << endl;  
  f << "<PPoints>" << endl;
  f << "<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>" << endl;
  f << "</PPoints>" << endl;
  f << "<PPointData Scalars=\"F_N\">" << endl;
  f << "<PDataArray type=\"Float32\" Name=\"F_N\"/>" << endl;
  f << "</PPointData>" << endl;
  f << "<PLines>" << endl;
  f << "<PDataArray Name=\"connectivity\" type=\"Int32\" format=\"ascii\">"
  	<< endl;
  f << "</PDataArray>" << endl;
  f << "<PDataArray Name=\"offsets\" type=\"Int32\" format=\"ascii\">" << endl;
  f << "</PDataArray>" << endl;
  f << "</PLines>" << endl;
  for (int i=0;i<m_nprocs;++i)
  {
    // Does this processor have to write down outputs ?
    if( PostProcessingWriter::m_bPPWindow[i] )
    {
      f << "<Piece Source=\"" << filename << "_" << i << ".vtp\">" << endl;
      f << "</Piece>" << endl; 
    }
  }
  f << "</PPolyData>" << endl;
  f << "</VTKFile>" << endl;  
  f.close();  
}




// ----------------------------------------------------------------------------
// Writes a pvtu file
void ParaviewPostProcessingWriter::writePVTU_Paraview( string const& filename,
  	list<string> const* pointVector,
	list<string> const* pointScalar,
	list<string> const* cellScalar )
{
  list<string>::const_iterator il;

  ofstream f( ( m_ParaviewFilename_dir + "/" + filename + ".pvtu" ).c_str(),
  	ios::out );
  f << "<?xml version=\"1.0\"?>" << endl; 
  f << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" "
    	<< "byte_order=\"LittleEndian\" ";
  if ( m_binary ) f << "compressor=\"vtkZLibDataCompressor\"";
  f << ">" << endl;
  f << "<PUnstructuredGrid GhostLevel=\"0\">" << endl;  
  f << "<PPoints>" << endl;
  f << "<PDataArray NumberOfComponents=\"3\" type=\"Float32\" format=\"ascii\">"
  	<< endl;
  f << "</PDataArray>" << endl;
  f << "</PPoints>" << endl;  
  f << "<PCells>" << endl;
  f << "<PDataArray Name=\"connectivity\" type=\"Int32\" format=\"ascii\">"
  	<< endl;
  f << "</PDataArray>" << endl;
  f << "<PDataArray Name=\"offsets\" type=\"Int32\" format=\"ascii\">" << endl;
  f << "</PDataArray>" << endl;
  f << "<PDataArray Name=\"types\" type=\"Int32\" format=\"ascii\">" << endl;
  f << "</PDataArray>" << endl;
  f << "</PCells>" << endl;    
  if ( pointVector->size() || pointScalar->size())
  {
    f << "<PPointData";
    if ( pointVector->size() )
    {
      il = pointVector->begin();
      f << " Vectors=\"" << *il;
      il++;
      for ( ;il!=pointVector->end();il++) f << "," << *il;
      f << "\"";      
    }
    if ( pointScalar->size() )
    {
      il = pointScalar->begin();
      f << " Scalars=\"" << *il;
      il++;
      for ( ;il!=pointScalar->end();il++) f << "," << *il;
      f << "\"";      
    }
    f << ">" << endl;
    for (il = pointVector->begin();il!=pointVector->end();il++)
    {    
      f << "<PDataArray Name=\"" << *il 
    	<< "\" NumberOfComponents=\"3\" type=\"Float32\""
    	<< " format=\"ascii\">" << endl;
      f << "</PDataArray>" << endl;
    }
    for (il = pointScalar->begin();il!=pointScalar->end();il++)
    {    
      f << "<PDataArray Name=\"" << *il 
    	<< "\" type=\"Float32\" format=\"ascii\">" << endl;
      f << "</PDataArray>" << endl;
    }    
    f << "</PPointData>" << endl;    
  }      
  if ( cellScalar->size() )
  {        
    f << "<PCellData Scalars=\"";
    il = cellScalar->begin();
    f << *il;
    il++;
    for ( ;il!=cellScalar->end();il++) f << "," << *il;
    f << "\">"; 
    for (il = cellScalar->begin();il!=cellScalar->end();il++)
    {     
      f << "<PDataArray Name=\"" << *il << "\" type=\"Float32\"" 
    	<< " format=\"ascii\">" << endl;
      f << "</PDataArray>" << endl;
    }
    f << "</PCellData>" << endl;  
  }  
  for (int i=0;i<m_nprocs;++i)
  {
    // Does this processor have to write down outputs ?
    if( PostProcessingWriter::m_bPPWindow[i] )
    {
      f << "<Piece Source=\"" << filename << "_" << i << ".vtu\">" << endl;
      f << "</Piece>" << endl; 
    }
  }
  f << "</PUnstructuredGrid>" << endl;
  f << "</VTKFile>" << endl;  
  f.close();  
}




// ----------------------------------------------------------------------------
// Writes components involved in a displacement or a contact error
void ParaviewPostProcessingWriter::writeErreurComponentsPostProcessing(
	string const& filename,
  	list<Component*> const& errcomposants )
{
  list<Component*>::const_iterator compo;

  ostringstream ossRK;
  ossRK << m_rank; 
    
  ofstream f( ( m_ParaviewFilename_dir + "/" + filename + "_" + ossRK.str()
  	+ ".vtu" ).c_str(), ios::out );  
  f << "<?xml version=\"1.0\"?>" << endl;
  f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
    	<< "byte_order=\"LittleEndian\" ";
  if ( m_binary ) f << "compressor=\"vtkZLibDataCompressor\"";
  f << ">" << endl;
  f << "<UnstructuredGrid>" << endl;
  int nbpts=0,nbcells=0;
  for (compo=errcomposants.begin();compo!=errcomposants.end();compo++)
  {
    nbpts += (*compo)->numberOfPoints_PARAVIEW();
    nbcells += (*compo)->numberOfCells_PARAVIEW();
  }
  f << "<Piece NumberOfPoints=\"" << nbpts << "\""
    	<< " NumberOfCells=\"" << nbcells << "\">" << endl;
  f << "<Points>" << endl;
  f << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    list<Point3> ppp;
    list<Point3>::iterator ilpp;    
    start_output_binary( sizeof_Float32, 3*nbpts ) ;
    for (compo=errcomposants.begin();compo!=errcomposants.end();compo++)
    {
      ppp = (*compo)->get_polygonsPts_PARAVIEW();
      for (ilpp=ppp.begin();ilpp!=ppp.end();ilpp++)
        for (int comp=0;comp<3;++comp)
	  write_double_binary( (*ilpp)[comp] ) ;
    }
    flush_binary( f, "writeErreurComponentsPostProcessing/Points" );      
  }
  else
    for (compo=errcomposants.begin();compo!=errcomposants.end();compo++)
      (*compo)->write_polygonsPts_PARAVIEW(f);
      //(*compo)->getRigidBody()->write_polygonsPts_PARAVIEW(f);
  f << "</DataArray>" << endl;
  f << "</Points>" << endl;
  list<int> connectivity,offsets,cellstype;
  list<int>::iterator ii;
  int firstpoint_globalnumber=0,last_offset=0;
  for (compo=errcomposants.begin();compo!=errcomposants.end();compo++)
    (*compo)->write_polygonsStr_PARAVIEW( connectivity,
    	offsets, cellstype, firstpoint_globalnumber, last_offset ); 
  f << "<Cells>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"connectivity\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(connectivity.size()) ) ;
    for (ii=connectivity.begin();ii!=connectivity.end();ii++)
      write_int_binary( *ii );
    flush_binary( f, "writeErreurComponentsPostProcessing/connectivity" );
  }
  else  
    for (ii=connectivity.begin();ii!=connectivity.end();ii++)
      f << *ii << " ";	
  f << endl;      
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"offsets\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(offsets.size()) ) ;
    for (ii=offsets.begin();ii!=offsets.end();ii++)
      write_int_binary( *ii );
    flush_binary( f, "writeErreurComponentsPostProcessing/offsets" );
  }
  else  
    for (ii=offsets.begin();ii!=offsets.end();ii++)
      f << *ii << " ";	
  f << endl; 
  f << "</DataArray>" << endl;
  f << "<DataArray type=\"Int32\" Name=\"types\" ";
  if ( m_binary ) f << "offset=\"" << OFFSET << "\" format=\"appended\">";
  else f << "format=\"ascii\">";
  f << endl;
  if ( m_binary )
  {
    start_output_binary( sizeof_Int32, int(cellstype.size()) ) ;
    for (ii=cellstype.begin();ii!=cellstype.end();ii++)
      write_int_binary( *ii );
    flush_binary( f, "writeErreurComponentsPostProcessing/types" );
  }
  else  
    for (ii=cellstype.begin();ii!=cellstype.end();ii++)
      f << *ii << " ";	
  f << endl;
  f << "</DataArray>" << endl;
  f << "</Cells>" << endl;
  f << "</Piece>" << endl;
  f << "</UnstructuredGrid>" << endl;
  if ( m_binary )
  {
    f << "<AppendedData encoding=\"raw\">" << endl << "    _" ;
    f.write( BUFFER, OFFSET ) ;
    delete [] BUFFER ; BUFFER = 0 ;
    ALLOCATED = 0 ;
    OFFSET = 0 ;
    f << endl << "</AppendedData>" << endl;    
  }  
  f << "</VTKFile>" << endl;	
  f.close();	    
}




// ----------------------------------------------------------------------------
// Delete all result files
void ParaviewPostProcessingWriter::clearResultFiles() const
{
  if ( m_rank == 0 ) 
  {
    string cmd = "bash " + GrainsExec::m_GRAINS_HOME 
     	+ "/Tools/ExecScripts/Paraview_clear.exec " + m_ParaviewFilename_dir +
	" " + m_ParaviewFilename;
    GrainsExec::m_return_syscmd = system( cmd.c_str() );
  }   
}




// ---------------------------------------------------------------------------
// Sets the initial cycle number
void ParaviewPostProcessingWriter::setInitialCycleNumber( const int& cycle0 ) 
{ 
  m_ParaviewCycleNumber = cycle0; 
  m_initialCycleNumber_forced = true ; 
}




// ----------------------------------------------------------------------------
// Writes a long ostringstream character by character (in case the
// classical fileOUT << oss.str() crashes)
void ParaviewPostProcessingWriter::writeBigOSS( ofstream& fileOUT, 
	ostringstream const& oss )
{
  string str = oss.str();      
  for ( string::iterator it=str.begin(); it!=str.end(); ++it)
    fileOUT << *it;  
}		




// ----------------------------------------------------------------------------
// Gets the post-processing writer type
string ParaviewPostProcessingWriter::getPostProcessingWriterType() const
{
  return ( "Paraview" );
}




// ----------------------------------------------------------------------------
// Methods to write binary data
void ParaviewPostProcessingWriter:: start_output_binary( int size, int number )
{
  int current_output_size = size*number ;
//   unsigned long ncomp = current_output_size + (current_output_size+999)/1000 
//   	+ 12 + sizeof_Int32 ;
  int ncomp = current_output_size + (current_output_size+999)/1000 
  	+ 12 + sizeof_Int32 ;	
  check_allocated_binary( ncomp ) ;
  CURRENT_LENGTH = store_int_binary(current_output_size) ;
}




void ParaviewPostProcessingWriter:: write_double_binary( double val )  
{
  *((float*)&(BUFFER[OFFSET])) = (float)val ;
  OFFSET += sizeof_Float32  ;
}




void ParaviewPostProcessingWriter:: write_int_binary( int val )  
{
  store_int_binary(val) ;
}




int ParaviewPostProcessingWriter:: store_int_binary( int val )  
{
  int result = OFFSET ;
  *((int*)&(BUFFER[OFFSET])) = val ;
  OFFSET += sizeof_Int32  ;
  return result ;
}




void ParaviewPostProcessingWriter:: check_allocated_binary( int size )  
{
  if ( OFFSET + size >= ALLOCATED ) 
  {
    int new_size = max( 2*ALLOCATED, (int)1024 ) ;
    new_size = max( new_size, 2*(OFFSET+size) ) ;
    new_size = 4 * ( new_size/4 +1 ) ; // alignment on 4 bytes
      
    char * new_buffer = new char [ new_size ] ;
    for( int i=0 ;i<OFFSET ;i++ ) new_buffer[i] = BUFFER[i] ;
    if( BUFFER!=0 ) delete [] BUFFER ;
    BUFFER = new_buffer ;
    ALLOCATED = new_size ;      
  }
}




void ParaviewPostProcessingWriter:: flush_binary( std::ofstream& file, 
	string const& calling )  
{
  compress_segment_binary( CURRENT_LENGTH, calling ) ;         
//  file << endl ;
}




void ParaviewPostProcessingWriter:: compress_segment_binary( int seg,
	string const& calling )  
{
   static int BlockSize = 32768 ;
   int size = (int)(*((int*)&BUFFER[seg])) ;
   
   int numFullBlocks = size / BlockSize;
   int lastBlockSize = size % BlockSize;
   int numBlocks = numFullBlocks + (lastBlockSize?1:0);

   int headerLength = numBlocks+3;

   int * CompressionHeader = new int[headerLength];
   CompressionHeader[0] = numBlocks;
   CompressionHeader[1] = BlockSize;
   CompressionHeader[2] = lastBlockSize;

   unsigned long encoded_buff_size = max(BlockSize,size)  ;
   unsigned char* encoded_buff = new unsigned char [ encoded_buff_size ] ;
   int encoded_offset = 0 ;
   for( int block=0 ; block<numBlocks ; block++ )
   {
      int buffer_start = seg + sizeof_Int32 + block*BlockSize ;
      int length = ( block+1<numBlocks || !lastBlockSize ? 
      	BlockSize : lastBlockSize ) ;
      unsigned char* to_encode = (unsigned char *)(&BUFFER[buffer_start]) ;
      unsigned char* encoded = &encoded_buff[encoded_offset] ;
      unsigned long ncomp = encoded_buff_size - encoded_offset ;

      if ( compress2( (Bytef*)encoded,
                   &ncomp,
                   (const Bytef*)to_encode,
                   length,
                   Z_DEFAULT_COMPRESSION) != Z_OK )
      {
         cout << "Zlib error while compressing data." << endl;
	 cout << "from " << calling << endl;
	 cout << "Details : block = " << block << "  numBlocks = " <<
	 	numBlocks << "  length = " << length << "  ncomp = "
		<< ncomp << endl;
	 exit(0);
      }
      CompressionHeader[3+block] = int(ncomp) ;
      encoded_offset += int(ncomp) ;      
   }
   
   OFFSET = seg ;
   check_allocated_binary( headerLength * sizeof_Int32 + encoded_offset ) ;
   
   for(int i=0 ; i<headerLength ; i++ )
      store_int_binary( CompressionHeader[i] ) ;     

   for(int i=0 ; i<encoded_offset ; i++ )
      BUFFER[OFFSET++] = encoded_buff[i] ;

   if( OFFSET%4 != 0 )
      OFFSET = 4*( OFFSET/4 +1 ) ; // Re-alignment
   
   delete [] CompressionHeader ;
   delete [] encoded_buff ;
}
// End of Methods to write binary data
// ----------------------------------------------------------------------------
