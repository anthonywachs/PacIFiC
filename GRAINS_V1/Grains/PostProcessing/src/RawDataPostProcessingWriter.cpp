#include "GrainsExec.hh"
#include "RawDataPostProcessingWriter.hh"
#include "Particle.hh"
#include "Obstacle.hh"
#include "GrainsMPIWrapper.hh"



// ----------------------------------------------------------------------------
// Default constructor
RawDataPostProcessingWriter::RawDataPostProcessingWriter()
{}




// ----------------------------------------------------------------------------
// Constructor with XML node, rank and number of processes as input parameters
RawDataPostProcessingWriter::RawDataPostProcessingWriter( DOMNode* dn,
    int const& rank_, int const& nbranks_, bool const& verbose )
  : PostProcessingWriter( dn, rank_, nbranks_ )
  , m_binary( false )
  , m_ndigits( 6 )  
{ 
  m_filerootname = ReaderXML::getNodeAttr_String( dn, "Name" );
  if ( ReaderXML::hasNodeAttr( dn, "WritingMode" ) )
  { 
    string sm_binary = ReaderXML::getNodeAttr_String( dn, "WritingMode" );
    if ( sm_binary == "Binary" ) m_binary = true;
  }  
 
  if ( m_rank == 0 && verbose )
  {
    cout << GrainsExec::m_shift9 << "Type = RawData" << endl;
    cout << GrainsExec::m_shift12 << "Output file name = " 
    	<< m_filerootname << endl;
    cout << GrainsExec::m_shift12 << "Writing mode = " 
    	<< ( m_binary ? "Binary" : "Text" ) << endl;	
  }
}




// ----------------------------------------------------------------------------
// Destructor
RawDataPostProcessingWriter::~RawDataPostProcessingWriter()
{}




// ----------------------------------------------------------------------------
// Initializes the post-processing writer
void RawDataPostProcessingWriter::PostProcessing_start(
    double const& time, 
    double const& dt,
    list<Particle*> const* particles,
    list<Particle*> const* inactiveparticles,
    list<Particle*> const* periodic_clones,
    vector<Particle*> const* referenceParticles,
    Obstacle* obstacle,
    LinkedCell const* LC,
    vector<Window> const& insert_windows )
{
  GrainsMPIWrapper const* wrapper = GrainsExec::getComm() ;
  size_t nb_total_part = GrainsExec::getTotalNumberPhysicalParticles() ;
  int type = 0;

  if ( m_rank == 0 )
  {
    ios_base::openmode mode = ios::app;
    if ( GrainsExec::m_ReloadType == "new" ) 
    {
      mode = ios::out;
      clearResultFiles();
    }
    prepareResultFiles( mode );
  }     

  // In parallel mode
  if ( wrapper )
  {
    vector<int>* types_Global = NULL;
    vector< vector<double> >* data_Global = NULL;    
    
    // Gather particles class from every proc on master proc
    types_Global = 
        wrapper->GatherParticlesClass_PostProcessing( *particles,
            nb_total_part );

    // Write down particles class only once at the begining
    if ( m_rank == 0 )
    {
      if ( m_binary )
      { 
        for (size_t i=0; i<nb_total_part; i++)
	{
	  type = (*types_Global)[i];
	  m_particle_class.write( reinterpret_cast<char*>( &type ), 
	  	sizeof(int) );
	}      
      }
      else
      {     
        for (size_t i=0; i<nb_total_part; i++)
          m_particle_class << (*types_Global)[i] << " " ;
        m_particle_class << endl ;
      }
    }
    if ( types_Global ) delete types_Global ;
    
    // Gather particle data ordered by particle ID on the master proc
    data_Global = wrapper->GatherParticleData_PostProcessing( 
    	  *particles, nb_total_part );

    // Write particle data
    if ( m_rank == 0 )
      if ( GrainsExec::m_ReloadType == "new" )
        one_output_MPI( time, nb_total_part, data_Global ) ;
    if ( data_Global ) delete data_Global ;
  }
  // In serial mode
  else
  {  
    // Extract the active particles that are not periodic clones
    // and create a map part ID - particle pointer such that we can write data
    // in increasing order of particle ID from 0 to nb_total_part-1
    map<size_t,Particle*> IDtoPart;
    list<Particle*>::const_iterator il;
    for (il=particles->begin(); il!=particles->end();il++)
      if ( (*il)->getTag() != 2 )
        IDtoPart.insert( pair<int,Particle*>( size_t((*il)->getID()), *il ) );
      
    map<size_t,Particle*>::const_iterator im; 
    if ( m_binary )
    {
      for (im=IDtoPart.begin();im!=IDtoPart.end();im++)
      {
	type = im->second->getGeometricType();
	m_particle_class.write( reinterpret_cast<char*>( &type ), sizeof(int) );
      }              
    }
    else
    { 
      for (im=IDtoPart.begin();im!=IDtoPart.end();im++)
        m_particle_class << im->second->getGeometricType() << " " ;
      m_particle_class << endl ;
    }

       
    // Write particle data
    one_output_Standard( time, nb_total_part, particles );    
  }
}




// ----------------------------------------------------------------------------
// Writes data
void RawDataPostProcessingWriter::PostProcessing( double const& time, 
    double const& dt,
    list<Particle*> const* particles,
    list<Particle*> const* inactiveparticles,
    list<Particle*> const* periodic_clones,
    vector<Particle*> const* referenceParticles,
    Obstacle* obstacle,
    LinkedCell const* LC )
{

  GrainsMPIWrapper const* wrapper = GrainsExec::getComm() ;
  size_t nb_total_part = GrainsExec::getTotalNumberPhysicalParticles() ;

  if ( wrapper )
  {
    vector< vector<double> >* data_Global = 
        wrapper->GatherParticleData_PostProcessing( *particles,
        nb_total_part );

    // Write data
    if ( m_rank == 0 )
      one_output_MPI( time, nb_total_part, data_Global ) ;

    if ( data_Global ) delete data_Global ;
  }
  else
    one_output_Standard( time, nb_total_part, particles );
}




// ----------------------------------------------------------------------------
// Finalizes writing data
void RawDataPostProcessingWriter::PostProcessing_end()
{
  if ( m_rank == 0 )
  {
    m_gc_coordinates_x.close();
    m_gc_coordinates_y.close();  
    m_gc_coordinates_z.close(); 
    m_translational_velocity_x.close();
    m_translational_velocity_y.close();   
    m_translational_velocity_z.close();
    m_angular_velocity_x.close();
    m_angular_velocity_y.close();   
    m_angular_velocity_z.close();
    m_coordination_number.close();
    m_particle_class.close();
  }
}




// ----------------------------------------------------------------------------
// Writes data in parallel mode at one physical time
void RawDataPostProcessingWriter::one_output_MPI( double const& time, 
    size_t& nb_total_part, 
    vector< vector<double> > const* data_Global )
{
  if ( m_binary )
  {
    double tt = time;
    m_gc_coordinates_x.write( reinterpret_cast<char*>( &tt ), 
    	sizeof(double) );
    m_gc_coordinates_y.write( reinterpret_cast<char*>( &tt ), 
    	sizeof(double) );
    m_gc_coordinates_z.write( reinterpret_cast<char*>( &tt ), 
    	sizeof(double) );
    m_translational_velocity_x.write( reinterpret_cast<char*>( &tt ), 
    	sizeof(double) );
    m_translational_velocity_y.write( reinterpret_cast<char*>( &tt ), 
    	sizeof(double) );
    m_translational_velocity_z.write( reinterpret_cast<char*>( &tt ), 
    	sizeof(double) );
    m_angular_velocity_x.write( reinterpret_cast<char*>( &tt ), 
    	sizeof(double) );
    m_angular_velocity_y.write( reinterpret_cast<char*>( &tt ), 
    	sizeof(double) );
    m_angular_velocity_z.write( reinterpret_cast<char*>( &tt ), 
    	sizeof(double) );
    m_angular_velocity_z.write( reinterpret_cast<char*>( &tt ), 
    	sizeof(double) );					  
  }
  else
  {
    string stime = GrainsExec::doubleToString( ios::scientific, 6, time );
    m_gc_coordinates_x << stime;
    m_gc_coordinates_y << stime;
    m_gc_coordinates_z << stime;
    m_translational_velocity_x << stime;
    m_translational_velocity_y << stime;
    m_translational_velocity_z << stime;
    m_angular_velocity_x << stime;
    m_angular_velocity_y << stime;
    m_angular_velocity_z << stime;
    m_angular_velocity_z << stime;
  }

  // Write in files: inactive particles are assigned 0 values and values are
  // written in increasing order of particle ID from 0 to nb_total_part-1
  if ( m_binary )
  {
    int coord = 0;
    double val = 0.;
    
    for (size_t i=0; i<nb_total_part; i++)
    {
      // Center of mass position
      val = (*data_Global)[i][0];
      m_gc_coordinates_x.write( reinterpret_cast<char*>( &val ), 
      	sizeof(double) ); 
      val = (*data_Global)[i][1];	
      m_gc_coordinates_y.write( reinterpret_cast<char*>( &val ), 
      	sizeof(double) );
      val = (*data_Global)[i][2];       
      m_gc_coordinates_z.write( reinterpret_cast<char*>( &val ), 
      	sizeof(double) ); 

      // Translational velocity
      val = (*data_Global)[i][3];      
      m_translational_velocity_x.write( reinterpret_cast<char*>( &val ), 
      	sizeof(double) );
      val = (*data_Global)[i][4];      
      m_translational_velocity_y.write( reinterpret_cast<char*>( &val ), 
      	sizeof(double) );	
      val = (*data_Global)[i][5];      
      m_translational_velocity_z.write( reinterpret_cast<char*>( &val ), 
      	sizeof(double) );	 
    
      // Angular velocity
      val = (*data_Global)[i][6];       
      m_angular_velocity_x.write( reinterpret_cast<char*>( &val ), 
      	sizeof(double) );
      val = (*data_Global)[i][7];       
      m_angular_velocity_y.write( reinterpret_cast<char*>( &val ), 
      	sizeof(double) );	
      val = (*data_Global)[i][8];       
      m_angular_velocity_z.write( reinterpret_cast<char*>( &val ), 
      	sizeof(double) );	

      // Number of contacts
      coord = int((*data_Global)[i][9]) ;
      m_coordination_number.write( reinterpret_cast<char*>( 
      	&coord ), sizeof(int) );
    }  
  }
  else
  { 
    for (size_t i=0; i<nb_total_part; i++)
    {
      // Center of mass position
      m_gc_coordinates_x << " " << GrainsExec::doubleToString( 
    	ios::scientific, m_ndigits, (*data_Global)[i][0] ) ;
      m_gc_coordinates_y << " " << GrainsExec::doubleToString( 
    	ios::scientific, m_ndigits, (*data_Global)[i][1] ) ;
      m_gc_coordinates_z << " " << GrainsExec::doubleToString( 
    	ios::scientific, m_ndigits, (*data_Global)[i][2] ) ;

      // Translational velocity
      m_translational_velocity_x << " " << GrainsExec::doubleToString( 
    	ios::scientific, m_ndigits, (*data_Global)[i][3] ) ;
      m_translational_velocity_y << " " << GrainsExec::doubleToString( 
    	ios::scientific, m_ndigits, (*data_Global)[i][4] ) ;
      m_translational_velocity_z << " " << GrainsExec::doubleToString( 
    	ios::scientific, m_ndigits, (*data_Global)[i][5] ) ;
    
      // Angular velocity
      m_angular_velocity_x << " " << GrainsExec::doubleToString( 
    	ios::scientific, m_ndigits, (*data_Global)[i][6] ) ;
      m_angular_velocity_y << " " << GrainsExec::doubleToString( 
    	ios::scientific, m_ndigits, (*data_Global)[i][7] ) ;
      m_angular_velocity_z << " " << GrainsExec::doubleToString( 
    	ios::scientific, m_ndigits, (*data_Global)[i][8] ) ;

      // Number of contacts
      m_coordination_number << " " << int((*data_Global)[i][9] ) ;
    }
  }
  
  if ( !m_binary )
  {
    m_gc_coordinates_x << endl ;
    m_gc_coordinates_y << endl ;
    m_gc_coordinates_z << endl ;
    m_translational_velocity_x << endl ;
    m_translational_velocity_y << endl ;
    m_translational_velocity_z << endl ;
    m_angular_velocity_x << endl ;
    m_angular_velocity_y << endl ;
    m_angular_velocity_z << endl ;
    m_coordination_number << endl ;
  }
}




// ----------------------------------------------------------------------------
// Delete all result files
void RawDataPostProcessingWriter::clearResultFiles() const
{
  if ( m_rank == 0 ) 
  {
    string cmd = "bash " + GrainsExec::m_GRAINS_HOME 
        + "/Tools/ExecScripts/Text_clear.exec " + m_filerootname;
    GrainsExec::m_return_syscmd = system( cmd.c_str() );
  }
}




// ----------------------------------------------------------------------------
// Creates output files and open streams
void RawDataPostProcessingWriter::prepareResultFiles( ios_base::openmode mode )
{
  string file;
  file = m_filerootname+"_position_x.dat";
  m_gc_coordinates_x.open( file.c_str(), mode | ( m_binary ? 
  	ios::binary : mode ) );
  file = m_filerootname+"_position_y.dat";
  m_gc_coordinates_y.open( file.c_str(), mode | ( m_binary ? 
  	ios::binary : mode ) );
  file = m_filerootname+"_position_z.dat";
  m_gc_coordinates_z.open( file.c_str(), mode | ( m_binary ? 
  	ios::binary : mode ) );

  file = m_filerootname+"_translational_velocity_x.dat";
  m_translational_velocity_x.open( file.c_str(), mode | ( m_binary ? 
  	ios::binary : mode ) );
  file = m_filerootname+"_translational_velocity_y.dat";
  m_translational_velocity_y.open( file.c_str(), mode | ( m_binary ? 
  	ios::binary : mode ) );
  file = m_filerootname+"_translational_velocity_z.dat";
  m_translational_velocity_z.open( file.c_str(), mode | ( m_binary ? 
  	ios::binary : mode ) ); 

  file = m_filerootname+"_angular_velocity_x.dat";
  m_angular_velocity_x.open( file.c_str(), mode | ( m_binary ? 
  	ios::binary : mode ) );
  file = m_filerootname+"_angular_velocity_y.dat";
  m_angular_velocity_y.open( file.c_str(), mode | ( m_binary ? 
  	ios::binary : mode ) );
  file = m_filerootname+"_angular_velocity_z.dat";
  m_angular_velocity_z.open( file.c_str(), mode | ( m_binary ? 
  	ios::binary : mode ) );
  
  file = m_filerootname+"_coordinationNumber.dat";
  m_coordination_number.open( file.c_str(), mode | ( m_binary ? 
  	ios::binary : mode ) );

  file = m_filerootname+"_particleType.dat";
  m_particle_class.open( file.c_str(), mode | ( m_binary ? 
  	ios::binary : mode ) );

}




// ----------------------------------------------------------------------------
// Writes data in serial mode at one physical time
void RawDataPostProcessingWriter::one_output_Standard( double const& time,
    size_t& nb_total_part, list<Particle*> const* particles )
{ 
  Point3* centre = NULL;
  Vector3* velT = NULL; 
  Vector3* velR = NULL; 
  double zero = 0., tt = time; 
  int coord = 0, izero = 0; 

  if ( m_binary )
  {
    m_gc_coordinates_x.write( reinterpret_cast<char*>( &tt ), 
    	sizeof(double) );
    m_gc_coordinates_y.write( reinterpret_cast<char*>( &tt ), 
    	sizeof(double) );
    m_gc_coordinates_z.write( reinterpret_cast<char*>( &tt ), 
    	sizeof(double) );
    m_translational_velocity_x.write( reinterpret_cast<char*>( &tt ), 
    	sizeof(double) );
    m_translational_velocity_y.write( reinterpret_cast<char*>( &tt ), 
    	sizeof(double) );
    m_translational_velocity_z.write( reinterpret_cast<char*>( &tt ), 
    	sizeof(double) );
    m_angular_velocity_x.write( reinterpret_cast<char*>( &tt ), 
    	sizeof(double) );
    m_angular_velocity_y.write( reinterpret_cast<char*>( &tt ), 
    	sizeof(double) );
    m_angular_velocity_z.write( reinterpret_cast<char*>( &tt ), 
    	sizeof(double) );
    m_angular_velocity_z.write( reinterpret_cast<char*>( &tt ), 
    	sizeof(double) );					  
  }
  else
  {
    string stime = GrainsExec::doubleToString( ios::scientific, 6, time );
    m_gc_coordinates_x << stime;
    m_gc_coordinates_y << stime;
    m_gc_coordinates_z << stime;
    m_translational_velocity_x << stime;
    m_translational_velocity_y << stime;
    m_translational_velocity_z << stime;
    m_angular_velocity_x << stime;
    m_angular_velocity_y << stime;
    m_angular_velocity_z << stime;
    m_angular_velocity_z << stime;
  }
  
  // Extract the active particles that are not periodic clones
  // and create a map part ID - particle pointer such that we can write data
  // in increasing order of particle ID from 0 to nb_total_part-1
  map<size_t,Particle*> IDtoPart;
  list<Particle*>::const_iterator il;
  Particle* pp = NULL;
  for (il=particles->begin(); il!=particles->end();il++)
    if ( (*il)->getTag() != 2 )
       IDtoPart.insert( pair<int,Particle*>( size_t((*il)->getID()), *il ) ); 
       
  // Write in files: inactive particles are assigned 0 values and values are
  // written in increasing order of particle ID from 0 to nb_total_part-1
  for (size_t i=0;i<nb_total_part;i++)
    if ( IDtoPart.count( i ) )
    {
      pp = IDtoPart[i];
      
      if ( m_binary )
      {
        // Center of mass position
        centre = const_cast<Point3*>(pp->getPosition());
	m_gc_coordinates_x.write( reinterpret_cast<char*>( &(*centre)[X] ), 
    		sizeof(double) );
	m_gc_coordinates_y.write( reinterpret_cast<char*>( &(*centre)[Y] ), 
    		sizeof(double) );
	m_gc_coordinates_z.write( reinterpret_cast<char*>( &(*centre)[Z] ), 
    		sizeof(double) );
		
        // Translational velocity
        velT = const_cast<Vector3*>(pp->getTranslationalVelocity());
        m_translational_velocity_x.write( reinterpret_cast<char*>( 
		&(*velT)[X] ), sizeof(double) );
        m_translational_velocity_y.write( reinterpret_cast<char*>( 
		&(*velT)[Y] ), sizeof(double) );
        m_translational_velocity_z.write( reinterpret_cast<char*>( 
		&(*velT)[Z] ), sizeof(double) );
		    
        // Angular velocity
        velR = const_cast<Vector3*>(pp->getAngularVelocity());
        m_angular_velocity_x.write( reinterpret_cast<char*>( 
		&(*velR)[X] ), sizeof(double) );
        m_angular_velocity_y.write( reinterpret_cast<char*>( 
		&(*velR)[Y] ), sizeof(double) );
        m_angular_velocity_z.write( reinterpret_cast<char*>( 
		&(*velR)[Z] ), sizeof(double) );
    
        // Number of contacts
        coord = pp->getCoordinationNumber(); 
	m_coordination_number.write( reinterpret_cast<char*>( 
		&coord ), sizeof(int) );
      }
      else
      {
        // Center of mass position
        centre = const_cast<Point3*>(pp->getPosition());
        m_gc_coordinates_x << " " << GrainsExec::doubleToString( 
    	ios::scientific, m_ndigits, (*centre)[X] ) ;
        m_gc_coordinates_y << " " << GrainsExec::doubleToString( 
    	ios::scientific, m_ndigits, (*centre)[Y] ) ;
        m_gc_coordinates_z << " " << GrainsExec::doubleToString( 
    	ios::scientific, m_ndigits, (*centre)[Z] ) ;

        // Translational velocity
        velT = const_cast<Vector3*>(pp->getTranslationalVelocity());
        m_translational_velocity_x << " " << GrainsExec::doubleToString( 
    	ios::scientific, m_ndigits, (*velT)[X] ) ;
        m_translational_velocity_y << " " << GrainsExec::doubleToString( 
    	ios::scientific, m_ndigits, (*velT)[Y] ) ;
        m_translational_velocity_z << " " << GrainsExec::doubleToString( 
    	ios::scientific, m_ndigits, (*velT)[Z] ) ;
    
        // Angular velocity
        velR = const_cast<Vector3*>(pp->getAngularVelocity());
        m_angular_velocity_x << " " << GrainsExec::doubleToString( 
    	ios::scientific, m_ndigits, (*velR)[X] ) ;
        m_angular_velocity_y << " " << GrainsExec::doubleToString( 
    	ios::scientific, m_ndigits, (*velR)[Y] ) ;
        m_angular_velocity_z << " " << GrainsExec::doubleToString( 
    	ios::scientific, m_ndigits, (*velR)[Z] ) ;
    
        // Number of contacts
        m_coordination_number << " " << pp->getCoordinationNumber(); 
      }         
    }
    else
    {
      if ( m_binary )
      {
        // Center of mass position
	m_gc_coordinates_x.write( reinterpret_cast<char*>( 
		&GrainsExec::m_defaultInactivePos[X] ), 
    		sizeof(double) );
	m_gc_coordinates_y.write( reinterpret_cast<char*>( 
		&GrainsExec::m_defaultInactivePos[Y] ), 
    		sizeof(double) );
	m_gc_coordinates_z.write( reinterpret_cast<char*>( &
		GrainsExec::m_defaultInactivePos[Z] ), 
    		sizeof(double) );
		
        // Translational velocity
        m_translational_velocity_x.write( reinterpret_cast<char*>( 
		&zero ), sizeof(double) );
        m_translational_velocity_y.write( reinterpret_cast<char*>( 
		&zero ), sizeof(double) );
        m_translational_velocity_z.write( reinterpret_cast<char*>( 
		&zero ), sizeof(double) );
		    
        // Angular velocity
        m_angular_velocity_x.write( reinterpret_cast<char*>( 
		&zero ), sizeof(double) );
        m_angular_velocity_y.write( reinterpret_cast<char*>( 
		&zero ), sizeof(double) );
        m_angular_velocity_z.write( reinterpret_cast<char*>( 
		&zero ), sizeof(double) );
    
        // Number of contacts
	m_coordination_number.write( reinterpret_cast<char*>( 
		&izero ), sizeof(int) );      
      }
      else
      {      
        // Center of mass position
        m_gc_coordinates_x << " " << GrainsExec::doubleToString( 
    	ios::scientific, m_ndigits, GrainsExec::m_defaultInactivePos[X] ) ;
        m_gc_coordinates_y << " " << GrainsExec::doubleToString( 
    	ios::scientific, m_ndigits, GrainsExec::m_defaultInactivePos[Y] ) ;
        m_gc_coordinates_z << " " << GrainsExec::doubleToString( 
    	ios::scientific, m_ndigits, GrainsExec::m_defaultInactivePos[Z] ) ;

        // Translational velocity
        m_translational_velocity_x << " " << GrainsExec::doubleToString( 
    	ios::scientific, m_ndigits, zero ) ;
        m_translational_velocity_y << " " << GrainsExec::doubleToString( 
    	ios::scientific, m_ndigits, zero ) ;
        m_translational_velocity_z << " " << GrainsExec::doubleToString( 
    	ios::scientific, m_ndigits, zero ) ;
    
        // Angular velocity
        m_angular_velocity_x << " " << GrainsExec::doubleToString( 
    	ios::scientific, m_ndigits, zero ) ;
        m_angular_velocity_y << " " << GrainsExec::doubleToString( 
    	ios::scientific, m_ndigits, zero ) ;
        m_angular_velocity_z << " " << GrainsExec::doubleToString( 
    	ios::scientific, m_ndigits, zero ) ;
    
        // Number of contacts
        m_coordination_number << " 0";
      }    
    }   
  
  m_gc_coordinates_x << endl;
  m_gc_coordinates_y << endl;  
  m_gc_coordinates_z << endl; 
  m_translational_velocity_x << endl;
  m_translational_velocity_y << endl;   
  m_translational_velocity_z << endl;
  m_angular_velocity_x << endl;
  m_angular_velocity_y << endl;   
  m_angular_velocity_z << endl;  
  m_coordination_number << endl;
}




// ----------------------------------------------------------------------------
// Gets the post-processing writer type
string RawDataPostProcessingWriter::getPostProcessingWriterType() const
{
  return ( "RawData" );
}
