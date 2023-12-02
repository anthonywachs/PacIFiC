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
    int const& rank_, int const& nbranks_ )
  : PostProcessingWriter( dn, rank_, nbranks_ )
{ 
  m_filerootname = ReaderXML::getNodeAttr_String( dn, "Name" );
 
  if ( m_rank == 0 )
  {
    cout << GrainsExec::m_shift9 << "Type = Text" << endl;
    cout << GrainsExec::m_shift12 << "Output file name = " 
    	<< m_filerootname << endl;
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
    list<Particle*> const* pwait,
    list<Particle*> const* periodic_clones,
    vector<Particle*> const* referenceParticles,
    Obstacle* obstacle,
    LinkedCell const* LC,
    vector<Window> const& insert_windows )
{
  GrainsMPIWrapper const* wrapper = GrainsExec::getComm() ;
  size_t nb_total_part = GrainsExec::getNumberParticlesOnAllProc() ;

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

  // IF PARALLEL
  if( wrapper )
  {
    vector< vector<double> >* cinematique_Global;

    // Réunir la cinematique des particles de tous les procs sur le master
    cinematique_Global = wrapper->GatherPositionVelocity_PostProcessing( 
    	  *particles, nb_total_part );

    // Ecrire les resultats contenus dans cinematique_Global
    if ( m_rank == 0 )
      if ( GrainsExec::m_ReloadType == "new" )
        one_output_MPI( time, nb_total_part, cinematique_Global ) ;

    // Gather particles class from every proc on master proc
    vector< vector<double> >* types_Global = 
        wrapper->GatherParticlesClass_PostProcessing( *particles,
            nb_total_part );

    // Write down particles class only once at the begining
    if ( m_rank == 0 )
    {
      for (size_t i=0; i<nb_total_part; i++)
        m_particle_class << (*types_Global)[0][i] << " " ;
      m_particle_class << endl ;
    }

    if ( cinematique_Global )
      delete cinematique_Global ;

    if ( types_Global )
      delete types_Global ;
  }
  else
    one_output_Standard( time, particles, pwait );
}




// ----------------------------------------------------------------------------
// Writes data
void RawDataPostProcessingWriter::PostProcessing( double const& time, 
    double const& dt,
    list<Particle*> const* particles,
    list<Particle*> const* pwait,
    list<Particle*> const* periodic_clones,
    vector<Particle*> const* referenceParticles,
    Obstacle* obstacle,
    LinkedCell const* LC )
{

  GrainsMPIWrapper const* wrapper = GrainsExec::getComm() ;
  size_t nb_total_part = GrainsExec::getNumberParticlesOnAllProc() ;

  if( wrapper )
  {
    vector< vector<double> >* cinematique_Global = 
        wrapper->GatherPositionVelocity_PostProcessing( *particles,
        nb_total_part );

    // Ecrire les résultats contenus dans cinematique_Global
    if ( m_rank == 0 )
      one_output_MPI( time, nb_total_part, cinematique_Global ) ;

    if ( cinematique_Global )
      delete cinematique_Global ;
  }
  else
    one_output_Standard( time, particles, pwait );
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
void RawDataPostProcessingWriter::one_output_MPI(double const& time, 
    size_t& nb_total_part, 
    vector< vector<double> > const* cinematique_Global)
{
  vector<double> InternalFeatures(9,0.);

  m_gc_coordinates_x << time;
  m_gc_coordinates_y << time;
  m_gc_coordinates_z << time;
  m_translational_velocity_x << time;
  m_translational_velocity_y << time;
  m_translational_velocity_z << time;
  m_angular_velocity_x << time;
  m_angular_velocity_y << time;
  m_angular_velocity_z << time;
  m_coordination_number << time;

  // Dans le cas d'une insertion, tant que la particle n'est pas insérée
  // la vitesse & la position sont nulles
  for (size_t i=0; i<nb_total_part; i++)
  {
    // Position du centre de gravité
    m_gc_coordinates_x << " " << (*cinematique_Global)[0][i] ;
    m_gc_coordinates_y << " " << (*cinematique_Global)[1][i] ;
    m_gc_coordinates_z << " " << (*cinematique_Global)[2][i] ;

    // Velocity translationnelle du centre de gravité
    m_translational_velocity_x << " " << (*cinematique_Global)[3][i] ;
    m_translational_velocity_y << " " << (*cinematique_Global)[4][i] ;
    m_translational_velocity_z << " " << (*cinematique_Global)[5][i] ;
    
    // Velocity de rotation du centre de gravité
    m_angular_velocity_x << " " << (*cinematique_Global)[6][i] ;
    m_angular_velocity_y << " " << (*cinematique_Global)[7][i] ;
    m_angular_velocity_z << " " << (*cinematique_Global)[8][i] ;

    // Coordination number
    m_coordination_number << " " << (*cinematique_Global)[9][i] ;
  }
  
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
  m_gc_coordinates_x.open(file.c_str(),mode);
  file = m_filerootname+"_position_y.dat";
  m_gc_coordinates_y.open(file.c_str(),mode);
  file = m_filerootname+"_position_z.dat";
  m_gc_coordinates_z.open(file.c_str(),mode);

  file = m_filerootname+"_translational_velocity_x.dat";
  m_translational_velocity_x.open(file.c_str(),mode);
  file = m_filerootname+"_translational_velocity_y.dat";
  m_translational_velocity_y.open(file.c_str(),mode);
  file = m_filerootname+"_translational_velocity_z.dat";
  m_translational_velocity_z.open(file.c_str(),mode); 

  file = m_filerootname+"_angular_velocity_x.dat";
  m_angular_velocity_x.open(file.c_str(),mode);
  file = m_filerootname+"_angular_velocity_y.dat";
  m_angular_velocity_y.open(file.c_str(),mode);
  file = m_filerootname+"_angular_velocity_z.dat";
  m_angular_velocity_z.open(file.c_str(),mode);
  
  file = m_filerootname+"_coordinationNumber.dat";
  m_coordination_number.open( file.c_str(), mode );

  file = m_filerootname+"_particleType.dat";
  m_particle_class.open( file.c_str(), mode );

}




// ----------------------------------------------------------------------------
// Writes data in serial mode at one physical time
void RawDataPostProcessingWriter::one_output_Standard(double const& time,
    list<Particle*> const* particles,
    list<Particle*> const* pwait)
{ 
  Point3 const* centre = NULL;
  Vector3 const* velT = NULL; 
  Vector3 const* velR = NULL;    

  m_gc_coordinates_x << time;
  m_gc_coordinates_y << time;  
  m_gc_coordinates_z << time; 
  m_translational_velocity_x << time;
  m_translational_velocity_y << time;   
  m_translational_velocity_z << time;
  m_angular_velocity_x << time;
  m_angular_velocity_y << time;   
  m_angular_velocity_z << time;
  m_coordination_number << time;

  list<Particle*>::const_iterator particle;
  for (particle=particles->begin(); particle!=particles->end();particle++)
    if ( (*particle)->getTag() != 2 )
    {
      // Center of mass position
      centre = (*particle)->getPosition();
      m_gc_coordinates_x << " " << (*centre)[X];
      m_gc_coordinates_y << " " << (*centre)[Y];
      m_gc_coordinates_z << " " << (*centre)[Z];

      // Translational velocity
      velT = (*particle)->getTranslationalVelocity();
      m_translational_velocity_x << " " << (*velT)[X];
      m_translational_velocity_y << " " << (*velT)[Y];
      m_translational_velocity_z << " " << (*velT)[Z]; 
    
      // Angular velocity
      velR = (*particle)->getAngularVelocity();
      m_angular_velocity_x << " " << (*velR)[X];
      m_angular_velocity_y << " " << (*velR)[Y];
      m_angular_velocity_z << " " << (*velR)[Z];
    
      // Number of contacts
      m_coordination_number << " " << (*particle)->getCoordinationNumber();
    }

  if( pwait )
    for (particle=pwait->begin(); particle!=pwait->end();particle++)
    {
      // Center of mass position
      centre = (*particle)->getPosition();
      m_gc_coordinates_x << " " << (*centre)[X];
      m_gc_coordinates_y << " " << (*centre)[Y];
      m_gc_coordinates_z << " " << (*centre)[Z];

      // Translational velocity
      velT = (*particle)->getTranslationalVelocity();
      m_translational_velocity_x << " " << (*velT)[X];
      m_translational_velocity_y << " " << (*velT)[Y];
      m_translational_velocity_z << " " << (*velT)[Z];

      // Angular velocity
      velR = (*particle)->getAngularVelocity();
      m_angular_velocity_x << " " << (*velR)[X];
      m_angular_velocity_y << " " << (*velR)[Y];
      m_angular_velocity_z << " " << (*velR)[Z];
      
      // Number of contacts
      m_coordination_number << " 0";
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
  return ( "Text" );
}
