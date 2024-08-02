#include "GrainsMPI.hh"
#include "ContactBuilderFactory.hh"
#include "LinkedCell.hh"
#include "ObstacleBuilderFactory.hh"
#include "ObstacleImposedVelocity.hh"
#include "GrainsBuilderFactory.hh"
#include "PostProcessingWriter.hh"
#include "PostProcessingWriterBuilderFactory.hh"
#include "RawDataPostProcessingWriter.hh"
#include "stdlib.h"


// ----------------------------------------------------------------------------
// Default constructor
GrainsMPI::GrainsMPI()
  : Grains()
{}




// ----------------------------------------------------------------------------
// Destructor
GrainsMPI::~GrainsMPI()
{
  delete m_wrapper; 
}




// ----------------------------------------------------------------------------
// Writes an initial message in the standard output only on the process ranked 0
void GrainsMPI::initialOutputMessage()
{
  int rankproc = 0;
  MPI_Comm_rank( MPI_COMM_WORLD, &rankproc );
  if ( rankproc == 0 )
    cout << "Grains3D MPI/Static uniform domain decomposition" << endl;
}




// ----------------------------------------------------------------------------
// Runs the simulation over the prescribed time interval
void GrainsMPI::Simulation( double time_interval )
{
  if ( m_processorIsActive )
  {
    double vmax = 0., vmean = 0. ;
    list<Particle*>* newBufPart = new list<Particle*>;
    bool forcestats = false;    

    // Timers
    SCT_insert_app( "ParticlesInsertion" );
    SCT_insert_app( "ComputeForces" );
    SCT_insert_app( "Move" );
    SCT_insert_app( "UpdateParticleActivity" );
    SCT_insert_app( "LinkUpdate" );
    SCT_insert_app( "OutputResults" );

    // Simulation: time marching algorithm
    while ( m_tend - m_time > 0.01 * m_dt )
    {
      try
      {
        m_time += m_dt;
	

        // Check whether data are output at this time
        m_lastTime_save = false;
        GrainsExec::m_output_data_at_this_time = false;
	GrainsExec::m_postprocess_forces_at_this_time = false;
        if ( m_timeSave != m_save.end() )
          if ( *m_timeSave - m_time < 0.01 * m_dt )
	  {
	    // Set the global data output boolean to true
	    GrainsExec::m_output_data_at_this_time = true;
	    if ( m_rank == 0 ) cout << endl << "Time = " << m_time << endl
	  	<< std::flush;

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
        if ( m_insertion_mode == IM_OVERTIME )
          insertParticle( m_insertion_order );	
        m_allcomponents.computeNumberParticles( m_wrapper );
	if ( GrainsExec::m_output_data_at_this_time )
          if ( m_rank == 0 )
	  {
	    if ( m_allcomponents.getNumberPhysicalParticlesToInsert() ) 
	    {
	      cout << "Number of active particles on all proc = " <<
	  	m_allcomponents.getNumberActiveParticlesOnAllProc() << endl;
	      cout << "Number of particles to insert = " <<
	  	m_allcomponents.getNumberPhysicalParticlesToInsert() 
		<< endl;
	    }		
	  }
        SCT_get_elapsed_time( "ParticlesInsertion" );

      
        // Move particles and obstacles
        // Update particle velocity over dt/2 and particle position over dt,
        // obstacle velocity and position over dt
        // v_i+1/2 = v_i + a_i * dt / 2
        // x_i+1 = x_i + v_i+1/2 * dt
        // Solve Newton's law and move particles
        SCT_set_start( "Move" );
        m_allcomponents.Move( m_time, 0.5 * m_dt, m_dt, m_dt, m_collision );

   
        // Update clone particle position and velocity
        m_wrapper->UpdateOrCreateClones_SendRecvLocal_GeoLoc( m_time,
		m_allcomponents.getActiveParticles(),
  		m_allcomponents.getParticlesInBufferzone(),
  		m_allcomponents.getCloneParticles(),
		m_allcomponents.getReferenceParticles(),
		m_collision, true, false );


        // Destroy out of domain clones
        m_collision->DestroyOutOfDomainClones( m_time,
		m_allcomponents.getActiveParticles(),
		m_allcomponents.getCloneParticles(),
		m_wrapper );      
        SCT_get_elapsed_time( "Move" );
      
      
        // Update particle activity
        SCT_set_start( "UpdateParticleActivity" );
        m_allcomponents.UpdateParticleActivity();
        SCT_get_elapsed_time( "UpdateParticleActivity" );


        // Update the particle & obstacles links with the grid
        SCT_set_start( "LinkUpdate" );
        m_collision->LinkUpdate( m_time, m_dt, 
      		m_allcomponents.getActiveParticles() );  
        m_allcomponents.updateParticleLists( m_time, newBufPart ); 


        // Create new clones
        m_wrapper->UpdateOrCreateClones_SendRecvLocal_GeoLoc( m_time,
		m_allcomponents.getActiveParticles(),
  		newBufPart,
  		m_allcomponents.getCloneParticles(),
		m_allcomponents.getReferenceParticles(),
		m_collision, false, false );
        SCT_get_elapsed_time( "LinkUpdate" );
      
      
        // Compute particle forces and acceleration
        // Compute f_i+1 and a_i+1 as a function of (x_i+1,v_i+1/2)
        SCT_set_start( "ComputeForces" );
        computeParticlesForceAndAcceleration();
        SCT_get_elapsed_time( "ComputeForces" );


        // Update particle velocity over dt/2
        // v_i+1 = v_i+1/2 + a_i+1 * dt / 2 
        SCT_set_start( "Move" );      
        m_allcomponents.advanceParticlesVelocity( m_time, 0.5 * m_dt );
        SCT_add_elapsed_time( "Move" );


        // Compute and write force & torque exerted on obstacles
        m_allcomponents.computeObstaclesLoad( m_time, m_dt, m_wrapper ); 
	m_allcomponents.outputObstaclesLoad( m_time, m_dt, false,
      		false, m_rank );
		

        // Compute and write force statistics
        if ( forcestats ) m_collision->outputForceStats( m_time, m_dt, m_rank, 
      		m_wrapper );		


        // Write postprocessing and reload files
        if ( GrainsExec::m_output_data_at_this_time )
        {
	  SCT_set_start( "OutputResults" );

          // Update clone particle velocity
          m_wrapper->UpdateOrCreateClones_SendRecvLocal_GeoLoc( m_time,
		m_allcomponents.getActiveParticles(),
  		m_allcomponents.getParticlesInBufferzone(),
  		m_allcomponents.getCloneParticles(),
		m_allcomponents.getReferenceParticles(),
		m_collision, true, true );
		
	  // Write time, track component max and mean velocity
	  m_allcomponents.ComputeMaxMeanVelocity( vmax, vmean, m_wrapper );
          if ( m_rank == 0 )
	  {
	    cout << "Component velocity : max = " << vmax
               << " average = " << vmean << endl;
            fVitMax << GrainsExec::doubleToString( ios::scientific, 6, m_time )
               << "\t" << GrainsExec::doubleToString( ios::scientific, 6,
               vmax ) << "\t" << GrainsExec::doubleToString( ios::scientific,
               6, vmean ) << endl;
	  }

	  // Display memory used by Grains
	  display_used_memory();
		
          // Summary of MPI particle comms
          if ( GrainsExec::m_MPI_verbose ) 
            m_wrapper->writeAndFlushMPIString( cout );		

	  // Write postprocessing files
          m_allcomponents.PostProcessing( m_time, m_dt, m_collision,
		m_rank, m_nprocs, m_wrapper );

	  // Write reload files
	  saveReload( m_time );

	  SCT_get_elapsed_time( "OutputResults" );
        }
      }
      catch (ContactError &errContact)
      {
        // Max overlap exceeded
        cout << endl;
        m_allcomponents.PostProcessingErreurComponents( "ContactError",
            errContact.getComponents() );
        errContact.Message( cout );
        m_error_occured = true;
        break;
      }
      catch (MotionError &errMotion)
      {
        // Particle motion over dt is too large
        cout << endl;
        m_allcomponents.PostProcessingErreurComponents( "MotionError",
            errMotion.getComponent() );
        errMotion.Message(cout);
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
}




// ----------------------------------------------------------------------------
// Reads data for MPI simulations and creates and sets the MPI wrapper
void GrainsMPI::readDomainDecomposition( DOMNode* root,
  	double const& lx, double const& ly, double const& lz )
{
  // Domain decomposition
  int nx, ny, nz = 1;
  DOMNode* decomp = ReaderXML::getNode( root, "DomainDecomposition" );
  nx = ReaderXML::getNodeAttr_Int( decomp, "NX" );
  ny = ReaderXML::getNodeAttr_Int( decomp, "NY" );
  if ( m_dimension == 3 ) nz = ReaderXML::getNodeAttr_Int( decomp, "NZ" ); 

  // MPI wrapper
  m_wrapper = new GrainsMPIWrapper( nx, ny, nz, m_periodicity[X], 
  	m_periodicity[Y], m_periodicity[Z], GrainsExec::m_shift3 ); 
  GrainsExec::setComm( m_wrapper );
  GrainsExec::m_MPI = true;
  m_processorIsActive = m_wrapper->isActive();
  m_rank = m_wrapper->get_rank();
  m_nprocs = m_wrapper->get_total_number_of_active_processes();
  DOMNode* mpiNode = ReaderXML::getNode( root, "Verbosity" );
  GrainsExec::m_MPI_verbose = 
  	ReaderXML::getNodeAttr_Int( mpiNode, "Level" );

  if ( m_processorIsActive )
  {      
    // Local domain geometric features
    App::set_local_domain_size( lx / m_wrapper->get_nb_procs_direction(0), 
  	ly / m_wrapper->get_nb_procs_direction(1), 
	lz / m_wrapper->get_nb_procs_direction(2) );
    App::set_local_domain_origin( m_wrapper->get_nb_procs_direction(),
  	m_wrapper->get_MPI_coordinates() );

    // MPI periodes
    if ( m_periodic )
      m_wrapper->setMPIperiodicVectors( lx, ly, lz );
    
    // Display the wrapper features
    m_wrapper->display( cout, GrainsExec::m_shift3 );	
  } 
}




// ----------------------------------------------------------------------------
// Attempts to insert a particle in the simulation
bool GrainsMPI::insertParticle( PullMode const& mode )
{
  static size_t insert_counter = 0;
  pair<bool,bool> insert( false, false );
  Vector3 vtrans, vrot ;
  Point3 position;
  Matrix mrot;
  Transform trot;
  Quaternion qrot;    
  int ptype = - 1;
  size_t npositions = m_insertion_position->size();
  Particle *particle = NULL;    
  
  if ( insert_counter == 0 )
  {
    // Return the particle to be inserted
    if ( m_rank == 0 ) 
    {
      particle = m_allcomponents.getParticleToInsert( mode );
      if ( particle ) ptype = particle->getGeometricType();
    }
    
    // Broadcast particle type
    ptype = m_wrapper->Broadcast_INT( ptype );
      
    if ( ptype > - 1 )
    {
      // Return the particle to be inserted
      particle = m_allcomponents.getParticleToInsert( ptype );
	
      // Initialisation of the centre of mass position of the particle
      if ( m_rank == 0 ) position = getInsertionPoint();
      position = m_wrapper->Broadcast_Point3( position );
      particle->setPosition( position );
      
      // Initialisation of the angular position of the particle
      // Rem: we compose to the right by a pure rotation as the particle
      // already has a non-zero position that we do not want to change (and
      // that we would change if we would compose to the left)
      if ( m_init_angpos == IAP_RANDOM )
      {
        if ( m_rank == 0 ) mrot = GrainsExec::RandomRotationMatrix( 
	  	m_dimension );
        mrot = m_wrapper->Broadcast_Matrix( mrot );
        trot.setBasis( mrot );
        particle->composePositionRightByTransform( trot );
      }            
      
      // Initialisation of the particle velocity
      if ( m_rank == 0 ) computeInitVit( vtrans, vrot );
      vtrans = m_wrapper->Broadcast_Vector3( vtrans ); 
      vrot = m_wrapper->Broadcast_Vector3( vrot ); 
      particle->setTranslationalVelocity( vtrans );
      particle->setAngularVelocity( vrot );      

      // Transform and quaternion
      particle->initialize_transformWithCrust_to_notComputed();
      qrot.setQuaternion( particle->getRigidBody()->getTransform()
		->getBasis() );
      particle->setQuaternionRotation( qrot );      

      // If insertion if successful, shift particle from wait to inserted
      // and initialize particle rotation quaternion from rotation matrix
      insert = m_collision->insertParticleParallel( m_time, particle,
      	m_allcomponents.getActiveParticles(),
	m_allcomponents.getCloneParticles(),
	m_allcomponents.getReferenceParticles(), 
	m_periodic, m_force_insertion,
	m_wrapper );

      // If no contact
      if ( !insert.second )
      {
        // If particle is in LinkedCell
	if ( insert.first )
	{
	  m_allcomponents.WaitToActive( true );
	  particle->InitializeForce( true );
	  particle->computeAcceleration( m_time );
        }
	// If not in LinkedCell
        else
          m_allcomponents.DeleteAndDestroyWait();

	if ( npositions ) m_insertion_position->erase( il_sp );
      }
    }
  }

  ++insert_counter;
  if ( insert_counter == m_insertion_frequency ) insert_counter = 0;

  return ( !insert.second && m_wrapper->max_UNSIGNED_INT( insert.first ) );
}




// ----------------------------------------------------------------------------
// Creates, inserts and links new particles in the simulation
void GrainsMPI::InsertCreateNewParticles()
{
  // IMPORTANT: for any system with N particles and M obstacles, particles are
  // numbered 1 to N and obstacles are numbered -1 to -M

  // Link all components with the grid
  m_allcomponents.Link( *m_collision );

  // In case of a restarted simulation, if the linked cell changed from the 
  // previous simulation, we need to check that all clones are there
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
  if ( error )  grainsAbort();

  // Insertion at initial time
  size_t nbPW = 0 ;
  bool b_insertion_BEFORE = true;  
  if ( m_insertion_mode == IM_INITIALTIME )
  {
    nbPW = m_allcomponents.getNumberPhysicalParticlesToInsert() ;
    for (size_t i=0;i<nbPW && b_insertion_BEFORE;++i)
      b_insertion_BEFORE = insertParticle( m_insertion_order );
    if ( !b_insertion_BEFORE )
    {
      cout << "Process " << m_rank << endl;
      cout << "Insertion issue with defined position method" << endl;
      cout << "Remaining number of particles to insert = " << 
      	m_allcomponents.getNumberPhysicalParticlesToInsert() << endl;
    }     
  }
  if ( !b_insertion_BEFORE ) grainsAbort();
  
  double volIN = m_allcomponents.getVolumeIn(), 
  	volOUT = m_allcomponents.getVolumeOut();
  volIN = m_wrapper->sum_DOUBLE( volIN );
  volOUT = m_wrapper->sum_DOUBLE( volOUT ); 
  
  if ( m_rank == 0 ) 
    cout << endl << "Volume des Particles IN  : " << volIN  << '\n'
      	<< "                     OUT : " << volOUT << '\n'
      	<< endl; 
}




// ----------------------------------------------------------------------------
// Sets particle initial positions from a file
size_t GrainsMPI::setPositionParticlesFromFile()
{
  size_t error = 0;
  
  // Note: only the master proc reads and stores positions  
  if ( m_rank == 0 ) error = GrainsMPI::setPositionParticlesFromFile();
  error = m_wrapper->Broadcast_UNSIGNED_INT( error );

  return ( error );
}




// ----------------------------------------------------------------------------
// Sets particle initial position with a structured array
size_t GrainsMPI::setPositionParticlesArray()
{
  size_t error = 0;
  
  // Checks that none of the structured array positions is exactly 
  // at a limit of the linked cell grid, otherwise shift by 1e-12
  m_collision->checkStructuredArrayPositionsMPI( m_InsertionArray, m_wrapper );
  
  // Note: only the master proc computes and stores positions  
  if ( m_rank == 0 ) error = Grains::setPositionParticlesArray();
  error = m_wrapper->Broadcast_UNSIGNED_INT( error );

  return ( error );    
}




// ----------------------------------------------------------------------------
// Sets the linked cell grid
void GrainsMPI::defineLinkedCell( double const& radius, string const& oshift )
{
  size_t error = m_collision->set( 2. * radius, 
  	m_wrapper->get_nb_procs_direction(),
  	m_wrapper->get_MPI_coordinates(), m_wrapper->get_MPI_neighbors(),
	oshift );
  if ( error ) grainsAbort();
}




// ----------------------------------------------------------------------------
// Emergency termination in case of an issue
void GrainsMPI::grainsAbort() const
{
  int error_code = 0;
  MPI_Abort( MPI_COMM_WORLD, error_code );
}




// ----------------------------------------------------------------------------
// Returns the maximum particle ID number
int GrainsMPI::getMaxParticleIDnumber() const
{
  int numMax = m_allcomponents.getMaxParticleIDnumber();
  int collective_numMax = m_wrapper->max_INT( numMax );

  return ( collective_numMax );
}




// ----------------------------------------------------------------------------
// Displays the memory used by the simulation
void GrainsMPI::display_used_memory() const
{
  m_wrapper->display_used_memory( cout );
}




// ----------------------------------------------------------------------------
// Synchronizes the PPWindow boolean relative to each sub-domain
void GrainsMPI::synchronize_PPWindow()
{
  vector<bool> b;
  b = PostProcessingWriter::get_PostProcessingWindow();
  
  for( int i=0; i<m_nprocs; i++ )
  {
    b[i] = m_wrapper->logical_and( b[i] );
    PostProcessingWriter::set_PostProcessingWindow( i, b[i] );
  }

  if ( m_rank == 0 )
  {
    cout << "\nProcessors that write Paraview particle files are : ";
    for( int i=0; i<m_nprocs; i++ )
      if ( b[i] > 0.5 ) cout << " " << i << endl;
  }
}




// ----------------------------------------------------------------------------
// Outputs timer summary */
void GrainsMPI::display_timer_summary()
{
  ostringstream oss;
  m_wrapper->timerSummary( oss );
  double cputime = CT_get_elapsed_time();
  oss << endl << "Full problem" << endl;
  write_elapsed_time_smhd( oss, cputime, "Computing time" );
  oss << "Mean number of particles on this sub-domain = " << 
	m_collision->getNbParticlesPerProcMean() << endl;     
  SCT_get_summary( oss, cputime );
  
  m_wrapper->writeStringPerProcess( cout, oss.str(), true, 
  	GrainsExec::m_shift0 ); 
}




// ----------------------------------------------------------------------------
// Returns a particle class among the classes of new particles to insert
Particle* GrainsMPI::getParticleClassForCreation( PullMode const& mode,
  	list< pair<Particle*,size_t> >& ParticleClassesForCreation,
	bool const& random_local )
{
  Particle* particleClass = NULL;
  list< pair<Particle*,size_t> >::iterator il;
  list<int> availableClasses;
  bool found = false;
  int i,i0;
  double v;
  
  switch ( mode ) 
  {
    case PM_ORDERED:
      for (il=ParticleClassesForCreation.begin(); 
      	il!=ParticleClassesForCreation.end() && !found;il++)
	if ( il->second != 0 )
	{
	  found = true;
	  particleClass = il->first;
	  --il->second;   
	}

      if ( !found )
      {
        cout << "No available particle in any class for creation " << endl;
	grainsAbort();
      }
      break;
      
    case PM_RANDOM:
      // We first find classes that still have particles available
      i=0;
      for (il=ParticleClassesForCreation.begin(); 
      	il!=ParticleClassesForCreation.end();il++,i++)
	if ( il->second != 0 ) availableClasses.push_back(i);
	
      // In local random mode, each process picks a random class
      if ( random_local )
      {
        v = double(random()) / double(INT_MAX);
        i0 = int(double(availableClasses.size()) * v);
      }      
      // In global random mode, the master process picks a random class and
      // sends it to the other processes
      else
      {
        if ( m_rank == 0 ) 
        {      
          v = double(random()) / double(INT_MAX);
          i0 = int(double(availableClasses.size()) * v);
        }
        i0 = m_wrapper->Broadcast_INT( i0 );
      } 

      il=ParticleClassesForCreation.begin();
      for (i=0;i<i0;i++) il++;
      particleClass = il->first;
      --il->second;          
      break;
      
    default:
      break;      
  }      

  return ( particleClass );
}




// ----------------------------------------------------------------------------
// Returns the number of insertion positions */
size_t GrainsMPI::getNbInsertionPositions() const
{
  size_t nbinsertionpos = m_wrapper->Broadcast_UNSIGNED_INT( 
	m_insertion_position->size() );
  return ( nbinsertionpos );
}




// ----------------------------------------------------------------------------
// Checks the clones when a simulation is reloaded */
void GrainsMPI::checkClonesReload()
{
  // If the linked cell size is smaller compared to that in the previous
  // simulation, we need to destroy the clones that are out of the 
  // linked cell   
  m_collision->DestroyOutOfDomainClones( m_time,
	m_allcomponents.getActiveParticles(),
	m_allcomponents.getCloneParticles(),
	m_wrapper ); 

  // If the linked cell size is larger compared to that in the previous
  // simulation or clones were not stored in the restart file, we need to 
  // create the missing clones
  m_wrapper->UpdateOrCreateClones_SendRecvLocal_GeoLoc( m_time,
	m_allcomponents.getActiveParticles(),
	m_allcomponents.getParticlesInBufferzone(),
	m_allcomponents.getCloneParticles(),
	m_allcomponents.getReferenceParticles(),
	m_collision, false, false );		
}
