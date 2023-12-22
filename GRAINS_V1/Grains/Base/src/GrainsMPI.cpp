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
  double vmax = 0., vmean = 0. ;

  // Timers
  SCT_insert_app( "ParticlesInsertion" );
  SCT_insert_app( "ComputeForces" );
  SCT_insert_app( "Move" );
  SCT_insert_app( "UpdateParticleActivity" );
  SCT_insert_app( "LinkUpdate" );
  SCT_insert_app( "OutputResults" );

  // Simulation: time marching algorithm
//  cout << "Time \t TO \tend \tParticles \tIn \tOut" << endl;
  while ( m_tend - m_time > 0.01 * m_dt )
  {
//     try
//     {
//       m_time += m_dt;
// 
// 
//       // Check whether data are output at this time
//       m_lastTime_save = false;
//       GrainsExec::m_output_data_at_this_time = false;
//       if ( m_timeSave != m_save.end() )
//         if ( *m_timeSave - m_time < 0.01 * m_dt )
// 	{
// 	  // Set the global data output boolean to true
// 	  GrainsExec::m_output_data_at_this_time = true;
// 
// 	  // Reset counter for force postprocessing
// 	  m_collision->resetPPForceIndex();
// 
// 	  // Next time of writing files
// 	  m_timeSave++;
// 
// 	  m_lastTime_save = true;
// 	}
// 
// 
//       // Insertion of particles
//       SCT_set_start( "ParticlesInsertion" );
//       if ( m_npwait_nm1 != m_allcomponents.getNumberInactiveParticles() )
//           cout << "\r                                              "
//                << "                 " << flush;
//       ostringstream oss;
//       oss.width(10);
//       oss << left << m_time;
//       cout << '\r' << oss.str() << "  \t" << m_tend << "\t\t\t"
//            << m_allcomponents.getNumberActiveParticlesOnProc() << '\t'
//            << m_allcomponents.getNumberInactiveParticles()    << flush;
//       m_npwait_nm1 = m_allcomponents.getNumberInactiveParticles();
//       if ( m_insertion_mode == IM_OVERTIME )
//         insertParticle( m_insertion_order );
//       SCT_get_elapsed_time( "ParticlesInsertion" );
// 
// 
//       if ( GrainsExec::m_TIScheme == "SecondOrderLeapFrog" )
//       {
//         // We implement the kick-drift-kick version of the LeapFrog scheme
// 
// 	// Move particles and obstacles
// 	// Update particle velocity over dt/2 and particle position over dt,
// 	// obstacle velocity and position over dt
// 	// v_i+1/2 = v_i + a_i * dt / 2
// 	// x_i+1 = x_i + v_i+1/2 * dt
//         moveParticlesAndObstacles( 0.5 * m_dt, m_dt, m_dt );
// 	
//         // Compute particle forces and acceleration
// 	// Compute f_i+1 and a_i+1 as a function of (x_i+1,v_i+1/2)
// 	SCT_set_start( "ComputeForces" );
//         computeParticlesForceAndAcceleration();
//         SCT_get_elapsed_time( "ComputeForces" );	
// 	
// 	// Update particle velocity over dt/2
// 	// v_i+1 = v_i+1/2 + a_i+1 * dt / 2 
// 	m_allcomponents.advanceParticlesVelocity( m_time, 0.5 * m_dt );
//       }
//       else
//       {
//         // Compute particle forces and acceleration
// 	// Compute f_i and a_i as a function of (x_i,v_i) 
// 	SCT_set_start( "ComputeForces" );
//         computeParticlesForceAndAcceleration();
//         SCT_get_elapsed_time( "ComputeForces" );
//       
//         // Move particles and obstacles
// 	// x_i+1 = x_i + g(v_i,v_i-1,a_i,dt)
// 	// v_i+1 = v_i + g(a_i,a_i-1,dt)
//         moveParticlesAndObstacles( m_dt, m_dt, m_dt );
//       }
// 
// 
//       // Write force & torque exerted on obstacles
//       m_allcomponents.outputObstaclesLoad( m_time, m_dt );
// 
// 
//       // Write postprocessing and reload files
//       if ( GrainsExec::m_output_data_at_this_time )
//       {
// 	SCT_set_start( "OutputResults" );
// 
// 	// Track component max and mean velocity
// 	m_allcomponents.ComputeMaxMeanVelocity( vmax, vmean );
//         cout << endl << "Component velocity : max = " << vmax
//                << " average = " << vmean << endl;
//         fVitMax << GrainsExec::doubleToString( ios::scientific, 6, m_time )
//                << "\t" << GrainsExec::doubleToString( ios::scientific, 6,
//                vmax ) << "\t" << GrainsExec::doubleToString( ios::scientific,
//                6, vmean ) << endl;
// 
// 	// Display memory used by Grains
// 	display_used_memory();
// 
// 	// Write reload files
// 	saveReload( m_time );
// 
// 	// Write postprocessing files
//         m_allcomponents.PostProcessing( m_time, m_dt, m_collision );
// 
// 	SCT_get_elapsed_time( "OutputResults" );
//       }
//     }
//     catch (ContactError &errContact)
//     {
//       // Max overlap exceeded
//       cout << endl;
//       m_allcomponents.PostProcessingErreurComponents( "ContactError",
//             errContact.getComponents() );
//       errContact.Message( cout );
//       m_error_occured = true;
//       break;
//     }
//     catch (DisplacementError &errDisplacement)
//     {
//       // Particle displacement over dt is too large
//       cout << endl;
//       m_allcomponents.PostProcessingErreurComponents( "DisplacementError",
//             errDisplacement.getComponent() );
//       errDisplacement.Message(cout);
//       m_error_occured = true;
//       break;
//     }
//     catch (SimulationError &errSimulation)
//     {
//       // Simulation error
//       cout << endl;
//       errSimulation.Message(cout);
//       m_error_occured = true;
//       break;
//     }
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
  m_rank = m_wrapper->get_rank_active();
  m_nprocs = m_wrapper->get_total_number_of_active_processes();
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
// Returns the full result file name
string GrainsMPI::fullResultFileName( string const& rootname ) const
{
  string fullname = rootname;
  ostringstream oss;
  oss << "_" << m_rank;
  fullname += oss.str()+".result";

  return ( fullname );
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
  for (int i=0;i<m_wrapper->get_total_number_of_active_processes();++i)
  {    
    if ( m_rank == i )
    {
      cout << endl;
      cout << "Processor " << m_rank << endl;
      m_wrapper->timerSummary();
      double cputime = CT_get_elapsed_time();
      cout << endl << "Full problem" << endl;
      write_elapsed_time_smhd( cout, cputime, "Computing time" );
      cout << "Mean number of particles on this sub-domain = " << 
	m_collision->getNbParticlesPerProcMean() << endl;     
      SCT_get_summary(cout,cputime);
    }
    m_wrapper->MPI_Barrier_ActivProc(); 
  }
}
