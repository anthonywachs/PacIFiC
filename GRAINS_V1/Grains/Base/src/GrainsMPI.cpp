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
  while ( m_tend - m_time > 0.01 * m_dt )
  {
    try
    {
      m_time += m_dt;


      // Check whether data are output at this time
      m_lastTime_save = false;
      GrainsExec::m_output_data_at_this_time = false;
      if ( m_timeSave != m_save.end() )
        if ( *m_timeSave - m_time < 0.01 * m_dt )
	{
	  // Set the global data output boolean to true
	  GrainsExec::m_output_data_at_this_time = true;
	  if ( m_rank == 0 ) cout << endl << "Time = " << m_time << endl;

	  // Reset counter for force postprocessing
	  m_collision->resetPPForceIndex();

	  // Next time of writing files
	  m_timeSave++;

	  m_lastTime_save = true;
	  
	}


      // Insertion of particles
      SCT_set_start( "ParticlesInsertion" );
      m_allcomponents.computeNumberParticles( m_wrapper );
      if ( GrainsExec::m_output_data_at_this_time )
        if ( m_rank == 0 )
	{
	  if ( m_allcomponents.getNumberActiveParticlesOnAllProc() !=
	  	m_allcomponents.getTotalNumberPhysicalParticles() ) 
	  {
	    cout << "Number of active particles on all proc = " <<
	  	m_allcomponents.getNumberActiveParticlesOnAllProc() << endl;
	    cout << "Number of inactive particles = " <<
	  	m_allcomponents.getNumberInactiveParticles() << endl;
	  }		
	}
      if ( m_insertion_mode == IM_OVERTIME )
        insertParticle( m_insertion_order );
      SCT_get_elapsed_time( "ParticlesInsertion" );

      
      // Move particles and obstacles
      // Update particle velocity over dt/2 and particle position over dt,
      // obstacle velocity and position over dt
      // v_i+1/2 = v_i + a_i * dt / 2
      // x_i+1 = x_i + v_i+1/2 * dt
      // Solve Newton's law and move particles
      SCT_set_start( "Move" );
      m_allcomponents.Move( m_time, 0.5 * m_dt, m_dt, m_dt );

   
      // Update particle position and velocity
      
      
      SCT_get_elapsed_time( "Move" );

      // Write force & torque exerted on obstacles
      m_allcomponents.outputObstaclesLoad( m_time, m_dt );


      // Write postprocessing and reload files
      if ( GrainsExec::m_output_data_at_this_time )
      {
	SCT_set_start( "OutputResults" );
		
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

	// Write reload files
	saveReload( m_time );

	// Write postprocessing files
        m_allcomponents.PostProcessing( m_time, m_dt, m_collision,
		m_rank, m_nprocs, m_wrapper );

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
    catch (DisplacementError &errDisplacement)
    {
      // Particle displacement over dt is too large
      cout << endl;
      m_allcomponents.PostProcessingErreurComponents( "DisplacementError",
            errDisplacement.getComponent() );
      errDisplacement.Message(cout);
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
// Attempts to insert a particle in the simulation
bool GrainsMPI::insertParticle( PullMode const& mode )
{
  static size_t insert_counter = 0;
  pair<bool,bool> insert(false,false);
  Vector3 vtrans, vrot ;
  Point3 position;
  Matrix mrot;
  Transform trot;    

  if ( insert_counter == 0 )
  {
    // In parallel, in case of an insertion via a position file or a structured
    // array, particles are created in this process if they are in the local 
    // domain. Therefore, each process has a different list of particles to 
    // insert. The random selection of particles is performed
    // in setPositionParticlesFromFile or setPositionParticlesArray, therefore
    // the list of particles to insert is scanned in an ordered way
    // In the case of insertion via a window, the procedure is similar to the 
    // serial insertion, all processes have the same list of particles to insert
    // and the random selection of particles is performed is performed in 
    //  m_allcomponents.getParticle( mode ) when mode == PM_RANDOM
    Particle *particle = NULL;
    if ( m_position != "" ) 
      particle = m_allcomponents.getParticle( PM_ORDERED );
    else particle = m_allcomponents.getParticle( mode );
    
    if ( particle )
    {
      if ( m_position == "" ) 
      {      
        // Initialisation of the particle velocity
        if ( m_rank == 0 ) computeInitVit( vtrans, vrot );
        vtrans = m_wrapper->Broadcast_Vector3( vtrans ); 
        vrot = m_wrapper->Broadcast_Vector3( vrot ); 
        particle->setTranslationalVelocity( vtrans );
        particle->setAngularVelocity( vrot );

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
      }
      else
      {
        // Initialisation of the particle velocity
        computeInitVit( vtrans, vrot );
        particle->setTranslationalVelocity( vtrans );
        particle->setAngularVelocity( vrot ); 
	
        // Initialisation of the angular position of the particle
        // Rem: we compose to the right by a pure rotation as the particle
        // already has a non-zero position that we do not want to change (and
        // that we would change if we would compose to the left)
        if ( m_init_angpos == IAP_RANDOM )
        {
          trot.setBasis( GrainsExec::RandomRotationMatrix( m_dimension ) );
          particle->composePositionRightByTransform( trot );
        }	     
      }

      particle->initialize_transformWithCrust_to_notComputed();

      // If insertion if successful, shift particle from wait to inserted
      // and initialize particle rotation quaternion from rotation matrix
      insert = m_collision->insertParticleParallel( particle,
      	m_allcomponents.getActiveParticles(), m_force_insertion,
	m_wrapper );
	
      // TO DO: periodic insertion !!
	
      // If no contact
      if ( !insert.second )
      {
        // If particle is in LinkedCell
	if ( insert.first )
	{
	  m_allcomponents.ShiftParticleOutIn();
	  Quaternion qrot;
	  qrot.setQuaternion( particle->getRigidBody()->getTransform()
		->getBasis() );
	  particle->setQuaternionRotation( qrot );
        }
	// If not in LinkedCell
        else
          m_allcomponents.DeleteAndDestroyWait();
      }
    }
  }

  ++insert_counter;
  if ( insert_counter == m_insertion_frequency ) insert_counter = 0;

  return ( !insert.second && insert.first );
}




// ----------------------------------------------------------------------------
// Creates, inserts and links new particles in the simulation
void GrainsMPI::InsertCreateNewParticles()
{
  // IMPORTANT: for any system with N particles and M obstacles, particles are
  // numbered 0 to N-1 and obstacles are numbered N to N+M-1
  // In the case of reload with additional insertion, it means that obstacles
  // are re-numbered

  int numPartMax = getMaxParticleIDnumber();
  list< pair<Particle*,int> >::iterator ipart;

  // New particles construction
  Component::setMaxIDnumber( numPartMax );
  if ( m_position != "" )
  {
    // From a structured array
    if ( m_position == "STRUCTURED" )
      setPositionParticlesArray( m_insertion_order );
    // From a file
    else setPositionParticlesFromFile( m_insertion_order );
  }
  else
  {
    for (ipart=m_newParticles.begin();ipart!=m_newParticles.end();ipart++)
    {
      int nbre = ipart->second;
      for (int ii=0; ii<nbre; ii++)
      {
        Particle* particle = ipart->first->createCloneCopy( true );
        m_allcomponents.AddParticle( particle );
      }
    }
  }

  // Obstacle renumbering
  int numInitObstacles = numPartMax;
  for (ipart=m_newParticles.begin();ipart!=m_newParticles.end();ipart++)
    numInitObstacles += ipart->second;
  list<SimpleObstacle*> lisObstaclesPrimaires =
    	m_allcomponents.getObstacles()->getObstacles();
  list<SimpleObstacle*>::iterator iobs;
  int ObstacleID = numInitObstacles;
  for (iobs=lisObstaclesPrimaires.begin();iobs!=lisObstaclesPrimaires.end();
    	iobs++)
  {
    (*iobs)->setID( ObstacleID );
    ++ObstacleID;
  }

  // Link all components with the grid
  m_allcomponents.Link( *m_collision );

  // Set particle positions from file or from a structured array
  if ( m_position != "" )
  {
    // From a structured array
    if ( m_position == "STRUCTURED" )
      setPositionParticlesArray( m_insertion_order );
    // From a file
    else setPositionParticlesFromFile( m_insertion_order );
  }

  // Insertion at initial time
  size_t nbPW = 0 ;
  bool b_insertion_BEFORE = true;  
  if ( m_insertion_mode == IM_INITIALTIME )
  {
    nbPW = m_allcomponents.getNumberInactiveParticles() ;
    for (size_t i=0;i<nbPW && b_insertion_BEFORE;++i)
      b_insertion_BEFORE = insertParticle( m_insertion_order );
    if ( !b_insertion_BEFORE )
    {
      cout << "Process " << m_rank << endl;
      cout << "Insertion issue with defined position method" << endl;
      cout << "Remaining number of particles to insert = " << 
      	m_allcomponents.getNumberInactiveParticles() << endl;
    }     

    // TO DO: create the MPI clones
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
void GrainsMPI::setPositionParticlesFromFile( const PullMode& mode )
{
  Particle* newPart = NULL;
  Particle* particleClass = NULL;
  int id;
  Point3 position;
  list< pair<Particle*,int> > newParticles_ = m_newParticles;
  list< pair<Particle*,int> >::iterator ipart;

  // Check that position file exists 
  ifstream filePos;
  size_t notok = 0;
  if ( m_rank == 0 )
  {
    filePos.open( m_position.c_str() );
    if ( !filePos.is_open() ) notok = 1;
  }
  notok = m_wrapper->sum_UNSIGNED_INT( notok );
  if ( notok )
  {
    if ( m_rank == 0 ) cout << "ERR : Position file does not exist" << 
    	m_position << endl << endl;  
    m_wrapper->MPI_Barrier_ActivProc();	
    grainsAbort();
  }

  // Check that the number of positions matches the number of particles
  // to insert
  size_t nwait = 0, npos = 0; 
  if ( m_rank == 0 )
  {
    for (ipart=m_newParticles.begin();ipart!=m_newParticles.end();ipart++)
      nwait += ipart->second;
    while ( filePos >> position ) ++npos;
    if ( nwait != npos ) notok = 1;
    filePos.close();
  }
  notok = m_wrapper->sum_UNSIGNED_INT( notok );
  if ( notok )
  {
    if ( m_rank == 0 )
      cout << "ERR: number of particles to insert is "
           << "different from the number of positions in the file " 
	   << m_position << " : " << nwait << " != " << npos << endl << endl;
    m_wrapper->MPI_Barrier_ActivProc();	
    grainsAbort();
  }
  
  // Initialisation of particle numbering
  id = getMaxParticleIDnumber();

  // All processes read the entire position file, but a particle is created
  // on this process if the position is inside the local domain of the process
  filePos.open( m_position.c_str() );
  while ( filePos >> position ) 
  {
    // Select a particle class
    particleClass = getParticleClassForCreation( mode, newParticles_, 
    	false );
       
    // Create the particle only if in the local domain
    if ( App::isInLocalDomain( &position ) )
    {
      newPart = new Particle( *particleClass );
      newPart->setPosition( position );
      newPart->setID( id );
      m_allcomponents.AddParticle( newPart ); 
    }
    id++;      	
  }
  filePos.close();
}




// ----------------------------------------------------------------------------
// Sets particle initial position with a structured array
void GrainsMPI::setPositionParticlesArray( const PullMode& mode )
{
  Particle* newPart = NULL; 
  Particle* particleClass = NULL;
  int id, k, l, m, kmin = 0, kmax = int(m_InsertionArray->NX) - 1, 
  	lmin = 0, lmax = int(m_InsertionArray->NY) - 1, 
	mmin = 0, mmax = int(m_InsertionArray->NZ) - 1,
	npartproc;
  Point3 position;
  list< pair<Particle*,int> > newParticles_ = m_newParticles;
  bool found = false, no_overlap = false ;

  // Test the number of particles to insert
  size_t ntotalinsert = 0;
  list< pair<Particle*,int> >::iterator ipart,ipart2;
  for (ipart=m_newParticles.begin();ipart!=m_newParticles.end();ipart++)
    ntotalinsert += ipart->second; 

  if ( ntotalinsert != m_InsertionArray->NX * m_InsertionArray->NY 
  	* m_InsertionArray->NZ )
  {
    if ( m_rank == 0 )
      cout << "ERR: number of particles to insert is different from"
    	<< " the number of positions in the structured array : " << 
	ntotalinsert << " != " <<
	m_InsertionArray->NX * m_InsertionArray->NY * m_InsertionArray->NZ 
	<< endl;
    m_wrapper->MPI_Barrier_ActivProc();	
    grainsAbort();	
  }


  // Uniform cell size of the structured array
  double deltax = 
  	( m_InsertionArray->box.ptB[X] - m_InsertionArray->box.ptA[X] ) 
  	/ double(m_InsertionArray->NX) ;
  double deltay = 
  	( m_InsertionArray->box.ptB[Y] - m_InsertionArray->box.ptA[Y] ) 
  	/ double(m_InsertionArray->NY) ;  
  double deltaz = 
  	( m_InsertionArray->box.ptB[Z] - m_InsertionArray->box.ptA[Z] ) 
  	/ double(m_InsertionArray->NZ) ;  


  // Coordinates of local domain
  Point3 localOrigin;
  App::get_local_domain_origin( localOrigin[X], localOrigin[Y], 
  	localOrigin[Z] );
  Vector3 localSize;
  App::get_local_domain_size( localSize[X], localSize[Y], localSize[Z] );
  Point3 MaxLocal = localOrigin + localSize;	
    

  // Test that no coordinate matches exactly the local domain limits
  // If any does, translate by 1e-12
  double geoshift = 1.e-12;
  vector<size_t> coorMatchLocLim( 3, 0 );
  for (k=0;k<int(m_InsertionArray->NX) && !coorMatchLocLim[0];++k)
  {
    position[X] = m_InsertionArray->box.ptA[X] + ( double(k) + 0.5 ) * deltax;
    if ( fabs( position[X] - localOrigin[X] ) < geoshift
    	|| fabs( position[X] - MaxLocal[X] ) < geoshift ) 
      coorMatchLocLim[0] = 1;
  }
  
  for (l=0;l<int(m_InsertionArray->NY) && !coorMatchLocLim[1];++l)
  {
    position[Y] = m_InsertionArray->box.ptA[Y] + ( double(l) + 0.5 ) * deltay;
    if ( fabs( position[Y] - localOrigin[Y] ) < geoshift
    	|| fabs( position[Y] - MaxLocal[Y] ) < geoshift ) 
      coorMatchLocLim[1] = 1;
  }  

  for (m=0;m<int(m_InsertionArray->NZ) && !coorMatchLocLim[2];++m)
  {
    position[Z] = m_InsertionArray->box.ptA[Z] + ( double(m) + 0.5 ) * deltaz;
    if ( fabs( position[Z] - localOrigin[Z] ) < geoshift
    	|| fabs( position[Z] - MaxLocal[Z] ) < geoshift ) 
      coorMatchLocLim[2] = 1;
  }
  
  coorMatchLocLim[0] = m_wrapper->max_UNSIGNED_INT( coorMatchLocLim[0] );
  if ( coorMatchLocLim[0] ) m_InsertionArray->box.ptA[X] += geoshift;
  coorMatchLocLim[1] = m_wrapper->max_UNSIGNED_INT( coorMatchLocLim[1] );
  if ( coorMatchLocLim[1] ) m_InsertionArray->box.ptA[Y] += geoshift;  
  coorMatchLocLim[2] = m_wrapper->max_UNSIGNED_INT( coorMatchLocLim[2] );
  if ( coorMatchLocLim[2] ) m_InsertionArray->box.ptA[Z] += geoshift;
  
  if ( ( coorMatchLocLim[0] || coorMatchLocLim[1] || coorMatchLocLim[2] )
  	&& m_rank == 0 )
  {
    cout << endl << "Warning: Structured array positions: some coordinates"
    	<< " exactly match local domain limits in the following directions:"
	<< endl;
    for (size_t i=0;i<3;++i)
      if ( coorMatchLocLim[i] )
        cout << "   * " << ( i == 0 ? "X" : i == 1 ? "Y" : "Z" ) << 
      	" automatix translation of " << geoshift << endl;	
  }	      
  m_wrapper->MPI_Barrier_ActivProc();	

  // Number of particles to insert on this process
  // We determine the min and max indices of particles located in this local
  // domain
  k = 0;
  found = false;
  while( !found && k < int(m_InsertionArray->NX) )
  {
    position[X] = m_InsertionArray->box.ptA[X] + ( double(k) + 0.5 ) * deltax;
    if ( position[X] > localOrigin[X] ) {kmin = k; found = true;}
    ++k;
  }
  if ( found )
  {
    k = kmax;
    found = false;
    while( !found && k >= 0 )
    {
      position[X] = m_InsertionArray->box.ptA[X] + ( double(k) + 0.5 ) * deltax;
      if ( position[X] < MaxLocal[X] ) {kmax = k; found = true;}
      --k;
    }
    if ( !found ) no_overlap = true ;
  }
  else no_overlap = true ;
  
  if ( !no_overlap )
  {    
    l = 0;
    found = false;  
    while( !found && l < int(m_InsertionArray->NY) )
    {
      position[Y] = m_InsertionArray->box.ptA[Y] + ( double(l) + 0.5 ) * deltay;
      if ( position[Y] > localOrigin[Y] ) {lmin = l; found = true;}
      ++l;
    }
    if ( found )
    {
      l = lmax; 
      found = false;
      while( !found && l >= 0 )
      {
        position[Y] = m_InsertionArray->box.ptA[Y] 
		+ ( double(l) + 0.5 ) * deltay;
        if ( position[Y] < MaxLocal[Y] ) {lmax = l; found = true;}
        --l;
      }
      if ( !found ) no_overlap = true ;
    }
    else no_overlap = true ;
  }
  
  if ( !no_overlap )
  {    
    m = 0;
    found = false; 
    while( !found && m < int(m_InsertionArray->NZ) )
    { 
      position[Z] = m_InsertionArray->box.ptA[Z] 
      	+ ( double(m) + 0.5 ) * deltaz;
      if ( position[Z] > localOrigin[Z] ) {mmin = m; found = true;}
      ++m;
    }
    if ( found )
    {
      m = mmax;
      found = false;
      while( !found && m >= 0 )
      {
        position[Z] = m_InsertionArray->box.ptA[Z] 
		+ ( double(m) + 0.5 ) * deltaz;
        if ( position[Z] < MaxLocal[Z] ) {mmax = m; found = true;}
        --m;
      }     
      if ( !found ) no_overlap = true ;
    }
    else no_overlap = true ;
  }  
  
  if ( no_overlap )
  {
    npartproc = kmin = lmin = mmin = 0;
    kmax = lmax = mmax = -1;
  }
  else npartproc = ( kmax - kmin + 1 ) * ( lmax - lmin + 1 )
  	* ( mmax - mmin + 1 );

        	
  // Number of particles of each class to insert in this local domain
  // Warning: works only if the number of particles per class is larger
  // than the number of processes
  if ( m_newParticles.size() == 1 )
    newParticles_.front().second = int(npartproc);
  else
    m_wrapper->distributeParticlesClassProc( m_newParticles,
    	newParticles_, npartproc, ntotalinsert );	
 
   
  // Initialisation of particle numbering
  size_t* vnbpartproc = m_wrapper->AllGather_UNSIGNED_INT( npartproc );  
  id = getMaxParticleIDnumber();
  for (k=0;k<m_rank;++k) id += int(vnbpartproc[k]);
  delete [] vnbpartproc;

  
  // Assign positions
  for (k=kmin;k<=kmax;++k)
    for (l=lmin;l<=lmax;++l)
      for (m=mmin;m<=mmax;++m)    
      {
	// Position
        position[X] = m_InsertionArray->box.ptA[X] 
		+ ( double(k) + 0.5 ) * deltax;
        position[Y] = m_InsertionArray->box.ptA[Y] 
		+ ( double(l) + 0.5 ) * deltay;	
        position[Z] = m_InsertionArray->box.ptA[Z] 
		+ ( double(m) + 0.5 ) * deltaz;

        // Select a particle class
        particleClass = getParticleClassForCreation( mode, newParticles_,
		true );

        // Create the particle
	newPart = new Particle( *particleClass );
        newPart->setPosition( position );
        newPart->setID( id );
        m_allcomponents.AddParticle( newPart );
	id++;   
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




// ----------------------------------------------------------------------------
// Returns a particle class among the classes of new particles to insert
Particle* GrainsMPI::getParticleClassForCreation( PullMode const& mode,
  	list< pair<Particle*,int> >& ParticleClassesForCreation,
	bool const& random_local )
{
  Particle* particleClass = NULL;
  list< pair<Particle*,int> >::iterator il;
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
