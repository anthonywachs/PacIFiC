/** 
# General DLMFD plugin 
*/

/** Definitions and global variables */
# include "DLMFD_Defs.h"

/** Fictitious domain implementation */
# include "DLMFD_Uzawa_velocity.h"

/** Paraview output functions */
# include "DLMFD_save_data_vtk.h"

/** Lambda criterion for visualizing wakes */
# include "lambda2.h"

/** Performance monitoring */
# include "DLMFD_Perf.h"



/** Generic granular solver init event: TO BE OVERLOADED */
//----------------------------------------------------------------------------
event GranularSolver_init (t < -1.)
//----------------------------------------------------------------------------
{} 




/** Generic granular solver predictor event: TO BE OVERLOADED */
//----------------------------------------------------------------------------
event GranularSolver_predictor (t < -1.)
//----------------------------------------------------------------------------
{} 




/** Generic granular solver velocity update event: TO BE OVERLOADED */
//----------------------------------------------------------------------------
event GranularSolver_updateVelocity (t < -1.)
//----------------------------------------------------------------------------
{}




/** Generic granular solver result saving event: TO BE OVERLOADED */
//----------------------------------------------------------------------------
event GranularSolver_saveResults (t < -1.)
//----------------------------------------------------------------------------
{}




/** Generic static embedded boundary cs and fs calculation event: TO BE 
OVERLOADED BY THE USER */
//----------------------------------------------------------------------------
event Compute_cs (t < -1.)
//----------------------------------------------------------------------------
{}




/** Overloading of the init event: initialize fluid and rigid bodies */
//----------------------------------------------------------------------------
event init (i = 0) 
//----------------------------------------------------------------------------
{
  if ( pid() == 0 )
  {
    printf( "==================================\n" );
    printf( "======   Basilisk + DLMFD   ======\n" );
    printf( "==================================\n" );        
  }   

  // Output basic fluid and geometric parameters
  if ( pid() == 0 )
  {
    printf( "Fluid density = %6.3e\n", FLUID_DENSITY );
    printf( "Fluid viscosity = %6.3e\n", FLUID_VISCOSITY );
    printf( "Space dimension = %d\n", dimension );    
    printf( "Domain size = %6.3e\n", L0 );
    printf( "\n" );            
  }  

# ifndef FLUID_DENSITY
#   define FLUID_DENSITY 1.
# endif
  
  // Initialize the density field
# if EMBED
    alpha = alphav;
    rho = rhov;
# else
    const scalar rhoc[] = FLUID_DENSITY;
    rho = rhoc;
# endif    

# ifndef FLUID_VISCOSITY
#   define FLUID_VISCOSITY 1.
# endif
  
  // Initialize the viscosity field
# if EMBED
    mu = muv;
# else
    const face vector muc[] = {FLUID_VISCOSITY, FLUID_VISCOSITY, 
    	FLUID_VISCOSITY};
    mu = muc;
# endif 
  
  // Initialize all DLMFD fields
  initialize_DLMFD_fields_to_zero();
  
  // Initiliaze the allDLMFDptscells stats to 0
  allDLMFDptscells.total_number_of_DLMFDpts = 0;  
  allDLMFDptscells.total_number_of_DLMFDcells = 0;

  // If new simulation: set fluid initial condition from user defined case file
  if ( ! restore ( FLUID_DUMP_FILENAME ) ) 
  {
    // Set the restarted simulation boolean to 0
    restarted_simu = 0;
  }
  else // Restart of a simulation 
  {
    // Set the restarted simulation boolean to 1
    restarted_simu = 1;

    // Read restart time
    if ( pid() == 0 ) printf( "Read t and dt from time restart file\n" );
    read_t_restart( DUMP_DIR, &trestart, &dtrestart, 
    	&imposed_periodicpressuredrop );  

    // Re-initialize the VTK writer
#   if PARAVIEW
      reinitialize_vtk_restart();	
#   endif	
  }

  
  // Initialize the granular solver
  if ( pid() == 0 ) 
  {  
    printf( "Granular solver initialization\n");
    printf( "------------------------------\n"); 
    printf( "Name : " );
  }    
  event( "GranularSolver_init" );
# if _MPI
    MPI_Barrier( MPI_COMM_WORLD );
# endif 


  // Creates boundary points and normal vectors of reference rigid bodies
  create_referencerigidbodies_boundary_geomfeatures( ReferenceRigidBodies, 
  	nbReferenceRigidBodies );
	
	
  // Print reference rigid bodies
  char outputshift[1]="";
  print_referencerigidbodies( ReferenceRigidBodies, nbReferenceRigidBodies, 
    	&outputshift[0] );	  

 
  // Check particle data files in case of restart
  if ( restarted_simu )
    check_particle_datafiles_restart( allRigidBodies, 
	nbRigidBodies, !RIGIDBODIES_AS_FIXED_OBSTACLES, trestart, dtrestart );


  // Initialize/open all DLMFD file pointers
  init_file_pointers( allRigidBodies, nbRigidBodies, pdata, 
  	!RIGIDBODIES_AS_FIXED_OBSTACLES, fdata, &converge, &cellvstime, 
	restarted_simu );      


  /* Construction of rigid bodies and their DLMFD features */
  if ( pid() == 0 ) printf( "Initial DLMFD RB construction\n" ); 
  DLMFD_construction();


  // Perform initial refinement around rigid bodies and write rigid body data
  int totalcell = 0;

  if ( restarted_simu ) // Restarted simulation
  {    
    totalcell = totalcells();
    if ( pid() == 0 ) printf( "Initial total cells = %d\n", totalcell );
    
    // If particles are lighter or neutrally buoyant and the flag 
    // B_SPLIT_EXPLICIT_ACCELERATION is on, we read the explication part of the 
    // particle split acceleration
#   if B_SPLIT_EXPLICIT_ACCELERATION
      read_explicit_splitAcceleration( DUMP_DIR, allRigidBodies, 
      	nbRigidBodies );
#   endif    
  }
  else // Simulation from t=0
  {
#   if ADAPTIVE
      if ( pid() == 0 )
      {
        printf( "\nInitial grid refinement in the rigid body volume\n" );
        printf( "------------------------------------------------\n" );      
      }
      
#     if INITIALGRIDADAPTIVE_NEWMETHOD == 0
            
        for (size_t k = 0; k < nbRigidBodies; k++) 
        {
          GeomParameter * gg;
          gg = &(allRigidBodies[k].g);    	  
    
          /* Perform initial refinement */
          totalcell = totalcells();
          if ( pid() == 0 )
	    printf( "# total cells before initial refinement of "
	  	"rigid body %d = %d\n", k, totalcell );

          coord position = gg->center;
          double radius = gg->radius;
          Point lpoint;
	
          /* Initial static refinement */
          lpoint = locate( position.x, position.y, position.z );
          int ilvl = 0;
          while( ilvl < MAXLEVEL ) 
          {
	    printf( "current level of refinement = %d, we want level %d\n",
	  	ilvl, MAXLEVEL );
	    printf ("refining ... \n");
	    if ( lpoint.level > -1 ) 
	    {
	      refine_cell( lpoint, all, 0, &tree->refined );
	      printf( "refined once on thread %d \n", pid() );
	    }
	    mpi_boundary_refine( all );
	    mpi_boundary_update( all );
      
	    lpoint = locate( position.x, position.y, position.z );
	    ilvl = lpoint.level;

#           if _MPI
	      mpi_all_reduce( ilvl, MPI_INT, MPI_MAX );
#           endif
          }
  
          radius = 1.2 * radius ;
          refine( sq( x - position.x ) + sq( y - position.y ) 
		+ sq( z - position.z ) - sq(radius) < 0. && level < MAXLEVEL );
	     
          // Compute the total number of cells after refinement 
          totalcell = totalcells();
          if ( pid() == 0 )
	    printf( "# total cells after initial refinement of rigid body %d = "
		"%d\n\n", k, totalcell );
        }
	      
#     else       
        // Create caches per rigid body boundary point
        Cache** stencil = (Cache **) calloc( nbRigidBodies, sizeof(Cache*) );
        for (size_t k = 0; k < nbRigidBodies; k++)
        { 
          stencil[k] = (Cache *) calloc( allRigidBodies[k].s.m, sizeof(Cache) );
	  for (int j = 0; j < allRigidBodies[k].s.m; j++)
	    initialize_and_allocate_Cache( &(stencil[k][j]) );
        }
      
        astats ss;
        int ic = 0, maxic = 100;
        double delta = L0 / (double)(1 << MAXLEVEL) ; 
        do 
        {
          ic++;
	
	  // Generate the 5^dim boundary point stencil
	  for (size_t k = 0; k < nbRigidBodies; k++)
            for (int j = 0; j < allRigidBodies[k].s.m; j++)
	    { 
              stencil[k][j].n = 0;
	      for (int ni=-2; ni<=2; ni++)
                for (int nj=-2; nj<=2; nj++) 
	        {
#                 if dimension < 3
                    Point point = locate( allRigidBodies[k].s.x[j] + ni * delta,
          		allRigidBodies[k].s.y[j] + nj * delta );
		    cache_append( &(stencil[k][j]), point, 0 );
#                 else
                    for (int nk=-2; nk<=2; nk++) 
		    {
                      Point point = locate( 
		      	allRigidBodies[k].s.bp[j].x + ni * delta,
          		allRigidBodies[k].s.bp[j].y + nj * delta,
			allRigidBodies[k].s.bp[j].z + nk * delta );
		      cache_append( &(stencil[k][j]), point, 0 );
		    }
#                 endif
                }
	    } 	
	
          // Assign the noisy distance function over cells that belong
	  // to each boundary point stencil
          foreach() DLM_Flag[] = 0; 	
	  for (size_t k = 0; k < nbRigidBodies; k++)
            for (int j = 0; j < allRigidBodies[k].s.m; j++) 	        
              foreach_cache(stencil[k][j]) 
                if ( point.level >= 0 ) 
	        {
                  coord dist;
                  dist.x = x - allRigidBodies[k].s.bp[j].x;
                  dist.y = y - allRigidBodies[k].s.bp[j].y;
#                 if dimension > 2
                    dist.z = z - allRigidBodies[k].s.bp[j].z;
#                 endif
#                 if dimension < 3
                    if ( fabs(dist.x) <= 2. * Delta 
		    		&& fabs(dist.y) <= 2. * Delta ) 
                      DLM_Flag[] = sq( dist.x + dist.y ) 
		      	/ sq( 2. * Delta ) * ( 2. + noise() );
#                 else
                    if ( fabs(dist.x) <= 2. * Delta 
		    		&& fabs(dist.y) <= 2. * Delta
          			&& fabs(dist.z) <= 2. * Delta )
                      DLM_Flag[] = sq( dist.x + dist.y + dist.z )
		      	/ cube( 2. * Delta ) * ( 2. + noise() );
#                 endif
                }

          // Run refinement using the noisy distance function
          ss = adapt_wavelet( {DLM_Flag}, (double[]) {1.e-16}, 
		maxlevel = MAXLEVEL, minlevel = LEVEL );

          totalcell = totalcells();
	  if ( pid() == 0 ) printf( "Refine initial mesh: step %d %d\n", ic,
		totalcell );
        } while ( ( ss.nf || ss.nc ) && ic < maxic );      
                        
        // Free all rigid body boundary point caches
        for (size_t k = 0; k < nbRigidBodies; k++)
        { 
          for (int j = 0; j < allRigidBodies[k].s.m; j++) 
	    free( stencil[k][j].p ); 
	  free( stencil[k] );
        }
        free(stencil);     
      
#     endif
#   endif
    
    // Save all rigid body data
    if ( !RIGIDBODIES_AS_FIXED_OBSTACLES )
      rigidbody_data( allRigidBodies, nbRigidBodies, t, i, pdata );  
  }

# if EMBED
    event( "Compute_cs" ); 
# endif 

  // Compute the AABB of the local domain
  compute_local_domain_AABB( &local_domain );
  
  // Simulation time interval
  maxtime = trestart + SIMUTIMEINTERVAL;
  if ( pid() == 0 ) 
    printf( "Simulation time interval: t_start = %f to t_end = %f\n\n", 
    	trestart, maxtime ); 

	
  // Initialize the field u_previoustime to compute x-velocity change
  foreach() u_previoustime[] = u.x[];
  
  
  // By default:
  // * do NOT DUMP the work fields used in the DLMFD sub-problem, the 
  // field u_previoustime and the pressure field
  // * DUMP the DLM_explicit field in the DLMFD sub-problem
  u_previoustime.nodump = true ;
  p.nodump = true ;
  DLM_Flag.nodump = true;
  DLM_FlagMesh.nodump = true ;
  foreach_dimension()
  {
    DLM_lambda.x.nodump = true ;
    DLM_Index.x.nodump = true ;
    DLM_PeriodicRefCenter.x.nodump = true ;
    DLM_r.x.nodump = true ;
    DLM_w.x.nodump = true ;
    DLM_v.x.nodump = true ;
    DLM_qu.x.nodump = true ;
    DLM_tu.x.nodump = true ;
#   if DLM_ALPHA_COUPLING
      DLM_explicit.x.nodump = false ;
#   endif 
  }
}



 
/** Logfile event to stop the simulation at maxtime */
//----------------------------------------------------------------------------
event logfile ( i=0; i++ ) 
//----------------------------------------------------------------------------
{
  // This is the condition that stops the simulation exactly at t=maxtime
  // We use -0.0001 * dt to avoid the problem of comparison of double that 
  // would exist if we would write if ( t > maxtime )
  if ( t - maxtime > - ROUNDDOUBLECOEF * dt ) 
  {
    // Close all DLMFD files
    close_file_pointers( nbRigidBodies, pdata, !RIGIDBODIES_AS_FIXED_OBSTACLES,
    	fdata, converge, cellvstime ); 

    // Write the dump time and time step for restart
    if ( pid() == 0 ) printf( "Write t and dt in time restart file\n" );
    save_t_dt_restart( DUMP_DIR, t, dt, imposed_periodicpressuredrop );  
 
    // Stop simulation
    return 1; 
  }
}




/** Output post-processing and restart data once */
//----------------------------------------------------------------------------
void do_output( char const* mess )
//----------------------------------------------------------------------------
{
  if ( pid() == 0 ) printf ( "Output data at %s %8.6f START\n", mess, t ); 

  scalar* dump_list = NULL;

# if LAMBDA2
    scalar l2[];
    lambda2( u, l2 );        
# endif

# if VORTICITY
#   if dimension == 2
      scalar omega[];
      vorticity( u, omega );
#   else
      vector omega[];      
      vorticity_3D( u, omega ); 
#   endif            
# endif

  // Paraview output for post-processing
# if PARAVIEW
    scalar* paraview_scalarlist = {PARAVIEW_SCALAR_LIST
#   if LAMBDA2
      , l2
#   endif
#   if VORTICITY
#     if dimension == 2
        , omega
#     endif            
#   endif         
    };
    vector* paraview_vectorlist = {PARAVIEW_VECTOR_LIST
#   if VORTICITY
      , omega
#   endif         
    };    
    save_data_vtk( paraview_scalarlist, paraview_vectorlist, allRigidBodies, 
    	nbRigidBodies, t );
# endif  
  
  // Bview output for post-processing   
# if BVIEW
    char dump_name[80] = "";
    char stime[80] = ""; 
    strcpy( dump_name, FLUID_DUMP_FILENAME );
    sprintf( stime, "%6f", t );
    strcat( dump_name, "_t" );
    strcat( dump_name, stime );     

    dump_list = (scalar *){BVIEW_LIST
#   if LAMBDA2
      , l2
#   endif     
#   if VORTICITY
      , omega
#   endif        
    };
    dump( dump_name, dump_list );  
# endif

  // Basilisk output for restart
# if DLM_ALPHA_COUPLING
    dump_list = (scalar *){u, g, DLM_explicit};
# else
    dump_list = (scalar *){u, g};    
# endif     
  dump( FLUID_DUMP_FILENAME, dump_list );
      
  // If particles are lighter or neutrally buoyant and the flag 
  // B_SPLIT_EXPLICIT_ACCELERATION is on, we save the explication part of the 
  // particle split acceleration for restart purposes
# if B_SPLIT_EXPLICIT_ACCELERATION
    if ( pid() == 0 ) 
      save_explicit_splitAcceleration( DUMP_DIR, allRigidBodies, 
      	nbRigidBodies );
# endif        

  // Granular solver output for both post-processing & restart
  event( "GranularSolver_saveResults" );
     
  if ( pid() == 0 ) printf ( "Output data at %s %8.6f END\n", mess, t );
}




/** Output post-processing and restart data at regular time intervals */
//----------------------------------------------------------------------------
event output_data (t += TINTERVALOUTPUT; 
 	t - ( maxtime < ROUNDDOUBLECOEF * TINTERVALOUTPUT ? 1.e20 : maxtime ) 
	< ROUNDDOUBLECOEF * dt )
//----------------------------------------------------------------------------	
{
  char mess[5] = "TIME";
  save_data_restart = true;

  // Freeing and reconstructing RB is required here as the remeshing step
  // might have induced some interpolation errors on the integer fields
  // DLM Flag, DLM FlagMesh and DLM_Index
# if DLMFD_INTERIORPOINTS || DLMFD_BOUNDARYPOINTS
    free_rigidbodies( allRigidBodies, nbRigidBodies, false );
    DLMFD_construction();  
# endif  

  do_output( mess );  
}




/** Output post-processing and restart data at the last time */
//----------------------------------------------------------------------------
event last_output_data (t = end) 
//----------------------------------------------------------------------------
{ 
  DLMFD_construction();

  char mess[10] = "LAST TIME";

  // Freeing and reconstructing RB is required here as the remeshing step
  // might have induced some interpolation errors on the integer fields
  // DLM Flag, DLM FlagMesh and DLM_Index
# if DLMFD_INTERIORPOINTS || DLMFD_BOUNDARYPOINTS
    free_rigidbodies( allRigidBodies, nbRigidBodies, false );
    DLMFD_construction();  
# endif 

  do_output( mess );    
  output_dlmfd_perf ( &DLMFD_UzawaTiming, &DLMFD_ConstructionTiming, i, 
  	&allDLMFDptscells );	       
}	




/** Writes the dump time and time step for restart, frees dynamic features 
of rigid bodies and predict the position & velocity of rigid bodies */
//----------------------------------------------------------------------------
event once_timestep_is_determined (i++)
//----------------------------------------------------------------------------
{
  if ( save_data_restart && i )
  {
    if ( pid() == 0 ) printf( "Write t and dt in time restart file\n" );
    save_t_dt_restart( DUMP_DIR, t, dt, imposed_periodicpressuredrop );        
    save_data_restart = false;   
  } 
  if ( i == 0 ) save_data_restart = false;

  if ( pid() == 0 )
  {
    printf( "\n****** ITER %d - TIME t=%8.5e to t+dt=%8.5e ******\n", 
    	i, t, t+dt );
    printf( "   Time step = %8.5e\n", dt );
  }
  
  // We free dynamic features of rigid bodies from the previous time step
  // or from initialization at the start of the current time step 
  // (run with -events to understand)
  free_rigidbodies( allRigidBodies, nbRigidBodies, true );  
  
  // In case of a periodic flow, we add the imposed pressure drop
# if IMPOSED_PERIODICFLOW
#   if IMPOSED_PERIODICFLOW_DIRECTION == 0 
      const face vector dp[] = { 
	- imposed_periodicpressuredrop / ( L0 * FLUID_DENSITY ), 0., 0. };
#   elif IMPOSED_PERIODICFLOW_DIRECTION == 1
      const face vector dp[] = { 0.,
	- imposed_periodicpressuredrop / ( L0 * FLUID_DENSITY ), 0. };  	
#   else 
      const face vector dp[] = { 0., 0., 
		- imposed_periodicpressuredrop / ( L0 * FLUID_DENSITY ) };  
#   endif
    a = dp;
# endif  
  
  /* Granular solver predictor */
  if ( pid() == 0 ) printf( "   GS predictor step: " ); 
  event( "GranularSolver_predictor" );  

  /* Construction of rigid bodies and their DLMFD features */
  if ( pid() == 0 ) printf( "   DLMFD RB construction\n" ); 
  DLMFD_construction();  
}




/** Frees dynamic features of rigid bodies at the very end */
//----------------------------------------------------------------------------
event cleanup (t = end)
//----------------------------------------------------------------------------
{
  // This is because we free dynamic features of rigid bodies from the previous 
  // time step at the start of the current time step (see above the event 
  // once_timestep_is_determined)
  // Since at the very end, there is no next time step, dynamic features 
  // of rigid bodies are not freed by start_timestep and hence are freed here
  free_rigidbodies( allRigidBodies, nbRigidBodies, true );
  free_rigidbodies( ReferenceRigidBodies, nbReferenceRigidBodies, true );   
  free_np_dep_arrays( nbParticles, allRigidBodies, RBnumToIndex, DLMFDtoGS_vel,
  	vpartbuf, pdata, !RIGIDBODIES_AS_FIXED_OBSTACLES, fdata );
}




/** Overloading of the viscous_term event: we add an explicit coupling term
that improves the coupling between the fluid problem and the DLMFD problem */
//----------------------------------------------------------------------------
event viscous_term (i++) 
//----------------------------------------------------------------------------
{
# if DLM_ALPHA_COUPLING
    foreach()
      foreach_dimension()
        u.x[] += - dt * DLM_explicit.x[] / ( FLUID_DENSITY * dlmfd_dv() );
# endif
}




/** Solves DLMFD velocity problem */
//----------------------------------------------------------------------------
void do_DLMFD( const int i )
//----------------------------------------------------------------------------
{
  /* Solve the DLMFD velocity problem by a Uzawa algorithm */ 
  DLMFD_Uzawa_velocity( i );

  /* Update velocity in the granular solver */
  if ( pid() == 0 && nbParticles )
  { 
    printf( "   GS velocity update in " );
    event( "GranularSolver_updateVelocity" );
  }    

  /* Save the forces acting on rigid bodies before adapting the mesh  */
  computeHydroForceTorque( allRigidBodies, nbRigidBodies, fdata, t + dt, dt, 
  	DLM_Flag, DLM_lambda, DLM_Index, FLUID_DENSITY, DLM_PeriodicRefCenter );
  
  /* Save all rigid body data */
  if ( !RIGIDBODIES_AS_FIXED_OBSTACLES )
    rigidbody_data( allRigidBodies, nbRigidBodies, t + dt, i, pdata );     
}




//----------------------------------------------------------------------------
event after_viscous_term (i++) 
//----------------------------------------------------------------------------
{ 
# if !DLMFD_PROB_AFTER_NAVIERSTOKES 
    do_DLMFD( i );
# endif 
}




/** Overloading of the end_timestep event: we plug the solution to the 
granular problem followed by the solution to the DLMFD problem */
//----------------------------------------------------------------------------
event end_timestep (i++) 
//----------------------------------------------------------------------------
{
  if ( pid() == 0 )
    printf( "   NS solver: MGu_niter = %d, MGp_niter = %d\n", mgu.i, mgp.i );

# if DLMFD_PROB_AFTER_NAVIERSTOKES 
    do_DLMFD( i );
# endif

  /* Compute periodic flow rate or adjust periodic pressure drop if 
  imposed periodic flow rate */
# if IMPOSED_PERIODICFLOW
    double flowrate = 0.;
#   if IMPOSED_PERIODICFLOW_TYPE == 0
#     if IMPOSED_PERIODICFLOW_DIRECTION == 0 
        flowrate = compute_flowrate_xperiodic( u );
#     elif IMPOSED_PERIODICFLOW_DIRECTION == 1
        flowrate = compute_flowrate_yperiodic( u );
#     else 
        flowrate = compute_flowrate_zperiodic( u );
#     endif 
      if ( pid() == 0 )
        printf( "   Periodic flow rate = %8.5e\n", flowrate );           
#   else
      double Q1 = - dt / ( L0 * FLUID_DENSITY );
      double deltaflowrate = 0.;                    
#     if IMPOSED_PERIODICFLOW_DIRECTION == 0 
        flowrate = compute_flowrate_xperiodic( u );
	deltaflowrate = imposed_periodicflowrate - flowrate;
	imposed_periodicpressuredrop += deltaflowrate / ( Q1 * L0 * L0 );
	foreach()
	  u.x[] += deltaflowrate / ( L0 * L0 );  
#     elif IMPOSED_PERIODICFLOW_DIRECTION == 1
        flowrate = compute_flowrate_yperiodic( u );
	deltaflowrate = imposed_periodicflowrate - flowrate;
	imposed_periodicpressuredrop += deltaflowrate / ( Q1 * L0 * L0 );
	foreach()
	  u.y[] += deltaflowrate / ( L0 * L0 ); 
#     else 
        flowrate = compute_flowrate_zperiodic( u );
	deltaflowrate = imposed_periodicflowrate - flowrate;		
	imposed_periodicpressuredrop += deltaflowrate / ( Q1 * L0 * L0 );
	foreach()
	  u.z[] += deltaflowrate / ( L0 * L0 );    
#     endif
      synchronize((scalar *){u});
      if ( pid() == 0 )
        printf( "   Periodic pressure drop = %8.5e %8.5e\n", 
		imposed_periodicpressuredrop, flowrate );
#   endif  
# endif
   
  
  /* Fluid velocity change over the time step */
  deltau = change( u.x, u_previoustime );
  if ( pid() == 0 )
    printf( "   Velocity change = %8.5e\n", deltau );       
}




/** Overloading of the adapt event: we refine the mesh with both the velocity
field and a phase indicator */
//----------------------------------------------------------------------------
event adapt (i++) 
//----------------------------------------------------------------------------
{
# if ADAPTIVE
    int totalcell = totalcells();
    if ( pid() == 0 )
    {
#     if dimension == 2  
        printf( "   Quadtree grid cells: " );
#     else
        printf( "   Octree grid cells: " );
#     endif
      printf( "total = %d, ", totalcell );
    }

# if LEVELDIFF_FLAG_U == 0
    astats s = adapt_wavelet( (scalar *){DLM_FlagMesh, u}, 
	(double[]){FLAG_ADAPT_CRIT, UX_ADAPT_CRIT, UY_ADAPT_CRIT, 
	UZ_ADAPT_CRIT}, maxlevel = MAXLEVEL, minlevel = LEVEL );	
# else	
    // Flag regions around rigid bodies
    flag_rigidbodies_with_boundarylayers( allRigidBodies, nbRigidBodies,
    	DLM_FlagMaxLev, BOUNDARY_LAYER_THICKNESS_COEF, &local_domain );       

    astats s = adapt_wavelet_spatial ( DLM_FlagMaxLev, MAXLEVEL,
	(scalar *){DLM_FlagMesh, u}, (double[]){FLAG_ADAPT_CRIT, UX_ADAPT_CRIT,
	UY_ADAPT_CRIT, UZ_ADAPT_CRIT}, MAXLEVEL-LEVELDIFF_FLAG_U, 
	minlevel = LEVEL );	
# endif	
	
# if EMBED
    event( "Compute_cs" ); 
# endif	

    if ( pid() == 0 ) 
      printf( "refined = %d, coarsened = %d\n", s.nf, s.nc );
# endif

  // Compute the AABB of the local domain
  compute_local_domain_AABB( &local_domain );	
}




/* Update fluid viscosity and density in case of embedded boundaries */
//----------------------------------------------------------------------------
event properties (i++) 
//----------------------------------------------------------------------------
{
# if EMBED
    foreach_face() 
    {
      muv.x[] = FLUID_VISCOSITY * fm.x[];
      alphav.x[] = fm.x[] / FLUID_DENSITY ;
    }
    foreach() rhov[] = FLUID_DENSITY * cm[];
# endif
}
