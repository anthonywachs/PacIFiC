/** 
# The general DLMFD plugin 
*/

/** File names definition and global variables */  
# ifndef FLUID_DUMP_FILENAME
#   define FLUID_DUMP_FILENAME "Savings/dump"
# endif
# ifndef DUMP_DIR
#   define DUMP_DIR "Savings"
# endif

# ifndef RESULT_DIR
#   define RESULT_DIR "Res"
# endif
# ifndef RESULT_RIGIDBODY_VP_ROOTFILENAME
#   define RESULT_RIGIDBODY_VP_ROOTFILENAME "rigidbody_data"
# endif
# ifndef RESULT_RIGIDBODY_HYDROFAT_ROOTFILENAME
#   define RESULT_RIGIDBODY_HYDROFAT_ROOTFILENAME "rigidbody_hydroFaT"
# endif
# ifndef RESULT_FLUID_ROOTFILENAME
#   define RESULT_FLUID_ROOTFILENAME "fluid"
# endif

# ifndef RESULT_DIR
#   define RESULT_DIR "Res"
# endif

# ifndef CONVERGE_UZAWA_FILENAME 
#   define CONVERGE_UZAWA_FILENAME "convergence_uzawa_velocity.dat"
# endif
# ifndef DLMFD_CELLS_FILENAME 
#   define DLMFD_CELLS_FILENAME "dlmfd_cells.dat"
# endif
# ifndef DLMFD_PERF_FILENAME 
#   define DLMFD_PERF_FILENAME "dlmfd_perf.dat"
# endif

# ifndef ROUNDDOUBLECOEF
#   define ROUNDDOUBLECOEF (1.e-4)
# endif

# ifndef TINTERVALOUTPUT
#   define TINTERVALOUTPUT (1.)
# endif

# ifndef SIMUTIMEINTERVAL
#   define SIMUTIMEINTERVAL (1.)
# endif

# ifndef FIGS_DIR
#   define FIGS_DIR "Figs"
# endif

# ifndef INITIALGRIDADAPTIVE_NEWMETHOD
#   define INITIALGRIDADAPTIVE_NEWMETHOD 1
# endif

# ifndef IMPOSED_PERIODICFLOW
#   define IMPOSED_PERIODICFLOW 0
# endif

# ifndef IMPOSED_PERIODICFLOW_TYPE
#   define IMPOSED_PERIODICFLOW_TYPE 0 // 0 for pressure and 1 for flow rate
# endif

# ifndef IMPOSED_PERIODICFLOW_DIRECTION
#   define IMPOSED_PERIODICFLOW_DIRECTION 0 
# endif

# ifndef PERIODICFLOWRATE_LEVEL
#   define PERIODICFLOWRATE_LEVEL MAXLEVEL 
# endif

# ifndef DLMFD_PROB_AFTER_NAVIERSTOKES
#   define DLMFD_PROB_AFTER_NAVIERSTOKES 0
# endif

# ifndef RIGIDBODIES_AS_FIXED_OBSTACLES
#   define RIGIDBODIES_AS_FIXED_OBSTACLES 0
# endif

# ifndef FLAG_ADAPT_CRIT
#   define FLAG_ADAPT_CRIT (1.E-9)
# endif

# ifndef UX_ADAPT_CRIT
#   define UX_ADAPT_CRIT (1.E-2)
# endif

# ifndef UY_ADAPT_CRIT
#   define UY_ADAPT_CRIT (1.E-2)
# endif

# ifndef UZ_ADAPT_CRIT
#   define UZ_ADAPT_CRIT (1.E-2)
# endif

# ifndef GRAVITY_VECTOR
#   if dimension == 2
#     define GRAVITY_VECTOR ((coord){0.,0.})
#   else 
#     define GRAVITY_VECTOR ((coord){0.,0.,0.})
#   endif
# endif

# ifndef LAMBDA2
#   define LAMBDA2 0
# endif

# ifndef VORTICITY
#   define VORTICITY 0
# endif

# ifndef PARAVIEW
#   define PARAVIEW 0
# endif

# ifndef PARAVIEW_SCALAR_LIST
#   define PARAVIEW_SCALAR_LIST p
# endif

# ifndef PARAVIEW_VECTOR_LIST
#   define PARAVIEW_VECTOR_LIST u
# endif

# ifndef BVIEW
#   define BVIEW 0
# endif

# ifndef BVIEW_LIST
#   define BVIEW_LIST u,p
# endif

# ifndef PARAVIEW_DLMFD_INTPTS
#   define PARAVIEW_DLMFD_INTPTS 0
# endif

# if PARAVIEW_DLMFD_INTPTS
#   ifndef PARAVIEW_DLMFD_INTPTS_FILENAME
#     define PARAVIEW_DLMFD_INTPTS_FILENAME "dlmfd_interior_points"
#   endif      
# endif

# ifndef PARAVIEW_DLMFD_BNDPTS
#   define PARAVIEW_DLMFD_BNDPTS 0
# endif

# if PARAVIEW_DLMFD_BNDPTS
#   ifndef PARAVIEW_DLMFD_BNDPTS_FILENAME
#     define PARAVIEW_DLMFD_BNDPTS_FILENAME "dlmfd_boundary_points"
#   endif      
# endif


double deltau;
int restarted_simu = 0;
char vtk_field_times_series[100000] = "";
# if PARAVIEW_DLMFD_BNDPTS
    char vtk_bndpts_times_series[100000] = "";      
# endif 
# if PARAVIEW_DLMFD_INTPTS
    char vtk_intpts_times_series[100000] = "";        
# endif 
int init_cycle_number = 0;
double maxtime = 0.;
double trestart = 0.;
double dtrestart = 0.;
bool save_data_restart = false;
scalar u_previoustime[];
double imposed_periodicpressuredrop = 0.;
double imposed_periodicflowrate = 0.;
size_t nbRigidBodies = 0;
size_t nbParticles = 0;
size_t NbObstacles = 0;


/** Fictitious domain implementation */
# include "DLMFD_Uzawa_velocity.h"

/** Paraview output functions */
# include "save_data_vtu.h"

/** Lambda criterion for visualizing wakes */
# include "lambda2.h"



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

# ifndef FLUID_DENSITY
#   define FLUID_DENSITY 1.
# endif
  
  /* Initialize the density field */
  const scalar rhoc[] = FLUID_DENSITY;
  rho = rhoc;

# ifndef FLUID_VISCOSITY
#   define FLUID_VISCOSITY 1.
# endif
  
  /* Initialize the viscosity field */
  const face vector muc[] = {FLUID_VISCOSITY, FLUID_VISCOSITY, FLUID_VISCOSITY};
  mu = muc;

  // Output basic fluid and geometric parameters
  if ( pid() == 0 )
  {
    printf( "Fluid density = %6.3e\n", FLUID_DENSITY );
    printf( "Fluid viscosity = %6.3e\n", FLUID_VISCOSITY );
    printf( "Space dimension = %d\n", dimension );    
    printf( "Domain size = %6.3e\n", L0 );
    printf( "\n" );            
  }  
  
  // Initialize all DLMFD fields
  initialize_DLMFD_fields_to_zero();

  // If new simulation: set fluid initial condition from user defined case file
  if ( ! restore ( file = FLUID_DUMP_FILENAME ) ) 
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
#   if Paraview
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


  // Initialize/open all DLMFD file pointers
  init_file_pointers( nbRigidBodies, pdata, !RIGIDBODIES_AS_FIXED_OBSTACLES, 
  	fdata, &converge, &cellvstime, restarted_simu );      


  // Perform initial refinement around rigid bodies and write rigid body data
  int totalcell = 0;

  if ( restarted_simu ) // Restarted simulation
  {    
    totalcell = totalcells();
    if ( pid() == 0 ) printf( "Initial total cells = %d\n", totalcell );
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
     
        // Create rigid body boundary points
        create_rigidbodies_boundary_points( allRigidBodies, nbRigidBodies );
      
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
		      	allRigidBodies[k].s.x[j] + ni * delta,
          		allRigidBodies[k].s.y[j] + nj * delta,
			allRigidBodies[k].s.z[j] + nk * delta );
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
                  dist.x = x - allRigidBodies[k].s.x[j];
                  dist.y = y - allRigidBodies[k].s.y[j];
#                 if dimension > 2
                    dist.z = z - allRigidBodies[k].s.z[j];
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
          ss = adapt_wavelet( {DLM_Flag}, (double[]) {1.e-30}, 
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
      
        // No need to free the rigid body boundary points, there will be freed
        // in the once_timestep_is_determined event by the free_rigidbodies 
	// function
      
#     endif
#   endif
    
    // Save all rigid body data
    if ( !RIGIDBODIES_AS_FIXED_OBSTACLES )
      rigidbody_data( allRigidBodies, nbRigidBodies, t, i, pdata );  
  }

  
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
# if BRINKMANN_DIRICHLET_PENALIZATION
    Xi.nodump = true;   
# endif  
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
#   if BRINKMANN_DIRICHLET_PENALIZATION
      Usolid.x.nodump = true ;   
#   endif   
  }
  
  
  /* Construction of rigid bodies and their DLMFD features */
  if ( pid() == 0 ) printf( "Initial DLMFD RB construction\n" ); 
  DLMFD_construction();    
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
    close_file_pointers( nbRigidBodies, pdata, fdata, converge, cellvstime ); 

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
    vector omega[];
#   if dimension == 2
      vorticity( u, omega );
#   else
      vorticity_3D( u, omega ); 
#   endif            
# endif

  // Paraview output for post-processing
# if PARAVIEW
    scalar* paraview_scalarlist = {PARAVIEW_SCALAR_LIST
#   if LAMBDA2
      , l2
#   endif         
    };
    vector* paraview_vectorlist = {PARAVIEW_VECTOR_LIST
#   if LAMBDA2
      , omega
#   endif         
    };    
    save_data( paraview_scalarlist, paraview_vectorlist, allRigidBodies, 
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
#   if LAMBDA2
      , omega
#   endif        
    };
    dump( file = dump_name, dump_list );  
# endif

  // Basilisk output for restart
# if DLM_alpha_coupling
    dump_list = (scalar *){u, g, DLM_explicit};
# else
    dump_list = (scalar *){u, g};    
# endif     
  dump( file = FLUID_DUMP_FILENAME, dump_list );    
   
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
    free_rigidbodies( allRigidBodies, nbRigidBodies );
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
    free_rigidbodies( allRigidBodies, nbRigidBodies );
    DLMFD_construction();  
# endif 

  do_output( mess );    
  output_dlmfd_perf ( dlmfd_globaltiming, i, allRigidBodies );       
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
  free_rigidbodies( allRigidBodies, nbRigidBodies );  
  
  // In case of a periodic flow, we add the imposed pressure drop
# if IMPOSED_PERIODICFLOW
#   if IMPOSED_PERIODICFLOW_DIRECTION == 0 
      const face vector dp[] = { 
		imposed_periodicpressuredrop / ( L0 * FLUID_DENSITY ), 0., 0. };  
#   elif IMPOSED_PERIODICFLOW_DIRECTION == 1
      const face vector dp[] = { 0.,
		imposed_periodicpressuredrop / ( L0 * FLUID_DENSITY ), 0. };  	
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
  free_rigidbodies( allRigidBodies, nbRigidBodies ); 
  free_np_dep_arrays( nbParticles, allRigidBodies, DLMFDtoGS_vel, vpartbuf, 
  	pdata, !RIGIDBODIES_AS_FIXED_OBSTACLES, fdata );
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
        flowrate = compute_flowrate_right( u, periodicflowrate_level );
#     elif IMPOSED_PERIODICFLOW_DIRECTION == 1
        flowrate = compute_flowrate_top( u, periodicflowrate_level );
#     else 
        flowrate = compute_flowrate_front( u, periodicflowrate_level );
#     endif 
      if ( pid() == 0 )
        printf( "   Periodic flow rate = %8.5e\n", flowrate );           
#   else
      double Q1 = - dt / ( L0 * FLUID_DENSITY );
      double deltaflowrate = 0.;                    
#     if IMPOSED_PERIODICFLOW_DIRECTION == 0 
        flowrate = compute_flowrate_right( u, PERIODICFLOWRATE_LEVEL );
	deltaflowrate = imposed_periodicflowrate - flowrate;
	imposed_periodicpressuredrop += deltaflowrate / ( Q1 * L0 * L0 );
	foreach()
	  u.x[] += deltaflowrate / ( L0 * L0 );  
#     elif IMPOSED_PERIODICFLOW_DIRECTION == 1
        flowrate = compute_flowrate_top( u, PERIODICFLOWRATE_LEVEL );
	deltaflowrate = imposed_periodicflowrate - flowrate;
	imposed_periodicpressuredrop += deltaflowrate / ( Q1 * L0 * L0 );
	foreach()
	  u.y[] += deltaflowrate / ( L0 * L0 ); 
#     else 
        flowrate = compute_flowrate_front( u, PERIODICFLOWRATE_LEVEL );
	deltaflowrate = imposed_periodicflowrate - flowrate;
	
//         // PID controller
//         double Kp = 1.e-1 / ( Q1 * L0 * L0 ) ;
//         double Ki = 2. * Kp / dt ;
//         double Kd = Kp * dt / 8.;
//       
//         static double deltaflowrate_nm1 = 0.;
//         static double deltaflowrate_nm2 = 0.;
//         static double imposed_periodicpressuredrop_nm1 = 0.;
//         static double dt_nm1 = 1.e-10; 
// 
// 	imposed_periodicpressuredrop = imposed_periodicpressuredrop_nm1 
// 		+ ( Kp + Ki * dt + Kd / dt ) * deltaflowrate
//       		- ( Kp + Kd * ( 1. / dt_nm1 + 1. / dt ) ) 
// 			* deltaflowrate_nm1 
// 		+ ( Kd / dt_nm1 ) * deltaflowrate_nm2;
// 
//         deltaflowrate_nm2 = deltaflowrate_nm1;
//         deltaflowrate_nm1 = deltaflowrate;
//         imposed_periodicpressuredrop_nm1 = imposed_periodicpressuredrop;
// 	dt_nm1 = dt;	
		
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

    astats s = adapt_wavelet( (scalar *){DLM_FlagMesh, u}, 
	(double[]){FLAG_ADAPT_CRIT, UX_ADAPT_CRIT, UY_ADAPT_CRIT, 
	UZ_ADAPT_CRIT}, maxlevel = MAXLEVEL, minlevel = LEVEL );		

    if ( pid() == 0 ) 
      printf( "refined = %d, coarsened = %d\n", s.nf, s.nc );
# endif
}
