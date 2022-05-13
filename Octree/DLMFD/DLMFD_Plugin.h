/** 
# The general DLMFD plugin 
*/

/** File names definition and global variables */  
# ifndef fluid_dump_filename
#   define fluid_dump_filename "Savings/dump"
# endif
# ifndef dump_dir
#   define dump_dir "Savings"
# endif

# ifndef result_dir
#   define result_dir "Res"
# endif
# ifndef result_particle_vp_rootfilename
#   define result_particle_vp_rootfilename "particle_data"
# endif
# ifndef result_particle_hydroFaT_rootfilename
#   define result_particle_hydroFaT_rootfilename "particle_hydroFaT"
# endif
# ifndef result_fluid_rootfilename
#   define result_fluid_rootfilename "fluid_basilisk"
# endif

# ifndef result_dir
#   define result_dir "Res"
# endif

# ifndef converge_uzawa_filename 
#   define converge_uzawa_filename "converge_uzawa.dat"
# endif
# ifndef dlmfd_cells_filename 
#   define dlmfd_cells_filename "dlmfd_cells.dat"
# endif
# ifndef dlmfd_perf_filename 
#   define dlmfd_perf_filename "dlmfd_perf.dat"
# endif

# ifndef RoundDoubleCoef
#   define RoundDoubleCoef (1.e-4)
# endif


double deltau;
int restarted_simu = 0;
char vtk_times_series[100000] = "";
int init_cycle_number = 0;
double maxtime = 0.;
double trestart = 0.;
double dtrestart = 0.;
bool save_data_restart = false;
scalar u_previoustime[];


/** Fictitious domain implementation */
# include "DLMFD_Uzawa_velocity.h"

/** Paraview output functions */
# include "save_data_vtu.h"




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




/** Overloading of the init event: initialize fluid and particles */
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

# ifndef rhoval
#   define rhoval 1.
# endif
  
  /* Initialize the density field */
  const scalar rhoc[] = rhoval;
  rho = rhoc;

# ifndef tval
#   define tval 1.
# endif
  
  /* Initialize the viscosity field */
  const face vector muc[] = {tval, tval, tval};
  mu = muc;

  // Output basic fluid and geometric parameters
  if ( pid() == 0 )
  {
    printf( "Fluid density = %6.3e\n", rhoval );
    printf( "Fluid viscosity = %6.3e\n", tval );
    printf( "Space dimension = %d\n", dimension );    
    printf( "Domain size = %6.3e\n", L0 );
    printf( "\n" );            
  }  
  
  // Initialize all DLMFD fields
  initialize_DLMFD_fields_to_zero();

  // If new simulation: set fluid initial condition from user defined case file
  if ( ! restore ( file = fluid_dump_filename ) ) 
  {
    // Set the restarted simulation boolean to 1
    restarted_simu = 0;
  }
  else // Restart of a simulation 
  {
    // Set the restarted simulation boolean to 1
    restarted_simu = 1;

    // Read restart time
    if ( pid() == 0 ) printf( "Read t and dt from time restart file\n" );
    read_t_restart( dump_dir, &trestart, &dtrestart );  

    // Re-initialize the VTK writer
#   if Paraview
      reinitialize_vtk_restart();	
#   endif	
  }


  // Initialize/open all DLMFD file pointers
  init_file_pointers( pdata, fdata, &converge, &cellvstime, restarted_simu );

  
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


  // Perform initial refinement around particles and write particle data
  int totalcell = 0;

  if ( restarted_simu ) // Restarted simulation
  {    
    totalcell = totalcells();
    if ( pid() == 0 ) printf( "Initial total cells = %d\n", totalcell );
  }
  else // Simulation from t=0
  {
#   if adaptive
      if ( pid() == 0 )
      {
        printf( "\nInitial grid refinement in the particle volume\n" );
        printf( "----------------------------------------------\n" );      
      }
      
      for (int k = 0; k < NPARTICLES; k++) 
      {
        GeomParameter * gg;
        gg = &(particles[k].g);    	  
    
        /* Perform initial refinement */
        totalcell = totalcells();
        if ( pid() == 0 )
	  printf( "# total cells before initial refinement of "
	  	"particle %d = %d\n", k, totalcell );

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
      
	  lpoint = locate (position.x, position.y, position.z);
	  ilvl = lpoint.level;

#         if _MPI
	    mpi_all_reduce( ilvl, MPI_INT, MPI_MAX );
#         endif
        }
  
        radius = 1.2 * radius ;
        refine( sq( x - position.x ) + sq( y - position.y ) 
		+ sq( z - position.z ) - sq(radius) < 0. && level < MAXLEVEL );
	     
        // Compute the total number of cells after refinement and 
        totalcell = totalcells();
        if ( pid() == 0 )
	  printf( "# total cells after initial refinement of particle %d = "
		"%d\n\n", k, totalcell );
      }
#   endif
    
    // Save all particles trajectories
    particle_data( particles, t, i, pdata );  
  }


  // Assign gravity to particles
# if DLM_Moving_particle
#   ifndef gravity_x
#     define gravity_x 0.
#   endif
#   ifndef gravity_y
#     define gravity_y 0.
#   endif
#   ifndef gravity_z
#     define gravity_z 0.
#   endif
    for (int k = 0; k < NPARTICLES; k++) 
    {    
      particles[k].gravity.x = gravity_x;
      particles[k].gravity.y = gravity_y;
      particles[k].gravity.z = gravity_z;
    }    
# endif
  
  // Simulation time interval
  maxtime = trestart + SimuTimeInterval;
  if ( pid() == 0 ) 
    printf( "Simulation time interval: t_start = %f to t_end = %f\n\n", 
    	trestart, maxtime ); 

	
  // Initialize the field u_previoustime to compute x-velocity change
  foreach() u_previoustime[] = u.x[];  
  
  
  // By default:
  // * do NOT DUMP the work fields used in the DLMFD sub-problem and the 
  // field u_previoustime
  // * DUMP the flagfield and DLM_explicit fields in the DLMFD sub-problem
  u_previoustime.nodump = true ;
  flagfield.nodump = false;
  flagfield_mailleur.nodump = true ;
  foreach_dimension()
  {
    DLM_lambda.x.nodump = true ;
    index_lambda.x.nodump = true ;
    DLM_periodic_shift.x.nodump = true ;
    DLM_r.x.nodump = true ;
    DLM_w.x.nodump = true ;
    DLM_v.x.nodump = true ;
    qu.x.nodump = true ;
    tu.x.nodump = true ;
#   if DLM_alpha_coupling
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
  if ( t - maxtime > - RoundDoubleCoef * dt ) 
  {
    // Close all DLMFD files
    close_file_pointers( pdata, fdata, converge, cellvstime ); 

    // Write the dump time and time step for restart
    if ( pid() == 0 ) printf( "Write t and dt in time restart file\n" );
    save_t_dt_restart( dump_dir, t, dt );  
    
    // Stop simulation
    return 1; 
  }
//  if ( t > maxtime ) return 1;   
}




/** Writes the dump time and time step for restart, and frees dynamic features 
of particles */
//----------------------------------------------------------------------------
event start_timestep (i++)
//----------------------------------------------------------------------------
{
  if ( save_data_restart && i )
  {
    if ( pid() == 0 ) printf( "Write t and dt in time restart file\n" );
    save_t_dt_restart( dump_dir, t, dt );        
    save_data_restart = false;
  } 
  if ( i == 0 ) save_data_restart = false;
  
  // We free dynamic features of particles from the previous time step
  // or from initialization at the start of the current time step 
  // (run with -events to understand)
  free_particles( particles, NPARTICLES ); 
}




/** Frees dynamic features of particles at the very end */
//----------------------------------------------------------------------------
event cleanup (t = end)
//----------------------------------------------------------------------------
{
  // This is because we free dynamic features of particles from the previous 
  // time step at the start of the current time step (see above start_timestep)
  // Since at the very end, there is no next time step, dynamic features 
  // of particles are not freed by start_timestep and hence are freed here
  free_particles( particles, NPARTICLES ); 
}




/** Overloading of the viscous_term event: we add an explicit coupling term
that improves the coupling between the fluid problem and the DLMFD problem */
//----------------------------------------------------------------------------
event viscous_term (i++) 
//----------------------------------------------------------------------------
{
# if DLM_alpha_coupling
    foreach()
      foreach_dimension()
        u.x[] += -dt*DLM_explicit.x[]/(rhoval*dlmfd_dv());

    boundary((scalar*){u});
# endif
}





/** Overloading of the end_timestep event: we plug the solution to the 
granular problem followed by the solution to the DLMFD problem */
//----------------------------------------------------------------------------
event end_timestep (i++) 
//----------------------------------------------------------------------------
{
  if ( pid() == 0 )
  {
    printf( "\n****** ITER %d - TIME t=%8.5e to t+dt=%8.5e ******\n", 
    	i, t, t+dt );
    printf( "   Time step = %8.5e\n", dt );
    printf( "   Navier-Stokes solver\n" ); 
    printf( "      MGu_niter = %d, MGp_niter = %d\n", mgu.i,
    	mgp.i );    
  }

  /* Granular solver predictor */
  if ( pid() == 0 ) printf( "   GS predictor step: " ); 
  event( "GranularSolver_predictor" );  
 
  /* Solve the DLMFD problem by a Uzawa algorithm */
  if ( pid() == 0 ) printf( "   DLMFD problem\n" ); 
  DLMFD_subproblem( particles, i, rhoval );

  /* Update velocity in the granular solver */
  if ( pid() == 0 ) printf( "   GS velocity update in " );
  event( "GranularSolver_updateVelocity" );    

  /* Save the forces acting on particles before adapting the mesh  */
  sumLambda( particles, fdata, t + dt, dt, flagfield, DLM_lambda, index_lambda, 
	rhoval, DLM_periodic_shift );
  
  /* Save all particles trajectories */
  particle_data( particles, t + dt, i, pdata );
  
  /* Fluid velocity change over the time step */
  deltau = change( u.x, u_previoustime );
  if ( pid() == 0 )
    printf( "   Velocity change = %8.5e\n", deltau );  
}




/** Overloading of the adapt event: we refine the mesh with both the velocity
field and a phase indicator */
# ifndef FlagAdaptCrit
#   define FlagAdaptCrit (1.E-9)
# endif

# ifndef UxAdaptCrit
#   define UxAdaptCrit (1.E-2)
# endif

# ifndef UyAdaptCrit
#   define UyAdaptCrit (1.E-2)
# endif

# ifndef UzAdaptCrit
#   define UzAdaptCrit (1.E-2)
# endif
//----------------------------------------------------------------------------
event adapt (i++) 
//----------------------------------------------------------------------------
{
# if adaptive
    int totalcell = totalcells();
    if ( pid() == 0 )
    {
#     if dimension == 2  
        printf( "   Quadtree grid\n" );
#     else
        printf( "   Octree grid\n" );
#     endif
      printf( "      Total number of cells = %d \n", totalcell );
    }

#   if DLM_Moving_particle
      astats s = adapt_wavelet ((scalar *){flagfield_mailleur, u}, 
	(double[]){FlagAdaptCrit, UxAdaptCrit, UyAdaptCrit, UzAdaptCrit}, 
	maxlevel = MAXLEVEL);
#   else
      astats s = adapt_wavelet ((scalar *){flagfield, u}, 
	(double[]){FlagAdaptCrit, UxAdaptCrit, UyAdaptCrit, UzAdaptCrit}, 
	maxlevel = MAXLEVEL, minlevel=LEVEL );
#   endif
    if ( pid() == 0 ) 
      printf( "      Refined %d cells, coarsened %d cells\n", s.nf, s.nc );
# endif
}
