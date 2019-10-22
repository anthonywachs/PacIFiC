/* File names definition and global variables */  
# ifndef fluid_dump_filename
# define fluid_dump_filename "Savings/dump"
# endif
# ifndef particle_dump_filename
# define particle_dump_filename "Savings/dump_particle"
# endif
# ifndef dump_dir
# define dump_dir "Savings"
# endif

# ifndef result_dir
# define result_dir "Res"
# endif
# ifndef result_particle_vp_rootfilename
# define result_particle_vp_rootfilename "particle_data"
# endif
# ifndef result_particle_hydroFaT_rootfilename
# define result_particle_hydroFaT_rootfilename "particle_hydroFaT"
# endif
# ifndef result_fluid_rootfilename
# define result_fluid_rootfilename "fluid_basilisk"
# endif

# ifndef result_dir
# define result_dir "Res"
# endif

# ifndef converge_uzawa_filename 
# define converge_uzawa_filename "converge_uzawa.dat"
# endif
# ifndef dlmfd_cells_filename 
# define dlmfd_cells_filename "dlmfd_cells.dat"
# endif
# ifndef dlmfd_perf_filename 
# define dlmfd_perf_filename "dlmfd_perf.dat"
# endif

# ifndef RoundDoubleCoef
# define RoundDoubleCoef (1.e-4)
# endif


double deltau;
int restarted_simu = 0;
char vtk_times_series[100000] = "";
int init_cycle_number = 0;
double maxtime = 0.;
scalar u_previoustime[];


/* Fictitious domain implementation */
# include "DLMFD_reverse_Uzawa.h"


/** Generic granular solver init event: TO BE OVERLOADED */
// --------------------------------------------------------
event GranularSolver_init (t < -1.)
{} 


/** Generic granular solver predictor event: TO BE OVERLOADED */
// -------------------------------------------------------------
event GranularSolver_predictor (t < -1.)
{} 


/** Generic granular solver velocity update event: TO BE OVERLOADED */
// -------------------------------------------------------------------
event GranularSolver_updateVelocity (t < -1.)
{}


/** Overloading of the init event: initialize fluid and particles */
// -----------------------------------------------------------------
event init (i = 0) {
#ifndef rhoval
# define rhoval 1.
#endif
  
  /* Initialize the density field */
  const scalar rhoc[] = rhoval;
  rho = rhoc;

#ifndef tval
# define tval 1.
#endif
  
  /* Initialize the viscosity field */
  const face vector muc[] = {tval, tval, tval};
  mu = muc;

  // If new simulation: set fluid initial condition and initialise data file 
  // pointers 
  if ( ! restore ( file = fluid_dump_filename ) ) 
  {
    // Set the restarted simulation boolean to 1
    restarted_simu = 0;
    
    // Initialise data file pointers
    init_file_pointers( pdata, fdata, 0 );

    /* Generic particle parameters, will be overwritten by the
       GranularSolver_init event if need be */
    for (int k = 0; k < NPARTICLES; k++) 
    {
      /* particle id */
      particles[k].pnum = k;
      
# if DLM_Moving_particle    
      /* density rho_s of the particle */
# ifndef rhosolid
# define rhosolid  HUGE
# endif
      particles[k].rho_s = rhosolid;
	
      /* DLMFD coupling factor */
      /* If explicit add mass, DLMFD_couplingFactor = 1 */
      /* otherwise DLMFD_couplingFactor = ( 1 - rhoval / rho_s ) */
    
      particles[k].DLMFD_couplingfactor = 1. ;

#if (!b_explicit_added_mass) 
      particles[k].DLMFD_couplingfactor -= rhoval / particles[k].rho_s ;
#endif
      
      if ( !particles[k].iswall && !particles[k].iscube ) 
      {
	GeomParameter gp = particles[k].g;
# if dimension == 2 
	/* Surface of the particle for a circle  */
	particles[k].Vp = pi*sq(gp.radius);

	/* total weight of the particle */
	particles[k].M = rhosolid*(particles[k].Vp);

	/* The inertia tensor is: */
	/* For a solid disk: */
	particles[k].Ip[0] = (particles[k].M)*sq(gp.radius)/4;  /* Ip[0] = Ixx */
	particles[k].Ip[1] = particles[k].Ip[0];                /* Ip[1] = Iyy */
	particles[k].Ip[2] = (particles[k].M)*sq(gp.radius)/2;  /* Ip[2] = Izz */
	particles[k].Ip[3] = 0.;                                /* Ip[3] = Ixy */
	particles[k].Ip[4] = 0.;                                /* Ip[4] = Ixz */
	particles[k].Ip[5] = 0.;                                /* Ip[5] = Iyz */
# endif
# if dimension == 3
	/* Volume  of the particle (sphere) */
	particles[k].Vp = 4.*pi*pow(gp.radius,3)/3.;

	/* total weight of the particle */
	particles[k].M = rhosolid*(particles[k].Vp);

	/* For a sphere it comes: I_xx = I_yy = I_zz = 2/5 M R^2 with R the radius */
	particles[k].Ip[0] = 2.*(particles[k].M)*sq(gp.radius)/5.; /* Ip[0] = Ixx */
	particles[k].Ip[1] = particles[k].Ip[0];                   /* Ip[1] = Iyy */
	particles[k].Ip[2] = particles[k].Ip[0];                   /* Ip[2] = Izz */
	particles[k].Ip[3] = 0.;                                   /* Ip[3] = Ixy */
	particles[k].Ip[4] = 0.;                                   /* Ip[4] = Ixz */
	particles[k].Ip[5] = 0.;                                   /* Ip[5] = Iyz */
# endif
      }
    
      if (particles[k].iscube) { /* To be done if needed */ }
#endif
    }
  }
  else // Restart of a simulation 
  {
      /* Set the restarted simulation boolean to 1 */
      restarted_simu = 1;
      
      /* Restore particle data */
      restore_particle( particles );

      /* Re-initialize file pointers */
      init_file_pointers( pdata, fdata, 1 );

      // Re-initialize the VTK writer
#if Paraview
      reinitialize_vtk_restart();	
#endif	
  }
  
  // Initialize the granular solver
  if ( pid() == 0 ) printf ("# granular solver initialization: "); 
  event ("GranularSolver_init");
      
  // Special case of cubes
  for (int k = 0; k < NPARTICLES; k++) 
  {
    GeomParameter * gg;
    gg = &(particles[k].g);
    if ( gg->ncorners == 8 ) 
    {
	particles[k].iscube = 1;
	compute_principal_vectors_Cubes( &(particles[k]) );
	printf ("length cube = %f\n", gg->radius);
    }
  }

  // Perform initial refinement around particles and write particle data
  int totalcell = 0;
  
  // Restarted simulation
  if ( restarted_simu )
  {    
    totalcell = totalcells();
    if ( pid() == 0 ) printf( "Initial total cells = %d\n", totalcell );
  }
  else // Simulation from t=0
  {
# if adaptive
    for (int k = 0; k < NPARTICLES; k++) 
    {
      GeomParameter * gg;
      gg = &(particles[k].g);    	  
    
      /* Perform initial refinement */
      totalcell = totalcells();
      if ( pid() == 0 )
	printf("# total cells before initial refinement = %d\n", totalcell);

      coord position = gg->center;
      double radius = gg->radius;
      Point lpoint;
	
      /* Initial static refinement */
      lpoint = locate( position.x, position.y, position.z );
      int ilvl = 0;
      while( ilvl < MAXLEVEL ) 
      {
	printf ("current level of refinement = %d, we want level %d\n",
	  	ilvl, MAXLEVEL );
	printf ("refining ... \n");
	if ( lpoint.level > -1 ) 
	{
	  refine_cell (lpoint, all, 0, &tree->refined);
	  printf ("refined once on thread %d \n", pid());
	}
	mpi_boundary_refine( all );
	mpi_boundary_update( all );
      
	lpoint = locate (position.x, position.y, position.z);
	ilvl = lpoint.level;

#if _MPI
	mpi_all_reduce (ilvl, MPI_INT, MPI_MAX);
#endif
      }
  
      radius = 1.2 * radius ;
      refine (sq(x - position.x) + sq(y - position.y) + sq(z - position.z) 
	      - sq(radius) < 0 && level < MAXLEVEL);
	     
      // Compute the total number of cells after refinement and 
      totalcell = totalcells();
      if ( pid() == 0 )
	printf ("# total cells after initial refinement = %d\n", totalcell);
    }
#endif
    
    // Save all particles trajectories
    particle_data( particles, t, i, pdata );  
  }

  // Assign gravity to particles
# ifndef gravity_x
# define gravity_x 0.
# endif
# ifndef gravity_y
# define gravity_y 0.
# endif
# ifndef gravity_z
# define gravity_z 0.
# endif
  for (int k = 0; k < NPARTICLES; k++) 
  {    
    particles[k].gravity.x = gravity_x;
    particles[k].gravity.y = gravity_y;
    particles[k].gravity.z = gravity_z;
  }    
  
  // Simulation time interval
  maxtime = trestart + SimuTimeInterval;
  fprintf( ferr, "Simulation time interval: t_start = %f to t_end = %f\n", 
    	trestart, maxtime ); 
	
  // Initialize the field u_previoustime to compute x-velocity change
  foreach() { u_previoustime[] = u.x[]; }	   
}




/** Overloading of the logfile event: we log the number of iterations of the 
multigrid solver for pressure and viscosity and stop the simulation at 
maxtime */
// -------------------------------------------------------------------------
event logfile ( i=0; i++ ) 
{
  deltau = change( u.x, u_previoustime );
  fprintf ( stderr, "%d %g %d %d %g %g \n", i, t, mgp.i, mgu.i, deltau, dt );

  // This is the condition that stops the simulation exactly at t=maxtime
  // We use -0.0001 * dt to avoid the problem of comparison of double that 
  // would exist if we would write if ( t > maxtime )
  if ( t - maxtime > - RoundDoubleCoef * dt ) return 1; 
//  if ( t > maxtime ) return 1; 
}




/** Overloading of the viscous_term event: we add an explicit coupling term
that improves the coupling between the fluid problem and the DLMFD problem */
// -------------------------------------------------------------------------- 
#if DLM_alpha_coupling
event viscous_term (i++) 
{
  foreach() {
    foreach_dimension() {
      u.x[] += -dt*DLM_explicit.x[]/(rhoval*dlmfd_dv());
    }
  }
  boundary((scalar*){u});
}
#endif




/** Overloading of the end_timestep event: we plug the solution to the 
granular problem followed by the solution to the DLMFD problem */
// -------------------------------------------------------------------------
event end_timestep (i++) 
{
  /* Granular solver predictor */
  if ( pid() == 0 ) printf( "Iter %d: ", i ); 
  event ("GranularSolver_predictor");  
 
  /* Solve the DLMFD problem by a Uzawa algorithm */
  DLMFD_subproblem( particles, i, rhoval );

  /* Update velocity in the granular solver */
  if ( pid() == 0 ) printf( "Iter %d: ", i );
  event ("GranularSolver_updateVelocity");    

  /* Save the forces acting on particles before adapting the mesh  */
  sumLambda( particles, fdata, t + dt, dt, flagfield, DLM_lambda, index_lambda, 
	rhoval, DLM_periodic_shift );

  /* Free particle's dlmfd points (we dont need them anymore) */
  free_particles( particles, NPARTICLES );
  
  /* Save all particles trajectories */
  particle_data( particles, t + dt, i, pdata );
}




/** Overloading of the adapt event: we refine the mesh with both the velocity
field and a phase indicator */
// -------------------------------------------------------------------------
# ifndef FlagAdaptCrit
# define FlagAdaptCrit (1.E-9)
# endif

# ifndef UxAdaptCrit
# define UxAdaptCrit (1.E-2)
# endif

# ifndef UyAdaptCrit
# define UyAdaptCrit (1.E-2)
# endif

# ifndef UzAdaptCrit
# define UzAdaptCrit (1.E-2)
# endif

event adapt (i++) {
  int totalcell = totalcells();
  if (pid() == 0)
    fprintf (stderr, "# Total number of cells: %d \n", totalcell);

#if adaptive
#if DLM_Moving_particle
  astats s = adapt_wavelet ((scalar *){flagfield_mailleur, u}, 
			    (double[]){FlagAdaptCrit, UxAdaptCrit, UyAdaptCrit, UzAdaptCrit}, 
			    maxlevel = MAXLEVEL);
#else
  astats s = adapt_wavelet ((scalar *){flagfield, u}, 
			    (double[]){FlagAdaptCrit, UxAdaptCrit, UyAdaptCrit, UzAdaptCrit}, 
			    maxlevel = MAXLEVEL);
#endif
  fprintf (ferr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
#endif
}
