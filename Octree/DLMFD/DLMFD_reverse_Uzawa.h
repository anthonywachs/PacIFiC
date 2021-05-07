/**
# Fictitious-domain method with distributed Lagrangian multipliers 

We are looking for a numerical approximation for the solution of the
following equations 
$$
\rho\left(\partial_t{\mathbf{u}} 
+ \mathbf{u}\cdot \mathbf{\nabla}\mathbf{u}\right) = 
-\mathbf{\nabla} p + \mathbf{\nabla} \cdot\left(2\mu\mathbf{D}\right) 
+ \rho\mathbf{a} + \mathbf{\lambda}~\text{over}~\Omega
$$
$$
\left(1-\frac{\rho}{\rho_s}\right)
\left(M\left(\partial_t{\mathbf{U}}-\mathbf{g}\right)\right)=
-\int_{P(t)} {\mathbf{\lambda}} dv~\text{over}~P(t)
$$
$$
\left(1-\frac{\rho}{\rho_s}\right)
\left(\mathbf{I}\partial_t{\mathbf{\omega}} 
+ \mathbf{\omega} \times\mathbf{I}\mathbf{\omega}\right) = 
-\int_{P(t)}\mathbf{r}\times{\mathbf{\lambda}} dv~\text{over}~P(t)
$$
$$
\mathbf{u}-\left(\mathbf{U}+\mathbf{\omega}\times \mathbf{r}\right)=
0~\text{over}~P(t)
$$
$$
\mathbf{\nabla}\cdot\mathbf{u} = \mathbf{0} ~\text{over}~\Omega
$$

From a numerical point of view, this is a complicated set of equations
to solve directly. One has to deal with three main difficulties:

  * an advection-diffusion equation
  * the imcompressibility constraint with the pressure $p$ as unknown 
  * the particle's solid body's constraint with the Lagrange
multipliers as unknow.

The classic strategy is to split the problem into subproblems and
solve them successively. We chose here a two-steps time-spliting.
*/

# define BGHOSTS 2
# define BSIZE 128

# ifndef DLM_Moving_particle 
#   define DLM_Moving_particle 0
# endif

# if ( ! DLM_Moving_particle ) 
#   ifdef TRANSLATION
#     undef TRANSLATION
#   endif
#   define TRANSLATION 0
#   ifdef ROTATION
#     undef ROTATION
#   endif
#   define ROTATION 0
# else
#   ifdef TRANSLATION
#     undef TRANSLATION
#   endif
#   define TRANSLATION 1
#   ifdef ROTATION
#     undef ROTATION
#   endif
#   define ROTATION 1
# endif

# ifndef NPARTICLES 
#   define NPARTICLES 1
# endif

# ifndef DLM_alpha_coupling
#   define DLM_alpha_coupling 0
# endif

# ifndef debugBD
#   define debugBD 0     // set 1 to desactivate boundary points
# endif
# ifndef debugInterior
#   define debugInterior 0  // set 1 to desactivate interior points
# endif

# ifndef DLMFD_OPT
#   define DLMFD_OPT 1     // use optimized version of DLMFD below
# endif

# ifndef PARTICLE_VERBOSE
#   define PARTICLE_VERBOSE 0     // print particle features
# endif


/** Functions and structures needed for the implementation of
    fictitious domain method. */
# include "dlmfd_functions.h"


/** Basilisk scalars and vectors needed for the implementation of the
    fictitious domain method. */
vector DLM_lambda[];
scalar flagfield[];
scalar flagfield_mailleur[];
vector index_lambda[];
vector DLM_periodic_shift[];
vector DLM_r[];
vector DLM_w[];
vector DLM_v[];
vector qu[];
vector tu[];
# if DLM_alpha_coupling
    vector DLM_explicit[];
# endif

particle particles[NPARTICLES] = {{0}};
double DLMFDtoGS_vel[NPARTICLES][6] = {{0}};
# if ( TRANSLATION && ROTATION )
#   if dimension == 3
#     define npartdata 6
#   else
#     define npartdata 3
#   endif     
# elif TRANSLATION
#   if dimension == 3 
#     define npartdata 3
#   else
#     define npartdata 2
#   endif 
# elif ROTATION
#   if dimension == 3 
#     define npartdata 3
#   else
#     define npartdata 1
#   endif 
# else
#   define npartdata 0  
# endif
double vpartbuf[npartdata*NPARTICLES];

FILE* pdata[NPARTICLES];
FILE* fdata[NPARTICLES];
FILE* converge = NULL;
FILE* cellvstime = NULL;
char converge_uzawa_filename_complete_name[80];
char dlmfd_cells_filename_complete_name[80]; 


/**
## First sub-problem: Navier Stokes problem

We are using the centred solver for this problem.
$$
\rho\left(\partial_t{\mathbf{u}}
+\mathbf{u}\cdot {\mathbf{\nabla}}\mathbf{u}\right) = 
-{\mathbf{\nabla}}p + {\mathbf{\nabla}}\cdot\left(2\mu\mathbf{D}\right) 
+ \rho\mathbf{a}~\text{over}~\Omega
$$
$$
{\mathbf{\nabla}}\cdot\mathbf{u} = \mathbf{0} ~\text{over}~\Omega. 
$$
*/

# include "dlmfd-ns-centered.h"


/* Adding this small macro because dv() breaks the compability with embed.h */
# if dimension == 1
#   define dlmfd_dv() (Delta)
# elif dimension == 2
#   define dlmfd_dv() (sq(Delta))
# else // dimension == 3
#   define dlmfd_dv() (cube(Delta))
# endif


/**
## Second sub-problem: Fictitious-domain problem
We are solving a fictitious problem given by:
$$
\rho \partial_t\mathbf{u} = {\mathbf{\lambda}}~\text{over}~\Omega
$$
$$
\left(1-\frac{\rho}{\rho_s}\right)
\left(M\left(\partial_t\mathbf{U}-\mathbf{g}\right)\right)=
-\int_{P(t)}{\mathbf{\lambda}} dv~\text{over}~P(t)
$$
$$
\left(1-\frac{\rho}{\rho_s}\right)
\left(\mathbf{I}\partial_t {\mathbf{\omega}} 
+ {\mathbf{\omega}} \times\mathbf{I}{\mathbf{\omega}}\right) = 
-\int_{P(t)}\mathbf{r}\times{\mathbf{\lambda}} dv~\text{over}~P(t)
$$
$$
\mathbf{u}-\left(\mathbf{U}+{\mathbf{\omega}}\times \mathbf{r}\right)=
0~\text{over}~P(t),
$$
with unknowns $\mathbf{u}, \mathbf{U}, {\mathbf{\omega}}$ and 
${\mathbf{\lambda}}$ being respectively the fluid velocity field, 
the particle's translational velocity, the particle's rotational velocity and 
the Lagrange multipliers. The particle occupies the domain $P(t)$ with density 
$\rho_s$, Inertia tensor $\mathbf{I}$ and mass $M$. The vector $\mathbf{g}$ is 
the gravity acceleration here.

This leads to a saddle-point problem which is solved with an iterative solver 
(Uzawa, conjugate gradient algorithm). It is implemented in the function below. 
*/


/* Timers and Timings */
timing dlmfd_globaltiming = {0.};

# include "dlmfd_perf.h"
# if DLMFD_OPT
#   include "dlmfd-fast.h"
# endif 




// Uzawa algorithm for the DLMFD problem
// Important comments:
// 1) only u, tu, DLM_lambda and DLM_w need to be synchronized with boundary, 
// all other fields are local to each thread
//----------------------------------------------------------------------------
void DLMFD_subproblem( particle * p, const int i, const double rho_f ) 
//----------------------------------------------------------------------------
{
  /* Timers and Timings */
  timer dlmfd_timer = timer_start ();
  double mpitimings[npe()];
  
  double DLM_alpha = 0., DLM_beta = 0.;
  double DLM_tol = 1.e-5, DLM_nr2 = 0., DLM_nr2_km1 = 0., DLM_wv = 0.;
  int ki = 0, DLM_maxiter = 200, allpts = 0, lm = 0, tcells = 0;
  coord ppshift = {0, 0, 0};
# if PARTICLE_VERBOSE
    char outputshift[6]="      ";
# endif
# if ( _MPI && DLM_Moving_particle )
    int counter = 0;
# endif 
# if !DLM_Moving_particle
    coord imposedU;
    coord imposedw;
# endif  


# if DLMFD_OPT
    // Fast loop structures to compute the DLMFD scalar products 
    BPFastLoop_LambdaMom LMloop;
    initialize_and_allocate_BPFastLoop_LambdaMom( &LMloop, DLMBLOCK );
    BPFastLoop_ResU RUloop;
    initialize_and_allocate_BPFastLoop_ResU( &RUloop, DLMBLOCK );
# endif  


  // Allocate and initialize particles
  allocate_and_init_particles( p, NPARTICLES, index_lambda, flagfield, 
  	flagfield_mailleur, DLM_periodic_shift );


  // Below we tranfer pointers in local arrays for ease of notation only  
# if debugInterior == 0
    Cache * Interior[NPARTICLES];
# endif

# if debugBD == 0
    Cache * Boundary[NPARTICLES];
    SolidBodyBoundary * sbm[NPARTICLES];
# endif
    GeomParameter * gci[NPARTICLES];

  
# if DLM_Moving_particle
#   if TRANSLATION
      coord * qU[NPARTICLES];
      coord * U[NPARTICLES];
      coord * tU[NPARTICLES];
#   endif
#   if ROTATION
      coord * qw[NPARTICLES];
      coord * w[NPARTICLES];
      coord * tw[NPARTICLES];
#   endif
# endif
  
  for (int k = 0; k < NPARTICLES; k++) 
  {
    gci[k] = &(p[k].g);
    
#   if debugInterior == 0
      Interior[k] = &(p[k].Interior);
#   endif
#   if debugBD == 0 
      Boundary[k] = &(p[k].reduced_domain);
      sbm[k] = &(p[k].s);
#   endif
    
#   if DLM_Moving_particle
#     if TRANSLATION
        qU[k] = &(p[k].qU);
        U[k] = &(p[k].U);
        tU[k] = &(p[k].tU);
#     endif
#     if ROTATION
        qw[k] = &(p[k].qw);
        w[k] = &(p[k].w);
        tw[k] = &(p[k].tw);
#     endif
#   endif
  }


  // In case 2 particles are too close, sort out their boundary points
  if ( NPARTICLES > 1 )
    remove_too_close_multipliers( p, index_lambda );


  // Tag the cells belonging to the stencil of the boundary-multipliers 
  // points, i.e. set flagfield to 1
# if debugBD == 0   
#   if DLM_Moving_particle
      reverse_fill_flagfield( p, flagfield_mailleur, index_lambda, 0, 
  	DLM_periodic_shift );
#   endif
  
    /* Consider fictitious domain's boundary points only if they are far
    enough of the domain's boundary (i.e a distance greater than
    2*Delta, with Delta being the smallest cell size) */
    double twodelta;

    foreach_level (depth()) 
    {
      twodelta = 2.*Delta;

#     if dimension == 2
        if (( x > (L0 - twodelta + X0)) || (x  < (twodelta + X0)) 
		|| ( y > (L0 - twodelta + Y0)) || (y < (twodelta + Y0)))
#     elif dimension == 3
        if ((x  > (L0 - twodelta + X0)) ||  (x  < (twodelta + X0)) 
      		|| (y > (L0 - twodelta + Y0)) || (y < (twodelta + Y0)) 
		|| (z > (L0 - twodelta +Z0)) || (z  < (twodelta + Z0)))
#     endif
	  if (index_lambda.x[] > -1.) 
	    index_lambda.x[] = -1.;
    }
    boundary((scalar*) {index_lambda});

    reverse_fill_flagfield( p, flagfield, index_lambda, 1, DLM_periodic_shift );
# endif


  // Create fictitious-domain's cache for the interior domain
# if debugInterior == 0 
    for (int k = 0; k < NPARTICLES; k++) 
    {
      switch( p[k].shape )
      {
        case SPHERE:
	  create_FD_Interior_Sphere( &p[k], index_lambda, DLM_periodic_shift, 
		flagfield );
	  break;
	  
	case CIRCULARCYLINDER2D:
	  create_FD_Interior_CircularCylinder2D( &p[k], index_lambda, 
	  	DLM_periodic_shift, flagfield );
	  break;
	  
	case CUBE:
	  create_FD_Interior_Cube( &p[k], index_lambda, DLM_periodic_shift );
	  break;
	  
	default:
          fprintf( stderr,"Unknown Rigid Body shape !!\n" );
      }
    }
# endif

# if PARTICLE_VERBOSE
    print_all_particles( p, &outputshift[0] );
# endif

# if DLMFD_OPT
    // Create the caches to loop over cells involved in DLMFD scalar products  
    // Cache for r, w, v and lambda
    Cache Traversal_rvwlambda;
    int bbb = 0;
    Point ppp;
    initialize_and_allocate_Cache( &Traversal_rvwlambda );
    // Interior points
    for (int m = 0; m < NPARTICLES; m++)
    {     
      bbb = 0; 
      foreach_cache((*Interior[m])) 
      {      
        if ( flagfield[] < 1 )
        {
          ppp.i = (Interior[m]->p)[bbb].i;
          ppp.j = (Interior[m]->p)[bbb].j;
          ppp.k = (Interior[m]->p)[bbb].k;		
          ppp.level = (Interior[m]->p)[bbb].level;	
	  cache_append( &Traversal_rvwlambda, ppp, 0);
        }
        bbb++;
      }
    }
    // Boundary points
    for (int m = 0; m < NPARTICLES; m++) 
    { 
      bbb = 0;
      foreach_cache((*Boundary[m])) 
      {   
        if ( index_lambda.x[] > -1 ) 
        {
          ppp.i = (Boundary[m]->p)[bbb].i;
          ppp.j = (Boundary[m]->p)[bbb].j;
          ppp.k = (Boundary[m]->p)[bbb].k;		
          ppp.level = (Boundary[m]->p)[bbb].level;	
	  cache_append( &Traversal_rvwlambda, ppp, 0);
        }
        bbb++;
      }
    } 
  
    // Cache for u, qu and tu
    Cache Traversal_uqutu;
    initialize_and_allocate_Cache( &Traversal_uqutu );
    // Interior points
    for (int m = 0; m < NPARTICLES; m++)
    {     
      bbb = 0; 
      foreach_cache((*Interior[m])) 
      {      
        if ( flagfield[] < 1 )
        {
          ppp.i = (Interior[m]->p)[bbb].i;
          ppp.j = (Interior[m]->p)[bbb].j;
          ppp.k = (Interior[m]->p)[bbb].k;		
          ppp.level = (Interior[m]->p)[bbb].level;	
	  cache_append( &Traversal_uqutu, ppp, 0);
        }
        bbb++;
      }
    }
    // Boundary points
    foreach_level(depth()) 
      if ( flagfield[] > 0.5 )
        cache_append( &Traversal_uqutu, point, 0);  
# endif  
      
  boundary ((scalar*) {DLM_periodic_shift});


# if DLMFD_OPT    
    // Only DLM_lambda has not been nullified at the end of the previous
    // call to DLMFD_subproblem as it is needed for the computation of the 
    // hydrodynamic force & torque for post-processing, so we nullify it here
    foreach() 
    {
      tcells ++; 
      foreach_dimension() 
        DLM_lambda.x[] = 0.;      
    }
# else
    foreach() 
    {
      tcells ++; 
      foreach_dimension()
      {
        DLM_lambda.x[] = 0.;
	DLM_r.x[] = 0.;
	DLM_w.x[] = 0.; 
	DLM_v.x[] = 0.; 
	qu.x[] = 0.; 
	tu.x[] = 0.;
      }
    }  
# endif  
  
# if _MPI
    mpi_all_reduce( tcells, MPI_INTEGER, MPI_SUM );
# endif


  // Statistics of number of Lagrange multiplier points  
  lm = total_dlmfd_multipliers (p, NPARTICLES);
  /* Get track of the number of multipliers for statistics */
  for (int k = 0; k < NPARTICLES; k++)
    p[k].tmultipliers += lm;

  allpts = total_dlmfd_cells (p, NPARTICLES);
  /* Get track of the number of cells involved in the dlmfd solver for 
  statistics */
  for (int k = 0; k < NPARTICLES; k++)
    p[k].tcells += allpts;

  if ( pid() == 0 )
  {
    printf( "      DLM points = %d, constrained cells = %d\n", lm, allpts );
    fprintf( cellvstime, "%d \t %d \t \t %d \t \t \t %d \n", i, lm, allpts, 
  	tcells );
    fflush( cellvstime );
  } 
  

  // Nullify the qU, tU, qw and tw vectors of all particles
# if DLM_Moving_particle
#   if TRANSLATION
      for (int k = 0; k < NPARTICLES; k++) 
        foreach_dimension()
        {
	  (*qU[k]).x = 0.; 
	  (*tU[k]).x = 0.;
        }
#   endif
#   if ROTATION
      for (int k = 0; k < NPARTICLES; k++) 
        foreach_dimension()
        {
	  (*qw[k]).x = 0.; 
	  (*tw[k]).x = 0.;
        }
#   endif
# endif



  /** # Iterative Uzawa/conjugate-gradient solver 
  Implementation of the iterative Uzawa/conjugate gradient solver for
  the fictitious domain problem. The code below follows closely the
  steps given in the annex of the paper "Wachs et al J.Eng. Math,
  2011". (minus all the sign errors)*/
 
  /* Initialization */
  /* ------------- */
  /* Compute qu = fu - (M_u^T)*lambda^0 with fu = rho*dV/dt*u^n+1/2 +
     alpha*lambda^n. */ 
  /* Note that (M_u^T)*lambda <=> <lambda,v>_P(t) (v is the test function
     for u)*/
  
  /* For moving particles add:  */
  /* qU = fU -(M_U^T)*lambda^0 with fU = (1 - rho_f /
     rho_s)*M*(U^{n+1/2}/dt + g) with M the particle's mass and g the
     gravity. */
  /* Note that (M_U^T)*lambda <=> - <lambda, V>_P(t) (V is the test
     function for U)*/
  /* qw = fw -(M_w^T)*lambda^0 with fw=(1-rho_f/rho_s)/dt*Ip*w^n+1/2
   * with Ip the Inertia tensor.*/
  /* Note that (M_w^T)*lambda <=> - <lambda, xi^r_GM>_P(t) (xi is the test
   * function for w)*/

  
  foreach() 
  {
    foreach_dimension() 
    {
      qu.x[] = rho_f*dlmfd_dv()*u.x[]/dt;

#     if DLM_alpha_coupling
        qu.x[] += DLM_explicit.x[];
        DLM_explicit.x[] = 0.;
#     endif
    }
  }
  
  /* Interior points qu = fu -(M_u^T)*lambda^0  */
  /* qu = fu -(M_u^T)*lambda^0 = fu -<lambda,v>_P(t) */

  /* Interior points qU = fU -(M_U^T)*lambda^0  */
  /* qU = fU -(M_U^T)*lambda^0 = fU  + <lambda,V>_P(t) */
  
  /* Interior points qw = fw -(M_w^T)*lambda^0  */
  /* qw = fw -(M_w^T)*lambda^0 = fw  + <lambda,xi^r_GM>_P(t) */
  
  /* For moving particles the fU and fw parts are added after the
     scalar product (for mpi purpose) */

# if debugInterior == 0
    for (int k = 0; k < NPARTICLES; k++) 
    {
      foreach_cache((*Interior[k])) 
      {
        if ( flagfield[] < 1 && (int)index_lambda.y[] == k ) 
        {
	  foreach_dimension() 
	    qu.x[] -= DLM_lambda.x[];
	
#         if DLM_Moving_particle
#           if TRANSLATION
              foreach_dimension()
	        (*qU[k]).x += DLM_lambda.x[];
#           endif
#           if ROTATION
	      /* Modify temporarily the particle center position for periodic 
	      boundary condition */
              foreach_dimension()
	        (*gci[k]).center.x += DLM_periodic_shift.x[];
	  
	      /* as a.(b^c) = c.(a^b) = b.(c^a) */
	      /* <lambda,xi^r_GM>_P(t)=<r_GM,lambda^xi>_P(t)
	      	=<xi,r_GM^lambda>_P(t)*/
#             if dimension == 3  
	        // lambda_z*r_y - lambda_y*r_z 
	        (*qw[k]).x += DLM_lambda.z[] * ( y - (*gci[k]).center.y ) 
			- DLM_lambda.y[] * ( z - (*gci[k]).center.z ); 
	        // lambda_x*r_z - lambda_z*r_x 
	        (*qw[k]).y += DLM_lambda.x[] * ( z - (*gci[k]).center.z ) 
			- DLM_lambda.z[] * ( x - (*gci[k]).center.x );
#             endif 
	      // lambda_y*r_x - lambda_x*r_y
	      (*qw[k]).z += DLM_lambda.y[] * ( x - (*gci[k]).center.x ) 
		- DLM_lambda.x[] * ( y - (*gci[k]).center.y ); 

	      foreach_dimension()
	        (*gci[k]).center.x -= DLM_periodic_shift.x[];
#           endif
#         endif
        }
      }
    }
# endif

  /* Boundary points qu = fu -(M_u^T)*lambda^0 */
  /* Boundary points qU = fU -(M_U^T)*lambda^0 */
  /* Boundary points qw = fw -(M_w^T)*lambda^0 */
  
# if debugBD == 0
    double weight = 0.;
    coord weightcellpos = {0., 0., 0.};
    coord lambdacellpos = {0., 0., 0.};
    coord lambdapos = {0., 0., 0.};
    coord sum = {0., 0., 0.};

    boundary ((scalar*) {DLM_lambda});

#   if DLMFD_OPT
      double* qux_ = NULL; 
      double* quy_ = NULL;
      double* quz_ = NULL;
#   endif
 
    for (int k = 0; k < NPARTICLES; k++) 
    {
      particle * pp = &p[k];
    
      foreach_cache ((*Boundary[k])) 
      {
#       if DLMFD_OPT
          qux_ = &(qu.x[]); 
	  quy_ = &(qu.y[]); 
#         if dimension == 3	  
	    quz_ = &(qu.z[]);
#         endif	  
#       endif  
	  
        weightcellpos.x = x; 
	weightcellpos.y = y; 
#       if dimension == 3
          weightcellpos.z = z;
#       endif

        foreach_dimension() 
	  sum.x = 0.; 

        foreach_neighbor() 
        {
	  if ( (int)index_lambda.x[] > -1 && level == depth() 
		&& is_leaf(cell) && (int)index_lambda.y[] == pp->pnum ) 
	  {
	    lambdacellpos.x = x; 
	    lambdacellpos.y = y; 
	    lambdapos.x = (*sbm[k]).x[(int)index_lambda.x[]];
	    lambdapos.y = (*sbm[k]).y[(int)index_lambda.x[]];	  
#           if dimension == 3
	      lambdacellpos.z = z;
	      lambdapos.z = (*sbm[k]).z[(int)index_lambda.x[]];
#           endif
	  
	    foreach_dimension()
	      ppshift.x = DLM_periodic_shift.x[];
	    
	    weight = reversed_weight( pp, weightcellpos, lambdacellpos, 
	  	lambdapos, Delta, ppshift );
	  
	    foreach_dimension()
	      sum.x += weight*DLM_lambda.x[];

#           if DLMFD_OPT	    
	      // Append to the fast loop for subsequent 
	      // M_u^T*w = <DLM_w, v>_P(t) computations over the iterative 
	      // process
	      if ( fabs(weight) > 1.e-8 )
#               if dimension == 3
	          append_BPFastLoop_LambdaMom( &LMloop, qux_, quy_, quz_, 
	      		&(DLM_w.x[]), &(DLM_w.y[]), &(DLM_w.z[]), weight );
#               else	      
	          append_BPFastLoop_LambdaMom_2D( &LMloop, qux_, quy_,
	      		&(DLM_w.x[]), &(DLM_w.y[]), weight );
#               endif			
#           endif		
	  }
        }

        // -= here as one fluid cell can be affected by multiples particle's
        // boundary multipliers 
	foreach_dimension() 
	  qu.x[] -= sum.x; 
      }
    }
    
#   if DLM_Moving_particle
      for (int k = 0; k < NPARTICLES; k++) 
      {
        particle * pp = &p[k];
    
        foreach_cache ((*Boundary[k])) 
        {
          if ( index_lambda.x[] > -1 && (int)index_lambda.y[] == pp->pnum ) 
          {
	    lambdapos.x = (*sbm[k]).x[(int)index_lambda.x[]];
	    lambdapos.y = (*sbm[k]).y[(int)index_lambda.x[]];
#           if dimension == 3
	      lambdapos.z = (*sbm[k]).z[(int)index_lambda.x[]];
#           endif

#           if TRANSLATION
	      foreach_dimension()
	        (*qU[k]).x += DLM_lambda.x[];
#           endif
#           if ROTATION
	      /* Modify temporarily the particle center position for periodic 
	      boundary condition */
	      foreach_dimension()
	        (*gci[k]).center.x += DLM_periodic_shift.x[];
	
	      /* as a.(b^c) = c.(a^b) = b.(c^a) */
	      // <lambda,xi^r_GM>_P(t)=<r_GM,lambda^xi>_P(t)
	      //	=<xi,r_GM^lambda>_P(t)
#             if dimension == 3
	        // lambda_z*r_y - lambda_y*r_z
	        (*qw[k]).x += 
			DLM_lambda.z[] * ( lambdapos.y - (*gci[k]).center.y ) 
			- DLM_lambda.y[] * ( lambdapos.z - (*gci[k]).center.z ); 
	        // lambda_x*r_z - lambda_z*r_x
	        (*qw[k]).y += 
			DLM_lambda.x[] * ( lambdapos.z - (*gci[k]).center.z ) 
			- DLM_lambda.z[] * ( lambdapos.x - (*gci[k]).center.x ); 
#             endif
	      // lambda_y*r_x - lambda_x*r_y
	      (*qw[k]).z += 
	      	DLM_lambda.y[] * ( lambdapos.x - (*gci[k]).center.x ) 
		- DLM_lambda.x[] * ( lambdapos.y - (*gci[k]).center.y );  

	      foreach_dimension()
	        (*gci[k]).center.x -= DLM_periodic_shift.x[];
#           endif
          }      
        }
      }
#   endif
# endif
  

# if DLM_Moving_particle
#   if _MPI /* _MPI Reduction */
      // Reduce to master
      // Pack data
      counter = 0;  
      for (int k = 0; k < NPARTICLES; k++) 
      {
#       if TRANSLATION
          vpartbuf[counter] = (*qU[k]).x;
          vpartbuf[counter+1] = (*qU[k]).y;
#         if dimension == 3	      
            vpartbuf[counter+2] = (*qU[k]).z;
            counter += 3;
#         else
            counter += 2;
#         endif    
#       endif
#       if ROTATION
#         if dimension == 3
            vpartbuf[counter] = (*qw[k]).x;
            vpartbuf[counter+1] = (*qw[k]).y;    
            vpartbuf[counter+2] = (*qw[k]).z;
            counter += 3; 
#         else
            vpartbuf[counter] = (*qw[k]).z;
	    counter += 1; 
#         endif 	    
#       endif
      }
  
      // Perform reduction on Master
      MPI_Reduce( pid() ? vpartbuf : MPI_IN_PLACE, vpartbuf, 
      	npartdata*NPARTICLES, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

      if ( pid() == 0 )
      {
        // Unpack data    
        counter = 0;
        for (int k = 0; k < NPARTICLES; k++) 
        {
#         if TRANSLATION
            (*qU[k]).x = vpartbuf[counter];
            (*qU[k]).y = vpartbuf[counter+1];
#           if dimension == 3	          
              (*qU[k]).z = vpartbuf[counter+2];      
              counter += 3;
#           else
              counter += 2;
#           endif	          
#         endif
#         if ROTATION
#           if dimension == 3
              (*qw[k]).x = vpartbuf[counter];
              (*qw[k]).y = vpartbuf[counter+1];    
              (*qw[k]).z = vpartbuf[counter+2];    
              counter += 3; 
#           else
              (*qw[k]).z = vpartbuf[counter];    
              counter += 1;
#           endif	      
#         endif
        }	   	
#   endif  /* end of _MPI Reduction */    
    
        // Perform the inversion (on master when in MPI)
        for (int k = 0; k < NPARTICLES; k++) 
        {
#         if TRANSLATION
            /* Add here fU to qU */
	    foreach_dimension()
              (*qU[k]).x += ( 1. - ( rho_f / p[k].rho_s ) ) *
    		p[k].M * ( p[k].gravity.x + p[k].addforce.x ) 
		+ p[k].DLMFD_couplingfactor * p[k].M * (*U[k]).x / dt ;
	
            /* Solution of M * DLMFD_couplingfactor * U / dt = qU */
            /* U = ( dt * qU ) / ( DLMFD_couplingfactor * M ) */
            foreach_dimension() 
	      (*U[k]).x = ( dt / ( p[k].DLMFD_couplingfactor * p[k].M ) ) 
	      	* (*qU[k]).x ;
#         endif
#         if ROTATION
            /* Add here fw to qw */
            /* The inertia tensor is */
            /*  Ixx  Ixy  Ixz */
            /*  Iyx  Iyy  Iyz */
            /*  Izx  Izy  Izz */ 
            /* with */
            /* Ip[0] = Ixx */
            /* Ip[1] = Iyy */
            /* Ip[2] = Izz */
            /* Ip[3] = Ixy */
            /* Ip[4] = Ixz */
            /* Ip[5] = Iyz */
#           if dimension == 3
              (*qw[k]).x += p[k].DLMFD_couplingfactor * 
	      		( (p[k].Ip[0]) * (*w[k]).x
    			- (p[k].Ip[3]) * (*w[k]).y 
			- (p[k].Ip[4]) * (*w[k]).z ) / dt;
              (*qw[k]).y += p[k].DLMFD_couplingfactor * 
	      		( - (p[k].Ip[3]) * (*w[k]).x
    			+ (p[k].Ip[1]) * (*w[k]).y 
			- (p[k].Ip[5]) * (*w[k]).z ) / dt;
#           endif
            (*qw[k]).z += p[k].DLMFD_couplingfactor * 
	    	( - (p[k].Ip[4]) * (*w[k]).x 
    		- (p[k].Ip[5]) * (*w[k]).y 
		+ (p[k].Ip[2]) * (*w[k]).z ) / dt;   

            /* Solution of Ip * DLMFD_couplingfactor * w / dt = qw */
            /* w = ( dt / ( DLMFD_couplingfactor ) * Ip_inv * qw 
            /* where Ip_inv is the inverse of Ip */        
#           if dimension == 3
              (*w[k]).x = ( dt / p[k].DLMFD_couplingfactor ) * 
    		( (p[k].Ip_inv)[0][0] * (*qw[k]).x 
		+ (p[k].Ip_inv)[0][1] * (*qw[k]).y 
    		+ (p[k].Ip_inv)[0][2] * (*qw[k]).z );
              (*w[k]).y = ( dt / p[k].DLMFD_couplingfactor ) * 
    		( (p[k].Ip_inv)[1][0] * (*qw[k]).x 
		+ (p[k].Ip_inv)[1][1] * (*qw[k]).y 
    		+ (p[k].Ip_inv)[1][2] * (*qw[k]).z );
              (*w[k]).z = ( dt / p[k].DLMFD_couplingfactor ) * 
    		( (p[k].Ip_inv)[2][0] * (*qw[k]).x 
		+ (p[k].Ip_inv)[2][1] * (*qw[k]).y 
    		+ (p[k].Ip_inv)[2][2] * (*qw[k]).z );
#           else
              (*w[k]).z = ( dt / p[k].DLMFD_couplingfactor ) * 
		(p[k].Ip_inv)[2][2] * (*qw[k]).z ;		
#           endif		
#         endif
        }  

#   if _MPI /* _MPI Broadcast */
        // Broadcast U and w
        // Pack data
        counter = 0;  
        for (int k = 0; k < NPARTICLES; k++) 
        {
#         if TRANSLATION
            vpartbuf[counter] = (*U[k]).x;
            vpartbuf[counter+1] = (*U[k]).y;    
#           if dimension == 3
              vpartbuf[counter+2] = (*U[k]).z;
              counter += 3;
#           else
              counter += 2;
#           endif	          
#         endif
#         if ROTATION
#           if dimension == 3
              vpartbuf[counter] = (*w[k]).x;
              vpartbuf[counter+1] = (*w[k]).y;    
              vpartbuf[counter+2] = (*w[k]).z;
              counter += 3;
#           else
              vpartbuf[counter] = (*w[k]).z;
	      counter += 1;
#           endif	       
#         endif
        }
      } /* End of "if pid() == 0" */    
  
      // Perform the broadcast from the master to the other processes
      MPI_Bcast( vpartbuf, npartdata*NPARTICLES, 
    	MPI_DOUBLE, 0, MPI_COMM_WORLD );

      // Unpack data    
      counter = 0;
      for (int k = 0; k < NPARTICLES; k++) 
      {
#       if TRANSLATION
          (*U[k]).x = vpartbuf[counter];
          (*U[k]).y = vpartbuf[counter+1];      
#         if dimension == 3
            (*U[k]).z = vpartbuf[counter+2];      
            counter += 3;
#         else	
            counter += 2;
#         endif        
#       endif
#       if ROTATION
#         if dimension == 3
            (*w[k]).x = vpartbuf[counter];
            (*w[k]).y = vpartbuf[counter+1];    
            (*w[k]).z = vpartbuf[counter+2];    
            counter += 3; 
#         else
            (*w[k]).z = vpartbuf[counter];    
            counter += 1; 	    
#         endif
#       endif
      }  
#   endif /* end of _MPI Broadcast */ 
# endif /* end of DLM_Moving_particle */

  
  /* Invert L*u^0 = qu with L=dv*rho/dt*Identity_Matrix */
  foreach()
    foreach_dimension()
      u.x[] = qu.x[] * dt / ( rho_f * dlmfd_dv() );


  /* Compute residual r^0 = G - M_u*u^0 - M_U*U^0 - M_w*w^0 */
  /* Note that M_u*u^0 <=> <alpha, u>_P(t) */
  /* Note that M_U*U^0 <=> - <alpha, U>_P(t) */
  /* Note that M_w*w^0 <=> - <alpha, w^r_GM>_P(t) */
  /* alpha is he test function for lambda the Lagrange multipliers */

  /* So r^0 = G -<alpha, u>_P(t) + <alpha, U>_P(t) + <alpha, w^r_GM>_P(t) */

  /* In the 2D case when rotation is enabled, the z-component of the
   * residual becomes involved in the algorithm, this is why the
   * foreach_dimension() iterators are being removed from this
   * point. Otherwise the z components would be discarded */


  /* Interior points: r^0 = G - M_u*u^0 - M_U*U^0 - M_w*w^0 */
  /* So r^0 = G -<alpha, u>_P(t) + <alpha, U>_P(t) + <alpha, w^r_GM>_P(t) */
# if debugInterior == 0
    for (int k = 0; k < NPARTICLES; k++) 
    {        
#     if !DLM_Moving_particle
        imposedU = (p[k]).imposedU;
        imposedw = (p[k]).imposedw;
#     endif

      foreach_cache((*Interior[k])) 
      {      
        if ((flagfield[]  < 1) && ((int)index_lambda.y[] == k)) 
        {
	  foreach_dimension() 
	    DLM_r.x[] = - u.x[];
 
#         if DLM_Moving_particle
#           if TRANSLATION
              foreach_dimension() 
	        DLM_r.x[] += (*U[k]).x;			
#           endif
#           if ROTATION
	      /* Modify temporarily the particle center position for periodic 
	      boundary condition */
	      foreach_dimension()
	        (*gci[k]).center.x += DLM_periodic_shift.x[];
	
#             if dimension == 3
	        // w_y*r_z - w_z*r_y
	        DLM_r.x[] += (*w[k]).y * ( z - (*gci[k]).center.z ) 
			- (*w[k]).z * ( y - (*gci[k]).center.y ); 
	        // w_z*r_x - w_x*r_z
	        DLM_r.y[] += (*w[k]).z * ( x - (*gci[k]).center.x ) 
			- (*w[k]).x * ( z - (*gci[k]).center.z ); 
#             endif
	      // w_x*r_y - w_y*r_x
	      DLM_r.z[] += (*w[k]).x * ( y - (*gci[k]).center.y ) 
			- (*w[k]).y * ( x - (*gci[k]).center.x ); 

	      foreach_dimension()
	        (*gci[k]).center.x -= DLM_periodic_shift.x[];
#           endif
#         else
	    foreach_dimension()   
	      DLM_r.x[] += imposedU.x;
	      
	    /* Modify temporarily the particle center position for periodic 
	    boundary condition */
	    foreach_dimension()
	      (*gci[k]).center.x += DLM_periodic_shift.x[];

#           if dimension == 3
	      DLM_r.x[] += imposedw.y * ( z - (*gci[k]).center.z ) 
		- imposedw.z * ( y - (*gci[k]).center.y);
	      DLM_r.y[] += imposedw.z * ( x - (*gci[k]).center.x ) 
		- imposedw.x * ( z - (*gci[k]).center.z);
#           endif		
	    DLM_r.z[] += imposedw.x * ( y - (*gci[k]).center.y ) 
		- imposedw.y * ( x - (*gci[k]).center.x);

	    foreach_dimension()
	      (*gci[k]).center.x -= DLM_periodic_shift.x[];	      
#         endif	
        }
      }
    }
# endif

  
  /* Boundary points: r^0 = G - M_u*u^0 - M_U*U^0 - M_w*w^0 */
  /* So r^0 = G -<alpha, u>_P(t) + <alpha, U>_P(t) + <alpha, w^r_GM>_P(t) */
# if debugBD == 0
    double testweight = 0.;
    boundary((scalar*){u});
#   if DLMFD_OPT
      int ndof = 0;
#   endif  

    for (int k = 0; k < NPARTICLES; k++) 
    {
      particle* pp = &p[k];
#     if !DLM_Moving_particle
        imposedU = (p[k]).imposedU;
        imposedw = (p[k]).imposedw;
#     endif      
      
      foreach_cache((*Boundary[k])) 
      {
#       if DLMFD_OPT
          ndof = 0;
#       endif      
        if (index_lambda.x[] > -1 && ((int)index_lambda.y[] == pp->pnum)) 
        {
	  lambdacellpos.x = x; 
	  lambdacellpos.y = y; 
	  lambdapos.x = (*sbm[k]).x[(int)index_lambda.x[]];
	  lambdapos.y = (*sbm[k]).y[(int)index_lambda.x[]];
#         if dimension == 3
	    lambdacellpos.z = z;
	    lambdapos.z = (*sbm[k]).z[(int)index_lambda.x[]];
#         endif

	  foreach_dimension() 
	    sum.x = 0.;
	    
	  testweight = 0.;
	
	  foreach_dimension()
	    ppshift.x = DLM_periodic_shift.x[];

	  foreach_neighbor() 
	  {
	    if ( level == depth() ) 
	    {
	      weightcellpos.x = x; 
	      weightcellpos.y = y; 
#             if dimension == 3
	        weightcellpos.z = z;
#             endif
	
	      weight = reversed_weight( pp, weightcellpos, lambdacellpos, 
	    	lambdapos, Delta, ppshift );

	      testweight += weight;
	      
	      foreach_dimension()
	        sum.x += weight * u.x[];

#             if DLMFD_OPT	      
	        // Append to fast loop for subsequent M_u*tu = <alpha, tu>_P(t)
	        // computations over the iterative process
	        if ( fabs(weight) > 1.e-8 )
	        {
#                 if dimension == 3
	            append_BPFastLoop_ResU_tu( &RUloop,
			&(tu.x[]), &(tu.y[]),  &(tu.z[]), weight );
#                 else
	            append_BPFastLoop_ResU_tu_2D( &RUloop,
			&(tu.x[]), &(tu.y[]),  weight );
#                 endif			
	          ++ndof;
	        }
#             endif	    	  
	    }
	  }

	  if ( testweight < 1. - 0.000001 || testweight > 1. + 0.000001 )
	    printf( "testweight = %f\n", testweight );
	
	  foreach_dimension() 
	    DLM_r.x[] = - sum.x; 

#         if DLMFD_OPT	
	    // Append to the fast loop for subsequent M_u*tu = <alpha, tu>_P(t)
	    // computations over the iterative process
#           if dimension == 3
	      append_BPFastLoop_ResU_v( &RUloop, &(DLM_v.x[]), &(DLM_v.y[]),
		&(DLM_v.z[]), ndof );
#           else
	      append_BPFastLoop_ResU_v_2D( &RUloop, &(DLM_v.x[]), &(DLM_v.y[]),
		ndof );
#           endif			
#         endif
			
#         if DLM_Moving_particle
#           if TRANSLATION
	      foreach_dimension()
	        DLM_r.x[] += (*U[k]).x;
#           endif
#           if ROTATION
	      /* Modify temporarily the particle center position for periodic 
	      boundary condition */
	      foreach_dimension()
	        (*gci[k]).center.x += DLM_periodic_shift.x[];

#             if dimension == 3		
	        // w_y*r_z - w_z*r_y
	        DLM_r.x[] +=  (*w[k]).y * ( lambdapos.z - (*gci[k]).center.z ) 
			- (*w[k]).z * ( lambdapos.y - (*gci[k]).center.y); 
	        // w_z*r_x - w_x*r_z
	        DLM_r.y[] +=  (*w[k]).z * ( lambdapos.x - (*gci[k]).center.x ) 
			- (*w[k]).x * ( lambdapos.z - (*gci[k]).center.z); 
#             endif
	      // w_x*r_y - w_y*r_x
	      DLM_r.z[] +=  (*w[k]).x * ( lambdapos.y - (*gci[k]).center.y ) 
		- (*w[k]).y * ( lambdapos.x - (*gci[k]).center.x ); 

	      foreach_dimension()
	        (*gci[k]).center.x -= DLM_periodic_shift.x[];	
#           endif
#         else
	    foreach_dimension()
	      DLM_r.x[] += imposedU.x;

	    /* Modify temporarily the particle center position for periodic 
	    boundary condition */
	    foreach_dimension()
	      (*gci[k]).center.x += DLM_periodic_shift.x[];
	
#           if dimension == 3	
	      DLM_r.x[] += 
	    	imposedw.y * ( lambdapos.z - (*gci[k]).center.z ) 
		- imposedw.z * ( lambdapos.y - (*gci[k]).center.y );
	      DLM_r.y[] +=  
	    	imposedw.z * ( lambdapos.x - (*gci[k]).center.x ) 
		- imposedw.x * ( lambdapos.z - (*gci[k]).center.z );
#           endif
	    DLM_r.z[] +=  
	    	imposedw.x * ( lambdapos.y - (*gci[k]).center.y ) 
		- imposedw.y * ( lambdapos.x - (*gci[k]).center.x );

	    foreach_dimension()
	      (*gci[k]).center.x -= DLM_periodic_shift.x[];	
#         endif
        }
      }
    }
# endif
  
  DLM_nr2 = 0.;
  
  /* set DLM_w = DLM_r */
# if DLMFD_OPT
    foreach_cache( Traversal_rvwlambda )
# else
    foreach() 
# endif
  {
    foreach_dimension() DLM_w.x[] = DLM_r.x[];  
#   if dimension == 3
      DLM_nr2 += sq(DLM_r.x[]) + sq(DLM_r.y[]) + sq(DLM_r.z[]);
# else
      DLM_nr2 += sq(DLM_r.x[]) + sq(DLM_r.y[]);
# endif      
  }

# if _MPI
    mpi_all_reduce( DLM_nr2, MPI_DOUBLE, MPI_SUM );
# endif

  
  /* Iterative loop */
  /* -------------- */  
  for (ki = 1; ki < DLM_maxiter && (sqrt(DLM_nr2) > DLM_tol || ki < 3); ki++) 
  {

    /* (1) Initialize and compute qu = (M^T)*DLM_w */
#   if DLMFD_OPT
      // At the first iteration, we need to nullify over all cells as at the 
      // initialization step we add the term fu = rho*dV/dt*u^n+1/2 to qu in
      // all cells, but from the 2nd iteration we can nullify over cells that are
      // involved in Lagrange multiplier stencils only
      if ( ki == 1 )
        foreach() 
	  foreach_dimension() 
	    qu.x[] = 0.;
      else    
        foreach_cache( Traversal_uqutu )
          foreach_dimension() 
	    qu.x[] = 0.; 
#   else
      foreach() 
        foreach_dimension()
        {
          qu.x[] = 0.; 
          tu.x[] = 0.; 
        }    
#   endif   
    
#   if DLM_Moving_particle
      for (int k = 0; k < NPARTICLES; k++) 
      {
#       if TRANSLATION
          foreach_dimension()
          {
            (*qU[k]).x = 0.;
	    (*tU[k]).x = 0.;
          }
#       endif
#       if ROTATION
#         if dimension == 3
            (*qw[k]).x = 0.; (*tw[k]).x = 0.;
            (*qw[k]).y = 0.; (*tw[k]).y = 0.;
#         endif
          (*qw[k]).z = 0.; (*tw[k]).z = 0.;
#       endif
      }
#   endif
     
    /* Interior points qu = (M_u^T)DLM_w =   <DLM_w, v>_P(t) */
    /* Interior points qU = (M_U^T)DLM_w = - <DLM_w, V>_P(t) */
    /* Interior points qw = (M_w^T)DLM_w = - <DLM_w, xi^r_GM>_P(t) */
    /* -<DLM_w, xi^r_GM>_P(t)=-<r_GM, DLM_w^xi>_P(t)=-<xi, r_GM^DLM_w>_P(t) */    
#   if debugInterior == 0
      for (int k = 0; k < NPARTICLES; k++) 
      {
        foreach_cache((*Interior[k])) 
        {
          if ((flagfield[]  < 1) && ((int)index_lambda.y[] == k)) 
          {
	    foreach_dimension() 
	      qu.x[] = DLM_w.x[];
	  
#           if DLM_Moving_particle
#             if TRANSLATION
                foreach_dimension() 
		  (*qU[k]).x += -DLM_w.x[];
#             endif
#             if ROTATION
	        /* Modify temporarily the particle center position for 
		periodic boundary condition */
	        foreach_dimension()
	          (*gci[k]).center.x += DLM_periodic_shift.x[];

#               if dimension == 3	    
	          // -(w_z*r_y - w_y*r_z)
	          (*qw[k]).x -=  DLM_w.z[] * ( y - (*gci[k]).center.y ) 
			- DLM_w.y[] * ( z - (*gci[k]).center.z ); 
	          // -(w_x*r_z - w_z*r_x)
	          (*qw[k]).y -=  DLM_w.x[] * ( z - (*gci[k]).center.z ) 
			- DLM_w.z[] * ( x - (*gci[k]).center.x ); 
#               endif 
	          // -(w_y*r_x - w_x*r_y)
	          (*qw[k]).z -=  DLM_w.y[] * ( x - (*gci[k]).center.x ) 
			- DLM_w.x[] * ( y - (*gci[k]).center.y );

	        foreach_dimension()
	          (*gci[k]).center.x -= DLM_periodic_shift.x[];	
#             endif
#           endif
          }
        }
      }
#   endif

  
    /* Boundary points qu = (M_u^T)DLM_w =   <DLM_w, v>_P(t) */
    /* Boundary points qU = (M_U^T)DLM_w = - <DLM_w, V>_P(t) */
    /* Boundary points qw = (M_w^T)DLM_w = - <DLM_w, xi^r_GM>_P(t) */
    /* -<DLM_w, xi^r_GM>_P(t)=-<r_GM, DLM_w^xi>_P(t)=-<xi, r_GM^DLM_w>_P(t) */  
#   if debugBD == 0
#     if DLMFD_OPT
        // Use of the fast loop for the computations of 
        // qu = (M_u^T)*DLM_w = <DLM_w, v>_P(t) over the boundary points
        boundary ((scalar*){DLM_w});
        int endloop = LMloop.n;  
        for (int k=0; k<endloop; ++k )
        {
          weight = LMloop.weight[k];
          *(LMloop.qux[k]) += weight * *(LMloop.dlmwx[k]);
          *(LMloop.quy[k]) += weight * *(LMloop.dlmwy[k]);    
#         if dimension == 3	
            *(LMloop.quz[k]) += weight * *(LMloop.dlmwz[k]); 
#         endif	       
        }
#     else
        boundary ((scalar*){DLM_w, qu});    
        for (int k = 0; k < NPARTICLES; k++) 
        {
          particle * pp = &p[k];
          foreach_cache((*Boundary[k])) 
          {
            weightcellpos.x = x; 
	    weightcellpos.y = y; 
#           if dimension == 3
              weightcellpos.z = z;
#           endif

            foreach_dimension() 
	      sum.x = 0; 

            foreach_neighbor() 
	    {
	      if ( (int)index_lambda.x[] > -1 && level == depth() 
		&& is_leaf(cell) && (int)index_lambda.y[] == pp->pnum ) 
	      {
	        lambdacellpos.x = x;
	        lambdacellpos.y = y;
	        lambdapos.x = (*sbm[k]).x[(int)index_lambda.x[]];
	        lambdapos.y = (*sbm[k]).y[(int)index_lambda.x[]];
#               if dimension == 3
	          lambdacellpos.z = z;
	          lambdapos.z = (*sbm[k]).z[(int)index_lambda.x[]];
#               endif

	        foreach_dimension()
	          ppshift.x = DLM_periodic_shift.x[];
	  
	        weight = reversed_weight( pp, weightcellpos, lambdacellpos, 
			lambdapos, Delta, ppshift );

	        foreach_dimension()
	          sum.x += weight*DLM_w.x[]; 
	      }
            }
      
            // += here as one fluid cell can be affected by multiples 
	    // particle's boundary multipliers
            foreach_dimension() 
	      qu.x[] += sum.x;
          }
        } 
#     endif

	
#     if DLM_Moving_particle
        for (int k = 0; k < NPARTICLES; k++) 
        {
          particle * pp = &p[k];
          foreach_cache((*Boundary[k])) 
          {
            if ( index_lambda.x[] > -1 && (int)index_lambda.y[] == pp->pnum ) 
            {
	      lambdapos.x = (*sbm[k]).x[(int)index_lambda.x[]];
	      lambdapos.y = (*sbm[k]).y[(int)index_lambda.x[]];	
#             if dimension == 3
	        lambdapos.z = (*sbm[k]).z[(int)index_lambda.x[]];
#             endif

#             if TRANSLATION
	        foreach_dimension() 
		  (*qU[k]).x += -DLM_w.x[];
#             endif
#             if ROTATION
	        /* Modify temporarily the particle center position for periodic 
	        boundary condition */
	        foreach_dimension()
	          (*gci[k]).center.x += DLM_periodic_shift.x[];

#               if dimension == 3	  
	          // -(w_z*r_y - w_y*r_z)
	          (*qw[k]).x -= DLM_w.z[] * ( lambdapos.y - (*gci[k]).center.y )
			- DLM_w.y[] * ( lambdapos.z - (*gci[k]).center.z ); 
	          // -(w_x*r_z - w_z*r_x)
	          (*qw[k]).y -= DLM_w.x[] * ( lambdapos.z - (*gci[k]).center.z )
			- DLM_w.z[] * ( lambdapos.x - (*gci[k]).center.x ); 
#               endif	        
		// -(w_y*r_x - w_x*r_y)
	        (*qw[k]).z -= DLM_w.y[] * ( lambdapos.x - (*gci[k]).center.x ) 
			- DLM_w.x[] * ( lambdapos.y - (*gci[k]).center.y ); 

	        foreach_dimension()
	          (*gci[k]).center.x -= DLM_periodic_shift.x[];	
#             endif
            }
          }
        }
#     endif  
#   endif
    

#   if DLM_Moving_particle
#     if _MPI /* _MPI Reduction */
        // Reduce to master
        // Pack data
        counter = 0;  
        for (int k = 0; k < NPARTICLES; k++) 
        {
#         if TRANSLATION
            vpartbuf[counter] = (*qU[k]).x;
            vpartbuf[counter+1] = (*qU[k]).y;
#           if dimension == 3	      
              vpartbuf[counter+2] = (*qU[k]).z;
              counter += 3;
#           else
              counter += 2;
#           endif    
#         endif
#         if ROTATION
#           if dimension == 3
              vpartbuf[counter] = (*qw[k]).x;
              vpartbuf[counter+1] = (*qw[k]).y;    
              vpartbuf[counter+2] = (*qw[k]).z;
              counter += 3; 
#           else
              vpartbuf[counter] = (*qw[k]).z;
	      counter += 1; 
#           endif 	    
#         endif
        }
  
        // Perform reduction on Master
        MPI_Reduce( pid() ? vpartbuf : MPI_IN_PLACE, vpartbuf, 
      		npartdata*NPARTICLES, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );	

        if ( pid() == 0 )
        {
          // Unpack data    
          counter = 0;
          for (int k = 0; k < NPARTICLES; k++) 
          {
#           if TRANSLATION
              (*qU[k]).x = vpartbuf[counter];
              (*qU[k]).y = vpartbuf[counter+1];
#             if dimension == 3	          
                (*qU[k]).z = vpartbuf[counter+2];      
                counter += 3;
#             else
                counter += 2;
#             endif	          
#           endif
#           if ROTATION
#             if dimension == 3
                (*qw[k]).x = vpartbuf[counter];
                (*qw[k]).y = vpartbuf[counter+1];    
                (*qw[k]).z = vpartbuf[counter+2];    
                counter += 3; 
#             else
                (*qw[k]).z = vpartbuf[counter];    
                counter += 1;
#             endif	      
#           endif
          } 	   	
#     endif  /* end of _MPI Reduction */    
    
          // Perform the inversion (on master when in MPI)
          for (int k = 0; k < NPARTICLES; k++) 
          {
#           if TRANSLATION
              /* Solution of M * DLMFD_couplingfactor * tU / dt = qU */
              /* tU = ( dt * qU ) / ( DLMFD_couplingfactor * M ) */
              foreach_dimension()              
                (*tU[k]).x = ( (*qU[k]).x * dt ) / 
	      		( p[k].DLMFD_couplingfactor * p[k].M );	   
#           endif
#           if ROTATION
              /* Solution of Ip * DLMFD_couplingfactor * tw / dt = qw */
              /* tw = ( dt / ( DLMFD_couplingfactor ) * Ip_inv * qw 
              /* where Ip_inv is the inverse of Ip */      
#             if dimension == 3              
	        (*tw[k]).x = ( dt / p[k].DLMFD_couplingfactor ) * 
    			( (p[k].Ip_inv)[0][0] * (*qw[k]).x 
			+ (p[k].Ip_inv)[0][1] * (*qw[k]).y 
    			+ (p[k].Ip_inv)[0][2] * (*qw[k]).z );
                (*tw[k]).y = ( dt / p[k].DLMFD_couplingfactor ) * 
    			( (p[k].Ip_inv)[1][0] * (*qw[k]).x 
			+ (p[k].Ip_inv)[1][1] * (*qw[k]).y 
    			+ (p[k].Ip_inv)[1][2] * (*qw[k]).z );
                (*tw[k]).z = ( dt / p[k].DLMFD_couplingfactor ) * 
    			( (p[k].Ip_inv)[2][0] * (*qw[k]).x 
			+ (p[k].Ip_inv)[2][1] * (*qw[k]).y 
    			+ (p[k].Ip_inv)[2][2] * (*qw[k]).z );
#             else
                (*tw[k]).z = ( dt / p[k].DLMFD_couplingfactor ) * 
    			(p[k].Ip_inv)[2][2] * (*qw[k]).z ;
#             endif			
#           endif
          }  

#     if _MPI /* _MPI Broadcast */
          // Broadcast tU and tw
          // Pack data
          counter = 0;  
          for (int k = 0; k < NPARTICLES; k++) 
          {
#           if TRANSLATION
              vpartbuf[counter] = (*tU[k]).x;
              vpartbuf[counter+1] = (*tU[k]).y;    
#             if dimension == 3
                vpartbuf[counter+2] = (*tU[k]).z;
                counter += 3;
#             else
                counter += 2;
#             endif	          
#           endif
#           if ROTATION
#             if dimension == 3
                vpartbuf[counter] = (*tw[k]).x;
                vpartbuf[counter+1] = (*tw[k]).y;    
                vpartbuf[counter+2] = (*tw[k]).z;
                counter += 3;
#             else
                vpartbuf[counter] = (*tw[k]).z;
	        counter += 1;
#             endif	       
#           endif
          }
        } /* End of "if pid() == 0" */    
  
        // Perform the broadcast from the master to the other processes
        MPI_Bcast( vpartbuf, npartdata*NPARTICLES, 
    		MPI_DOUBLE, 0, MPI_COMM_WORLD );

        // Unpack data    
        counter = 0;
        for (int k = 0; k < NPARTICLES; k++) 
        {
#         if TRANSLATION
            (*tU[k]).x = vpartbuf[counter];
            (*tU[k]).y = vpartbuf[counter+1];      
#           if dimension == 3
              (*tU[k]).z = vpartbuf[counter+2];      
              counter += 3;
#           else	
              counter += 2;
#           endif        
#         endif
#         if ROTATION
#           if dimension == 3
              (*tw[k]).x = vpartbuf[counter];
              (*tw[k]).y = vpartbuf[counter+1];    
              (*tw[k]).z = vpartbuf[counter+2];    
              counter += 3; 
#           else
              (*tw[k]).z = vpartbuf[counter];    
              counter += 1; 	    
#           endif
#         endif
        }  
#     endif /* end of _MPI Broadcast */ 
#   endif /* end of DLM_Moving_particle */


    /* (2) Invert L*t = qu with L=rho*dV/dt*I */
#   if DLMFD_OPT
      foreach_cache( Traversal_uqutu ) 
#   else
      foreach()
#   endif
        foreach_dimension() 
	  tu.x[] = qu.x[] * dt / ( rho_f * dlmfd_dv() );

         
    /* (3) Compute residual y = M*t, the residual vector 
       y = <alpha, tu-(tU+tw^GM)>_P(t) 
       Sign error in eq (A.9) in J.Eng. Math, 2011 
       Written -M*t, should be +M*t */    

#   if debugInterior == 0
      /* Interior points: y = M*t  */
      /* Interior points: y = M_u*tu + M_U*tU + M_w*tw */
      /* So y = <alpha, tu>_P(t) - <alpha, tU>_P(t) - <alpha, tw^r_GM>_P(t) */
      for (int k = 0; k < NPARTICLES; k++) 
      {
        foreach_cache((*Interior[k])) 
        {
	  if ((flagfield[]  < 1) && ((int)index_lambda.y[] == k)) 
	  {
	    foreach_dimension() 
	      DLM_v.x[] = tu.x[];

#           if DLM_Moving_particle
#             if TRANSLATION
	        foreach_dimension() 
		  DLM_v.x[] -= (*tU[k]).x;
#             endif
#             if ROTATION
	        /* Modify temporarily the particle center position for periodic 
	        boundary condition */
	        foreach_dimension()
	          (*gci[k]).center.x += DLM_periodic_shift.x[];

#               if dimension == 3	
	          // -(t_y*r_z - t_z*r_y)
	          DLM_v.x[] -= (*tw[k]).y * (z - (*gci[k]).center.z ) 
	  		- (*tw[k]).z * ( y - (*gci[k]).center.y ); 
	          // -(t_z*r_x - t_x*r_z)
	          DLM_v.y[] -= (*tw[k]).z * ( x - (*gci[k]).center.x ) 
	  		- (*tw[k]).x * ( z - (*gci[k]).center.z );
#               endif			 
	        // -(t_x*r_y - t_y*r_x)
	        DLM_v.z[] -= (*tw[k]).x * ( y - (*gci[k]).center.y ) 
	  		- (*tw[k]).y * ( x - (*gci[k]).center.x ); 

	        foreach_dimension()
	          (*gci[k]).center.x -= DLM_periodic_shift.x[];	
#             endif
#           endif
	  }
        }
      } 
#   endif
    
#   if debugBD == 0
      /* Boundary points: y = M*t */
      /* Boundary points: y = M_u*tu + M_U*tU + M_w*tw */
      /* So y = <alpha, tu>_P(t) - <alpha, tU>_P(t) - <alpha, tw^r_GM>_P(t) */
      boundary ((scalar*){tu});
#     if DLMFD_OPT
        // Use of the fast loop for the computations of 
        // M_u * tu = <alpha, tu>_P(t) over the boundary points
        int endloopv = RUloop.nv, pos = 0, endlooptu = 0;  
        for (int k=0; k<endloopv; ++k )
        {
          sum.x = 0; sum.y = 0; sum.z = 0;
          endlooptu = RUloop.ndof[k] + pos;
          for (int l=pos;l<endlooptu;++l)
          {
            weight = RUloop.weight[l];
	    sum.x += weight * *(RUloop.tux[l]);
	    sum.y += weight * *(RUloop.tuy[l]);
	    sum.z += weight * *(RUloop.tuz[l]);
          }
          *(RUloop.dlmvx[k]) = sum.x;
          *(RUloop.dlmvy[k]) = sum.y;      
          *(RUloop.dlmvz[k]) = sum.z; 
          pos += RUloop.ndof[k];              
        }
#     else
        for (int k = 0; k < NPARTICLES; k++) 
        {
          particle * pp = &p[k];
      
          foreach_cache((*Boundary[k])) 
          {
	    if (index_lambda.x[] > -1 && ((int)index_lambda.y[] == pp->pnum)) 
	    {
	      lambdacellpos.x = x; 
	      lambdacellpos.y = y;
	      lambdapos.x = (*sbm[k]).x[(int)index_lambda.x[]];
	      lambdapos.y = (*sbm[k]).y[(int)index_lambda.x[]];
#             if dimension == 3
	        lambdacellpos.z = z;
	        lambdapos.z = (*sbm[k]).z[(int)index_lambda.x[]];
#             endif

	      foreach_dimension() 
	        sum.x = 0; 
		
	      double testweight = 0.;

	      foreach_dimension()
	        ppshift.x = DLM_periodic_shift.x[];

	      foreach_neighbor() 
	      {
	        if ( level == depth() ) 
	        {
	          weightcellpos.x = x; 
		  weightcellpos.y = y; 
#                 if dimension == 3
	            weightcellpos.z = z;
#                 endif
	     
	          weight = reversed_weight( pp, weightcellpos, lambdacellpos, 
	      	  	lambdapos, Delta, ppshift );
	          testweight += weight;
	      	      	      
	          foreach_dimension() 
		    sum.x += weight * tu.x[];
	        }
	      }

	      if ( testweight < 1. - 0.000001 || testweight > 1. + 0.000001 )
	        printf("testweight = %f\n",testweight);
	    
	      foreach_dimension() 
	        DLM_v.x[] = sum.x;
	    }
          }
        }
#     endif
    
#     if DLM_Moving_particle
        for (int k = 0; k < NPARTICLES; k++) 
        {
          particle * pp = &p[k];
      
          foreach_cache((*Boundary[k])) 
          {
	    if (index_lambda.x[] > -1 && ((int)index_lambda.y[] == pp->pnum)) 
	    {    
	      lambdapos.x = (*sbm[k]).x[(int)index_lambda.x[]];
	      lambdapos.y = (*sbm[k]).y[(int)index_lambda.x[]];
#             if dimension == 3
	        lambdapos.z = (*sbm[k]).z[(int)index_lambda.x[]];
#             endif

#             if TRANSLATION
	        foreach_dimension() 
		  DLM_v.x[] -= (*tU[k]).x;
#             endif
#             if ROTATION
	        /* Modify temporarily the particle center position for periodic 
	        boundary condition */
	        foreach_dimension()
	          (*gci[k]).center.x += DLM_periodic_shift.x[];

#               if dimension == 3	
	          // -(t_y*r_z - t_z*r_y)
	          DLM_v.x[] -= (*tw[k]).y * ( lambdapos.z - (*gci[k]).center.z )
	  		- (*tw[k]).z * ( lambdapos.y - (*gci[k]).center.y ); 
	          // -(t_z*r_x - t_x*r_z)
	          DLM_v.y[] -= (*tw[k]).z * ( lambdapos.x - (*gci[k]).center.x )
	  		- (*tw[k]).x * ( lambdapos.z - (*gci[k]).center.z );
#               endif			 
	        // -(t_x*r_y - t_y*r_x)
	        DLM_v.z[] -= (*tw[k]).x * ( lambdapos.y - (*gci[k]).center.y ) 
	  		- (*tw[k]).y * ( lambdapos.x - (*gci[k]).center.x ); 

	        foreach_dimension()
	          (*gci[k]).center.x -= DLM_periodic_shift.x[];	
#             endif
	    }
          }
        }
#      endif
#   endif
     
    /* (4) Compute alpha = r^(k-1)*r^(k-1) / DLM_w.y */
    DLM_nr2_km1 = DLM_nr2;
    DLM_wv = 0.;

#   if DLMFD_OPT
      foreach_cache( Traversal_rvwlambda )
#   else
      foreach()
#   endif
#     if dimension == 3
        DLM_wv += DLM_w.x[]*DLM_v.x[] + DLM_w.y[]*DLM_v.y[] 
      		+ DLM_w.z[]*DLM_v.z[];
#     else
        DLM_wv += DLM_w.x[]*DLM_v.x[] + DLM_w.y[]*DLM_v.y[];
#     endif		
    
#   if _MPI
      mpi_all_reduce( DLM_wv, MPI_DOUBLE, MPI_SUM );
#   endif

    DLM_alpha = DLM_nr2_km1 / DLM_wv;
        
    /* (5) Update lambda = lambda - alpha*DLM_w and r = r - alpha*y */
#   if DLMFD_OPT
      foreach_cache( Traversal_rvwlambda )
#   else
      foreach()
#   endif
        foreach_dimension()
	{
	  DLM_lambda.x[] -= DLM_alpha * DLM_w.x[];
          DLM_r.x[] -= DLM_alpha * DLM_v.x[];
        }  

    /* (7) Update u = u + alpha*t */
#   if DLMFD_OPT
      foreach_cache( Traversal_uqutu )
#   else
      foreach()
#   endif
        foreach_dimension()
          u.x[] += DLM_alpha * tu.x[];

#   if DLM_Moving_particle
      for (int k = 0; k < NPARTICLES; k++) 
      {
#       if TRANSLATION
          foreach_dimension() 
	    (*U[k]).x += DLM_alpha * (*tU[k]).x;
#       endif
#       if ROTATION
#         if dimension == 3   
            (*w[k]).x += DLM_alpha * (*tw[k]).x;
            (*w[k]).y += DLM_alpha * (*tw[k]).y;
            (*w[k]).z += DLM_alpha * (*tw[k]).z;
#         else
            (*w[k]).z += DLM_alpha * (*tw[k]).z;
#         endif
#       endif
      }
#   endif
   
    /* (8) Compute beta = nr2^k / nr2^(k-1) */
    DLM_nr2 = 0.;
    
    /* (6) compute norme ||r^k||^2 =  DLM_nr2 */
#   if DLMFD_OPT
      foreach_cache( Traversal_rvwlambda )
#   else
      foreach()
#   endif
#     if dimension == 3
        DLM_nr2 += sq(DLM_r.x[]) + sq(DLM_r.y[]) + sq(DLM_r.z[]);
#     else
        DLM_nr2 += sq(DLM_r.x[]) + sq(DLM_r.y[]);
#     endif	
    
#   if _MPI
      mpi_all_reduce( DLM_nr2, MPI_DOUBLE, MPI_SUM );
#   endif
    
    DLM_beta = DLM_nr2 / DLM_nr2_km1;

    /* (9) Update DLM_w = r + beta*DLM_w */
    /* Sign error in eq (A.15) of J.Eng. Math, 2011 */
    /* Written -beta*DLM_w, should be +beta*DLM_w */
#   if DLMFD_OPT
      foreach_cache( Traversal_rvwlambda )
#   else
      foreach()
#   endif
        foreach_dimension()
          DLM_w.x[] = DLM_r.x[] + DLM_beta * DLM_w.x[];

  }  /* End of Iterative loop */




  // Once algorithm has converged
  if ( pid() == 0 )
  {
    printf( "      niter = %d residual = %8.5e\n", ki, sqrt(DLM_nr2) );
    fprintf( converge,"%d \t %d \t \t %10.8e\n", i, ki, sqrt(DLM_nr2) );
    fflush( converge );
  }

  
  /* Compute the explicit term here */
# if DLM_alpha_coupling 
    boundary ((scalar*) {DLM_lambda});

    for (int k = 0; k < NPARTICLES; k++) 
    {
#     if debugBD == 0
        particle * pp = &p[k];
    
        foreach_cache ((*Boundary[k])) 
        {
          weightcellpos.x = x; 
	  weightcellpos.y = y;
#         if dimension == 3
            weightcellpos.z = z;
#         endif

          foreach_dimension() 
	    sum.x = 0.; 

          foreach_neighbor() 
          {
	    if (((int)index_lambda.x[] > -1) && (level == depth()) && 
		is_leaf(cell) && ((int)index_lambda.y[]) == k) 
	    {
	      lambdacellpos.x = x;
	      lambdacellpos.y = y;
	      lambdapos.x = (*sbm[k]).x[(int)index_lambda.x[]];
	      lambdapos.y = (*sbm[k]).y[(int)index_lambda.x[]];	
#             if dimension == 3
	        lambdacellpos.z = z;
	        lambdapos.z = (*sbm[k]).z[(int)index_lambda.x[]];
#             endif

	      foreach_dimension()
	        ppshift.x = DLM_periodic_shift.x[];

	      weight = reversed_weight (pp, weightcellpos, lambdacellpos, 
	  	lambdapos, Delta, ppshift);
	  
	      foreach_dimension()
	        sum.x += weight*DLM_lambda.x[];
	    }
          }

          foreach_dimension()
	    DLM_explicit.x[] += sum.x;
        }
#     endif
    
#     if debugInterior == 0
        foreach_cache ((*Interior[k])) 
        {
          if ((flagfield[]  < 1) && ((int)index_lambda.y[] == k))
	    foreach_dimension() 
	      DLM_explicit.x[] = DLM_lambda.x[];
        }
#     endif
    }
    boundary((scalar *) {DLM_explicit});
# endif


# if DLMFD_OPT
    // Free fast loops
    free_BPFastLoop_LambdaMom( &LMloop );
    free_BPFastLoop_ResU( &RUloop );

    // Reset DLM work vectors to 0 for next call using the global caches
    foreach_cache( Traversal_rvwlambda ) 
      foreach_dimension() 
      {
        DLM_r.x[] = 0.; 
        DLM_w.x[] = 0.; 
        DLM_v.x[] = 0.; 
      }
    boundary((scalar *) {DLM_w});      
    
    // We do not need to nullify qu as it will be assigned on all cells
    // at the next call of DLMFD_subproblem  
    foreach_cache( Traversal_uqutu ) 
      foreach_dimension()
        tu.x[] = 0.;      
    boundary((scalar *) {tu});  

  
    // Free DLMFD global caches
    free( Traversal_rvwlambda.p );
    free( Traversal_uqutu.p );  
# endif

  
  // To guarantee that u is the same on all sub-domains
  boundary ((scalar*) {u});


  /* Timers and Timings */
  /* give 1 as number of cell so that timer_timing does not compute the total
     number of cells
     the total number of cell is used to compute the speed of the solver 
     which does not make sens for a single step here, we will compute it 
     globally at the end of the run */
  timing dlmfd_timing = timer_timing (dlmfd_timer, 1, 1, mpitimings);

  dlmfd_globaltiming.cpu += dlmfd_timing.cpu;
  dlmfd_globaltiming.real += dlmfd_timing.real;
  dlmfd_globaltiming.min += dlmfd_timing.min;
  dlmfd_globaltiming.avg += dlmfd_timing.avg;
  dlmfd_globaltiming.max += dlmfd_timing.max;
}




/** Initialize all DLMFD fields */
//----------------------------------------------------------------------------
void initialize_DLMFD_fields_to_zero( void )
//----------------------------------------------------------------------------
{
  foreach()
  {
    flagfield[] = 0.;
    flagfield_mailleur[] = 0.;
    foreach_dimension()
    {
      DLM_lambda.x[] = 0. ;
      index_lambda.x[] = 0. ;
      DLM_periodic_shift.x[] = 0. ;
      DLM_r.x[] = 0. ;
      DLM_w.x[] = 0. ;
      DLM_v.x[] = 0. ;
      qu.x[] = 0. ;
      tu.x[] = 0. ;
#     if DLM_alpha_coupling
        DLM_explicit.x[] = 0. ;
#     endif    
    }
  }
}
