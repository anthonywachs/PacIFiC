/**
# General functions for the DLMFD implementation.
*/

/** Define NSDF, the number of significant digits after the decimal point
to output data in files in text mode. NSDF cannot be lower than 3, larger 
than 15 and 9 (for formatting reasons). */
# ifndef NSDF
#   define NSDF 7
# else 
#   if ( NSDF == 9 )
#     undef NSDF
#     define NSDF 8
#   else
#     if ( NSDF < 3 )
#       undef NSDF
#       define NSDF 3
#     else
#       if ( NSDF > 15 )
#         undef NSDF
#         define NSDF 15
#       endif
#     endif
#   endif
# endif




/** Define the factor alpha (generally between 1 and 2) that is involved 
in the inter-boundary point distance on the rigid body surface. */
# ifndef INTERBPCOEF
#   define INTERBPCOEF 2.
# endif 




/** Different rigid body shapes supported */      
enum RigidBodyShape {
  SPHERE,
  CIRCULARCYLINDER2D,
  CUBE,
  TETRAHEDRON
};


 

/** Structure for the coordinates of a rigid body boundary. */
typedef struct {
  double *x;
  double *y;
  double *z;
  int m;
  int nm;
} SolidBodyBoundary;




/** Additional geometric parameters for polygons/polyhedrons */
typedef struct {
  int allPoints, allFaces;
  double ** cornersCoord;
  long int ** cornersIndex;
  long int * numPointsOnFaces;
  
  /* Special Cube: 3 principal vectors */
  coord u1, v1, w1;
  coord mins, maxs;
  
} PolyGeomParameter;





/** Rigid body geometric parameters */
typedef struct {
  coord center;
  double radius;
  int ncorners;
  PolyGeomParameter* pgp;  
} GeomParameter;




/** Rigid body parameters for the toy granular solver */
typedef struct {
  double kn, en, vzero, wished_ratio;
  coord normalvector;
  GeomParameter gnm1;  
} ToyGSParameter;




/** Set of parameters describing a rigid body (also named particle) */
typedef struct {
  size_t pnum;
  enum RigidBodyShape shape;  
  SolidBodyBoundary s;
  GeomParameter g;
  double M, Ip[6], rho_s, Vp, DLMFD_couplingfactor, RotMat[3][3];  
# if DLM_Moving_particle
    ToyGSParameter *toygsp;
    coord gravity;
    double Ip_inv[3][3];
    coord addforce;    
#   if TRANSLATION
      coord U, Unm1, qU, tU;
#   endif
#   if ROTATION
      coord w, wnm1, qw, tw;
#   endif
# else
    coord imposedU, imposedw;    
# endif
  Cache Interior;
  Cache reduced_domain;
  long tcells, tmultipliers;
} particle;




# include "CircularCylinder2D.h"
# include "Sphere.h"
# include "Cube.h"
# include "Tetrahedron.h"


/** Allocates memory for m points in the SolidBodyBoundary structure. */
//----------------------------------------------------------------------------
void allocate_SolidBodyBoundary( SolidBodyBoundary* sbm, const int m ) 
//----------------------------------------------------------------------------
{
  sbm->x = (double*) calloc( m, sizeof(double)); 
  sbm->y = (double*) calloc( m, sizeof(double));
# if dimension == 3  
    sbm->z = (double*) calloc( m, sizeof(double));
# else
    sbm->z = NULL;
# endif

  sbm->m = m;
}




/** Re-allocates memory for m points in the SolidBodyBoundary structure. */
//----------------------------------------------------------------------------
void reallocate_SolidBodyBoundary( SolidBodyBoundary* sbm, const int m ) 
//----------------------------------------------------------------------------
{
  sbm->x = (double*) realloc (sbm->x, m*sizeof(double)); 
  sbm->y = (double*) realloc (sbm->y, m*sizeof(double));  
# if dimension == 3 
    sbm->z = (double*) realloc (sbm->z, m*sizeof(double));
# endif    
  
  sbm->m = m;
}




/** Frees memory associated to the points in the SolidBodyBoundary structure. */
//----------------------------------------------------------------------------
void free_SolidBodyBoundary( SolidBodyBoundary* sbm ) 
//----------------------------------------------------------------------------
{
  free( sbm->x ); sbm->x = NULL;
  free( sbm->y ); sbm->y = NULL;
# if dimension == 3 
    free( sbm->z ); sbm->z = NULL;
# endif 
}




/** Allocates an initial Cache structure */
//----------------------------------------------------------------------------
void initialize_and_allocate_Cache( Cache* p ) 
//----------------------------------------------------------------------------
{
  p->n = 0;
  p->nm = 1;
  p->p = (Index *) calloc( 1, sizeof(Index) );
}




/** Frees the particle data that were dynamically allocated */
//----------------------------------------------------------------------------
void free_particles( particle* pp, const int n ) 
//----------------------------------------------------------------------------
{
  for (size_t k=0;k<n;k++) 
  {        
    // Free the boundary point coordinate arrays
    SolidBodyBoundary* sbm = &(pp[k].s);
    free_SolidBodyBoundary( sbm );
    
    // Free the caches 
    Cache* c = &(pp[k].Interior);
    free( c->p );
    c->p = NULL;
    c = &(pp[k].reduced_domain);
    free( c->p );
    c->p = NULL;    
    
    // Free the toy granular solver parameter structure
# if DLM_Moving_particle
    if ( pp[k].toygsp )
    {
      free( pp[k].toygsp );
      pp[k].toygsp = NULL;
    }     
# endif

    // Free the additional geometric features of the particle
    switch ( pp[k].shape )
    {
      case SPHERE:
        free_Sphere( &(pp[k].g) );
	break;
	  
      case CIRCULARCYLINDER2D:
        free_CircularCylinder2D( &(pp[k].g) );
	break;
	  
      case CUBE:
        free_Cube( &(pp[k].g) );
	break;

      case TETRAHEDRON:
        free_Tetrahedron( &(pp[k].g) );
	break;
	  
      default:
        fprintf( stderr,"Unknown Rigid Body shape !!\n" );
    }                               
  }
}




/** Prints data of a particle */
//----------------------------------------------------------------------------
void print_particle( particle const* pp, char const* poshift )
//----------------------------------------------------------------------------
{
  if ( pid() == 0 ) 
  {
    printf( "%sNumber = %lu\n", poshift, pp->pnum ); 
    printf( "%sShape = ", poshift );
    switch ( pp->shape )
    {
      case SPHERE:
        printf( "SPHERE" );
        break;
	  
      case CIRCULARCYLINDER2D:
        printf( "CIRCULARCYLINDER2D" );
        break;
	  
      case CUBE:
        printf( "CUBE" );
        break;
	
      case TETRAHEDRON:
        printf( "TETRAHEDRON" );
	break;	
	  
      default:
        fprintf( stderr,"Unknown Rigid Body shape !!\n" );
    }
    printf( "\n" );
    printf( "%sCenter of mass = %e %e", poshift, pp->g.center.x, 
    	pp->g.center.y );
#   if dimension == 3
      printf( " %e", pp->g.center.z );
#   endif
    printf( "\n" );
    printf( "%sRadius = %e\n", poshift, pp->g.radius );  
    printf( "%sMass = %e\n", poshift, pp->M ); 
    printf( "%sVolume = %e\n", poshift, pp->Vp );     
    printf( "%sDensity = %e\n", poshift, pp->rho_s ); 
#   if DLM_Moving_particle
#     if TRANSLATION
        printf( "%sTranslational velocity = %e %e", poshift, pp->U.x, pp->U.y );
#       if dimension == 3
          printf( " %e", pp->U.z );
#       endif
        printf( "\n" );
#     endif
#     if ROTATION
        printf( "%sAngular velocity = ", poshift );
#       if dimension == 3
          printf( "%e %e", pp->w.x, pp->w.y );
#       endif
        printf( " %e\n", pp->w.z );
#     endif
#   else
      printf( "%sImposed translational velocity = %e %e", 
    	poshift, pp->imposedU.x, pp->imposedU.y );
#     if dimension == 3
        printf( " %e", pp->imposedU.z );
#     endif
      printf( "\n" );
      printf( "%sImposed angular velocity = ", poshift );
#     if dimension == 3
        printf( "%e %e", pp->imposedw.x, pp->imposedw.y );
#     endif
      printf( " %e\n", pp->imposedw.z );    
#   endif
  }
  int intpts = pp->Interior.n;
  int bdpts = pp->reduced_domain.n;  
# if _MPI
    mpi_all_reduce( intpts, MPI_INTEGER, MPI_SUM );
    mpi_all_reduce( bdpts, MPI_INTEGER, MPI_SUM );
# endif
  if ( pid() == 0 )
  {   
    printf( "%sNumber of interior DLM points = %d\n", poshift, intpts ); 
    printf( "%sNumber of boundary DLM points = %d\n", poshift, pp->s.m ); 
    printf( "%sNumber of cells in boundary DLM point stencils = %d\n", poshift, 
    	bdpts );
  }
} 




/** Prints all particle data */
//----------------------------------------------------------------------------
void print_all_particles( particle const* allpart, char const* oshift )
//----------------------------------------------------------------------------
{
  if ( pid() == 0 ) printf( "%sTotal number of particles = %d\n", oshift, 
  	NPARTICLES );
  char poshift[10]="   ";
  strcat( poshift, oshift );
  for (size_t k=0;k<NPARTICLES;k++)
  { 
    if ( pid() == 0 ) printf( "%sParticle %lu\n", oshift, k );    
    print_particle( &(allpart[k]), &poshift[0] );
  }
} 




/** Ceates the scalar field where the cells contain the
indices of the dlmfd points. If there is no dlmfd point it is
tagged -1. */
//----------------------------------------------------------------------------
void create_index_lambda_scalar (const SolidBodyBoundary dlm_bd, 
	vector Index_lambda, const int kk) 
//----------------------------------------------------------------------------
{  
  Point lpoint;
  int i;
  Cache *fdlocal;
   
  fdlocal = (Cache *){calloc(dlm_bd.m, sizeof(Cache))};
  
  for (i = 0; i < dlm_bd.m; i++) {

#if dimension == 2 
    lpoint = locate(dlm_bd.x[i], dlm_bd.y[i]);
#elif dimension == 3
    lpoint = locate(dlm_bd.x[i], dlm_bd.y[i], dlm_bd.z[i]);
#endif    
 /*  */
    if ((lpoint.level)  == depth()) {
	/* Create a cache for each fictitious domain's boundary point */
	cache_append(&fdlocal[i], lpoint, 0);
	
      	foreach_cache(fdlocal[i]) {
	  /* Tag cell only if it was not tagged by another particle */
	  if (Index_lambda.x[] < 0){
	    Index_lambda.x[] = i;
	    Index_lambda.y[] = kk;
	    if (level != depth()) {
	      printf("On thread %d, point dlmfd %d at (%f, %f, %f) is in a"
	      	" cell that has not the maximum level of refinnement %d, it "
		"is on level %d \n", pid(), i, x, y, z, depth(), level);
	    }
	  }
	}
    }
  }
   
  for (i = 0; i < dlm_bd.m; i++) {
    free(fdlocal[i].p);
  }
  
  free (fdlocal);
  
  boundary ((scalar*) {Index_lambda});
}





#include "DLMFD_Weight_functions.h"

/** Returns the weight associated to a Lagrange 
multipliers for a cell within a foreach() loop (ie for a "real" or local cell) 
associated to a Lagrange multiplier in it's neighborhood. The neighboor can be 
in another process (i.e be a ghost cell for the current thread). */
//----------------------------------------------------------------------------
double reversed_weight (particle * pp, const coord weightcellpos, 
	const coord lambdacellpos, const coord lambdapos, const double delta, 
	const coord pshift) 
//----------------------------------------------------------------------------
{
  /* this function has to be embedded within a double foreach_cache() and 
  foreach_neighbor() loop  */
  coord rel = {0, 0, 0};
  coord relnl = {0, 0, 0};
  int NCX = 0, CX = 0, weight_id = 0;
  size_t goflag = 0;
  double weight = 0.;
  GeomParameter gcb = pp->g; 
  
  /* compute relative vector X_boundary - X_local = rel from the cell 
  (containning the boundary) position to the boundary's (analytical) position 
  i.e. vector goes from first argument to the second */
  compute_relative_vector (lambdacellpos, lambdapos, &rel);

  /* reset dials integers */
  NCX = 0; CX = 0; weight_id = 0; goflag = 0;

  /* assign first normal (gives the quadrant) */ 
  assign_dial (rel, &CX);
  
  GeomParameter gcbdum;
  gcbdum = gcb;
	
  foreach_dimension()
    gcbdum.center.x += pshift.x;

  /* assign fictitious-boundary's normal (use boundary's analytical position)*/
  assign_dial_fd_boundary (pp, lambdapos, gcbdum, delta, &NCX);
  
  /* compute relative vector from the neighbor cell to the local cell */
  compute_relative_vector (weightcellpos, lambdacellpos , &relnl);

  /* Assign weight ids. */
  assign_weight_id_quad_outward (NCX, CX, relnl, delta, &weight_id, &goflag);
    
  /* Compute weight */
  if ( goflag == 1 ) 
    weight = compute_weight_Quad (weight_id, lambdapos, lambdacellpos, 
    	NCX, CX, delta);
  
  return weight;
}




//----------------------------------------------------------------------------
void remove_too_close_multipliers(particle * p, vector index_lambda) 
//----------------------------------------------------------------------------
{   
  for (int k = 0; k < NPARTICLES; k++) {

    SolidBodyBoundary dlm_lambda_to_desactivate;
    int allocated = 1;
    allocate_SolidBodyBoundary (&dlm_lambda_to_desactivate, allocated*BSIZE);
    int countalloc = 0;
    int other_part;
    int is_in_other_particle;
    foreach_level(depth()) {
      
      int direct_neigh = 0;  is_in_other_particle = 0; other_part = -1;
      if (((int)index_lambda.x[] > -1)  && level == depth() && is_leaf(cell) 
      	&& (p[k].pnum == (int)index_lambda.y[])) {

	/* check here if this multiplier is not in another particle's
	   domain, if yes desactive it */
	
	for (int l = 0; l < NPARTICLES; l++) {
	  if (l != p[k].pnum) {
      	    particle * other_particle = &(p[l]);
      	    /* printf("thread %d, this is particle %zu checkin on particle 
	    %zu iteration %d\n", pid(), p[k].pnum, other_particle->pnum, l); */

      	    /* Check particle's type */
      	    if ( other_particle->shape == CUBE ) {
      	      compute_principal_vectors_Cube( other_particle );
      	      coord u = other_particle->g.pgp->u1;
      	      coord v = other_particle->g.pgp->v1;
      	      coord w = other_particle->g.pgp->w1;
      	      coord mins = other_particle->g.pgp->mins;
      	      coord maxs = other_particle->g.pgp->maxs;

      	      /* current cell's position */
      	      coord checkpt = {x, y, z};
      	      if ( is_in_Cube( &u, &v, &w, &mins, &maxs, &checkpt ) ) {
      		is_in_other_particle = 1; other_part =  other_particle->pnum;
		break;
      	      }
      	    }
      	    /* if sphere */
      	    else {
      	      GeomParameter gp = other_particle->g;
      	      if (is_in_Sphere (x, y, z, gp)) {
      		is_in_other_particle = 1; other_part =  other_particle->pnum;
		break;
      	      }
      	    }
      	  }
	  
	}

	if (is_in_other_particle) {
	  /* printf("thread %d, this is particle %zu which has a cell 
	  in particle %d\n", pid(), p[k].pnum, other_part); */
	  index_lambda.x[] = -1; index_lambda.y[] = other_part;
	}
	
	/* Check if two (or more) Lagrange multipliers (from different
	   particles) are in the same neighborhoods (i.e a 5x5 stencil
	   in 2D), if yes desactivate them. Otherwise the cells of
	   their stencils may be located in each other's domain */
	direct_neigh = 0;
	foreach_neighbor() {
	  if (((int)index_lambda.x[] > -1)  && level == depth() 
	  	&& is_leaf(cell) && (p[k].pnum != (int)index_lambda.y[])) {

	    direct_neigh = 1;
	    /* We may catch and identical multiplier more than once here */
	    /* Add the multiplier(s) indices to a list because MPI imposes 
	    that we can't modify the cell within a foreach_neighbor loop */
	    if (countalloc >= allocated*BSIZE) {
	      allocated ++;
	      reallocate_SolidBodyBoundary (&dlm_lambda_to_desactivate, 
	      	allocated*BSIZE);
	    }
	    dlm_lambda_to_desactivate.x[countalloc] = x;
	    dlm_lambda_to_desactivate.y[countalloc] = y;
#if dimension == 3
	    dlm_lambda_to_desactivate.z[countalloc] = z;
#endif
	    countalloc ++;
	  }
	}
	/* Desactivate the local multiplier here  */
	if (direct_neigh) {
	index_lambda.x[] = -1; index_lambda.y[] = -1;
	}
      }
    }
    /* if (countalloc > 0) */
    /*   printf("there is %d points to desactivate on thread %d\n", 
    countalloc, pid()); */

    Cache c = {0};
    Point lpoint;
    
#if _MPI
    int size = npe();
    int counts[size];
    
    /* Each thread tells the root how many multipliers it holds  */
    MPI_Gather(&countalloc, 1, MPI_INT, &counts, 1, MPI_INT, 0, 
    	MPI_COMM_WORLD);

    /* if (pid() == 0) */
    /*   for (int i = 0; i < size; i++) */
    /* 	printf("particle %d, this is root receiving %d elements 
    by thread %d\n",k, counts[i], i); */

    /* Displacements in the receive buffer for MPI_GATHERV */
    int disps[size];
     /* Displacement for the first chunk of data - 0 */
    for (int i = 0; i < size; i++)
      disps[i] = (i > 0) ? (disps[i-1] + counts[i-1]) : 0;

    double * alldatax = NULL;
    double * alldatay = NULL;
#if dimension == 3 
    double * alldataz = NULL;
#endif

    int m = 0;
    /* Threads 0 allocates and compute the total number of multipliers */
    if (pid() == 0) {
      /* disps[size-1] + counts[size-1] == total number of elements */
      /* printf("thread %d : disps[size-1]+counts[size-1]  = %d\n", pid(), 
      disps[size-1]+counts[size-1]); */
      m = disps[size-1] + counts[size-1];
      alldatax = (double*) calloc (m, sizeof(double));
      alldatay = (double*) calloc (m, sizeof(double));
#if dimension == 3
      alldataz = (double*) calloc (m, sizeof(double));
#endif
    }

    /* Gather on thread 0 the coordinates of the multipliers */ 
    MPI_Gatherv (&dlm_lambda_to_desactivate.x[0], countalloc, MPI_DOUBLE, 
    	&alldatax[0], counts, disps, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gatherv (&dlm_lambda_to_desactivate.y[0], countalloc, MPI_DOUBLE, 
    	&alldatay[0], counts, disps, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#if dimension == 3 
    MPI_Gatherv (&dlm_lambda_to_desactivate.z[0], countalloc, MPI_DOUBLE, 
    	&alldataz[0], counts, disps, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    
    /* send the total number of multipliers to all threads */
    mpi_all_reduce (m, MPI_INT, MPI_MAX);

    /* allocate now alldatax, alldatay, alldataz on threads !=0 */
    if (pid() != 0) {
      alldatax = (double*) calloc (m, sizeof(double));
      alldatay = (double*) calloc (m, sizeof(double));
#if dimension == 3
      alldataz = (double*) calloc (m, sizeof(double));
#endif
    }

  
    /* Now broadcast alldatax,alldatay,alldataz from thread 0 
    to other threads */
    MPI_Bcast (alldatax, m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (alldatay, m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#if dimension == 3 
    MPI_Bcast (alldataz, m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    
/*     for (int i = 0; i < m; i++) { */
/* #if dimension == 2  */
/*       printf("thread %d: particle %d, coord of the multiplier %d to be 
	removed (%g,%g,%g)\n", pid(), k, i, alldatax[i], alldatay[i], 0.); */
/* #elif dimension ==3 */
/*       printf("thread %d: particle %d, coord of the multiplier %d to be 
	removed (%g,%g,%g)\n", pid(), k, i, alldatax[i], alldatay[i], 
	alldataz[i]); */
/* #endif */
/*     } */
    
    for (int i = 0; i < m; i++) {
#if dimension == 2 
      lpoint = locate(alldatax[i], alldatay[i]); 
#elif dimension == 3
      lpoint = locate(alldatax[i], alldatay[i], alldataz[i]); 
#endif
      if (lpoint.level > -1)  cache_append(&c, lpoint, 0);
    }

    free(alldatax);
    free(alldatay);
#if dimension == 3 
    free(alldataz);
#endif

#elif _MPI == 0

    for (int i = 0; i < countalloc; i++) {
#if dimension == 2 
      lpoint = locate(dlm_lambda_to_desactivate.x[i], 
      	dlm_lambda_to_desactivate.y[i]); 
#elif dimension == 3
      lpoint = locate(dlm_lambda_to_desactivate.x[i], 
      	dlm_lambda_to_desactivate.y[i], dlm_lambda_to_desactivate.z[i]); 
#endif
      if (lpoint.level > -1)  cache_append(&c, lpoint, 0);
    }

#endif   /* End if _MPI */

    /* Finally desactivate the multiplicators */
    int counter = 0;
    foreach_cache(c) {
      if (index_lambda.x[] > -1 && index_lambda.y[] > -1) {
	index_lambda.x[] = -1;
	index_lambda.y[] = -1;
	counter ++;
      }
    }
    
    /* printf("thread %d: removed %d multipliers \n", pid(), counter); */

    free(c.p);
    free_SolidBodyBoundary(&dlm_lambda_to_desactivate);
  } /* End loop NPARTICLES */

  boundary((scalar*) {index_lambda});
}




/** Tgs cells associated to a Lagrange multiplier on the boundary of the 
fictitious domain. */
//----------------------------------------------------------------------------
void reverse_fill_flagfield (particle * p, scalar f, vector index_lambda, 
	const int cacheflag, vector pshift) 
//----------------------------------------------------------------------------
{
  for (int k = 0; k < NPARTICLES; k++) 
  {  
    coord rel = {0., 0., 0.};
    coord relnl = {0., 0., 0.};
    int NCX, CX, weight_id;
    size_t goflag = 0;
    coord lambdacellpos = {0., 0., 0.};
    coord lambdapos = {0., 0., 0.};
    coord localcellpos = {0., 0., 0.};
 
    SolidBodyBoundary dlm_lambda = p[k].s;
    GeomParameter gcb = p[k].g;
    Cache * reduced_domain = &(p[k].reduced_domain);
    
    foreach_level(depth()) {      
      localcellpos.x = x;
      localcellpos.y = y;
#if dimension == 3 
      localcellpos.z = z;
#endif
    
      goflag = 0;       
      // Check if there is a Lagrange multiplier in the neigborhood of this 
      // cell;
      foreach_neighbor() {
	if (((int)index_lambda.x[] > -1)  && level == depth() 
		&& is_leaf(cell) && (p[k].pnum == (int)index_lambda.y[])) {
	  lambdacellpos.x = x; lambdacellpos.y = y; 
	  lambdapos.x = dlm_lambda.x[(int)index_lambda.x[]];
	  lambdapos.y = dlm_lambda.y[(int)index_lambda.x[]];
	
#if dimension == 3
	  lambdacellpos.z = z;
	  lambdapos.z = dlm_lambda.z[(int)index_lambda.x[]];
#endif
	  /* compute relative vector X_boundary - X_local = rel from the 
	  cell (containning the boundary) position to the boundary's 
	  (analytical) position */
	  compute_relative_vector (lambdacellpos, lambdapos, &rel);

	  /* reset dials integers */
	  NCX = 0; CX = 0; weight_id = 0; 

	  /* assign first normal (gives the quadrant) */ 
	  assign_dial (rel, &CX);

	  GeomParameter gcbdum;
	  gcbdum = gcb;
	
	  foreach_dimension()
	    gcbdum.center.x += pshift.x[];
	  
	  /* assign fictitious-boundary's normal (use boundary's position) */
	  assign_dial_fd_boundary (&p[k], lambdapos, gcbdum, Delta, &NCX);	  
	
	  /* compute relative vector neighbor to the local cell */
	  compute_relative_vector (localcellpos, lambdacellpos , &relnl);
	
	  /* Assign weight ids. */
	  assign_weight_id_quad_outward (NCX, CX, relnl, Delta, &weight_id, 
	  	&goflag);
	} // end if (lambda.x[] > -1)
      } // end foreach_neigboor loop

      /* If the cell belongs to the stencil of a Lagrange multiplier tag it 
      and create a reduced domain to optimize the stencil-traversal cost */
      if (goflag == 1) {
	f[] = 1;
	if (cacheflag == 1)
	  cache_append (reduced_domain, point, 0);
      }
    }
  
    boundary ({f, index_lambda.x, index_lambda.y});
  }// loop on particles id 
}




/** Initializes the particles and the scalar/vector fields needed to the 
method */
//----------------------------------------------------------------------------
void allocate_and_init_particles (particle * p, const int n, vector e, 
	scalar g, scalar h, vector pshift)
//----------------------------------------------------------------------------
{  
  foreach() {
    e.x[] = -1;// index_lambda.x => lambda indices
    e.y[] = -1;// particle structure/fd number;
    e.z[] = -1;// constrained cells (-1:not constrained)
    g[] = 0;   // flagfield
    h[] = 0;   // flagfield_mailleur

    /* Reset the shift vector for the periodic case */
    if (Period.x || Period.y || Period.z)
      foreach_dimension() {
	pshift.x[] = 0.;
      }
  }
  
  boundary ({g, h});
  boundary ((scalar *){e, pshift});


  for (int k = 0; k < n; k++) {
    Cache * c = NULL;

#if debugBD == 0
    GeomParameter gci = p[k].g;
    int m = 0;
    int lN = 0;
    
    switch( p[k].shape )
    {
      case SPHERE:
        compute_nboundary_Sphere( gci, &m );
        allocate_SolidBodyBoundary( &(p[k].s), m );
        create_FD_Boundary_Sphere( gci, &(p[k].s), m, pshift );
	break;
	  
      case CIRCULARCYLINDER2D:
        compute_nboundary_CircularCylinder2D( gci, &m );
        allocate_SolidBodyBoundary( &(p[k].s), m );
        create_FD_Boundary_CircularCylinder2D( gci, &(p[k].s), m, pshift );
	break;
	  
      case CUBE:
	compute_nboundary_Cube( &gci, &m, &lN );
        allocate_SolidBodyBoundary( &(p[k].s), m );
        create_FD_Boundary_Cube( &gci, &(p[k].s), m, lN, pshift );
	break;
	
      case TETRAHEDRON:
	compute_nboundary_Tetrahedron( &gci, &m, &lN );
        allocate_SolidBodyBoundary( &(p[k].s), m );
        create_FD_Boundary_Tetrahedron( &gci, &(p[k].s), m, lN, pshift );
	break;	
	  
      default:
        fprintf( stderr, "Unknown Rigid Body shape !!\n" );
    }    

    create_index_lambda_scalar ((p[k].s), e, k);
    c = &(p[k].reduced_domain);
    initialize_and_allocate_Cache(c);
#endif     
#if debugInterior == 0
    c = &(p[k].Interior);
    initialize_and_allocate_Cache(c);
#endif
  }
}




/** Writes headers in particle data output files */
//----------------------------------------------------------------------------
void writer_headers( FILE * pdata, FILE * sl ) 
//----------------------------------------------------------------------------
{
  // Comments: 
  // 1) we assume that initial time is always 0 when the simulation is 
  // not restarted
  // 2) we assume all hydro forces and torques are 0 at t=0
  // 3) we position the headers over the 1st line as a function of NSD, the 
  // number of significant digits after the decimal point

#if _MPI
  if ( pid() == 0 )
#endif
  {
    /* Write header for particles positions and velocities */
#if dimension == 2
    // Position and velocity file
    if ( NSDF > 9 )
      fprintf( pdata, "#time\t\t\tposition.x\t\tposition.y"
    	"\t\tU.x\t\t\tU.y\t\t\tw.z\n" );
    else
      fprintf( pdata, "#time\t\tposition.x\tposition.y"
    	"\tU.x\t\tU.y\t\tw.z\n" );

    // Hydro force and torque file
    if ( NSDF > 9 )
      fprintf( sl, "#time\t\t\tF.x\t\t\tF.y\t\t\tT.z\n" );
    else
      fprintf( sl, "#time\t\tF.x\t\tF.y\t\tT.z\n" );
    for (int i=0; i<3; ++i) fprintf( sl, "%.*e\t", NSDF, 0. ); 
    fprintf( sl, "%.*e\n", NSDF, 0. );   
#elif dimension == 3
    // Position and velocity file
    if ( NSDF > 9 )
      fprintf( pdata, "#time\t\t\tposition.x\t\tposition.y\t\tposition.z"
    	"\t\tU.x\t\t\tU.y\t\t\tU.z"
	"\t\t\tw.x\t\t\tw.y\t\t\tw.z\n" );
    else
      fprintf( pdata, "#time\t\tposition.x\tposition.y\tposition.z"
    	"\tU.x\t\tU.y\t\tU.z"
	"\t\tw.x\t\tw.y\t\tw.z\n" );

    // Hydro force and torque file
    if ( NSDF > 9 )
      fprintf( sl, "#time\t\t\tF.x\t\t\tF.y\t\t\tF.z\t\t\tT.x\t\t\tT.y"
      	"\t\t\tT.z\n" );
    else
      fprintf( sl, "#time\t\tF.x\t\tF.y\t\tF.z\t\tT.x\t\tT.y\t\tT.z\n" );
    for (int i=0; i<6; ++i) fprintf( sl, "%.*e\t", NSDF, 0. ); 
    fprintf( sl, "%.*e\n", NSDF, 0. );
#endif
    
    // Flush the buffers such that files are updated dynamically
    fflush( sl );
    fflush( pdata );
  }
}




/** Writes particle data in files */
//----------------------------------------------------------------------------
void particle_data( particle * p, const double t, const int i, FILE ** pdata ) 
//----------------------------------------------------------------------------
{  
  for (int k = 0; k < NPARTICLES; k++) 
  {
    GeomParameter* GCi = &(p[k].g);
#   if DLM_Moving_particle
#     if TRANSLATION 
        coord* U = &(p[k].U);
#     endif
#     if ROTATION
        coord* w =  &(p[k].w);
#     endif
#   endif
#   if dimension == 2
      if ( pid() == 0 ) 
      {
        fprintf( pdata[k], "%.*e\t", NSDF, t );
        fprintf( pdata[k], "%.*e\t", NSDF, (*GCi).center.x );
        fprintf( pdata[k], "%.*e\t", NSDF, (*GCi).center.y );         
        fprintf( pdata[k], "%.*e\t%.*e\t%.*e\n"     
#       if DLM_Moving_particle
#         if TRANSLATION
	    , NSDF, (*U).x, NSDF, (*U).y 
#         elif TRANSLATION == 0
	    , NSDF, 0., NSDF, 0.
#         endif   
#         if ROTATION
	    , NSDF, (*w).z
#         elif ROTATION == 0
	    , NSDF, 0.
#         endif
#       elif DLM_Moving_particle == 0
	  , NSDF, 0., NSDF, 0., NSDF, 0.
#       endif
        );      
        fflush(pdata[k]);
      }
#   elif dimension == 3
      if ( pid() == 0 ) 
      {
        fprintf( pdata[k], "%.*e\t", NSDF, t );
        fprintf( pdata[k], "%.*e\t", NSDF, (*GCi).center.x );
        fprintf( pdata[k], "%.*e\t", NSDF, (*GCi).center.y );      
        fprintf( pdata[k], "%.*e\t", NSDF, (*GCi).center.z );      
        fprintf( pdata[k], "%.*e\t%.*e\t%.*e\t%.*e\t%.*e\t%.*e\n"     
#       if DLM_Moving_particle
#         if TRANSLATION
	    , NSDF, (*U).x, NSDF, (*U).y, NSDF, (*U).z
#         elif TRANSLATION == 0
	    , NSDF, 0., NSDF, 0., NSDF, 0.
#         endif   
#         if ROTATION
	    , NSDF, (*w).x, NSDF, (*w).y, NSDF, (*w).z
#         elif ROTATION == 0
	    , NSDF, 0., NSDF, 0., NSDF, 0.
#         endif
#       elif DLM_Moving_particle == 0
	  , NSDF, 0., NSDF, 0., NSDF, 0., NSDF, 0., NSDF, 0., NSDF, 0.
#       endif
        );      
        fflush(pdata[k]);      	
      }
#   endif
  }
}




/** Compute hydrodynamic force & torque and write the values in files */
//----------------------------------------------------------------------------
coord sumLambda( particle* p, FILE** sl, const double t, const double dt, 
	scalar Flagfield, vector Lambda, vector Index_lambda, 
	const double rho_f, vector pshift ) 
//----------------------------------------------------------------------------
{
  /* Compute hydrodynamic force & torque and write to a file */
  /* Fh = <lambda,V>_P + (rho_f/rho_s)MdU/dt */
  /* Mh = <lambda,xi^GM>_P + (rho_f/rho_s).( Idom/dt + om x (I.om)) */

  /* The inertia tensor is: */
  /* Ip[0] = Ixx */
  /* Ip[1] = Iyy */
  /* Ip[2] = Izz */
  /* Ip[3] = Ixy */
  /* Ip[4] = Ixz */
  /* Ip[5] = Iyz */  

  coord lambdasumint;
  coord lambdasumboundary;
  coord lambdasum;
  coord crossLambdaSumInt;
  coord crossLambdaSumBoundary;
  coord crossLambdaSum; 
  Cache * Interior[NPARTICLES];
  Cache * Boundary[NPARTICLES];
  GeomParameter * gci;
  SolidBodyBoundary * sbm;
  
  /* Loop over all particles */
  for (int k = 0; k < NPARTICLES; k++) {

    foreach_dimension() {
      lambdasumint.x = 0.;
      lambdasumboundary.x = 0.;
      lambdasum.x = 0;
      
      crossLambdaSumInt.x = 0.;
      crossLambdaSumBoundary.x = 0.;
      crossLambdaSum.x = 0.;
    }

#if dimension ==2 
    lambdasumint.z = 0.;
    lambdasumboundary.z = 0.;
    lambdasum.z = 0;
      
    crossLambdaSumInt.z = 0.;
    crossLambdaSumBoundary.z = 0.;
    crossLambdaSum.z = 0.;
#endif

    /* Interior domain of the particle */
    Interior[k] = &(p[k].Interior);

    /* Boundary domain of the particle */
    Boundary[k] = &(p[k].reduced_domain);
    gci = &(p[k].g);
    sbm = &(p[k].s);
    coord lambdapos;
#if DLM_Moving_particle
    double rho_s = p[k].rho_s;
#if TRANSLATION
    double M = p[k].M;
    coord Unm1 = p[k].Unm1;
    coord U = p[k].U;
#endif
#if ROTATION
    coord wnm1 = p[k].wnm1;
    coord w = p[k].w;
    coord Iom, omIom, Idomdt;
#endif
#endif
    
    // Particule's interior multipliers
    foreach_cache((*Interior[k])) {
      if (Flagfield[] < 1 && k == Index_lambda.y[]) {

	/* Compute Fh = <lambda,V>_P */
	/* For moving particles the additional term +
	   (rho_f/rho_s).M.dU/dt is added after the mpi_calls below */
	
	foreach_dimension() {
	  lambdasumint.x += Lambda.x[];
	}

	/* Compute Mh = <lambda,xi^GM>_P */  
	/* For moving particles, the additional term +
	   (rho_f/rho_s).(I.dom/dt + om ^ (I.om)) is added after mpi
	   calls below */
	
	/* Modify temporerly the particle center position for periodic 
	boundary condition */
	foreach_dimension()
	  (*gci).center.x += pshift.x[];

	crossLambdaSumInt.x += (Lambda.z[]*(y - (*gci).center.y) 
		- Lambda.y[]*(z - (*gci).center.z));
	crossLambdaSumInt.y += (Lambda.x[]*(z - (*gci).center.z) 
		- Lambda.z[]*(x - (*gci).center.x));
	crossLambdaSumInt.z += (Lambda.y[]*(x - (*gci).center.x) 
		- Lambda.x[]*(y - (*gci).center.y));
	foreach_dimension()
	  (*gci).center.x -= pshift.x[];

      }
    }//end foreach_cache()

    // Particle's boundary multipliers
    foreach_cache((*Boundary[k])) {
      if ((Index_lambda.x[] > -1) && (p[k].pnum == (int)Index_lambda.y[])) {
	lambdapos.x = (*sbm).x[(int)Index_lambda.x[]];
	lambdapos.y = (*sbm).y[(int)Index_lambda.x[]];
	lambdapos.z = 0.;
#if dimension == 3
	lambdapos.z = (*sbm).z[(int)Index_lambda.x[]];
#endif

	/* Compute Fh = <lambda,V>_P */
	foreach_dimension() {
	  lambdasumboundary.x += Lambda.x[];
	}


	/* Modify temporerly the particle center position for periodic 
	boundary condition */
	foreach_dimension()
	  (*gci).center.x += pshift.x[];
	
	/* Compute Mh = <lambda,xi^GM>_P */
        crossLambdaSumBoundary.x += (Lambda.z[]*(lambdapos.y - (*gci).center.y)
		 - Lambda.y[]*(lambdapos.z - (*gci).center.z));
        crossLambdaSumBoundary.y += (Lambda.x[]*(lambdapos.z - (*gci).center.z)
		 - Lambda.z[]*(lambdapos.x - (*gci).center.x));
        crossLambdaSumBoundary.z += (Lambda.y[]*(lambdapos.x - (*gci).center.x)
		 - Lambda.x[]*(lambdapos.y - (*gci).center.y));

	foreach_dimension()
	  (*gci).center.x -= pshift.x[];
      }
    }
 
#if _MPI
    foreach_dimension() {
      mpi_all_reduce (lambdasumint.x, MPI_DOUBLE, MPI_SUM);
      mpi_all_reduce (lambdasumboundary.x, MPI_DOUBLE, MPI_SUM);
      mpi_all_reduce (crossLambdaSumInt.x, MPI_DOUBLE, MPI_SUM);
      mpi_all_reduce (crossLambdaSumBoundary.x, MPI_DOUBLE, MPI_SUM);
    }
#if dimension == 2
    mpi_all_reduce (crossLambdaSumInt.z, MPI_DOUBLE, MPI_SUM);
    mpi_all_reduce (crossLambdaSumBoundary.z, MPI_DOUBLE, MPI_SUM);
#endif
#endif

    /* The Force term (rho_f/rho_s)MdU/dt is added here for
       translating particles */
    
#if DLM_Moving_particle
#if TRANSLATION
    foreach_dimension()
      lambdasumint.x += (rho_f/rho_s)*M*(U.x - Unm1.x)/dt;
#endif
    
    /* The Torque term (rho_f/rho_s).(Idom/dt + om ^ (I.om)) is added
	 here for rotating particles */
    
#if ROTATION
    /* I.om term */
    Iom.x = p[k].Ip[0]*w.x + p[k].Ip[3]*w.y + p[k].Ip[4]*w.z;
    Iom.y = p[k].Ip[3]*w.x + p[k].Ip[1]*w.y + p[k].Ip[5]*w.z;
    Iom.z = p[k].Ip[4]*w.x + p[k].Ip[5]*w.y + p[k].Ip[2]*w.z;

    /* om^(I.om) term */
    omIom.x = w.y*Iom.z - w.z*Iom.y;
    omIom.y = w.z*Iom.x - w.x*Iom.z;
    omIom.z = w.x*Iom.y - w.y*Iom.x;

    /* Idom/dt term; */
    Idomdt.x = p[k].Ip[0]*(w.x - wnm1.x) + p[k].Ip[3]*(w.y - wnm1.y) 
    	+ p[k].Ip[4]*(w.z - wnm1.z);
    Idomdt.y = p[k].Ip[3]*(w.x - wnm1.x) + p[k].Ip[1]*(w.y - wnm1.y) 
    	+ p[k].Ip[5]*(w.z - wnm1.z);
    Idomdt.z = p[k].Ip[4]*(w.x - wnm1.x) + p[k].Ip[5]*(w.y - wnm1.y) 
    	+ p[k].Ip[2]*(w.z - wnm1.z);
	
    crossLambdaSumInt.x += (rho_f/rho_s)*((Idomdt.x)/dt + omIom.x);
    crossLambdaSumInt.y += (rho_f/rho_s)*((Idomdt.y)/dt + omIom.y);
    crossLambdaSumInt.z += (rho_f/rho_s)*((Idomdt.z)/dt + omIom.z);
#endif	  
#endif

    
    foreach_dimension() {
      lambdasum.x = lambdasumint.x + lambdasumboundary.x;
      crossLambdaSum.x = crossLambdaSumInt.x + crossLambdaSumBoundary.x;
    }
    

#if dimension == 2
    crossLambdaSum.z = crossLambdaSumInt.z + crossLambdaSumBoundary.z;
      
    if( pid() == 0 ) 
    {
      fprintf( sl[k], "%.*e\t", NSDF, t );
      fprintf( sl[k], "%.*e\t", NSDF, lambdasum.x );      
      fprintf( sl[k], "%.*e\t", NSDF, lambdasum.y );          
      fprintf( sl[k], "%.*e\n", NSDF, crossLambdaSum.z );      
      fflush( sl[k] );
    }
#elif dimension == 3
    if( pid() == 0 ) 
    {
      fprintf( sl[k], "%.*e\t", NSDF, t );
      fprintf( sl[k], "%.*e\t", NSDF, lambdasum.x );      
      fprintf( sl[k], "%.*e\t", NSDF, lambdasum.y );      
      fprintf( sl[k], "%.*e\t", NSDF, lambdasum.z );      
      fprintf( sl[k], "%.*e\t", NSDF, crossLambdaSum.x );      
      fprintf( sl[k], "%.*e\t", NSDF, crossLambdaSum.y );      
      fprintf( sl[k], "%.*e\n", NSDF, crossLambdaSum.z );      
      fflush( sl[k] );
    }
#endif
  }
  
  return lambdasum;
}




/** Initialize/open all DLMFD file pointers */
//----------------------------------------------------------------------------
void init_file_pointers( FILE** p, FILE** d, FILE** UzawaCV, FILE** CVT, 
	const size_t rflag ) 
//----------------------------------------------------------------------------
{
  char name[80] = "";
  char name2[80] = "";
  char suffix[80] = "";
  char buffer[80] = "";  
# if _MPI
    if ( pid() == 0 )
# endif
    {
      // Particle data
      for (int k = 0; k < NPARTICLES; k++) 
      {
        sprintf( suffix, "_%d.dat", k );

        strcpy( name, result_dir );
        strcat( name, "/" );
        strcat( name, result_particle_vp_rootfilename );
        strcat( name, suffix );
      
        strcpy( name2, result_dir );
        strcat( name2, "/" );
        strcat( name2, result_particle_hydroFaT_rootfilename );
        strcat( name2, suffix );      

        if ( !rflag ) 
        {
	  p[k] = fopen( name,  "w" ); 
	  d[k] = fopen( name2, "w" );
	
	  // Write headers in these files
	  writer_headers( p[k],  d[k] );
        }
        else 
        {
	  p[k] = fopen( name,  "a" );
	  d[k] = fopen( name2, "a" );
        }
      }
    
      // Uzawa convergence
      char converge_uzawa_filename_complete_name[80];
      strcpy( buffer, result_dir );
      strcat( buffer, "/" );
      strcat( buffer, converge_uzawa_filename );
      strcpy( converge_uzawa_filename_complete_name, buffer );
      if ( !rflag )
      {
        *UzawaCV = fopen( converge_uzawa_filename_complete_name, "w" ); 
        fprintf( *UzawaCV, "# Iter \t Uzawa Iter \t ||u-u_imposed||\n" );
      }
      else
        *UzawaCV = fopen( converge_uzawa_filename_complete_name, "a" );   
    
      // Cells, contrained cells and nb of multipliers 
      char dlmfd_cells_filename_complete_name[80]; 
      strcpy( buffer, result_dir );
      strcat( buffer, "/" );
      strcat( buffer, dlmfd_cells_filename );
      strcpy( dlmfd_cells_filename_complete_name, buffer );
      if ( !rflag )
      {
        *CVT = fopen ( dlmfd_cells_filename_complete_name, "w" ); 
        fprintf ( *CVT,"# Iter \t LagMult \t ConstrainedCells \t "
    		"TotalCells\n" );
      }
      else
        *CVT = fopen( dlmfd_cells_filename_complete_name, "a" );          
    }    
}




/** Close all DLMFD files */
//----------------------------------------------------------------------------
void close_file_pointers( FILE** p, FILE** d, FILE* UzawaCV, FILE* CVT ) 
//----------------------------------------------------------------------------
{ 
# if _MPI
    if ( pid() == 0 )
# endif
    {
      // Particle data
      for (int k = 0; k < NPARTICLES; k++) 
      {
        fclose( p[k] ); 
        fclose( d[k] );
      }
    
      // Uzawa convergence
      fclose( UzawaCV );  
    
      // Cells, contrained cells and nb of multipliers 
      fclose( CVT );               
    }    
}




/** Compute vorticity in 3D and store it in field omega */
//----------------------------------------------------------------------------
void vorticity_3D( const vector u, vector omega ) 
//----------------------------------------------------------------------------
{  
  foreach()
    foreach_dimension()
      omega.x[] = ( (u.z[0,1,0] - u.z[0,-1,0]) 
      	- (u.y[0,0,1] - u.y[0,0,-1]) )/2.*Delta;
  
  boundary((scalar *){omega});
}




/** Computes the flow rate on the right boundary and writes the value
in a file */
//----------------------------------------------------------------------------
void compute_flowrate_right( FILE * pf, const vector u, const int level ) 
//----------------------------------------------------------------------------
{
  coord velo = {0., 0., 0.};

  foreach_dimension()
    velo.x = 0.;

#if !adaptive 
  Cache intDomain = {0};
  Point lpoint;
  foreach() {
    if (x > (L0 - Delta)) { 
	lpoint = locate(x, y, z);
	cache_append(&intDomain, lpoint, 0);
      }
  }
		
  cache_shrink(&intDomain);
    
  foreach_cache(intDomain){
    foreach_dimension()
      velo.x += sq(Delta)*u.x[];
  }
  
  free(intDomain.p);
#endif

#if adaptive
  double hh = L0/pow(2,level);
  double xi = 0., yj = 0., zval = 0.;
  int ii = 0, jj = 0;
  
  foreach_level(level) {
    zval = L0 - 0.5*Delta + Z0;
  }

  /* printf("zz = %f, h = %f\n", zval, hh); */
  
  for (ii = 0; ii < pow(2,level); ii++) {
    xi = 0.5*hh + ii*hh + X0;
    for (jj = 0; jj < pow(2,level); jj++) {
      yj = 0.5*hh + jj*hh + Y0;
      coord uu = {0., 0., 0.};
      foreach_dimension() {
	uu.x = interpolate (u.x, xi, yj, zval);
	if (uu.x < nodata)
	  velo.x += sq(hh)*uu.x;
      }
      /* printf("thread %d, u = (%f,%f,%f)\n",pid(), 
      	interpolate(u.x, xi, yj, zz), interpolate(u.y, xi, yj, zz), 
	interpolate(u.z, xi, yj, zz)); */
    }
  }  
#endif

  /* printf("thread %d, velo = (%f,%f,%f)\n",pid(), velo.x, velo.y, velo.z); */


#if _MPI
  foreach_dimension ()
    mpi_all_reduce (velo.x, MPI_DOUBLE, MPI_SUM);
#endif
  
  fprintf(pf, "%f\t %20.18f\t %20.18f\t %20.18f\n", t, velo.x, velo.y, velo.z);
  fflush(pf);
}




/** Gets the pressure at a point and writes the value in a file */
//----------------------------------------------------------------------------
void pressure_law (scalar pres, FILE * ppf, const coord hv, const double t) 
//----------------------------------------------------------------------------
{
  double plaw = 0.;
  
  plaw = interpolate (pres, hv.x, hv.y, hv.z);

#if _MPI
  MPI_Barrier(MPI_COMM_WORLD);
  if (plaw == nodata)
    plaw = 0.;
  MPI_Barrier(MPI_COMM_WORLD);

  if (plaw != nodata)
    mpi_all_reduce(plaw, MPI_DOUBLE, MPI_SUM);
#endif
 
  if (pid() == 0) {
    fprintf(ppf, "%20.18f\t %20.18f\n", t, plaw);
    fflush(ppf);
  }
}



/** Inverts a 3 x 3 matrix and stores it in a double** */
//----------------------------------------------------------------------------
void inverse3by3matrix( double Matrix[3][3], double** inversedMatrix ) 
//----------------------------------------------------------------------------
{
  double block1 = Matrix[1][1] * Matrix[2][2] 
  	- Matrix[1][2] * Matrix[2][1];
  double block2 = Matrix[1][2] * Matrix[2][0] 
  	- Matrix[1][0] * Matrix[2][2];
  double block3 = Matrix[1][0] * Matrix[2][1] 
  	- Matrix[1][1] * Matrix[2][0];

  double determinant = Matrix[0][0] * block1 + Matrix[0][1] * block2 
  	+ Matrix[0][2] * block3;

  /* m[i][j], *(*(m + i) + j) */
  
  *(*(inversedMatrix + 0) + 0) = block1 / determinant;
  *(*(inversedMatrix + 0) + 1) = ( Matrix[0][2] * Matrix[2][1] 
  	- Matrix[0][1] * Matrix[2][2] ) / determinant;
  *(*(inversedMatrix + 0) + 2) = ( Matrix[0][1] * Matrix[1][2] 
  	- Matrix[0][2] * Matrix[1][1]) / determinant; 
  *(*(inversedMatrix + 1) + 0) = block2 / determinant;  
  *(*(inversedMatrix + 1) + 1) = ( Matrix[0][0] * Matrix[2][2] 
  	- Matrix[0][2] * Matrix[2][0] ) / determinant;
  *(*(inversedMatrix + 1) + 2) = ( Matrix[0][2] * Matrix[1][0] 
  	- Matrix[0][0] * Matrix[1][2] ) / determinant;
  *(*(inversedMatrix + 2) + 0) = block3 / determinant;
  *(*(inversedMatrix + 2) + 1) = ( Matrix[0][1] * Matrix[2][0] 
  	- Matrix[0][0] * Matrix[2][1]) / determinant;
  *(*(inversedMatrix + 2) + 2) = ( Matrix[0][0] * Matrix[1][1] 
  	- Matrix[0][1] * Matrix[1][0] ) / determinant ;
}




/** Inverts a 3 x 3 matrix and stores it in a double[][] */
//----------------------------------------------------------------------------
void inverse3by3matrix__( double Matrix[3][3], double inversedMatrix[3][3] ) 
//----------------------------------------------------------------------------
{
  double block1 = Matrix[1][1] * Matrix[2][2] 
  	- Matrix[1][2] * Matrix[2][1];
  double block2 = Matrix[1][2] * Matrix[2][0] 
  	- Matrix[1][0] * Matrix[2][2];
  double block3 = Matrix[1][0] * Matrix[2][1] 
  	- Matrix[1][1] * Matrix[2][0];

  double determinant = Matrix[0][0] * block1 + Matrix[0][1] * block2 
  	+ Matrix[0][2] * block3;
  
  inversedMatrix[0][0] = block1 / determinant;
  inversedMatrix[0][1] = ( Matrix[0][2] * Matrix[2][1] 
  	- Matrix[0][1] * Matrix[2][2]) / determinant;
  inversedMatrix[0][2] = ( Matrix[0][1] * Matrix[1][2] 
  	- Matrix[0][2] * Matrix[1][1] ) / determinant; 
  inversedMatrix[1][0] = block2 / determinant;
  inversedMatrix[1][1] = ( Matrix[0][0] * Matrix[2][2] 
  	- Matrix[0][2] * Matrix[2][0] ) / determinant;
  inversedMatrix[1][2] = ( Matrix[0][2] * Matrix[1][0] 
  	- Matrix[0][0] * Matrix[1][2] ) / determinant;
  inversedMatrix[2][0] = block3 / determinant;
  inversedMatrix[2][1] = ( Matrix[0][1] * Matrix[2][0] 
  	- Matrix[0][0] * Matrix[2][1] ) / determinant;
  inversedMatrix[2][2] = ( Matrix[0][0] * Matrix[1][1] 
  	- Matrix[0][1] * Matrix[1][0]) / determinant ;
}




#if DLM_Moving_particle
/** Computes the inverse of the moment of inertia matrix of a rigid body */
//----------------------------------------------------------------------------
void compute_inv_inertia( particle *p )
//----------------------------------------------------------------------------
{
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

  // Transfer Ip to a 3x3 matrix
  double Imat[3][3]; 
  Imat[0][0] = p->Ip[0];
  Imat[1][1] = p->Ip[1];
  Imat[2][2] = p->Ip[2];

  Imat[0][1] = p->Ip[3];
  Imat[0][2] = p->Ip[4];
  Imat[1][2] = p->Ip[5];
    
  Imat[1][0] = Imat[0][1]; 
  Imat[2][0] = Imat[0][2];
  Imat[2][1] = Imat[1][2];
  
  // Invert Imat and copy the result in p->Ip_inv
  // Ip_inv is a 2D 3 by 3 array
  inverse3by3matrix__( Imat, p->Ip_inv );  
}
#endif




/** Computes and returns the total number of cells of the grid */
//----------------------------------------------------------------------------
int totalcells() 
//----------------------------------------------------------------------------
{
  int t = 0;
  
  foreach() t++;
  
# if _MPI
    mpi_all_reduce( t, MPI_INTEGER, MPI_SUM );
# endif
  
  return t;
}




/** Computes and returns the total number of cells related to Distributed
Lagrange multiplier points for all rigid bodies */
//----------------------------------------------------------------------------
int total_dlmfd_cells( particle* allpart, const int np ) 
//----------------------------------------------------------------------------
{
  int apts = 0;
  
  for (int k = 0; k < np; k++) 
  {
#   if debugInterior == 0
      apts += allpart[k].Interior.n;
#   endif
#   if debugBD == 0
      apts += allpart[k].reduced_domain.n;
#   endif
  }
  
# if _MPI
    mpi_all_reduce( apts, MPI_INTEGER, MPI_SUM );
# endif

  return apts;
}




/** Computes and returns the total number of Ditributed Lagrange multiplier 
points for all rigid bodies */
//----------------------------------------------------------------------------
int total_dlmfd_multipliers( particle* allpart, const int np ) 
//----------------------------------------------------------------------------
{
  int apts = 0;
  
# if debugInterior == 0
    for (int k = 0; k < np; k++) 
      apts += allpart[k].Interior.n;
  
#   if _MPI
      mpi_all_reduce (apts, MPI_INTEGER, MPI_SUM);
#   endif
# endif
  
# if debugBD == 0
    for (int k = 0; k < np; k++) 
    {
      SolidBodyBoundary * bla;
      bla = &(allpart[k].s);
      apts += bla->m;
    }
# endif
  
  return apts;
}




/** Return the smallest grid size */
//----------------------------------------------------------------------------
double getDelta() 
//----------------------------------------------------------------------------
{
  int ii = 0;
  double dd = 0.;
  
  foreach_level(depth()) 
  {
    ii++;
    dd = Delta;
    if ( ii > 0 ) break;
  }

# if _MPI
    mpi_all_reduce( dd, MPI_DOUBLE, MPI_MIN );
# endif
  
  return dd;
}




/** Save the time and time step in a file */
//----------------------------------------------------------------------------
void save_t_dt_restart( char* dirname, double time, double deltat )
//----------------------------------------------------------------------------
{
  char dump_name[80] = "";
  strcpy( dump_name, dirname );
  strcat( dump_name, "/t_dt_restart.res" );
  FILE* ft = fopen( dump_name, "w" );
  fprintf ( ft, "%10e %10e", time, deltat );
  fclose( ft );  
}




/** Read the restart time and time step from a file */
//----------------------------------------------------------------------------
void read_t_restart( char* dirname, double* time, double* deltat )
//----------------------------------------------------------------------------
{
  char dump_name[80] = "";
  strcpy( dump_name, dirname );
  strcat( dump_name, "/t_dt_restart.res" );
  FILE* ft = fopen( dump_name, "r" );
  fscanf ( ft, "%lf %lf", time, deltat );
  fclose( ft );  
}
