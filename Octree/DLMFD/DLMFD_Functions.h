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
  TETRAHEDRON,
  OCTAHEDRON,
  DODECAHEDRON,
  ICOSAHEDRON
};


 

/** Structure for the coordinates of a rigid body boundary. */
typedef struct {
  double* x;
  double* y;
  double* z;
  int m;
  int nm;
} SolidBodyBoundary;




/** Additional geometric parameters for polygons/polyhedrons */
typedef struct {
  int allPoints, allFaces;
  double** cornersCoord;
  long int** cornersIndex;
  long int* numPointsOnFaces;
} PolyGeomParameter;





/** Rigid body geometric parameters */
typedef struct {
  coord center;
  coord* perclonecenters;  
  double radius;
  int ncorners;
  int nperclones;
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
  char tag[3];
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




/** Force synchronization of a field in parallel by setting the dirty flag
of the fields to true */
//----------------------------------------------------------------------------
trace void synchronize( scalar* list )
//----------------------------------------------------------------------------
{
  for (scalar s in list)
    s.dirty = true;
  boundary(list);
}




# include "CircularCylinder2D.h"
# include "Sphere.h"
# include "Cube.h"
# include "Tetrahedron.h"
# include "Octahedron.h"
# include "Dodecahedron.h"
# include "Icosahedron.h"


/** Allocates memory for m points in the SolidBodyBoundary structure. */
//----------------------------------------------------------------------------
void allocate_SolidBodyBoundary( SolidBodyBoundary* sbm, const int m ) 
//----------------------------------------------------------------------------
{
  sbm->x = (double*) calloc( m, sizeof(double) ); 
  sbm->y = (double*) calloc( m, sizeof(double) );
# if dimension == 3  
    sbm->z = (double*) calloc( m, sizeof(double) );
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
  sbm->x = (double*) realloc( sbm->x, m * sizeof(double) ); 
  sbm->y = (double*) realloc( sbm->y, m * sizeof(double) );  
# if dimension == 3 
    sbm->z = (double*) realloc( sbm->z, m * sizeof(double) );
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

  sbm->m = 0;  
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
void free_particles( particle* allparticles, const size_t n ) 
//----------------------------------------------------------------------------
{
  for (size_t k=0;k<n;k++) 
  {        
    // Free the boundary point coordinate arrays
    SolidBodyBoundary* sbm = &(allparticles[k].s);
    free_SolidBodyBoundary( sbm );
    
    // Free the periodic clones position vector
    if ( allparticles[k].g.nperclones ) 
      free( allparticles[k].g.perclonecenters ); 
    
    // Free the caches 
    Cache* c = &(allparticles[k].Interior);
    free( c->p );
    c->p = NULL;
    c = &(allparticles[k].reduced_domain);
    free( c->p );
    c->p = NULL;    
    
    // Free the toy granular solver parameter structure
#   if DLM_Moving_particle
      if ( allparticles[k].toygsp )
      {
        free( allparticles[k].toygsp );
        allparticles[k].toygsp = NULL;
      }     
#   endif

    // Free the additional geometric features of the particle
    switch ( allparticles[k].shape )
    {
      case SPHERE:
        free_Sphere( &(allparticles[k].g) );
	break;
	  
      case CIRCULARCYLINDER2D:
        free_CircularCylinder2D( &(allparticles[k].g) );
	break;
	  
      case CUBE:
        free_Cube( &(allparticles[k].g) );
	break;

      case TETRAHEDRON:
        free_Tetrahedron( &(allparticles[k].g) );
	break;
	
      case OCTAHEDRON:
	free_Octahedron( &(allparticles[k].g) );
	break;
	
      case ICOSAHEDRON:
	free_Icosahedron( &(allparticles[k].g) );
	break;

      case DODECAHEDRON:
	free_Dodecahedron( &(allparticles[k].g) );
	break;	
	  
      default:
        fprintf( stderr,"Unknown Rigid Body shape !!\n" );
    }                               
  }
}




/** Prints data of a particle */
//----------------------------------------------------------------------------
void print_particle( particle const* p, char const* poshift )
//----------------------------------------------------------------------------
{
  if ( pid() == 0 ) 
  {
    printf( "%sNumber = %lu\n", poshift, p->pnum ); 
    printf( "%sTag = %s\n", poshift, p->tag ); 
    printf( "%sShape = ", poshift );
    switch ( p->shape )
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
	
      case OCTAHEDRON:
        printf( "OCTAHEDRON" );
	break;
	
      case ICOSAHEDRON:
        printf( "ICOSAHEDRON" );
	break;

      case DODECAHEDRON:
        printf( "DODECAHEDRON" );
	break;	
	  
      default:
        fprintf( stderr,"Unknown Rigid Body shape !!\n" );
    }
    printf( "\n" );
    printf( "%sCenter of mass position = %e %e", poshift, p->g.center.x, 
    	p->g.center.y );
#   if dimension == 3
      printf( " %e", p->g.center.z );
#   endif
    printf( "\n" );
    if ( p->g.nperclones )
    {
      printf( "%sNumber of periodic clones = %d\n", poshift, p->g.nperclones );
      printf( "%sPeriodic clones center of mass positions\n", poshift );
      for (int j=0; j < p->g.nperclones; j++)
      {
        printf( "%s   %e %e", poshift, p->g.perclonecenters[j].x, 
		 p->g.perclonecenters[j].y );
#       if dimension == 3
          printf( " %e", p->g.perclonecenters[j].z );
#       endif
        printf( "\n" );          
      }      
    }
    printf( "%sRadius = %e\n", poshift, p->g.radius );  
    printf( "%sMass = %e\n", poshift, p->M ); 
    printf( "%sVolume = %e\n", poshift, p->Vp );     
    printf( "%sDensity = %e\n", poshift, p->rho_s ); 
#   if DLM_Moving_particle
#     if dimension == 3
        printf( "%sInertia tensor\n", poshift );
        printf( "%s   Ixx = %e\n", poshift, p->Ip[0] );
        printf( "%s   Iyy = %e\n", poshift, p->Ip[1] );	  
        printf( "%s   Izz = %e\n", poshift, p->Ip[2] );	  
        printf( "%s   Ixy = %e\n", poshift, p->Ip[3] );	  
        printf( "%s   Ixz = %e\n", poshift, p->Ip[4] );	  
        printf( "%s   Iyz = %e\n", poshift, p->Ip[5] );
#     else
        printf( "%s   Inertia tensor component Izz = %e\n", poshift, 
		p->Ip[2] );
#     endif             
#     if TRANSLATION
        printf( "%sTranslational velocity = %e %e", poshift, p->U.x, p->U.y );
#       if dimension == 3
          printf( " %e", p->U.z );
#       endif
        printf( "\n" );
#     endif
#     if ROTATION
        printf( "%sAngular velocity = ", poshift );
#       if dimension == 3
          printf( "%e %e", p->w.x, p->w.y );
#       endif
        printf( " %e\n", p->w.z );
#     endif
#   else
      printf( "%sImposed translational velocity = %e %e", 
    	poshift, p->imposedU.x, p->imposedU.y );
#     if dimension == 3
        printf( " %e", p->imposedU.z );
#     endif
      printf( "\n" );
      printf( "%sImposed angular velocity = ", poshift );
#     if dimension == 3
        printf( "%e %e", p->imposedw.x, p->imposedw.y );
#     endif
      printf( " %e\n", p->imposedw.z );    
#   endif
  }
  int intpts = p->Interior.n;
  int bdpts = p->reduced_domain.n;  
# if _MPI
    mpi_all_reduce( intpts, MPI_INT, MPI_SUM );
    mpi_all_reduce( bdpts, MPI_INT, MPI_SUM );
# endif
  if ( pid() == 0 )
  {   
    printf( "%sNumber of interior DLM points = %d\n", poshift, intpts ); 
    printf( "%sNumber of boundary DLM points = %d\n", poshift, p->s.m ); 
    printf( "%sNumber of cells in boundary DLM point stencils = %d\n", poshift, 
    	bdpts );
  }
}  




/** Prints all particle data */
//----------------------------------------------------------------------------
void print_all_particles( particle const* allparticles, const size_t n,
	char const* oshift )
//----------------------------------------------------------------------------
{
  if ( pid() == 0 ) printf( "%sTotal number of particles = %lu\n", oshift, 
  	NPARTICLES );
  char poshift[20]="   ";
  strcat( poshift, oshift );
  for (size_t k=0;k<n;k++)
  { 
    if ( pid() == 0 ) printf( "%sParticle %lu\n", oshift, k );    
    print_particle( &(allparticles[k]), &poshift[0] );
  }
} 




/** Tags the cells that contain a DLMFD boundary point, i.e. assign the point 
number to the x component of the index field and the particle number to the y 
component of the index field. If there is no DLMFD boundary point, index.x is
set to -1 */
//----------------------------------------------------------------------------
void fill_DLM_Index( const SolidBodyBoundary dlm_bd, vector Index, 
	const size_t kk ) 
//----------------------------------------------------------------------------
{  
  Point lpoint;
  int i;
  Cache* fdlocal;
   
  fdlocal = (Cache *){calloc(dlm_bd.m, sizeof(Cache))};
  
  for (i=0;i<dlm_bd.m;i++) 
  {
    # if dimension == 2 
        lpoint = locate( dlm_bd.x[i], dlm_bd.y[i] );
    # elif dimension == 3
        lpoint = locate( dlm_bd.x[i], dlm_bd.y[i], dlm_bd.z[i] );
    # endif    

    if ( (lpoint.level) == depth() ) 
    {
      /* Create a cache for each fictitious domain's boundary point */
      cache_append( &fdlocal[i], lpoint, 0 );
	
      foreach_cache(fdlocal[i]) 
      {
	/* Tag cell only if it was not tagged by another particle */
	if ( Index.x[] < 0 )
	{
	  Index.x[] = i;
	  Index.y[] = kk;
	  if ( level != depth() ) 
	  {
	    printf( "On thread %d, point dlmfd %d at (%f, %f, %f) is in a"
		" cell that has not the maximum level of refinnement %d, it "
		"is on level %d \n", pid(), i, x, y, z, depth(), level );
	  }
	}
      }
    }
  }
   
  for (i=0;i<dlm_bd.m;i++)
    free(fdlocal[i].p);
  
  free (fdlocal);
  
  synchronize((scalar*) {Index});
}





#include "DLMFD_Weight_functions.h"

/** Returns the weight associated to a cell if it belongs to the 3^dim stencil
associated to a Lagrange multiplier, otherwise return 0. The cell that contains 
the Lagrange multiplier can be on another process (i.e be a ghost cell for the 
current process). This function has to be embedded within a double 
foreach_cache() and foreach_neighbor() loop. The foreach_cache() loops on all 
non-ghost cells on this process that have been tagged to belong to a Lagrange 
multiplier stencil and the foreach_neighbor() loops on the neighbors of the 
non-ghost cells in a 5^dim stencil. */
//----------------------------------------------------------------------------
double reversed_weight( particle* pp, const coord weightcellpos, 
	const coord lambdacellpos, const coord lambdapos, const double delta ) 
//----------------------------------------------------------------------------
{
  coord rel = {0, 0, 0};
  coord relnl = {0, 0, 0};
  int NCX = 0, CX = 0, weight_id = 0;
  size_t goflag = 0;
  double weight = 0.;
  GeomParameter gcb = pp->g; 
  
  /* Compute relative vector from the cell (containning the boundary) position 
  to the boundary's (analytical) position */
  foreach_dimension() rel.x = lambdapos.x - lambdacellpos.x;

  /* In a periodic case, the position of the Lagrange multipliers on the
  boundary are those of the master rigid body, i.e., they are not shifted for
  each periodic clone. Therefore we may have a case where the cell and the
  position of the Lagrange multiplier are on opposite side of the domain. The
  solution is to subtract the periodic domain length L0 when this happens */
  foreach_dimension()
  {
    if ( rel.x > L0/2. ) rel.x = rel.x - L0;
    else if ( rel.x < -L0/2. ) rel.x = rel.x + L0;
  }  

  /* Reset dial integers */
  NCX = 0; CX = 0; weight_id = 0; goflag = 0;

  /* Assign quadrant number CX defining relative position of the Lagrange point 
  with respect to the center of the cell it belongs to */ 
  assign_dial( rel, &CX );
  
  GeomParameter gcbdum;
  gcbdum = gcb;

  /* Assign quadrant number NCX given by the direction of the normal over 
  the geometric boundary of the rigid body */ 
  assign_dial_fd_boundary( pp, lambdapos, gcbdum, delta, &NCX );
  
  /* Compute relative vector from the cell to the cell that contains the 
  Lagrange multiplier */
  foreach_dimension() relnl.x = lambdacellpos.x - weightcellpos.x;  

  /* Assign weight id if this cell is flagged ( goflag = 1 ) */
  assign_weight_id_quad_outward( NCX, CX, relnl, delta, &weight_id, &goflag );
    
  /* If the cell belongs to the 3^dim stencil of this Lagrange multiplier, then 
  compute the weight of this cell in the quadratic interpolation
  Otherwise return 0 */
  if ( goflag == 1 ) 
    weight = compute_weight_Quad( weight_id, lambdapos, lambdacellpos, 
    	NCX, CX, delta );
  
  return ( weight );
}




//----------------------------------------------------------------------------
void remove_too_close_multipliers( particle* p, vector DLM_Index) 
//----------------------------------------------------------------------------
{   
  for (size_t k = 0; k < NPARTICLES; k++) {

    SolidBodyBoundary dlm_lambda_to_desactivate;
    int allocated = 1;
    allocate_SolidBodyBoundary (&dlm_lambda_to_desactivate, allocated*BSIZE);
    int countalloc = 0;
    int other_part;
    int is_in_other_particle;
    foreach_level(depth()) {
      
      int direct_neigh = 0;  is_in_other_particle = 0; other_part = -1;
      if (((int)DLM_Index.x[] > -1)  && level == depth() && is_leaf(cell) 
      	&& (p[k].pnum == (int)DLM_Index.y[])) {

	/* check here if this multiplier is not in another particle's
	   domain, if yes desactive it */
	
	for (size_t l = 0; l < NPARTICLES; l++) {
	  if (l != p[k].pnum) {
      	    particle * other_particle = &(p[l]);
      	    /* printf("thread %d, this is particle %zu checkin on particle 
	    %zu iteration %d\n", pid(), p[k].pnum, other_particle->pnum, l); */

      	    /* Check particle's type */
      	    if ( other_particle->shape == CUBE ) {
//       	      compute_principal_vectors_Cube( other_particle );
//       	      coord u = other_particle->g.pgp->u1;
//       	      coord v = other_particle->g.pgp->v1;
//       	      coord w = other_particle->g.pgp->w1;
//       	      coord mins = other_particle->g.pgp->mins;
//       	      coord maxs = other_particle->g.pgp->maxs;
// 
//       	      /* current cell's position */
//       	      coord checkpt = {x, y, z};
//       	      if ( is_in_Cube( &u, &v, &w, &mins, &maxs, &checkpt ) ) {
//       		is_in_other_particle = 1; other_part =  other_particle->pnum;
// 		break;
//       	      }
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
	  DLM_Index.x[] = -1; DLM_Index.y[] = other_part;
	}
	
	/* Check if two (or more) Lagrange multipliers (from different
	   particles) are in the same neighborhoods (i.e a 5x5 stencil
	   in 2D), if yes desactivate them. Otherwise the cells of
	   their stencils may be located in each other's domain */
	direct_neigh = 0;
	foreach_neighbor() {
	  if (((int)DLM_Index.x[] > -1)  && level == depth() 
	  	&& is_leaf(cell) && (p[k].pnum != (int)DLM_Index.y[])) {

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
	DLM_Index.x[] = -1; DLM_Index.y[] = -1;
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
      if (DLM_Index.x[] > -1 && DLM_Index.y[] > -1) {
	DLM_Index.x[] = -1;
	DLM_Index.y[] = -1;
	counter ++;
      }
    }
    
    /* printf("thread %d: removed %d multipliers \n", pid(), counter); */

    free(c.p);
    free_SolidBodyBoundary(&dlm_lambda_to_desactivate);
  } /* End loop NPARTICLES */

  synchronize((scalar*) {DLM_Index});
}




/** Tag cells that belong to a 3^dim stencil associated to a Lagrange multiplier
point of the rigid body boundary. */
//----------------------------------------------------------------------------
void reverse_fill_DLM_Flag( particle* allparticles, const size_t n, scalar Flag, 
	vector Index, const int cacheflag ) 
//----------------------------------------------------------------------------
{
  for (size_t k = 0; k < n; k++) 
  {  
    coord rel = {0., 0., 0.};
    coord relnl = {0., 0., 0.};
    int NCX, CX, weight_id;
    size_t goflag = 0;
    coord lambdacellpos = {0., 0., 0.};
    coord lambdapos = {0., 0., 0.};
    coord localcellpos = {0., 0., 0.};
 
    SolidBodyBoundary dlm_lambda = allparticles[k].s;
    GeomParameter gcb = allparticles[k].g;
    Cache* reduced_domain = &(allparticles[k].reduced_domain);
        
    foreach_level(depth()) 
    {      
      localcellpos.x = x;
      localcellpos.y = y;
#     if dimension == 3 
        localcellpos.z = z;
#     endif

      // IMPORTANT REMARK: we initialize the goflag variable to 0 outside 
      // the foreach_neighbor() because a cell might belong to 2 or more 
      // different Lagrange multiplier stencils but we want to tag it and add 
      // it to the reduced domain cache once only. Later, when we compute the 
      // contribution of this cell to the different stencils in 
      // DLM_Uzawa_velocity with the function reversed_weight, the contribution
      // of this cell to each stencil will be properly computed      
      goflag = 0; 
             
      // Check if there is a Lagrange multiplier in the neigborhood of this 
      // cell
      foreach_neighbor() 
      {
	if ( (int)Index.x[] > -1 && level == depth() 
		&& is_leaf(cell) && allparticles[k].pnum == (int)Index.y[] ) 
	{
	  lambdacellpos.x = x; 
	  lambdacellpos.y = y; 
	  lambdapos.x = dlm_lambda.x[(int)Index.x[]];
	  lambdapos.y = dlm_lambda.y[(int)Index.x[]];	
#         if dimension == 3
	    lambdacellpos.z = z;
	    lambdapos.z = dlm_lambda.z[(int)Index.x[]];
#         endif

          /* Compute relative vector from the cell (containning the boundary) 
	  position to the boundary's (analytical) position */
          foreach_dimension() rel.x = lambdapos.x - lambdacellpos.x;

          /* In a periodic case, the position of the Lagrange multipliers on 
	  the boundary are those of the master rigid body, i.e., they are not 
	  shifted for each periodic clone. Therefore we may have a case where 
	  the cell and the position of the Lagrange multiplier are on opposite 
	  side of the domain. The solution is to subtract the periodic domain 
	  length L0 when this happens */
          foreach_dimension()
          {
            if ( rel.x > L0/2. ) rel.x = rel.x - L0;
            else if ( rel.x < -L0/2. ) rel.x = rel.x + L0;
          }
	  
          /* Reset dial integers */
	  NCX = 0; CX = 0; weight_id = 0;

          /* Assign quadrant number CX defining relative position of the 
	  Lagrange point with respect to the center of the cell it belongs to */ 
          assign_dial( rel, &CX );

	  GeomParameter gcbdum;
	  gcbdum = gcb;

          /* Assign quadrant number NCX given by the direction of the normal 
	  over the geometric boundary of the rigid body */ 
          assign_dial_fd_boundary( &allparticles[k], lambdapos, gcbdum, Delta, 
	  	&NCX );

          /* Compute relative vector from the cell to the cell that contains 
	  the Lagrange multiplier */
          foreach_dimension() relnl.x = lambdacellpos.x - localcellpos.x;
	
	  /* Assign weight id if this cell is flagged ( goflag = 1 ) */
	  assign_weight_id_quad_outward( NCX, CX, relnl, Delta, &weight_id, 
	  	&goflag );
	} // end if (lambda.x[] > -1)
      } // end foreach_neigboor loop

      /* If the cell belongs to at least one stencil of a Lagrange multiplier 
      tag it and add it the reduced domain of this particle to optimize the 
      stencil-traversal cost */
      if ( goflag == 1 ) 
      {
	Flag[] = 1;
	if ( cacheflag == 1 )
	  cache_append( reduced_domain, point, 0 );
      }
    }
  
    synchronize({Flag});    
  } // end loop on particles id 
}




/** Creates DLM/FD boundary points of a given particle. Sets the
PeriodicRefCenter field only if setPeriodicRefCenter is true */
//----------------------------------------------------------------------------
void create_boundary_points( particle* p, vector* pPeriodicRefCenter,
	const bool setPeriodicRefCenter )
//----------------------------------------------------------------------------
{  
  GeomParameter gci = p->g;
  int m = 0;
  int lN = 0;
    
  switch( p->shape )
  {
    case SPHERE:
      compute_nboundary_Sphere( gci, &m );
      allocate_SolidBodyBoundary( &(p->s), m );
      create_FD_Boundary_Sphere( gci, &(p->s), m, pPeriodicRefCenter, 
      		setPeriodicRefCenter );
      break;
	  
    case CIRCULARCYLINDER2D:
      compute_nboundary_CircularCylinder2D( gci, &m );
      allocate_SolidBodyBoundary( &(p->s), m );
      create_FD_Boundary_CircularCylinder2D( gci, &(p->s), m, 
      		pPeriodicRefCenter, setPeriodicRefCenter );
      break;
	  
    case CUBE:
      compute_nboundary_Cube( &gci, &m, &lN );
      allocate_SolidBodyBoundary( &(p->s), m );
      create_FD_Boundary_Cube( &gci, &(p->s), m, lN, pPeriodicRefCenter, 
		setPeriodicRefCenter );
      break;
	
    case TETRAHEDRON:
      compute_nboundary_Tetrahedron( &gci, &m, &lN );
      allocate_SolidBodyBoundary( &(p->s), m );
      create_FD_Boundary_Tetrahedron( &gci, &(p->s), m, lN, pPeriodicRefCenter, 
		setPeriodicRefCenter );
      break;
	
    case OCTAHEDRON:
      compute_nboundary_Octahedron( &gci, &m, &lN );
      allocate_SolidBodyBoundary( &(p->s), m );
      create_FD_Boundary_Octahedron( &gci, &(p->s), m, lN, pPeriodicRefCenter, 
		setPeriodicRefCenter );
      break;

    case ICOSAHEDRON:
      compute_nboundary_Icosahedron( &gci, &m, &lN );
      allocate_SolidBodyBoundary( &(p->s), m );
      create_FD_Boundary_Icosahedron( &gci, &(p->s), m, lN, pPeriodicRefCenter, 
		setPeriodicRefCenter );
      break;	

    case DODECAHEDRON:     
      compute_nboundary_Dodecahedron( &gci, &m, &lN );	
      allocate_SolidBodyBoundary( &(p->s), m );
      create_FD_Boundary_Dodecahedron( &gci, &(p->s), m, lN, pPeriodicRefCenter, 
		setPeriodicRefCenter );
      break;	
	  
    default:
      fprintf( stderr, "Unknown Rigid Body shape !!\n" );
  }
}




/** Initializes the particles and the scalar/vector fields needed to the 
method */
//----------------------------------------------------------------------------
void allocate_and_init_particles( particle* allparticles, const size_t n, 
	vector Index, scalar Flag, scalar FlagMesh, vector PeriodicRefCenter )
//----------------------------------------------------------------------------
{  
  /* Initialize fields */
  foreach() 
  {
    Index.x[] = -1;   // DLM/FD boundary point index
    Index.y[] = -1;   // particle number
    Index.z[] = -1;   // cell is constrained ? (-1:not constrained)
    Flag[] = 0;       // Flag
    FlagMesh[] = 0;   // DLM_FlagMesh
    if ( Period.x || Period.y || Period.z )
      foreach_dimension() PeriodicRefCenter.x[] = 0.;
  }
  
  synchronize({Flag, FlagMesh});
  synchronize((scalar *){Index, PeriodicRefCenter});


  for (size_t k = 0; k < n; k++) 
  {
    Cache* c = NULL;

#   if debugBD == 0
      create_boundary_points( &(allparticles[k]), &PeriodicRefCenter, true );
      fill_DLM_Index( (allparticles[k].s), Index, k );
      c = &(allparticles[k].reduced_domain);
      initialize_and_allocate_Cache( c );
#   endif     
#   if debugInterior == 0
      c = &(allparticles[k].Interior);
      initialize_and_allocate_Cache( c );
#   endif
  }
  
  synchronize((scalar*) {Index});
}




/** Creates boundary points of all particles but does not set the 
PeriodicRefCenter field */
//----------------------------------------------------------------------------
void create_particles_boundary_points( particle* allparticles, const size_t n )
//----------------------------------------------------------------------------
{  
# if debugBD == 0
    for (size_t k = 0; k < n; k++) 
      create_boundary_points( &(allparticles[k]), NULL, false );    
# endif
}




/** Frees boundary points of all particles  */
//----------------------------------------------------------------------------
void free_particles_boundary_points( particle* allparticles, const size_t n )
//----------------------------------------------------------------------------
{  
# if debugBD == 0
    for (size_t k = 0; k < n; k++) 
      free_SolidBodyBoundary( &(allparticles[k].s) );   
# endif
}




/** Writes headers in particle data output files */
//----------------------------------------------------------------------------
void writer_headers( FILE* pdata, FILE* sl ) 
//----------------------------------------------------------------------------
{
  // Comments: 
  // 1) we assume that initial time is always 0 when the simulation is 
  // not restarted
  // 2) we assume all hydro forces and torques are 0 at t=0
  // 3) we position the headers over the 1st line as a function of NSD, the 
  // number of significant digits after the decimal point

# if _MPI
    if ( pid() == 0 )
# endif
    {
      /* Write header for particles positions and velocities */
#     if dimension == 2
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
#     elif dimension == 3
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
#     endif
    
      // Flush the buffers such that files are updated dynamically
      fflush( sl );
      fflush( pdata );
    }
}




/** Writes particles data in files */
//----------------------------------------------------------------------------
void particle_data( particle* allparticles, const size_t n, const double t, 
	const int i, FILE** pdata ) 
//----------------------------------------------------------------------------
{  
  for (size_t k = 0; k < n; k++) 
  {
    GeomParameter* GCi = &(allparticles[k].g);
#   if DLM_Moving_particle
#     if TRANSLATION 
        coord* U = &(allparticles[k].U);
#     endif
#     if ROTATION
        coord* w =  &(allparticles[k].w);
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
void sumLambda( particle* allparticles, const size_t n, FILE** sl, 
	const double t, const double dt, 
	scalar Flag, vector lambda, vector Index, 
	const double rho_f, vector prefcenter ) 
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
  Cache* Interior[NPARTICLES];
  Cache* Boundary[NPARTICLES];
  SolidBodyBoundary * sbm;
  
  /* Loop over all particles */
  for (size_t k = 0; k < n; k++) 
  {
    foreach_dimension() 
    {
      lambdasumint.x = 0.;
      lambdasumboundary.x = 0.;
      lambdasum.x = 0;
      
      crossLambdaSumInt.x = 0.;
      crossLambdaSumBoundary.x = 0.;
      crossLambdaSum.x = 0.;
    }

#   if dimension == 2 
      lambdasumint.z = 0.;
      lambdasumboundary.z = 0.;
      lambdasum.z = 0;
      
      crossLambdaSumInt.z = 0.;
      crossLambdaSumBoundary.z = 0.;
      crossLambdaSum.z = 0.;
#   endif

    /* Interior domain of the particle */
    Interior[k] = &(allparticles[k].Interior);

    /* Boundary domain of the particle */
    Boundary[k] = &(allparticles[k].reduced_domain);
    sbm = &(allparticles[k].s);
    coord lambdapos;
#   if DLM_Moving_particle
      double rho_s = allparticles[k].rho_s;
#     if TRANSLATION
        double M = allparticles[k].M;
        coord Unm1 = allparticles[k].Unm1;
        coord U = allparticles[k].U;
#     endif
#     if ROTATION
        coord wnm1 = allparticles[k].wnm1;
        coord w = allparticles[k].w;
        coord Iom, omIom, Idomdt;
#     endif
#   endif
    
    // Particle's interior multipliers
    foreach_cache((*Interior[k])) 
    {
      if ( Flag[] < 1 && k == Index.y[] ) 
      {
	/* Compute Fh = <lambda,V>_P */
	/* For moving particles the additional term +
	   (rho_f/rho_s).M.dU/dt is added after the mpi_calls below */
	
	foreach_dimension() 
	  lambdasumint.x += lambda.x[];

	/* Compute Mh = <lambda,xi^GM>_P */  
	/* For moving particles, the additional term +
	   (rho_f/rho_s).(I.dom/dt + om ^ (I.om)) is added after mpi
	   calls below */
#       if dimension == 3
	  crossLambdaSumInt.x += ( lambda.z[] * ( y - prefcenter.y[] )
		- lambda.y[] * ( z - prefcenter.z[] ) );
	  crossLambdaSumInt.y += ( lambda.x[] * ( z - prefcenter.z[] )
		- lambda.z[] * ( x - prefcenter.x[] ) );
#       endif
	crossLambdaSumInt.z += ( lambda.y[] * ( x - prefcenter.x[] )
		- lambda.x[] * ( y - prefcenter.y[] ) );
      }
    }

    // Particle's boundary multipliers
    foreach_cache((*Boundary[k])) 
    {
      if ( Index.x[] > -1 && allparticles[k].pnum == (int)Index.y[] ) 
      {
	lambdapos.x = (*sbm).x[(int)Index.x[]];
	lambdapos.y = (*sbm).y[(int)Index.x[]];
	lambdapos.z = 0.;
#       if dimension == 3
	  lambdapos.z = (*sbm).z[(int)Index.x[]];
#       endif

	/* Compute Fh = <lambda,V>_P */
	foreach_dimension()
	  lambdasumboundary.x += lambda.x[];

	/* Modify temporerly the particle center position for periodic 
	boundary condition */
#       if dimension == 3
          crossLambdaSumBoundary.x += 
		( lambda.z[] * ( lambdapos.y - prefcenter.y[] )
		- lambda.y[] * ( lambdapos.z - prefcenter.z[] ) );
          crossLambdaSumBoundary.y += 
		( lambda.x[] * ( lambdapos.z - prefcenter.z[] )
		- lambda.z[] * ( lambdapos.x - prefcenter.x[] ) );
#       endif		
        crossLambdaSumBoundary.z += 
		( lambda.y[] * ( lambdapos.x - prefcenter.x[] )
		- lambda.x[] * ( lambdapos.y - prefcenter.y[] ) );
      }
    }
 
#   if _MPI
#     if dimension == 3
        foreach_dimension() 
        {
          mpi_all_reduce( lambdasumint.x, MPI_DOUBLE, MPI_SUM );
          mpi_all_reduce( lambdasumboundary.x, MPI_DOUBLE, MPI_SUM );
          mpi_all_reduce( crossLambdaSumInt.x, MPI_DOUBLE, MPI_SUM );
          mpi_all_reduce( crossLambdaSumBoundary.x, MPI_DOUBLE, MPI_SUM );
        }
#     else
        mpi_all_reduce( lambdasumint.x, MPI_DOUBLE, MPI_SUM );
        mpi_all_reduce( lambdasumboundary.x, MPI_DOUBLE, MPI_SUM );
        mpi_all_reduce( lambdasumint.y, MPI_DOUBLE, MPI_SUM );
        mpi_all_reduce( lambdasumboundary.y, MPI_DOUBLE, MPI_SUM );
        mpi_all_reduce( crossLambdaSumInt.z, MPI_DOUBLE, MPI_SUM );
        mpi_all_reduce( crossLambdaSumBoundary.z, MPI_DOUBLE, MPI_SUM );
#     endif
#   endif

    /* The Force term (rho_f/rho_s)MdU/dt is added here for
       translating particles */
    
#   if DLM_Moving_particle

    /* The Force term (rho_f/rho_s)MdU/dt is added here for
       translating particles */

#     if TRANSLATION
        foreach_dimension()
          lambdasumint.x += ( rho_f / rho_s ) * M * ( U.x - Unm1.x ) / dt;
#     endif
    
    /* The Torque term (rho_f/rho_s).(Idom/dt + om ^ (I.om)) is added
	 here for rotating particles */
    
#     if ROTATION
        /* I.om term */
#       if dimension == 3
          Iom.x = allparticles[k].Ip[0] * w.x 
	  	+ allparticles[k].Ip[3] * w.y + allparticles[k].Ip[4] * w.z;
          Iom.y = allparticles[k].Ip[3] * w.x 
	  	+ allparticles[k].Ip[1] * w.y + allparticles[k].Ip[5] * w.z;
          Iom.z = allparticles[k].Ip[4] * w.x 
	  	+ allparticles[k].Ip[5] * w.y + allparticles[k].Ip[2] * w.z;
#       else
          Iom.z = allparticles[k].Ip[2] * w.z;
#       endif	 

        /* om^(I.om) term */
#       if dimension == 3
          omIom.x = w.y * Iom.z - w.z * Iom.y;
          omIom.y = w.z * Iom.x - w.x * Iom.z;
          omIom.z = w.x * Iom.y - w.y * Iom.x;
#       else
          omIom.z = 0.;
#       endif	  

        /* Idom/dt term */
#       if dimension == 3	
          Idomdt.x = allparticles[k].Ip[0] * ( w.x - wnm1.x ) 
	  	+ allparticles[k].Ip[3] * ( w.y - wnm1.y ) 
    		+ allparticles[k].Ip[4] * ( w.z - wnm1.z );
          Idomdt.y = allparticles[k].Ip[3] * ( w.x - wnm1.x ) 
	  	+ allparticles[k].Ip[1] * ( w.y - wnm1.y ) 
    		+ allparticles[k].Ip[5] * ( w.z - wnm1.z );
          Idomdt.z = allparticles[k].Ip[4] * ( w.x - wnm1.x ) 
	  	+ allparticles[k].Ip[5] * ( w.y - wnm1.y ) 
    		+ allparticles[k].Ip[2] * ( w.z - wnm1.z );
#       else
          Idomdt.z = allparticles[k].Ip[2] * ( w.z - wnm1.z );
#       endif		
	
#       if dimension == 3
          crossLambdaSumInt.x += ( rho_f / rho_s ) * ( ( Idomdt.x ) / dt 
		+ omIom.x );
          crossLambdaSumInt.y += ( rho_f / rho_s ) * ( ( Idomdt.y ) / dt 
		+ omIom.y );
#       endif		
        crossLambdaSumInt.z += ( rho_f / rho_s ) * ( ( Idomdt.z ) / dt 
		+ omIom.z );
#     endif	  
#   endif

#   if dimension == 3    
      foreach_dimension() 
      {
        lambdasum.x = lambdasumint.x + lambdasumboundary.x;
        crossLambdaSum.x = crossLambdaSumInt.x + crossLambdaSumBoundary.x;
      }
#   else
      lambdasum.x = lambdasumint.x + lambdasumboundary.x;
      lambdasum.y = lambdasumint.y + lambdasumboundary.y;            
      crossLambdaSum.z = crossLambdaSumInt.z + crossLambdaSumBoundary.z;
#   endif
      
    if( pid() == 0 ) 
#     if dimension == 3  
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
#     else
      {
        fprintf( sl[k], "%.*e\t", NSDF, t );
        fprintf( sl[k], "%.*e\t", NSDF, lambdasum.x );      
        fprintf( sl[k], "%.*e\t", NSDF, lambdasum.y );          
        fprintf( sl[k], "%.*e\n", NSDF, crossLambdaSum.z );      
        fflush( sl[k] );
      }
#endif
  }
}




/** Initialize/open all DLMFD file pointers */
//----------------------------------------------------------------------------
void init_file_pointers( const size_t n, FILE** p, FILE** d, FILE** UzawaCV, 
	FILE** CVT, const size_t rflag ) 
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
      for (size_t k = 0; k < n; k++) 
      {
        sprintf( suffix, "_%lu.dat", k );

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
void close_file_pointers( const size_t n, FILE** p, FILE** d, FILE* UzawaCV, 
	FILE* CVT ) 
//----------------------------------------------------------------------------
{ 
# if _MPI
    if ( pid() == 0 )
# endif
    {
      // Particle data
      for (size_t k = 0; k < n; k++) 
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
  
  synchronize((scalar *){omega});
}




/** Computes the flow rate on the right boundary (x+ direction) */
//----------------------------------------------------------------------------
double compute_flowrate_right( const vector u, const int level ) 
//----------------------------------------------------------------------------
{
  double flowrate = 0.;

# if !adaptive 
    Cache intDomain = {0};
    Point lpoint;
    foreach() 
      if ( x > L0 - Delta ) 
      { 
	lpoint = locate( x, y, z );
	cache_append( &intDomain, lpoint, 0 );
      }
		
    cache_shrink( &intDomain );
    
    foreach_cache(intDomain)
      flowrate += sq(Delta) * u.x[];
  
    free( intDomain.p );
# endif

# if adaptive
    double hh = L0 / pow(2,level);
    double zi = 0., yj = 0., xval = 0., uinter = 0.;
    int ii = 0, jj = 0;
  
    foreach_level(level) xval = L0 - 0.5 * Delta + X0;
  
    for (ii = 0; ii < pow(2,level); ii++) 
    {
      zi = 0.5 * hh + ii * hh + Z0;
      for (jj = 0; jj < pow(2,level); jj++) 
      {
        yj = 0.5 * hh + jj * hh + Y0;
        uinter = interpolate( u.x, xval, yj, zi );
	if ( uinter < nodata )
	  flowrate += sq(hh) * uinter;
      }
    }  
# endif

# if _MPI
    mpi_all_reduce( flowrate, MPI_DOUBLE, MPI_SUM );
# endif
 
  return flowrate;
}




/** Computes the flow rate on the top boundary (y+ direction) */
//----------------------------------------------------------------------------
double compute_flowrate_top( const vector u, const int level ) 
//----------------------------------------------------------------------------
{
  double flowrate = 0.;

# if !adaptive 
    Cache intDomain = {0};
    Point lpoint;
    foreach() 
      if ( y > L0 - Delta ) 
      { 
	lpoint = locate( x, y, z );
	cache_append( &intDomain, lpoint, 0 );
      }
		
    cache_shrink( &intDomain );
    
    foreach_cache(intDomain)
      flowrate += sq(Delta) * u.y[];
  
    free( intDomain.p );
# endif

# if adaptive
    double hh = L0 / pow(2,level);
    double zi = 0., yval = 0., xj = 0., uinter = 0.;
    int ii = 0, jj = 0;
  
    foreach_level(level) yval = L0 - 0.5 * Delta + Y0;
  
    for (ii = 0; ii < pow(2,level); ii++) 
    {
      zi = 0.5 * hh + ii * hh + Z0;
      for (jj = 0; jj < pow(2,level); jj++) 
      {
        xj = 0.5 * hh + jj * hh + X0;
        uinter = interpolate( u.y, xj, yval, zi );
	if ( uinter < nodata )
	  flowrate += sq(hh) * uinter;
      }
    }  
# endif

# if _MPI
    mpi_all_reduce( flowrate, MPI_DOUBLE, MPI_SUM );
# endif
 
  return flowrate;
}




/** Computes the flow rate on the front boundary (z+ direction) */
//----------------------------------------------------------------------------
double compute_flowrate_front( const vector u, const int level ) 
//----------------------------------------------------------------------------
{
  double flowrate = 0.;

# if !adaptive 
    Cache intDomain = {0};
    Point lpoint;
    foreach() 
      if ( z > L0 - Delta ) 
      { 
	lpoint = locate( x, y, z );
	cache_append( &intDomain, lpoint, 0 );
      }
		
    cache_shrink( &intDomain );
    
    foreach_cache(intDomain)
      flowrate += sq(Delta) * u.z[];
  
    free( intDomain.p );
# endif

# if adaptive
    double hh = L0 / pow(2,level);
    double xi = 0., yj = 0., zval = 0., uinter = 0.;
    int ii = 0, jj = 0;
  
    foreach_level(level) zval = L0 - 0.5 * Delta + Z0;
  
    for (ii = 0; ii < pow(2,level); ii++) 
    {
      xi = 0.5 * hh + ii * hh + X0;
      for (jj = 0; jj < pow(2,level); jj++) 
      {
        yj = 0.5 * hh + jj * hh + Y0;
        uinter = interpolate ( u.z, xi, yj, zval );
	if ( uinter < nodata )
	  flowrate += sq(hh) * uinter;
      }
    }  
# endif

# if _MPI
    mpi_all_reduce( flowrate, MPI_DOUBLE, MPI_SUM );
# endif
 
  return flowrate;
}




/** Gets the pressure at a point and writes the value in a file */
//----------------------------------------------------------------------------
void pressure_at_point( scalar pres, FILE* ppf, const coord hv, 
	const double t )
//----------------------------------------------------------------------------
{
  double plaw = 0.;
  
  plaw = interpolate (pres, hv.x, hv.y, hv.z);

# if _MPI
    MPI_Barrier(MPI_COMM_WORLD);
    if ( plaw == nodata ) plaw = 0.;
    MPI_Barrier(MPI_COMM_WORLD);

    if ( plaw != nodata )
      mpi_all_reduce(plaw, MPI_DOUBLE, MPI_SUM);
#  endif
 
  if ( pid() == 0 ) 
  {
    fprintf( ppf, "%20.18f\t %20.18f\n", t, plaw );
    fflush( ppf );
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
void compute_inv_inertia( particle* p )
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
  foreach(reduction(+:t)) t++;  
  
  return t;
}




/** Computes and returns the total number of cells related to Distributed
Lagrange multiplier points for all rigid bodies */
//----------------------------------------------------------------------------
int total_dlmfd_cells( particle* allparticles, const size_t np ) 
//----------------------------------------------------------------------------
{
  int apts = 0;
  
  for (size_t k = 0; k < np; k++) 
  {
#   if debugInterior == 0
      apts += allparticles[k].Interior.n;
#   endif
#   if debugBD == 0
      apts += allparticles[k].reduced_domain.n;
#   endif
  }

# if _MPI
    mpi_all_reduce( apts, MPI_INT, MPI_SUM );
# endif

  return apts;
}




/** Computes and returns the total number of Ditributed Lagrange multiplier 
points for all rigid bodies */
//----------------------------------------------------------------------------
int total_dlmfd_multipliers( particle* allparticles, const size_t np ) 
//----------------------------------------------------------------------------
{
  int apts = 0;
  
# if debugInterior == 0
    for (size_t k = 0; k < np; k++) 
      apts += allparticles[k].Interior.n;
  
#   if _MPI
      mpi_all_reduce (apts, MPI_INT, MPI_SUM);
#   endif
# endif
  
# if debugBD == 0
    for (int k = 0; k < np; k++) 
    {
      SolidBodyBoundary * bla;
      bla = &(allparticles[k].s);
      apts += bla->m;
    }
# endif
  
  return apts;
}




/** Save the time and time step in a file */
//----------------------------------------------------------------------------
void save_t_dt_restart( char* dirname, double time, double deltat, double ppd )
//----------------------------------------------------------------------------
{
  char dump_name[80] = "";
  strcpy( dump_name, dirname );
  strcat( dump_name, "/t_dt_restart.res" );
  FILE* ft = fopen( dump_name, "w" );
  fprintf ( ft, "%.10e %.10e %.10e", time, deltat, ppd );
  fclose( ft );  
}




/** Read the restart time and time step from a file */
//----------------------------------------------------------------------------
void read_t_restart( char* dirname, double* time, double* deltat, double* ppd )
//----------------------------------------------------------------------------
{
  char dump_name[80] = "";
  strcpy( dump_name, dirname );
  strcat( dump_name, "/t_dt_restart.res" );
  FILE* ft = fopen( dump_name, "r" );
  fscanf ( ft, "%lf %lf %lf", time, deltat, ppd );
  fclose( ft );  
}




/** Allocate number of particles dependent arrays */
//----------------------------------------------------------------------------
void allocate_np_dep_arrays( const size_t npart, particle** particles_, 
	double*** DLMFDtoGS_vel_, double** vpartbuf_, FILE*** pdata_, 
	FILE*** fdata_ )
//----------------------------------------------------------------------------
{
  *particles_ = (particle*) calloc( npart, sizeof(particle) );
  *DLMFDtoGS_vel_ = (double**) calloc( npart, sizeof(double*) );
  for (size_t k=0;k<npart;++k)
    (*DLMFDtoGS_vel_)[k] = (double*) calloc( 6, sizeof(double) );
  *vpartbuf_ = (double*) calloc( npartdata, npart * sizeof(double) );
  *pdata_ = (FILE**) calloc( npart, sizeof( FILE* ) );
  *fdata_ = (FILE**) calloc( npart, sizeof( FILE* ) );    
}




/** Free number of particles dependent arrays */
//----------------------------------------------------------------------------
void free_np_dep_arrays( const size_t npart, particle* particles_, 
	double** DLMFDtoGS_vel_, double* vpartbuf_, FILE** pdata_, 
	FILE** fdata_ )
//----------------------------------------------------------------------------
{
  free( particles_ );
  for (size_t k=0;k<npart;++k) free( DLMFDtoGS_vel_[k] );
  free( DLMFDtoGS_vel_ );
  free( vpartbuf_ );
  free( pdata_ );
  free( fdata_ );  
}
