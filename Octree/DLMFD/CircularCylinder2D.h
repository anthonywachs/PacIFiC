/** 
# Set of functions for a 2D circular cylinder 
*/


/** Tests whether a point lies inside the 2D circular cylinder */
//----------------------------------------------------------------------------
bool is_in_CircularCylinder2D( const double x, const double y, 
	const GeomParameter gp )
//----------------------------------------------------------------------------
{
  return ( sqrt( sq( x - gp.center.x ) + sq( y - gp.center.y ) ) < gp.radius );
}




/** Computes the number of boundary points on the surface of the 2D circular 
cylinder */
//----------------------------------------------------------------------------
void compute_nboundary_CircularCylinder2D( const GeomParameter gcp, int* nb ) 
//----------------------------------------------------------------------------
{
  coord pos[6];
  Cache poscache = {0};
  Point lpoint;
   
  pos[0].x = gcp.center.x + gcp.radius + X0;
  pos[0].y = gcp.center.y + Y0;

  pos[1].x = gcp.center.x - gcp.radius + X0;
  pos[1].y = gcp.center.y + Y0;

  pos[2].x = gcp.center.x + X0;
  pos[2].y = gcp.center.y + gcp.radius + Y0;

  pos[3].x = gcp.center.x + X0;
  pos[3].y = gcp.center.y - gcp.radius + Y0;
  
  *nb = 0;
  lpoint.level = -1;
  int ip = 0;
  
  /* We assume that around the circular cylinder at least one point has to be 
     at maximum level of refinement. */
  for (ip = 0; ip < dimension*2; ip++) 
  {
    lpoint = locate( pos[ip].x, pos[ip].y );
    if ( lpoint.level == depth() )
      break;
  }
  
  /** Multiple threads can have a different point at the same level of
      refinement. This is fine, they will do the same computation for
      *nb. This problem does not exist in serial but the code below
      still works. */
  if ( lpoint.level == depth() ) 
  {           
    /** Only this thread creates the Cache ... */
    cache_append( &poscache, lpoint, 0 );

    /** and only this thread computes the number of boundary points ... */
    foreach_cache(poscache) 
      *nb += (int)( floor( gcp.radius * 2. * pi / ( sqrt(2.) * Delta ) ) );

    /** and finally, this thread destroys the cache. */    
    free( poscache.p );
  }

#if _MPI
  MPI_Barrier( MPI_COMM_WORLD );
  mpi_all_reduce ( *nb, MPI_INT, MPI_MAX );
#endif
      
  if( *nb == 0 )
    printf( "nboundary = 0: No boundary points !!!\n" );
}




/** Creates boundary points on the surface of the 2D circular cylinder */
//----------------------------------------------------------------------------
void create_FD_Boundary_CircularCylinder2D( GeomParameter gcp, 
	SolidBodyBoundary* dlm_bd, const int nsphere, vector pshift ) 
//----------------------------------------------------------------------------
{

  int m = nsphere;
  coord pos;
  pos.z = 0.;
  Point lpoint;
  
  coord shift = {0., 0., 0.};
  coord ori = {X0, Y0, Z0};
 
  int i;
  double theta[m];
  double radius  = gcp.radius;
  coord fact_theta[m];
  
  for (i = 0; i < m; i++) 
  {
    theta[i] = (double)(i) * 2. * pi / (double)(m); 
    fact_theta[i].x = cos( theta[i] );
    fact_theta[i].y = sin( theta[i] );
    
    /* Assign positions x,y on a circle's boundary */ 
    pos.x = radius * fact_theta[i].x + gcp.center.x + X0;
    pos.y = radius * fact_theta[i].y + gcp.center.y + Y0;

    /* Check if the point falls outside of the domain */
    foreach_dimension() 
    {
      if ( Period.x ) 
      {
	shift.x = 0.;
	if ( pos.x > L0 + ori.x ) 
	{
	  pos.x -= L0;
	  shift.x = -L0;
	} 
	if ( pos.x < 0. + ori.x ) 
	{
	  pos.x += L0;
	  shift.x = L0;
	}
      }
      dlm_bd->x[i] = pos.x; 
    }
   
    Cache poscache = {0};
    lpoint = locate( pos.x, pos.y );

    if ( lpoint.level > -1 ) 
    {	
      cache_append( &poscache, lpoint, 0 );

      foreach_cache(poscache) 
	foreach_dimension()
	  pshift.x[] = shift.x;
      
      free( poscache.p );
    }
  }
  
  boundary((scalar*){pshift});
}




/** Finds cells lying inside the 2D circular cylinder */
//----------------------------------------------------------------------------
void create_FD_Interior_CircularCylinder2D( particle* p, vector Index_lambda, 
	vector shift, scalar flag )
//----------------------------------------------------------------------------
{
  GeomParameter gci = p->g;
  Cache* fd = &(p->Interior);
  
  /** Create the cache for the interior points */
  foreach() 
  {
    if ((sq(x - (gci.center.x + X0)) + sq(y - (gci.center.y + Y0))) 
	< sq(gci.radius)) 
    {
      cache_append (fd, point, 0);
      
      /* tag cell with the number of the particle */
      if ((int)Index_lambda.y[] == -1)
	Index_lambda.y[] = p->pnum; 
    }

    if (Period.x) {
      if ((sq(x - L0 - (gci.center.x + X0)) + sq(y - (gci.center.y + Y0))) 
	< sq(gci.radius)) 
      {
    	cache_append (fd, point, 0);
	
    	/* tag cell with the number of the particle */
    	if ((int)Index_lambda.y[] == -1) 
    	  Index_lambda.y[] = p->pnum;
	if ((int)Index_lambda.y[] == p->pnum)
	  shift.x[] = L0;
      }
      if ((sq(x + L0 - (gci.center.x + X0)) + sq(y - (gci.center.y + Y0))) 
      	< sq(gci.radius)) 
      {
    	cache_append (fd, point, 0);

	/* tag cell with the number of the particle */
    	if ((int)Index_lambda.y[] == -1) 
    	  Index_lambda.y[] = p->pnum;
	if ((int)Index_lambda.y[] == p->pnum)
	  shift.x[] = -L0;
      }
    }
   
    if (Period.y) {
      if ((sq(x - (gci.center.x + X0)) + sq(y - L0 - (gci.center.y + Y0))) 
      	< sq(gci.radius)) 
      {
    	cache_append (fd, point, 0);
	
	/* tag cell with the number of the particle */
	if ((int)Index_lambda.y[] == -1) {
	  Index_lambda.y[] = p->pnum;
	  if (flag[] < 1)
	    shift.y[] = L0;
	}
      }
      if ((sq(x - (gci.center.x + X0)) + sq(y + L0 - (gci.center.y + Y0))) 
      	< sq(gci.radius)) 
      {
    	cache_append (fd, point, 0);
	
	/* tag cell with the number of the particle */
	if ((int)Index_lambda.y[] == -1) {
	  Index_lambda.y[] = p->pnum;
	  if (flag[] < 1)
	    shift.y[] = -L0;
	}
      }
    }        
  }
  cache_shrink( fd );
}




/** Reads geometric parameters of the 2D circular cylinder */
//----------------------------------------------------------------------------
void update_CircularCylinder2D( GeomParameter* gcp ) 
//----------------------------------------------------------------------------
{    
  // TO DO
}




/** Frees the geometric parameters of the 2D circular cylinder */
//----------------------------------------------------------------------------
void free_CircularCylinder2D( GeomParameter* gcp ) 
//----------------------------------------------------------------------------
{  
  // Nothing to do
}
