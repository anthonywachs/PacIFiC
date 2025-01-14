/** 
# Set of functions for a 2D circular cylinder 
*/


/** Tests whether a point lies inside the 2D circular cylinder */
//----------------------------------------------------------------------------
bool is_in_CircularCylinder2D_clone( const double x, const double y, 
	const GeomParameter gp )
//----------------------------------------------------------------------------
{
  return ( sqrt( sq( x - gp.center.x ) + sq( y - gp.center.y ) ) < gp.radius );
}




/** Tests whether a point lies inside the 2D circular cylinder or any of its 
periodic clones */
//----------------------------------------------------------------------------
bool is_in_CircularCylinder2D( const double x, const double y, 
	const GeomParameter gp )
//----------------------------------------------------------------------------
{
  // Check if it is in the master particle
  bool status = is_in_CircularCylinder2D_clone( x, y, gp );

  // Check if it is in any clone particle
  if ( gp.nperclones && !status )
    for (int i = 0; i < gp.nperclones && !status; i++)
    {
      GeomParameter clones = gp;
      clones.center = gp.perclonecenters[i];
      status = is_in_CircularCylinder2D_clone( x, y, clones );
    }

  return status;
}




/** Tests whether a point lies inside the 2D circular cylinder or any of its 
periodic clones and assign the proper center of mass coordinates associated to
this point */
//----------------------------------------------------------------------------
bool in_which_CircularCylinder2D( double x1, double y1, 
	const GeomParameter gp, vector PeriodicRefCenter, 
	const bool setPeriodicRefCenter )
//----------------------------------------------------------------------------
{
  // Check if it is in the master particle
  bool status = is_in_CircularCylinder2D_clone( x1, y1, gp );
  if ( status )
  {
    Cache poscache = {0};
    Point lpoint;
    lpoint = locate( x1, y1 );

    if ( lpoint.level > -1 )
    {
      cache_append( &poscache, lpoint, 0 );
      foreach_cache(poscache)
        foreach_dimension()
          PeriodicRefCenter.x[] = gp.center.x;
      free( poscache.p );
    }
  }

  //  Check if it is in any clone particle
  if ( gp.nperclones && !status )
    for (int i = 0; i < gp.nperclones && !status; i++) 
    {
      GeomParameter clones = gp;
      clones.center = gp.perclonecenters[i];
      status = is_in_CircularCylinder2D_clone( x1, y1, clones );
      if ( status && setPeriodicRefCenter )
      {
        Cache poscache = {0};
	Point lpoint;
	lpoint = locate( x1, y1 );
	
	if ( lpoint.level > -1 )
	{
	  cache_append( &poscache, lpoint, 0 );
	  foreach_cache(poscache)
	    foreach_dimension()
	      PeriodicRefCenter.x[] = clones.center.x;
          free( poscache.p );
        }
      }
    }

  return status;
}




/** Computes the number of boundary points on the perimeter of the 2D circular 
cylinder */
//----------------------------------------------------------------------------
void compute_nboundary_CircularCylinder2D( const GeomParameter gcp, int* nb ) 
//----------------------------------------------------------------------------
{
  double delta = L0 / (double)(1 << MAXLEVEL) ; 

  *nb = (int)( floor( 2. * pi * gcp.radius
	/ ( INTERBPCOEF * delta ) ) );
      
  if( *nb == 0 )
    printf( "nboundary = 0: No boundary points !!!\n" );
}




/** Creates boundary points on the surface of the 2D circular cylinder */
//----------------------------------------------------------------------------
void create_FD_Boundary_CircularCylinder2D( GeomParameter gcp,
	SolidBodyBoundary* dlm_bd, const int nsphere, 
	vector* pPeriodicRefCenter, const bool setPeriodicRefCenter  ) 
//----------------------------------------------------------------------------
{
  int m = nsphere;
  coord pos; pos.z = 0.;
  Point lpoint;
  
  coord shift = {0., 0., 0.};
  coord ori = {X0, Y0, 0.};
 
  int i;
  double theta[m];
  double radius = gcp.radius;
  coord fact_theta[m];
  
  for (i = 0; i < m; i++) 
  {
    theta[i] = (double)(i) * 2. * pi / (double)(m); 
    fact_theta[i].x = cos( theta[i] );
    fact_theta[i].y = sin( theta[i] );
    
    /* Assign positions x, y on the circular cylinder boundary */ 
    pos.x = radius * fact_theta[i].x + gcp.center.x + X0;
    pos.y = radius * fact_theta[i].y + gcp.center.y + Y0;

    /* Check if the point falls outside of the domain */    
    foreach_dimension()
    {
      shift.x = 0.;      
      if ( Period.x )
      {
	if ( pos.x > L0 + ori.x )
	{  
	  pos.x -= L0;
	  shift.x = - L0;
	}
        else if ( pos.x < 0. + ori.x )
	{
	  pos.x += L0;
	  shift.x = L0;
	}
      }
    }
   
    if ( setPeriodicRefCenter )
    {
      // Setting the periodic clone center vector field
      Cache poscache = {0};
      lpoint = locate( pos.x, pos.y );

      if ( lpoint.level > -1 )
      {
        cache_append( &poscache, lpoint, 0 );
        foreach_cache(poscache)
          foreach_dimension()
            pPeriodicRefCenter->x[] = gcp.center.x + shift.x;
        free( poscache.p );
      }
    }

    dlm_bd->x[i] = pos.x;
    dlm_bd->y[i] = pos.y;
  }

  if ( setPeriodicRefCenter ) synchronize((scalar*){pPeriodicRefCenter->x,
  	pPeriodicRefCenter->y, pPeriodicRefCenter->z});
}




/** Finds cells lying inside the 2D circular cylinder */
//----------------------------------------------------------------------------
void create_FD_Interior_CircularCylinder2D( particle* p, vector Index,
	vector PeriodicRefCenter )
//----------------------------------------------------------------------------
{
  GeomParameter gci = p->g;
  Cache* fd = &(p->Interior);

  /** Create the cache for the interior points */
  foreach()
  {
    bool in_clone;
    in_clone = in_which_CircularCylinder2D( x, y, gci, PeriodicRefCenter, 
    	true );
    if ( in_clone )
    {
      cache_append( fd, point, 0 );
      if ( (int)Index.y[] == -1 )
        Index.y[] = p->pnum;
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
