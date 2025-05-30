/** 
# Set of functions for a 2D circular cylinder 
*/


/** Tests whether a point lies inside the 2D circular cylinder */
//----------------------------------------------------------------------------
bool is_in_CircularCylinder2D_clone( const double x, const double y, 
	GeomParameter const* gcp )
//----------------------------------------------------------------------------
{
  return ( sqrt( sq( x - gcp->center.x ) + sq( y - gcp->center.y ) ) 
  	< gcp->radius );
}




/** Tests whether a point lies inside the 2D circular cylinder or any of its 
periodic clones */
//----------------------------------------------------------------------------
bool is_in_CircularCylinder2D( const double x, const double y, 
	GeomParameter const* gcp )
//----------------------------------------------------------------------------
{
  // Check if it is in the master rigid body
  bool status = is_in_CircularCylinder2D_clone( x, y, gcp );

  // Check if it is in any clone rigid body
  if ( gcp->nperclones && !status )
    for (int i = 0; i < gcp->nperclones && !status; i++)
    {
      GeomParameter clone = *gcp;
      clone.center = gcp->perclonecenters[i];
      status = is_in_CircularCylinder2D_clone( x, y, &clone );
    }

  return ( status );
}




/** Tests whether a point lies inside the 2D circular cylinder or any of its 
periodic clones and assign the proper center of mass coordinates associated to
this point */
//----------------------------------------------------------------------------
bool in_which_CircularCylinder2D( double x1, double y1, 
	GeomParameter const* gcp, vector* pPeriodicRefCenter, 
	const bool setPeriodicRefCenter )
//----------------------------------------------------------------------------
{
  // Check if it is in the master rigid body
  bool status = is_in_CircularCylinder2D_clone( x1, y1, gcp );
  if ( status && setPeriodicRefCenter )
    foreach_point( x1, y1 )
      foreach_dimension()
        pPeriodicRefCenter->x[] = gcp->center.x;

  //  Check if it is in any clone rigid body
  if ( gcp->nperclones && !status )
    for (int i = 0; i < gcp->nperclones && !status; i++) 
    {
      GeomParameter clone = *gcp;
      clone.center = gcp->perclonecenters[i];
      status = is_in_CircularCylinder2D_clone( x1, y1, &clone );
      if ( status && setPeriodicRefCenter )
        foreach_point( x1, y1 )
          foreach_dimension()
	    pPeriodicRefCenter->x[] = clone.center.x;
    }

  return ( status );
}




/** Computes the number of boundary points on the perimeter of the 2D circular 
cylinder */
//----------------------------------------------------------------------------
void compute_nboundary_CircularCylinder2D( GeomParameter const* gcp, int* nb ) 
//----------------------------------------------------------------------------
{
  double delta = L0 / (double)(1 << MAXLEVEL) ; 

  *nb = (int)( floor( 2. * pi * gcp->radius
	/ ( INTERBPCOEF * delta ) ) );
      
  if( *nb == 0 )
    printf( "nboundary = 0: No boundary points !!!\n" );
}




/** Creates boundary points on the surface of the 2D circular cylinder */
//----------------------------------------------------------------------------
void create_FD_Boundary_CircularCylinder2D( GeomParameter const* gcp,
	RigidBodyBoundary* dlm_bd, const int nsphere, 
	vector* pPeriodicRefCenter, const bool setPeriodicRefCenter  ) 
//----------------------------------------------------------------------------
{
  int m = nsphere;
  coord pos; pos.z = 0.;  
  coord shift = {0., 0., 0.};
  coord ori = {X0, Y0, 0.};
 
  int i;
  double theta[m];
  double radius = gcp->radius;
  coord fact_theta[m];
  
  for (i = 0; i < m; i++) 
  {
    theta[i] = (double)(i) * 2. * pi / (double)(m); 
    fact_theta[i].x = cos( theta[i] );
    fact_theta[i].y = sin( theta[i] );
    
    /* Assign positions x, y on the circular cylinder boundary */ 
    pos.x = radius * fact_theta[i].x + gcp->center.x + X0;
    pos.y = radius * fact_theta[i].y + gcp->center.y + Y0;

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
      foreach_point( pos.x, pos.y )
        foreach_dimension()
	  pPeriodicRefCenter->x[] = gcp->center.x + shift.x;
    }

    dlm_bd->x[i] = pos.x;
    dlm_bd->y[i] = pos.y;
  }

  if ( setPeriodicRefCenter ) synchronize((scalar*){pPeriodicRefCenter->x,
  	pPeriodicRefCenter->y});
}




/** Finds cells lying inside the 2D circular cylinder */
//----------------------------------------------------------------------------
void create_FD_Interior_CircularCylinder2D( RigidBody* p, vector Index,
	vector PeriodicRefCenter )
//----------------------------------------------------------------------------
{
  GeomParameter const* gcp = &(p->g);
  Cache* fd = &(p->Interior);
  Point ppp;

  /** Create the cache for the interior points */
  foreach(serial)
    if ( in_which_CircularCylinder2D( x, y, gcp, &PeriodicRefCenter, true ) )
      if ( (int)Index.y[] == -1 )
      {
	ppp.i = point.i;
        ppp.j = point.j;		
        ppp.level = point.level;
	cache_append( fd, ppp, 0 );
        Index.y[] = p->pnum;
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
