/**
# Set of functions for a sphere
*/


/** Tests whether a point lies inside the sphere */
//----------------------------------------------------------------------------
bool is_in_Sphere_clone( const double x, const double y, const double z,
	GeomParameter const* gcp )
//----------------------------------------------------------------------------
{
  return ( sqrt( sq( x - gcp->center.x ) + sq( y - gcp->center.y )
    	+ sq( z - gcp->center.z ) ) < gcp->radius );
}




/** Tests whether a point lies inside the sphere or any of its 
periodic clones */
//----------------------------------------------------------------------------
bool is_in_Sphere( const double x, const double y, const double z,
	GeomParameter const* gcp )
//----------------------------------------------------------------------------
{
  // Check if it is in the master rigid body
  bool status = is_in_Sphere_clone( x, y, z, gcp );

  // Check if it is in any clone rigid body
  if ( gcp->nperclones && !status )
    for (int i = 0; i < gcp->nperclones && !status; i++)
    {
      GeomParameter clone = *gcp;
      clone.center = gcp->perclonecenters[i];
      status = is_in_Sphere_clone( x, y, z, &clone );
    }

  return ( status );
}




/** Tests whether a point lies inside the sphere or any of its 
periodic clones and assign the proper center of mass coordinates associated to
this point */
//----------------------------------------------------------------------------
bool in_which_Sphere( double x1, double y1, double z1,
	GeomParameter const* gcp, vector* pPeriodicRefCenter, 
	const bool setPeriodicRefCenter )
//----------------------------------------------------------------------------
{
  // Check if it is in the master rigid body
  bool status = is_in_Sphere_clone( x1, y1, z1, gcp );
  if ( status && setPeriodicRefCenter )
    foreach_point( x1, y1, z1 )
      foreach_dimension()
        pPeriodicRefCenter->x[] = gcp->center.x;

  //  Check if it is in any clone rigid body
  if ( gcp->nperclones && !status )
    for (int i = 0; i < gcp->nperclones && !status; i++) 
    {
      GeomParameter clone = *gcp;
      clone.center = gcp->perclonecenters[i];
      status = is_in_Sphere_clone( x1, y1, z1, &clone );
      if ( status && setPeriodicRefCenter )
        foreach_point( x1, y1, z1 )
          foreach_dimension()
            pPeriodicRefCenter->x[] = clone.center.x;
    }

  return ( status );
}



/** Computes the number of boundary points on the surface of the sphere */
//----------------------------------------------------------------------------
void compute_nboundary_Sphere( GeomParameter const* gcp, int* nb )
//----------------------------------------------------------------------------
{
  double delta = L0 / (double)(1 << MAXLEVEL) ;  

  *nb = (int)( floor( pow( 3.809 * gcp->radius
	/ ( INTERBPCOEF * delta ), 2.) ) );

  if ( *nb == 0 )
    printf( "nboundary = 0: No boundary points !!!\n" );
}




/** Creates boundary points on the surface of the sphere */
//----------------------------------------------------------------------------
void create_FD_Boundary_Sphere( GeomParameter const* gcp,
	RigidBodyBoundary* dlm_bd, const int nsphere, 
	vector* pPeriodicRefCenter, const bool setPeriodicRefCenter )
//----------------------------------------------------------------------------
{
  coord pos;
  coord shift = {0., 0., 0.};
  coord ori = {X0, Y0, Z0};

  /* Lay out points on a sphere with a Z oriented spiral scheme */
  /* This portion of code is recovered from Peligriff and adapted */
  /* More information of this method can be found at: */
  /* Saff, E. & Kuijlaars, A. Distributing many points on a sphere The
     Mathematical Intelligencer, 1997 */

  double delta = L0 / (double)(1 << MAXLEVEL) ;  
  double spiral_spacing_correction = INTERBPCOEF * delta / sqrt(3.) ;
  double hydro_radius = gcp->radius;
  size_t k;

  /* Number of points */
  size_t NSpiral = nsphere;

  /* Spiral points construction */
  double hk, thetak, phik, phikm1 = 0., TwoPi = 2. * pi,
    Cspiral = 3.6 / sqrt( NSpiral ), dphi = 0.;

  for (k = 0; k < NSpiral; ++k)
  {
    hk = - 1. + 2. * (double)(k) / ( NSpiral - 1. );
    thetak = acos( hk ) ;
    if ( k == 0 )
    {
      phik = 0.;
      thetak = pi - 0.5 * spiral_spacing_correction / hydro_radius ;
    }
    else if ( k == NSpiral - 1 )
    {
      phik = phikm1 + 1. * dphi;
      if ( phik > TwoPi ) phik -= TwoPi ;
      thetak = 0.5 * spiral_spacing_correction / hydro_radius ;
    }
    else
    {
      dphi = Cspiral / sqrt( 1. - hk * hk ) ;
      phik = phikm1 + dphi ;
      if ( phik > TwoPi ) phik -= TwoPi ;
    }

    phikm1 = phik ;

    if ( k == 1 ) thetak -= 0.4 * spiral_spacing_correction / hydro_radius ;

    if ( k == NSpiral - 2 )
    {
      phik -= 0.1 * dphi ;
      thetak += 0.25 * spiral_spacing_correction / hydro_radius ;
    }

    pos.x = hydro_radius * cos( phik ) * sin( thetak ) + gcp->center.x;
    pos.y = hydro_radius * sin( phik ) * sin( thetak ) + gcp->center.y;
    pos.z = hydro_radius * cos( thetak ) + gcp->center.z;

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
      foreach_point( pos.x, pos.y, pos.z )
        foreach_dimension()
	  pPeriodicRefCenter->x[] = gcp->center.x + shift.x;
    }

    foreach_dimension() 
      dlm_bd->x[k] = pos.x;
  }

  if ( setPeriodicRefCenter ) synchronize((scalar*){pPeriodicRefCenter->x,
  	pPeriodicRefCenter->y, pPeriodicRefCenter->z});
}




/** Finds cells lying inside the sphere */
//----------------------------------------------------------------------------
void create_FD_Interior_Sphere( RigidBody* p, vector Index,
	vector PeriodicRefCenter )
//----------------------------------------------------------------------------
{
  GeomParameter* gcp = &(p->g);
  Cache* fd = &(p->Interior);
  Point ppp;
  
  /** Create the cache for the interior points */
  foreach(serial)
    if ( in_which_Sphere( x, y, z, gcp, &PeriodicRefCenter, true ) )
      if ( (int)Index.y[] == -1 )
      {
	ppp.i = point.i;
        ppp.j = point.j;
        ppp.k = point.k;			
        ppp.level = point.level;
	cache_append( fd, ppp, 0 );
        Index.y[] = p->pnum;
      } 

  cache_shrink( fd );
}




/** Reads geometric parameters of the sphere */
//----------------------------------------------------------------------------
void update_Sphere( GeomParameter* gcp )
//----------------------------------------------------------------------------
{
  // We already have all parameters for the sphere but the input array of
  // characters contains again a "1", the gravity center coordinates and
  // a "0", hence we need to read five tokens but we do not do
  // anything with them
  for (size_t i=0;i<5;++i)
    strtok( NULL, " " );
}




/** Frees the geometric parameters of the sphere */
//----------------------------------------------------------------------------
void free_Sphere( GeomParameter* gcp )
//----------------------------------------------------------------------------
{
  // Nothing to do
}
