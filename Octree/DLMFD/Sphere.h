/** 
# Set of functions for a sphere 
*/


/** Tests whether a point lies inside the sphere */
//----------------------------------------------------------------------------
bool is_in_Sphere( const double x, const double y, const double z, 
	const GeomParameter gp )
//----------------------------------------------------------------------------
{
  return ( sqrt( sq( x - gp.center.x ) + sq( y - gp.center.y ) 
    	+ sq( z - gp.center.z ) ) < gp.radius );
}




/** Computes the number of boundary points on the surface of the sphere */
//----------------------------------------------------------------------------
void compute_nboundary_Sphere( const GeomParameter gcp, int* nb ) 
//----------------------------------------------------------------------------
{
  coord pos[6];
  Cache poscache = {0};
  Point lpoint;
   
  pos[0].x = gcp.center.x + gcp.radius + X0;
  pos[0].y = gcp.center.y + Y0;
  pos[0].z = gcp.center.z + Z0;  

  pos[1].x = gcp.center.x - gcp.radius + X0;
  pos[1].y = gcp.center.y + Y0;
  pos[1].z = gcp.center.z + Z0;  

  pos[2].x = gcp.center.x + X0;
  pos[2].y = gcp.center.y + gcp.radius + Y0;
  pos[2].z = gcp.center.z + Z0;
  
  pos[3].x = gcp.center.x + X0;
  pos[3].y = gcp.center.y - gcp.radius + Y0;
  pos[3].z = gcp.center.z + Z0;

  pos[4].x = gcp.center.x + X0;
  pos[4].y = gcp.center.y + Y0;  
  pos[4].z = gcp.center.z + gcp.radius + Z0;

  pos[5].x = gcp.center.x + X0;
  pos[5].y = gcp.center.y + Y0;
  pos[5].z = gcp.center.z - gcp.radius + Z0;
  
  *nb = 0;
  lpoint.level = -1;
  int ip = 0;
  
  /* We assume that around the sphere at least one point has to be at
     maximum level of refinement. */
  for (ip = 0; ip < 6; ip++) 
  {
    lpoint = locate( pos[ip].x, pos[ip].y, pos[ip].z );
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
      *nb += (int)( floor( pow( 3.809 * gcp.radius 
      		/ ( INTERBPCOEF * Delta ), 2.) ) );
      
    /** and finally, this thread destroys the cache. */    
    free( poscache.p );
  }

# if _MPI
    MPI_Barrier( MPI_COMM_WORLD );
    mpi_all_reduce( *nb, MPI_INT, MPI_MAX );
# endif
      
  if ( *nb == 0 )
    printf( "nboundary = 0: No boundary points !!!\n" );
}




/** Creates boundary points on the surface of the sphere */
//----------------------------------------------------------------------------
void create_FD_Boundary_Sphere( GeomParameter gcp, 
	SolidBodyBoundary * dlm_bd, const int nsphere, vector pshift ) 
//----------------------------------------------------------------------------
{
  coord pos;
  Point lpoint;
  
  coord shift = {0., 0., 0.};
  coord ori = {X0, Y0, Z0};
 
  /* Lay out points on a sphere with a Z oriented spiral scheme */
  /* This portion of code is recovered from Peligriff and adapted */
  /* More information of this method can be found at: */
  /* Saff, E. & Kuijlaars, A. Distributing many points on a sphere The
     Mathematical Intelligencer, 1997 */
  
  double spiral_spacing_correction = 0.;
  foreach_level(depth()) 
     spiral_spacing_correction = INTERBPCOEF * Delta / sqrt(3.) ; 
  
  double hydro_radius = gcp.radius;
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
      
    pos.x = hydro_radius * cos( phik ) * sin( thetak ) + gcp.center.x + X0;
    pos.y = hydro_radius * sin( phik ) * sin( thetak ) + gcp.center.y + Y0;
    pos.z = hydro_radius * cos( thetak ) + gcp.center.z + Z0;
        
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
    }
        
    Cache poscache = {0};
    lpoint = locate( pos.x, pos.y, pos.z );

    if ( lpoint.level > -1 ) 
    {	
      cache_append( &poscache, lpoint, 0 );

      foreach_cache(poscache)
	foreach_dimension()
	  pshift.x[] = shift.x;
      
      free( poscache.p );
    }
          
    dlm_bd->x[k] = pos.x;
    dlm_bd->y[k] = pos.y;
    dlm_bd->z[k] = pos.z;
  }

//  boundary((scalar*){pshift});
}




/** Finds cells lying inside the sphere */
//----------------------------------------------------------------------------
void create_FD_Interior_Sphere( particle * p, vector Index_lambda, 
	vector shift, scalar flag ) 
//----------------------------------------------------------------------------
{
  GeomParameter gci = p->g;
  Cache* fd = &(p->Interior);
  
  /** Create the cache for the interior points */
  foreach() 
  {
    if ((sq(x - (gci.center.x + X0)) + sq(y - (gci.center.y + Y0)) 
    	+ sq(z - (gci.center.z + Z0))) < sq(gci.radius)) 
    {
      cache_append (fd, point, 0);
      /* tagg cell with the number of the particle */
      if ((int)Index_lambda.y[] == -1)
	Index_lambda.y[] = p->pnum;
    }
    if (Period.x) {
      if ((sq(x - L0 - (gci.center.x + X0) ) + sq(y - (gci.center.y + Y0)) 
      	+ sq(z - (gci.center.z + Z0))) < sq(gci.radius)) 
      {
    	cache_append (fd, point, 0);
    	/* tagg cell with the number of the particle */
    	if ((int)Index_lambda.y[] == -1) 
    	  Index_lambda.y[] = p->pnum;
	if ((int)Index_lambda.y[] == p->pnum)
	  shift.x[] = L0;
      }
	if ((sq(x + L0 - (gci.center.x + X0) ) + sq(y - (gci.center.y + Y0)) 
	+ sq(z - (gci.center.z + Z0))) < sq(gci.radius)) 
      {
    	cache_append (fd, point, 0);
    	/* tagg cell with the number of the particle */
    	if ((int)Index_lambda.y[] == -1) 
    	  Index_lambda.y[] = p->pnum;
	if ((int)Index_lambda.y[] == p->pnum)
	  shift.x[] = -L0;
      }
    }
   
    if (Period.y) {
      if ((sq(x - (gci.center.x + X0)) + sq(y - L0 - (gci.center.y + Y0)) 
      	+ sq(z - (gci.center.z + Z0))) < sq(gci.radius)) 
      {
    	cache_append (fd, point, 0);
	/* tagg cell with the number of the particle */
	if ((int)Index_lambda.y[] == -1) {
	  Index_lambda.y[] = p->pnum;
	  if (flag[] < 1)
	    shift.y[] = L0;
	}
      }
      if ((sq(x - (gci.center.x + X0)) + sq(y + L0 - (gci.center.y + Y0)) 
      	+ sq(z - (gci.center.z + Z0))) < sq(gci.radius)) 
      {
    	cache_append (fd, point, 0);
	/* tagg cell with the number of the particle */
	if ((int)Index_lambda.y[] == -1) {
	  Index_lambda.y[] = p->pnum;
	  if (flag[] < 1)
	    shift.y[] = -L0;
	}
      }
    }

    if (Period.z) {
      if ((sq(x - (gci.center.x + X0)) + sq(y - (gci.center.y + Y0)) 
      	+ sq(z - L0 - (gci.center.z + Z0))) < sq(gci.radius)) 
      {
    	cache_append (fd, point, 0);
	/* tagg cell with the number of the particle */
	if ((int)Index_lambda.y[] == -1) {
	  Index_lambda.y[] = p->pnum;
	  if (flag[] < 1)
	    shift.z[] = L0;
	}
      }
      if ((sq(x - (gci.center.x + X0)) + sq(y - (gci.center.y + Y0)) 
      	+ sq(z + L0 - (gci.center.z + Z0))) < sq(gci.radius)) 
      {
    	cache_append (fd, point, 0);
	/* tagg cell with the number of the particle */
	if ((int)Index_lambda.y[] == -1) {
	  Index_lambda.y[] = p->pnum;
	  if (flag[] < 1)
	    shift.z[] = -L0;
	}
      }
    }    
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
