/** 
# Set of functions for a 3D circular cylinder 
*/


# include "Polyhedron.h"

/** Tests whether a point lies inside the 3D circular cylinder */
//----------------------------------------------------------------------------
bool is_in_CircularCylinder3D_clone( const double x, const double y, 
	const double z, GeomParameter const* gcp )
//----------------------------------------------------------------------------
{
  bool isin = false; 
  double dot = ( ( x - gcp->cgp->BottomCenter.x ) * gcp->cgp->BottomToTopVec.x
  	+ ( y - gcp->cgp->BottomCenter.y ) * gcp->cgp->BottomToTopVec.y
  	+ ( z - gcp->cgp->BottomCenter.z ) * gcp->cgp->BottomToTopVec.z ) 
	/ gcp->cgp->height ;

  if ( dot < gcp->cgp->height && dot > 0. )
    if ( ( x - gcp->cgp->BottomCenter.x ) * ( x - gcp->cgp->BottomCenter.x )
  	+ ( y - gcp->cgp->BottomCenter.y ) * ( y - gcp->cgp->BottomCenter.y )
  	+ ( z - gcp->cgp->BottomCenter.z ) * ( z - gcp->cgp->BottomCenter.z ) 
	- sq( dot ) < sq( gcp->cgp->radius ) )
      isin = true ;
		  
  return ( isin );
}




/** Tests whether a point lies inside the 3D circular cylinder or any of its 
periodic clones */
//----------------------------------------------------------------------------
bool is_in_CircularCylinder3D( const double x1, const double y1, 
	const double z1, GeomParameter const* gcp )
//----------------------------------------------------------------------------
{
  // Check if it is in the master rigid body
  bool status = is_in_CircularCylinder3D_clone( x1, y1, z1, gcp );

  double x2, y2, z2;

  // Check if it is in any clone rigid body
  if ( gcp->nperclones && !status )
    for (int i = 0; i < gcp->nperclones && !status; i++)
    {
      GeomParameter clone = *gcp;
      clone.center = gcp->perclonecenters[i];
      x2 = x1 + gcp->center.x - clone.center.x;
      y2 = y1 + gcp->center.y - clone.center.y;
      z2 = z1 + gcp->center.z - clone.center.z;
      status = is_in_CircularCylinder3D_clone( x2, y2, z2, gcp );
    }

  return ( status );
}




/** Tests whether a point lies inside the 3D circular cylinder or any of its 
periodic clones and assign the proper center of mass coordinates associated to
this point */
//----------------------------------------------------------------------------
bool in_which_CircularCylinder3D( double x1, double y1, double z1, 
	GeomParameter const* gcp, vector* pPeriodicRefCenter, 
	const bool setPeriodicRefCenter )
//----------------------------------------------------------------------------
{
  // Check if it is in the master rigid body
  bool status = is_in_CircularCylinder3D_clone( x1, y1, z1, gcp );    
  if ( status && setPeriodicRefCenter )
    foreach_point( x1, y1, z1 )
      foreach_dimension()
        pPeriodicRefCenter->x[] = gcp->center.x;

  double x2, y2, z2;

  // Check if it is in any clone rigid body
  if ( gcp->nperclones && !status )
    for (int i = 0; i < gcp->nperclones && !status; i++) 
    {
      GeomParameter clone = *gcp;
      clone.center = gcp->perclonecenters[i];
      x2 = x1 + gcp->center.x - clone.center.x;
      y2 = y1 + gcp->center.y - clone.center.y;
      z2 = z1 + gcp->center.z - clone.center.z;
      status = is_in_CircularCylinder3D_clone( x2, y2, z2, gcp );
      if ( status && setPeriodicRefCenter )
        foreach_point( x1, y1, z1 )
          foreach_dimension()
            pPeriodicRefCenter->x[] = clone.center.x;
    }

  return ( status );
}




/** Computes the number of boundary points on the perimeter of the 3D circular 
cylinder */
//----------------------------------------------------------------------------
void compute_nboundary_CircularCylinder3D( GeomParameter const* gcp, int* nb ) 
//----------------------------------------------------------------------------
{
  double delta = L0 / (double)(1 << MAXLEVEL) ;
  double spacing = INTERBPCOEF * delta;
  int npts_height = (int)( 2. * gcp->cgp->height / ( sqrt(3.) * spacing ) ) 
  	+ 1 ;     
  int npts_local_radius = (int)( 2. * pi * gcp->cgp->radius / spacing );
    
  *nb = npts_height * npts_local_radius + 2;
  
  size_t npts_radius = (size_t)( gcp->cgp->radius / spacing ) + 1 ;
  double delta_radius = gcp->cgp->radius / ( (double)(npts_radius) - 1. ) ;
  for (size_t i=1;i<npts_radius-1;++i)
  {
    double local_radius = (double)(i) * delta_radius ;
    *nb += 2 * (size_t)( 2. * pi * local_radius / spacing ) ;    
  }  
        
  if( *nb == 0 )
    printf( "nboundary = 0: No boundary points !!!\n" );
}




/** Creates boundary points on the surface of the 3D circular cylinder */
//----------------------------------------------------------------------------
void create_FD_Boundary_CircularCylinder3D( GeomParameter const* gcp,
	RigidBodyBoundary* dlm_bd, const int nsphere, 
	vector* pPeriodicRefCenter, const bool setPeriodicRefCenter  ) 
//----------------------------------------------------------------------------
{
  double delta = L0 / (double)(1 << MAXLEVEL) ;
  double spacing = INTERBPCOEF * delta;
  coord pos, unit_axial, n_cross_rad;
  int isb = 0;
  
  foreach_dimension() 
    unit_axial.x = gcp->cgp->BottomToTopVec.x / gcp->cgp->height;
  n_cross_rad.x = unit_axial.y * gcp->cgp->RadialRefVec.z 
  	- unit_axial.z * gcp->cgp->RadialRefVec.y; 
  n_cross_rad.y = unit_axial.z * gcp->cgp->RadialRefVec.x 
  	- unit_axial.x * gcp->cgp->RadialRefVec.z;    
  n_cross_rad.z = unit_axial.x * gcp->cgp->RadialRefVec.y 
  	- unit_axial.y * gcp->cgp->RadialRefVec.x; 
    
  // Cylinder height (diamond meshing)
  size_t npts_height = (size_t)( 2. * gcp->cgp->height / ( sqrt(3.) * spacing ) ) 
  	+ 1 ;
  size_t npts_local_radius = (size_t)( 2. * pi * gcp->cgp->radius / spacing );
  double bin, local_angle, dangle = pi / (double)(npts_local_radius),
	delta_height = gcp->cgp->height / ( (double)(npts_height) - 1. ) ;  
  for (size_t i=0;i<npts_height;++i)
  {
    // odd or even
    if ( i % 2 == 0 ) bin = 0.;
    else bin = 1.;

    for (size_t j=0;j<npts_local_radius;++j)
    {
      local_angle = ( 2. * (double)(j) + pi * bin ) * dangle ;      
      
      foreach_dimension() 
        pos.x = cos( local_angle ) * gcp->cgp->RadialRefVec.x
      		+ sin( local_angle ) * n_cross_rad.x
		+ (double)(i) * delta_height * unit_axial.x
		+ gcp->cgp->BottomCenter.x;                

      periodic_correction( gcp, &pos, pPeriodicRefCenter, 
		setPeriodicRefCenter );
		
      foreach_dimension() dlm_bd->x[isb] = pos.x;
      isb++;
    }
  }

  // Bottom and top centers
  foreach_dimension() pos.x = gcp->cgp->BottomCenter.x;
  periodic_correction( gcp, &pos, pPeriodicRefCenter, 
		setPeriodicRefCenter );
  foreach_dimension() dlm_bd->x[isb] = pos.x;
  isb++;
  		  
  foreach_dimension() pos.x = gcp->cgp->TopCenter.x;
  periodic_correction( gcp, &pos, pPeriodicRefCenter, 
		setPeriodicRefCenter );
  foreach_dimension() dlm_bd->x[isb] = pos.x;
  isb++;  

  // Bottom and top disks in concentric circles 
  size_t npts_radius = (size_t)( gcp->cgp->radius / spacing ) + 1 ;
  double delta_radius = gcp->cgp->radius / ( (double)(npts_radius) - 1. ) ;
  for (size_t i=1;i<npts_radius-1;++i)
  {
    double local_radius = (double)(i) * delta_radius ;
    double local_radius_ratio = local_radius / gcp->cgp->radius ;
    npts_local_radius = (size_t)( 2. * pi * local_radius / spacing ) ;
      
    for (size_t j=0;j<npts_local_radius;++j)
    {      
      local_angle = 2. * pi * (double)(j) / (double)(npts_local_radius) ;
      
      foreach_dimension() 
        pos.x = local_radius_ratio * ( 
			cos( local_angle ) * gcp->cgp->RadialRefVec.x
			+ sin( local_angle ) * n_cross_rad.x );     	
      
      // Bottom disk
      foreach_dimension() 
        pos.x += gcp->cgp->BottomCenter.x;
      periodic_correction( gcp, &pos, pPeriodicRefCenter, 
		setPeriodicRefCenter );
      foreach_dimension() dlm_bd->x[isb] = pos.x;
      isb++;      
      
      // Top disk
      foreach_dimension() 
        pos.x += gcp->cgp->BottomToTopVec.x;
      periodic_correction( gcp, &pos, pPeriodicRefCenter, 
		setPeriodicRefCenter );
      foreach_dimension() dlm_bd->x[isb] = pos.x;
      isb++;      
    }
  }
  
  if ( setPeriodicRefCenter ) synchronize((scalar*){pPeriodicRefCenter->x,
  	pPeriodicRefCenter->y, pPeriodicRefCenter->z});      
}




/** Finds cells lying inside the 3D circular cylinder */
//----------------------------------------------------------------------------
void create_FD_Interior_CircularCylinder3D( RigidBody* p, vector Index,
	vector PeriodicRefCenter )
//----------------------------------------------------------------------------
{
  GeomParameter const* gcp = &(p->g);  
  Cache* fd = &(p->Interior);
  Point ppp;

  /** Create the cache of the interior points of the polyhedron */
  foreach(serial)
    if ( in_which_CircularCylinder3D( x, y, z, gcp, &PeriodicRefCenter, true ) )
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




/** Reads geometric parameters of the 3D circular cylinder */
//----------------------------------------------------------------------------
void update_CircularCylinder3D( GeomParameter* gcp ) 
//----------------------------------------------------------------------------
{    
  char* token = NULL;

  // Read number of points, check that it is 3
  size_t np = 0;
  token = strtok(NULL, " " );
  sscanf( token, "%lu", &np );
  if ( np != 3 )
    printf ("Error in number of points in update_CircularCylinder3D\n");
    
  // Allocate the CylGeomParameter structure
  gcp->cgp = (CylGeomParameter*) malloc( sizeof(CylGeomParameter) );

  // Read the bottom disk center
  foreach_dimension()
  {
    token = strtok(NULL, " " );
    sscanf( token, "%lf", &(gcp->cgp->BottomCenter.x) );
  }
  
  // Read the point to compute the radial reference vector
  foreach_dimension()
  {
    token = strtok(NULL, " " );
    sscanf( token, "%lf", &(gcp->cgp->TopCenter.x) );
    gcp->cgp->RadialRefVec.x = gcp->cgp->TopCenter.x - gcp->cgp->BottomCenter.x;
  }

  // Read the top disk center
  foreach_dimension()
  {
    token = strtok(NULL, " " );
    sscanf( token, "%lf", &(gcp->cgp->TopCenter.x) );
  }
  
  // We already have all parameters for the 3D circular cylinder but the input 
  // array of characters contains an additional "0", hence we need to read one 
  // token but we do not do anything with it
  strtok( NULL, " " );
  
  // Compute the bottom to top vector
  foreach_dimension() 
    gcp->cgp->BottomToTopVec.x = gcp->cgp->TopCenter.x 
    		- gcp->cgp->BottomCenter.x;
    
  // Compute the radius and the height
  gcp->cgp->radius = sqrt( sq( gcp->cgp->RadialRefVec.x ) 
  	+ sq( gcp->cgp->RadialRefVec.y )
  	+ sq( gcp->cgp->RadialRefVec.z ) );
  gcp->cgp->height = sqrt( sq( gcp->cgp->BottomToTopVec.x ) 
  	+ sq( gcp->cgp->BottomToTopVec.y )
  	+ sq( gcp->cgp->BottomToTopVec.z ) );	  
}




/** Frees the geometric parameters of the 3D circular cylinder */
//----------------------------------------------------------------------------
void free_CircularCylinder3D( GeomParameter* gcp ) 
//----------------------------------------------------------------------------
{  
  // Free the CylGeomParameter structure
  free( gcp->cgp );
  gcp->cgp = NULL;
}
