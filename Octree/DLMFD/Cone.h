/** 
# Set of functions for a truncated cone 
*/


# include "Polyhedron.h"

/** Tests whether a point lies inside the truncated cone */
//----------------------------------------------------------------------------
bool is_in_Cone_clone( const double x, const double y, 
	const double z, GeomParameter const* gcp )
//----------------------------------------------------------------------------
{
  bool isin = false; 

  coord vec;
  vec.x = x - gcp->tcgp->BottomCenter.x;
  vec.y = y - gcp->tcgp->BottomCenter.y;
  vec.z = z - gcp->tcgp->BottomCenter.z; 
  
  double proj = vec.x * gcp->tcgp->BottomToTopVec.x
	+ vec.y * gcp->tcgp->BottomToTopVec.y 
	+ vec.z * gcp->tcgp->BottomToTopVec.z; 
  proj /= gcp->tcgp->height;
	
  if ( proj > 0. && proj < gcp->tcgp->height )
  {
    foreach_dimension() 
      vec.x -= proj * gcp->tcgp->BottomToTopVec.x / gcp->tcgp->height;
      
    double local_radius = - gcp->tcgp->BottomRadius
    	* proj / gcp->tcgp->height + gcp->tcgp->BottomRadius;
	
    if ( sqrt( sq( vec.x ) + sq( vec.y ) + sq( vec.z ) ) <= local_radius )
      isin = true;	  
  }	 
		  
  return ( isin );
}




/** Tests whether a point lies inside the truncated cone or any of its 
periodic clones */
//----------------------------------------------------------------------------
bool is_in_Cone( const double x1, const double y1, 
	const double z1, GeomParameter const* gcp )
//----------------------------------------------------------------------------
{
  // Check if it is in the master rigid body
  bool status = is_in_Cone_clone( x1, y1, z1, gcp );

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
      status = is_in_Cone_clone( x2, y2, z2, gcp );
    }

  return ( status );
}




/** Tests whether a point lies inside the truncated cone or any of its 
periodic clones and assign the proper center of mass coordinates associated to
this point */
//----------------------------------------------------------------------------
bool in_which_Cone( double x1, double y1, double z1, 
	GeomParameter const* gcp, vector* pPeriodicRefCenter, 
	const bool setPeriodicRefCenter )
//----------------------------------------------------------------------------
{
  // Check if it is in the master rigid body
  bool status = is_in_Cone_clone( x1, y1, z1, gcp );    
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
      status = is_in_Cone_clone( x2, y2, z2, gcp );
      if ( status && setPeriodicRefCenter )
        foreach_point( x1, y1, z1 )
          foreach_dimension()
            pPeriodicRefCenter->x[] = clone.center.x;
    }

  return ( status );
}




/** Computes the number of boundary points on the perimeter of the truncated 
cone */
//----------------------------------------------------------------------------
void compute_nboundary_Cone( GeomParameter const* gcp, int* nb ) 
//----------------------------------------------------------------------------
{
  double delta = L0 / (double)(1 << MAXLEVEL) ;
  double spacing = INTERBPCOEF * delta;
  
  *nb = 2;
 
  size_t npts_radius = (size_t)( gcp->tcgp->BottomRadius / spacing ) + 1 ;
  double delta_radius = gcp->tcgp->BottomRadius 
  	/ ( (double)(npts_radius) - 1. ) ;
  for (size_t i=1;i<npts_radius;++i)
  {
    double local_radius = (double)(i) * delta_radius ;
    *nb += (size_t)( 2. * pi * local_radius / spacing ) ;    
  }
  
  coord pos;
  foreach_dimension()
    pos.x = gcp->tcgp->BottomToTopVec.x - gcp->tcgp->BottomRadialRefVec.x;
  double inclined_height = sqrt( sq( pos.x ) + sq( pos.y ) + sq( pos.z ) );
  size_t npts_height = (size_t)( 2. * inclined_height / ( sqrt(3.) * spacing ) )
  	+ 1;
  double delta_height = gcp->tcgp->height / ( (double)(npts_height) - 1. ) ;	
  for (size_t i=1;i<npts_height-1;++i)
  {
    double local_radius = - gcp->tcgp->BottomRadius
    	* (double)(i) * delta_height / gcp->tcgp->height 
	+ gcp->tcgp->BottomRadius;
    *nb += (size_t)( 2. * pi * local_radius / spacing ) ;
  }    
        
  if( *nb == 0 )
    printf( "nboundary = 0: No boundary points !!!\n" );
}




/** Creates boundary points on the surface of the truncated cone */
//----------------------------------------------------------------------------
void create_FD_Boundary_Cone( GeomParameter const* gcp,
	RigidBodyBoundary* dlm_bd, const int nsphere, 
	vector* pPeriodicRefCenter, const bool setPeriodicRefCenter  ) 
//----------------------------------------------------------------------------
{
  double delta = L0 / (double)(1 << MAXLEVEL) ;
  double spacing = INTERBPCOEF * delta, local_angle, local_radius,
  	local_radius_ratio, delta_radius, inclined_height, bin, dangle, 
	delta_height;
  int isb = 0;
  coord pos, unit_axial, n_cross_rad;
  size_t npts_local_radius, npts_radius, npts_height;
  
  foreach_dimension() 
    unit_axial.x = gcp->tcgp->BottomToTopVec.x / gcp->tcgp->height;  
  
  // Bottom center
  foreach_dimension() pos.x = gcp->tcgp->BottomCenter.x;
  periodic_correction( gcp, &pos, pPeriodicRefCenter, 
		setPeriodicRefCenter );
  foreach_dimension() dlm_bd->x[isb] = pos.x;
  isb++;
  
  // Bottom disk in concentric circles
  n_cross_rad.x = unit_axial.y * gcp->tcgp->BottomRadialRefVec.z 
  	- unit_axial.z * gcp->tcgp->BottomRadialRefVec.y; 
  n_cross_rad.y = unit_axial.z * gcp->tcgp->BottomRadialRefVec.x 
  	- unit_axial.x * gcp->tcgp->BottomRadialRefVec.z;    
  n_cross_rad.z = unit_axial.x * gcp->tcgp->BottomRadialRefVec.y 
  	- unit_axial.y * gcp->tcgp->BottomRadialRefVec.x;   
  npts_radius = (size_t)( gcp->tcgp->BottomRadius / spacing ) + 1 ;
  delta_radius = gcp->tcgp->BottomRadius 
  	/ ( (double)(npts_radius) - 1. ) ;
  for (size_t i=1;i<npts_radius;++i)
  {
    local_radius = (double)(i) * delta_radius ;
    local_radius_ratio = local_radius / gcp->tcgp->BottomRadius ;
    npts_local_radius = (size_t)( 2. * pi * local_radius / spacing ) ;
      
    for (size_t j=0;j<npts_local_radius;++j)
    {      
      local_angle = 2. * pi * (double)(j) / (double)(npts_local_radius) ;
      
      foreach_dimension() 
        pos.x = local_radius_ratio * ( 
			cos( local_angle ) * gcp->tcgp->BottomRadialRefVec.x
			+ sin( local_angle ) * n_cross_rad.x );     	
      
      foreach_dimension() 
        pos.x += gcp->tcgp->BottomCenter.x;
      periodic_correction( gcp, &pos, pPeriodicRefCenter, 
		setPeriodicRefCenter );
      foreach_dimension() dlm_bd->x[isb] = pos.x;
      isb++;          
    }
  }
  
  // Top center
  foreach_dimension() pos.x = gcp->tcgp->TopCenter.x;
  periodic_correction( gcp, &pos, pPeriodicRefCenter, 
		setPeriodicRefCenter );
  foreach_dimension() dlm_bd->x[isb] = pos.x;
  isb++;

  // Lateral surface
  foreach_dimension()
    pos.x = gcp->tcgp->BottomToTopVec.x - gcp->tcgp->BottomRadialRefVec.x;
  inclined_height = sqrt( sq( pos.x ) + sq( pos.y ) + sq( pos.z ) );
  npts_height = (size_t)( 2. * inclined_height / ( sqrt(3.) * spacing ) )
  	+ 1;
  delta_height = gcp->tcgp->height / ( (double)(npts_height) - 1. ) ;	
  for (size_t i=1;i<npts_height-1;++i)
  {
    local_radius = - gcp->tcgp->BottomRadius
    	* (double)(i) * delta_height / gcp->tcgp->height 
	+ gcp->tcgp->BottomRadius;
    npts_local_radius = (size_t)( 2. * pi * local_radius / spacing ) ;
    dangle = pi / (double)(npts_local_radius);
     
    // odd or even
    if ( i % 2 == 0 ) bin = 0.;
    else bin = 1.;
    
    for (size_t j=0;j<npts_local_radius;++j)
    {
      local_angle = ( 2. * (double)(j) + bin ) * dangle ;      
      
      foreach_dimension() 
        pos.x = cos( local_angle ) * gcp->tcgp->BottomRadialRefVec.x 
		* local_radius / gcp->tcgp->BottomRadius 
      		+ sin( local_angle ) * n_cross_rad.x * local_radius 
			/ gcp->tcgp->BottomRadius
		+ (double)(i) * delta_height * unit_axial.x
		+ gcp->tcgp->BottomCenter.x;                

      periodic_correction( gcp, &pos, pPeriodicRefCenter, 
		setPeriodicRefCenter );
		
      foreach_dimension() dlm_bd->x[isb] = pos.x;
      isb++;
    }         		    	
  }           
  
  if ( setPeriodicRefCenter ) synchronize((scalar*){pPeriodicRefCenter->x,
  	pPeriodicRefCenter->y, pPeriodicRefCenter->z});      
}




/** Finds cells lying inside the truncated cone */
//----------------------------------------------------------------------------
void create_FD_Interior_Cone( RigidBody* p, vector Index,
	vector PeriodicRefCenter )
//----------------------------------------------------------------------------
{
  GeomParameter const* gcp = &(p->g);  
  Cache* fd = &(p->Interior);
  Point ppp;

  /** Create the cache of the interior points of the polyhedron */
  foreach(serial)
    if ( in_which_Cone( x, y, z, gcp, &PeriodicRefCenter, true ) )
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




/** Reads geometric parameters of the truncated cone */
//----------------------------------------------------------------------------
void update_Cone( GeomParameter* gcp ) 
//----------------------------------------------------------------------------
{    
  char* token = NULL;

  // Read number of points, check that it is 4
  size_t np = 0;
  token = strtok(NULL, " " );
  sscanf( token, "%lu", &np );
  if ( np != 3 )
    printf ("Error in number of points in update_Cone\n");
    
  // Allocate the CylGeomParameter structure
  gcp->tcgp = (TruncConeGeomParameter*) malloc( 
  	sizeof(TruncConeGeomParameter) );

  // Read the bottom disk center
  foreach_dimension()
  {
    token = strtok(NULL, " " );
    sscanf( token, "%lf", &(gcp->tcgp->BottomCenter.x) );
  }
  
  // Read the point to compute the bottom radial reference vector
  foreach_dimension()
  {
    token = strtok(NULL, " " );
    sscanf( token, "%lf", &(gcp->tcgp->BottomRadialRefVec.x) );
    gcp->tcgp->BottomRadialRefVec.x -= gcp->tcgp->BottomCenter.x;
  }

  // For a cone, the top center corresponds to the tip
  foreach_dimension()
  {
    token = strtok(NULL, " " );
    sscanf( token, "%lf", &(gcp->tcgp->TopCenter.x) );
  }
  
  // The top radial reference vector is a zero vector for a cone
  foreach_dimension()
    gcp->tcgp->TopRadialRefVec.x = 0.; 
  
  // We already have all parameters for the 3D circular cylinder but the input 
  // array of characters contains an additional "0", hence we need to read one 
  // token but we do not do anything with it
  strtok( NULL, " " );
  
  // Compute the bottom to top vector
  foreach_dimension() 
    gcp->tcgp->BottomToTopVec.x = gcp->tcgp->TopCenter.x 
    		- gcp->tcgp->BottomCenter.x;
    
  // Compute the bottom radius and the height
  gcp->tcgp->BottomRadius = sqrt( sq( gcp->tcgp->BottomRadialRefVec.x ) 
  	+ sq( gcp->tcgp->BottomRadialRefVec.y )
  	+ sq( gcp->tcgp->BottomRadialRefVec.z ) );
  gcp->tcgp->height = sqrt( sq( gcp->tcgp->BottomToTopVec.x ) 
  	+ sq( gcp->tcgp->BottomToTopVec.y )
  	+ sq( gcp->tcgp->BottomToTopVec.z ) );
	
  // The top radius is zero for a cone
  gcp->tcgp->TopRadius = 0.;				  
}




/** Frees the geometric parameters of the truncated cone */
//----------------------------------------------------------------------------
void free_Cone( GeomParameter* gcp ) 
//----------------------------------------------------------------------------
{  
  // Free the CylGeomParameter structure
  free( gcp->tcgp );
  gcp->tcgp = NULL;
}
