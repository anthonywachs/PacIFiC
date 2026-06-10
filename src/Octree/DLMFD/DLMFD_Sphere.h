/**
# Set of functions for a sphere
*/

# include "DLMFD_foreach_region_plusplus.h"


/** Tests whether a point lies inside the sphere */
//----------------------------------------------------------------------------
bool is_in_Sphere_geomtest( const double x, const double y, const double z,
	const coord center, const double radius )
//----------------------------------------------------------------------------
{
  return ( sqrt( sq( x - center.x ) + sq( y - center.y )
    	+ sq( z - center.z ) ) < radius );
}




/** Tests whether a point lies inside the sphere or any of its 
periodic clones */
//----------------------------------------------------------------------------
bool is_in_Sphere( const double x, const double y, const double z,
	GeomParameter const* gcp )
//----------------------------------------------------------------------------
{
  // Check if it is in the master rigid body
  bool status = is_in_Sphere_geomtest( x, y, z, gcp->center, gcp->radius );

  // Check if it is in any clone rigid body
  if ( gcp->nperclones && !status )
    for (int i = 0; i < gcp->nperclones && !status; i++)
      status = is_in_Sphere_geomtest( x, y, z, gcp->perclonecenters[i], 
      	gcp->radius );

  return ( status );
}




/** Flag boundary layer around the sphere */
//----------------------------------------------------------------------------
void flag_boundarylayer_Sphere( scalar flag_maxlevel, double const dcoef, 
	RigidBody const* p, AABB const* ld )
//----------------------------------------------------------------------------
{
  GeomParameter const* gcp = &(p->g); 
  AABB ExpBBox;
  double delta = L0 / (double)(1 << MAXLEVEL) ;
  double ext_radius = gcp->radius + dcoef * delta; 
  
  foreach_dimension()
  {
    ExpBBox.min.x = gcp->BBox.min.x - dcoef * delta;
    ExpBBox.max.x = gcp->BBox.max.x + dcoef * delta;
  } 
      
  // Loops over cells in the bounding box of the expanded sphere
  if ( intersect( ld, &ExpBBox ) )
    foreach_region_plus_plus( ExpBBox.min, ExpBBox.max ) 
      if ( is_leaf(cell) )
        if ( flag_maxlevel[] == 0. )
        {    
          if ( sqrt( sq( x - gcp->center.x ) + sq( y - gcp->center.y )
    		+ sq( z - gcp->center.z ) ) < ext_radius )
	    flag_maxlevel[] = 1.;
        }
	
  // Loops over cells in the bounding box of its clones
  AABB cloneBBox;
  coord shift;
  for (size_t i = 0; i < gcp->nperclones; i++)
  {
    foreach_dimension() shift.x = gcp->perclonecenters[i].x - gcp->center.x; 
    assign_shifted_BBox( &cloneBBox, &ExpBBox, shift );
    if ( intersect( ld, &cloneBBox ) )
      foreach_region_plus_plus(cloneBBox.min, cloneBBox.max) 
        if ( is_leaf(cell) )     
          if ( flag_maxlevel[] == 0. )
          {    
            if ( sqrt( sq( x - gcp->perclonecenters[i].x ) 
	    	+ sq( y - gcp->perclonecenters[i].y )
    		+ sq( z - gcp->perclonecenters[i].z ) ) < ext_radius )
	      flag_maxlevel[] = 1.;
          }
  }	
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




/** Creates boundary points and normal vectors of the reference sphere */
//----------------------------------------------------------------------------
void create_referenceRB_boundary_geomfeatures_Sphere( 
	GeomParameter const* gcp,
	RigidBodyBoundary* dlm_bd, const int nsphere )
//----------------------------------------------------------------------------
{
  coord pos;

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

    pos.x = hydro_radius * cos( phik ) * sin( thetak );
    pos.y = hydro_radius * sin( phik ) * sin( thetak );
    pos.z = hydro_radius * cos( thetak );

    foreach_dimension()
    { 
      dlm_bd->bp[k].x = pos.x;
      // We arbitrary set the norm of the normal vector to 0.25 * radius
      dlm_bd->outwardnormalvector[k].x = ( pos.x - gcp->center.x ) / 4.;
    }
  }
}



/** Finds cells lying inside the sphere */
//----------------------------------------------------------------------------
void create_FD_Interior_Sphere( RigidBody* p, vector Index,
	vector PeriodicRefCenter, AABB const* ld )
//----------------------------------------------------------------------------
{
  GeomParameter* gcp = &(p->g);
  Cache* fd = &(p->Interior);
  Point ppp;

  // Loops over cells in the bounding box of the sphere
  if ( intersect( ld, &(gcp->BBox) ) )
    foreach_region_plus_plus(gcp->BBox.min, gcp->BBox.max) 
      if ( is_leaf(cell) ) 
        if ( is_in_Sphere_geomtest( x, y, z, gcp->center, gcp->radius ) )
          if ( (int)Index.y[] == -1 )
          {
            foreach_dimension() PeriodicRefCenter.x[] = gcp->center.x;
	    ppp.i = point.i;
            ppp.j = point.j;
            ppp.k = point.k;			
            ppp.level = point.level;
	    cache_append( fd, ppp, 0 );
            Index.y[] = p->pnum;
          } 
	
  // Loops over cells in the bounding box of its clones
  AABB cloneBBox;
  coord shift;
  for (size_t i = 0; i < gcp->nperclones; i++)
  {
    foreach_dimension() shift.x = gcp->perclonecenters[i].x - gcp->center.x; 
    assign_shifted_BBox( &cloneBBox, &(gcp->BBox), shift );
    if ( intersect( ld, &cloneBBox ) )
      foreach_region_plus_plus(cloneBBox.min, cloneBBox.max) 
        if ( is_leaf(cell) )     
          if ( is_in_Sphere_geomtest( x, y, z, gcp->perclonecenters[i], 
		gcp->radius ) )
            if ( (int)Index.y[] == -1 )
            {
              foreach_dimension() 
	        PeriodicRefCenter.x[] = gcp->perclonecenters[i].x;
	      ppp.i = point.i;
              ppp.j = point.j;
              ppp.k = point.k;			
              ppp.level = point.level;
	      cache_append( fd, ppp, 0 );
              Index.y[] = p->pnum;
            }
  }   		     

  cache_shrink( fd );  
}




/** Reads geometric parameters of the sphere */
//----------------------------------------------------------------------------
void read_reference_Sphere( GeomParameter* gcp, const double RotMat[3][3] )
//----------------------------------------------------------------------------
{
  // We already have all parameters for the sphere but the input array of
  // characters contains again a "1", the gravity center coordinates and
  // a "0", hence we need to read five tokens but we do not do
  // anything with them
  for (size_t i=0;i<5;++i)
    strtok( NULL, " " );      
}




/** Update geometric parameters with the reference rigid body */
//----------------------------------------------------------------------------
void update_Sphere_from_RBRef( GeomParameter* gcp, RigidBody const* RBRef )
//----------------------------------------------------------------------------
{
  // Nothing to do  
}




/** Frees the geometric parameters of the sphere */
//----------------------------------------------------------------------------
void free_Sphere( GeomParameter* gcp )
//----------------------------------------------------------------------------
{
  // Nothing to do
}
