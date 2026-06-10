/** 
# Set of functions for a 2D circular cylinder 
*/

# include "DLMFD_foreach_region_plusplus.h"


/** Tests whether a point lies inside the 2D circular cylinder */
//----------------------------------------------------------------------------
bool is_in_CircularCylinder2D_geomtest( const double x, const double y, 
	const coord center, const double radius )
//----------------------------------------------------------------------------
{
  return ( sqrt( sq( x - center.x ) + sq( y - center.y ) ) < radius );
}




/** Tests whether a point lies inside the 2D circular cylinder or any of its 
periodic clones */
//----------------------------------------------------------------------------
bool is_in_CircularCylinder2D( const double x, const double y, 
	GeomParameter const* gcp )
//----------------------------------------------------------------------------
{
  // Check if it is in the master rigid body
  bool status = is_in_CircularCylinder2D_geomtest( x, y, gcp->center, 
  	gcp->radius );

  // Check if it is in any clone rigid body
  if ( gcp->nperclones && !status )
    for (int i = 0; i < gcp->nperclones && !status; i++)
      status = is_in_CircularCylinder2D_geomtest( x, y, gcp->perclonecenters[i], 
      	gcp->radius );

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




/** Flag boundary layer around the  2D circular cylinder */
//----------------------------------------------------------------------------
void flag_boundarylayer_CircularCylinder2D( scalar flag_maxlevel, 
	double const dcoef, RigidBody const* p, AABB const* ld )
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

  // Loops over cells in the bounding box of the 2D circular cylinder      
  if ( intersect( ld, &ExpBBox ) )
    foreach_region_plus_plus( ExpBBox.min, ExpBBox.max ) 
      if ( is_leaf(cell) )
        if ( flag_maxlevel[] == 0. )
        {    
          if ( sqrt( sq( x - gcp->center.x ) + sq( y - gcp->center.y ) ) 
		< ext_radius )
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
	    	+ sq( y - gcp->perclonecenters[i].y ) ) < ext_radius )
	      flag_maxlevel[] = 1.;
          }
  }	
}




/** Creates boundary points and normal vectors of the 2D circular cylinder */
//----------------------------------------------------------------------------
void create_referenceRB_boundary_geomfeatures_CircularCylinder2D( 
	GeomParameter const* gcp,
	RigidBodyBoundary* dlm_bd, const int nsphere ) 
//----------------------------------------------------------------------------
{
  int k, m = nsphere;  
  double thetak, radius = gcp->radius;
  
  for (k = 0; k < m; k++) 
  {
    thetak = (double)(k) * 2. * pi / (double)(m); 
    dlm_bd->bp[k].x = radius * cos( thetak );
    dlm_bd->bp[k].y = radius * sin( thetak );
    dlm_bd->bp[k].z = 0.;
    // We arbitrary set the norm of the normal vector to 0.25 * radius
    dlm_bd->outwardnormalvector[k].x = ( dlm_bd->bp[k].x - gcp->center.x ) / 4.;
    dlm_bd->outwardnormalvector[k].y = ( dlm_bd->bp[k].y - gcp->center.y ) / 4.;
    dlm_bd->outwardnormalvector[k].z = 0.;         
  }
}




/** Finds cells lying inside the 2D circular cylinder */
//----------------------------------------------------------------------------
void create_FD_Interior_CircularCylinder2D( RigidBody* p, vector Index,
	vector PeriodicRefCenter, AABB const* ld )
//----------------------------------------------------------------------------
{
  GeomParameter const* gcp = &(p->g);
  Cache* fd = &(p->Interior);
  Point ppp;

  // Loops over cells in the bounding box of the 2D circular cylinder
  if ( intersect( ld, &(gcp->BBox) ) )  
    foreach_region_plus_plus(gcp->BBox.min, gcp->BBox.max) 
      if ( is_leaf(cell) ) 
        if ( is_in_CircularCylinder2D_geomtest( x, y, gcp->center, 
		gcp->radius ) )
          if ( (int)Index.y[] == -1 )
          {
            foreach_dimension() PeriodicRefCenter.x[] = gcp->center.x;
	    ppp.i = point.i;
            ppp.j = point.j;			
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
          if ( is_in_CircularCylinder2D_geomtest( x, y, 
	  	gcp->perclonecenters[i], gcp->radius ) )
            if ( (int)Index.y[] == -1 )
            {
              foreach_dimension() 
	        PeriodicRefCenter.x[] = gcp->perclonecenters[i].x;
	      ppp.i = point.i;
              ppp.j = point.j;		
              ppp.level = point.level;
	      cache_append( fd, ppp, 0 );
              Index.y[] = p->pnum;
            }
  }

  cache_shrink( fd );
}




/** Reads geometric parameters of the 2D circular cylinder */
//----------------------------------------------------------------------------
void read_reference_CircularCylinder2D( GeomParameter* gcp, 
	const double RotMat[3][3] ) 
//----------------------------------------------------------------------------
{    
  // TO DO
}




/** Update geometric parameters with the reference rigid body */
//----------------------------------------------------------------------------
void update_CircularCylinder2D_from_RBRef( GeomParameter* gcp, 
	RigidBody const* RBRef )
//----------------------------------------------------------------------------
{
  // Nothing to do  
}



/** Frees the geometric parameters of the 2D circular cylinder */
//----------------------------------------------------------------------------
void free_CircularCylinder2D( GeomParameter* gcp ) 
//----------------------------------------------------------------------------
{  
  // Nothing to do
}
