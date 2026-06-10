/** 
# Set of functions for a cone 
*/

# include "DLMFD_foreach_region_plusplus.h"


/** Tests whether a point lies inside the cone */
//----------------------------------------------------------------------------
bool is_in_Cone_geomtest( const double x, const double y, 
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
	
  if ( proj >= 0. && proj <= gcp->tcgp->height )
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




/** Tests whether a point lies inside the cone or any of its periodic clones */
//----------------------------------------------------------------------------
bool is_in_Cone( const double x1, const double y1, 
	const double z1, GeomParameter const* gcp )
//----------------------------------------------------------------------------
{
  // Check if it is in the master rigid body
  bool status = is_in_Cone_geomtest( x1, y1, z1, gcp );

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
      status = is_in_Cone_geomtest( x2, y2, z2, gcp );
    }

  return ( status );
}




/** Flag boundary layer around the cone */
//----------------------------------------------------------------------------
void flag_boundarylayer_Cone( scalar flag_maxlevel, double const dcoef, 
	RigidBody const* p, AABB const* ld )
//----------------------------------------------------------------------------
{
  GeomParameter const* gcp = &(p->g); 
  AABB ExpBBox;
  coord vec;
  double delta = L0 / (double)(1 << MAXLEVEL), x2, y2, z2 ;
  double swell = 1. + dcoef * delta * ( 1. / gcp->tcgp->height
    	+ sqrt( sq( 1. / gcp->tcgp->height ) 
		+ sq( 1. / gcp->tcgp->BottomRadius) ) );
  
  coord bottomToTopVec;  
  foreach_dimension()
    bottomToTopVec.x = swell * gcp->tcgp->BottomToTopVec.x;
  coord bottomCenter;
  double vecnorm = sqrt( sq( gcp->tcgp->BottomToTopVec.x ) 
  	+ sq( gcp->tcgp->BottomToTopVec.y ) 
	+ sq( gcp->tcgp->BottomToTopVec.z ) );
  foreach_dimension()
    bottomCenter.x = gcp->tcgp->BottomCenter.x - dcoef * delta
    	* gcp->tcgp->BottomToTopVec.x / vecnorm;
  double bottomRadius = swell * gcp->tcgp->BottomRadius;
  double height = swell * gcp->tcgp->height;   
       
  foreach_dimension()
  {
    ExpBBox.min.x = gcp->BBox.min.x - dcoef * delta;
    ExpBBox.max.x = gcp->BBox.max.x + dcoef * delta;
  } 
      
  // Loops over cells in the bounding box of the expanded cone
  if ( intersect( ld, &ExpBBox ) )  
    foreach_region_plus_plus( ExpBBox.min, ExpBBox.max ) 
      if ( is_leaf(cell) )
        if ( flag_maxlevel[] == 0. )
        {    
          vec.x = x - bottomCenter.x;
          vec.y = y - bottomCenter.y;
          vec.z = z - bottomCenter.z; 
  
          double proj = vec.x * bottomToTopVec.x
		+ vec.y * bottomToTopVec.y 
		+ vec.z * bottomToTopVec.z; 
          proj /= height;
	
          if ( proj >= 0. && proj <= height )
          {
            foreach_dimension() 
              vec.x -= proj * bottomToTopVec.x / height;
      
            double local_radius = - bottomRadius * proj / height 
	  	+ bottomRadius;
	
            if ( sqrt( sq( vec.x ) + sq( vec.y ) + sq( vec.z ) ) 
	    		<= local_radius )
              flag_maxlevel[] = 1.;	  
          }
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
            x2 = x - shift.x;
            y2 = y - shift.y;
            z2 = z - shift.z;        

            vec.x = x2 - bottomCenter.x;
            vec.y = y2 - bottomCenter.y;
            vec.z = z2 - bottomCenter.z; 
  
            double proj = vec.x * bottomToTopVec.x
		+ vec.y * bottomToTopVec.y 
		+ vec.z * bottomToTopVec.z; 
            proj /= height;
	
            if ( proj >= 0. && proj <= height )
            {
              foreach_dimension() 
                vec.x -= proj * bottomToTopVec.x / height;
      
              double local_radius = - bottomRadius * proj / height 
	  	+ bottomRadius;
	
              if ( sqrt( sq( vec.x ) + sq( vec.y ) + sq( vec.z ) ) 
	    		<= local_radius )
                flag_maxlevel[] = 1.;	  
            }
          }
  }	
}




/** Computes the number of boundary points on the perimeter of the cone */
//----------------------------------------------------------------------------
void compute_nboundary_Cone( GeomParameter const* gcp, int* nb ) 
//----------------------------------------------------------------------------
{
  double delta = L0 / (double)(1 << MAXLEVEL) ;
  double spacing = INTERBPCOEF * delta;
  
  *nb = 2;
 
  size_t npts_radius = (size_t)( gcp->tcgp->BottomRadius / spacing ) + 1, n ;
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
    n = (size_t)( 2. * pi * local_radius / spacing ); 	
    *nb += n ? n : 1 ;
  }    
        
  if( *nb == 0 )
    printf( "nboundary = 0: No boundary points !!!\n" );
}




/** Creates boundary points and normal vectors of the reference cone */
//----------------------------------------------------------------------------
void create_referenceRB_boundary_geomfeatures_Cone( GeomParameter const* gcp,
	RigidBodyBoundary* dlm_bd, const int m ) 
//----------------------------------------------------------------------------
{
  double delta = L0 / (double)(1 << MAXLEVEL) ;
  double spacing = INTERBPCOEF * delta, local_angle, local_radius,
  	local_radius_ratio, delta_radius, inclined_height, bin, dangle, 
	delta_height, vecnorm, div = 4. * ( sqrt( 1 + sq( gcp->tcgp->BottomRadius )
		/ sq( gcp->tcgp->height ) ) - gcp->tcgp->BottomRadius
		/ gcp->tcgp->height );
  int isb = 0;
  coord pos, unit_axial, n_cross_rad, bottom_normal;
  size_t npts_local_radius, npts_radius, npts_height;

  // Note: we arbitrary set the vecnorm of the normal vector to 0.25 * radius

  foreach_dimension() 
    unit_axial.x = gcp->tcgp->BottomToTopVec.x / gcp->tcgp->height;  

  foreach_dimension() 
    bottom_normal.x = ( gcp->tcgp->BottomCenter.x - gcp->center.x ) 
    	* ( gcp->tcgp->BottomRadius / gcp->tcgp->height );    

  // Bottom center
  foreach_dimension() 
  {
    dlm_bd->bp[isb].x = gcp->tcgp->BottomCenter.x;
    dlm_bd->outwardnormalvector[isb].x = bottom_normal.x;
  }  
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
      { 
        pos.x = local_radius_ratio * ( 
			cos( local_angle ) * gcp->tcgp->BottomRadialRefVec.x
			+ sin( local_angle ) * n_cross_rad.x );     	 
        dlm_bd->bp[isb].x = pos.x + gcp->tcgp->BottomCenter.x;
      }
      
      if ( i == npts_radius - 1 )
      {
        vecnorm = 0.;
	foreach_dimension()
	{
	  dlm_bd->outwardnormalvector[isb].x = bottom_normal.x
	  	+ ( cos( local_angle ) * gcp->tcgp->BottomRadialRefVec.x
			+ sin( local_angle ) * n_cross_rad.x ) / div;
	  vecnorm += sq( dlm_bd->outwardnormalvector[isb].x );
        }
        vecnorm = sqrt( vecnorm );
        foreach_dimension() 
          dlm_bd->outwardnormalvector[isb].x *= 0.25 * gcp->tcgp->BottomRadius / vecnorm;
      }		     
      else
        foreach_dimension()
	  dlm_bd->outwardnormalvector[isb].x = bottom_normal.x;      
      
      isb++;          
    }
  }

  // Top center
  foreach_dimension() 
  {
    dlm_bd->bp[isb].x = gcp->tcgp->TopCenter.x;
    dlm_bd->outwardnormalvector[isb].x = ( gcp->tcgp->TopCenter.x - gcp->center.x )
    	* (  gcp->tcgp->BottomRadius / ( 3. * gcp->tcgp->height ) );
  }  
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
    npts_local_radius = npts_local_radius ? npts_local_radius : 1;
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
		
      vecnorm = 0.;
      foreach_dimension() 
      {
        dlm_bd->bp[isb].x = pos.x;
	dlm_bd->outwardnormalvector[isb].x = 
		( cos( local_angle ) * gcp->tcgp->BottomRadialRefVec.x  
      		+ sin( local_angle ) * n_cross_rad.x ) / gcp->tcgp->BottomRadius
		+ ( gcp->tcgp->BottomRadius / gcp->tcgp->height ) 
			* unit_axial.x ;
	vecnorm += sq( dlm_bd->outwardnormalvector[isb].x );
      }
      vecnorm = sqrt( vecnorm );
      foreach_dimension() 
        dlm_bd->outwardnormalvector[isb].x *= 0.25 * gcp->tcgp->BottomRadius / vecnorm;
            
      isb++;
    }         		    	
  }
}




/** Finds cells lying inside the cone */
//----------------------------------------------------------------------------
void create_FD_Interior_Cone( RigidBody* p, vector Index,
	vector PeriodicRefCenter, AABB const* ld )
//----------------------------------------------------------------------------
{
  GeomParameter const* gcp = &(p->g);  
  Cache* fd = &(p->Interior);
  Point ppp;
  double x2, y2, z2;  

  // Loops over cells in the bounding box of the sphere
  if ( intersect( ld, &(gcp->BBox) ) )
    foreach_region_plus_plus(gcp->BBox.min, gcp->BBox.max) 
      if ( is_leaf(cell) ) 
        if ( is_in_Cone_geomtest( x, y, z, gcp ) )
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
        {    
          x2 = x - shift.x;
          y2 = y - shift.y;
          z2 = z - shift.z;        
	  if ( is_in_Cone_geomtest( x2, y2, z2, gcp ) )
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
  }

  cache_shrink( fd );
}




/** Reads geometric parameters of the cone */
//----------------------------------------------------------------------------
void read_reference_Cone( GeomParameter* gcp, const double RotMat[3][3] ) 
//----------------------------------------------------------------------------
{    
  char* token = NULL;
  coord v;
  
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


  // In case the reference rigid body was sent by the granular solver with 
  // a non zero center of mass and/or a non-zero identity angular position
  // we need to reset all corners to the neutral reference position
  // Bottom disk center
  // Translation    
  foreach_dimension() v.x = gcp->tcgp->BottomCenter.x - gcp->center.x;
  // Rotation
  matTransposedCoordDotProduct( RotMat, v, &(gcp->tcgp->BottomCenter) );
  
  // Bottom radial reference vector 
  foreach_dimension() v.x = gcp->tcgp->BottomRadialRefVec.x; 
  // Rotation
  matTransposedCoordDotProduct( RotMat, v, &(gcp->tcgp->BottomRadialRefVec) );
  
  // Top disk center
  // Translation    
  foreach_dimension() v.x = gcp->tcgp->TopCenter.x - gcp->center.x;
  // Rotation
  matTransposedCoordDotProduct( RotMat, v, &(gcp->tcgp->TopCenter) );

  
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




/** Update geometric parameters with the reference rigid body */
//----------------------------------------------------------------------------
void update_Cone_from_RBRef( GeomParameter* gcp, 
	RigidBody const* RBRef, const double RotMat[3][3] ) 
//----------------------------------------------------------------------------
{        
  // Allocate the CylGeomParameter structure
  gcp->tcgp = (TruncConeGeomParameter*) malloc( 
  	sizeof(TruncConeGeomParameter) );

  // Bottom disk center
  // Rotation
  matCoordDotProduct( RotMat, RBRef->g.tcgp->BottomCenter, 
    	&(gcp->tcgp->BottomCenter) );	
  // Translation
  foreach_dimension() gcp->tcgp->BottomCenter.x += gcp->center.x;

  // Bottom radial reference vector 
  // Rotation
  matCoordDotProduct( RotMat, RBRef->g.tcgp->BottomRadialRefVec, 
    	&(gcp->tcgp->BottomRadialRefVec) );
	
  // Top disk center
  // Rotation
  matCoordDotProduct( RotMat, RBRef->g.tcgp->TopCenter, 
    	&(gcp->tcgp->TopCenter) );	
  // Translation
  foreach_dimension() gcp->tcgp->TopCenter.x += gcp->center.x;
  
  // Top radial reference vector 
  foreach_dimension() gcp->tcgp->TopRadialRefVec.x = 0.; 	 
  
  // Compute the bottom to top vector
  foreach_dimension() 
    gcp->tcgp->BottomToTopVec.x = gcp->tcgp->TopCenter.x 
    		- gcp->tcgp->BottomCenter.x;
    
  // Assign the bottom radius, top radius and the height
  gcp->tcgp->BottomRadius = RBRef->g.tcgp->BottomRadius;
  gcp->tcgp->TopRadius = 0.;	
  gcp->tcgp->height = RBRef->g.tcgp->height;			  
}




/** Frees the geometric parameters of the cone */
//----------------------------------------------------------------------------
void free_Cone( GeomParameter* gcp ) 
//----------------------------------------------------------------------------
{  
  // Free the CylGeomParameter structure
  free( gcp->tcgp );
  gcp->tcgp = NULL;
}
