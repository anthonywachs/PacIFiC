/** 
# Set of functions for a 3D circular cylinder 
*/

# include "DLMFD_foreach_region_plusplus.h"


/** Tests whether a point lies inside the 3D circular cylinder */
//----------------------------------------------------------------------------
bool is_in_CircularCylinder3D_geomtest( const double x, const double y, 
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
  bool status = is_in_CircularCylinder3D_geomtest( x1, y1, z1, gcp );

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
      status = is_in_CircularCylinder3D_geomtest( x2, y2, z2, gcp );
    }

  return ( status );
}




/** Flag boundary layer around the 3D circular cylinder */
//----------------------------------------------------------------------------
void flag_boundarylayer_CircularCylinder3D( scalar flag_maxlevel, 
	double const dcoef, RigidBody const* p, AABB const* ld )
//----------------------------------------------------------------------------
{
  GeomParameter const* gcp = &(p->g); 
  AABB ExpBBox;
  double delta = L0 / (double)(1 << MAXLEVEL) ;
  double swellheight = 1. + 2. * dcoef * delta / gcp->cgp->height, x2, y2, z2;
  
  coord bottomToTopVec;  
  foreach_dimension()
    bottomToTopVec.x = swellheight * gcp->cgp->BottomToTopVec.x;
  coord bottomCenter;
  double norm = sqrt( sq( gcp->cgp->BottomToTopVec.x ) 
  	+ sq( gcp->cgp->BottomToTopVec.y ) 
	+ sq( gcp->cgp->BottomToTopVec.z ) );
  foreach_dimension()
    bottomCenter.x = gcp->cgp->BottomCenter.x - dcoef * delta
    	* gcp->cgp->BottomToTopVec.x / norm;
  double height = swellheight * gcp->cgp->height; 
  double radius = gcp->cgp->radius + dcoef * delta; 
       
  foreach_dimension()
  {
    ExpBBox.min.x = gcp->BBox.min.x - dcoef * delta;
    ExpBBox.max.x = gcp->BBox.max.x + dcoef * delta;
  } 
      
  // Loops over cells in the bounding box of the expanded 3D circular cylinder
  if ( intersect( ld, &ExpBBox ) )  
    foreach_region_plus_plus( ExpBBox.min, ExpBBox.max ) 
      if ( is_leaf(cell) )
        if ( flag_maxlevel[] == 0. )
        {    
          double dot = ( ( x - bottomCenter.x ) * bottomToTopVec.x
  		+ ( y - bottomCenter.y ) * bottomToTopVec.y
  		+ ( z - bottomCenter.z ) * bottomToTopVec.z ) / height ;

          if ( dot < height && dot > 0. )
            if ( ( x - bottomCenter.x ) * ( x - bottomCenter.x )
  		+ ( y - bottomCenter.y ) * ( y - bottomCenter.y )
  		+ ( z - bottomCenter.z ) * ( z - bottomCenter.z ) - sq( dot ) 
		< sq( radius ) )
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
            x2 = x - shift.x;
            y2 = y - shift.y;
            z2 = z - shift.z;        

            double dot = ( ( x2 - bottomCenter.x ) * bottomToTopVec.x
  		+ ( y2 - bottomCenter.y ) * bottomToTopVec.y
  		+ ( z2 - bottomCenter.z ) * bottomToTopVec.z ) / height ;

            if ( dot < height && dot > 0. )
              if ( ( x2 - bottomCenter.x ) * ( x2 - bottomCenter.x )
  		+ ( y2 - bottomCenter.y ) * ( y2 - bottomCenter.y )
  		+ ( z2 - bottomCenter.z ) * ( z2 - bottomCenter.z ) - sq( dot ) 
		< sq( radius ) )
              flag_maxlevel[] = 1.;
          }
  }	
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




/** Creates boundary points and normal vectors of the reference 3D circular 
cylinder */
//----------------------------------------------------------------------------
void create_referenceRB_boundary_geomfeatures_CircularCylinder3D( 
	GeomParameter const* gcp, RigidBodyBoundary* dlm_bd, const int m ) 
//----------------------------------------------------------------------------
{
  double delta = L0 / (double)(1 << MAXLEVEL) ;
  double spacing = INTERBPCOEF * delta;
  coord pos, unit_axial, n_cross_rad, top_normal, bottom_normal;
  int isb = 0;

  // Note: we arbitrary set the norm of the normal vector to 0.25 * radius
  
  foreach_dimension() 
    unit_axial.x = gcp->cgp->BottomToTopVec.x / gcp->cgp->height;
  n_cross_rad.x = unit_axial.y * gcp->cgp->RadialRefVec.z 
  	- unit_axial.z * gcp->cgp->RadialRefVec.y; 
  n_cross_rad.y = unit_axial.z * gcp->cgp->RadialRefVec.x 
  	- unit_axial.x * gcp->cgp->RadialRefVec.z;    
  n_cross_rad.z = unit_axial.x * gcp->cgp->RadialRefVec.y 
  	- unit_axial.y * gcp->cgp->RadialRefVec.x; 

  foreach_dimension() 
  {
    bottom_normal.x = ( gcp->cgp->BottomCenter.x - gcp->center.x ) 
    	* ( 0.5 * gcp->cgp->radius / gcp->cgp->height );
    top_normal.x = ( gcp->cgp->TopCenter.x - gcp->center.x ) 
    	* ( 0.5 * gcp->cgp->radius / gcp->cgp->height );
  }    
    
  // Cylinder height (diamond meshing)
  size_t npts_height = (size_t)( 2. * gcp->cgp->height / 
  	( sqrt(3.) * spacing ) ) + 1 ;
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
      local_angle = ( 2. * (double)(j) + bin ) * dangle ;      
      
      foreach_dimension() 
        pos.x = cos( local_angle ) * gcp->cgp->RadialRefVec.x
      		+ sin( local_angle ) * n_cross_rad.x
		+ (double)(i) * delta_height * unit_axial.x
		+ gcp->cgp->BottomCenter.x;                
		
      foreach_dimension() 
        dlm_bd->bp[isb].x = pos.x;
      
      if ( i == 0 )
        foreach_dimension() 
          dlm_bd->outwardnormalvector[isb].x = ( ( cos( local_angle ) 
	  	* gcp->cgp->RadialRefVec.x
      		+ sin( local_angle ) * n_cross_rad.x ) / 4.
		+ bottom_normal.x ) / sqrt(2.) ;
      else if ( i == npts_height - 1 )
        foreach_dimension() 
          dlm_bd->outwardnormalvector[isb].x = ( ( cos( local_angle ) 
	  	* gcp->cgp->RadialRefVec.x
      		+ sin( local_angle ) * n_cross_rad.x ) / 4.
		+ top_normal.x ) / sqrt(2.) ;
      else
        foreach_dimension() 
          dlm_bd->outwardnormalvector[isb].x = ( cos( local_angle ) 
	  	* gcp->cgp->RadialRefVec.x
      		+ sin( local_angle ) * n_cross_rad.x ) / 4.;
      isb++;
    }
  }

  // Bottom and top centers
  foreach_dimension() 
  {
    dlm_bd->bp[isb].x = gcp->cgp->BottomCenter.x;
    dlm_bd->outwardnormalvector[isb].x = bottom_normal.x;
  }
  isb++;
  		  
  foreach_dimension() 
  {
    dlm_bd->bp[isb].x = gcp->cgp->TopCenter.x;
    dlm_bd->outwardnormalvector[isb].x = top_normal.x;
  }    
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
      foreach_dimension() 
      {
        dlm_bd->bp[isb].x = pos.x;
	dlm_bd->outwardnormalvector[isb].x = bottom_normal.x;
      }
      isb++;      
      
      // Top disk
      foreach_dimension() 
        pos.x += gcp->cgp->BottomToTopVec.x;
      foreach_dimension() 
      {
        dlm_bd->bp[isb].x = pos.x;
        dlm_bd->outwardnormalvector[isb].x = top_normal.x;
      }	
      isb++;      
    }
  }
}




/** Finds cells lying inside the 3D circular cylinder */
//----------------------------------------------------------------------------
void create_FD_Interior_CircularCylinder3D( RigidBody* p, vector Index,
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
        if ( is_in_CircularCylinder3D_geomtest( x, y, z, gcp ) )
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
	  if ( is_in_CircularCylinder3D_geomtest( x2, y2, z2, gcp ) )
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




/** Reads geometric parameters of the 3D circular cylinder */
//----------------------------------------------------------------------------
void read_reference_CircularCylinder3D( GeomParameter* gcp, 
	const double RotMat[3][3] ) 
//----------------------------------------------------------------------------
{    
  char* token = NULL;
  coord v;

  // Read number of points, check that it is 3
  size_t np = 0;
  token = strtok(NULL, " " );
  sscanf( token, "%lu", &np );
  if ( np != 3 )
    printf ("Error in number of points in read_reference_CircularCylinder3D\n");
    
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
  
  
  // In case the reference rigid body was sent by the granular solver with 
  // a non zero center of mass and/or a non-zero identity angular position
  // we need to reset all corners to the neutral reference position
  // Bottom disk center
  // Translation    
  foreach_dimension() v.x = gcp->cgp->BottomCenter.x - gcp->center.x;
  // Rotation
  matTransposedCoordDotProduct( RotMat, v, &(gcp->cgp->BottomCenter) );
  
  // Radial reference vector 
  foreach_dimension() v.x = gcp->cgp->RadialRefVec.x; 
  // Rotation
  matTransposedCoordDotProduct( RotMat, v, &(gcp->cgp->RadialRefVec) );
  
  // Top disk center
  // Translation    
  foreach_dimension() v.x = gcp->cgp->TopCenter.x - gcp->center.x;
  // Rotation
  matTransposedCoordDotProduct( RotMat, v, &(gcp->cgp->TopCenter) );
  
  
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




/** Update geometric parameters with the reference rigid body */
//----------------------------------------------------------------------------
void update_CircularCylinder3D_from_RBRef( GeomParameter* gcp, 
	RigidBody const* RBRef, const double RotMat[3][3] ) 
//----------------------------------------------------------------------------
{        
  // Allocate the CylGeomParameter structure
  gcp->cgp = (CylGeomParameter*) malloc( sizeof(CylGeomParameter) );

  // Bottom disk center
  // Rotation
  matCoordDotProduct( RotMat, RBRef->g.cgp->BottomCenter, 
    	&(gcp->cgp->BottomCenter) );	
  // Translation
  foreach_dimension() gcp->cgp->BottomCenter.x += gcp->center.x;

  // Radial reference vector 
  // Rotation
  matCoordDotProduct( RotMat, RBRef->g.cgp->RadialRefVec, 
    	&(gcp->cgp->RadialRefVec) );

  // Top disk center
  // Rotation
  matCoordDotProduct( RotMat, RBRef->g.cgp->TopCenter, 
    	&(gcp->cgp->TopCenter) );	
  // Translation
  foreach_dimension() gcp->cgp->TopCenter.x += gcp->center.x; 	 

  // Compute the bottom to top vector
  foreach_dimension() 
    gcp->cgp->BottomToTopVec.x = gcp->cgp->TopCenter.x 
    		- gcp->cgp->BottomCenter.x;

  // Assign the bottom radius, top radius and the height
  gcp->cgp->radius = RBRef->g.cgp->radius;	
  gcp->cgp->height = RBRef->g.cgp->height;  
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
