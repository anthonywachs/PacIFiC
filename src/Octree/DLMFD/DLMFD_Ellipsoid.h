/** 
# Set of functions for an ellipsoid
*/

# include "DLMFD_foreach_region_plusplus.h"


/** Tests whether a point lies inside the ellipsoid */
//----------------------------------------------------------------------------
bool is_in_Ellipsoid_geomtest( const double x, const double y, 
	const double z, RigidBody const* p )
//----------------------------------------------------------------------------
{
  bool isin = false; 
  coord v, vr;

  // Translation   
  v.x = x - p->g.center.x;
  v.y = y - p->g.center.y;  
  v.z = z - p->g.center.z;  
  
  // Rotation
  matTransposedCoordDotProduct( p->RotMat, v, &vr );
  
  // Compute potential
  isin = sq( vr.x / p->g.elgp->a ) + sq( vr.y / p->g.elgp->b )
  	+ sq( vr.z / p->g.elgp->c ) - 1. <= 0.;
  		  
  return ( isin );
}




/** Tests whether a point lies inside the ellipsoid or any of its 
periodic clones */
//----------------------------------------------------------------------------
bool is_in_Ellipsoid( const double x1, const double y1, 
	const double z1, RigidBody const* p )
//----------------------------------------------------------------------------
{
  // Check if it is in the master rigid body
  bool status = is_in_Ellipsoid_geomtest( x1, y1, z1, p );

  double x2, y2, z2;

  // Check if it is in any clone rigid body
  if ( p->g.nperclones && !status )
    for (int i = 0; i < p->g.nperclones && !status; i++)
    {
      x2 = x1 + p->g.center.x - p->g.perclonecenters[i].x;
      y2 = y1 + p->g.center.y - p->g.perclonecenters[i].y;
      z2 = z1 + p->g.center.z - p->g.perclonecenters[i].z;
      status = is_in_Ellipsoid_geomtest( x2, y2, z2, p );
    }

  return ( status );
}




/** Flag boundary layer around the ellipsoid */
//----------------------------------------------------------------------------
void flag_boundarylayer_Ellipsoid( scalar flag_maxlevel, double const dcoef, 
	RigidBody const* p, AABB const* ld )
//----------------------------------------------------------------------------
{
  GeomParameter const* gcp = &(p->g); 
  AABB ExpBBox;
  coord v, vr;
  double delta = L0 / (double)(1 << MAXLEVEL), x2, y2, z2 ;
  double a = gcp->elgp->a + dcoef * delta,
  	b = gcp->elgp->b + dcoef * delta,
	c = gcp->elgp->c + dcoef * delta; 
  
  foreach_dimension()
  {
    ExpBBox.min.x = gcp->BBox.min.x - dcoef * delta;
    ExpBBox.max.x = gcp->BBox.max.x + dcoef * delta;
  } 
      
  // Loops over cells in the bounding box of the expanded ellipsoid
  if ( intersect( ld, &ExpBBox ) )  
    foreach_region_plus_plus( ExpBBox.min, ExpBBox.max ) 
      if ( is_leaf(cell) )
        if ( flag_maxlevel[] == 0. )
        {    
          v.x = x - gcp->center.x;
	  v.y = y - gcp->center.y;  
	  v.z = z - gcp->center.z;  
  
          matTransposedCoordDotProduct( p->RotMat, v, &vr );
  
          if ( sq( vr.x / a ) + sq( vr.y / b ) + sq( vr.z / c ) - 1. <= 0. )
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

            v.x = x2 - gcp->center.x;
	    v.y = y2 - gcp->center.y;  
	    v.z = z2 - gcp->center.z;  
  
            matTransposedCoordDotProduct( p->RotMat, v, &vr );
  
            if ( sq( vr.x / a ) + sq( vr.y / b ) + sq( vr.z / c ) - 1. <= 0. )
	      flag_maxlevel[] = 1.;
          }
  }	
}




/** Computes the distance between 2 points (xip1,y,0) and (xi,y,0) that lie
on the ellipsoid surface */
//----------------------------------------------------------------------------
double fx( double xip1, double xi, double a, double b, double l )
//----------------------------------------------------------------------------
{
  return( sq( xip1 - xi ) + sq( b ) * sq( sqrt( 1. - sq( xip1 / a )  ) 
  	- sqrt( 1. - sq( xi / a ) ) ) - sq( l ) ); 
}




/** Distributes points over the x-axis of the ellipsoid surface that points 
that lie on the surface at z=0 are equidistant by roughly spacing */
//----------------------------------------------------------------------------
double* xaxis_points_distribution( GeomParameter const* gcp, size_t* totalnx,
	double* spacing )
//----------------------------------------------------------------------------	
{
  *totalnx = 0;  
  double delta = L0 / (double)(1 << MAXLEVEL) ;
  *spacing = INTERBPCOEF * delta;
  double refs = *spacing;  
  double lb, ub, mid, sol = - 1.;
  double a = gcp->elgp->a, b = gcp->elgp->b;
  size_t maxnx = (size_t)( 2. * ( a + b ) / delta ), ninterior = 0, i;
  double* vx = (double*) calloc( maxnx, sizeof(double) );

  // Distribution of points along the x-axis such that points on the surface
  // are approximately distant by spacing
  vx[0] = - a;
  for (i=1;i<maxnx && sol < 0.;++i)
  {
    // Knowing x[i-1], we find x[i] by solving a non-linear equation via the 
    // bi-section method
    lb = vx[i-1];
    
    // At lb, fx is negative, so we look for the first sign change before
    // running the bi-section method to set ub
    ub = lb;    
    while ( fx( ub, lb, a, b, refs ) < 0. ) ub += 1.e-4 * refs;
    
    // Bi-section method
    while ( ( ub - lb ) > 1.e-8 * refs )
    {
      mid = 0.5 * ( lb + ub );
      if ( copysign( 1., fx( mid, vx[i-1], a, b, refs ) ) ==
      	copysign( 1., fx( ub, vx[i-1], a, b, refs ) ) )
        ub = mid;
      else
        lb = mid; 	
    }
    sol = 0.5 * ( lb + ub );
    vx[i] = sol;
    ++ninterior;     
  }
  *totalnx = 2 * ( ninterior - 1 ) + 3;
  
  // Bi-section method to find optimal spacing
  double lbs = 0.5 * refs, ubs = refs, fubs = vx[ninterior], mids, fmids;
  while ( ( ubs - lbs ) > 1.e-6 * refs )
  {
    mids = 0.5 * ( ubs + lbs );
     
    // The function to find the root of is vx[ninterior]
    vx[0] = - a;
    for (i=1;i<ninterior+1;++i)
    {
      // Knowing x[i-1], we find x[i] by solving a non-linear equation via the 
      // bi-section method
      lb = vx[i-1];
    
      // At lb, fx is negative, so we look for the first sign change before
      // running the bi-section method to set ub
      ub = lb;      
      while ( fx( ub, lb, a, b, mids ) < 0. ) ub += 1.e-4 * mids;
    
      // Bi-section method
      while ( ( ub - lb ) > 1.e-8 * refs )
      {
        mid = 0.5 * ( lb + ub );
        if ( copysign( 1., fx( mid, vx[i-1], a, b, mids ) ) ==
      	  copysign( 1., fx( ub, vx[i-1], a, b, mids ) ) )
          ub = mid;
        else
          lb = mid; 	
      }
      sol = 0.5 * ( lb + ub );
      vx[i] = sol;    
    }
    
    fmids = vx[ninterior];
    if ( copysign( 1., fmids ) == copysign( 1., fubs ) )
    {
      ubs = mids;
      fubs = fmids;
    }
    else
      lbs = mids;
  }

  // Computes point distribution with the optimal spacing  
  *spacing = 0.5 * ( ubs + lbs );
  vx[0] = - a;
  for (i=1;i<ninterior+1;++i)
  {
    // Knowing x[i-1], we find x[i] by solving a non-linear equation via the 
    // bi-section method
    lb = vx[i-1];
    
    // At lb, fx is negative, so we look for the first sign change before
    // running the bi-section method to set ub
    ub = lb;      
    while ( fx( ub, lb, a, b, *spacing ) < 0. ) ub += 1.e-4 * *spacing;
    
    // Bi-section method
    while ( ( ub - lb ) > 1.e-8 * refs )
    {
      mid = 0.5 * ( lb + ub );
      if ( copysign( 1., fx( mid, vx[i-1], a, b, *spacing ) ) ==
      	  copysign( 1., fx( ub, vx[i-1], a, b, *spacing ) ) )
        ub = mid;
      else
        lb = mid; 	
    }
    sol = 0.5 * ( lb + ub );
    vx[i] = sol;    
  }
           
  // Points mirroring over positive x
  vx[ninterior] = 0.;
  for (i=ninterior+1;i<*totalnx;++i) vx[i] = - vx[2*ninterior-i]; 
  
  return( vx ); 
} 




/** Computes the number of boundary points on the perimeter of the 3D circular 
cylinder */
//----------------------------------------------------------------------------
void compute_nboundary_Ellipsoid( GeomParameter const* gcp, int* nb ) 
//----------------------------------------------------------------------------
{
  size_t totalnx = 0, i;
  double spacing = 0., a = gcp->elgp->a, b = gcp->elgp->b, rr;
  double* vx = xaxis_points_distribution( gcp, &totalnx, &spacing );
    
  // Since we impose b=c, each perimeter at a given x is a circle, we equally
  // distribute point over each circle  
  for (i=1;i<totalnx-1;++i)
  {
    rr = b * sqrt( 1. - sq( vx[i] / a ) );    
    *nb += (size_t)( 2. * pi * rr / spacing ); 
  }
  *nb += 2;
            
  if ( *nb == 0 )
    printf( "nboundary = 0: No boundary points !!!\n" );    
  free( vx ); vx = NULL; 
}




/** Creates boundary points and normal vectors of the reference 3D circular 
cylinder */
//----------------------------------------------------------------------------
void create_referenceRB_boundary_geomfeatures_Ellipsoid( 
	GeomParameter const* gcp, RigidBodyBoundary* dlm_bd, const int m ) 
//----------------------------------------------------------------------------
{  
  size_t totalnx = 0, i, ntheta = 0, j;
  int isb = 0; 
  double spacing = 0., a = gcp->elgp->a, b = gcp->elgp->b, local_radius, 
  	dangle, local_angle, bin, vecnorm;
  double* vx = xaxis_points_distribution( gcp, &totalnx, &spacing );  

  // x=-a tip
  dlm_bd->bp[isb].x = - a;
  dlm_bd->bp[isb].y = 0.;  
  dlm_bd->bp[isb].z = 0.;
  dlm_bd->outwardnormalvector[isb].x = - 0.25 * a;
  dlm_bd->outwardnormalvector[isb].y = 0.;  
  dlm_bd->outwardnormalvector[isb].z = 0.;  
  ++isb;

  // Points over each circular perimeter at the determined x positions 
  for (i=1;i<totalnx-1;++i)
  {
    local_radius = b * sqrt( 1. - sq( vx[i] / a ) );
    ntheta = (size_t)( 2. * pi * local_radius / spacing );
    dangle = pi / (double)(ntheta);
    
    // odd or even
    if ( i % 2 == 0 ) bin = 0.;
    else bin = 1.;    
    
    for (j=0;j<ntheta;++j)
    {
      local_angle = ( 2. * (double)(j) + bin ) * dangle ;
      
      dlm_bd->bp[isb].x = vx[i];
      dlm_bd->bp[isb].y = cos( local_angle ) * local_radius ;  
      dlm_bd->bp[isb].z = sin( local_angle ) * local_radius;
      dlm_bd->outwardnormalvector[isb].x = 2. * dlm_bd->bp[isb].x / sq( a );
      dlm_bd->outwardnormalvector[isb].y = 2. * dlm_bd->bp[isb].y / sq( b );  
      dlm_bd->outwardnormalvector[isb].z = 2. * dlm_bd->bp[isb].z / sq( b );
      vecnorm = sqrt( sq( dlm_bd->outwardnormalvector[isb].x ) 
      	+ sq( dlm_bd->outwardnormalvector[isb].y )
      	+ sq( dlm_bd->outwardnormalvector[isb].z ) );
      foreach_dimension() 
        dlm_bd->outwardnormalvector[isb].x *= 0.25 * a / vecnorm;	 
      ++isb;      
    }    
  }

  // x=+a tip
  dlm_bd->bp[isb].x = a;
  dlm_bd->bp[isb].y = 0.;  
  dlm_bd->bp[isb].z = 0.;
  dlm_bd->outwardnormalvector[isb].x = 0.25 * a;
  dlm_bd->outwardnormalvector[isb].y = 0.;  
  dlm_bd->outwardnormalvector[isb].z = 0.; 
    
  free( vx ); vx = NULL;     
}




/** Finds cells lying inside the ellipsoid */
//----------------------------------------------------------------------------
void create_FD_Interior_Ellipsoid( RigidBody* p, vector Index,
	vector PeriodicRefCenter, AABB const* ld )
//----------------------------------------------------------------------------
{
  GeomParameter const* gcp = &(p->g);  
  Cache* fd = &(p->Interior);
  Point ppp;

  // Loops over cells in the bounding box of the sphere
  if ( intersect( ld, &(gcp->BBox) ) )
    foreach_region_plus_plus(gcp->BBox.min, gcp->BBox.max) 
      if ( is_leaf(cell) ) 
        if ( is_in_Ellipsoid_geomtest( x, y, z, p ) )
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

  double x2, y2, z2;

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
	  if ( is_in_Ellipsoid_geomtest( x2, y2, z2, p ) )
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




/** Reads geometric parameters of the ellipsoid */
//----------------------------------------------------------------------------
void read_reference_Ellipsoid( GeomParameter* gcp, 
	const double RotMat[3][3] ) 
//----------------------------------------------------------------------------
{    
  char* token = NULL;
  coord v,w,vr,wr;

  // Read number of points, check that it is 2
  size_t np = 0;
  token = strtok(NULL, " " );
  sscanf( token, "%lu", &np );
  if ( np != 2 )
    printf ("Error in number of points in read_reference_Ellipsoid\n");	
    
  // Allocate the EllipsoidGeomParameter structure
  gcp->elgp = (EllipsoidGeomParameter*) malloc( 
  	sizeof(EllipsoidGeomParameter) );
	
  // Read the point that contains the three semi-axis
  foreach_dimension()
  {
    token = strtok(NULL, " " );
    sscanf( token, "%lf", &(v.x) );
  }
  
  // Read the point that contains the two power
  foreach_dimension()
  {
    token = strtok(NULL, " " );
    sscanf( token, "%lf", &(w.x) );
  }
  
  // In case the reference rigid body was sent by the granular solver with 
  // a non zero center of mass and/or a non-zero identity angular position
  // we need to reset it to its neutral reference position
  // Translation    
  foreach_dimension() v.x -= gcp->center.x;
  // Rotation
  matTransposedCoordDotProduct( RotMat, v, &vr );
  gcp->elgp->a = vr.x;
  gcp->elgp->b = vr.y;
  gcp->elgp->c = vr.z;
  
  // Translation    
  foreach_dimension() w.x -= gcp->center.x;
  // Rotation
  matTransposedCoordDotProduct( RotMat, w, &wr );
  gcp->elgp->n1 = wr.x;
  gcp->elgp->n2 = wr.y; 
  
  // We check that the rigid body is an ellipsoid
  if ( fabs( gcp->elgp->n1 - 2. ) > 1.e-8 || fabs( gcp->elgp->n2 - 2. )> 1.e-8 )
    printf ("Error: superquadric is not an ellipsoid in "
    	"read_reference_Ellipsoid\n");
  else
  {
    gcp->elgp->n1 = 2.;
    gcp->elgp->n2 = 2.;    
  }
  
  // For now, only ellipsoids with equal semi-axis in y and z are supported
  if ( fabs( gcp->elgp->b - gcp->elgp->c ) > 1.e-8 * gcp->radius )
    printf ("Error: ellipsoid does not have equal y and z semi-axis in "
    	"read_reference_Ellipsoid\n");	
}




/** Update geometric parameters with the reference rigid body */
//----------------------------------------------------------------------------
void update_Ellipsoid_from_RBRef( GeomParameter* gcp, 
	RigidBody const* RBRef, const double RotMat[3][3] ) 
//----------------------------------------------------------------------------
{        
  // Allocate the EllipsoidGeomParameter structure
  gcp->elgp = (EllipsoidGeomParameter*) malloc( 
  	sizeof(EllipsoidGeomParameter) );
	
  // Parameters
  gcp->elgp->a = RBRef->g.elgp->a;
  gcp->elgp->b = RBRef->g.elgp->b;
  gcp->elgp->c = RBRef->g.elgp->c;
  gcp->elgp->n1 = RBRef->g.elgp->n1;
  gcp->elgp->n2 = RBRef->g.elgp->n2;   
}




/** Frees the geometric parameters of the ellipsoid */
//----------------------------------------------------------------------------
void free_Ellipsoid( GeomParameter* gcp ) 
//----------------------------------------------------------------------------
{  
  // Free the EllipsoidGeomParameter structure
  free( gcp->elgp );
  gcp->elgp = NULL;
}
