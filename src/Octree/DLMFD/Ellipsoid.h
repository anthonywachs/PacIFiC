/** 
# Set of functions for an ellipsoid
*/

# include "foreach_region_plusplus.h"


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




double fx( double xip1, double xi, double a, double b, double l )
{
  return( sq( xip1 - xi ) + sq( b ) * sq( sqrt( 1. - sq( xip1 / a )  ) 
  	- sqrt( 1. - sq( xi / a ) ) ) - sq( l ) ); 
}

/** Computes the number of boundary points on the perimeter of the 3D circular 
cylinder */
//----------------------------------------------------------------------------
void compute_nboundary_Ellipsoid( GeomParameter const* gcp, int* nb ) 
//----------------------------------------------------------------------------
{
  double delta = L0 / (double)(1 << MAXLEVEL) ;
  double spacing = INTERBPCOEF * delta, lb, ub, mid, sol = - 1.;
  double a = gcp->elgp->a, b = gcp->elgp->b, shift, rr, vx1;
  size_t maxnx = (size_t)( 2. * ( a + b ) / delta ), ninterior = 0, 
  	totalnx = 0, i;
  double* vx = (double*) calloc( maxnx, sizeof(double) );
  bool compress = true; 

  // Distribution of points along the x-axis such that points on the surface
  // are approximately distant by spacing
  vx[0] = - a;
  for (i=1;i<maxnx && sol < 0.;++i)
  {
    // Knowing x[i-1], we find x[i] by solving a non-linear equation via the 
    // bi-section method
    lb = vx[i-1];
    ub = vx[i-1] + 1.2 * spacing;
    while ( ( ub - lb ) > 1.e-8 * gcp->radius )
    {
      mid = 0.5 * ( lb + ub );
      if ( copysign( 1., fx( mid, vx[i-1], a, b, spacing ) ) ==
      	copysign( 1., fx( ub, vx[i-1], a, b, spacing ) ) )
        ub = mid;
      else
        lb = mid; 	
    }
    sol = 0.5 * ( lb + ub );
    vx[i] = sol;
    ++ninterior;     
  }
  totalnx = 2 * ( ninterior - 1 ) + 3;

  // We need to correct the x positions such that x=0 is a point and points 
  // are symmetric wrt x=0
  // We first check whether we need to stretch or compress the x positions
  // We test whether compressing would work, if not we stretch
  // After compression, the inter-point distance on the ellipsoid surface should
  // not be less than 0.9 the prescribed distance
  shift = vx[ninterior] / (double)( ninterior );
  vx1 = vx[1] - shift;
  if ( - a - vx1 > 0. ) compress = false;
  else if ( sqrt( sq( - a - vx1 ) 
  	+ sq( b ) * ( 1. - sq( vx1 / a ) ) ) < 0.9 * spacing ) 
	compress = false;  
    
  if ( compress )
  {
    shift = vx[ninterior] / (double)( ninterior ); 
    for (i=1;i<ninterior+1;++i) vx[i] -= (double)(i) * shift;
  }
  else
  {
    --ninterior;
    totalnx -= 2;
    shift = fabs( vx[ninterior] ) / (double)( ninterior );
    for (i=1;i<ninterior+1;++i) vx[i] += (double)(i) * shift;        
  } 
  for (i=ninterior+1;i<totalnx;++i) vx[i] = - vx[2*ninterior-i];
    
  // Since we impose b=c, each perimeter at a given x is a circle, we equally
  // distribute point over each circle  
  for (i=1;i<totalnx-1;++i)
  {
    rr = b * sqrt( 1. - sq( vx[i] / a ) );    
    *nb += (size_t)( 2. * pi * rr / spacing ); 
  }
  *nb += 2;   
            
  if( *nb == 0 )
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
  double delta = L0 / (double)(1 << MAXLEVEL) ;
  double spacing = INTERBPCOEF * delta, lb, ub, mid, sol = - 1.;
  double a = gcp->elgp->a, b = gcp->elgp->b, shift, local_radius, dangle, 
  	local_angle, bin, norm, vx1;
  size_t maxnx = (size_t)( 2. * ( a + b ) / delta ), ninterior = 0, 
  	totalnx = 0, i, ntheta = 0, j;
  double* vx = (double*) calloc( maxnx, sizeof(double) );
  int isb = 0; 
  bool compress = true; 
    
  // Distribution of points along the x-axis such that points on the surface
  // are approximately distant by spacing
  vx[0] = - a;
  for (i=1;i<maxnx && sol < 0.;++i)
  {
    lb = vx[i-1];
    ub = vx[i-1] + 1.2 * spacing;
    while ( ( ub - lb ) > 1.e-8 * gcp->radius )
    {
      mid = 0.5 * ( lb + ub );
      if ( copysign( 1., fx( mid, vx[i-1], a, b, spacing ) ) ==
      	copysign( 1., fx( ub, vx[i-1], a, b, spacing ) ) )
        ub = mid;
      else
        lb = mid; 	
    }
    sol = 0.5 * ( lb + ub );
    vx[i] = sol;
    ++ninterior;     
  }
  totalnx = 2 * ( ninterior - 1 ) + 3;

  // We need to correct the x positions such that x=0 is a point and points 
  // are symmetric wrt x=0
  // We first check whether we need to stretch or compress the x positions
  // We test whether compressing would work, if not we stretch
  // After compression, the inter-point distance on the ellipsoid surface should
  // not be less than 0.9 the prescribed distance
  shift = vx[ninterior] / (double)( ninterior );
  vx1 = vx[1] - shift;
  if ( - a - vx1 > 0. ) compress = false;
  else if ( sqrt( sq( - a - vx1 ) 
  	+ sq( b ) * ( 1. - sq( vx1 / a ) ) ) < 0.9 * spacing ) 
	compress = false;  
    
  if ( compress )
  {
    shift = vx[ninterior] / (double)( ninterior ); 
    for (i=1;i<ninterior+1;++i) vx[i] -= (double)(i) * shift;
  }
  else
  {
    --ninterior;
    totalnx -= 2;
    shift = fabs( vx[ninterior] ) / (double)( ninterior );
    for (i=1;i<ninterior+1;++i) vx[i] += (double)(i) * shift;        
  } 
  for (i=ninterior+1;i<totalnx;++i) vx[i] = - vx[2*ninterior-i];

  // x=-a tip
  dlm_bd->bp[isb].x = - a;
  dlm_bd->bp[isb].y = 0.;  
  dlm_bd->bp[isb].z = 0.;
  dlm_bd->normal[isb].x = - 0.25 * a;
  dlm_bd->normal[isb].y = 0.;  
  dlm_bd->normal[isb].z = 0.;  
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
      dlm_bd->normal[isb].x = 2. * dlm_bd->bp[isb].x / sq( a );
      dlm_bd->normal[isb].y = 2. * dlm_bd->bp[isb].y / sq( b );  
      dlm_bd->normal[isb].z = 2. * dlm_bd->bp[isb].z / sq( b );     
      norm = sqrt( sq( dlm_bd->normal[isb].x ) + sq( dlm_bd->normal[isb].y )
      	+ sq( dlm_bd->normal[isb].z ) );
      foreach_dimension() dlm_bd->normal[isb].x *= 0.25 * a / norm;	 
      ++isb;      
    }    
  }

  // x=+a tip
  dlm_bd->bp[isb].x = a;
  dlm_bd->bp[isb].y = 0.;  
  dlm_bd->bp[isb].z = 0.;
  dlm_bd->normal[isb].x = 0.25 * a;
  dlm_bd->normal[isb].y = 0.;  
  dlm_bd->normal[isb].z = 0.; 
    
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
