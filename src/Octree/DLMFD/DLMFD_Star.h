/** 
# Set of functions for a 6-branch star 
*/

# include "DLMFD_Polyhedron.h"


/** Tests whether a point lies inside the 6-branch star */
//----------------------------------------------------------------------------
bool is_in_Star_geomtest( const double x, const double y, 
	const double z, RigidBody const* p, const double explength )
//----------------------------------------------------------------------------
{
  bool isin = false; 
  coord v, vr, bounds;
  bounds.x = p->g.sbsgp->armLength + p->g.sbsgp->armWidth 
  	* cos( M_PI / 6. ) + explength;
  bounds.y = p->g.sbsgp->armWidth / 2. + explength;	
  bounds.z = p->g.sbsgp->depth / 2. + explength;	
  double wx, wy;
  
  // The neutral position of the 6-branch star is depth in z and
  // one branch aligned with x

  // Translation   
  v.x = x - p->g.center.x;
  v.y = y - p->g.center.y;  
  v.z = z - p->g.center.z;  
  
  // Rotation
  matTransposedCoordDotProduct( p->RotMat, v, &vr ); 
  
  if ( vr.z >= - bounds.z && vr.z <= bounds.z )
  {
    // x-aligned branch
    if ( vr.x >= - bounds.x && vr.x <= bounds.x )
      if ( vr.y >= - bounds.y && vr.y <= bounds.y ) isin = true;
	
    // +pi/3 branch
    // We rotate by -pi/3 and do the same check as the x-aligned branch
    if ( !isin )
    {
      wx = 0.5 * vr.x + 0.5 * sqrt( 3. ) * vr.y;
      wy = - 0.5 * sqrt( 3. ) * vr.x + 0.5 * vr.y;
      if ( wx >= - bounds.x && wx <= bounds.x )
        if ( wy >= - bounds.y && wy <= bounds.y ) isin = true;        
    }
        
    // +2*pi/3 branch 
    // We rotate by -2*pi/3 and do the same check as the x-aligned branch
    if ( !isin )
    {
      wx = - 0.5 * vr.x + 0.5 * sqrt( 3. ) * vr.y;
      wy = - 0.5 * sqrt( 3. ) * vr.x - 0.5 * vr.y;
      if ( wx >= - bounds.x && wx <= bounds.x )
        if ( wy >= - bounds.y && wy <= bounds.y ) isin = true;        
    }       
  } 
  	 		  
  return ( isin );
}




/** Tests whether a point lies inside the 6-branch star or any of its periodic 
clones */
//----------------------------------------------------------------------------
bool is_in_Star( const double x1, const double y1, 
	const double z1, RigidBody const* p )
//----------------------------------------------------------------------------
{
  // Check if it is in the master rigid body
  bool status = is_in_Star_geomtest( x1, y1, z1, p, 0. );

  double x2, y2, z2;

  // Check if it is in any clone rigid body
  if ( p->g.nperclones && !status )
    for (int i = 0; i < p->g.nperclones && !status; i++)
    {
      x2 = x1 + p->g.center.x - p->g.perclonecenters[i].x;
      y2 = y1 + p->g.center.y - p->g.perclonecenters[i].y;
      z2 = z1 + p->g.center.z - p->g.perclonecenters[i].z;
      status = is_in_Star_geomtest( x2, y2, z2, p, 0. );
    }

  return ( status );
}




/** Flag boundary layer around the 6-branch star */
//----------------------------------------------------------------------------
void flag_boundarylayer_Star( scalar flag_maxlevel, double const dcoef, 
	RigidBody const* p, AABB const* ld )
//----------------------------------------------------------------------------
{
  GeomParameter const* gcp = &(p->g); 
  AABB ExpBBox;
  double delta = L0 / (double)(1 << MAXLEVEL), x2, y2, z2 ; 

  foreach_dimension()
  {
    ExpBBox.min.x = gcp->BBox.min.x - dcoef * delta;
    ExpBBox.max.x = gcp->BBox.max.x + dcoef * delta;
  } 
      
  // Loops over cells in the bounding box of the expanded polyhedron
  if ( intersect( ld, &ExpBBox ) )  
    foreach_region_plus_plus( ExpBBox.min, ExpBBox.max ) 
      if ( is_leaf(cell) )
        if ( flag_maxlevel[] == 0. )
          if ( is_in_Star_geomtest( x, y, z, p, dcoef * delta ) )
	    flag_maxlevel[] = 1.;
	
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

            if ( is_in_Star_geomtest( x2, y2, z2, p, dcoef * delta ) )
	      flag_maxlevel[] = 1.;
          }
  }	
}




/** Computes the number of boundary points on the perimeter of the 6-branch 
star */
//----------------------------------------------------------------------------
void compute_nboundary_Star( GeomParameter const* gcp, int* nb, int* lN, 
	int* lH, int* lL, int* lB ) 
//----------------------------------------------------------------------------
{
  double delta = L0 / (double)(1 << MAXLEVEL);
  double spacing = INTERBPCOEF * delta;
  
  // Core hexagon
  /* Number of points on the hexagonal prism edge (i.e. arm width) */
  *lN = floor( gcp->sbsgp->armWidth / spacing ) + 1 ;
  
  /* Number of points on the hexagonal prism depth */
  *lH = floor( gcp->sbsgp->depth / spacing ) + 1;
  
  /* Number of points required over the 12 hexa edges */
  *nb += 12 * ( *lN - 2 );
  
  /* Number of points required over the 6 height edges */
  *nb += 6 * ( *lH - 2 ); 

  /* Number of points required over the inner edges */
  *nb += 12 * ( *lN - 2 );      

  /* Number of points required over the 2 hexagonal faces */
  *nb += 12 * ( *lN - 2 ) * ( *lN - 3 ) / 2;  

  /* Number of points required for the 12 corners and 2 central points */
  *nb += 12 + 2;
  

  // Six box branches
  /* Number of points on the branch (we do not +1 to avoid duplicates with
  the core hexagon edge points */ 
  *lL = floor( gcp->sbsgp->armLength / spacing ) ;
  
  /* Number of points required over each branch */
  *lB = 2 * ( *lN - 2 ) * ( *lL - 1 ) + 2 * ( *lH - 2 ) * ( *lL - 1 )
  	+ ( *lN - 2 ) * ( *lH - 2 ) + 4 * ( *lL - 1)
	+ 2 * ( *lN - 2 ) + 2 * ( *lH - 2 ) + 4;    
    
  /* Number of points required over the 6 branches */
  *nb += 6 * *lB;
        
  if ( *nb == 0 ) printf( "nboundary = 0: No boundary points !!!\n" );
}




/** Creates boundary points and normal vectors of the reference 6-branch star */
//----------------------------------------------------------------------------
void create_referenceRB_boundary_geomfeatures_Star( GeomParameter const* gcp,
	RigidBodyBoundary* dlm_bd, const int m, const int lN, const int lH,
	const int lL, const int lB ) 
//----------------------------------------------------------------------------
{
  double ndist = gcp->sbsgp->armWidth * cos( M_PI / 6. ), angle, 
	surfnormvecnorm = 0., sq2o2 = sqrt( 2. ) / 2., sq3o3 = sqrt( 3. ) / 3.;
  int iref, nc = 12, i1, i2, npoints, isb = 0, itb = 0;
  coord gc_to_center_face, surfnormvec, dir1, dir2, refcorner, pt1;

  // Note: we arbitrary set the vecnorm of the normal vector to 0.25 * radius

  // We first redefine the corners and face connectivity of the hexagonal prism
  // to recycle the code used for that shape for the core hexagonal prism of
  // the 6-branch star
  coord* cornersCoord = (coord*) calloc( 12, sizeof(coord) );
  int** cornersIndex = (int**) malloc( 8 * sizeof(int*) );
  for (size_t i=0;i<6;++i) cornersIndex[i] = (int*) malloc( 4 * sizeof(int) );
  cornersIndex[6] = (int*) malloc( 6 * sizeof(int) );
  cornersIndex[7] = (int*) malloc( 6 * sizeof(int) );

  cornersIndex[0][0] = 6; 
  cornersIndex[0][1] = 0; 
  cornersIndex[0][2] = 1;
  cornersIndex[0][3] = 7;

  cornersIndex[1][0] = 7; 
  cornersIndex[1][1] = 1; 
  cornersIndex[1][2] = 2;
  cornersIndex[1][3] = 8;

  cornersIndex[2][0] = 8; 
  cornersIndex[2][1] = 2; 
  cornersIndex[2][2] = 3;
  cornersIndex[2][3] = 9;
  
  cornersIndex[3][0] = 9; 
  cornersIndex[3][1] = 3; 
  cornersIndex[3][2] = 4;
  cornersIndex[3][3] = 10;
  
  cornersIndex[4][0] = 10; 
  cornersIndex[4][1] = 4; 
  cornersIndex[4][2] = 6;
  cornersIndex[4][3] = 11; 
   
  cornersIndex[5][0] = 11; 
  cornersIndex[5][1] = 5; 
  cornersIndex[5][2] = 0;
  cornersIndex[5][3] = 6; 

  cornersIndex[6][0] = 0; 
  cornersIndex[6][1] = 5; 
  cornersIndex[6][2] = 4;
  cornersIndex[6][3] = 3;
  cornersIndex[6][4] = 2;
  cornersIndex[6][5] = 1;  

  cornersIndex[7][0] = 6; 
  cornersIndex[7][1] = 7; 
  cornersIndex[7][2] = 8;
  cornersIndex[7][3] = 9;
  cornersIndex[7][4] = 10;
  cornersIndex[7][5] = 11;  
  
  for (size_t i=0;i<6;++i)
  {
    angle = M_PI / 6. + (double)i * M_PI / 3.;
    cornersCoord[i].x = gcp->sbsgp->armWidth * cos( angle); 
    cornersCoord[i].y = gcp->sbsgp->armWidth * sin( angle);
    cornersCoord[i].z = - 0.5 * gcp->sbsgp->depth;        
  }    
  for (size_t i=0;i<6;++i)
  {
    angle = M_PI / 6. + (double)i * M_PI / 3.;
    cornersCoord[i+6].x = gcp->sbsgp->armWidth * cos( angle); 
    cornersCoord[i+6].y = gcp->sbsgp->armWidth * sin( angle);
    cornersCoord[i+6].z = 0.5 * gcp->sbsgp->depth;        
  }  


  // We distribute points on the core hexagonal prism
  /* Normal at the corners */
  coord* corner_normals = (coord*) calloc( nc, sizeof(coord) );
  for (size_t k=0;k<nc;++k)
  {
    foreach_dimension()
      corner_normals[k].x = cornersCoord[k].x - gcp->center.x;
    surfnormvecnorm = sqrt( sq( corner_normals[k].x ) 
    	+ sq( corner_normals[k].y ) );    
    corner_normals[k].z = corner_normals[k].z > 0. ? 
    	surfnormvecnorm : - surfnormvecnorm;
    surfnormvecnorm = 0.;
    foreach_dimension() surfnormvecnorm += sq( corner_normals[k].x );
    surfnormvecnorm = sqrt( surfnormvecnorm );
    foreach_dimension() corner_normals[k].x *= 0.25 * gcp->radius 
    	/ surfnormvecnorm;       
  }  

  /* Points over edges */
  surfnormvec.x = surfnormvec.y = 0.; 
  for (size_t i=6;i<8;++i)
  {
    for (size_t k = 0; k < 6; ++k)
    {
      i1 = cornersIndex[i][k];
      i2 = cornersIndex[i][(k+1) % 6];
    
      surfnormvec.z = ( cornersCoord[i1].z > 0. ? 1. : -1. ) 
      	* 0.25 * gcp->radius; 
    
      distribute_points_edge( gcp, cornersCoord[i1], cornersCoord[i2], dlm_bd, 
      	lN, isb, surfnormvec );
      isb += lN - 2; 
    }     
  }
  
  for (size_t i=0;i<6;++i)
  {
    i1 = i;
    i2 = i + 6;	
    
    foreach_dimension() 
      surfnormvec.x = 0.5 * ( corner_normals[i1].x + corner_normals[i2].x );
    surfnormvecnorm = 0.;
    foreach_dimension() surfnormvecnorm += sq( surfnormvec.x );
    surfnormvecnorm = sqrt( surfnormvecnorm );
    foreach_dimension() surfnormvec.x *= 0.25 * gcp->radius / surfnormvecnorm;
    
    distribute_points_edge( gcp, cornersCoord[i1], cornersCoord[i2], dlm_bd,
    	lH, isb, surfnormvec );
    isb += lH - 2; 
  }
  
  /* Points over hexagonal faces */
  for (size_t i=6;i<8;++i)
  {
    npoints = 6;
    foreach_dimension() refcorner.x = 0. ;

    for (int j = 0; j < npoints; j++)
    {
      iref = cornersIndex[i][j];
      foreach_dimension() 
        refcorner.x += cornersCoord[iref].x / npoints;
    }
    
    foreach_dimension() 
      gc_to_center_face.x = refcorner.x - gcp->center.x;    

    // Add central point on each hexagonal face
    foreach_dimension()
      dlm_bd->bp[isb].x = refcorner.x;    

    for (int k = 0; k < npoints; k++)
    {
      i1 = cornersIndex[i][k];
      i2 = cornersIndex[i][(k+1) % npoints];

      foreach_dimension() 
      {
        dir1.x = cornersCoord[i1].x;
	pt1.x = dir1.x;
	dir2.x = cornersCoord[i2].x;
      }

      foreach_dimension()
      {
        dir1.x -= refcorner.x;
        dir2.x -= refcorner.x;
        dir1.x /= ( lN - 1 );
        dir2.x /= ( lN - 1 );
      }
		
      if ( !k )
      {
        // Compute surfnormvec vector of the face
	VecVecCrossProduct( dir1, dir2, &surfnormvec );
        if ( VecVecDotProduct( gc_to_center_face, surfnormvec ) < 0. )
          foreach_dimension() surfnormvec.x *= -1.;
        surfnormvecnorm = 0.;
        foreach_dimension() surfnormvecnorm += sq( surfnormvec.x );
        surfnormvecnorm = sqrt( surfnormvecnorm );
        foreach_dimension() surfnormvec.x *= 0.25 * gcp->radius 
		/ surfnormvecnorm;
	
	// Set surfnormvec vector of the central point
	foreach_dimension()
          dlm_bd->outwardnormalvector[isb].x = surfnormvec.x;
	isb++;
      }		

      // Insert points on the innerface edges
      distribute_points_edge( gcp, refcorner, pt1, dlm_bd, lN, isb, 
      	surfnormvec );
      isb += lN - 2;

      // Insert points on the innerface triangles
      for (int ii = 1; ii <= lN-2; ii++)
      {
        for (int jj = 1; jj <= lN-2 - ii; jj++)
        {
          foreach_dimension()
	  {
	    dlm_bd->bp[isb].x = refcorner.x + (double) ii * dir1.x
		      + (double) jj * dir2.x;
	    dlm_bd->outwardnormalvector[isb].x = surfnormvec.x ;	      
	  }
          isb++;
        }
      }
    }    
  }      

  /* Add the 12 corners points */
  for (size_t i = 0; i < nc; i++)
  {
    foreach_dimension()
    {
      dlm_bd->bp[isb].x = cornersCoord[i].x;
      dlm_bd->outwardnormalvector[isb].x = corner_normals[i].x;
    }
    isb++;
  }


  // We distribute points on the six branches
  // Template branch aligned with +x
  coord* branch_points = (coord*) calloc( lB, sizeof(coord) );  
  coord* branch_normals = (coord*) calloc( lB, sizeof(coord) );

  // 2 xy faces
  refcorner.x = ndist; 
  refcorner.y = - 0.5 * gcp->sbsgp->armWidth;    
  dir1.x = gcp->sbsgp->armLength / lL; dir1.y = 0.; dir1.z = 0.;
  dir2.x = 0.; dir2.y = gcp->sbsgp->armWidth / ( lN - 1 ); dir2.z = 0.;
  surfnormvec.x = 0.; surfnormvec.y = 0.; 
  for ( size_t k=0; k<2; ++k)
  {
    refcorner.z = pow( -1., (double)k ) * 0.5 * gcp->sbsgp->depth;
    surfnormvec.z = pow( -1., (double)k ) * 0.25 * gcp->radius;
    for (int ii = 1; ii <= lL-1; ii++)
    {
      for (int jj = 1; jj <= lN-2; jj++)
      {
        foreach_dimension()
	{
	  branch_points[itb].x = refcorner.x + (double) ii * dir1.x
      		+ (double) jj * dir2.x;
	  branch_normals[itb].x = surfnormvec.x ;	
	}
      	itb++;
      }
    }           
  }
  
  // 2 xz faces
  refcorner.x = ndist; 
  refcorner.z = - 0.5 * gcp->sbsgp->depth;    
  dir1.x = gcp->sbsgp->armLength / lL; dir1.y = 0.; dir1.z = 0.;
  dir2.x = 0.; dir2.y = 0.; dir2.z = gcp->sbsgp->depth / ( lH - 1 );
  surfnormvec.x = 0.; surfnormvec.z = 0.; 
  for ( size_t k=0; k<2; ++k)
  {
    refcorner.y = pow( -1., (double)k ) * 0.5 * gcp->sbsgp->armWidth;
    surfnormvec.y = pow( -1., (double)k ) * 0.25 * gcp->radius;
    for (int ii = 1; ii <= lL-1; ii++)
    {
      for (int jj = 1; jj <= lH-2; jj++)
      {
        foreach_dimension()
	{
	  branch_points[itb].x = refcorner.x + (double) ii * dir1.x
      		+ (double) jj * dir2.x;
	  branch_normals[itb].x = surfnormvec.x ;	
	}
      	itb++;
      }
    }           
  }
  
  // The tip yz face
  refcorner.x = ndist + gcp->sbsgp->armLength;
  refcorner.y = - 0.5 * gcp->sbsgp->armWidth;   
  refcorner.z = - 0.5 * gcp->sbsgp->depth;    
  dir1.x = 0.; dir1.y = gcp->sbsgp->armWidth / ( lN - 1 ); dir1.z = 0.;
  dir2.x = 0.; dir2.y = 0.; dir2.z = gcp->sbsgp->depth / ( lH - 1 );
  surfnormvec.x = 0.25 * gcp->radius; surfnormvec.y = 0.; surfnormvec.z = 0.; 
  for (int ii = 1; ii <= lN-2; ii++)
  {
    for (int jj = 1; jj <= lH-2; jj++)
    {
      foreach_dimension()
      {
	branch_points[itb].x = refcorner.x + (double) ii * dir1.x
      		+ (double) jj * dir2.x;
	branch_normals[itb].x = surfnormvec.x ;	
      }
      itb++;
    }           
  }
  
  // The 4 +x edges
  refcorner.x = ndist; 
  surfnormvec.x = 0.;  
  dir1.x = gcp->sbsgp->armLength / lL; dir1.y = 0.; dir1.z = 0.; 
  
  // -y,-z
  refcorner.y = - 0.5 * gcp->sbsgp->armWidth;   
  refcorner.z = - 0.5 * gcp->sbsgp->depth;
  surfnormvec.y = - sq2o2 * 0.25 * gcp->radius;   
  surfnormvec.z = - sq2o2 * 0.25 * gcp->radius; 
  for (int ii = 1; ii <= lL-1; ii++)
  {
    foreach_dimension()
    {
      branch_points[itb].x = refcorner.x + (double) ii * dir1.x;
      branch_normals[itb].x = surfnormvec.x ;	
    }
    itb++;
  }
  
  // -y,+z
  refcorner.y = - 0.5 * gcp->sbsgp->armWidth;   
  refcorner.z = 0.5 * gcp->sbsgp->depth;
  surfnormvec.y = - sq2o2 * 0.25 * gcp->radius;   
  surfnormvec.z = sq2o2 * 0.25 * gcp->radius; 
  for (int ii = 1; ii <= lL-1; ii++)
  {
    foreach_dimension()
    {
      branch_points[itb].x = refcorner.x + (double) ii * dir1.x;
      branch_normals[itb].x = surfnormvec.x ;	
    }
    itb++;
  } 
  
  // +y,-z
  refcorner.y = 0.5 * gcp->sbsgp->armWidth;   
  refcorner.z = - 0.5 * gcp->sbsgp->depth;
  surfnormvec.y = sq2o2 * 0.25 * gcp->radius;   
  surfnormvec.z = - sq2o2 * 0.25 * gcp->radius; 
  for (int ii = 1; ii <= lL-1; ii++)
  {
    foreach_dimension()
    {
      branch_points[itb].x = refcorner.x + (double) ii * dir1.x;
      branch_normals[itb].x = surfnormvec.x ;	
    }
    itb++;
  }  
  
  // +y,+z
  refcorner.y = 0.5 * gcp->sbsgp->armWidth;   
  refcorner.z = 0.5 * gcp->sbsgp->depth;
  surfnormvec.y = sq2o2 * 0.25 * gcp->radius;   
  surfnormvec.z = sq2o2 * 0.25 * gcp->radius; 
  for (int ii = 1; ii <= lL-1; ii++)
  {
    foreach_dimension()
    {
      branch_points[itb].x = refcorner.x + (double) ii * dir1.x;
      branch_normals[itb].x = surfnormvec.x ;	
    }
    itb++;
  } 
  
  // The 2 y edges
  refcorner.x = ndist + gcp->sbsgp->armLength;
  refcorner.y = - 0.5 * gcp->sbsgp->armWidth;  
  surfnormvec.x = sq2o2 * 0.25 * gcp->radius;     
  surfnormvec.y = 0.;  
  dir1.x = 0.; dir1.y = gcp->sbsgp->armWidth / ( lN - 1 ); dir1.z = 0.; 
  
  // -z 
  refcorner.z = - 0.5 * gcp->sbsgp->depth; 
  surfnormvec.z = - sq2o2 * 0.25 * gcp->radius; 
  for (int ii = 1; ii <= lN-2; ii++)
  {
    foreach_dimension()
    {
      branch_points[itb].x = refcorner.x + (double) ii * dir1.x;
      branch_normals[itb].x = surfnormvec.x ;	
    }
    itb++;
  }
  
  // +z 
  refcorner.z = 0.5 * gcp->sbsgp->depth; 
  surfnormvec.z = sq2o2 * 0.25 * gcp->radius; 
  for (int ii = 1; ii <= lN-2; ii++)
  {
    foreach_dimension()
    {
      branch_points[itb].x = refcorner.x + (double) ii * dir1.x;
      branch_normals[itb].x = surfnormvec.x ;	
    }
    itb++;
  }
  
  // The 2 z edges
  refcorner.x = ndist + gcp->sbsgp->armLength;
  refcorner.z = - 0.5 * gcp->sbsgp->depth;  
  surfnormvec.x = sq2o2 * 0.25 * gcp->radius;     
  surfnormvec.z = 0.;  
  dir1.x = 0.; dir1.y = 0.; dir1.z = gcp->sbsgp->depth / ( lH - 1 );
  
  // -y 
  refcorner.y = - 0.5 * gcp->sbsgp->armWidth; 
  surfnormvec.y = - sq2o2 * 0.25 * gcp->radius; 
  for (int ii = 1; ii <= lH-2; ii++)
  {
    foreach_dimension()
    {
      branch_points[itb].x = refcorner.x + (double) ii * dir1.x;
      branch_normals[itb].x = surfnormvec.x ;	
    }
    itb++;
  }
  
  // +y 
  refcorner.y = 0.5 * gcp->sbsgp->armWidth; 
  surfnormvec.y = sq2o2 * 0.25 * gcp->radius; 
  for (int ii = 1; ii <= lH-2; ii++)
  {
    foreach_dimension()
    {
      branch_points[itb].x = refcorner.x + (double) ii * dir1.x;
      branch_normals[itb].x = surfnormvec.x ;	
    }
    itb++;
  } 
  
  // The 4 tip corners
  refcorner.x = ndist + gcp->sbsgp->armLength; 
  surfnormvec.x = sq3o3 * 0.25 * gcp->radius;    
  
  // -y,-z
  refcorner.y = - 0.5 * gcp->sbsgp->armWidth;   
  refcorner.z = - 0.5 * gcp->sbsgp->depth;
  surfnormvec.y = - sq3o3 * 0.25 * gcp->radius;   
  surfnormvec.z = - sq3o3 * 0.25 * gcp->radius; 
  foreach_dimension()
  {
    branch_points[itb].x = refcorner.x ;
    branch_normals[itb].x = surfnormvec.x ;	
  }
  itb++;
  
  // -y,+z
  refcorner.y = - 0.5 * gcp->sbsgp->armWidth;   
  refcorner.z = 0.5 * gcp->sbsgp->depth;
  surfnormvec.y = - sq2o2 * 0.25 * gcp->radius;   
  surfnormvec.z = sq2o2 * 0.25 * gcp->radius; 
  foreach_dimension()
  {
    branch_points[itb].x = refcorner.x ;
    branch_normals[itb].x = surfnormvec.x ;	
  }
  itb++; 
  
  // +y,-z
  refcorner.y = 0.5 * gcp->sbsgp->armWidth;   
  refcorner.z = - 0.5 * gcp->sbsgp->depth;
  surfnormvec.y = sq2o2 * 0.25 * gcp->radius;   
  surfnormvec.z = - sq2o2 * 0.25 * gcp->radius; 
  foreach_dimension()
  {
    branch_points[itb].x = refcorner.x ;
    branch_normals[itb].x = surfnormvec.x ;	
  }
  itb++;  
  
  // +y,+z
  refcorner.y = 0.5 * gcp->sbsgp->armWidth;   
  refcorner.z = 0.5 * gcp->sbsgp->depth;
  surfnormvec.y = sq2o2 * 0.25 * gcp->radius;   
  surfnormvec.z = sq2o2 * 0.25 * gcp->radius; 
  foreach_dimension()
  {
    branch_points[itb].x = refcorner.x ;
    branch_normals[itb].x = surfnormvec.x ;	
  }
 
  // Distribution by i * PI / 3 rotation
  double rotang[6] = { 0., M_PI / 3., 2. * M_PI / 3., 3. * M_PI / 3.,
  	4. * M_PI / 3., 5. * M_PI / 3.};
  double MatRot[3][3];
  MatRot[0][2] = MatRot[1][2] = MatRot[2][0] = MatRot[2][1] = 0.;
  MatRot[2][2] = 1.;
  for (size_t k=0; k<6; ++k)
  {
    // Assign rotation matrix
    MatRot[0][0] = cos( rotang[k] );
    MatRot[0][1] = - sin( rotang[k] );
    MatRot[1][0] = sin( rotang[k] );
    MatRot[1][1] = cos( rotang[k] );
    
    // Copy rotated points from template branch
    for (size_t i=0; i<lB; ++i)
    {
      matCoordDotProduct( MatRot, branch_points[i], &(dlm_bd->bp[isb]) );
      matCoordDotProduct( MatRot, branch_normals[i], 
      	&(dlm_bd->outwardnormalvector[isb]) );
      isb++;
    }               
  }


  // We free locally allocated memory
  free( cornersCoord );
  cornersCoord = NULL;
  int* in = NULL;
  for (size_t i=0;i<8;++i)
  {
    in = &(cornersIndex[i][0]);
    free( in );
    in = NULL;
  }
  free( cornersIndex );
  cornersIndex = NULL; 
  free( corner_normals );
  corner_normals = NULL;
  free( branch_points );
  branch_points = NULL;
  free( branch_normals );
  branch_normals = NULL;       
}




/** Finds cells lying inside the 6-branch star */
//----------------------------------------------------------------------------
void create_FD_Interior_Star( RigidBody* p, vector Index,
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
        if ( is_in_Star_geomtest( x, y, z, p, 0. ) )
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
	  if ( is_in_Star_geomtest( x2, y2, z2, p, 0. ) )
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




/** Reads geometric parameters of the 6-branch star */
//----------------------------------------------------------------------------
void read_reference_Star( GeomParameter* gcp, const double RotMat[3][3] ) 
//----------------------------------------------------------------------------
{    
  char* token = NULL;

  // Allocate the SixBranchStarGeomParameter structure
  gcp->sbsgp = (SixBranchStarGeomParameter*) malloc( 
  	sizeof(SixBranchStarGeomParameter) );
  
  // Read the arm width, arm length and depth
  token = strtok( NULL, " " );
  sscanf( token, "%lf", &(gcp->sbsgp->armWidth) );
  token = strtok( NULL, " " );
  sscanf( token, "%lf", &(gcp->sbsgp->armLength) );    
  token = strtok( NULL, " " );
  sscanf( token, "%lf", &(gcp->sbsgp->depth) );
}




/** Update geometric parameters with the reference rigid body */
//----------------------------------------------------------------------------
void update_Star_from_RBRef( GeomParameter* gcp, 
	RigidBody const* RBRef, const double RotMat[3][3] ) 
//----------------------------------------------------------------------------
{        
  // Allocate the SixBranchStarGeomParameter structure
  gcp->sbsgp = (SixBranchStarGeomParameter*) malloc( 
  	sizeof(SixBranchStarGeomParameter) );
    
  // Assign the arm width, arm length and depth
  gcp->sbsgp->armWidth = RBRef->g.sbsgp->armWidth;
  gcp->sbsgp->armLength = RBRef->g.sbsgp->armLength;
  gcp->sbsgp->depth = RBRef->g.sbsgp->depth;			  
}




/** Frees the geometric parameters of the 6-branch star */
//----------------------------------------------------------------------------
void free_Star( GeomParameter* gcp ) 
//----------------------------------------------------------------------------
{  
  // Free the SixBranchStarGeomParameter structure
  free( gcp->sbsgp );
  gcp->sbsgp = NULL;
}
