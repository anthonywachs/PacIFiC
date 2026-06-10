/**
# Set of functions for a box
*/

# include "DLMFD_Polyhedron.h"


/** Computes the number of boundary points on the surface of the box */
//----------------------------------------------------------------------------
void compute_nboundary_Box( GeomParameter const* gcp, int* nb )
//----------------------------------------------------------------------------
{
  double delta = L0 / (double)(1 << MAXLEVEL) ; 
  int lNdir[3];   

  /* We compute the cube edge in x */    
  double lx = sqrt( 
  	sq( gcp->pgp->cornersCoord[0].x - gcp->pgp->cornersCoord[1].x ) 
  	+ sq( gcp->pgp->cornersCoord[0].y - gcp->pgp->cornersCoord[1].y ) 
  	+ sq( gcp->pgp->cornersCoord[0].z - gcp->pgp->cornersCoord[1].z ) );
	
  /* We compute the number of intervals on the cube edge in x */
  lNdir[0] = floor( lx / ( INTERBPCOEF * delta ) ) + 1;
  
  /* We compute the cube edge in y */    
  double ly = sqrt( 
  	sq( gcp->pgp->cornersCoord[0].x - gcp->pgp->cornersCoord[4].x ) 
  	+ sq( gcp->pgp->cornersCoord[0].y - gcp->pgp->cornersCoord[4].y ) 
  	+ sq( gcp->pgp->cornersCoord[0].z - gcp->pgp->cornersCoord[4].z ) );
	
  /* We compute the number of intervals on the cube edge in y */
  lNdir[1] = floor( ly / ( INTERBPCOEF * delta ) ) + 1;  
  
  /* We compute the cube edge in z */    
  double lz = sqrt( 
  	sq( gcp->pgp->cornersCoord[0].x - gcp->pgp->cornersCoord[3].x ) 
  	+ sq( gcp->pgp->cornersCoord[0].y - gcp->pgp->cornersCoord[3].y ) 
  	+ sq( gcp->pgp->cornersCoord[0].z - gcp->pgp->cornersCoord[3].z ) );
	
  /* We compute the number of intervals on the cube edge in x */
  lNdir[2] = floor( lz / ( INTERBPCOEF * delta ) ) + 1;  

  /* Number of points required for the 12 edges of the cube */
  *nb = ( lNdir[0] - 2 ) * 4 + ( lNdir[1] - 2 ) * 4 + ( lNdir[2] - 2 ) * 4;
  
  /* Number of points required for the 6 faces of the cube */
  *nb += 2 * ( lNdir[0] - 2 ) * ( lNdir[1] - 2 )
  	+ 2 * ( lNdir[0] - 2 ) * ( lNdir[2] - 2 )
	+ 2 * ( lNdir[1] - 2 ) * ( lNdir[2] - 2 );
      
  /* Number of points required for the 8 corners */
  *nb += 8;

  if ( *nb == 0 )
    fprintf( stderr,"nboundary = 0: No boundary points for the"
    	" cube/square !!!\n" );
}




/** Creates boundary points and normal vectors of the reference box */
//----------------------------------------------------------------------------
void create_referenceRB_boundary_geomfeatures_Box( GeomParameter const* gcp,
	RigidBodyBoundary* dlm_bd, const int m )
//----------------------------------------------------------------------------
{
  int nfaces = gcp->pgp->nfaces, nc = gcp->ncorners;
  int iref, i1, i2, i3, isb = 0, npoints, ndir1, ndir2;
  coord gc_to_center_face, surfnormvec, refcorner, dir1, dir2, p4;
  double delta = L0 / (double)(1 << MAXLEVEL), surfnormvecnorm = 0.,
  	cornercomp = 0.25 * gcp->radius / sqrt( 3. ) ;    

  // Note: we arbitrary set the norm of the surface normal vector to 0.25 *
  // circumscribed radius

  /* Normal at the corners */
  coord* corner_normals = (coord*) calloc( nc, sizeof(coord) );
  for (size_t k=0;k<nc;++k)
  {
    foreach_dimension()
      corner_normals[k].x = gcp->pgp->cornersCoord[k].x - gcp->center.x;
    foreach_dimension() 
      corner_normals[k].x = corner_normals[k].x > 0. ? cornercomp : 
      	- cornercomp;     
  }


  /* Add first interior points on surfaces */
  for (int i = 0; i < nfaces; i++)
  {
    npoints = gcp->pgp->numPointsOnFaces[i];

    iref = gcp->pgp->cornersIndex[i][0];
    i1 = gcp->pgp->cornersIndex[i][1];
    i2 = gcp->pgp->cornersIndex[i][npoints-1];
    i3 = gcp->pgp->cornersIndex[i][2];    

    foreach_dimension() 
    {
      refcorner.x = gcp->pgp->cornersCoord[iref].x;
      dir1.x = gcp->pgp->cornersCoord[i1].x;
      dir2.x = gcp->pgp->cornersCoord[i2].x;
      p4.x = gcp->pgp->cornersCoord[i3].x;      
    }
	
    foreach_dimension() 
      gc_to_center_face.x = refcorner.x + dir1.x + dir2.x + p4.x 
      	- gcp->center.x;

    foreach_dimension()
    {
      dir1.x -= refcorner.x;
      dir2.x -= refcorner.x;
    }
    
    ndir1 = floor( sqrt( sq( dir1.x ) + sq( dir1.y ) + sq( dir1.z ) ) 
    	/ ( INTERBPCOEF * delta ) );
    ndir2 = floor( sqrt( sq( dir2.x ) + sq( dir2.y ) + sq( dir2.z ) ) 
    	/ ( INTERBPCOEF * delta ) );

    foreach_dimension()
    {
      dir1.x /= ndir1;
      dir2.x /= ndir2;
    }

    VecVecCrossProduct( dir1, dir2, &surfnormvec );
    if ( VecVecDotProduct( gc_to_center_face, surfnormvec ) < 0. )
      foreach_dimension() surfnormvec.x *= -1.;
    surfnormvecnorm = 0.;
    foreach_dimension() surfnormvecnorm += sq( surfnormvec.x );
    surfnormvecnorm = sqrt( surfnormvecnorm );
    foreach_dimension() surfnormvec.x *= 0.25 * gcp->radius / surfnormvecnorm;
    
    for (int ii = 1; ii <= ndir1-1; ii++)
    {
      for (int jj = 1; jj <= ndir2-1; jj++)
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

  // We have 8 corner points for the box
  int allindextable[8][8] = {{0}};
  int j1,jm1;

  /* Add points on the edges without the corners */
  for (int i = 0; i < nfaces; i++)
  {
    npoints = gcp->pgp->numPointsOnFaces[i];

    for (int j = 0; j < npoints; j++)
    {
      jm1 = gcp->pgp->cornersIndex[i][j];
      j1 = gcp->pgp->cornersIndex[i][(j+1) % npoints];

      foreach_dimension() 
        surfnormvec.x = 0.5 * ( corner_normals[jm1].x + corner_normals[j1].x );
      surfnormvecnorm = 0.;
      foreach_dimension() surfnormvecnorm += sq( surfnormvec.x );
      surfnormvecnorm = sqrt( surfnormvecnorm );
      foreach_dimension() surfnormvec.x *= 0.25 * gcp->radius / surfnormvecnorm; 

      if ( jm1 > j1 )
      {
	if ( allindextable[jm1][j1] == 0 )
	{
          ndir1 = floor( sqrt( sq( gcp->pgp->cornersCoord[jm1].x 
	  			- gcp->pgp->cornersCoord[j1].x ) 
		+ sq( gcp->pgp->cornersCoord[jm1].y 
				- gcp->pgp->cornersCoord[j1].y )
	  	+ sq( gcp->pgp->cornersCoord[jm1].z 
				- gcp->pgp->cornersCoord[j1].z ) ) 
			/ ( INTERBPCOEF * delta ) ) + 1;
	  distribute_points_edge( gcp, gcp->pgp->cornersCoord[jm1], 
	  	gcp->pgp->cornersCoord[j1], dlm_bd, ndir1, isb, surfnormvec );
	  allindextable[jm1][j1] = 1;
	  isb += ndir1 - 2;
	}
      }
      else
      {
	if ( allindextable[j1][jm1] == 0 )
	{
          ndir1 = floor( sqrt( sq( gcp->pgp->cornersCoord[j1].x 
	  			- gcp->pgp->cornersCoord[jm1].x ) 
	  	+ sq( gcp->pgp->cornersCoord[j1].y 
				- gcp->pgp->cornersCoord[jm1].y )
	  	+ sq( gcp->pgp->cornersCoord[j1].z 
				- gcp->pgp->cornersCoord[jm1].z ) ) 
			/ ( INTERBPCOEF * delta ) ) + 1;		
	  distribute_points_edge( gcp, gcp->pgp->cornersCoord[j1], 
	  	gcp->pgp->cornersCoord[jm1], dlm_bd, ndir1, isb, surfnormvec );
	  allindextable[j1][jm1] = 1;
	  isb += ndir1 - 2;
	}
      }
    }
  }

  /* Add the final 8 corners points */
  for (int i = 0; i  < nc; i++)
  {
    foreach_dimension()
    {
      dlm_bd->bp[isb].x = gcp->pgp->cornersCoord[i].x;
      dlm_bd->outwardnormalvector[isb].x = corner_normals[i].x;
    }

    isb++;
  }
  
  free( corner_normals ); corner_normals = NULL;  
}




/** Reads geometric parameters of the box */
//----------------------------------------------------------------------------
void read_reference_Box( GeomParameter* gcp, const double RotMat[3][3] )
//----------------------------------------------------------------------------
{
  char* token = NULL;

  // Read number of corners, check that it is 8
  size_t nc = 0;
  token = strtok(NULL, " " );
  sscanf( token, "%lu", &nc );
  if ( nc != 8 )
    printf ("Error in number of corners in update_Box\n");

  // Allocate the PolyGeomParameter structure
  gcp->pgp = (PolyGeomParameter*) malloc( sizeof(PolyGeomParameter) );
  gcp->pgp->ncorners = nc;

  // Allocate the array of corner coordinates
  gcp->pgp->cornersCoord = (coord*) malloc( nc * sizeof(coord) );

  // Read the point/corner coordinates
  for (size_t i=0;i<nc;++i)
    foreach_dimension()
    {
      token = strtok(NULL, " " );
      sscanf( token, "%lf", &(gcp->pgp->cornersCoord[i].x) );
    }

  // Read number of faces, check that it is 6
  size_t nf = 0;
  token = strtok(NULL, " " );
  sscanf( token, "%lu", &nf );
  if ( nf != 6 )
    printf ("Error in number of faces in update_Box\n");
  gcp->pgp->nfaces = nf;

  // Allocate the array of number of points/corners on each face
  gcp->pgp->numPointsOnFaces = (long int*) malloc( nf * sizeof(long int) );

  // Allocate the array of point/corner indices on each face
  gcp->pgp->cornersIndex = (long int**) malloc( nf * sizeof(long int*) );

  // Read the face indices
  long int nppf = 0;
  for (size_t i=0;i<nf;++i)
  {
    // Read the number of points/corners on the face, check that it is 4
    token = strtok(NULL, " " );
    sscanf( token, "%ld", &nppf );
    if ( nppf != 4 )
      printf ("Error in number of corners per face in update_Box\n");
    gcp->pgp->numPointsOnFaces[i] = nppf;

    // Allocate the point/corner index vector on the face
    gcp->pgp->cornersIndex[i] = (long int*) malloc( nppf * sizeof(long int) );

    // Read the point/corner indices
    for (size_t j=0;j<4;++j)
    {
      token = strtok(NULL, " " );
      sscanf( token, "%ld", &(gcp->pgp->cornersIndex[i][j]));
    }
  }
  
  
  // In case the reference rigid body was sent by the granular solver with 
  // a non zero center of mass and/or a non-zero identity angular position
  // we need to reset all corners to the neutral reference position
  coord v;
  for (size_t i=0;i<nc;++i)
  {
    // Translation    
    foreach_dimension()
      v.x = gcp->pgp->cornersCoord[i].x - gcp->center.x;

    // Rotation
    matTransposedCoordDotProduct( RotMat, v, &(gcp->pgp->cornersCoord[i]) );
  }
  
# if LEVELDIFF_FLAG_U
    // Allocate the array of corner coordinates of the expanded cube
    gcp->pgp->cornersCoordExp = (coord*) malloc( nc * sizeof(coord) ); 
    
    // Compute corner coordinates of the expanded cube
    double delta = L0 / (double)(1 << MAXLEVEL), 
    	width = BOUNDARY_LAYER_THICKNESS_COEF * delta;
    for (size_t i=0;i<nc;++i)
      foreach_dimension()
        gcp->pgp->cornersCoordExp[i].x = gcp->pgp->cornersCoord[i].x
		+ copysign( width, gcp->pgp->cornersCoord[i].x );
# endif    
}
