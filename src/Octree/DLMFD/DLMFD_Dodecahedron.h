/**
# Set of functions for a dodecahedron
*/

# include "DLMFD_Polyhedron.h"

/** Computes the number of boundary points on the surface of the dodecahedron */
//----------------------------------------------------------------------------
void compute_nboundary_Dodecahedron( GeomParameter const* gcp, int* nb, 
	int* lN )
//----------------------------------------------------------------------------
{
  double delta = L0 / (double)(1 << MAXLEVEL) ; 
  
  /* Grains sends the dodecahedron circumscribed radius, so to get the
  dodecahedron edge length we multiply by sqrt(3.)/3.*(sqrt(5.)-1) */
  double lengthedge = gcp->radius * sqrt(3.) / 3. * ( sqrt(5.) - 1. );  

  /* We compute the number of intervals on the dodecahedron edge */
  *lN = floor( lengthedge / ( INTERBPCOEF * delta ) );

  /* The number of points on a dodecahedron edge is the number of intervals 
  + 1 */
  *lN += 1;

  /* Number of points required for the 30 edges and 12 * 5 inner edges of the 
  dodecahedron */
  *nb += ( *lN - 2 ) * ( 30 + 12 * 5 );
  
  /* Number of points required for the 12 faces of the dodecahedron */
  *nb += 12 * ( *lN - 2 ) * ( *lN - 3 ) / 2 * 5;
  
  /* Number of points required for the 20 corners and 12 central points */
  *nb += 20 + 12;

  if ( *nb == 0 )
    fprintf( stderr,"nboundary = 0: No boundary points for the"
    	" Dodecahedron !!!\n" );
}




/** Creates boundary points and normal vectors of the reference dodecahedron */
//----------------------------------------------------------------------------
void create_referenceRB_boundary_geomfeatures_Dodecahedron( 
	GeomParameter const* gcp, RigidBodyBoundary* dlm_bd, const int m, 
	const int lN )
//----------------------------------------------------------------------------
{
  int nfaces = gcp->pgp->nfaces, nc = gcp->ncorners;
  int iref, i1, i2, isb = 0, npoints;
  coord gc_to_center_face, surfnormvec, refcorner, dir1, dir2, pt1;
  double surfnormvecnorm = 0.;
  
  // Note: we arbitrary set the norm of the surface normal vector to 0.25 *
  // circumscribed radius

  /* Normal at the corners */
  coord* corner_normals = (coord*) calloc( nc, sizeof(coord) );
  for (size_t k=0;k<nc;++k)
  {
    foreach_dimension()
      corner_normals[k].x = gcp->pgp->cornersCoord[k].x - gcp->center.x;
    surfnormvecnorm = 0.;
    foreach_dimension() surfnormvecnorm += sq( corner_normals[k].x );
    surfnormvecnorm = sqrt( surfnormvecnorm );
    foreach_dimension() corner_normals[k].x *= 0.25 * gcp->radius 
    	/ surfnormvecnorm;       
  }  


  /* Add first interior points on surfaces */
  for (int i = 0; i < nfaces; i++)
  {
    npoints = gcp->pgp->numPointsOnFaces[i];
    foreach_dimension() refcorner.x = 0. ;

    for (int j = 0; j < npoints; j++)
    {
      iref = gcp->pgp->cornersIndex[i][j];
      foreach_dimension() 
       refcorner.x += gcp->pgp->cornersCoord[iref].x / npoints;
    }
    
    foreach_dimension() 
      gc_to_center_face.x = refcorner.x - gcp->center.x;    

    // Add central points on each surface of the Dodecahedron (12 points)
    foreach_dimension()
      dlm_bd->bp[isb].x = refcorner.x;

    for (int k = 0; k < npoints; k++)
    {
      i1 = gcp->pgp->cornersIndex[i][k];
      i2 = gcp->pgp->cornersIndex[i][(k+1) % npoints];

      foreach_dimension() 
      {
        dir1.x = gcp->pgp->cornersCoord[i1].x;
	pt1.x = dir1.x;
	dir2.x = gcp->pgp->cornersCoord[i2].x;
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

  // We have 20 corner points for dodecahedron
  int allindextable[20][20] = {{0}};
  int j1, jm1;

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
	  distribute_points_edge( gcp, gcp->pgp->cornersCoord[jm1], 
	  	gcp->pgp->cornersCoord[j1], dlm_bd, lN, isb, surfnormvec );
	  allindextable[jm1][j1] = 1;
	  isb += lN - 2;
	}
      }
      else
      {
	if ( allindextable[j1][jm1] == 0 )
	{
	  distribute_points_edge( gcp, gcp->pgp->cornersCoord[j1], 
	  	gcp->pgp->cornersCoord[jm1], dlm_bd, lN, isb, surfnormvec );
	  allindextable[j1][jm1] = 1;
	  isb += lN - 2;
	}
      }
    }
  }

  /* Add the final 20 corners points */
  for (size_t i = 0; i < nc; i++)
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



/** Reads geometric parameters of the Dodecahedron */
//----------------------------------------------------------------------------
void read_reference_Dodecahedron( GeomParameter* gcp, 
	const double RotMat[3][3] )
//----------------------------------------------------------------------------
{
  char* token = NULL;

  // Read number of corners, check that it is 20
  size_t nc = 0;
  token = strtok(NULL, " " );
  sscanf( token, "%lu", &nc );
  if ( nc != 20 )
    printf ("Error in number of corners in update_Dodecahedron \n");

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

  // Read number of faces, check that it is 12
  size_t nf = 0;
  token = strtok(NULL, " " );
  sscanf( token, "%lu", &nf );
  if ( nf != 12 )
    printf ("Error in number of faces in update_Dodecahedron\n");
  gcp->pgp->nfaces = nf;

  // Allocate the array of number of points/corners on each face
  gcp->pgp->numPointsOnFaces = (long int*) malloc( nf * sizeof(long int) );

  // Allocate the array of point/corner indices on each face
  gcp->pgp->cornersIndex = (long int**) malloc( nf * sizeof(long int*) );

  // Read the face indices
  long int nppf = 0;
  for (size_t i=0;i<nf;++i)
  {
    // Read the number of points/corners on the face, check that it is 5
    token = strtok(NULL, " " );
    sscanf( token, "%ld", &nppf );
    if ( nppf != 5 )
      printf ("Error in number of corners per face in update_Dodecahedron\n");
    gcp->pgp->numPointsOnFaces[i] = nppf;

    // Allocate the point/corner index vector on the face
    gcp->pgp->cornersIndex[i] = (long int*) malloc( nppf * sizeof(long int) );

    // Read the point/corner indices : attention 5 points on a face
    for (size_t j=0;j<5;++j)
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
    // Allocate the array of corner coordinates of the expanded dodecahedron
    gcp->pgp->cornersCoordExp = (coord*) malloc( nc * sizeof(coord) ); 
    
    // Compute corner coordinates of the expanded dodecahedron
    double delta = L0 / (double)(1 << MAXLEVEL), 
    	width = BOUNDARY_LAYER_THICKNESS_COEF * delta,
	ext = sqrt( 15. / ( 5. + 2. * sqrt( 5. ) ) ) * width, 
	vecnorm = 0.;
    for (size_t i=0;i<nc;++i)
    {	
      vecnorm = sqrt( sq( gcp->pgp->cornersCoord[i].x ) 
      		+ sq( gcp->pgp->cornersCoord[i].y ) 
		+ sq( gcp->pgp->cornersCoord[i].z ) ); 
      foreach_dimension()
        gcp->pgp->cornersCoordExp[i].x = ( 1. + ext / vecnorm )
		* gcp->pgp->cornersCoord[i].x ;		
    }          
# endif  
}
