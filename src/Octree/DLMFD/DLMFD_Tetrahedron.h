/**
# Set of functions for a tetrahedron
*/

# include "DLMFD_Polyhedron.h"

/** Computes the number of boundary points on the surface of the tetrahedron */
//----------------------------------------------------------------------------
void compute_nboundary_Tetrahedron( GeomParameter const* gcp, int* nb, int* lN )
//----------------------------------------------------------------------------
{
  double delta = L0 / (double)(1 << MAXLEVEL) ;  
  
  /* Grains sends the tetrahedron circumscribed radius, so to get the 
  tetrahedron edge length we multiply by sqrt(8/3) */
  double lengthedge = gcp->radius * sqrt (8.) / sqrt(3.);  

  /* We compute the number of intervals on the cube edge */
  *lN = floor( lengthedge / ( INTERBPCOEF * delta ) );

  /* The number of points on a cube edge is the number of intervals + 1 */
  *lN += 1;
  
  /* Number of points required for the 6 edges of the tetrahedron */
  *nb += ( *lN - 2 ) * 6;
  
  /* Number of points required for the 4 faces of the tetrahedron inside a 
  triangle */
  *nb += 4 * ( *lN - 2 ) * ( *lN - 3 ) / 2;
  
  /* Number of points required for the 4 corners */
  *nb += 4;  

  if ( *nb == 0 )
    fprintf( stderr,"nboundary = 0: No boundary points for the"
    	" Tetrahedron !!!\n" );
}




/** Creates boundary points and normal vectors of the reference tetrahedron */
//----------------------------------------------------------------------------
void create_referenceRB_boundary_geomfeatures_Tetrahedron( 
	GeomParameter const* gcp,
	RigidBodyBoundary* dlm_bd, const int m, const int lN )
//----------------------------------------------------------------------------
{
  int nfaces = gcp->pgp->nfaces, nc = gcp->ncorners;
  int iref, i1, i2, isb = 0, npoints;
  coord gc_to_center_face, surfnormvec, refcorner, dir1, dir2;
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
    
    iref = gcp->pgp->cornersIndex[i][0];
    i1 = gcp->pgp->cornersIndex[i][1];
    i2 = gcp->pgp->cornersIndex[i][npoints-1];

    foreach_dimension() 
    {
      refcorner.x = gcp->pgp->cornersCoord[iref].x;
      dir1.x = gcp->pgp->cornersCoord[i1].x;
      dir2.x = gcp->pgp->cornersCoord[i2].x;    
    }

    foreach_dimension() 
      gc_to_center_face.x = refcorner.x + dir1.x + dir2.x - gcp->center.x;

    foreach_dimension()
    {
      dir1.x -= refcorner.x;
      dir2.x -= refcorner.x;
      dir1.x /= (lN-1);
      dir2.x /= (lN-1);
    }
    
    VecVecCrossProduct( dir1, dir2, &surfnormvec );
    if ( VecVecDotProduct( gc_to_center_face, surfnormvec ) < 0. )
      foreach_dimension() surfnormvec.x *= -1.;
    surfnormvecnorm = 0.;
    foreach_dimension() surfnormvecnorm += sq( surfnormvec.x );
    surfnormvecnorm = sqrt( surfnormvecnorm );
    foreach_dimension() surfnormvec.x *= 0.25 * gcp->radius / surfnormvecnorm;
    
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

  // We have 4 corner points for the tetrahedron
  int allindextable[4][4] = {{0}};
  int j1,jm1;


  /* Add points on the edges without the corners*/
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

  /* Add the final 4 corners points */
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




// Reads geometric parameters of the tetrahedron
//----------------------------------------------------------------------------
void read_reference_Tetrahedron( GeomParameter* gcp, const double RotMat[3][3] )
//----------------------------------------------------------------------------
{
  char* token = NULL;

  // Read number of corners, check that it is 4
  size_t nc = 0;
  token = strtok(NULL, " " );
  sscanf( token, "%lu", &nc );
  if ( nc != 4 )
    printf ("Error in number of corners in update_Tetrahedron \n");

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

  // Read number of faces, check that it is 4
  size_t nf = 0;
  token = strtok(NULL, " " );
  sscanf( token, "%lu", &nf );
  if ( nf != 4 )
    printf ("Error in number of faces in update_Tetrahedron\n");
  gcp->pgp->nfaces = nf;

  // Allocate the array of number of points/corners on each face
  gcp->pgp->numPointsOnFaces = (long int*) malloc( nf * sizeof(long int) );

  // Allocate the array of point/corner indices on each face
  gcp->pgp->cornersIndex = (long int**) malloc( nf * sizeof(long int*) );

  // Read the face indices
  long int nppf = 0;
  for (size_t i=0;i<nf;++i)
  {
    // Read the number of points/corners on the face, check that it is 3
    token = strtok(NULL, " " );
    sscanf( token, "%ld", &nppf );
    if ( nppf != 3 )
      printf ("Error in number of corners per face in update_Tetrahedron\n");
    gcp->pgp->numPointsOnFaces[i] = nppf;

    // Allocate the point/corner index vector on the face
    gcp->pgp->cornersIndex[i] = (long int*) malloc( nppf * sizeof(long int) );

    // Read the point/corner indices
    for (size_t j=0;j<3;++j)
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
    // Allocate the array of corner coordinates of the expanded tetrahedron
    gcp->pgp->cornersCoordExp = (coord*) malloc( nc * sizeof(coord) ); 
    
    // Compute corner coordinates of the expanded tetrahedron
    double delta = L0 / (double)(1 << MAXLEVEL), 
    	width = BOUNDARY_LAYER_THICKNESS_COEF * delta,
	ext = 3. * width, vecnorm = 0.;
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
