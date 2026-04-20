/**
# Set of functions for a dodecahedron
*/

# include "Polyhedron.h"

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
  int nfaces = gcp->pgp->allFaces, nc = gcp->ncorners;
  int iref, i1, i2, isb = 0, npoints;
  double orig_x, orig_y, orig_z;
  coord pos, gc_to_center_face, normal;
  double norm = 0.;
  
  // Note: we arbitrary set the norm of the normal vector to 0.25 *
  // circumscribed radius

  /* Normal at the corners */
  coord* corner_normals = (coord*) calloc( nc, sizeof(coord) );
  for (size_t k=0;k<nc;++k)
  {
    corner_normals[k].x = gcp->pgp->cornersCoord[k][0] - gcp->center.x;
    corner_normals[k].y = gcp->pgp->cornersCoord[k][1] - gcp->center.y;    
    corner_normals[k].z = gcp->pgp->cornersCoord[k][2] - gcp->center.z; 
    norm = 0.;
    foreach_dimension() norm += sq( corner_normals[k].x );
    norm = sqrt( norm );
    foreach_dimension() corner_normals[k].x *= 0.25 * gcp->radius / norm;       
  }  


  /* Add first interior points on surfaces */
  for (int i = 0; i < nfaces; i++)
  {
    npoints = gcp->pgp->numPointsOnFaces[i];
    orig_x = 0.;
    orig_y = 0.;
    orig_z = 0.;

    for (int j = 0; j < npoints; j++)
    {
      iref = gcp->pgp->cornersIndex[i][j];
      orig_x += gcp->pgp->cornersCoord[iref][0] / npoints;
      orig_y += gcp->pgp->cornersCoord[iref][1] / npoints;
      orig_z += gcp->pgp->cornersCoord[iref][2] / npoints;
    }

    coord refcorner = { orig_x, orig_y, orig_z } ;
    
    foreach_dimension() 
      gc_to_center_face.x = refcorner.x - gcp->center.x;    

    // Add central points on each surface of the Dodecahedron (12 points)
    foreach_dimension()
      dlm_bd->bp[isb].x = refcorner.x;

    for (int k = 0; k < npoints; k++)
    {
      i1 = gcp->pgp->cornersIndex[i][k];
      i2 = gcp->pgp->cornersIndex[i][(k+1) % npoints];

      coord dir1 = {gcp->pgp->cornersCoord[i1][0],
    		gcp->pgp->cornersCoord[i1][1],
    		gcp->pgp->cornersCoord[i1][2]};
      coord pt1;
      foreach_dimension() pt1.x = dir1.x;		

      coord dir2 = {gcp->pgp->cornersCoord[i2][0],
   	 	gcp->pgp->cornersCoord[i2][1],
   	 	gcp->pgp->cornersCoord[i2][2]};

      foreach_dimension()
      {
        dir1.x -= refcorner.x;
        dir2.x -= refcorner.x;
        dir1.x /= ( lN - 1 );
        dir2.x /= ( lN - 1 );
      }
		
      if ( !k )
      {
        // Compute normal vector of the face
	VecVecCrossProduct( dir1, dir2, &normal );
        if ( VecVecDotProduct( gc_to_center_face, normal ) < 0. )
          foreach_dimension() normal.x *= -1.;
        norm = 0.;
        foreach_dimension() norm += sq( normal.x );
        norm = sqrt( norm );
        foreach_dimension() normal.x *= 0.25 * gcp->radius / norm;
	
	// Set normal vector of the central point
	foreach_dimension()
          dlm_bd->normal[isb].x = normal.x;
	isb++;
      }		

      // Insert points on the innerface edges
      distribute_points_edge( gcp, refcorner, pt1, dlm_bd, lN, isb, normal );
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
	    dlm_bd->normal[isb].x = normal.x ;	      
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
        normal.x = 0.5 * ( corner_normals[jm1].x + corner_normals[j1].x );
      norm = 0.;
      foreach_dimension() norm += sq( normal.x );
      norm = sqrt( norm );
      foreach_dimension() normal.x *= 0.25 * gcp->radius / norm;

      if ( jm1 > j1 )
      {
	if ( allindextable[jm1][j1] == 0 )
	{
	  coord c1 = {gcp->pgp->cornersCoord[jm1][0],
	  	gcp->pgp->cornersCoord[jm1][1],
	  	gcp->pgp->cornersCoord[jm1][2]};
	  coord c2 = {gcp->pgp->cornersCoord[j1][0],
	  	gcp->pgp->cornersCoord[j1][1],
	  	gcp->pgp->cornersCoord[j1][2]};
	  distribute_points_edge( gcp, c1, c2, dlm_bd, lN, isb, normal );
	  allindextable[jm1][j1] = 1;
	  isb += lN - 2;
	}
      }
      else
      {
	if ( allindextable[j1][jm1] == 0 )
	{
	  coord c1 = {gcp->pgp->cornersCoord[j1][0],
	  	gcp->pgp->cornersCoord[j1][1],
	  	gcp->pgp->cornersCoord[j1][2]};
	  coord c2 = {gcp->pgp->cornersCoord[jm1][0],
	  	gcp->pgp->cornersCoord[jm1][1],
	  	gcp->pgp->cornersCoord[jm1][2]};
	  distribute_points_edge( gcp, c1, c2, dlm_bd, lN, isb, normal );
	  allindextable[j1][jm1] = 1;
	  isb += lN - 2;
	}
      }
    }
  }

  /* Add the final 20 corners points */
  for (size_t i = 0; i < nc; i++)
  {
    pos.x = gcp->pgp->cornersCoord[i][0];
    pos.y = gcp->pgp->cornersCoord[i][1];
    pos.z = gcp->pgp->cornersCoord[i][2];

    foreach_dimension()
    {
      dlm_bd->bp[isb].x = pos.x;
      dlm_bd->normal[isb].x = corner_normals[i].x;
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
  gcp->pgp->allPoints = nc;

  // Allocate the array of corner coordinates
  gcp->pgp->cornersCoord = (double**) malloc( nc * sizeof(double*) );
  for (size_t i=0;i<nc;i++)
    gcp->pgp->cornersCoord[i] = (double*) malloc( 3 * sizeof(double) );

  // Read the point/corner coordinates
  for (size_t i=0;i<nc;++i)
    for (size_t j=0;j<3;++j)
    {
      token = strtok(NULL, " " );
      sscanf( token, "%lf", &(gcp->pgp->cornersCoord[i][j]) );
    }

  // Read number of faces, check that it is 12
  size_t nf = 0;
  token = strtok(NULL, " " );
  sscanf( token, "%lu", &nf );
  if ( nf != 12 )
    printf ("Error in number of faces in update_Dodecahedron\n");
  gcp->pgp->allFaces = nf;

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
  double v[3];
  for (size_t i=0;i<nc;++i)
  {
    // Translation    
    v[0] = gcp->pgp->cornersCoord[i][0] - gcp->center.x;
    v[1] = gcp->pgp->cornersCoord[i][1] - gcp->center.y;
    v[2] = gcp->pgp->cornersCoord[i][2] - gcp->center.z;

    // Rotation
    matTransposedVecDotProduct( RotMat, v, gcp->pgp->cornersCoord[i] );
  } 
}
