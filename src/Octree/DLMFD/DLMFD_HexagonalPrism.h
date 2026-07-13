/**
# Set of functions for a hexagonal prism
In the reference position:
* Width is in the z direction
* Vertices are numbered as follows: 6 vertices with negative z followed by 6
vertices with positive z
* Faces are numbered as follows: 0-5 side faces, 6 bottom, 7 top 
*/

# include "DLMFD_Polyhedron.h"

/** Computes the number of boundary points on the surface of the hexagonal 
prism */
//----------------------------------------------------------------------------
void compute_nboundary_HexagonalPrism( GeomParameter const* gcp, int* nb, 
	int* lN, int* lH )
//----------------------------------------------------------------------------
{
  double delta = L0 / (double)(1 << MAXLEVEL), lengthedge, height ;  
  *nb = 0;
  
  /* Compute edge length and height */
  lengthedge = sqrt( 
  	sq( gcp->pgp->cornersCoord[1].x - gcp->pgp->cornersCoord[0].x )
  	+ sq( gcp->pgp->cornersCoord[1].y - gcp->pgp->cornersCoord[0].y )
	+ sq( gcp->pgp->cornersCoord[1].z - gcp->pgp->cornersCoord[0].z ) );
  height = sqrt( sq( gcp->pgp->cornersCoord[6].x - gcp->pgp->cornersCoord[0].x )
  	+ sq( gcp->pgp->cornersCoord[6].y - gcp->pgp->cornersCoord[0].y )
	+ sq( gcp->pgp->cornersCoord[6].z - gcp->pgp->cornersCoord[0].z ) );

  /* We compute the number of intervals on the hexagonal prism edge */
  *lN = floor( lengthedge / ( INTERBPCOEF * delta ) ); 
  
  /* The number of points on a hexagonal prism edge is the number of intervals 
  + 1 */
  *lN += 1;
  
  /* We compute the number of intervals on the hexagonal prism height */
  *lH = floor( height / ( INTERBPCOEF * delta ) ); 
  
  /* The number of points on a hexagonal prism height is the number of 
  intervals + 1 */
  *lH += 1;     

  /* Number of points required over the 12 hexa edges */
  *nb += 12 * ( *lN - 2 );
  
  /* Number of points required over the 6 height edges */
  *nb += 6 * ( *lH - 2 ); 

  /* Number of points required over the inner edges */
  *nb += 12 * ( *lN - 2 ); 
  
  /* Number of points required over the 6 rectangular faces */
  *nb += 6 * ( *lH - 2 ) * ( *lN - 2 );     

  /* Number of points required over the 2 hexagonal faces */
  *nb += 12 * ( *lN - 2 ) * ( *lN - 3 ) / 2;  

  /* Number of points required for the 12 corners and 2 central points */
  *nb += 12 + 2;    

  if ( *nb == 0 )
    fprintf( stderr,"nboundary = 0: No boundary points for the"
    	" HexagonalPrism !!!\n" );
}




/** Creates boundary points and normal vectors of the reference hexagonal 
prism */
//----------------------------------------------------------------------------
void create_referenceRB_boundary_geomfeatures_HexagonalPrism( 
	GeomParameter const* gcp, RigidBodyBoundary* dlm_bd, const int m, 
	const int lN, const int lH )
//----------------------------------------------------------------------------
{
  int nc = gcp->ncorners, npoints;
  int iref, isb = 0, i1, i2, ndir1, ndir2;  
  coord gc_to_center_face, surfnormvec, refcorner, dir1, dir2, pt1;
  double surfnormvecnorm = 0., delta = L0 / (double)(1 << MAXLEVEL);
  
  // Note: we arbitrary set the norm of the surface normal vector to 0.25 *
  // circumscribed radius
  
  /* Normal at the corners */
  coord* corner_normals = (coord*) calloc( nc, sizeof(coord) );
  for (size_t k=0;k<nc;++k)
  {
    foreach_dimension()
      corner_normals[k].x = gcp->pgp->cornersCoord[k].x - gcp->center.x;
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
  for (size_t i=6;i<8;++i)
  {
    for (size_t k = 0; k < 6; ++k)
    {
      i1 = gcp->pgp->cornersIndex[i][k];
      i2 = gcp->pgp->cornersIndex[i][(k+1) % 6];		
    
      foreach_dimension() 
        surfnormvec.x = 0.5 * ( corner_normals[i1].x + corner_normals[i2].x );
      surfnormvecnorm = 0.;
      foreach_dimension() surfnormvecnorm += sq( surfnormvec.x );
      surfnormvecnorm = sqrt( surfnormvecnorm );
      foreach_dimension() surfnormvec.x *= 0.25 * gcp->radius 
      	/ surfnormvecnorm;
    
      distribute_points_edge( gcp, gcp->pgp->cornersCoord[i1], 
      	gcp->pgp->cornersCoord[i2], dlm_bd, lN, isb, surfnormvec );
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
    
    distribute_points_edge( gcp, gcp->pgp->cornersCoord[i1], 
    	gcp->pgp->cornersCoord[i2], dlm_bd, lH, isb, surfnormvec );
    isb += lH - 2; 
  }        


  /* Points over lateral faces */
  for (int i = 0; i < 6; i++)
  {
    i1 = gcp->pgp->cornersIndex[i][0];
    i2 = gcp->pgp->cornersIndex[i][1];    
    foreach_dimension()
    {
      refcorner.x = gcp->pgp->cornersCoord[i1].x;
      dir1.x = gcp->pgp->cornersCoord[i2].x - refcorner.x;
    }
    
    surfnormvecnorm = 0.;
    foreach_dimension() surfnormvecnorm += sq( dir1.x );
    surfnormvecnorm = sqrt( surfnormvecnorm );
    ndir1 = floor( surfnormvecnorm / ( INTERBPCOEF * delta ) );
    foreach_dimension() dir1.x /= ndir1;    
        
    i2 = gcp->pgp->cornersIndex[i][3];
    foreach_dimension()
      dir2.x = gcp->pgp->cornersCoord[i2].x - refcorner.x;
    surfnormvecnorm = 0.;
    foreach_dimension() surfnormvecnorm += sq( dir2.x );
    surfnormvecnorm = sqrt( surfnormvecnorm );
    ndir2 = floor( surfnormvecnorm / ( INTERBPCOEF * delta ) );
    foreach_dimension() dir2.x /= ndir2;     
    
    foreach_dimension() surfnormvec.x = 0.;
    for (size_t k=0;k<4;++k)
      foreach_dimension() 
        surfnormvec.x += corner_normals[gcp->pgp->cornersIndex[i][k]].x;
    foreach_dimension() surfnormvec.x /= 4.;
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
  
  
  /* Points over hexagonal faces */
  for (size_t i=6;i<8;++i)
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

    // Add central point on each hexagonal face
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


  /* Add the 12 corners points */
  for (size_t i = 0; i < nc; i++)
  {
    foreach_dimension()
    {
      dlm_bd->bp[isb].x = gcp->pgp->cornersCoord[i].x;
      dlm_bd->outwardnormalvector[isb].x = corner_normals[i].x;
    }

    isb++;
  }
  
  free( corner_normals );
  corner_normals = NULL;
}




/** Reads geometric parameters of the hexagonal prism */
//----------------------------------------------------------------------------
void read_reference_HexagonalPrism( GeomParameter* gcp, 
	const double RotMat[3][3] )
//----------------------------------------------------------------------------
{
  char* token = NULL;

  // Read number of corners, check that it is 12
  size_t nc = 0;
  token = strtok(NULL, " " );
  sscanf( token, "%lu", &nc );
  if ( nc != 12 )
    printf ("Error in number of corners in update_HexagonalPrism \n");

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

  // Read number of faces, check that it is 8
  size_t nf = 0;
  token = strtok(NULL, " " );
  sscanf( token, "%lu", &nf );
  if ( nf != 8 )
    printf ("Error in number of faces in update_HexagonalPrism\n");
  gcp->pgp->nfaces = nf;

  // Allocate the array of number of points/corners on each face
  gcp->pgp->numPointsOnFaces = (long int*) malloc( nf * sizeof(long int) );

  // Allocate the array of point/corner indices on each face
  gcp->pgp->cornersIndex = (long int**) malloc( nf * sizeof(long int*) );

  // Read the face indices
  long int nppf = 0;
  for (size_t i=0;i<nf;++i)
  {
    // Read the number of points/corners on the face
    token = strtok(NULL, " " );
    sscanf( token, "%ld", &nppf );
    gcp->pgp->numPointsOnFaces[i] = nppf;

    // Allocate the point/corner index vector on the face
    gcp->pgp->cornersIndex[i] = (long int*) malloc( nppf * sizeof(long int) );

    // Read the point/corner indices
    for (size_t j=0;j<nppf;++j)
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
    // Allocate the array of corner coordinates of the expanded hexagonal prism
    gcp->pgp->cornersCoordExp = (coord*) malloc( nc * sizeof(coord) ); 
    
    // Compute corner coordinates of the expanded hexagonal prism
    double delta = L0 / (double)(1 << MAXLEVEL), 
    	width = BOUNDARY_LAYER_THICKNESS_COEF * delta,
	hext = 2. * width / sqrt( 3. ), vecnorm = 0.;
    for (size_t i=0;i<nc;++i)
    {	
      vecnorm = sqrt( sq( gcp->pgp->cornersCoord[i].x ) 
      		+ sq( gcp->pgp->cornersCoord[i].y ) ); 
      gcp->pgp->cornersCoordExp[i].x = ( 1. + hext / vecnorm )
		* gcp->pgp->cornersCoord[i].x ;
      gcp->pgp->cornersCoordExp[i].y = ( 1. + hext / vecnorm )
		* gcp->pgp->cornersCoord[i].y ;
      gcp->pgp->cornersCoordExp[i].z = gcp->pgp->cornersCoord[i].z
      		+ ( i < 6 ? - width : width ) ;
    }          
# endif  
}
