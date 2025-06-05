/**
# Set of functions for a box
*/

# include "Polyhedron.h"

/** Computes the number of boundary points on the surface of the box */
//----------------------------------------------------------------------------
void compute_nboundary_Box( GeomParameter const* gcp, int* nb )
//----------------------------------------------------------------------------
{
  double delta = L0 / (double)(1 << MAXLEVEL) ; 
  int lNdir[3];   

  /* We compute the cube edge in x */    
  double lx = sqrt( 
  	sq( gcp->pgp->cornersCoord[0][0] - gcp->pgp->cornersCoord[1][0] ) 
  	+ sq( gcp->pgp->cornersCoord[0][1] - gcp->pgp->cornersCoord[1][1] ) 
  	+ sq( gcp->pgp->cornersCoord[0][2] - gcp->pgp->cornersCoord[1][2] ) );
	
  /* We compute the number of intervals on the cube edge in x */
  lNdir[0] = floor( lx / ( INTERBPCOEF * delta ) ) + 1;
  
  /* We compute the cube edge in y */    
  double ly = sqrt( 
  	sq( gcp->pgp->cornersCoord[0][0] - gcp->pgp->cornersCoord[4][0] ) 
  	+ sq( gcp->pgp->cornersCoord[0][1] - gcp->pgp->cornersCoord[4][1] ) 
  	+ sq( gcp->pgp->cornersCoord[0][2] - gcp->pgp->cornersCoord[4][2] ) );
	
  /* We compute the number of intervals on the cube edge in y */
  lNdir[1] = floor( ly / ( INTERBPCOEF * delta ) ) + 1;  
  
  /* We compute the cube edge in z */    
  double lz = sqrt( 
  	sq( gcp->pgp->cornersCoord[0][0] - gcp->pgp->cornersCoord[3][0] ) 
  	+ sq( gcp->pgp->cornersCoord[0][1] - gcp->pgp->cornersCoord[3][1] ) 
  	+ sq( gcp->pgp->cornersCoord[0][2] - gcp->pgp->cornersCoord[3][2] ) );
	
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




/** Creates boundary points on the surface of the box */
//----------------------------------------------------------------------------
void create_FD_Boundary_Box( GeomParameter const* gcp,
	RigidBodyBoundary* dlm_bd, const int m, 
	vector* pPeriodicRefCenter, const bool setPeriodicRefCenter )
//----------------------------------------------------------------------------
{
  int nfaces = gcp->pgp->allFaces;
  int iref, i1, i2, ichoice = 0, isb = 0, npoints, ndir1, ndir2;
  coord pos;
  double delta = L0 / (double)(1 << MAXLEVEL) ;    

  /* Add first interior points on surfaces */
  for (int i = 0; i < nfaces; i++)
  {
    npoints = gcp->pgp->numPointsOnFaces[i];

    iref = gcp->pgp->cornersIndex[i][ichoice];
    i1 = gcp->pgp->cornersIndex[i][ichoice + 1];
    i2 = gcp->pgp->cornersIndex[i][npoints-1];

    coord refcorner = {gcp->pgp->cornersCoord[iref][0],
    	gcp->pgp->cornersCoord[iref][1],
    	gcp->pgp->cornersCoord[iref][2]} ;

    coord dir1 = {gcp->pgp->cornersCoord[i1][0],
    	gcp->pgp->cornersCoord[i1][1],
    	gcp->pgp->cornersCoord[i1][2]};

    coord dir2 = {gcp->pgp->cornersCoord[i2][0],
    	gcp->pgp->cornersCoord[i2][1],
    	gcp->pgp->cornersCoord[i2][2]};

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
    
    for (int ii = 1; ii <= ndir1-1; ii++)
    {
      for (int jj = 1; jj <= ndir2-1; jj++)
      {
        pos.x = refcorner.x + (double) ii * dir1.x
      		+ (double) jj * dir2.x;
        pos.y = refcorner.y + (double) ii * dir1.y
      		+ (double) jj * dir2.y;
        pos.z = refcorner.z + (double) ii * dir1.z
      		+ (double) jj * dir2.z;
		
        periodic_correction( gcp, &pos, pPeriodicRefCenter, 
		setPeriodicRefCenter );

        foreach_dimension()
	  dlm_bd->x[isb] = pos.x;

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
    i1 = gcp->pgp->cornersIndex[i][1];

    for (int j = 1; j < npoints; j++)
    {
      jm1 = gcp->pgp->cornersIndex[i][j-1];
      j1 = gcp->pgp->cornersIndex[i][j];

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
          ndir1 = floor( sqrt( sq( c1.x - c2.x ) + sq( c1.y - c2.y )
	  	+ sq( c1.z - c2.z ) ) / ( INTERBPCOEF * delta ) ) + 1;
	  distribute_points_edge( gcp, c1, c2, dlm_bd, ndir1, isb, 
	  	pPeriodicRefCenter, setPeriodicRefCenter );
	  allindextable[jm1][j1] = 1;
	  isb += ndir1 - 2;
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
          ndir1 = floor( sqrt( sq( c1.x - c2.x ) + sq( c1.y - c2.y )
	  	+ sq( c1.z - c2.z ) ) / ( INTERBPCOEF * delta ) ) + 1;		
	  distribute_points_edge( gcp, c1, c2, dlm_bd, ndir1, isb, 
	  	pPeriodicRefCenter, setPeriodicRefCenter );
	  allindextable[j1][jm1] = 1;
	  isb += ndir1 - 2;
	}
      }
    }
  }

  /* Add the final 8 corners points */
  for (int i = 0; i  < gcp->ncorners; i++)
  {
    pos.x = gcp->pgp->cornersCoord[i][0];
    pos.y = gcp->pgp->cornersCoord[i][1];
    pos.z = gcp->pgp->cornersCoord[i][2];
    
    periodic_correction( gcp, &pos, pPeriodicRefCenter, setPeriodicRefCenter );

    foreach_dimension()
      dlm_bd->x[isb] = pos.x;

    isb++;
  }
  
  if ( setPeriodicRefCenter ) synchronize((scalar*){pPeriodicRefCenter->x,
  	pPeriodicRefCenter->y, pPeriodicRefCenter->z});  
}




/** Reads geometric parameters of the box */
//----------------------------------------------------------------------------
void update_Box( GeomParameter* gcp )
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

  // Read number of faces, check that it is 6
  size_t nf = 0;
  token = strtok(NULL, " " );
  sscanf( token, "%lu", &nf );
  if ( nf != 6 )
    printf ("Error in number of faces in update_Box\n");
  gcp->pgp->allFaces = nf;

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
}




/** Frees the geometric parameters of the box */
//----------------------------------------------------------------------------
void free_Box( GeomParameter* gcp )
//----------------------------------------------------------------------------
{
  // Free the point/corner coordinate array
  double* cc = NULL;
  for (size_t i=0;i<gcp->pgp->allPoints;++i)
  {
    cc = &(gcp->pgp->cornersCoord[i][0]);
    free( cc );
    cc = NULL;
  }
  free( gcp->pgp->cornersCoord );
  gcp->pgp->cornersCoord = NULL;
  gcp->pgp->allPoints = 0;

  // Free the point/corner arrays
  long int* in = NULL;
  for (size_t i=0;i<gcp->pgp->allFaces;++i)
  {
    in = &(gcp->pgp->cornersIndex[i][0]);
    free( in );
    in = NULL;
  }
  free( gcp->pgp->cornersIndex );
  gcp->pgp->cornersIndex = NULL;
  free( gcp->pgp->numPointsOnFaces );
  gcp->pgp->numPointsOnFaces = NULL;
  gcp->pgp->allFaces = 0;

  // Free the PolyGeomParameter structure
  free( gcp->pgp );
  gcp->pgp = NULL;
}
