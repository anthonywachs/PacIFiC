/**
# Set of functions for a dodecahedron
*/

/** Periodic correction */
//----------------------------------------------------------------------------
void periodic_correction_Dodeca( GeomParameter* gcp, coord* pos, 
	vector* pPeriodicRefCenter, const bool setPeriodicRefCenter )
//----------------------------------------------------------------------------
{
  coord ori = {X0, Y0, Z0};
  coord shift = {0., 0., 0.};
  Point lpoint;

  foreach_dimension()
  {
    if ( Period.x )
    {
      shift.x = 0.;
      if ( pos->x > L0 + ori.x )
      {
        pos->x -= L0;
        shift.x = -L0;
      }
      if ( pos->x < 0. + ori.x )
      {
        pos->x += L0;
        shift.x = L0;
      }
    }
  }

  if ( setPeriodicRefCenter )
  {
    Cache poscache = {0};
    lpoint = locate( pos->x, pos->y, pos->z );

    if ( lpoint.level > -1 )
    {
      cache_append( &poscache, lpoint, 0 );
      foreach_cache(poscache)
        foreach_dimension()
          pPeriodicRefCenter->x[] = gcp->center.x + shift.x;
      free( poscache.p );
    }
  }
}




/** Distributes points on an edge of the dodecahedron */
//----------------------------------------------------------------------------
void distribute_points_edge_Dodecahedron(GeomParameter* gcp, 
	coord const corner1, coord const corner2, SolidBodyBoundary* dlm_bd, 
	int const lN, int const istart, vector* pPeriodicRefCenter, 
	const bool setPeriodicRefCenter)
//----------------------------------------------------------------------------
{
  if ( lN > 0 )
  {
    coord dinc;
    coord pos;

    foreach_dimension()
      dinc.x = (corner2.x - corner1.x)/(lN-1);

    for (int i = 1; i <= lN-2; i++)
    {
      pos.x = corner1.x + (double)i * dinc.x;
      pos.y = corner1.y + (double)i * dinc.y;
      pos.z = corner1.z + (double)i * dinc.z;

      periodic_correction( gcp, &pos, pPeriodicRefCenter, 
      	setPeriodicRefCenter );
      
      dlm_bd->x[istart + i -1] = pos.x;
      dlm_bd->y[istart + i -1] = pos.y;
      dlm_bd->z[istart + i -1] = pos.z;
    }
  }
}




/** Determines on which side of a plane a point is */
//----------------------------------------------------------------------------
int determ_posi_plane_Dodecahedron( coord* pointOne, coord* pointTwo,
	coord* pointThree, coord* pointCheck )
//----------------------------------------------------------------------------
{
  // Equation of the plane: ax + by + cz + k = 0
  double x1 = pointOne->x;
  double y1 = pointOne->y;
  double z1 = pointOne->z;
  double x2 = pointTwo->x;
  double y2 = pointTwo->y;
  double z2 = pointTwo->z;
  double x3 = pointThree->x;
  double y3 = pointThree->y;
  double z3 = pointThree->z;
  double xc = pointCheck->x;
  double yc = pointCheck->y;
  double zc = pointCheck->z;

  // Compute coefficients a, b, c
  double cross_a = (y2 - y1) * (z3 - z2) - (z2 - z1) * (y3 - y2);
  double cross_b = (z2 - z1) * (x3 - x2) - (x2 - x1) * (z3 - z2);
  double cross_c = (x2 - x1) * (y3 - y2) - (y2 - y1) * (x3 - x2);

  // Compute the value of k
  double k =  - (cross_a*x1 + cross_b*y1 + cross_c*z1);

  // Plugin the tested point (xc, yc, zc)
  double check_value  = cross_a*xc + cross_b*yc + cross_c*zc + k;

  // Determine the position of the test point relative to the plan
  // created by three surface points
  int retVal;
  if ( fabs(check_value) > 1.e-13 )
    retVal = (k*check_value > 0) - (k*check_value < 0);
  else
    retVal = 0;

  return retVal;
}



/** Tests whether a point lies inside the dodecahedron */
//----------------------------------------------------------------------------
bool is_in_Dodeca_clone( double x1, double y1, double z1,
	const GeomParameter gp )
//----------------------------------------------------------------------------
{
  GeomParameter const* gcp = &gp;

  // Coordinates of the checkpoint, relative to the center of mass
  coord checkpt;
  checkpt.x = x1 - gcp->center.x;
  checkpt.y = y1 - gcp->center.y;
  checkpt.z = z1 - gcp->center.z;

  int nfaces = gcp->pgp->allFaces;
  int iref, i1, i2, ichoice;

  ichoice = 0;
  int npoints;

  int* position = (int *) calloc( nfaces, sizeof(int) );

  for (int i = 0; i < nfaces; i++)
  {
    npoints = gcp->pgp->numPointsOnFaces[i];

    iref = gcp->pgp->cornersIndex[i][ichoice];
    i1 = gcp->pgp->cornersIndex[i][ichoice + 1];
    i2 = gcp->pgp->cornersIndex[i][npoints - 1];

    // Compute the coordinate of the three points relative to the
    // center of mass
    coord refcorner = {gcp->pgp->cornersCoord[iref][0] - gcp->center.x,
    	gcp->pgp->cornersCoord[iref][1] - gcp->center.y,
    	gcp->pgp->cornersCoord[iref][2] - gcp->center.z};

    coord cornerTwo = {gcp->pgp->cornersCoord[i1][0] - gcp->center.x,
    	gcp->pgp->cornersCoord[i1][1] - gcp->center.y,
    	gcp->pgp->cornersCoord[i1][2] - gcp->center.z};

    coord cornerThree = {gcp->pgp->cornersCoord[i2][0] - gcp->center.x,
    	gcp->pgp->cornersCoord[i2][1] - gcp->center.y,
    	gcp->pgp->cornersCoord[i2][2] - gcp->center.z};

    // Judge if the point lies at the same side with origin
    // compared to the plane generated by refcorner, cornerTwo
    // and cornerThree
    position[i] = determ_posi_plane_Dodecahedron( &refcorner, &cornerTwo,
	        &cornerThree, &checkpt );
  }

  bool isin=false;
  int all_position = 0;
  for (int i = 0; i < nfaces; i++) all_position += position[i] ;

  // The test point lies inside the polyhedron only if it lies on the same side
  // as the geometry center relative to all faces of the polyhedron
  if ( all_position == nfaces ) isin = true;

  free( position );

  return isin;
}




/** Tests whether a point lies inside the dodecahedron or any of its 
periodic clones */
//----------------------------------------------------------------------------
bool is_in_Dodecahedron( double x1, double y1, double z1,
	const GeomParameter gp )
//----------------------------------------------------------------------------
{
  //Check if it is in the master particle
  bool status = is_in_Dodeca_clone( x1, y1, z1, gp );
  
  double x2, y2, z2;

  // Check if it is in any periodic clone particle
  if ( gp.nperclones && !status )
    for (int i = 0; i < gp.nperclones && !status; i++)
    {
      GeomParameter clones = gp;
      clones.center = gp.perclonecenters[i];
      x2 = x1 + gp.center.x - clones.center.x;
      y2 = y1 + gp.center.y - clones.center.y;
      z2 = z1 + gp.center.z - clones.center.z;
      status = is_in_Dodeca_clone( x2, y2, z2, gp );
    }
  
  return status;
}




/** Tests whether a point lies inside the dodecahedron or any of its 
periodic clones and assign the proper center of mass coordinates associated to
this point */
//----------------------------------------------------------------------------
bool in_which_Dodecahedron(double x1, double y1, double z1,
	const GeomParameter gp, vector* pPeriodicRefCenter, 
	const bool setPeriodicRefCenter )
//----------------------------------------------------------------------------
{
  // Check if it is in the master particle
  bool status = is_in_Dodeca_clone(x1, y1, z1, gp);
  if ( status && setPeriodicRefCenter )
  {
    Cache poscache = {0};
    Point lpoint;
    lpoint = locate( x1, y1, z1 );

    if ( lpoint.level > -1 )
    {
      cache_append( &poscache, lpoint, 0 );
      foreach_cache(poscache)
        foreach_dimension()
	  pPeriodicRefCenter->x[] = gp.center.x;
      free( poscache.p );
    }
  }

  double x2, y2, z2;

  // Check if it is in any clone particle
  if ( gp.nperclones && !status )
    for (int i = 0; i < gp.nperclones && !status; i++) 
    {
      GeomParameter clones = gp;
      clones.center = gp.perclonecenters[i];
      x2 = x1 + gp.center.x - clones.center.x;
      y2 = y1 + gp.center.y - clones.center.y;
      z2 = z1 + gp.center.z - clones.center.z;
      status = is_in_Dodeca_clone(x2, y2, z2, gp);
      if ( status && setPeriodicRefCenter )
      {
        Cache poscache = {0};
        Point lpoint;
	lpoint = locate( x1, y1, z1 );

        if ( lpoint.level > -1 )
	{
	  cache_append( &poscache, lpoint, 0 );
	  foreach_cache(poscache)
	    foreach_dimension()
              pPeriodicRefCenter->x[] = clones.center.x;
          free( poscache.p );
	}
      }
    }

  return status;
}




/** Computes the number of boundary points on the surface of the dodecahedron */
//----------------------------------------------------------------------------
void compute_nboundary_Dodecahedron( GeomParameter* gcp, int* nb, int* lN )
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









/** Creates boundary points on the surface of the dodecahedron */
//----------------------------------------------------------------------------
void create_FD_Boundary_Dodecahedron( GeomParameter* gcp,
	SolidBodyBoundary* dlm_bd, const int m, const int lN, 
	vector* pPeriodicRefCenter, const bool setPeriodicRefCenter )
//----------------------------------------------------------------------------
{
  int nfaces = gcp->pgp->allFaces;
  int iref, i1, i2, isb = 0, npoints;
  double orig_x, orig_y, orig_z;
  coord pos;

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

    // Add central points on each surface of the Dodecahedron (12 points)
    pos.x = refcorner.x;
    pos.y = refcorner.y;
    pos.z = refcorner.z;
    periodic_correction_Dodeca( gcp, &pos, pPeriodicRefCenter, 
    	setPeriodicRefCenter );

    dlm_bd->x[isb] = pos.x;
    dlm_bd->y[isb] = pos.y;
    dlm_bd->z[isb] = pos.z;
    isb++;

    for (int k = 0; k < npoints; k++)
    {
      i1 = gcp->pgp->cornersIndex[i][k];

      if ( k == npoints - 1 )
        i2 = gcp->pgp->cornersIndex[i][0];
      else
        i2 = gcp->pgp->cornersIndex[i][k+1];

      coord dir1 = {gcp->pgp->cornersCoord[i1][0],
    		gcp->pgp->cornersCoord[i1][1],
    		gcp->pgp->cornersCoord[i1][2]};

      coord dir2 = {gcp->pgp->cornersCoord[i2][0],
   	 	gcp->pgp->cornersCoord[i2][1],
   	 	gcp->pgp->cornersCoord[i2][2]};

      // Insert points on the innerface edges
      distribute_points_edge_Dodecahedron( gcp, refcorner, dir1, dlm_bd, lN, 
      	isb, pPeriodicRefCenter, setPeriodicRefCenter );
      isb += lN - 2;

      foreach_dimension()
      {
        dir1.x -= refcorner.x;
        dir2.x -= refcorner.x;
        dir1.x /= ( lN - 1 );
        dir2.x /= ( lN - 1 );
      }

      for (int ii = 1; ii <= lN-2; ii++)
      {
        for (int jj = 1; jj <= lN-2 - ii; jj++)
        {
          pos.x = refcorner.x + (double) ii * dir1.x
		+ (double) jj * dir2.x;
          pos.y = refcorner.y + (double) ii * dir1.y
		+ (double) jj * dir2.y;
          pos.z = refcorner.z + (double) ii * dir1.z
		+ (double) jj * dir2.z;

          periodic_correction_Dodeca( gcp, &pos, pPeriodicRefCenter, 
    		setPeriodicRefCenter );

          dlm_bd->x[isb] = pos.x;
          dlm_bd->y[isb] = pos.y;
          dlm_bd->z[isb] = pos.z;

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
	  distribute_points_edge_Dodecahedron( gcp, c1, c2, dlm_bd, lN, isb, 
	  	pPeriodicRefCenter, setPeriodicRefCenter );
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
	  distribute_points_edge_Dodecahedron(gcp, c1, c2, dlm_bd, lN, isb, 
	  	pPeriodicRefCenter, setPeriodicRefCenter );
	  allindextable[j1][jm1] = 1;
	  isb += lN - 2;
	}
      }
    }
  }

  /* Add the final 20 corners points */
  for (int i = 0; i  < gcp->ncorners; i++)
  {
    pos.x = gcp->pgp->cornersCoord[i][0];
    pos.y = gcp->pgp->cornersCoord[i][1];
    pos.z = gcp->pgp->cornersCoord[i][2];
    
    periodic_correction_Dodeca( gcp, &pos, pPeriodicRefCenter, 
    	setPeriodicRefCenter );

    dlm_bd->x[isb] = pos.x;
    dlm_bd->y[isb] = pos.y;
    dlm_bd->z[isb] = pos.z;

    isb++;
  }
  
  if ( setPeriodicRefCenter ) synchronize((scalar*){pPeriodicRefCenter->x,
  	pPeriodicRefCenter->y, pPeriodicRefCenter->z});   
}




/** Finds cells lying inside the Dodecahedron */
//----------------------------------------------------------------------------
void create_FD_Interior_Dodecahedron( particle* p, vector Index, 
	vector PeriodicRefCenter )
//----------------------------------------------------------------------------
{
  Cache* c;

  /** Create the cache of the interior points for a dodecahedron */
  c = &(p->Interior);

  GeomParameter gp = p->g;
  foreach()
  {
    if ( in_which_Dodecahedron( x, y, z, gp, &PeriodicRefCenter, true ) )
    {
      cache_append( c, point, 0 );
      /* tag cell with the number of the particle */
      if ( (int)Index.y[] == -1 )
        Index.y[] = p->pnum;
    }
  }

  cache_shrink( c );
}




/** Reads geometric parameters of the Dodecahedron */
//----------------------------------------------------------------------------
void update_Dodecahedron( GeomParameter* gcp)
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
}




/** Frees the geometric parameters of the Dodecahedron */
//----------------------------------------------------------------------------
void free_Dodecahedron( GeomParameter* gcp )
//----------------------------------------------------------------------------
{
  // Free the point/corner coordinate array
  double* cc = NULL;
  for (int i=0; i < gcp->pgp->allPoints; ++i)
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
  for (int i=0;i < gcp->pgp->allFaces; ++i)
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
