/**
# Set of functions for a polyhedron
*/

/** Periodic correction */
//----------------------------------------------------------------------------
void periodic_correction( GeomParameter const* gcp, coord* pos, 
	vector* pPeriodicRefCenter, const bool setPeriodicRefCenter )
//----------------------------------------------------------------------------
{
  coord ori = {X0, Y0, Z0};
  coord shift = {0., 0., 0.};

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
    foreach_point( pos->x, pos->y, pos->z )
      foreach_dimension()
        pPeriodicRefCenter->x[] = gcp->center.x + shift.x;
}




/** Distributes points over an edge of the polyhedron */
//----------------------------------------------------------------------------
void distribute_points_edge( GeomParameter const* gcp, coord const corner1, 
	coord const corner2, RigidBodyBoundary* dlm_bd, int const lN, 
	int const istart, vector* pPeriodicRefCenter, 
	const bool setPeriodicRefCenter )
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
      foreach_dimension()
        pos.x = corner1.x + (double)i * dinc.x;

      periodic_correction( gcp, &pos, pPeriodicRefCenter, 
      	setPeriodicRefCenter );
      
      foreach_dimension()
        dlm_bd->x[istart + i -1] = pos.x;
    }
  }
}




/** Determines on which side of a plane a point lies */
//----------------------------------------------------------------------------
int determ_posi_plane( coord* pointOne, coord* pointTwo,
	coord* pointThree, coord* pointCheck )
//----------------------------------------------------------------------------
{
  // Equation of the plan: ax + by + cz + k = 0
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
  if ( fabs( check_value ) > 1.e-13 )
    retVal = (k*check_value > 0) - (k*check_value < 0);
  else
    retVal = 0;

  return ( retVal );
}




/** Tests whether a point lies inside the polyhedron */
//----------------------------------------------------------------------------
bool is_in_Polyhedron_clone( double x1, double y1, double z1,
	GeomParameter const* gcp )
//----------------------------------------------------------------------------
{
  // Coordinates of the checkpoint, relative to the center of mass
  coord checkpt;
  checkpt.x = x1 - gcp->center.x;
  checkpt.y = y1 - gcp->center.y;
  checkpt.z = z1 - gcp->center.z;

  int nfaces = gcp->pgp->allFaces;
  int iref, i1, i2, ichoice;

  ichoice = 0;
  int npoints;

  int* position = (int *) calloc(nfaces, sizeof(int));

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
    	gcp->pgp->cornersCoord[iref][2] - gcp->center.z} ;

    coord cornerTwo = {gcp->pgp->cornersCoord[i1][0] - gcp->center.x,
    	gcp->pgp->cornersCoord[i1][1] - gcp->center.y,
    	gcp->pgp->cornersCoord[i1][2] - gcp->center.z};

    coord cornerThree = {gcp->pgp->cornersCoord[i2][0] - gcp->center.x,
    	gcp->pgp->cornersCoord[i2][1] - gcp->center.y,
    	gcp->pgp->cornersCoord[i2][2] - gcp->center.z};

    // Judge if the point lies at the same side with origin
    // compared to the plane generated by refcorner, cornerTwo
    // and cornerThree
    position[i] = determ_posi_plane( &refcorner, &cornerTwo,
	        &cornerThree, &checkpt );
  }

  bool isin = false;
  int all_position = 0;
  for (int i = 0; i < nfaces; i++) all_position += position[i] ;

  // The test point lies inside the polyhedron only if it lies on the same side
  // as the center of mass relative to all faces of the polyhedron
  if ( all_position == nfaces ) isin = true;

  free( position );

  return ( isin );
}




/** Tests whether a point lies inside the polyhedron or any of its 
periodic clones */
//----------------------------------------------------------------------------
bool is_in_Polyhedron( double x1, double y1, double z1,
	GeomParameter const* gcp )
//----------------------------------------------------------------------------
{
  // Check if it is in the master rigid body
  bool status = is_in_Polyhedron_clone( x1, y1, z1, gcp );
  
  double x2, y2, z2;

  // Check if it is in any periodic clone rigid body
  if ( gcp->nperclones && !status )
    for (int i = 0; i < gcp->nperclones && !status; i++)
    {
      GeomParameter clone = *gcp;
      clone.center = gcp->perclonecenters[i];
      x2 = x1 + gcp->center.x - clone.center.x;
      y2 = y1 + gcp->center.y - clone.center.y;
      z2 = z1 + gcp->center.z - clone.center.z;
      status = is_in_Polyhedron_clone( x2, y2, z2, gcp );
    }

  return ( status );
}




/** Tests whether a point lies inside the polyhedron or any of its 
periodic clones and assign the proper center of mass coordinates associated to
this point */
//----------------------------------------------------------------------------
bool in_which_Polyhedron( double x1, double y1, double z1,
	GeomParameter const* gcp, vector* pPeriodicRefCenter, 
	const bool setPeriodicRefCenter )
//----------------------------------------------------------------------------
{
  // Check if it is in the master rigid body
  bool status = is_in_Polyhedron_clone( x1, y1, z1, gcp );    
  if ( status && setPeriodicRefCenter )
    foreach_point( x1, y1, z1 )
      foreach_dimension()
        pPeriodicRefCenter->x[] = gcp->center.x;

  double x2, y2, z2;

  // Check if it is in any clone rigid body
  if ( gcp->nperclones && !status )
    for (int i = 0; i < gcp->nperclones && !status; i++) 
    {
      GeomParameter clone = *gcp;
      clone.center = gcp->perclonecenters[i];
      x2 = x1 + gcp->center.x - clone.center.x;
      y2 = y1 + gcp->center.y - clone.center.y;
      z2 = z1 + gcp->center.z - clone.center.z;
      status = is_in_Polyhedron_clone( x2, y2, z2, gcp );
      if ( status && setPeriodicRefCenter )
        foreach_point( x1, y1, z1 )
          foreach_dimension()
            pPeriodicRefCenter->x[] = clone.center.x;
    }

  return ( status );
}




/** Finds cells lying inside the polyhedron */
//----------------------------------------------------------------------------
void create_FD_Interior_Polyhedron( RigidBody* p, vector Index, 
	vector PeriodicRefCenter )
//----------------------------------------------------------------------------
{
  GeomParameter const* gcp = &(p->g);  
  Cache* fd = &(p->Interior);
  Point ppp;

  /** Create the cache of the interior points of the polyhedron */
  foreach(serial)
    if ( in_which_Polyhedron( x, y, z, gcp, &PeriodicRefCenter, true ) )
      if ( (int)Index.y[] == -1 )
      {
	ppp.i = point.i;
        ppp.j = point.j;
        ppp.k = point.k;			
        ppp.level = point.level;
	cache_append( fd, ppp, 0 );
	Index.y[] = p->pnum;
      }

  cache_shrink( fd );
}
