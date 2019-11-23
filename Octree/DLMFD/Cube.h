/** Set of functions for the cube as fictitious domain */
void compute_nboundary_Cube_v2( GeomParameter* gcp, int* nb, int* lN ) 
{
  Cache poscache = {0};
  Point lpoint;
  *nb = 0;
  coord pos = {0., 0., 0.};
  int ip = 0;

  while ( *nb == 0 ) 
  {
    pos.x = gcp->cornersCoord[ip][0];
    pos.y = gcp->cornersCoord[ip][1];
    
#if dimension == 3
    pos.z = gcp->cornersCoord[ip][2];
    lpoint = locate( pos.x, pos.y, pos.z );
#elif dimension ==2
    lpoint = locate( pos.x, pos.y );
#endif
   
    /** Only one thread has the point in its domain (works in serial
	too). */
    if ( lpoint.level > -1 ) 
    {    
      /** Only this thread creates the Cache ... */
      cache_append(&poscache, lpoint, 0);

      /** and only this thread computes the number of boundary points
	  ... */
    
      /* Grains sends the cube circumscribed radius, so to get the cube edge
      length we multiply by 2/sqrt(3) */
      double lengthedge = 2. * gcp->radius / sqrt(3.);
    
      /* We compute the number of intervals on the cube edge */ 
      foreach_cache (poscache) 
      {
	*lN = floor( lengthedge / ( INTERBPCOEF * Delta ) );

//         double actual_INTERBPCOEF = ( lengthedge / *lN ) / Delta;
//         fprintf( stderr, "actual inter-BP coef = %6.4f\n", actual_INTERBPCOEF );
	
        /* The numberof points on a cube edge is the number of intervals + 1 */
        *lN += 1;      
      }
      
#if dimension == 2
      /* number of points required for the 4 edges of the square */
      *nb += (*lN-2)*4;
      /* number of points required for the face */
      *nb += (*lN-2)*(*lN-2);
      /* number of points required for the 4 corners */
      *nb += 4;
#elif dimension == 3
      /* number of points required for the 12 edges of the cube */
      *nb += (*lN-2)*12;
      /* number of points required for the 6 faces of the cube */
      *nb += 6*(*lN-2)*(*lN-2);
      /* number of points required for the 8 corners */
      *nb += 8;
#endif
      
      /** and finally, this thread destroys the cache. */
      free(poscache.p);
    }
  
#if _MPI
    MPI_Barrier(MPI_COMM_WORLD);
    mpi_all_reduce(*nb, MPI_INT, MPI_MAX);
    mpi_all_reduce(*lN, MPI_INT, MPI_MAX);
#endif
    if ( ip < gcp->ncorners )
      ip++;
    else
      break;
  }
  
  if ( *nb == 0 )
    fprintf( stderr,"nboundary = 0: No boundary points for the"
    	" cube/square !!!\n" );
}




void distribute_points_edge_Cube_v2 (coord const corner1, coord const corner2, 
	SolidBodyBoundary * dlm_bd, int const lN, int const istart) 
{
# if dimension == 3 
  if ( lN > 0 ) 
  {
    coord dinc; 

    foreach_dimension()
      dinc.x = (corner2.x - corner1.x)/(lN-1);

    for (int i = 1; i <= lN-2; i++) 
    {
      dlm_bd->x[istart + i -1] = corner1.x + (double)i * dinc.x; 
      dlm_bd->y[istart + i -1] = corner1.y + (double)i * dinc.y;
      dlm_bd->z[istart + i -1] = corner1.z + (double)i * dinc.z;
    }
  }
#endif
}




bool is_it_in_cube_v2 (coord * u, coord * v, coord * w, coord * mins, 
	coord * maxs, coord * checkpt) 
{
  
  /* a point x with coord (x,y,z) lies in the rectangle if the 3
     following conditions are satisfied */
  /* 1- u.p0 <= u.x <= u.p3 */
  /* 2- v.p0 <= v.x <= v.p1 */
  /* 3- w.p0 <= w.w <= w.p4 */
  double x = checkpt->x;
  double y = checkpt->y;
  double z = checkpt->z;
  double checkval = 0.;
  coord u1 = *u;
  coord v1 = *v;
  coord w1 = *w;
  coord cpt = {0., 0., 0.};
  bool isin = false;
  
  checkval = u1.x*x + u1.y*y + u1.z*z;
  cpt.x = checkval;

  checkval = v1.x*x + v1.y*y + v1.z*z;
  cpt.y = checkval; 

  checkval = w1.x*x + w1.y*y + w1.z*z;
  cpt.z = checkval;

  if ((cpt.x >= maxs->x) && (cpt.x <= mins->x) && (cpt.y >= maxs->y) 
  	&& (cpt.y <= mins->y) && (cpt.z >= maxs->z) && (cpt.z <= mins->z) )
    isin = true;

  return isin;
}




void compute_principal_vectors_Cubes (particle * p) 
{
  GeomParameter * gcp = &(p->g);
  int nfaces = gcp->allFaces;
  int npoints;

  /* get the 3 directions u1, u2, u3 of the cube with such that */
  /* u1 = corner0 - corner3; */
  /* u2 = corner0 - corner1; */
  /* u3 = corner0 - corner4; */
  coord corner0 = {0, 0, 0}, corner1 = {0, 0, 0}, 
  	corner3 = {0, 0, 0}, corner4 = {0, 0, 0};
  coord u1, v1, w1;
  coord mins, maxs;
  
  /* find the coordinates of these 4 corners */
  for (int i = 0; i < nfaces; i++) 
  {
    npoints = gcp->numPointsOnFaces[i];
    
    for (int j = 0; j < npoints; j++) 
    {
      /* printf("index = %lu\n", gcp->cornersIndex[i][j]); */
      if (gcp->cornersIndex[i][j] == 0) 
      {
	long int ii = gcp->cornersIndex[i][j];
	corner0.x = gcp->cornersCoord[ii][0];
	corner0.y = gcp->cornersCoord[ii][1];
	corner0.z = gcp->cornersCoord[ii][2];
      }
      
      if (gcp->cornersIndex[i][j] == 1) 
      {
	long int ii = gcp->cornersIndex[i][j];
	corner1.x = gcp->cornersCoord[ii][0];
	corner1.y = gcp->cornersCoord[ii][1];
	corner1.z = gcp->cornersCoord[ii][2];
      }
      
      if (gcp->cornersIndex[i][j] == 3) 
      {
	long int ii = gcp->cornersIndex[i][j];
	corner3.x = gcp->cornersCoord[ii][0];
	corner3.y = gcp->cornersCoord[ii][1];
	corner3.z = gcp->cornersCoord[ii][2];
      }
      
      if (gcp->cornersIndex[i][j] == 4) 
      {
	long int ii = gcp->cornersIndex[i][j];
	corner4.x = gcp->cornersCoord[ii][0];
	corner4.y = gcp->cornersCoord[ii][1];
	corner4.z = gcp->cornersCoord[ii][2];
      }
    }  
  }

  foreach_dimension() 
  {
    u1.x = corner0.x - corner3.x;
    v1.x = corner0.x - corner1.x;
    w1.x = corner0.x - corner4.x;
    mins.x = 0.;
    maxs.x = 0.;
  }
  
  gcp->u1 = u1;
  gcp->v1 = v1;
  gcp->w1 = w1;
  
  double minval = 0., maxval = 0.;

  foreach_dimension() {
    minval += u1.x*corner0.x;
    maxval += u1.x*corner3.x;
  }

  mins.x = (minval); maxs.x = (maxval);

  minval = 0; maxval = 0;
  foreach_dimension() 
  {
    minval += v1.x*corner0.x;
    maxval += v1.x*corner1.x;
  }
  
  mins.y = (minval); maxs.y = (maxval);

  minval = 0; maxval = 0;
  foreach_dimension() 
  {
    minval += w1.x*corner0.x;
    maxval += w1.x*corner4.x;
  }
  
  mins.z = (minval); maxs.z = (maxval);


  gcp->mins = mins;
  gcp->maxs = maxs;  
}




void create_FD_Boundary_Cube_v2 (GeomParameter * gcp, 
	SolidBodyBoundary * dlm_bd, const int m, const int lN, vector pshift) 
{
  int nfaces = gcp->allFaces;
  int iref, i1, i2, ichoice;

  ichoice = 0;
  int isb = 0;
  int npoints;

  /* Add first interrior points on surfaces */
  for (int i = 0; i < nfaces; i++) 
  {
    npoints = gcp->numPointsOnFaces[i];
    
    iref = gcp->cornersIndex[i][ichoice];
    i1 = gcp->cornersIndex[i][ichoice + 1];
    i2 = gcp->cornersIndex[i][npoints-1];

    coord refcorner = {gcp->cornersCoord[iref][0], gcp->cornersCoord[iref][1],
    	gcp->cornersCoord[iref][2]} ; 

    coord dir1 = {gcp->cornersCoord[i1][0], gcp->cornersCoord[i1][1],
    	gcp->cornersCoord[i1][2]};

    coord dir2 = {gcp->cornersCoord[i2][0], gcp->cornersCoord[i2][1],
    	gcp->cornersCoord[i2][2]};
    
    foreach_dimension() 
    {
      dir1.x -= refcorner.x;
      dir2.x -= refcorner.x;
      dir1.x /= (lN-1);
      dir2.x /= (lN-1);
    }
       
    for (int ii = 1; ii <= lN-2; ii++) 
    {
      for (int jj = 1; jj <= lN-2; jj++) 
      { 
	dlm_bd->x[isb] = refcorner.x + (double) ii * dir1.x 
		+ (double) jj * dir2.x;

	dlm_bd->y[isb] = refcorner.y + (double) ii * dir1.y 
		+ (double) jj * dir2.y;

	dlm_bd->z[isb] = refcorner.z + (double) ii * dir1.z 
		+ (double) jj * dir2.z;
	isb++;
      }
    }
  }

  int allindextable[8][8] = {{0}};
  int j1,jm1;

  
  /* Add points on the edges without the corners*/
  for (int i = 0; i < nfaces; i++) 
  {
    npoints = gcp->numPointsOnFaces[i];
    i1 = gcp->cornersIndex[i][1];

    for (int j = 1; j < npoints; j++) 
    {
      jm1 = gcp->cornersIndex[i][j-1];
      j1 = gcp->cornersIndex[i][j];
      
      if (jm1 > j1) 
      {
	if (allindextable[jm1][j1] == 0) 
	{
	  coord c1 = {gcp->cornersCoord[jm1][0], gcp->cornersCoord[jm1][1], 
	  	gcp->cornersCoord[jm1][2]};
	  coord c2 = {gcp->cornersCoord[j1][0], gcp->cornersCoord[j1][1], 
	  	gcp->cornersCoord[j1][2]};
	  distribute_points_edge_Cube_v2 (c1, c2, dlm_bd, lN, isb);
	  allindextable[jm1][j1] = 1;
	  isb +=lN-2;
	}
      }      
      else 
      {
	if (allindextable[j1][jm1] == 0) 
	{
	  coord c1 = {gcp->cornersCoord[j1][0], gcp->cornersCoord[j1][1], 
	  	gcp->cornersCoord[j1][2]};
	  coord c2 = {gcp->cornersCoord[jm1][0], gcp->cornersCoord[jm1][1], 
	  	gcp->cornersCoord[jm1][2]};
	  distribute_points_edge_Cube_v2 (c1, c2, dlm_bd, lN, isb);
	  allindextable[j1][jm1] = 1;
	  isb +=lN-2;
	}
      }
    }   
  }
  
  /* Add the final 8 corners points */
  for (int i = 0; i  < gcp->ncorners; i++) 
  {
    dlm_bd->x[isb] = gcp->cornersCoord[i][0];
    dlm_bd->y[isb] = gcp->cornersCoord[i][1];
    dlm_bd->z[isb] = gcp->cornersCoord[i][2];
    isb++;
  }
}




void create_FD_Interior_Cube_v2 (particle *p, vector Index_lambda, 
	vector pshift) 
{
  Cache * c;
  
  /** Create the cache of the interior points for a cube*/
  c = &(p->Interior);

  /* compute the 3 principal vector of the cube */
  compute_principal_vectors_Cubes (p);

  /* a point x with coord (x,y,z) lies in the cube if the 3
     following conditions are satisfied */
  /* 1- u.p0 <= u.x <= u.p3 */
  /* 2- v.p0 <= v.x <= v.p1 */
  /* 3- w.p0 <= w.w <= w.p4  */
  coord checkpt;
  coord u1 = p->g.u1;
  coord v1 = p->g.v1;
  coord w1 = p->g.w1;
  coord mins = p->g.mins;
  coord maxs = p->g.maxs;

  /* Min/Max coordinates for the AABB (Axed-Aligned-Bounding-Box) */
  coord mincoord = {HUGE, HUGE, HUGE};
  coord maxcoord = {-HUGE, -HUGE, -HUGE};

  double ** table = p->g.cornersCoord;
  for (int ii = 0; ii < p->g.ncorners; ii++) 
  {
    if (mincoord.x > table[ii][0])
      mincoord.x = table[ii][0];

    if (mincoord.y > table[ii][1])
      mincoord.y = table[ii][1];

    if (mincoord.z > table[ii][2])
      mincoord.z = table[ii][2];

    if (maxcoord.x < table[ii][0])
      maxcoord.x = table[ii][0];

    if (maxcoord.y < table[ii][1])
      maxcoord.y = table[ii][1];

    if (maxcoord.z < table[ii][2])
      maxcoord.z = table[ii][2];
  }

  foreach() 
  {
    checkpt.x = x;
    checkpt.y = y;
    checkpt.z = z;

    /* Check only if the point is in the AABB (Axed-Aligned-Bounding-Box) */
    if ((x > mincoord.x) && (x < maxcoord.x)) 
      if ((y > mincoord.y) && (y < maxcoord.y))
	if ((z > mincoord.z) && (z < maxcoord.z))

	  /* If yes: check if it is inside the cube now */
    	  if (is_it_in_cube_v2 (&u1, &v1, &w1, &mins, &maxs, &checkpt)) 
	  {
	    cache_append (c, point, 0);
	    /* tagg cell with the number of the particle */
	    if ((int)Index_lambda.y[] == -1)
	      Index_lambda.y[] = p->pnum;
	  }
  }
 
  cache_shrink (c);    
}
