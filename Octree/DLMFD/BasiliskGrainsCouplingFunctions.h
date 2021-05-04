/** 
# Helper functions for interfacing Basilisk/Grains3D 
*/
void printBasiliskDataStructure (struct BasiliskDataStructure * b) 
{
  printf("Grains nparticles = %lu\n", b->particlessize);
  
  /* Geometric aspect */
  printf("Grains radius = %f\n", b->rayon);
  printf("Grains center.x = %f\n", b->centreX);
  printf("Grains center.y = %f\n", b->centreY);
  printf("Grains center.z = %f\n", b->centreZ);
  printf("Grains ncorners = %d\n", b->ncorners);


  /* Coordinates of the corners, number of faces and their numerotation */

  printf("Grains allPoints = %lu\n", b->allPoints);
  printf("Grains allFaces = %lu\n", b->allFaces);

  for (int i = 0; i < (b->allPoints); i ++) {
# if dimension == 2
    printf ("Grains Coordinates of point %d is (%f,%f)\n", i, 
    	b->PointsCoord[i][0],b->PointsCoord[i][1]);
    
# elif dimension ==3
    printf ("Grains Coordinates of point %d is (%f,%f,%f)\n", i, 
    	b->PointsCoord[i][0], b->PointsCoord[i][1], b->PointsCoord[i][2]);
# endif
  } 
   
  /* Velocities */
  printf("Grains U.x = %f\n", b->vitesseTX);
  printf("Grains U.y = %f\n", b->vitesseTY);
  printf("Grains U.z = %f\n", b->vitesseTZ);
  
  printf("Grains w.x = %f\n", b->vitesseRX);
  printf("Grains w.y = %f\n", b->vitesseRY);
  printf("Grains w.z = %f\n", b->vitesseRZ);

  /* Physical properties */
  printf ("Grains mass = %f\n", b->masse);
  printf ("Grains rho_s = %f\n", b->masseVol);
  for (int k = 0; k < 6; k++) {
      printf ("Grains inertia = %g\n", b->inertie[k]);
    }
}




void UpdateParticlesBasilisk( struct BasiliskDataStructure * b, particle * p, 
	const int m, bool explicit_added_mass_, double rhoval_ )
{
  int r;
  for (int k = 0; k < m; k++) {
    
#if _MPI

    MPI_Bcast (&(b[k].particlessize), 1, MPI_INT, 0, MPI_COMM_WORLD);

    /* Geometric aspect */
    MPI_Bcast (&(b[k].centreX), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&(b[k].centreY), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&(b[k].centreZ), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&(b[k].rayon), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&(b[k].ncorners), 1, MPI_INT, 0, MPI_COMM_WORLD);

    /* Coordinates of the corners, number of faces and their numerotation */
    MPI_Bcast (&(b[k].allPoints), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&(b[k].allFaces), 1, MPI_INT, 0, MPI_COMM_WORLD);

    r = b[k].allPoints;
  
    /* Allocate the structure in the other threads (Grains is serial) */
    if (pid() > 0) {
      b[k].PointsCoord = (double **) malloc (r*sizeof(double *));
#if dimension == 2 
      for (int i = 0; i < r; i++) 
	b[k].PointsCoord[i] = (double *) malloc (2 * sizeof(double));
#elif dimension == 3
      for (int i = 0; i < r; i++) 
	b[k].PointsCoord[i] = (double *) malloc (3 * sizeof(double));
#endif
    }

    for (int i = 0; i < r; i++) {
      MPI_Bcast (&(b[k].PointsCoord[i][0]), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast (&(b[k].PointsCoord[i][1]), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#if dimension ==3
      MPI_Bcast (&(b[k].PointsCoord[i][2]), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    }

    /* Allocate the structures in the other threads (Grains is serial) */
    r = b[k].allFaces;
    if (pid() > 0) {
      b[k].FacesIndex = (long int **) malloc (r * sizeof(long int *));
      b[k].numPointsOnFaces = (long int *) malloc (r * sizeof(long int));
    }

    /* Broadcast-allocate-broadcast */
    for (int i = 0; i < r; i++) {
      MPI_Bcast (&(b[k].numPointsOnFaces[i]), 1, MPI_LONG, 0, MPI_COMM_WORLD);
      int rr = b[k].numPointsOnFaces[i];
      if (pid() > 0)
	b[k].FacesIndex[i] = (long int *)malloc(rr * sizeof(long int));
      MPI_Barrier (MPI_COMM_WORLD);

      for (int j = 0; j < rr; j++)
	MPI_Bcast (&(b[k].FacesIndex[i][j]), 1, MPI_LONG, 0, MPI_COMM_WORLD);
    }

    /* Velocities */
    MPI_Bcast (&(b[k].vitesseTX), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&(b[k].vitesseTY), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&(b[k].vitesseTZ), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&(b[k].vitesseRX), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&(b[k].vitesseRY), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&(b[k].vitesseRZ), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* Physical properties */
    MPI_Bcast (&(b[k].masseVol), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&(b[k].masse), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (int i = 0; i < 6; i++)
      MPI_Bcast (&(b[k].inertie[i]), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

 
    /* Geometric aspect */ 
    GeomParameter gg; 
    coord cc = {0., 0., 0.};
    cc.x = b[k].centreX;
    cc.y = b[k].centreY;
    cc.z = b[k].centreZ;
    gg.center = cc;
    gg.radius = b[k].rayon;
    gg.ncorners = b[k].ncorners;
    gg.pgp = (PolyGeomParameter*) malloc (sizeof(PolyGeomParameter));
    gg.pgp->allPoints = b[k].allPoints;
    gg.pgp->allFaces = b[k].allFaces;

    /* We need the coordinates of the corners for a cube */

    /* Allocate on all threads */
    r = gg.ncorners;
    gg.pgp->cornersCoord = (double **) malloc (r*sizeof(double *));
    
#if dimension == 2 
    for (int i = 0; i < r; i++) 
      gg.pgp->cornersCoord[i] = (double *) malloc (2 * sizeof(double));
#elif dimension == 3
    for (int i = 0; i < r; i++) 
      gg.pgp->cornersCoord[i] = (double *) malloc (3 * sizeof(double));
#endif

    int ndim = 0;
    for (int i = 0; i < r; i++) {
#if dimension == 2 
      ndim  = 2; 
#elif dimension == 3
      ndim = 3;
#endif
      for (int j = 0; j < ndim; j++) {
	  gg.pgp->cornersCoord[i][j] = b[k].PointsCoord[i][j];
      }
    }

    /* We need the indices of the corners for a cube */
    /* Allocate */
    r = gg.pgp->allFaces;
    
    gg.pgp->cornersIndex = (long int **) malloc (r * sizeof(long int *));
    gg.pgp->numPointsOnFaces = (long int *) malloc (r * sizeof(long int));
    
    for (int i = 0; i < r; i++) {
      int rr = b[k].numPointsOnFaces[i];
      gg.pgp->cornersIndex[i] = (long int *) malloc (rr * sizeof(long int));
      gg.pgp->numPointsOnFaces[i] = rr;
      for (int j = 0; j < rr; j++) {
	gg.pgp->cornersIndex[i][j] = b[k].FacesIndex[i][j];
      }
    }
    
    p[k].g = gg;

    // Rigid body shape: only cube or sphere so far
    if ( gg.ncorners == 8 ) 
    {
      p[k].shape = CUBE;
      compute_principal_vectors_Cubes( &(p[k]) );
    }
    else
      p[k].shape = SPHERE;
    
#if DLM_Moving_particle    
    /* Velocities */
    cc.x = b[k].vitesseTX;
    cc.y = b[k].vitesseTY;
    cc.z = b[k].vitesseTZ;

#if TRANSLATION
    /* Save the previous translational velocity before updating */
    p[k].Unm1 = p[k].U;
    p[k].U = cc;
#endif
    cc.x = b[k].vitesseRX;
    cc.y = b[k].vitesseRY;
    cc.z = b[k].vitesseRZ;

#if ROTATION
    /* Save the previous angular velocity before updating */
    p[k].wnm1 = p[k].w;
    p[k].w = cc;
#endif

    /* Particle number */
    p[k].pnum = k;
    
    /* Physical properties */
    p[k].M = b[k].masse;
    p[k].rho_s = b[k].masseVol;
    p[k].Vp = (p[k].M)/(p[k].rho_s);

    /* Inertia tensor: Grains stores them as */
    /* inertie[0] = Ixx; */
    /* inertie[1] = Ixy; */
    /* inertie[2] = Ixz; */
    /* inertie[3] = Iyy; */
    /* inertie[4] = Iyz; */
    /* inertie[5] = Izz; */

    /* Basilisk stores these as */
    /* Ip[0] = Ixx */
    /* Ip[1] = Iyy */
    /* Ip[2] = Izz */
    /* Ip[3] = Ixy */
    /* Ip[4] = Ixz */
    /* Ip[5] = Iyz */
      
    /* Ixx */
    p[k].Ip[0] = b[k].inertie[0];
    /* Iyy */
    p[k].Ip[1] = b[k].inertie[3];
    /* Izz */
    p[k].Ip[2] = b[k].inertie[5];
    /* Ixy */
    p[k].Ip[3] = b[k].inertie[1];
    /* Ixz */
    p[k].Ip[4] = b[k].inertie[2];
    /* Iyz */
    p[k].Ip[5] = b[k].inertie[4];
    
    // DLMFD coupling factor
    // If explicit add mass, DLMFD_couplingFactor = 1
    // otherwise DLMFD_couplingFactor = ( 1 - rhoval / rho_s )
    p[k].DLMFD_couplingfactor = 1. ;
    if ( !explicit_added_mass_ ) 
      p[k].DLMFD_couplingfactor -= rhoval_ / p[k].rho_s ;
      
    // Compute the inverse of the moment of inertia matrix
    compute_inv_inertia( &(p[k]) );           
#endif
  }
}




void UpdateBasiliskStructure (struct BasiliskDataStructure * b, particle * p, 
	const int m) 
{
  coord U = {0., 0., 0.};
  coord w = {0., 0., 0.};
  
  for (int k = 0; k < m; k++) {
#if DLM_Moving_particle
#if TRANSLATION
    U = p[k].U;
#endif
#if ROTATION
    w = p[k].w;
#endif
#endif
    b[k].vitesseTX = U.x;
    b[k].vitesseTY = U.y;
    b[k].vitesseTZ = U.z;
    b[k].vitesseRX = w.x;
    b[k].vitesseRY = w.y;
    b[k].vitesseRZ = w.z;   
  }
}




void UpdateDLMFDtoGS_vel( double arrayv[][6], particle* p, 
	const int m ) 
{
  coord U = {0., 0., 0.};
  coord w = {0., 0., 0.};
  
  for (size_t k=0;k<m;k++) 
  {
#   if DLM_Moving_particle
#     if TRANSLATION
        U = p[k].U;
#     endif
#     if ROTATION
        w = p[k].w;
#     endif
#   else
      U = p[k].imposedU;
      w = p[k].imposedw;
#   endif
    arrayv[k][0] = U.x;
    arrayv[k][1] = U.y;
    arrayv[k][2] = U.z;
    arrayv[k][3] = w.x;
    arrayv[k][4] = w.y;
    arrayv[k][5] = w.z;   
  }
}




void unallocateBasiliskDataStructure(struct BasiliskDataStructure * b, 
	const int m) 
{

  double * c;
  
  for (int k = 0; k < m; k++) {
    /* Unallocate fields of BasiliskDataStructure that are no longer needed 
    (we transferred everything on the particles structure */
    /* printf("addresse of b[k].PointsCoord = %p\n", 
    	(void *) b[k].PointsCoord); */
    
    for (signed long int i = 0; i < b[k].allPoints; i++) {
      c = b[k].PointsCoord[i];
      /* printf("addresse of b[k].PointsCoord[i] = %p\n",
      	(void *) b[k].PointsCoord[i]); */
      /* printf("c = %p\n",(void *) c); */
      if(c)
	free(c);
    }

    /* printf("ok free(b[k].PointsCoord[i]) on thread %d\n",pid()); */
    free(b[k].PointsCoord);

    /* printf("ok free(b[k].PointsCoord) on thread %d\n",pid()); */
    
    for (signed long int i = 0; i < b[k].allFaces; i++) {
      free(b[k].FacesIndex[i]);
    }
    /* printf("ok free(b[k].FaceIndex[i]) on thread %d\n",pid()); */
    free(b[k].FacesIndex);
    free(b[k].numPointsOnFaces);
  }

}




char* UpdateParticlesBasilisk2( char* pstr, const int pstrsize,
	particle* allpart, const int npart_, bool explicit_added_mass_, 
	double rhoval_ )
{
# if _MPI
    // Broadcast the size of the array of characters
    int sstr = pstrsize;
    MPI_Bcast( &sstr, 1, MPI_INT, 0, MPI_COMM_WORLD );
    
    // Allocate the array of characters of other processes than 0
    if ( pid() != 0 )
      pstr = (char*) malloc( sstr * sizeof(char) ); 
    
    // Broadcast the array of characters
    MPI_Bcast( pstr, sstr, MPI_CHAR, 0, MPI_COMM_WORLD );
# endif

  char* token = NULL;
  
  // Parse the array of character coming from Grains3D
  token = strtok( pstr, " " );  

  // First entry is the number of particles
  int np = 0;
  sscanf( token, "%d", &np );
  if ( np != npart_ )
    printf ("Error in number of particles in UpdateParticlesBasilisk2\n");
  
  // Read the parsed array of character for each particle
  double Ux = 0., Uy = 0., Uz = 0., omx = 0., omy = 0., omz = 0., rhop = 0., 
  	massp  = 0., Ixx = 0., Ixy = 0., Ixz = 0., Iyy = 0., Iyz = 0., Izz = 0.,
	gx = 0., gy = 0., gz = 0., radiusp = 0.;
  int ncornersp = 0;
  for (size_t k = 0; k < npart_; k++) 
  { 
    GeomParameter* gg = &(allpart[k].g);
    gg->pgp = NULL;
#   if DLM_Moving_particle    
      allpart[k].toygsp = NULL;
#   endif
    
    // Read the particle number but assign k
    token = strtok( NULL, " " );
    allpart[k].pnum = k;

    // Read the particle's number of corners or tag
    token = strtok( NULL, " " );
    sscanf( token, "%d", &ncornersp ); 
    
    // Read the particle type: standard, periodic or obstacle (not used for now)
    token = strtok( NULL, " " );
    
    // Read Ux
    token = strtok( NULL, " " );
    sscanf( token, "%lf", &Ux );
    
    // Read Uy
    token = strtok( NULL, " " );
    sscanf( token, "%lf", &Uy );    
    
#   if dimension == 3
      // Read Uz
      token = strtok( NULL, " " );
      sscanf( token, "%lf", &Uz ); 
      
      // Read omega_x
      token = strtok( NULL, " " );
      sscanf( token, "%lf", &omx ); 
      
      // Read omega_y
      token = strtok( NULL, " " );
      sscanf( token, "%lf", &omy );                 
#   endif 

    // Read omega_z
    token = strtok( NULL, " " );
    sscanf( token, "%lf", &omz );
    
    // Read density
    token = strtok( NULL, " " );
    sscanf( token, "%lf", &rhop );
    
    // Read mass
    token = strtok( NULL, " " );
    sscanf( token, "%lf", &massp );
    
#   if dimension == 3
      // Read Ixx
      token = strtok( NULL, " " );
      sscanf( token, "%lf", &Ixx ); 

      // Read Ixy
      token = strtok( NULL, " " );
      sscanf( token, "%lf", &Ixy );
      
      // Read Ixz
      token = strtok( NULL, " " );
      sscanf( token, "%lf", &Ixz );       

      // Read Iyy
      token = strtok( NULL, " " );
      sscanf( token, "%lf", &Iyy ); 

      // Read Iyz
      token = strtok( NULL, " " );
      sscanf( token, "%lf", &Iyz );
#   endif  
      
    // Read Izz
    token = strtok( NULL, " " );
    sscanf( token, "%lf", &Izz );                 

    // Read gx
    token = strtok( NULL, " " );
    sscanf( token, "%lf", &gx );
    
    // Read gy
    token = strtok( NULL, " " );
    sscanf( token, "%lf", &gy );
    
#   if dimension == 3
      // Read gz
      token = strtok( NULL, " " );
      sscanf( token, "%lf", &gz );
#   endif                    

    // Read radius
    token = strtok( NULL, " " );
    sscanf( token, "%lf", &radiusp );


    // Assign the values read to the particle data
#   if DLM_Moving_particle 
      // Save previous velocity before updating
#     if TRANSLATION
        allpart[k].Unm1 = allpart[k].U;
#     endif 
#     if ROTATION   	 
        allpart[k].wnm1 = allpart[k].w; 
#     endif

#     if TRANSLATION
        allpart[k].U.x = Ux;
        allpart[k].U.y = Uy;	
#       if dimension == 3
          allpart[k].U.z = Uz;
#       else
          allpart[k].U.z = 0.;  	  
#       endif       
#     endif
#     if ROTATION
        allpart[k].w.z = omz;
#       if dimension == 3
          allpart[k].w.x = omx;
          allpart[k].w.y = omy;	
#       else
          allpart[k].w.x = 0.;
          allpart[k].w.y = 0.;	    
#       endif 	
#     endif
#   else
      allpart[k].imposedU.x = Ux;
      allpart[k].imposedU.y = Uy;	
#     if dimension == 3
        allpart[k].imposedU.z = Uz;
        allpart[k].imposedw.x = omx;
        allpart[k].imposedw.y = omy;
#     else
        allpart[k].imposedU.z = 0.;
        allpart[k].imposedw.x = 0.;
        allpart[k].imposedw.y = 0.;		  
#     endif
      allpart[k].imposedw.z = omz; 
# endif
    allpart[k].rho_s = rhop;
    allpart[k].M = massp;
    allpart[k].Vp = (allpart[k].M)/(allpart[k].rho_s); 
    /* Inertia tensor: Grains stores them as */
    /* inertie[0] = Ixx; */
    /* inertie[1] = Ixy; */
    /* inertie[2] = Ixz; */
    /* inertie[3] = Iyy; */
    /* inertie[4] = Iyz; */
    /* inertie[5] = Izz; */
    /* Basilisk stores these as */
    /* Ip[0] = Ixx */
    /* Ip[1] = Iyy */
    /* Ip[2] = Izz */
    /* Ip[3] = Ixy */
    /* Ip[4] = Ixz */
    /* Ip[5] = Iyz */
#   if dimension == 3
      allpart[k].Ip[0] = Ixx;
      allpart[k].Ip[1] = Iyy;
      allpart[k].Ip[3] = Ixy;
      allpart[k].Ip[4] = Ixz;
      allpart[k].Ip[5] = Iyz;      
#   else
      allpart[k].Ip[0] = 0.;
      allpart[k].Ip[1] = 0.;
      allpart[k].Ip[3] = 0.;
      allpart[k].Ip[4] = 0.;
      allpart[k].Ip[5] = 0.;
#   endif
    allpart[k].Ip[2] = Izz;
    gg->center.x = gx;
    gg->center.y = gy;
#   if dimension == 3
      gg->center.z = gz;
#   else
      gg->center.z = 0.;      
#   endif
    gg->ncorners = ncornersp;
    gg->radius = radiusp;             
    
    // DLMFD coupling factor
    // If explicit add mass, DLMFD_couplingFactor = 1
    // otherwise DLMFD_couplingFactor = ( 1 - rhoval / rho_s )
    allpart[k].DLMFD_couplingfactor = 1. ;
    if ( !explicit_added_mass_ ) 
      allpart[k].DLMFD_couplingfactor -= rhoval_ / allpart[k].rho_s ;

#   if DLM_Moving_particle      
      // Compute the inverse of the moment of inertia matrix
      compute_inv_inertia( &(allpart[k]) );
#   endif
    
    // Read the additional geometric features of the particle
    switch ( ncornersp )
    {
#     if dimension == 3
        case 1: 
          allpart[k].shape = SPHERE;
	  update_Sphere( gg ); 
          break;
	  
        // For now, we assume that all 8-corner polyhedrons are cubes
	case 8: 
          allpart[k].shape = CUBE;
	  update_Cube( gg );
	  compute_principal_vectors_Cubes( &(allpart[k]) ); 
          break;	  
#     else
        case 1: 
          allpart[k].shape = CIRCULARCYLINDER2D;
	  update_CircularCylinder2D( gg );
          break;
#     endif	  	  
	        
      default:
        fprintf( stderr,"Unknown ncorners in UpdateParticlesBasilisk2!!\n" );
    }                               
  }
  
  return ( pstr );         
}
