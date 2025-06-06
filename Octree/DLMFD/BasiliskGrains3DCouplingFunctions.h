/** 
# Helper functions for interfacing Basilisk/Grains3D 
*/


/** Transfers particles velocity into a 2D array to be sent to the 
granular solver */
//----------------------------------------------------------------------------
void UpdateDLMFDtoGS_vel( double** arrayv, RigidBody* allrb, 
	const int npart )
//---------------------------------------------------------------------------- 
{
  coord U = {0., 0., 0.};
  coord w = {0., 0., 0.};
  
  // We only transfer the particles velocity, i.e. the npart first elements
  // of the allrb array
  for (size_t k=0;k<npart;k++) 
  {    
#   if TRANSLATION
      U = allrb[k].U;
#   endif
#   if ROTATION
      w = allrb[k].w;
#   endif
 
    arrayv[k][0] = U.x;
    arrayv[k][1] = U.y;
    arrayv[k][2] = U.z;
    arrayv[k][3] = w.x;
    arrayv[k][4] = w.y;
    arrayv[k][5] = w.z;   
  }
}




/** Updates rigid bodies through parsing and reading a C string coming
from the granular solver */
//----------------------------------------------------------------------------
char* UpdateParticlesBasilisk( char* pstr, const int pstrsize,
	RigidBody* allrb, const size_t nrb_, 
	bool fluidCorrectedAcceleration_, double rho_f_ )
//----------------------------------------------------------------------------
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

  // First entry is the number of rigid bodies
  size_t np = 0;
  sscanf( token, "%lu", &np );
  if ( np != nrb_ )
    printf ("Error in number of rigid bodies in UpdateParticlesBasilisk\n");
  
  // Read the parsed array of character for each rigid body
  double Ux = 0., Uy = 0., Uz = 0., omx = 0., omy = 0., omz = 0., rhop = 0., 
  	massp  = 0., Ixx = 0., Ixy = 0., Ixz = 0., Iyy = 0., Iyz = 0., Izz = 0.,
	gx = 0., gy = 0., gz = 0., radiusp = 0.,
	MRxx = 0., MRxy = 0., MRxz = 0., MRyx = 0., MRyy = 0., MRyz = 0.,
	MRzx = 0., MRzy = 0., MRzz = 0.;
  int ncornersp = 0;
  char RBTag[3] = "";
  char particleDefaultTag[] = "P";
  char periodicParticleDefaultTag[] = "PP"; 
  char obstacleDefaultTag[] = "O";       
  int nperclonesp = 0;
  double* vecx = NULL;
  double* vecy = NULL;
  double* vecz = NULL;  
  
  for (size_t k = 0; k < nrb_; k++) 
  { 
    nperclonesp = 0;
    GeomParameter* gg = &(allrb[k].g);
    gg->pgp = NULL;
    gg->cgp = NULL;    
    allrb[k].toygsp = NULL;
    
    // Read the rigid body number but assign k
    token = strtok( NULL, " " );
    allrb[k].pnum = k;

    // Read the rigid body's number of corners or code
    token = strtok( NULL, " " );
    sscanf( token, "%d", &ncornersp ); 
    
    // Read the rigid body type: particle, periodic particle or obstacle 
    token = strtok( NULL, " " );
    sscanf( token, "%s", RBTag );
    strcpy( allrb[k].typetag, RBTag );     
    if ( !strcmp( RBTag, particleDefaultTag ) )
      allrb[k].type = PARTICLE;
    else if ( !strcmp( RBTag, periodicParticleDefaultTag ) )
      allrb[k].type = PERIODICPARTICLE; 
    else if ( !strcmp( RBTag, obstacleDefaultTag ) )
      allrb[k].type = OBSTACLE;
    else
      printf( "Warning: Unknown rigid body type in the string sent by Grains3D"
      	"\n" );          
            
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

#   if dimension == 3
      // Read MRxx
      token = strtok( NULL, " " );
      sscanf( token, "%lf", &MRxx ); 

      // Read MRxy
      token = strtok( NULL, " " );
      sscanf( token, "%lf", &MRxy );
      
      // Read MRxz
      token = strtok( NULL, " " );
      sscanf( token, "%lf", &MRxz );       

      // Read MRyx
      token = strtok( NULL, " " );
      sscanf( token, "%lf", &MRyx );
      
      // Read MRyy
      token = strtok( NULL, " " );
      sscanf( token, "%lf", &MRyy );      
      
      // Read MRyz
      token = strtok( NULL, " " );
      sscanf( token, "%lf", &MRyz );      
      
      // Read MRz
      token = strtok( NULL, " " );
      sscanf( token, "%lf", &MRzx );
      
      // Read MRzy
      token = strtok( NULL, " " );
      sscanf( token, "%lf", &MRzy );       
#   endif

    // Read MRzz
    token = strtok( NULL, " " );
    sscanf( token, "%lf", &MRzz );  

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

    // If RBTag is "PP", read periodic clone vectors
    if ( allrb[k].type == PERIODICPARTICLE )
    {
      // Read the number of clones
      token = strtok( NULL, " " );
      sscanf( token, "%d", &nperclonesp);
      vecx = (double*) malloc( nperclonesp * sizeof(double) );
      vecy = (double*) malloc( nperclonesp * sizeof(double) );
      vecz = (double*) malloc( nperclonesp * sizeof(double) );
      
      // Read the periodic clone vectors
      for (size_t j = 0; j < nperclonesp; j++)
      {
        // Read vecx
	token = strtok( NULL, " " );
	sscanf( token, "%lf", &vecx[j] );

	// Read vecy
	token = strtok( NULL, " " );
	sscanf( token, "%lf", &vecy[j] );

#       if dimension == 3
	  // Read vecz
	 token = strtok( NULL, " " );
	 sscanf( token, "%lf", &vecz[j] );
#        endif
      }                 
    }

    // Read radius
    token = strtok( NULL, " " );
    sscanf( token, "%lf", &radiusp );

    // Assign the values read to the rigid body data
    if ( allrb[k].type != OBSTACLE )
    {
      // Save previous velocity before updating
#     if TRANSLATION
        allrb[k].Unm1 = allrb[k].U;
	foreach_dimension() allrb[k].imposedU.x = 0.;
#     endif 
#     if ROTATION   	 
        allrb[k].wnm1 = allrb[k].w;
	foreach_dimension() allrb[k].imposedw.x = 0.;
#     endif

#     if TRANSLATION
        allrb[k].U.x = Ux;
        allrb[k].U.y = Uy;	
#       if dimension == 3
          allrb[k].U.z = Uz;
#       else
          allrb[k].U.z = 0.;  	  
#       endif       
#     endif
#     if ROTATION
        allrb[k].w.z = omz;
#       if dimension == 3
          allrb[k].w.x = omx;
          allrb[k].w.y = omy;	
#       else
          allrb[k].w.x = 0.;
          allrb[k].w.y = 0.;	    
#       endif 	
#     endif
    }
    else
    {
#     if TRANSLATION
        allrb[k].Unm1 = allrb[k].imposedU;
	foreach_dimension() allrb[k].U.x = 0.;
#     endif 
#     if ROTATION   	 
        allrb[k].wnm1 = allrb[k].imposedw;
	foreach_dimension() allrb[k].w.x = 0.; 
#     endif 
   
#     if TRANSLATION
        allrb[k].imposedU.x = Ux;
        allrb[k].imposedU.y = Uy;	
#       if dimension == 3
          allrb[k].imposedU.z = Uz;
#       else
          allrb[k].imposedU.z = 0.;  	  
#       endif       
#     endif
#     if ROTATION
        allrb[k].imposedw.z = omz;
#       if dimension == 3
          allrb[k].imposedw.x = omx;
          allrb[k].imposedw.y = omy;	
#       else
          allrb[k].imposedw.x = 0.;
          allrb[k].imposedw.y = 0.;	    
#       endif 	
#     endif
    }
    
    allrb[k].rho_s = rhop;
    allrb[k].M = massp;
    allrb[k].Vp = (allrb[k].M)/(allrb[k].rho_s); 
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
      allrb[k].Ip[0] = Ixx;
      allrb[k].Ip[1] = Iyy;
      allrb[k].Ip[3] = Ixy;
      allrb[k].Ip[4] = Ixz;
      allrb[k].Ip[5] = Iyz;      
#   else
      allrb[k].Ip[0] = 0.;
      allrb[k].Ip[1] = 0.;
      allrb[k].Ip[3] = 0.;
      allrb[k].Ip[4] = 0.;
      allrb[k].Ip[5] = 0.;
#   endif
    allrb[k].Ip[2] = Izz;
#   if dimension == 3
      allrb[k].RotMat[0][0] = MRxx;
      allrb[k].RotMat[0][1] = MRxy;
      allrb[k].RotMat[0][2] = MRxz;
      allrb[k].RotMat[1][0] = MRyx;
      allrb[k].RotMat[1][1] = MRyy;
      allrb[k].RotMat[1][2] = MRyz;
      allrb[k].RotMat[2][0] = MRzx;
      allrb[k].RotMat[2][1] = MRzy;                        
#   else
      allrb[k].RotMat[0][0] = 0.;
      allrb[k].RotMat[0][1] = 0.;
      allrb[k].RotMat[0][2] = 0.;
      allrb[k].RotMat[1][0] = 0.;
      allrb[k].RotMat[1][1] = 0.;
      allrb[k].RotMat[1][2] = 0.;
      allrb[k].RotMat[2][0] = 0.;
      allrb[k].RotMat[2][1] = 0.;
#   endif
    allrb[k].RotMat[2][2] = MRzz;        
    gg->center.x = gx;
    gg->center.y = gy;
#   if dimension == 3
      gg->center.z = gz;
#   else
      gg->center.z = 0.;      
#   endif
    gg->ncorners = ncornersp;
    if ( gg->ncorners == 666 ) gg->ncorners = 8;
    else if ( gg->ncorners == 777 ) gg->ncorners = 0;
    gg->radius = radiusp; 
    gg->nperclones = nperclonesp;    
    if ( nperclonesp )
    {
	gg->perclonecenters = (coord*) malloc( nperclonesp * sizeof(coord) );
	for (int j=0; j < nperclonesp; j++)
 	{
          gg->perclonecenters[j].x = gx + vecx[j];
 	  gg->perclonecenters[j].y = gy + vecy[j];
#         if dimension == 3
	    gg->perclonecenters[j].z = gz + vecz[j];
#         else
	    gg->perclonecenters[j].z = 0.;
#         endif
   	}

	// Free the vecx pointers
	free(vecx);
	free(vecy);
	free(vecz);      
    }            
    
    // DLMFD coupling factor
    // If fluidCorrectedAcceleration_ == true, DLMFD_couplingFactor = 
    //   ( 1 - FLUID_DENSITY / rho_s )
    // otherwise DLMFD_couplingFactor = 1
    allrb[k].DLMFD_couplingfactor = 1. ;
    if ( fluidCorrectedAcceleration_ ) 
      allrb[k].DLMFD_couplingfactor -= rho_f_ / allrb[k].rho_s ;

    // In case rigid bodies are treated as fixed obstacles
    if ( RIGIDBODIES_AS_FIXED_OBSTACLES )
    {
      strcpy( allrb[k].typetag, obstacleDefaultTag ); 
      allrb[k].type = OBSTACLE; 
      allrb[k].imposedU.x = 0.;
      allrb[k].imposedU.y = 0.;	
#     if dimension == 3
        allrb[k].imposedU.z = 0.;
        allrb[k].imposedw.x = 0.;
        allrb[k].imposedw.y = 0.;
#     else
        allrb[k].imposedU.z = 0.;
        allrb[k].imposedw.x = 0.;
        allrb[k].imposedw.y = 0.;		  
#     endif
      allrb[k].imposedw.z = 0.; 
#     if TRANSLATION
        foreach_dimension() 
	{
	  allrb[k].U.x = 0.;
	  allrb[k].Unm1.x = 0.;
	}
#     endif 
#     if ROTATION   	 
        foreach_dimension() 
	{
	  allrb[k].w.x = 0.;
	  allrb[k].wnm1.x = 0.; 
	}
#     endif           
    }
      
    // Compute the inverse of the moment of inertia matrix
    if ( allrb[k].type != OBSTACLE )
      compute_inv_inertia( &(allrb[k]) );
    
    // Read the additional geometric features of the rigid body
    // Note that the C function strtok keeps track of the pointer to 
    // the last C string, which explains why we do not need to pass any
    // parameter to the functions below 
    switch ( ncornersp )
    {
#     if dimension == 3
        case 1: 
          allrb[k].shape = SPHERE;
	  update_Sphere( gg ); 
          break;

        // For now, we assume that all 4-corner polyhedrons are tetrahedrons
	case 4: 
          allrb[k].shape = TETRAHEDRON;
	  update_Tetrahedron( gg );
          break;
	  
        // For now, we assume that all 8-corner polyhedrons are cubes
	case 8: 
          allrb[k].shape = CUBE;
	  update_Cube( gg );
          break;
	  
       // For now, we assume that all 12-corner polyhedrons are icosahedrons
       case 12: 
         allrb[k].shape = ICOSAHEDRON;
	 update_Icosahedron( gg );
         break;  
          
       // For now, we assume that all 20-corner polyhedrons are dodecahedrons
       case 20: 
         allrb[k].shape = DODECAHEDRON;
	 update_Dodecahedron( gg );
         break;
	 
       case 666: 
         allrb[k].shape = BOX;
	 update_Box( gg );
         break;	
	 
       case 777: 
         allrb[k].shape = CIRCULARCYLINDER3D;
	 update_CircularCylinder3D( gg );
         break;
	 
       case 8888: 
         allrb[k].shape = TRUNCATEDCONE;
	 update_TruncatedCone( gg );
         break;	 	        	  
#     else
        case 1: 
          allrb[k].shape = CIRCULARCYLINDER2D;
	  update_CircularCylinder2D( gg );
          break;
#     endif	  	  
	        
      default:
        fprintf( stderr, "Unknown ncorners in UpdateParticlesBasilisk!!\n" );
    }                               
  }
  
  return ( pstr );         
}
