/** 
# Helper functions for interfacing Basilisk/Grains3D 
*/


/** Transfers particles velocity into a 2D array to be sent to the 
granular solver */
//----------------------------------------------------------------------------
void UpdateDLMFDtoGS_vel( double arrayv[][6], particle* p, 
	const int m )
//---------------------------------------------------------------------------- 
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




/** Updates particles through parsing and reading a C string coming
from the granular solver */
//----------------------------------------------------------------------------
char* UpdateParticlesBasilisk( char* pstr, const int pstrsize,
	particle* allpart, const int npart_, bool fluidCorrectedAcceleration_, 
	double rhoval_ )
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

  // First entry is the number of particles
  int np = 0;
  sscanf( token, "%d", &np );
  if ( np != npart_ )
    printf ("Error in number of particles in UpdateParticlesBasilisk\n");
  
  // Read the parsed array of character for each particle
  double Ux = 0., Uy = 0., Uz = 0., omx = 0., omy = 0., omz = 0., rhop = 0., 
  	massp  = 0., Ixx = 0., Ixy = 0., Ixz = 0., Iyy = 0., Iyz = 0., Izz = 0.,
	gx = 0., gy = 0., gz = 0., radiusp = 0.,
	MRxx = 0., MRxy = 0., MRxz = 0., MRyx = 0., MRyy = 0., MRyz = 0.,
	MRzx = 0., MRzy = 0., MRzz = 0.;
  int ncornersp = 0;
  char particleTag[3] = "";
  char particleDefaultTagPer[] = "PP";    
  int nperclonesp = 0;
  double* vecx = NULL;
  double* vecy = NULL;
  double* vecz = NULL;  
  
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
    sscanf( token, "%s", particleTag );
        
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

    // If particleTag is "PP", read periodic clone vectors
    if ( !strcmp( particleTag, particleDefaultTagPer ) )
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
#   if dimension == 3
      allpart[k].RotMat[0][0] = MRxx;
      allpart[k].RotMat[0][1] = MRxy;
      allpart[k].RotMat[0][2] = MRxz;
      allpart[k].RotMat[1][0] = MRyx;
      allpart[k].RotMat[1][1] = MRyy;
      allpart[k].RotMat[1][2] = MRyz;
      allpart[k].RotMat[2][0] = MRzx;
      allpart[k].RotMat[2][1] = MRzy;                        
#   else
      allpart[k].RotMat[0][0] = 0.;
      allpart[k].RotMat[0][1] = 0.;
      allpart[k].RotMat[0][2] = 0.;
      allpart[k].RotMat[1][0] = 0.;
      allpart[k].RotMat[1][1] = 0.;
      allpart[k].RotMat[1][2] = 0.;
      allpart[k].RotMat[2][0] = 0.;
      allpart[k].RotMat[2][1] = 0.;
#   endif
    allpart[k].RotMat[2][2] = MRzz;        
    strcpy( allpart[k].tag, particleTag ); 
    gg->center.x = gx;
    gg->center.y = gy;
#   if dimension == 3
      gg->center.z = gz;
#   else
      gg->center.z = 0.;      
#   endif
    gg->ncorners = ncornersp;
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
    //   ( 1 - rhoval / rho_s )
    // otherwise DLMFD_couplingFactor = 1
    allpart[k].DLMFD_couplingfactor = 1. ;
    if ( fluidCorrectedAcceleration_ ) 
      allpart[k].DLMFD_couplingfactor -= rhoval_ / allpart[k].rho_s ;

#   if DLM_Moving_particle      
      // Compute the inverse of the moment of inertia matrix
      compute_inv_inertia( &(allpart[k]) );
#   endif
    
    // Read the additional geometric features of the particle
    // Note that the C function strtok keeps track of the pointer to 
    // the last C string, which explains why we do not need to pass any
    // parameter to the functions below 
    switch ( ncornersp )
    {
#     if dimension == 3
        case 1: 
          allpart[k].shape = SPHERE;
	  update_Sphere( gg ); 
          break;

        // For now, we assume that all 4-corner polyhedrons are tetrahedrons
	case 4: 
          allpart[k].shape = TETRAHEDRON;
	  update_Tetrahedron( gg );
          break;
	  
        // For now, we assume that all 8-corner polyhedrons are cubes
	case 8: 
          allpart[k].shape = CUBE;
	  update_Cube( gg );
          break;
	  
       // For now, we assume that all 12-corner polyhedrons are icosahedrons
       case 12: 
         allpart[k].shape = ICOSAHEDRON;
	 update_Icosahedron( gg );
         break;  
          
       // For now, we assume that all 20-corner polyhedrons are dodecahedrons
       case 20: 
         allpart[k].shape = DODECAHEDRON;
	 update_Dodecahedron( gg );
         break;      	  
#     else
        case 1: 
          allpart[k].shape = CIRCULARCYLINDER2D;
	  update_CircularCylinder2D( gg );
          break;
#     endif	  	  
	        
      default:
        fprintf( stderr, "Unknown ncorners in UpdateParticlesBasilisk!!\n" );
    }                               
  }
  
  return ( pstr );         
}
