/** 
# Helper functions for interfacing Basilisk/Grains3D 
*/


/** Transfers particles velocity into a 2D array to be sent to the 
granular solver */
//----------------------------------------------------------------------------
void UpdateDLMFDtoGS_vel( double** arrayv, RigidBody* allrbs, 
	const int npart )
//---------------------------------------------------------------------------- 
{
  coord U = {0., 0., 0.};
  coord w = {0., 0., 0.};
  
  // We only transfer the particles velocity, i.e. the npart first elements
  // of the allrbs array
  for (size_t k=0;k<npart;k++) 
  {    
#   if TRANSLATION
      U = allrbs[k].U;
#   endif
#   if ROTATION
      w = allrbs[k].w;
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
	RigidBody* allrbs, const size_t nrb_, RigidBody const* allrefrbs,
	double rho_f_, bool init_ )
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
  double Ux = 0., Uy = 0., Uz = 0., omx = 0., omy = 0., omz = 0., 
	gx = 0., gy = 0., gz = 0., 
	MRxx = 0., MRxy = 0., MRxz = 0., MRyx = 0., MRyy = 0., MRyz = 0.,
	MRzx = 0., MRzy = 0., MRzz = 0.,
	bbminx, bbminy, bbminz, bbmaxx, bbmaxy, bbmaxz;
  int pnum = 0;
  size_t gtype = 0;
  char RBTag[3] = "";
  char particleDefaultTag[] = "P";
  char periodicParticleDefaultTag[] = "PP"; 
  char obstacleDefaultTag[] = "O";       
  int nperclonesp = 0;
  double* vecx = NULL;
  double* vecy = NULL;
  double* vecz = NULL; 
  RigidBody const* RBRef = NULL; 
  
  for (size_t k = 0; k < nrb_; k++) 
  { 
    nperclonesp = 0;
    GeomParameter* gg = &(allrbs[k].g);
    gg->pgp = NULL;
    gg->cgp = NULL;    
    allrbs[k].toygsp = NULL;
    
    // Read the rigid body number
    token = strtok( NULL, " " );
    sscanf( token, "%u", &pnum );     
    allrbs[k].pnum = pnum;

    // Read the rigid body geometric type
    token = strtok( NULL, " " );
    sscanf( token, "%lu", &gtype );         
    
    // Read the rigid body type: particle, periodic particle or obstacle 
    token = strtok( NULL, " " );
    sscanf( token, "%s", RBTag );
    strcpy( allrbs[k].typetag, RBTag );     
    if ( !strcmp( RBTag, particleDefaultTag ) )
      allrbs[k].type = PARTICLE;
    else if ( !strcmp( RBTag, periodicParticleDefaultTag ) )
      allrbs[k].type = PERIODICPARTICLE; 
    else if ( !strcmp( RBTag, obstacleDefaultTag ) )
      allrbs[k].type = OBSTACLE;
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

    // Read bbminx
    token = strtok( NULL, " " );
    sscanf( token, "%lf", &bbminx );
    
    // Read bbminy
    token = strtok( NULL, " " );
    sscanf( token, "%lf", &bbminy );
    
#   if dimension == 3
      // Read bbminz
      token = strtok( NULL, " " );
      sscanf( token, "%lf", &bbminz );
#   endif

    // Read bbmaxx
    token = strtok( NULL, " " );
    sscanf( token, "%lf", &bbmaxx );
    
    // Read bbmaxy
    token = strtok( NULL, " " );
    sscanf( token, "%lf", &bbmaxy );
    
#   if dimension == 3
      // Read bbmaxz
      token = strtok( NULL, " " );
      sscanf( token, "%lf", &bbmaxz );
#   endif    
            
    // If RBTag is "PP", read periodic clone vectors
    if ( allrbs[k].type == PERIODICPARTICLE )
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

    // Assign the values and coming from the reference rigid body read to the 
    // rigid body data 
    RBRef = &(allrefrbs[gtype]);
    allrbs[k].geomType = gtype;       
    allrbs[k].rho_s = RBRef->rho_s;
    allrbs[k].M = RBRef->M;
    allrbs[k].Vp = RBRef->Vp; 
#   if dimension == 3
      allrbs[k].RotMat[0][0] = MRxx;
      allrbs[k].RotMat[0][1] = MRxy;
      allrbs[k].RotMat[0][2] = MRxz;
      allrbs[k].RotMat[1][0] = MRyx;
      allrbs[k].RotMat[1][1] = MRyy;
      allrbs[k].RotMat[1][2] = MRyz;
      allrbs[k].RotMat[2][0] = MRzx;
      allrbs[k].RotMat[2][1] = MRzy;                        
#   else
      allrbs[k].RotMat[0][0] = 0.;
      allrbs[k].RotMat[0][1] = 0.;
      allrbs[k].RotMat[0][2] = 0.;
      allrbs[k].RotMat[1][0] = 0.;
      allrbs[k].RotMat[1][1] = 0.;
      allrbs[k].RotMat[1][2] = 0.;
      allrbs[k].RotMat[2][0] = 0.;
      allrbs[k].RotMat[2][1] = 0.;
#   endif
    allrbs[k].RotMat[2][2] = MRzz;        
    gg->center.x = gx;
    gg->center.y = gy;
#   if dimension == 3
      gg->center.z = gz;
#   else
      gg->center.z = 0.;      
#   endif
    gg->BBox.min.x = bbminx;
    gg->BBox.min.y = bbminy;    
    gg->BBox.min.z = bbminz;    
    gg->BBox.max.x = bbmaxx;
    gg->BBox.max.y = bbmaxy;    
    gg->BBox.max.z = bbmaxz;        
    if ( allrbs[k].type != OBSTACLE )
    {
#     if TRANSLATION
        allrbs[k].U.x = Ux;
        allrbs[k].U.y = Uy;	
#       if dimension == 3
          allrbs[k].U.z = Uz;
#       else
          allrbs[k].U.z = 0.;  	  
#       endif       
#     endif
#     if ROTATION
        allrbs[k].w.z = omz;
#       if dimension == 3
          allrbs[k].w.x = omx;
          allrbs[k].w.y = omy;	
#       else
          allrbs[k].w.x = 0.;
          allrbs[k].w.y = 0.;	    
#       endif 	
#     endif

#     if TRANSLATION
        if ( init_ ) allrbs[k].Unm1 = allrbs[k].U;
	foreach_dimension() allrbs[k].imposedU.x = 0.;
#     endif 
#     if ROTATION   	 
        if ( init_ ) 
	{
	  allrbs[k].wnm1 = allrbs[k].w;
#         if B_SPLIT_EXPLICIT_ACCELERATION
#           if dimension == 3
              allrbs[k].Iwnm1.x = allrbs[k].Ip[0] * omx 
	  	+ allrbs[k].Ip[3] * omy + allrbs[k].Ip[4] * omz;
              allrbs[k].Iwnm1.y = allrbs[k].Ip[3] * omx 
	  	+ allrbs[k].Ip[1] * omy + allrbs[k].Ip[5] * omz;
              allrbs[k].Iwnm1.z = allrbs[k].Ip[4] * omx 
	  	+ allrbs[k].Ip[5] * omy + allrbs[k].Ip[2] * omz;
#           else
              allrbs[k].Iwnm1.z = allrbs[k].Ip[2] * omz;
#           endif
#         endif		  
	}
	foreach_dimension() allrbs[k].imposedw.x = 0.;
#     endif
    }
    else
    {   
#     if TRANSLATION
        allrbs[k].imposedU.x = Ux;
        allrbs[k].imposedU.y = Uy;	
#       if dimension == 3
          allrbs[k].imposedU.z = Uz;
#       else
          allrbs[k].imposedU.z = 0.;  	  
#       endif       
#     endif
#     if ROTATION
        allrbs[k].imposedw.z = omz;
#       if dimension == 3
          allrbs[k].imposedw.x = omx;
          allrbs[k].imposedw.y = omy;	
#       else
          allrbs[k].imposedw.x = 0.;
          allrbs[k].imposedw.y = 0.;	    
#       endif 	
#     endif

#     if TRANSLATION
        if ( init_ ) allrbs[k].Unm1 = allrbs[k].imposedU;
	foreach_dimension() allrbs[k].U.x = 0.;
#     endif 
#     if ROTATION   	 
        if ( init_ ) allrbs[k].wnm1 = allrbs[k].imposedw;
	foreach_dimension() allrbs[k].w.x = 0.; 
#     endif 
    }
     
    gg->ncorners = RBRef->g.ncorners;
    gg->radius = RBRef->g.radius; 
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
    
    allrbs[k].DLMFD_couplingfactor = RBRef->DLMFD_couplingfactor ;

    // In case rigid bodies are treated as fixed obstacles
    if ( RIGIDBODIES_AS_FIXED_OBSTACLES )
    {
      strcpy( allrbs[k].typetag, obstacleDefaultTag ); 
      allrbs[k].type = OBSTACLE; 
      allrbs[k].imposedU.x = 0.;
      allrbs[k].imposedU.y = 0.;	
#     if dimension == 3
        allrbs[k].imposedU.z = 0.;
        allrbs[k].imposedw.x = 0.;
        allrbs[k].imposedw.y = 0.;
#     else
        allrbs[k].imposedU.z = 0.;
        allrbs[k].imposedw.x = 0.;
        allrbs[k].imposedw.y = 0.;		  
#     endif
      allrbs[k].imposedw.z = 0.; 
#     if TRANSLATION
        foreach_dimension() 
	{
	  allrbs[k].U.x = 0.;
	  allrbs[k].Unm1.x = 0.;
	}
#     endif 
#     if ROTATION   	 
        foreach_dimension() 
	{
	  allrbs[k].w.x = 0.;
	  allrbs[k].wnm1.x = 0.; 
	}
#     endif           
    }
      
    // Compute the inverse of the moment of inertia matrix
    if ( allrbs[k].type != OBSTACLE )
      compute_inertia_inv_inertia( &(allrbs[k]), RBRef );
    
    // Set the additional geometric features of the rigid body
    allrbs[k].shape = RBRef->shape;
    switch ( allrbs[k].shape )
    {
      case SPHERE: 
	update_Sphere_from_RBRef( gg, RBRef ); 
        break;

      case TETRAHEDRON: 
	update_Polyhedron_from_RBRef( gg, RBRef, allrbs[k].RotMat );
        break;
	  
      case CUBE: 
	update_Polyhedron_from_RBRef( gg, RBRef, allrbs[k].RotMat );
        break;
	
      case OCTAHEDRON: 
	update_Polyhedron_from_RBRef( gg, RBRef, allrbs[k].RotMat );
        break;	
	  
      case ICOSAHEDRON: 
	update_Polyhedron_from_RBRef( gg, RBRef, allrbs[k].RotMat );
        break;  
          
      case DODECAHEDRON: 
	update_Polyhedron_from_RBRef( gg, RBRef, allrbs[k].RotMat );
        break;
	 
      case BOX: 
	update_Polyhedron_from_RBRef( gg, RBRef, allrbs[k].RotMat );
        break;	
	 
      case CIRCULARCYLINDER3D: 
	update_CircularCylinder3D_from_RBRef( gg, RBRef, allrbs[k].RotMat );
        break;

      case CONE: 
	update_Cone_from_RBRef( gg, RBRef, allrbs[k].RotMat );
        break;
	 
      case TRUNCATEDCONE: 
        update_TruncatedCone_from_RBRef( gg, RBRef, allrbs[k].RotMat );
        break;	 	        	  

      case CIRCULARCYLINDER2D: 
	update_CircularCylinder2D_from_RBRef( gg, RBRef );
        break;  	  

      case ELLIPSOID: 
	update_Ellipsoid_from_RBRef( gg, RBRef, allrbs[k].RotMat );
        break;
	
      case HEXAGONALPRISM:
	update_Polyhedron_from_RBRef( gg, RBRef, allrbs[k].RotMat );
        break;
	
      case SIXBRANCHSTAR:
	update_Star_from_RBRef( gg, RBRef, allrbs[k].RotMat );
        break;	      	 
	        
      default:
        printf( "Unknown shape in UpdateParticlesBasilisk!!\n" );
    }                               
  }
  
  return ( pstr );        
}




/** Create the array of reference rigid bodies */
//----------------------------------------------------------------------------
char* CreateReferenceRBBasilisk( char* pstr, const int pstrsize,
	RigidBody* allrefrbs, const size_t nrefrb_, double rho_f_ )
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

  // First entry is the number of reference rigid bodies
  size_t np = 0;
  sscanf( token, "%lu", &np );
  if ( np != nrefrb_ )
    printf ("Error in number of reference rigid bodies in "
    	"CreateReferenceRBBasilisk\n");
  
  // Read the parsed array of character for each rigid body
  double rhop = 0., massp  = 0., Ixx = 0., Ixy = 0., Ixz = 0., Iyy = 0., 
  	Iyz = 0., Izz = 0., gx = 0., gy = 0., gz = 0., radiusp = 0.,
        MRxx = 0., MRxy = 0., MRxz = 0., MRyx = 0., MRyy = 0., MRyz = 0.,
 	MRzx = 0., MRzy = 0., MRzz = 0.;
  int nshapecode = 0;
  size_t geomType = 0;
  
  for (size_t k = 0; k < nrefrb_; k++) 
  { 
    GeomParameter* gg = &(allrefrbs[k].g);
    gg->pgp = NULL;
    gg->cgp = NULL;    
    allrefrbs[k].toygsp = NULL;
    
    // Read the rigid body geometric type
    token = strtok( NULL, " " );
    sscanf( token, "%lu", &geomType );     
    allrefrbs[k].geomType = geomType;
    if ( k != geomType )
      printf ("Error in geometric type numbering of reference rigid bodies in "
    	"CreateReferenceRBBasilisk\n");    

    // Read the rigid body's shape code
    token = strtok( NULL, " " );
    sscanf( token, "%d", &nshapecode ); 
    
    // Set the rigid body type to reference
    allrefrbs[k].type = REFERENCERIGIDBODY; 
    strcpy( allrefrbs[k].typetag, "REF" );                 
    
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

    // Read radius
    token = strtok( NULL, " " );
    sscanf( token, "%lf", &radiusp );

    // Assign the values read to the rigid body data    
    allrefrbs[k].rho_s = rhop;
    allrefrbs[k].M = massp;
    allrefrbs[k].Vp = (allrefrbs[k].M)/(allrefrbs[k].rho_s); 
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
      allrefrbs[k].Ip[0] = Ixx;
      allrefrbs[k].Ip[1] = Iyy;
      allrefrbs[k].Ip[3] = Ixy;
      allrefrbs[k].Ip[4] = Ixz;
      allrefrbs[k].Ip[5] = Iyz;      
#   else
      allrefrbs[k].Ip[0] = 0.;
      allrefrbs[k].Ip[1] = 0.;
      allrefrbs[k].Ip[3] = 0.;
      allrefrbs[k].Ip[4] = 0.;
      allrefrbs[k].Ip[5] = 0.;
#   endif
    allrefrbs[k].Ip[2] = Izz;
#   if dimension == 3
      allrefrbs[k].RotMat[0][0] = MRxx;
      allrefrbs[k].RotMat[0][1] = MRxy;
      allrefrbs[k].RotMat[0][2] = MRxz;
      allrefrbs[k].RotMat[1][0] = MRyx;
      allrefrbs[k].RotMat[1][1] = MRyy;
      allrefrbs[k].RotMat[1][2] = MRyz;
      allrefrbs[k].RotMat[2][0] = MRzx;
      allrefrbs[k].RotMat[2][1] = MRzy;                        
#   else
      allrefrbs[k].RotMat[0][0] = 0.;
      allrefrbs[k].RotMat[0][1] = 0.;
      allrefrbs[k].RotMat[0][2] = 0.;
      allrefrbs[k].RotMat[1][0] = 0.;
      allrefrbs[k].RotMat[1][1] = 0.;
      allrefrbs[k].RotMat[1][2] = 0.;
      allrefrbs[k].RotMat[2][0] = 0.;
      allrefrbs[k].RotMat[2][1] = 0.;
#   endif
    allrefrbs[k].RotMat[2][2] = MRzz;      
    gg->center.x = gx;
    gg->center.y = gy;
#   if dimension == 3
      gg->center.z = gz;
#   else
      gg->center.z = 0.;      
#   endif    
      
    if ( nshapecode < 10000 ) // Disc, sphere, polygon and polyhedron
      gg->ncorners = (int)( (double)nshapecode / 100. );
    else if ( nshapecode == 10000 ) // 2D Box, not implemented
      gg->ncorners = 4;
    else if ( nshapecode == 10001 ) // 3D Box
      gg->ncorners = 8;       
    else gg->ncorners = 0; // Any other shape
    gg->radius = radiusp; 
            
    
    // DLMFD coupling factor
    // If B_SPLIT_EXPLICIT_ACCELERATION == false, DLMFD_couplingFactor = 
    //   ( 1 - FLUID_DENSITY / rho_s )
    // otherwise DLMFD_couplingFactor = 1
    allrefrbs[k].DLMFD_couplingfactor = 1. ;
#   if !B_SPLIT_EXPLICIT_ACCELERATION
      allrefrbs[k].DLMFD_couplingfactor -= rho_f_ / allrefrbs[k].rho_s ;
#   endif
    
    // Read the additional geometric features of the rigid body
    // Note that the C function strtok keeps track of the pointer to 
    // the last C string, which explains why we do not need to pass any
    // parameter to the functions below 
    switch ( nshapecode )
    {
#     if dimension == 3
        case 1: 
          allrefrbs[k].shape = SPHERE;
	  read_reference_Sphere( gg, allrefrbs[k].RotMat ); 
          break;

        // For now, we assume that all 4-corner/4-face polyhedrons are regular 
	// tetrahedrons
	case 404: 
          allrefrbs[k].shape = TETRAHEDRON;
	  read_reference_Tetrahedron( gg, allrefrbs[k].RotMat );
          break;
	  
        // For now, we assume that all 6-corner/8-face polyhedrons are regular 
	// octahedrons
	case 608: 
          allrefrbs[k].shape = OCTAHEDRON;
	  read_reference_Octahedron( gg, allrefrbs[k].RotMat );
          break;	  
	  
        // For now, we assume that all 8-corner/6-face polyhedrons are cubes
	case 806: 
          allrefrbs[k].shape = CUBE;
	  read_reference_Cube( gg, allrefrbs[k].RotMat );
          break;
	  
        // For now, we assume that all 12-corner/20-face polyhedrons are regular
	// icosahedrons 
        case 1220: 
	  allrefrbs[k].shape = ICOSAHEDRON;
	  read_reference_Icosahedron( gg, allrefrbs[k].RotMat );
          break;
	  
        // For now, we assume that all 12-corner/8-face polyhedrons are 
	// hexagonal prism (8 faces) 
	case 1208:
	  allrefrbs[k].shape = HEXAGONALPRISM;
	  read_reference_HexagonalPrism( gg, allrefrbs[k].RotMat );
          break;  
          
        // For now, we assume that all 20-corner/12-face polyhedrons are regular
	// dodecahedrons
        case 2012: 
          allrefrbs[k].shape = DODECAHEDRON;
	  read_reference_Dodecahedron( gg, allrefrbs[k].RotMat );
          break;

        case 10005: 
          allrefrbs[k].shape = ELLIPSOID;
	  read_reference_Ellipsoid( gg, allrefrbs[k].RotMat );
          break;
	 
        case 10001: 
          allrefrbs[k].shape = BOX;
	  read_reference_Box( gg, allrefrbs[k].RotMat );
          break;	
	 
        case 10002: 
          allrefrbs[k].shape = CIRCULARCYLINDER3D;
	  read_reference_CircularCylinder3D( gg, allrefrbs[k].RotMat );
          break;

        case 10003: 
          allrefrbs[k].shape = CONE;
	  read_reference_Cone( gg, allrefrbs[k].RotMat );
          break;
	 
        case 10004: 
          allrefrbs[k].shape = TRUNCATEDCONE;
	  read_reference_TruncatedCone( gg, allrefrbs[k].RotMat );
          break;
	 
        case 1000003: 
          allrefrbs[k].shape = SIXBRANCHSTAR;
	  read_reference_Star( gg, allrefrbs[k].RotMat );
          break;	  	 	        	  
#     else
        case 2: 
          allrefrbs[k].shape = CIRCULARCYLINDER2D;
	  read_reference_CircularCylinder2D( gg, allrefrbs[k].RotMat );
          break;
#     endif	  	  
	        
      default:
        fprintf( stderr, "Unknown ncorners in CreateReferenceRBBasilisk!!\n" );
    }
    
    // Once we defined all rigid bodies in their neutral angular position and
    // centered at the origin, we reset the center of mass position and rotation
    // matrix to identity
    foreach_dimension() gg->center.x = 0.;
    for (size_t i=0;i<3;++i)
      for (size_t j=0;j<3;++j) 
        if ( i == j )
          allrefrbs[k].RotMat[i][j] = 1.;
	else
	  allrefrbs[k].RotMat[i][j] = 0.;
  }

  return ( pstr );         
}
