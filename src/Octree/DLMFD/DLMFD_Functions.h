/**
# General functions for the DLMFD implementation.
*/

/** Define NSDF, the number of significant digits after the decimal point
to output data in files in text mode. NSDF cannot be lower than 3, larger 
than 15 and 9 (for formatting reasons). */
# ifndef NSDF
#   define NSDF 7
# else 
#   if ( NSDF == 9 )
#     undef NSDF
#     define NSDF 8
#   else
#     if ( NSDF < 3 )
#       undef NSDF
#       define NSDF 3
#     else
#       if ( NSDF > 15 )
#         undef NSDF
#         define NSDF 15
#       endif
#     endif
#   endif
# endif




/** Define the factor alpha (generally between 1 and 2) that is involved 
in the inter-boundary point distance on the rigid body surface. */
# ifndef INTERBPCOEF
#   define INTERBPCOEF 2.
# endif 




/** Different rigid body shapes supported */      
enum RigidBodyShape {
  SPHERE,
  CIRCULARCYLINDER2D,
  CUBE,
  TETRAHEDRON,
  OCTAHEDRON,
  DODECAHEDRON,
  ICOSAHEDRON,
  BOX,
  CIRCULARCYLINDER3D,
  CONE,
  TRUNCATEDCONE,
  ELLIPSOID
};




/** Different rigid body shapes supported */      
enum RigidBodyType {
  PARTICLE,
  PERIODICPARTICLE,
  OBSTACLE, 
  REFERENCERIGIDBODY
};




/** Structure for an axis-aligned bounding box */
typedef struct {
  coord min;
  coord max;
} AABB;


 

/** Structure for the rigid body boundary points */
typedef struct {
  coord* bp;
  coord* normal;
  bool* deactivated;
  int m;
} RigidBodyBoundary;




/** Additional geometric parameters for polygons/polyhedrons */
typedef struct {
  int allPoints, allFaces;
  double** cornersCoord;
  long int** cornersIndex;
  long int* numPointsOnFaces;
} PolyGeomParameter;




/** Additional geometric parameters for 3D cylinders */
typedef struct {
  coord BottomCenter;
  coord TopCenter;
  coord BottomToTopVec;
  coord RadialRefVec;
  double radius;
  double height;
} CylGeomParameter;




/** Additional geometric parameters for full or truncated cones */
typedef struct {
  coord BottomCenter;
  coord TopCenter;
  coord BottomToTopVec;
  coord BottomRadialRefVec;
  coord TopRadialRefVec;  
  double BottomRadius;
  double TopRadius;  
  double height;
} TruncConeGeomParameter;




/** Additional geometric parameters for ellpsoids */
typedef struct { 
  double a,b,c;
  double n1, n2;
} EllipsoidGeomParameter;




/** Rigid body geometric parameters */
typedef struct {
  coord center;
  AABB BBox;
  coord* perclonecenters;  
  double radius;
  int ncorners;
  int nperclones;
  PolyGeomParameter* pgp;
  CylGeomParameter* cgp; 
  TruncConeGeomParameter* tcgp;
  EllipsoidGeomParameter* elgp;  
} GeomParameter;




/** Rigid body parameters for the toy granular solver */
typedef struct {
  double kn, en, vzero, wished_ratio;
  coord normalvector;
  GeomParameter gnm1;  
} ToyGSParameter;




/** Set of parameters describing a rigid body */
typedef struct {
  size_t pnum;
  size_t geomType;
  char typetag[4];
  enum RigidBodyType type;
  enum RigidBodyShape shape;  
  RigidBodyBoundary s;
  GeomParameter g;
  double M, Ip[6], rho_s, Vp, DLMFD_couplingfactor, RotMat[3][3];  
  ToyGSParameter *toygsp;
  double Ip_inv[3][3];
  coord addforce;    
# if TRANSLATION
    coord U, Unm1, splitUacc, qU, tU, imposedU;    
# endif
# if ROTATION
    coord w, wnm1, Iwnm1, splitwacc, qw, tw, imposedw;
# endif
  Cache Interior;
  Cache Boundary;
} RigidBody;




/** Total number of DLMFD cells and points for statistics */
typedef struct {
  size_t total_number_of_DLMFDcells; 
  size_t total_number_of_DLMFDpts;  
} DLMFDptscells;




# define DYNARRAYBLOCK 128

/** Structure for a dynamic array of unsigned integers */
typedef struct {
  size_t n;
  size_t nm;   
  size_t* elem; 
} dynUIarray;




/** Initialize and allocate a dynUIarray */
//----------------------------------------------------------------------------
void initialize_and_allocate_dynUIarray( dynUIarray* a, const int nm ) 
//----------------------------------------------------------------------------
{
  a->elem = (size_t*) calloc( nm, sizeof(size_t) );
  a->n = 0;
  a->nm = nm;
}




/** Free a dynUIarray */
//----------------------------------------------------------------------------
void free_dynUIarray( dynUIarray* a ) 
//----------------------------------------------------------------------------
{
  free(a->elem); a->elem = NULL;
  a->n = 0;
  a->nm = 0;
}




/** Append a new element to a dynUIarray */
//----------------------------------------------------------------------------
void append_dynUIarray( dynUIarray* a, size_t newelem )
//----------------------------------------------------------------------------
{
  if ( a->n >= a->nm )
  {
    a->nm += DYNARRAYBLOCK;
    a->elem = (size_t*) realloc( a->elem, a->nm*sizeof(size_t) );
  }      
  a->elem[a->n] = newelem;
  a->n += 1;    
}




/** Structure for a dynamic array of pointers to double */
typedef struct {
  size_t n;
  size_t nm;   
  double** elem; 
} dynPDBarray;




/** Initialize and allocate a dynPDBarray */
//----------------------------------------------------------------------------
void initialize_and_allocate_dynPDBarray( dynPDBarray* a, const int nm ) 
//----------------------------------------------------------------------------
{
  a->elem = (double**) calloc( nm, sizeof(double*) );
  a->n = 0;
  a->nm = nm;
}




/** Free a dynPDBarray */
//----------------------------------------------------------------------------
void free_dynPDBarray( dynPDBarray* a ) 
//----------------------------------------------------------------------------
{
  free(a->elem); a->elem = NULL;
  a->n = 0;
  a->nm = 0;
}




/** Append a new element to a dynPDBarray */
//----------------------------------------------------------------------------
void append_dynPDBarray( dynPDBarray* a, double* newelem )
//----------------------------------------------------------------------------
{
  if ( a->n >= a->nm )
  {
    a->nm += DYNARRAYBLOCK;
    a->elem = (double**) realloc( a->elem, a->nm*sizeof(double*) );
  }      
  a->elem[a->n] = newelem;
  a->n += 1;    
}




/** Force synchronization of a field in parallel by setting the dirty flag
of the field to true */
//----------------------------------------------------------------------------
trace void synchronize( scalar* list )
//----------------------------------------------------------------------------
{
  for (scalar s in list)
    s.dirty = true;
  boundary(list);
}




/** Allocates memory for m points in the RigidBodyBoundary structure. */
//----------------------------------------------------------------------------
void allocate_RigidBodyBoundary( RigidBodyBoundary* sbm, const int m ) 
//----------------------------------------------------------------------------
{
  sbm->bp = (coord*) calloc( m, sizeof(coord) );
  sbm->normal = (coord*) calloc( m, sizeof(coord) );   
  sbm->deactivated = (bool*) calloc( m, sizeof(bool) );
  sbm->m = m;
}




/** Re-allocates memory for m points in the RigidBodyBoundary structure. */
//----------------------------------------------------------------------------
void reallocate_RigidBodyBoundary( RigidBodyBoundary* sbm, const int m ) 
//----------------------------------------------------------------------------
{
  sbm->bp = (coord*) realloc( sbm->bp, m * sizeof(coord) ); 
  sbm->normal = (coord*) realloc( sbm->normal, m * sizeof(coord) );   
  sbm->deactivated = (bool*) realloc( sbm->deactivated, m * sizeof(bool) );       
  sbm->m = m;
}




/** Frees memory associated to the points in the RigidBodyBoundary structure. */
//----------------------------------------------------------------------------
void free_RigidBodyBoundary( RigidBodyBoundary* sbm ) 
//----------------------------------------------------------------------------
{
  free( sbm->bp ); sbm->bp = NULL;
  free( sbm->normal ); sbm->normal = NULL;  
  free( sbm->deactivated ); sbm->deactivated = NULL;
  sbm->m = 0;  
}




/** Allocates an initial Cache structure */
//----------------------------------------------------------------------------
void initialize_and_allocate_Cache( Cache* p ) 
//----------------------------------------------------------------------------
{
  p->n = 0;
  p->nm = 1;
  p->p = (Index*) calloc( 1, sizeof(Index) );
}




/** Inverts a 3 x 3 matrix and stores it in a double** */
//----------------------------------------------------------------------------
void inverse3by3matrix( const double Matrix[3][3], double** inversedMatrix ) 
//----------------------------------------------------------------------------
{
  double block1 = Matrix[1][1] * Matrix[2][2] 
  	- Matrix[1][2] * Matrix[2][1];
  double block2 = Matrix[1][2] * Matrix[2][0] 
  	- Matrix[1][0] * Matrix[2][2];
  double block3 = Matrix[1][0] * Matrix[2][1] 
  	- Matrix[1][1] * Matrix[2][0];

  double determinant = Matrix[0][0] * block1 + Matrix[0][1] * block2 
  	+ Matrix[0][2] * block3;

  /* m[i][j], *(*(m + i) + j) */
  
  *(*(inversedMatrix + 0) + 0) = block1 / determinant;
  *(*(inversedMatrix + 0) + 1) = ( Matrix[0][2] * Matrix[2][1] 
  	- Matrix[0][1] * Matrix[2][2] ) / determinant;
  *(*(inversedMatrix + 0) + 2) = ( Matrix[0][1] * Matrix[1][2] 
  	- Matrix[0][2] * Matrix[1][1]) / determinant; 
  *(*(inversedMatrix + 1) + 0) = block2 / determinant;  
  *(*(inversedMatrix + 1) + 1) = ( Matrix[0][0] * Matrix[2][2] 
  	- Matrix[0][2] * Matrix[2][0] ) / determinant;
  *(*(inversedMatrix + 1) + 2) = ( Matrix[0][2] * Matrix[1][0] 
  	- Matrix[0][0] * Matrix[1][2] ) / determinant;
  *(*(inversedMatrix + 2) + 0) = block3 / determinant;
  *(*(inversedMatrix + 2) + 1) = ( Matrix[0][1] * Matrix[2][0] 
  	- Matrix[0][0] * Matrix[2][1]) / determinant;
  *(*(inversedMatrix + 2) + 2) = ( Matrix[0][0] * Matrix[1][1] 
  	- Matrix[0][1] * Matrix[1][0] ) / determinant ;
}




/** Inverts a 3 x 3 matrix and stores it in a double[][] */
//----------------------------------------------------------------------------
void inverse3by3matrix__( const double Matrix[3][3], 
	double inversedMatrix[3][3] ) 
//----------------------------------------------------------------------------
{
  double block1 = Matrix[1][1] * Matrix[2][2] 
  	- Matrix[1][2] * Matrix[2][1];
  double block2 = Matrix[1][2] * Matrix[2][0] 
  	- Matrix[1][0] * Matrix[2][2];
  double block3 = Matrix[1][0] * Matrix[2][1] 
  	- Matrix[1][1] * Matrix[2][0];

  double determinant = Matrix[0][0] * block1 + Matrix[0][1] * block2 
  	+ Matrix[0][2] * block3;
  
  inversedMatrix[0][0] = block1 / determinant;
  inversedMatrix[0][1] = ( Matrix[0][2] * Matrix[2][1] 
  	- Matrix[0][1] * Matrix[2][2]) / determinant;
  inversedMatrix[0][2] = ( Matrix[0][1] * Matrix[1][2] 
  	- Matrix[0][2] * Matrix[1][1] ) / determinant; 
  inversedMatrix[1][0] = block2 / determinant;
  inversedMatrix[1][1] = ( Matrix[0][0] * Matrix[2][2] 
  	- Matrix[0][2] * Matrix[2][0] ) / determinant;
  inversedMatrix[1][2] = ( Matrix[0][2] * Matrix[1][0] 
  	- Matrix[0][0] * Matrix[1][2] ) / determinant;
  inversedMatrix[2][0] = block3 / determinant;
  inversedMatrix[2][1] = ( Matrix[0][1] * Matrix[2][0] 
  	- Matrix[0][0] * Matrix[2][1] ) / determinant;
  inversedMatrix[2][2] = ( Matrix[0][0] * Matrix[1][1] 
  	- Matrix[0][1] * Matrix[1][0]) / determinant ;
}




/** 3 x 3 matrix - vector dot product where the matrix is stored as a 
double[][] and the vector as a coord, the result is stored as a coord */
//----------------------------------------------------------------------------
void matCoordDotProduct( const double Matrix[3][3], const coord v, coord* res ) 
//----------------------------------------------------------------------------
{
  (*res).x = Matrix[0][0] * v.x + Matrix[0][1] * v.y + Matrix[0][2] * v.z;
  (*res).y = Matrix[1][0] * v.x + Matrix[1][1] * v.y + Matrix[1][2] * v.z;
  (*res).z = Matrix[2][0] * v.x + Matrix[2][1] * v.y + Matrix[2][2] * v.z;
}




/** 3 x 3 matrix - vector dot product where the matrix is stored as a 
double[][] and the vector as a double*, the result is stored as a double* */
//----------------------------------------------------------------------------
void matVecDotProduct( const double Matrix[3][3], double const* v, double* res ) 
//----------------------------------------------------------------------------
{
  res[0] = Matrix[0][0] * v[0] + Matrix[0][1] * v[1] + Matrix[0][2] * v[2];
  res[1] = Matrix[1][0] * v[0] + Matrix[1][1] * v[1] + Matrix[1][2] * v[2];
  res[2] = Matrix[2][0] * v[0] + Matrix[2][1] * v[1] + Matrix[2][2] * v[2];
}




/** 3 x 3 matrix - vector dot product where the matrix is stored as a 
double[][] and the vector as a coord, the result is stored as a coord */
//----------------------------------------------------------------------------
void matTransposedCoordDotProduct( const double Matrix[3][3], const coord v, 
	coord* res ) 
//----------------------------------------------------------------------------
{
  (*res).x = Matrix[0][0] * v.x + Matrix[1][0] * v.y + Matrix[2][0] * v.z;
  (*res).y = Matrix[0][1] * v.x + Matrix[1][1] * v.y + Matrix[2][1] * v.z;
  (*res).z = Matrix[0][2] * v.x + Matrix[1][2] * v.y + Matrix[2][2] * v.z;
}




/** 3 x 3 matrix transposed - vector dot product where the matrix is stored as 
a double[][] and the vector as a coord, the result is stored as a double* */
//----------------------------------------------------------------------------
void matTransposedVecDotProduct( const double Matrix[3][3], double const* v, 
	double* res ) 
//----------------------------------------------------------------------------
{
  res[0] = Matrix[0][0] * v[0] + Matrix[1][0] * v[1] + Matrix[2][0] * v[2];
  res[1] = Matrix[0][1] * v[0] + Matrix[1][1] * v[1] + Matrix[2][1] * v[2];
  res[2] = Matrix[0][2] * v[0] + Matrix[1][2] * v[1] + Matrix[2][2] * v[2];
}




/** 3 x 3 matrix - 3 x 3 matrix dot product where both matrices are stored as
double[][], the result is stored as a double[3][3] */
//----------------------------------------------------------------------------
void matMatDotProduct( const double A[3][3], const double B[3][3], 
	double Res[3][3] ) 
//----------------------------------------------------------------------------
{
  for (size_t i=0;i<3;++i)
    for (size_t j=0;j<3;++j)
    {
      Res[i][j] = 0.;
      for (size_t k=0;k<3;++k) Res[i][j] += A[i][k] * B[k][j];
    }    
}




/** 3 x 3 matrix - 3 x 3 matrix transposed dot product where both matrices are 
stored as double[][], the result is stored as a double[3][3] */
//----------------------------------------------------------------------------
void matMatTransposedDotProduct( const double A[3][3], const double B[3][3], 
	double Res[3][3] ) 
//----------------------------------------------------------------------------
{
  for (size_t i=0;i<3;++i)
    for (size_t j=0;j<3;++j)
    {
      Res[i][j] = 0.;
      for (size_t k=0;k<3;++k) Res[i][j] += A[i][k] * B[j][k];
    }    
}




/** Vec - Vec cross product where both vectors are stored as coord, the 
result is stored as a coord */
//----------------------------------------------------------------------------
void VecVecCrossProduct( coord const u, coord const v, coord* res )
//----------------------------------------------------------------------------
{
  (*res).x = u.y * v.z - u.z * v.y;
  (*res).y = u.z * v.x - u.x * v.z;
  (*res).z = u.x * v.y - u.y * v.x; 
}




/** Vec - Vec dot product where both vectors are stored as coord */
//----------------------------------------------------------------------------
double VecVecDotProduct( coord const u, coord const v )
//----------------------------------------------------------------------------
{
  double res = 0.;
  
  foreach_dimension() res += u.x * v.x; 
  
  return ( res );
}




/** Computes the inverse of the moment of inertia matrix of a rigid body */
//----------------------------------------------------------------------------
void compute_inv_inertia( RigidBody* p )
//----------------------------------------------------------------------------
{
  /* The inertia tensor is */
  /*  Ixx  Ixy  Ixz */
  /*  Iyx  Iyy  Iyz */
  /*  Izx  Izy  Izz */ 
  /* with */
  /* Ip[0] = Ixx */
  /* Ip[1] = Iyy */
  /* Ip[2] = Izz */
  /* Ip[3] = Ixy */
  /* Ip[4] = Ixz */
  /* Ip[5] = Iyz */

  // Transfer Ip to a 3x3 matrix
  double Imat[3][3]; 
  Imat[0][0] = p->Ip[0];
  Imat[1][1] = p->Ip[1];
  Imat[2][2] = p->Ip[2];

  Imat[0][1] = p->Ip[3];
  Imat[0][2] = p->Ip[4];
  Imat[1][2] = p->Ip[5];
    
  Imat[1][0] = Imat[0][1]; 
  Imat[2][0] = Imat[0][2];
  Imat[2][1] = Imat[1][2];
  
  // Invert Imat and copy the result in p->Ip_inv
  // Ip_inv is a 2D 3 by 3 array
  inverse3by3matrix__( Imat, p->Ip_inv );  
}




/** Computes the inertia matrix and the inverse of the moment of inertia matrix
 of a rigid body in the world frame of reference */
//----------------------------------------------------------------------------
void compute_inertia_inv_inertia( RigidBody* p, RigidBody const* rbref )
//----------------------------------------------------------------------------
{
  /* The inertia tensor is */
  /*  Ixx  Ixy  Ixz */
  /*  Iyx  Iyy  Iyz */
  /*  Izx  Izy  Izz */ 
  /* with */
  /* Ip[0] = Ixx */
  /* Ip[1] = Iyy */
  /* Ip[2] = Izz */
  /* Ip[3] = Ixy */
  /* Ip[4] = Ixz */
  /* Ip[5] = Iyz */
  
  // Transfer I in the body frame of reference to a 3x3 matrix
  double Ib[3][3]; 
  Ib[0][0] = rbref->Ip[0];
  Ib[1][1] = rbref->Ip[1];
  Ib[2][2] = rbref->Ip[2];

  Ib[0][1] = rbref->Ip[3];
  Ib[0][2] = rbref->Ip[4];
  Ib[1][2] = rbref->Ip[5];
    
  Ib[1][0] = Ib[0][1]; 
  Ib[2][0] = Ib[0][2];
  Ib[2][1] = Ib[1][2];
  
  // Compute I = RotMat * Ib * RotMat^t in the world frame of reference
  double MIb[3][3];
  matMatDotProduct( p->RotMat, Ib, MIb );
  double MIbMt[3][3];  
  matMatTransposedDotProduct( MIb, p->RotMat, MIbMt );
  
  // Transfer to p->I
  p->Ip[0] = MIbMt[0][0];
  p->Ip[1] = MIbMt[1][1];
  p->Ip[2] = MIbMt[2][2];
  p->Ip[3] = MIbMt[0][1];
  p->Ip[4] = MIbMt[0][2];
  p->Ip[5] = MIbMt[1][2];
		   
  // Invert MIbMt and copy the result in p->Ip_inv
  // Ip_inv is a 2D 3 by 3 array
  inverse3by3matrix__( MIbMt, p->Ip_inv );  
}




/** Assigns the shifted coordinates of a AABB to another AABB */
//----------------------------------------------------------------------------
void assign_shifted_BBox( AABB* shifted, AABB const* original, 
	coord const shift )
//----------------------------------------------------------------------------
{
  foreach_dimension()
  {
    shifted->min.x = original->min.x + shift.x;
    shifted->max.x = original->max.x + shift.x;    
  }
}




/** Returns whether two AABB intersect */
//----------------------------------------------------------------------------
bool intersect( AABB const* A, AABB const* B )
//----------------------------------------------------------------------------
{
  return ( ( A->min.x <= B->max.x && A->max.x >= B->min.x ) &&
	( A->min.y <= B->max.y && A->max.y >= B->min.y ) &&
	( A->min.z <= B->max.z && A->max.z >= B->min.z ) );
}




/** Computes the local domain AABB */
//----------------------------------------------------------------------------
void compute_local_domain_AABB( AABB* A )
//----------------------------------------------------------------------------
{
  foreach_dimension()
  {
    A->min.x = 1.e20;
    A->max.x = -1.e20;
  }
  foreach(serial)
  {
    A->min.x = min( x, A->min.x );
    A->min.y = min( y, A->min.y );    
    A->max.x = max( x, A->max.x );
    A->max.y = max( y, A->max.y );    
#   if dimension == 3
      A->max.z = max( z, A->max.z );
      A->min.z = min( z, A->min.z ); 
#   endif           
  }
}

 
# include "CircularCylinder2D.h"
# include "Sphere.h"
# include "Cube.h"
# include "Tetrahedron.h"
# include "Octahedron.h"
# include "Dodecahedron.h"
# include "Icosahedron.h"
# include "Box.h"
# include "CircularCylinder3D.h"
# include "TruncatedCone.h"
# include "Cone.h"
# include "Ellipsoid.h"

/** Frees the rigid body data that were dynamically allocated */
//----------------------------------------------------------------------------
void free_rigidbodies( RigidBody* allrbs, const size_t nrb, bool full_free ) 
//----------------------------------------------------------------------------
{
  for (size_t k=0;k<nrb;k++) 
  {            
    // Free the boundary point coordinate arrays
    RigidBodyBoundary* sbm = &(allrbs[k].s);
    free_RigidBodyBoundary( sbm );

    // Free the caches 
    Cache* c = &(allrbs[k].Interior);
    free( c->p );
    c->p = NULL;
    c = &(allrbs[k].Boundary);
    free( c->p );
    c->p = NULL;    

    if ( full_free )
    {
      // Free the periodic clones position vector
      if ( allrbs[k].g.nperclones )
      { 
        free( allrbs[k].g.perclonecenters );
	allrbs[k].g.perclonecenters = NULL;
      }

      // Free the toy granular solver parameter structure
      if ( allrbs[k].toygsp )
      {
        free( allrbs[k].toygsp );
        allrbs[k].toygsp = NULL;
      }     

      // Free the additional geometric features of the rigid body
      switch ( allrbs[k].shape )
      {
        case SPHERE:
          free_Sphere( &(allrbs[k].g) );
	  break;
	  
        case CIRCULARCYLINDER2D:
          free_CircularCylinder2D( &(allrbs[k].g) );
	  break;
	  
        case CUBE:
          free_Polyhedron( &(allrbs[k].g) );
	  break;

        case TETRAHEDRON:
          free_Polyhedron( &(allrbs[k].g) );
	  break;
	
        case OCTAHEDRON:
	  free_Polyhedron( &(allrbs[k].g) );
	  break;
	
        case ICOSAHEDRON:
	  free_Polyhedron( &(allrbs[k].g) );
	  break;

        case DODECAHEDRON:
	  free_Polyhedron( &(allrbs[k].g) );
	  break;
	  
        case BOX:
	  free_Polyhedron( &(allrbs[k].g) );
	  break;
	  
        case CIRCULARCYLINDER3D:
	  free_CircularCylinder3D( &(allrbs[k].g) );
	  break;

        case CONE:
	  free_Cone( &(allrbs[k].g) );
	  break;
	  
        case TRUNCATEDCONE:
	  free_TruncatedCone( &(allrbs[k].g) );
	  break;
	  
        case ELLIPSOID:
	  free_Ellipsoid( &(allrbs[k].g) );
	  break;	  	  	  	  	
	  
        default:
          fprintf( stderr,"Unknown Rigid Body shape !!\n" );
      }
    }                                      
  }
}



 
/** Prints data of a rigid body */
//----------------------------------------------------------------------------
void print_rigidbody( RigidBody const* p, char const* poshift )
//----------------------------------------------------------------------------
{
  if ( pid() == 0 ) 
  {
    printf( "%sNumber = %lu\n", poshift, p->pnum ); 
    printf( "%sTag = %s\n", poshift, p->typetag ); 
    printf( "%sGeometric type = %lu\n", poshift, p->geomType );     
    printf( "%sShape = ", poshift );
    switch ( p->shape )
    {
      case SPHERE:
        printf( "SPHERE" );
        break;
	  
      case CIRCULARCYLINDER2D:
        printf( "CIRCULARCYLINDER2D" );
        break;
	  
      case CUBE:
        printf( "CUBE" );
        break;
	
      case TETRAHEDRON:
        printf( "TETRAHEDRON" );
	break;
	
      case OCTAHEDRON:
        printf( "OCTAHEDRON" );
	break;
	
      case ICOSAHEDRON:
        printf( "ICOSAHEDRON" );
	break;

      case DODECAHEDRON:
        printf( "DODECAHEDRON" );
	break;
	
      case BOX:
        printf( "BOX" );
	break;
	
      case CIRCULARCYLINDER3D:
        printf( "CIRCULARCYLINDER3D" );
	break;

      case CONE:
        printf( "CONE" );
	break;	
	
      case TRUNCATEDCONE:
        printf( "TRUNCATEDCONE" );
	break;
	
      case ELLIPSOID:
        printf( "ELLIPSOID" );
	break;					
	  
      default:
        fprintf( stderr,"Unknown Rigid Body shape !!\n" );
    }
    printf( "\n" );
    printf( "%sCenter of mass position = %e %e", poshift, p->g.center.x, 
    	p->g.center.y );
#   if dimension == 3
      printf( " %e", p->g.center.z );
#   endif
    printf( "\n" );
    if ( p->g.nperclones )
    {
      printf( "%sNumber of periodic clones = %d\n", poshift, p->g.nperclones );
      printf( "%sPeriodic clones center of mass positions\n", poshift );
      for (int j=0; j < p->g.nperclones; j++)
      {
        printf( "%s   %e %e", poshift, p->g.perclonecenters[j].x, 
		 p->g.perclonecenters[j].y );
#       if dimension == 3
          printf( " %e", p->g.perclonecenters[j].z );
#       endif
        printf( "\n" );          
      }      
    }
    printf( "%sRadius = %e\n", poshift, p->g.radius );  
    printf( "%sMass = %e\n", poshift, p->M ); 
    printf( "%sVolume = %e\n", poshift, p->Vp );     
    printf( "%sDensity = %e\n", poshift, p->rho_s ); 
    if ( p->type != OBSTACLE )
    {
#     if dimension == 3
        printf( "%sInertia tensor\n", poshift );
        printf( "%s   Ixx = %e\n", poshift, p->Ip[0] );
        printf( "%s   Iyy = %e\n", poshift, p->Ip[1] );	  
        printf( "%s   Izz = %e\n", poshift, p->Ip[2] );	  
        printf( "%s   Ixy = %e\n", poshift, p->Ip[3] );	  
        printf( "%s   Ixz = %e\n", poshift, p->Ip[4] );	  
        printf( "%s   Iyz = %e\n", poshift, p->Ip[5] );
#     else
        printf( "%s   Inertia tensor component Izz = %e\n", poshift, 
		p->Ip[2] );
#     endif             
#     if TRANSLATION
        printf( "%sTranslational velocity = %e %e", poshift, p->U.x, p->U.y );
#       if dimension == 3
          printf( " %e", p->U.z );
#       endif
        printf( "\n" );
#     endif
#     if ROTATION
        printf( "%sAngular velocity = ", poshift );
#       if dimension == 3
          printf( "%e %e", p->w.x, p->w.y );
#       endif
        printf( " %e\n", p->w.z );
#     endif
    }
    else
    {
      printf( "%sImposed translational velocity = %e %e", 
    	poshift, p->imposedU.x, p->imposedU.y );
#     if dimension == 3
        printf( " %e", p->imposedU.z );
#     endif
      printf( "\n" );
      printf( "%sImposed angular velocity = ", poshift );
#     if dimension == 3
        printf( "%e %e", p->imposedw.x, p->imposedw.y );
#     endif
      printf( " %e\n", p->imposedw.z );
    }    
  }
  int intpts = p->Interior.n;
  int bdpts = p->Boundary.n;  
# if _MPI
    mpi_all_reduce( intpts, MPI_INT, MPI_SUM );
    mpi_all_reduce( bdpts, MPI_INT, MPI_SUM );
# endif
  if ( pid() == 0 )
  {   
    printf( "%sNumber of interior DLM points = %d\n", poshift, intpts ); 
    printf( "%sNumber of boundary DLM points = %d\n", poshift, p->s.m ); 
    printf( "%sNumber of cells in boundary DLM point stencils = %d\n", poshift, 
    	bdpts );
  }
}  



 
/** Prints all RigidBody* data */
//----------------------------------------------------------------------------
void print_all_rigidbodies( RigidBody const* allrbs, const size_t nrb,
	const size_t npart, const size_t nobs, char const* oshift )
//----------------------------------------------------------------------------
{
  if ( pid() == 0 ) printf( "%sNumber of rigid bodies / particles / obstacles ="
  	" %lu %lu %lu\n", oshift, nrb, npart, nobs );
  char poshift[4]="   ";
  strcat( poshift, oshift );
  for (size_t k=0;k<nrb;k++)
  { 
    if ( pid() == 0 ) printf( "%sRigid body %lu\n", oshift, allrbs[k].pnum );    
    print_rigidbody( &(allrbs[k]), &poshift[0] );
  }
}



 
/** Prints data of a reference rigid body */
//----------------------------------------------------------------------------
void print_referencerigidbody( RigidBody const* p, char const* poshift )
//----------------------------------------------------------------------------
{
  if ( pid() == 0 ) 
  {
    printf( "%sGeometric type = %lu\n", poshift, p->geomType ); 
    printf( "%sTag = %s\n", poshift, p->typetag ); 
    printf( "%sShape = ", poshift );
    switch ( p->shape )
    {
      case SPHERE:
        printf( "SPHERE" );
        break;
	  
      case CIRCULARCYLINDER2D:
        printf( "CIRCULARCYLINDER2D" );
        break;
	  
      case CUBE:
        printf( "CUBE" );
        break;
	
      case TETRAHEDRON:
        printf( "TETRAHEDRON" );
	break;
	
      case OCTAHEDRON:
        printf( "OCTAHEDRON" );
	break;
	
      case ICOSAHEDRON:
        printf( "ICOSAHEDRON" );
	break;

      case DODECAHEDRON:
        printf( "DODECAHEDRON" );
	break;
	
      case BOX:
        printf( "BOX" );
	break;
	
      case CIRCULARCYLINDER3D:
        printf( "CIRCULARCYLINDER3D" );
	break;

      case CONE:
        printf( "CONE" );
	break;	
	
      case TRUNCATEDCONE:
        printf( "TRUNCATEDCONE" );
	break;				

      case ELLIPSOID:
        printf( "ELLIPSOID" );
	break;	
	  
      default:
        fprintf( stderr,"Unknown Rigid Body shape !!\n" );
    }
    printf( "\n" );
    printf( "%sRadius = %e\n", poshift, p->g.radius );  
    printf( "%sMass = %e\n", poshift, p->M ); 
    printf( "%sVolume = %e\n", poshift, p->Vp );     
    printf( "%sDensity = %e\n", poshift, p->rho_s ); 
#   if dimension == 3
        printf( "%sInertia tensor\n", poshift );
        printf( "%s   Ixx = %e\n", poshift, p->Ip[0] );
        printf( "%s   Iyy = %e\n", poshift, p->Ip[1] );	  
        printf( "%s   Izz = %e\n", poshift, p->Ip[2] );	  
        printf( "%s   Ixy = %e\n", poshift, p->Ip[3] );	  
        printf( "%s   Ixz = %e\n", poshift, p->Ip[4] );	  
        printf( "%s   Iyz = %e\n", poshift, p->Ip[5] );
#   else
      printf( "%s   Inertia tensor component Izz = %e\n", poshift, 
		p->Ip[2] );
#   endif
  }             
}  



 
/** Prints reference RigidBody* data */
//----------------------------------------------------------------------------
void print_referencerigidbodies( RigidBody const* allrefrbs, 
	const size_t nrefrb, char const* oshift )
//----------------------------------------------------------------------------
{
  if ( pid() == 0 ) printf( "%sNumber of reference rigid bodies = %lu \n",
  	oshift, nrefrb );
  char poshift[4]="   ";
  strcat( poshift, oshift );
  for (size_t k=0;k<nrefrb;k++)
  { 
    if ( pid() == 0 ) printf( "%sReference rigid body %lu\n", oshift, 
    	allrefrbs[k].geomType );    
    print_referencerigidbody( &(allrefrbs[k]), &poshift[0] );
  }
} 




/** Tags the cells that contain a DLMFD boundary point, i.e. assign the point 
number to the x component of the index field and the rigid body number to the y 
component of the index field. If there is no DLMFD boundary point, index.x is
set to -1 */
//----------------------------------------------------------------------------
void fill_DLM_Index( size_t const* rbnumToIndex, const RigidBodyBoundary dlm_bd,
	vector Index, const size_t pnum, dynUIarray* deactivatedBPindices_,
	dynPDBarray* deactivatedIndexFieldValues_,
	bool* at_least_one_deactivated_ ) 
//----------------------------------------------------------------------------
{    
  for (int i=0;i<dlm_bd.m;i++) 
  {
    # if dimension == 2 
        foreach_point( serial, dlm_bd.bp[i].x, dlm_bd.bp[i].y )
    # elif dimension == 3
        foreach_point( serial, dlm_bd.bp[i].x, dlm_bd.bp[i].y, dlm_bd.bp[i].z )
    # endif    

    if ( level == depth() ) 
    {      
	/* Tag cell only if it was not tagged by another rigid body, else
	store the indices of the BP to be deactivated in both rigid bodies */
	if ( Index.x[] < 0 )
	{
	  Index.x[] = i;
	  Index.y[] = pnum;
	}
	else
	{
	  append_dynUIarray( deactivatedBPindices_, 
	  	rbnumToIndex[(size_t)Index.y[]] );
	  append_dynUIarray( deactivatedBPindices_, (size_t)Index.x[] );
	  append_dynUIarray( deactivatedBPindices_, rbnumToIndex[pnum] );
	  append_dynUIarray( deactivatedBPindices_, i );	  	  	
	  append_dynPDBarray( deactivatedIndexFieldValues_, &(Index.x[]) );
	  append_dynPDBarray( deactivatedIndexFieldValues_, &(Index.y[]) );
	  *at_least_one_deactivated_ = true;
	}
    }
    else  
      printf( "On thread %d, point dlmfd %d of RB %lu at (%f, %f, %f) is in a"
	" cell that has not the maximum level of refinement %d, it "
	"is on level %d \n", pid(), i, pnum, dlm_bd.bp[i].x, dlm_bd.bp[i].y, 
        # if dimension == 3 
            dlm_bd.bp[i].z, 
        # endif	
	depth(), level );    
  }    
}



 
/** Returns whether a point lies inside a rigid body */
//----------------------------------------------------------------------------
bool is_in_rigidbody( RigidBody const* p, double x, double y, double z )
//----------------------------------------------------------------------------
{
  bool is_in = false;
  GeomParameter const* gcp = &(p->g);
  switch ( p->shape )
  {
    case SPHERE:
      is_in = is_in_Sphere( x, y, z, gcp );
      break;
	  
    case CIRCULARCYLINDER2D:
      is_in = is_in_CircularCylinder2D( x, y, gcp );        
      break;
	  
    case CUBE:
      is_in = is_in_Polyhedron( x, y, z, gcp );        
      break;
	
    case TETRAHEDRON:
      is_in = is_in_Polyhedron( x, y, z, gcp );             
      break;
	
    case OCTAHEDRON:
      is_in = is_in_Polyhedron( x, y, z, gcp );             
      break;
	
    case ICOSAHEDRON:
      is_in = is_in_Polyhedron( x, y, z, gcp );             
      break;

    case DODECAHEDRON:
      is_in = is_in_Polyhedron( x, y, z, gcp );             
      break;
      
    case BOX:
      is_in = is_in_Polyhedron( x, y, z, gcp );             
      break;
      
    case CIRCULARCYLINDER3D:
      is_in = is_in_CircularCylinder3D( x, y, z, gcp );             
      break; 

    case CONE:
      is_in = is_in_Cone( x, y, z, gcp );             
      break; 
      
    case TRUNCATEDCONE:
      is_in = is_in_TruncatedCone( x, y, z, gcp );             
      break;
      
    case ELLIPSOID:
      is_in = is_in_Ellipsoid( x, y, z, p );             
      break;                        
	  
    default:
      fprintf( stderr,"Unknown Rigid Body shape !!\n" );
  }

  return ( is_in );
}




/** Tags the grid along rigid body boundaries two cell layers into the fluid */
//----------------------------------------------------------------------------
void fill_FlagMesh( scalar Flag, vector const Index, RigidBody const* allrbs,
	const size_t nrb, size_t const* rbnumToIndex )
//----------------------------------------------------------------------------
{
  foreach() Flag[] = 0; 

  bool goflag = false;  
  foreach_level(MAXLEVEL,serial) 
  {      
    goflag = false;
    if ( (int)Index.x[] > -1 ) goflag = true;
    else
    {
      double xc = x;
      double yc = y;
#     if dimension == 3
        double  zc = z;
#     else
        double zc = 0.;
#     endif      
      foreach_neighbor()
        if ( !goflag ) 
          if ( is_leaf(cell) )
	    if ( (int)Index.x[] > -1 )
	      if ( !is_in_rigidbody( &(allrbs[rbnumToIndex[(size_t)Index.y[]]]),
	      	 xc, yc, zc ) ) 
	        goflag = true;
    }

    if ( goflag ) Flag[] = 1;
  }
  
  synchronize((scalar*) {Flag});    
}  




# include "DLMFD_Weight_functions.h"

/** Returns the weight associated to a cell if it belongs to the 3^dim stencil
associated to a Lagrange multiplier, otherwise return 0. The cell that contains 
the Lagrange multiplier can be on another process (i.e be a ghost cell for the 
current process). This function has to be embedded within a double 
foreach_cache() and foreach_neighbor() loop. The foreach_cache() loops on all 
non-ghost cells on this process that have been tagged to belong to a Lagrange 
multiplier stencil and the foreach_neighbor() loops on the neighbors of the 
non-ghost cells in a 5^dim stencil. */
//----------------------------------------------------------------------------
double reversed_weight( RigidBody* pp, const coord weightcellpos, 
	const coord lambdacellpos, const coord lambdapos, const double delta,
	const int CX, const int NCX ) 
//----------------------------------------------------------------------------
{
  coord relnl;
  int weight_id = 0;
  size_t goflag = 0;
  double weight = 0.;
	
  /* Compute relative vector from the cell to the cell that contains the 
  Lagrange multiplier */
  foreach_dimension() relnl.x = lambdacellpos.x - weightcellpos.x;  

  /* Assign weight id if this cell is flagged ( goflag = 1 ) */
  assign_weight_id_quad_outward( NCX, CX, relnl, delta, &weight_id, &goflag );
    
  /* If the cell belongs to the 3^dim stencil of this Lagrange multiplier, then 
  compute the weight of this cell in the quadratic interpolation
  Otherwise return 0 */
  if ( goflag == 1 ) 
    weight = compute_weight_Quad( weight_id, lambdapos, lambdacellpos, 
    	NCX, CX, delta );	
  
  return ( weight );
}




/** Deactivate boundary points that are too close either to a rigid wall
domain boundary or to another rigid body */
//----------------------------------------------------------------------------
void deactivate_critical_boundary_points( RigidBody* allrbs, const size_t nrb, 
	size_t const* rbnumToIndex, vector Index, 
	dynUIarray* deactivatedBPindices_,
	dynPDBarray* deactivatedIndexFieldValues_,
	bool* at_least_one_deactivated_ ) 
//----------------------------------------------------------------------------
{   
  bool domain_has_rigid_walls = !Period.x || !Period.y || !Period.z;
  double critical_distance = 2. * L0 / (double)(1 << MAXLEVEL) ;
  bool deactivate = false;
  size_t RBid0 = 0, RBid1 = 0, PTid0 = 0, PTid1 = 0, indexp = 0, indexbp;
  double x0, y0, z0;  

  // Distance to rigid walls
  if ( domain_has_rigid_walls )
    foreach_level(MAXLEVEL,serial)
    {      
      if ( (int)Index.x[] > -1 )
      {
        indexp = rbnumToIndex[(size_t)Index.y[]];
	indexbp = (size_t)Index.x[];
	deactivate = false;
	if ( !Period.x )
          if ( allrbs[indexp].s.bp[indexbp].x 
	  		> L0 - critical_distance + X0 
		|| allrbs[indexp].s.bp[indexbp].x  
	  		< critical_distance + X0 )
	    deactivate = true;
	    	    	    
        if ( !Period.y )
          if ( allrbs[indexp].s.bp[indexbp].y 
	  		> L0 - critical_distance + Y0 
		|| allrbs[indexp].s.bp[indexbp].y 
			< critical_distance + Y0 )
	    deactivate = true;

#       if dimension == 3
          if ( !Period.z )
            if ( allrbs[indexp].s.bp[indexbp].z 
	    		> L0 - critical_distance + Z0 
		|| allrbs[indexp].s.bp[indexbp].z 
			< critical_distance + Z0 )
	      deactivate = true;
#       endif 

        if ( deactivate )
	{
    	  append_dynUIarray( deactivatedBPindices_, indexp );
	  append_dynUIarray( deactivatedBPindices_, indexbp );
	  append_dynPDBarray( deactivatedIndexFieldValues_, &(Index.x[]) );
	  append_dynPDBarray( deactivatedIndexFieldValues_, &(Index.y[]) );
	  *at_least_one_deactivated_ = true;
	} 
      }     
    }

  // Proximity of rigid bodies
  if ( nrb > 1 )
  {
    foreach_level(MAXLEVEL,serial)
      if ( (int)Index.x[] > -1 )
      {
	RBid0 = rbnumToIndex[(size_t)Index.y[]];
	PTid0 = (size_t)Index.x[];
        x0 = allrbs[RBid0].s.bp[PTid0].x; 
	y0 = allrbs[RBid0].s.bp[PTid0].y; 
#       if dimension == 3	
	  z0 = allrbs[RBid0].s.bp[PTid0].z;
#       endif 	  
	deactivate = false;
	
        foreach_neighbor()
          if ( is_leaf(cell) )
	    if ( (int)Index.x[] > - 1 )
	    {
	      RBid1 = rbnumToIndex[(size_t)Index.y[]];
	      if ( RBid1 != RBid0 )
	      {
		PTid1 = (size_t)Index.x[];
// 		if ( sqrt( 
// 			sq( allrbs[RBid1].s.x[PTid1] - x0 ) 
// 			+ sq( allrbs[RBid1].s.y[PTid1] - y0 )
// #                 if dimension == 3	
// 	            	+ sq( allrbs[RBid1].s.z[PTid1] - z0 )
// #                 endif
//                   ) < critical_distance )
// 		  deactivate = true;
		  
		if ( abs( allrbs[RBid1].s.bp[PTid1].x - x0 ) < critical_distance
			|| abs( allrbs[RBid1].s.bp[PTid1].y - y0 ) 
				< critical_distance
#                 if dimension == 3	
	            || abs( allrbs[RBid1].s.bp[PTid1].z - z0 ) 
		    	< critical_distance
#                 endif
                  ) deactivate = true;		  		 
              }
	    }
	          
        if ( deactivate ) 
	{
          append_dynUIarray( deactivatedBPindices_, RBid0 );
	  append_dynUIarray( deactivatedBPindices_, PTid0 );
	  append_dynPDBarray( deactivatedIndexFieldValues_, &(Index.x[]) );
	  append_dynPDBarray( deactivatedIndexFieldValues_, &(Index.y[]) );
          *at_least_one_deactivated_ = true;
	}
      }
  }

#   if  _MPI
      int local = *at_least_one_deactivated_, global;
      MPI_Allreduce( &local, &global, 
      	1, MPI_INT, MPI_SUM, MPI_COMM_WORLD ); 
      if ( global > 0 ) *at_least_one_deactivated_ = true;
#   endif
   
  if ( *at_least_one_deactivated_ )
  {
    size_t total_size = 0;
    size_t* alldeactivatedBPindices = NULL;

    // In MPI concatenate the arrays of deactivated boundary points
#   if  _MPI
      int* deactivated_count = NULL;
      int* deactivated_displace = NULL;
      int nb_ranks;
      MPI_Comm_size( MPI_COMM_WORLD, &nb_ranks );
      
      deactivated_count = (int*) calloc( nb_ranks, sizeof(int) );
      deactivated_displace = (int*) calloc( nb_ranks, sizeof(int) );
      
      int nn = (int)deactivatedBPindices_->n;
      MPI_Allgather( &nn, 1, MPI_INT, deactivated_count, 1, MPI_INT, 
      	MPI_COMM_WORLD );

      int temp = 0;
      for (size_t i = 0; i < nb_ranks; i++)
      {
         deactivated_displace[i] = temp;
         temp += deactivated_count[i];
      }
	  		
      MPI_Allreduce( &(deactivatedBPindices_->n), &total_size, 
      	1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD ); 
   
      alldeactivatedBPindices = (size_t*) calloc( total_size, 
    	sizeof(size_t) );
	
      MPI_Allgatherv( deactivatedBPindices_->elem, deactivatedBPindices_->n, 
    	MPI_UNSIGNED_LONG, alldeactivatedBPindices, deactivated_count, 
	deactivated_displace, MPI_UNSIGNED_LONG, MPI_COMM_WORLD );	
#   else
      total_size = deactivatedBPindices_->n;
      alldeactivatedBPindices = deactivatedBPindices_->elem;
#   endif

    // Deactivate boundary points in the rigid bodies
    for (size_t i=0;i<total_size;i+=2)
      allrbs[alldeactivatedBPindices[i]].s.deactivated[
      	alldeactivatedBPindices[i+1]] = true;
 
#   if  _MPI
      free( deactivated_count );
      free( deactivated_displace );
      free( alldeactivatedBPindices );
#   endif 

    // Set the stored Index values to - 1
    for (size_t i=0;i<deactivatedIndexFieldValues_->n;++i)
      *(deactivatedIndexFieldValues_->elem[i]) = - 1;     

    synchronize((scalar*) {Index});
  }
}



 
/** Tag cells that belong to a 3^dim stencil associated to a Lagrange multiplier
point of the rigid body boundary. */
//----------------------------------------------------------------------------
void reverse_fill_DLM_Flag( RigidBody* allrbs, const size_t nrb, 
	size_t const* rbnumToIndex, scalar Flag, vector CX_NCX, 
	const vector Index ) 
//----------------------------------------------------------------------------
{
  coord rel = {0., 0., 0.};
  coord relnl = {0., 0., 0.};
  int NCX, CX, weight_id, pindex, nrbID = 0;
  int rbID[20], bpnum; 
  size_t goflag = 0;
  coord lambdacellpos = {0., 0., 0.};
  coord lambdapos = {0., 0., 0.};
  coord localcellpos = {0., 0., 0.}; 
  coord normal = {0., 0., 0.}; 
  GeomParameter const* gcp = NULL;


  // Assign values to CX_NCX 
  foreach_level(MAXLEVEL,serial) 
    if ( (int)Index.x[] > -1 )
    {
      bpnum = (int)Index.x[];
      pindex = rbnumToIndex[(int)Index.y[]];
      gcp = &(allrbs[pindex].g);       
      lambdacellpos.x = x;      
      lambdacellpos.y = y;
      lambdapos.x = allrbs[pindex].s.bp[bpnum].x;
      lambdapos.y = allrbs[pindex].s.bp[bpnum].y;
      normal.x = allrbs[pindex].s.normal[bpnum].x;
      normal.y = allrbs[pindex].s.normal[bpnum].y;            
#     if dimension == 3 
        lambdacellpos.z = z;
	lambdapos.z = allrbs[pindex].s.bp[bpnum].z;
	normal.z = allrbs[pindex].s.normal[bpnum].z; 	
#     endif

      /* Compute relative vector from the cell (containning the boundary) 
      position to the boundary's (analytical) position */
      foreach_dimension() rel.x = lambdapos.x - lambdacellpos.x;

      /* In a periodic case, the position of the Lagrange multipliers on 
      the boundary are those of the master rigid body, i.e., they are not 
      shifted for each periodic clone. Therefore we may have a case where 
      the cell and the position of the Lagrange multiplier are on opposite 
      side of the domain. The solution is to subtract/add the periodic 
      domain length L0 when this happens */
      foreach_dimension()
      {
        if ( rel.x > L0/2. ) rel.x = rel.x - L0;
        else if ( rel.x < -L0/2. ) rel.x = rel.x + L0;
      }
      
      /* Reset dial integers */
      NCX = 0; CX = 0; 
       
      /* Assign quadrant number CX defining relative position of the 
      Lagrange point with respect to the center of the cell it belongs to */
      assign_dial( rel, &CX );
      CX_NCX.x[] = CX;
       
      /* Assign quadrant number NCX given by the direction of the normal 
      over the geometric boundary of the rigid body */ 
      assign_dial_fd_boundary( &allrbs[pindex], lambdapos, gcp, Delta, normal,
      	&NCX );
      CX_NCX.y[] = NCX;                      
    }

  synchronize((scalar*){CX_NCX.x, CX_NCX.y}); 
	

  // Tag cells that belong to a 3^dim stencil
  foreach_level(MAXLEVEL,serial) 
  {      
    localcellpos.x = x;
    localcellpos.y = y;
#   if dimension == 3 
      localcellpos.z = z;
#   endif
    nrbID = 0;

    // Check if there is a Lagrange multiplier in the neigborhood of this 
    // cell
    foreach_neighbor() 
    {
      if ( (int)Index.x[] > -1 && level == depth() && is_leaf(cell) ) 
      {		
        goflag = 0;
	pindex = rbnumToIndex[(int)Index.y[]]; 

        /* Compute relative vector from the cell to the cell that contains 
	the Lagrange multiplier */
        relnl.x = x - localcellpos.x;
	relnl.y = y - localcellpos.y;
#       if dimension == 3 
          relnl.z = z - localcellpos.z;
#       endif
	
	/* Assign weight id if this cell is flagged ( goflag = 1 ) */
	assign_weight_id_quad_outward( (int)CX_NCX.y[], (int)CX_NCX.x[], relnl,
		Delta, &weight_id, &goflag );		
		
	if ( goflag == 1 )
	{
	  bool already_inserted = false;
	  for (int m=0;m<nrbID;++m)
	    if ( rbID[m] == pindex ) already_inserted = true;
	  if ( !already_inserted ) { rbID[nrbID] = pindex; ++nrbID; }  
	}		
      } // end if (Index.x[] > -1)
    } // end foreach_neigboor loop
    
    /* If the cell belongs to at least one stencil of a Lagrange multiplier 
    of a rigid body tag it and add it the reduced domain of all rigid body
    stencils to optimize the stencil-traversal cost */
    if ( nrbID ) 
    {
      Flag[] = 1;
      Point ppp;
      ppp.i = point.i;
      ppp.j = point.j;
      ppp.k = point.k;			
      ppp.level = point.level;	  
      for (int m=0;m<nrbID;++m)	
	cache_append( &(allrbs[rbID[m]].Boundary), ppp, 0 );
    }    
  }

  synchronize({Flag});       
}




/** Creates DLM/FD boundary points of a given rigid body. Sets the
PeriodicRefCenter field only if setPeriodicRefCenter is true */
//----------------------------------------------------------------------------
void create_boundary_points( RigidBody* p, RigidBody const* refrb,
	vector* pPeriodicRefCenter, 
	const bool setPeriodicRefCenter )
//----------------------------------------------------------------------------
{   
  coord pos, rotp;
  coord shift = {0., 0., 0.};
  coord ori = {X0, Y0, Z0};
  size_t m = refrb->s.m; 
  
  allocate_RigidBodyBoundary( &(p->s), m );
  
  for (size_t i=0;i<m;++i)
  {
    // Rotate point of the reference rigid body
    matCoordDotProduct( p->RotMat, refrb->s.bp[i], &rotp );
    
    // Translate the rotated point   
    foreach_dimension()
      pos.x = p->g.center.x + rotp.x;
      
    /* Check if the point falls outside of the domain in case of periodicity */
    foreach_dimension()
    {
      shift.x = 0.;      
      if ( Period.x )
      {
	if ( pos.x > L0 + ori.x )
	{  
	  pos.x -= L0;
	  shift.x = - L0;
	}
        else if ( pos.x < 0. + ori.x )
	{
	  pos.x += L0;
	  shift.x = L0;
	}
      }
    }
    
    if ( setPeriodicRefCenter )
    {
      // Setting the periodic clone center vector field
      foreach_point( pos.x, pos.y, pos.z )
        foreach_dimension()
	  pPeriodicRefCenter->x[] = p->g.center.x + shift.x;
    }

    foreach_dimension() 
      p->s.bp[i].x = pos.x; 
      
    // Rotate normal vector of the reference rigid body
    matCoordDotProduct( p->RotMat, refrb->s.normal[i], &(p->s.normal[i]) );
  }      
}





/** Creates DLM/FD boundary points of a given rigid body. Sets the
PeriodicRefCenter field only if setPeriodicRefCenter is true */
//----------------------------------------------------------------------------
void create_referencerigidbody_boundary_geomfeatures( RigidBody* p )
//----------------------------------------------------------------------------
{  
  GeomParameter gci = p->g;
  int m = 0;
  int lN = 0;
    
  switch( p->shape )
  {
    case SPHERE:
      compute_nboundary_Sphere( &gci, &m );
      allocate_RigidBodyBoundary( &(p->s), m );
      create_referenceRB_boundary_geomfeatures_Sphere( &gci, &(p->s), m );
      break;
	  
    case CIRCULARCYLINDER2D:
      compute_nboundary_CircularCylinder2D( &gci, &m );
      allocate_RigidBodyBoundary( &(p->s), m );
      create_referenceRB_boundary_geomfeatures_CircularCylinder2D( &gci, 
      	&(p->s), m );
      break;
	  
    case CUBE:
      compute_nboundary_Cube( &gci, &m, &lN );
      allocate_RigidBodyBoundary( &(p->s), m );
      create_referenceRB_boundary_geomfeatures_Cube( &gci, &(p->s), m, lN );
      break;
	
    case TETRAHEDRON:
      compute_nboundary_Tetrahedron( &gci, &m, &lN );
      allocate_RigidBodyBoundary( &(p->s), m );
      create_referenceRB_boundary_geomfeatures_Tetrahedron( &gci, &(p->s), m, 
      	lN );
      break;
	
    case OCTAHEDRON:
      compute_nboundary_Octahedron( &gci, &m, &lN );
      allocate_RigidBodyBoundary( &(p->s), m );
      create_referenceRB_boundary_geomfeatures_Octahedron( &gci, &(p->s), m, 
      	lN );      
      break;

    case ICOSAHEDRON:
      compute_nboundary_Icosahedron( &gci, &m, &lN );
      allocate_RigidBodyBoundary( &(p->s), m );
      create_referenceRB_boundary_geomfeatures_Icosahedron( &gci, &(p->s), m, 
      	lN );
      break;	

    case DODECAHEDRON:     
      compute_nboundary_Dodecahedron( &gci, &m, &lN );	
      allocate_RigidBodyBoundary( &(p->s), m );
      create_referenceRB_boundary_geomfeatures_Dodecahedron( &gci, &(p->s), m, 
      	lN );      
      break;
      
    case BOX:     
      compute_nboundary_Box( &gci, &m );	
      allocate_RigidBodyBoundary( &(p->s), m );
      create_referenceRB_boundary_geomfeatures_Box( &gci, &(p->s), m );
      break;
      
    case CIRCULARCYLINDER3D:     
      compute_nboundary_CircularCylinder3D( &gci, &m );	
      allocate_RigidBodyBoundary( &(p->s), m );
      create_referenceRB_boundary_geomfeatures_CircularCylinder3D( &gci, 
      	&(p->s), m );
      break; 

    case CONE:     
      compute_nboundary_Cone( &gci, &m );	
      allocate_RigidBodyBoundary( &(p->s), m );
      create_referenceRB_boundary_geomfeatures_Cone( &gci, &(p->s), m );
      break;
      
    case TRUNCATEDCONE:     
      compute_nboundary_TruncatedCone( &gci, &m );	
      allocate_RigidBodyBoundary( &(p->s), m );
      create_referenceRB_boundary_geomfeatures_TruncatedCone( &gci, 
      	&(p->s), m );
      break;
      
    case ELLIPSOID:     
      compute_nboundary_Ellipsoid( &gci, &m );	
      allocate_RigidBodyBoundary( &(p->s), m );
      create_referenceRB_boundary_geomfeatures_Ellipsoid( &gci, 
      	&(p->s), m );
      break;                       	
	  
    default:
      fprintf( stderr, "Unknown Rigid Body shape !!\n" );
  }
}


 
/** Initializes the rigid bodies and the scalar/vector fields needed to the 
method */
//----------------------------------------------------------------------------
void allocate_and_init_rigidbodies( RigidBody* allrbs, const size_t nrb,
	size_t const* rbnumToIndex, RigidBody const* allrefrbs,
	vector Index, vector CX_NCX, scalar Flag, scalar FlagMesh, 
	vector PeriodicRefCenter, dynUIarray* deactivatedBPindices_,
	dynPDBarray* deactivatedIndexFieldValues_,
	bool* at_least_one_deactivated_ )
//----------------------------------------------------------------------------
{  
  /* Initialize fields */
  foreach() 
  {
    Index.x[] = - 1;   // DLM/FD boundary point index
    Index.y[] = - 1;   // rigid body number
    CX_NCX.x[] = - 1;  // Relative position the Lagrange point with respect to 
    		       // the center of the cell it belongs to
    CX_NCX.y[] = - 1;  // Quadrant number given by the direction of the normal 
		       // over the geometric boundary of the rigid body
    Flag[] = 0;        // Flag
    FlagMesh[] = 0;    // FlagMesh
    foreach_dimension() PeriodicRefCenter.x[] = 0.;
  }

  for (size_t k = 0; k < nrb; k++) 
  {
    Cache* c = NULL;

#   if DLMFD_BOUNDARYPOINTS
      create_boundary_points( &(allrbs[k]), &(allrefrbs[allrbs[k].geomType]), 
      		&PeriodicRefCenter, true );
      fill_DLM_Index( rbnumToIndex, (allrbs[k].s), Index, allrbs[k].pnum, 
      	deactivatedBPindices_, deactivatedIndexFieldValues_, 
	at_least_one_deactivated_ );
      c = &(allrbs[k].Boundary);
      initialize_and_allocate_Cache( c );
#   endif     
#   if DLMFD_INTERIORPOINTS
      c = &(allrbs[k].Interior);
      initialize_and_allocate_Cache( c );
#   endif
  }

  synchronize((scalar*){PeriodicRefCenter.x,
  	PeriodicRefCenter.y, PeriodicRefCenter.z});  
  synchronize((scalar*){Index});
}




/** Creates boundary points and normal vectors of reference rigid bodies */
//----------------------------------------------------------------------------
void create_referencerigidbodies_boundary_geomfeatures( RigidBody* allrefrbs, 
	const size_t nrefrb )
//----------------------------------------------------------------------------
{  
# if DLMFD_BOUNDARYPOINTS
    for (size_t k = 0; k < nrefrb; k++) 
      create_referencerigidbody_boundary_geomfeatures( &(allrefrbs[k]) );    
# endif
}




/** Frees boundary points of all rigid bodies */
//----------------------------------------------------------------------------
void free_rigidbodies_boundary_points( RigidBody* allrbs, const size_t nrb )
//----------------------------------------------------------------------------
{  
# if DLMFD_BOUNDARYPOINTS
    for (size_t k = 0; k < nrb; k++) 
      free_RigidBodyBoundary( &(allrbs[k].s) );   
# endif
}




/** Writes headers in rigid body data output files */
//----------------------------------------------------------------------------
void writer_headers( FILE* pdata, const bool pdata_is_open, FILE* sl ) 
//----------------------------------------------------------------------------
{
  // Comments: 
  // 1) we assume that initial time is always 0 when the simulation is 
  // not restarted
  // 2) we assume all hydro forces and torques are 0 at t=0
  // 3) we position the headers over the 1st line as a function of NSD, the 
  // number of significant digits after the decimal point

# if _MPI
    if ( pid() == 0 )
# endif
    {
      /* Write header for rigid body positions and velocities */
#     if dimension == 2
        // Position and velocity file
	if ( pdata_is_open )
	{
          if ( NSDF > 9 )
            fprintf( pdata, "#time\t\t\tposition.x\t\tposition.y"
    		"\t\tU.x\t\t\tU.y\t\t\tw.z\n" );
          else
            fprintf( pdata, "#time\t\tposition.x\tposition.y"
    		"\tU.x\t\tU.y\t\tw.z\n" );
	}

        // Hydro force and torque file
        if ( NSDF > 9 )
          fprintf( sl, "#time\t\t\tF.x\t\t\tF.y\t\t\tT.z\n" );
        else
          fprintf( sl, "#time\t\tF.x\t\tF.y\t\tT.z\n" );
        for (int i=0; i<3; ++i) fprintf( sl, "%.*e\t", NSDF, 0. ); 
        fprintf( sl, "%.*e\n", NSDF, 0. );   
#     elif dimension == 3
        // Position and velocity file
        if ( pdata_is_open )
	{
	  if ( NSDF > 9 )
            fprintf( pdata, "#time\t\t\tposition.x\t\tposition.y\t\tposition.z"
    		"\t\tU.x\t\t\tU.y\t\t\tU.z"
		"\t\t\tw.x\t\t\tw.y\t\t\tw.z\n" );
          else
            fprintf( pdata, "#time\t\tposition.x\tposition.y\tposition.z"
    		"\tU.x\t\tU.y\t\tU.z"
		"\t\tw.x\t\tw.y\t\tw.z\n" );
	}

        // Hydro force and torque file
        if ( NSDF > 9 )
          fprintf( sl, "#time\t\t\tF.x\t\t\tF.y\t\t\tF.z\t\t\tT.x\t\t\tT.y"
      		"\t\t\tT.z\n" );
        else
          fprintf( sl, "#time\t\tF.x\t\tF.y\t\tF.z\t\tT.x\t\tT.y\t\tT.z\n" );
        for (int i=0; i<6; ++i) fprintf( sl, "%.*e\t", NSDF, 0. ); 
        fprintf( sl, "%.*e\n", NSDF, 0. );
#     endif
    
      // Flush the buffers such that files are updated dynamically
      fflush( sl );
      fflush( pdata );
    }
}




/** Writes rigid body data in files */
//----------------------------------------------------------------------------
void rigidbody_data( RigidBody* allrbs, const size_t nrb, const double t, 
	const int i, FILE** pdata ) 
//----------------------------------------------------------------------------
{  
  for (size_t k = 0; k < nrb; k++) 
  {
    GeomParameter* GCi = &(allrbs[k].g);
    coord* U;
    coord* w;    
    if ( allrbs[k].type != OBSTACLE )
    {
#     if TRANSLATION 
        U = &(allrbs[k].U);
#     endif
#     if ROTATION
        w = &(allrbs[k].w);
#     endif
    }
    else
    {
#     if TRANSLATION 
        U = &(allrbs[k].imposedU);
#     endif
#     if ROTATION
        w = &(allrbs[k].imposedw);
#     endif
    }
    
#   if dimension == 2
      if ( pid() == 0 ) 
      {
        fprintf( pdata[k], "%.*e\t", NSDF, t );
        fprintf( pdata[k], "%.*e\t", NSDF, (*GCi).center.x );
        fprintf( pdata[k], "%.*e\t", NSDF, (*GCi).center.y );         
        fprintf( pdata[k], "%.*e\t%.*e\t%.*e\n"     
#       if TRANSLATION
	  , NSDF, (*U).x, NSDF, (*U).y 
#       elif TRANSLATION == 0
	  , NSDF, 0., NSDF, 0.
#       endif   
#       if ROTATION
	  , NSDF, (*w).z
#       elif ROTATION == 0
	  , NSDF, 0.
#       endif
        );      
        fflush(pdata[k]);
      }
#   elif dimension == 3
      if ( pid() == 0 ) 
      {
        fprintf( pdata[k], "%.*e\t", NSDF, t );
        fprintf( pdata[k], "%.*e\t", NSDF, (*GCi).center.x );
        fprintf( pdata[k], "%.*e\t", NSDF, (*GCi).center.y );      
        fprintf( pdata[k], "%.*e\t", NSDF, (*GCi).center.z );      
        fprintf( pdata[k], "%.*e\t%.*e\t%.*e\t%.*e\t%.*e\t%.*e\n"     
#       if TRANSLATION
	  , NSDF, (*U).x, NSDF, (*U).y, NSDF, (*U).z
#       elif TRANSLATION == 0
	  , NSDF, 0., NSDF, 0., NSDF, 0.
#       endif   
#       if ROTATION
	  , NSDF, (*w).x, NSDF, (*w).y, NSDF, (*w).z
#       elif ROTATION == 0
	  , NSDF, 0., NSDF, 0., NSDF, 0.
#       endif
        );      
        fflush(pdata[k]);      	
      }
#   endif
  }
}




/** Compute hydrodynamic force & torque and write the values in files */
//----------------------------------------------------------------------------
void computeHydroForceTorque( RigidBody* allrbs, const size_t nrb, FILE** sl, 
	const double t, const double dt, 
	scalar Flag, vector lambda, vector Index, 
	const double rho_f, vector prefcenter ) 
//----------------------------------------------------------------------------
{
  /* Compute hydrodynamic force & torque and write to a file */
  /* Fh = <lambda,V>_P + (rho_f/rho_s)MdU/dt */
  /* Mh = <lambda,xi^GM>_P + (rho_f/rho_s).( Idom/dt + om x (I.om)) */

  /* The inertia tensor is: */
  /* Ip[0] = Ixx */
  /* Ip[1] = Iyy */
  /* Ip[2] = Izz */
  /* Ip[3] = Ixy */
  /* Ip[4] = Ixz */
  /* Ip[5] = Iyz */  

  coord lambdasumint;
  coord lambdasumboundary;
  coord lambdasum;
  coord crossLambdaSumInt;
  coord crossLambdaSumBoundary;
  coord crossLambdaSum; 
  Cache* Interior[nbRigidBodies];
  Cache* Boundary[nbRigidBodies];
  RigidBodyBoundary * sbm;
  
  /* Loop over all rigid bodies */
  for (size_t k = 0; k < nrb; k++) 
  {
    foreach_dimension() 
    {
      lambdasumint.x = 0.;
      lambdasumboundary.x = 0.;
      lambdasum.x = 0;
      
      crossLambdaSumInt.x = 0.;
      crossLambdaSumBoundary.x = 0.;
      crossLambdaSum.x = 0.;
    }

#   if dimension == 2 
      lambdasumint.z = 0.;
      lambdasumboundary.z = 0.;
      lambdasum.z = 0;
      
      crossLambdaSumInt.z = 0.;
      crossLambdaSumBoundary.z = 0.;
      crossLambdaSum.z = 0.;
#   endif

    /* Interior domain of the RigidBody* */
    Interior[k] = &(allrbs[k].Interior);

    /* Boundary domain of the RigidBody* */
    Boundary[k] = &(allrbs[k].Boundary);
    sbm = &(allrbs[k].s);
    coord lambdapos;
    double rho_s = allrbs[k].rho_s;
    if ( allrbs[k].type == OBSTACLE ) rho_s = rho_f;
#   if TRANSLATION
      double M = allrbs[k].M;
      coord Unm1 = allrbs[k].Unm1;
      coord U = allrbs[k].U;
#   endif
#   if ROTATION
      coord wnm1 = allrbs[k].wnm1;
      coord w = allrbs[k].w;
      coord Iom, omIom, Idomdt;
#   endif
    
    // Particle's interior multipliers
    foreach_cache((*Interior[k])) 
    {
      if ( Flag[] < 1 && k == Index.y[] ) 
      {
	/* Compute Fh = <lambda,V>_P */
	/* For moving rigid body the additional term +
	   (rho_f/rho_s).M.dU/dt is added after the mpi_calls below */
	
	foreach_dimension() 
	  lambdasumint.x += lambda.x[];

	/* Compute Mh = <lambda,xi^GM>_P */  
	/* For moving rigid body, the additional term +
	   (rho_f/rho_s).(I.dom/dt + om ^ (I.om)) is added after mpi
	   calls below */
#       if dimension == 3
	  crossLambdaSumInt.x += ( lambda.z[] * ( y - prefcenter.y[] )
		- lambda.y[] * ( z - prefcenter.z[] ) );
	  crossLambdaSumInt.y += ( lambda.x[] * ( z - prefcenter.z[] )
		- lambda.z[] * ( x - prefcenter.x[] ) );
#       endif
	crossLambdaSumInt.z += ( lambda.y[] * ( x - prefcenter.x[] )
		- lambda.x[] * ( y - prefcenter.y[] ) );
      }
    }

    // Particle's boundary multipliers
    foreach_cache((*Boundary[k])) 
    {
      if ( Index.x[] > -1 && allrbs[k].pnum == (int)Index.y[] ) 
      {
	lambdapos.x = (*sbm).bp[(int)Index.x[]].x;
	lambdapos.y = (*sbm).bp[(int)Index.x[]].y;
	lambdapos.z = 0.;
#       if dimension == 3
	  lambdapos.z = (*sbm).bp[(int)Index.x[]].z;
#       endif

	/* Compute Fh = <lambda,V>_P */
	foreach_dimension()
	  lambdasumboundary.x += lambda.x[];

	/* Modify temporerly the RigidBody* center position for periodic 
	boundary condition */
#       if dimension == 3
          crossLambdaSumBoundary.x += 
		( lambda.z[] * ( lambdapos.y - prefcenter.y[] )
		- lambda.y[] * ( lambdapos.z - prefcenter.z[] ) );
          crossLambdaSumBoundary.y += 
		( lambda.x[] * ( lambdapos.z - prefcenter.z[] )
		- lambda.z[] * ( lambdapos.x - prefcenter.x[] ) );
#       endif		
        crossLambdaSumBoundary.z += 
		( lambda.y[] * ( lambdapos.x - prefcenter.x[] )
		- lambda.x[] * ( lambdapos.y - prefcenter.y[] ) );
      }
    }
 
#   if _MPI
#     if dimension == 3
        foreach_dimension() 
        {
          mpi_all_reduce( lambdasumint.x, MPI_DOUBLE, MPI_SUM );
          mpi_all_reduce( lambdasumboundary.x, MPI_DOUBLE, MPI_SUM );
          mpi_all_reduce( crossLambdaSumInt.x, MPI_DOUBLE, MPI_SUM );
          mpi_all_reduce( crossLambdaSumBoundary.x, MPI_DOUBLE, MPI_SUM );
        }
#     else
        mpi_all_reduce( lambdasumint.x, MPI_DOUBLE, MPI_SUM );
        mpi_all_reduce( lambdasumboundary.x, MPI_DOUBLE, MPI_SUM );
        mpi_all_reduce( lambdasumint.y, MPI_DOUBLE, MPI_SUM );
        mpi_all_reduce( lambdasumboundary.y, MPI_DOUBLE, MPI_SUM );
        mpi_all_reduce( crossLambdaSumInt.z, MPI_DOUBLE, MPI_SUM );
        mpi_all_reduce( crossLambdaSumBoundary.z, MPI_DOUBLE, MPI_SUM );
#     endif
#   endif


  /* The force term (rho_f/rho_s)MdU/dt and torque term 
  (rho_f/rho_s)(Idw/dt + w x Iw) are added for all rigid bodies */

#   if TRANSLATION
      foreach_dimension()
        lambdasumint.x += ( rho_f / rho_s ) * M * ( U.x - Unm1.x ) / dt;
#   endif
    
#   if ROTATION
      /* I.om term */
#     if dimension == 3
        Iom.x = allrbs[k].Ip[0] * w.x 
	  	+ allrbs[k].Ip[3] * w.y + allrbs[k].Ip[4] * w.z;
        Iom.y = allrbs[k].Ip[3] * w.x 
	  	+ allrbs[k].Ip[1] * w.y + allrbs[k].Ip[5] * w.z;
        Iom.z = allrbs[k].Ip[4] * w.x 
	  	+ allrbs[k].Ip[5] * w.y + allrbs[k].Ip[2] * w.z;
#     else
        Iom.z = allrbs[k].Ip[2] * w.z;
#     endif	 

      /* om^(I.om) term */
#     if dimension == 3
        omIom.x = w.y * Iom.z - w.z * Iom.y;
        omIom.y = w.z * Iom.x - w.x * Iom.z;
        omIom.z = w.x * Iom.y - w.y * Iom.x;
#     else
        omIom.z = 0.;
#     endif	  

      /* Idom/dt term */
#     if dimension == 3	
        Idomdt.x = allrbs[k].Ip[0] * ( w.x - wnm1.x ) 
	  	+ allrbs[k].Ip[3] * ( w.y - wnm1.y ) 
    		+ allrbs[k].Ip[4] * ( w.z - wnm1.z );
        Idomdt.y = allrbs[k].Ip[3] * ( w.x - wnm1.x ) 
	  	+ allrbs[k].Ip[1] * ( w.y - wnm1.y ) 
    		+ allrbs[k].Ip[5] * ( w.z - wnm1.z );
        Idomdt.z = allrbs[k].Ip[4] * ( w.x - wnm1.x ) 
	  	+ allrbs[k].Ip[5] * ( w.y - wnm1.y ) 
    		+ allrbs[k].Ip[2] * ( w.z - wnm1.z );
#     else
        Idomdt.z = allrbs[k].Ip[2] * ( w.z - wnm1.z );
#     endif		
	
#     if dimension == 3
        crossLambdaSumInt.x += ( rho_f / rho_s ) * ( ( Idomdt.x ) / dt 
		+ omIom.x );
        crossLambdaSumInt.y += ( rho_f / rho_s ) * ( ( Idomdt.y ) / dt 
		+ omIom.y );
#     endif		
      crossLambdaSumInt.z += ( rho_f / rho_s ) * ( ( Idomdt.z ) / dt 
		+ omIom.z );
#   endif	  

#   if dimension == 3    
      foreach_dimension() 
      {
        lambdasum.x = lambdasumint.x + lambdasumboundary.x;
        crossLambdaSum.x = crossLambdaSumInt.x + crossLambdaSumBoundary.x;
      }
#   else
      lambdasum.x = lambdasumint.x + lambdasumboundary.x;
      lambdasum.y = lambdasumint.y + lambdasumboundary.y;            
      crossLambdaSum.z = crossLambdaSumInt.z + crossLambdaSumBoundary.z;
#   endif
      
    if( pid() == 0 ) 
#     if dimension == 3  
      {
        fprintf( sl[k], "%.*e\t", NSDF, t );
        fprintf( sl[k], "%.*e\t", NSDF, lambdasum.x );      
        fprintf( sl[k], "%.*e\t", NSDF, lambdasum.y );      
        fprintf( sl[k], "%.*e\t", NSDF, lambdasum.z );      
        fprintf( sl[k], "%.*e\t", NSDF, crossLambdaSum.x );      
        fprintf( sl[k], "%.*e\t", NSDF, crossLambdaSum.y );      
        fprintf( sl[k], "%.*e\n", NSDF, crossLambdaSum.z );      
        fflush( sl[k] );
      }
#     else
      {
        fprintf( sl[k], "%.*e\t", NSDF, t );
        fprintf( sl[k], "%.*e\t", NSDF, lambdasum.x );      
        fprintf( sl[k], "%.*e\t", NSDF, lambdasum.y );          
        fprintf( sl[k], "%.*e\n", NSDF, crossLambdaSum.z );      
        fflush( sl[k] );
      }
#endif
  }
}




/** Initialize/open all DLMFD file pointers */
//----------------------------------------------------------------------------
void init_file_pointers( RigidBody const* allrbs, const size_t nrb, FILE** p, 
	const bool open_pdata, FILE** d, FILE** UzawaCV, FILE** CVT, 
	const size_t rflag ) 
//----------------------------------------------------------------------------
{
  char name[80] = "";
  char suffix[80] = "";
  char buffer[80] = "";  
# if _MPI
    if ( pid() == 0 )
# endif
    {
      // Particle data
      for (size_t k = 0; k < nrb; k++) 
      {
        sprintf( suffix, "_%lu.dat", allrbs[k].pnum );

        if ( open_pdata )
	{
          strcpy( name, RESULT_DIR );
          strcat( name, "/" );
          strcat( name, RESULT_RIGIDBODY_VP_ROOTFILENAME );
          strcat( name, suffix );
	  if ( !rflag ) p[k] = fopen( name,  "w" ); 
	  else p[k] = fopen( name,  "a" );
	}
		 	       
        strcpy( name, RESULT_DIR );
        strcat( name, "/" );
        strcat( name, RESULT_RIGIDBODY_HYDROFAT_ROOTFILENAME );
        strcat( name, suffix );      
        if ( !rflag ) d[k] = fopen( name, "w" );
	else d[k] = fopen( name, "a" );
	
	// Write headers in these files
	if ( !rflag ) 
	  writer_headers( open_pdata ? p[k] : NULL, open_pdata, d[k] );
      }
    
      // Uzawa convergence
      char CONVERGE_UZAWA_FILENAME_complete_name[80];
      strcpy( buffer, RESULT_DIR );
      strcat( buffer, "/" );
      strcat( buffer, CONVERGE_UZAWA_FILENAME );
      strcpy( CONVERGE_UZAWA_FILENAME_complete_name, buffer );
      if ( !rflag )
      {
        *UzawaCV = fopen( CONVERGE_UZAWA_FILENAME_complete_name, "w" ); 
        fprintf( *UzawaCV, "# Iter \t Uzawa Iter \t ||u-u_imposed||\n" );
      }
      else
        *UzawaCV = fopen( CONVERGE_UZAWA_FILENAME_complete_name, "a" );   
    
      // Cells, contrained cells and nb of multipliers 
      char DLMFD_CELLS_FILENAME_complete_name[80]; 
      strcpy( buffer, RESULT_DIR );
      strcat( buffer, "/" );
      strcat( buffer, DLMFD_CELLS_FILENAME );
      strcpy( DLMFD_CELLS_FILENAME_complete_name, buffer );
      if ( !rflag )
      {
        *CVT = fopen ( DLMFD_CELLS_FILENAME_complete_name, "w" ); 
        fprintf ( *CVT,"# Iter \t DLMFDPts \t DLMFDCells \t "
    		"TotalCells\n" );
      }
      else
        *CVT = fopen( DLMFD_CELLS_FILENAME_complete_name, "a" );          
    }    
}




/** Close all DLMFD files */
//----------------------------------------------------------------------------
void close_file_pointers( const size_t nrb, FILE** p, const bool pdata_is_open,
	FILE** d, FILE* UzawaCV, FILE* CVT ) 
//----------------------------------------------------------------------------
{ 
# if _MPI
    if ( pid() == 0 )
# endif
    {
      // Rigid body data
      if ( pdata_is_open ) for (size_t k = 0; k < nrb; k++) fclose( p[k] );
      for (size_t k = 0; k < nrb; k++) fclose( d[k] );
    
      // Uzawa convergence
      fclose( UzawaCV );  
    
      // Cells, contrained cells and nb of multipliers 
      fclose( CVT );               
    }    
}




/** Compute vorticity in 3D and store it in field omega */
//----------------------------------------------------------------------------
void vorticity_3D( const vector u, vector omega ) 
//----------------------------------------------------------------------------
{  
  foreach()
    foreach_dimension()
      omega.x[] = ( ( u.z[0,1,0] - u.z[0,-1,0] )  
      	- ( u.y[0,0,1] - u.y[0,0,-1] ) ) / ( 2. * Delta );
	
  synchronize((scalar*){omega});
}




/** Computes the flow rate on the right boundary (x+ direction) */
//----------------------------------------------------------------------------
double compute_flowrate_right( const vector u, const int level ) 
//----------------------------------------------------------------------------
{
  double flowrate = 0.;
# if ADAPTIVE
    double hh = L0 / pow(2,level);
# else
    double hh = Delta;
# endif 
  double zi = 0., yj = 0., xval = X0 + 0.999999 * L0, uinter = 0.;
  int ii = 0, jj = 0;
    
  for (ii = 0; ii < pow(2,level); ii++) 
  {
    zi = 0.5 * hh + ii * hh + Z0;
    for (jj = 0; jj < pow(2,level); jj++) 
    {
      yj = 0.5 * hh + jj * hh + Y0;
      uinter = interpolate( u.x, xval, yj, zi );
      flowrate += sq(hh) * uinter;
    }
  }      
 
  return flowrate;
}




/** Computes the flow rate on the top boundary (y+ direction) */
//----------------------------------------------------------------------------
double compute_flowrate_top( const vector u, const int level ) 
//----------------------------------------------------------------------------
{
  double flowrate = 0.;
# if ADAPTIVE
    double hh = L0 / pow(2,level);
# else
    double hh = Delta;
# endif 
  double zi = 0., yval = Y0 + 0.999999 * L0, xj = 0., uinter = 0.;
  int ii = 0, jj = 0;
  
  for (ii = 0; ii < pow(2,level); ii++) 
  {
    zi = 0.5 * hh + ii * hh + Z0;
    for (jj = 0; jj < pow(2,level); jj++) 
    {
      xj = 0.5 * hh + jj * hh + X0;
      uinter = interpolate( u.y, xj, yval, zi );
      flowrate += sq(hh) * uinter;
    }
  }  
 
  return flowrate;
}




/** Computes the flow rate on the front boundary (z+ direction) */
//----------------------------------------------------------------------------
double compute_flowrate_front( const vector u, const int level ) 
//----------------------------------------------------------------------------
{
  double flowrate = 0.;  
# if ADAPTIVE
    double hh = L0 / pow(2,level);
# else
    double hh = Delta;
# endif    
  double xi = 0., yj = 0., zval = Z0 + 0.999999 * L0, uinter = 0.;
  int ii = 0, jj = 0;

  for (ii = 0; ii < pow(2,level); ii++) 
  {
    xi = 0.5 * hh + ii * hh + X0;
    for (jj = 0; jj < pow(2,level); jj++) 
    {
      yj = 0.5 * hh + jj * hh + Y0;
      uinter = interpolate ( u.z, xi, yj, zval );
      flowrate += sq(hh) * uinter;
    }
  }
 
  return flowrate;
}




/** Gets the pressure at a point and writes the value in a file */
//----------------------------------------------------------------------------
void pressure_at_point( scalar pres, FILE* ppf, const coord hv, 
	const double t )
//----------------------------------------------------------------------------
{
  double plaw = 0.;
  
  plaw = interpolate (pres, hv.x, hv.y, hv.z);

# if _MPI
    MPI_Barrier(MPI_COMM_WORLD);
    if ( plaw == nodata ) plaw = 0.;
    MPI_Barrier(MPI_COMM_WORLD);

    if ( plaw != nodata )
      mpi_all_reduce(plaw, MPI_DOUBLE, MPI_SUM);
#  endif
 
  if ( pid() == 0 ) 
  {
    fprintf( ppf, "%20.18f\t %20.18f\n", t, plaw );
    fflush( ppf );
  }
}




/** Computes and returns the total number of cells of the grid */
//----------------------------------------------------------------------------
int totalcells() 
//----------------------------------------------------------------------------
{
  int t = 0;
  foreach(reduction(+:t)) t++;  
  
  return t;
}




/** Computes and returns the total number of cells related to Distributed
Lagrange multiplier points for all rigid bodies */
//----------------------------------------------------------------------------
int total_dlmfd_cells( RigidBody const* allrbs, const size_t nrb ) 
//----------------------------------------------------------------------------
{
  int apts = 0;
  
  for (size_t k = 0; k < nrb; k++) 
  {
#   if DLMFD_INTERIORPOINTS
      apts += allrbs[k].Interior.n;
#   endif
#   if DLMFD_BOUNDARYPOINTS
      apts += allrbs[k].Boundary.n;
#   endif
  }

# if _MPI
    mpi_all_reduce( apts, MPI_INT, MPI_SUM );
# endif

  return apts;
}




/** Computes and returns the total number of Ditributed Lagrange multiplier 
points for all rigid bodies */
//----------------------------------------------------------------------------
int total_dlmfd_multipliers( RigidBody const* allrbs, const size_t nrb ) 
//----------------------------------------------------------------------------
{
  int apts = 0;
  
# if DLMFD_INTERIORPOINTS
    for (size_t k = 0; k < nrb; k++) 
      apts += allrbs[k].Interior.n;
  
#   if _MPI
      mpi_all_reduce( apts, MPI_INT, MPI_SUM );
#   endif
# endif
  
# if DLMFD_BOUNDARYPOINTS
    for (int k = 0; k < nrb; k++) 
    {
      RigidBodyBoundary const* sbb = &(allrbs[k].s);
      for (size_t j = 0; j < sbb->m; j++)
	if ( sbb->deactivated[j] == 0 ) ++apts;      
    }
# endif
  
  return apts;
}




/** Save the time and time step in a file */
//----------------------------------------------------------------------------
void save_t_dt_restart( char* dirname, double time, double deltat, double ppd )
//----------------------------------------------------------------------------
{
  if ( pid() == 0 )
  {
    char dump_name[80] = "";
    strcpy( dump_name, dirname );
    strcat( dump_name, "/t_dt_restart.res" );
    FILE* ft = fopen( dump_name, "w" );
    fprintf ( ft, "%.10e %.10e %.10e", time, deltat, ppd );
    fclose( ft );
  }  
}




/** Read the restart time and time step from a file */
//----------------------------------------------------------------------------
void read_t_restart( char* dirname, double* time, double* deltat, double* ppd )
//----------------------------------------------------------------------------
{
  char dump_name[80] = "";
  strcpy( dump_name, dirname );
  strcat( dump_name, "/t_dt_restart.res" );
  FILE* ft = fopen( dump_name, "r" );
  int result = fscanf ( ft, "%lf %lf %lf", time, deltat, ppd );
  if ( result != 3 ) 
    if ( pid() == 0 )
      printf( "WARNING: problem when reading data in file %s\n", dump_name );
  fclose( ft );  
}




/** Allocate number of rigid body dependent arrays */
//----------------------------------------------------------------------------
void allocate_np_dep_arrays( const size_t nrb, const size_t npart,
	const size_t nrbdata_, RigidBody** allrb_, size_t** RBnumToIndex_,
	double*** DLMFDtoGS_vel_, 
	double** vpartbuf_, FILE*** pdata_, const bool openpdata_, 
	FILE*** fdata_ )
//----------------------------------------------------------------------------
{
  *allrb_ = (RigidBody*) calloc( nrb, sizeof(RigidBody) );
  *RBnumToIndex_ = (size_t*) calloc( nrb + 1, sizeof(size_t) );
  if ( npart )
  {
    *DLMFDtoGS_vel_ = (double**) calloc( npart, sizeof(double*) );
    for (size_t k=0;k<npart;++k)
      (*DLMFDtoGS_vel_)[k] = (double*) calloc( 6, sizeof(double) );
  }
  else
    *DLMFDtoGS_vel_ = NULL;
  *vpartbuf_ = (double*) calloc( nrbdata_ * nrb, sizeof(double) );
  if ( openpdata_ ) *pdata_ = (FILE**) calloc( nrb, sizeof( FILE* ) );
  *fdata_ = (FILE**) calloc( nrb, sizeof( FILE* ) );    
}




/** Free number of rigid body dependent arrays */
//----------------------------------------------------------------------------
void free_np_dep_arrays( const size_t npart, RigidBody* allrb_, 
	size_t* RBnumToIndex_, double** DLMFDtoGS_vel_, double* vpartbuf_, 
	FILE** pdata_, const bool pdata_is_open, FILE** fdata_ )
//----------------------------------------------------------------------------
{
  free( allrb_ );
  free( RBnumToIndex_ );
  if ( npart )
  {
    for (size_t k=0;k<npart;++k) free( DLMFDtoGS_vel_[k] );
    free( DLMFDtoGS_vel_ );
  }
  free( vpartbuf_ );
  if ( pdata_is_open ) free( pdata_ );
  free( fdata_ );  
}




/** Compute and store velocity and explicit part of the split acceleration 
at previous times */
//----------------------------------------------------------------------------
void update_previoustimes_velocity_splitAcceleration( RigidBody* allrbs, 
	const size_t nrb, const double dt ) 
//----------------------------------------------------------------------------
{  
  for (size_t k = 0; k < nrb; k++) 
  {
    if ( allrbs[k].type != OBSTACLE )
    {
      // Save previous velocity before updating
#     if TRANSLATION
#       if B_SPLIT_EXPLICIT_ACCELERATION
          foreach_dimension()
	    allrbs[k].splitUacc.x = ( FLUID_DENSITY / allrbs[k].rho_s )
	    	* allrbs[k].M * ( allrbs[k].U.x
		- allrbs[k].Unm1.x ) / dt;
#       endif 
        allrbs[k].Unm1 = allrbs[k].U;
#     endif 
#     if ROTATION
#       if B_SPLIT_EXPLICIT_ACCELERATION
          coord Iom;
#         if dimension == 3
            Iom.x = allrbs[k].Ip[0] * allrbs[k].w.x 
	  	+ allrbs[k].Ip[3] * allrbs[k].w.y 
		+ allrbs[k].Ip[4] * allrbs[k].w.z;
            Iom.y = allrbs[k].Ip[3] * allrbs[k].w.x 
	  	+ allrbs[k].Ip[1] * allrbs[k].w.y 
		+ allrbs[k].Ip[5] * allrbs[k].w.z;
            Iom.z = allrbs[k].Ip[4] * allrbs[k].w.x 
	  	+ allrbs[k].Ip[5] * allrbs[k].w.y 
		+ allrbs[k].Ip[2] * allrbs[k].w.z;
#         else
            Iom.z = allrbs[k].Ip[2] * allrbs[k].w.z;
#         endif
          foreach_dimension()
	    allrbs[k].splitwacc.x = ( FLUID_DENSITY / allrbs[k].rho_s )
	    	* ( Iom.x - allrbs[k].Iwnm1.x ) / dt;
	  allrbs[k].Iwnm1 = Iom;
#       endif    	 
        allrbs[k].wnm1 = allrbs[k].w;
#     endif
    }
    else
    {
#     if TRANSLATION
        allrbs[k].Unm1 = allrbs[k].imposedU;
#     endif 
#     if ROTATION   	 
        allrbs[k].wnm1 = allrbs[k].imposedw;
#     endif 
    }
  }
}




/** Save the explicit part of the split acceleration of all particles for 
restart purposes in a file */
//----------------------------------------------------------------------------
void save_explicit_splitAcceleration( char* dirname, RigidBody const* allrbs, 
	const size_t nrb )
//----------------------------------------------------------------------------
{
  char dump_name[160] = "";
  strcpy( dump_name, dirname );
  strcat( dump_name, "/particle_explicit_acceleration.res" );
  FILE* ft = fopen( dump_name, "w" );

  for (size_t k = 0; k < nrb; k++)
  { 
    fprintf ( ft, "%.10e", allrbs[k].splitUacc.x );
    fprintf ( ft, " %.10e", allrbs[k].splitUacc.y );
#   if dimension == 3
      fprintf ( ft, " %.10e", allrbs[k].splitUacc.z );
      fprintf ( ft, " %.10e", allrbs[k].splitwacc.x );
      fprintf ( ft, " %.10e", allrbs[k].splitwacc.y );      
#   endif
    fprintf ( ft, " %.10e\n", allrbs[k].splitwacc.z );                        
  }  

  fclose( ft );  
}




/** Read the explicit part of the split acceleration of all particles for 
restart purposes from a file */
//----------------------------------------------------------------------------
void read_explicit_splitAcceleration( char* dirname, RigidBody* allrbs, 
	const size_t nrb )
//----------------------------------------------------------------------------
{
  char dump_name[160] = "";
  strcpy( dump_name, dirname );
  strcat( dump_name, "/particle_explicit_acceleration.res" );
  FILE* ft = fopen( dump_name, "r" );
  int result = 0;

  for (size_t k = 0; k < nrb; k++)
  { 
    result = 0;
    result += fscanf ( ft, "%lf", &(allrbs[k].splitUacc.x) );
    result += fscanf ( ft, "%lf", &(allrbs[k].splitUacc.y) );
#   if dimension == 3
      result += fscanf ( ft, "%lf", &(allrbs[k].splitUacc.z) );
      result += fscanf ( ft, "%lf", &(allrbs[k].splitwacc.x) );
      result += fscanf ( ft, "%lf", &(allrbs[k].splitwacc.y) );      
#   endif
    result += fscanf ( ft, "%lf", &(allrbs[k].splitwacc.z) ); 
#   if dimension == 3
      if ( result != 6 )    
#   else
      if ( result != 3 )            
#   endif 
        if ( pid() == 0 )
	  printf( "WARNING: problem when reading data in file %s\n", 
	  	dump_name );                          
  }  

  fclose( ft ); 
}




/** Adapt_wavelet algorithm with a maxlevel per field component */
//----------------------------------------------------------------------------
trace
astats adapt_wavelet_multimaxlevel (scalar * slist,       // list of scalars
		      double * max,         // tolerance for each scalar
		      int * maxlevel,         // maximum level of refinement
		      int minlevel = 1,     // minimum level of refinement
		      scalar * list = all)  // list of fields to update
//----------------------------------------------------------------------------
{
  scalar * ilist = list;
  
  if (is_constant(cm)) {
    if (list == NULL || list == all)
      list = list_copy (all);
    boundary (list);
    restriction (slist);
  }
  else {
    if (list == NULL || list == all) {
      list = list_copy ({cm, fm});
      for (scalar s in all)
	list = list_add (list, s);
    }
    boundary (list);
    scalar * listr = list_concat (slist, {cm});
    restriction (listr);
    free (listr);
  }

  astats st = {0, 0};
  scalar * listc = NULL;
  for (scalar s in list)
    listc = list_add_depend (listc, s);

  // refinement
  if (minlevel < 1)
    minlevel = 1;
  tree->refined.n = 0;
  static const int refined = 1 << user, too_fine = 1 << (user + 1);
  foreach_cell() {
    if (is_active(cell)) {
      static const int too_coarse = 1 << (user + 2);
      if (is_leaf (cell)) {
	if (cell.flags & too_coarse) {
	  cell.flags &= ~too_coarse;
	  refine_cell (point, listc, refined, &tree->refined);
	  st.nf++;
	}
	continue;
      }
      else { // !is_leaf (cell)
	if (cell.flags & refined) {
	  // cell has already been refined, skip its children
	  cell.flags &= ~too_coarse;
	  continue;
	}
	// check whether the cell or any of its children is local
	bool local = is_local(cell);
	if (!local)
	  foreach_child()
	    if (is_local(cell)) {
	      local = true; break;
	    }
	if (local) {
	  int i = 0;
	  static const int just_fine = 1 << (user + 3);
	  for (scalar s in slist) {
	    double emax = max[i], sc[1 << dimension];
	    int c = 0, smaxlevel = maxlevel[i++];
	    foreach_child()
	      sc[c++] = s[];
	    s.prolongation (point, s);
	    c = 0;
	    foreach_child() {
	      double e = fabs(sc[c] - s[]);
	      if (e > emax && level < smaxlevel) {
		cell.flags &= ~too_fine;
		cell.flags |= too_coarse;
	      }
	      else if ((e <= emax/1.5 || level > smaxlevel) &&
		       !(cell.flags & (too_coarse|just_fine))) {
		if (level >= minlevel)
		  cell.flags |= too_fine;
	      }
	      else if (!(cell.flags & too_coarse)) {
		cell.flags &= ~too_fine;
		cell.flags |= just_fine;
	      }
	      s[] = sc[c++];
	    }
	  }
	  foreach_child() {
	    cell.flags &= ~just_fine;
	    if (!is_leaf(cell)) {
	      cell.flags &= ~too_coarse;
	      if (level >= MAXLEVEL)
		cell.flags |= too_fine;
	    }
	    else if (!is_active(cell))
	      cell.flags &= ~too_coarse;
	  }
	}
      }
    }
    else // inactive cell
      continue;
  }
  mpi_boundary_refine (listc);
  
  // coarsening
  // the loop below is only necessary to ensure symmetry of 2:1 constraint
  for (int l = depth(); l >= 0; l--) {
    foreach_cell()
      if (!is_boundary(cell)) {
	if (level == l) {
	  if (!is_leaf(cell)) {
	    if (cell.flags & refined)
	      // cell was refined previously, unset the flag
	      cell.flags &= ~(refined|too_fine);
	    else if (cell.flags & too_fine) {
	      if (is_local(cell) && coarsen_cell (point, listc))
		st.nc++;
	      cell.flags &= ~too_fine; // do not coarsen parent
	    }
	  }
	  if (cell.flags & too_fine)
	    cell.flags &= ~too_fine;
	  else if (level > 0 && (aparent(0).flags & too_fine))
	    aparent(0).flags &= ~too_fine;
	  continue;
	}
	else if (is_leaf(cell))
	  continue;
      }
    mpi_boundary_coarsen (l, too_fine);
  }
  free (listc);

  mpi_all_reduce (st.nf, MPI_INT, MPI_SUM);
  mpi_all_reduce (st.nc, MPI_INT, MPI_SUM);
  if (st.nc || st.nf)
    mpi_boundary_update (list);

  if (list != ilist)
    free (list);
  
  return st;
}
