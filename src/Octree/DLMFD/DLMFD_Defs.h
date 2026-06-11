/** 
# Definitions and global variables
*/


//----------------------------------------------------------------------------
// DEFINTIONS
//----------------------------------------------------------------------------
# ifndef FLUID_DUMP_FILENAME
#   define FLUID_DUMP_FILENAME "Savings/dump"
# endif
# ifndef DUMP_DIR
#   define DUMP_DIR "Savings"
# endif

# ifndef RESULT_DIR
#   define RESULT_DIR "Res"
# endif
# ifndef RESULT_RIGIDBODY_VP_ROOTFILENAME
#   define RESULT_RIGIDBODY_VP_ROOTFILENAME "rigidbody_data"
# endif
# ifndef RESULT_RIGIDBODY_HYDROFAT_ROOTFILENAME
#   define RESULT_RIGIDBODY_HYDROFAT_ROOTFILENAME "rigidbody_hydroFaT"
# endif
# ifndef RESULT_FLUID_ROOTFILENAME
#   define RESULT_FLUID_ROOTFILENAME "fluid"
# endif

# ifndef RESULT_DIR
#   define RESULT_DIR "Res"
# endif

# ifndef CONVERGE_UZAWA_FILENAME 
#   define CONVERGE_UZAWA_FILENAME "convergence_uzawa_velocity.dat"
# endif
# ifndef DLMFD_CELLS_FILENAME 
#   define DLMFD_CELLS_FILENAME "dlmfd_cells.dat"
# endif
# ifndef DLMFD_PERF_FILENAME 
#   define DLMFD_PERF_FILENAME "dlmfd_perf.dat"
# endif

# ifndef ROUNDDOUBLECOEF
#   define ROUNDDOUBLECOEF (1.e-4)
# endif

# ifndef TINTERVALOUTPUT
#   define TINTERVALOUTPUT (1.)
# endif

# ifndef SIMUTIMEINTERVAL
#   define SIMUTIMEINTERVAL (1.)
# endif

# ifndef FIGS_DIR
#   define FIGS_DIR "Figs"
# endif

# ifndef INITIALGRIDADAPTIVE_NEWMETHOD
#   define INITIALGRIDADAPTIVE_NEWMETHOD 1
# endif

# ifndef IMPOSED_PERIODICFLOW
#   define IMPOSED_PERIODICFLOW 0
# endif

# ifndef IMPOSED_PERIODICFLOW_TYPE
#   define IMPOSED_PERIODICFLOW_TYPE 0 // 0 for pressure and 1 for flow rate
# endif

# ifndef IMPOSED_PERIODICFLOW_DIRECTION
#   define IMPOSED_PERIODICFLOW_DIRECTION 0 
# endif

# ifndef DLMFD_PROB_AFTER_NAVIERSTOKES
#   define DLMFD_PROB_AFTER_NAVIERSTOKES 0
# endif

# ifndef RIGIDBODIES_AS_FIXED_OBSTACLES
#   define RIGIDBODIES_AS_FIXED_OBSTACLES 0
# endif

# ifndef LEVELDIFF_FLAG_U
#   define LEVELDIFF_FLAG_U 0
# endif

# ifndef FLAG_ADAPT_CRIT
#   define FLAG_ADAPT_CRIT (1.E-16)
# endif

# ifndef UX_ADAPT_CRIT
#   define UX_ADAPT_CRIT (1.E-2)
# endif

# ifndef UY_ADAPT_CRIT
#   define UY_ADAPT_CRIT (1.E-2)
# endif

# ifndef UZ_ADAPT_CRIT
#   define UZ_ADAPT_CRIT (1.E-2)
# endif

# ifndef GRAVITY_VECTOR
#   if dimension == 2
#     define GRAVITY_VECTOR ((coord){0.,0.})
#   else 
#     define GRAVITY_VECTOR ((coord){0.,0.,0.})
#   endif
# endif

# ifndef LAMBDA2
#   define LAMBDA2 0
# endif

# ifndef VORTICITY
#   define VORTICITY 0
# endif

# ifndef PARAVIEW_VTU
#   define PARAVIEW_VTU 0
# endif

# ifndef PARAVIEW_HTG
#   define PARAVIEW_HTG 0
# endif

# ifndef PARAVIEW_SCALAR_LIST
#   define PARAVIEW_SCALAR_LIST p
# endif

# ifndef PARAVIEW_VECTOR_LIST
#   define PARAVIEW_VECTOR_LIST u
# endif

# ifndef BVIEW
#   define BVIEW 0
# endif

# ifndef BVIEW_LIST
#   define BVIEW_LIST u,p
# endif

# ifndef PARAVIEW_DLMFD_INTPTS
#   define PARAVIEW_DLMFD_INTPTS 0
# endif

# if PARAVIEW_DLMFD_INTPTS
#   ifndef PARAVIEW_DLMFD_INTPTS_FILENAME
#     define PARAVIEW_DLMFD_INTPTS_FILENAME "dlmfd_interior_points"
#   endif      
# endif

# ifndef PARAVIEW_DLMFD_BNDPTS
#   define PARAVIEW_DLMFD_BNDPTS 0
# endif

# if PARAVIEW_DLMFD_BNDPTS
#   ifndef PARAVIEW_DLMFD_BNDPTS_FILENAME
#     define PARAVIEW_DLMFD_BNDPTS_FILENAME "dlmfd_boundary_points"
#   endif      
# endif

# if PARAVIEW_VTU || PARAVIEW_HTG || PARAVIEW_DLMFD_INTPTS \
	|| PARAVIEW_DLMFD_BNDPTS
#   define PARAVIEW 1
# else 
#   define PARAVIEW 0
# endif

# define BGHOSTS 2
# define BSIZE 128

# ifndef TRANSLATION
#   define TRANSLATION 1
# endif

# ifndef ROTATION
#   define ROTATION 1
# endif

# ifndef DLM_ALPHA_COUPLING
#   define DLM_ALPHA_COUPLING 0
# endif

# ifndef DLM_UZAWA_TOL
#   define DLM_UZAWA_TOL 1.e-5
# endif

# ifndef DLMFD_BOUNDARYPOINTS
#   define DLMFD_BOUNDARYPOINTS 1
# endif

# ifndef DLMFD_INTERIORPOINTS
#   define DLMFD_INTERIORPOINTS 1
# endif

# ifndef DLMFD_OPT
#   define DLMFD_OPT 1     // use optimized version of DLMFD algorithm
# endif

# ifndef RIGIDBODY_VERBOSE
#   define RIGIDBODY_VERBOSE 0     // print rigid body features
# endif

# if ( TRANSLATION && ROTATION )
#   if dimension == 3
#     define NRBDATA 6
#   else
#     define NRBDATA 3
#   endif     
# elif TRANSLATION
#   if dimension == 3 
#     define NRBDATA 3
#   else
#     define NRBDATA 2
#   endif 
# elif ROTATION
#   if dimension == 3 
#     define NRBDATA 3
#   else
#     define NRBDATA 1
#   endif 
# else
#   define NRBDATA 0  
# endif

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

/** Define the width of the boundary layer with multiple max levels */
# if LEVELDIFF_FLAG_U
#   ifndef BOUNDARY_LAYER_THICKNESS_COEF
#     define BOUNDARY_LAYER_THICKNESS_COEF 4. 
#   endif
# endif  

# ifndef PARAVIEW_DATATYPE_DOUBLE
#   define PARAVIEW_DATATYPE_DOUBLE 0 // 1 for double and 0 for float
# endif

# if ( PARAVIEW_DATATYPE_DOUBLE == 1 )
#   define PARAVIEW_DATATYPE double
#   define PARAVIEW_DATANAME "Float64"
# else
#   define PARAVIEW_DATATYPE float
#   define PARAVIEW_DATANAME "Float32"
# endif

# ifndef PARAVIEW_BINFILE
#   define PARAVIEW_BINFILE 1
# endif
# ifndef PARAVIEW_VTU_MPIIO_WRITER 
#   define PARAVIEW_VTU_MPIIO_WRITER 0
# endif

# if PARAVIEW_VTU_MPIIO_WRITER && !_MPI
#   undef PARAVIEW_VTU_MPIIO_WRITER
#   define PARAVIEW_VTU_MPIIO_WRITER 0
# endif

# define DYNARRAYBLOCK 128

# ifndef DLMBLOCK 
#   define DLMBLOCK 5000
# endif

# if dimension == 2
#   define NDOFSTENCIL 9
# else
#   define NDOFSTENCIL 27
# endif

/** Split explicit acceleration treatment in case of particles are lighter 
than the fluid or neutrally buoyant */
# ifndef B_SPLIT_EXPLICIT_ACCELERATION
#   define B_SPLIT_EXPLICIT_ACCELERATION 0
# endif




//----------------------------------------------------------------------------
// STRUCTURES AND ENUMERATIONS
//----------------------------------------------------------------------------
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
  ELLIPSOID,
  HEXAGONALPRISM
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
  coord* outwardnormalvector;
  bool* deactivated;
  int m;
} RigidBodyBoundary;


/** Additional geometric parameters for polygons/polyhedrons */
typedef struct {
  int ncorners, nfaces;
  coord* cornersCoord;
# if LEVELDIFF_FLAG_U
    coord* cornersCoordExp;  
# endif  
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


/** Structure for a dynamic array of unsigned integers */
typedef struct {
  size_t n;
  size_t nm;   
  size_t* elem; 
} dynUIarray;


/** Structure for a dynamic array of pointers to double */
typedef struct {
  size_t n;
  size_t nm;   
  double** elem; 
} dynPDBarray;


/** Structure and routines to store pointers and coefficients from the
initialization of the Uzawa algorithm and use them over the iterative process
to compute faster M_u^T*w = <DLM_w, v>_P(t) over the boundary points */
typedef struct {
  int n;
  int nm;
  double** qux;
  double** quy;
  double** dlmwx;
  double** dlmwy;
# if dimension == 3  
    double** quz;
    double** dlmwz;
# endif   
  double* weight; 
} BPFastLoop_LambdaMom;


/** Structure and routines to store pointers and coefficients from the
initialization of the Uzawa algorithm and use them over the iterative process
to compute faster M_u*tu = <alpha, tu>_P(t) over the boundary points */
typedef struct {
  int nv;
  int nvm;
  double** dlmvx;
  double** dlmvy;
  int* ndof;
  int ntu;
  int ntum;      
  double** tux;
  double** tuy;
# if dimension == 3 
    double** dlmvz; 
    double** tuz;
# endif 
  double* weight; 
} BPFastLoop_ResU;




//----------------------------------------------------------------------------
// BASILISK FIELDS
//----------------------------------------------------------------------------
vector DLM_lambda[];
scalar DLM_Flag[];
scalar DLM_FlagMesh[];
# if LEVELDIFF_FLAG_U
    scalar DLM_FlagMaxLev[];
# endif    
vector DLM_Index[];
vector DLM_CX_NCX[];
vector DLM_PeriodicRefCenter[];
vector DLM_r[];
vector DLM_w[];
vector DLM_v[];
vector DLM_qu[];
vector DLM_tu[];
# if DLM_ALPHA_COUPLING
    vector DLM_explicit[];
# endif
# if EMBED
    scalar rhov[];
    face vector muv[];
    face vector alphav[];
# endif




//----------------------------------------------------------------------------
// STANDARD GLOBAL VARIABLES
//----------------------------------------------------------------------------
double deltau;
int restarted_simu = 0;
# if PARAVIEW_VTU
    char vtu_field_times_series[100000] = "";
# endif
# if PARAVIEW_HTG
    char htg_field_times_series[100000] = "";
# endif
# if PARAVIEW_DLMFD_BNDPTS
    char vtk_bndpts_times_series[100000] = "";      
# endif 
# if PARAVIEW_DLMFD_INTPTS
    char vtk_intpts_times_series[100000] = "";        
# endif 
int init_cycle_number = 0;
double maxtime = 0.;
double trestart = 0.;
double dtrestart = 0.;
bool save_data_restart = false;
scalar u_previoustime[];
double imposed_periodicpressuredrop = 0.;
double imposed_periodicflowrate = 0.;
size_t nbRigidBodies = 0;
size_t nbParticles = 0;
size_t NbObstacles = 0;
size_t nbReferenceRigidBodies = 0; 
RigidBody* allRigidBodies = NULL;
RigidBody* ReferenceRigidBodies = NULL;
size_t* RBnumToIndex;
double** DLMFDtoGS_vel = NULL;
double* vpartbuf = NULL;
FILE** pdata = NULL;
FILE** fdata = NULL;
FILE* converge = NULL;
FILE* cellvstime = NULL;
dynUIarray deactivatedBPindices;
dynPDBarray deactivatedIndexFieldValues;
AABB local_domain;




//----------------------------------------------------------------------------
// TIMING AND STATISTICS
//----------------------------------------------------------------------------
# include "utils.h"
timing DLMFD_UzawaTiming = {0.};
timing DLMFD_ConstructionTiming = {0.};
DLMFDptscells allDLMFDptscells;
