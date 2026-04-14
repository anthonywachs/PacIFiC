/**
# Triangulated mesh of a capsule

In this file, we implement a Lagrangian mesh for the front-tracking method,
meant to track the position and compute the stresses of an elasitc membrane.
*/

#ifndef dimension
  #define dimension 3
#endif


#ifndef MULT_GRID
  #define MULT_GRID 0
#endif

/*The stencil type is chosen as 3 if not 5 */
#ifndef STENCIL_TYPE
  #define STENCIL_TYPE 3
#endif

/* Repulsive lubrication force to avoid overlapping */
#ifndef LUBR_FORCE
  #define LUBR_FORCE 0
#endif

#ifndef LUBR_VEL
  #define LUBR_VEL 0
#endif

/*Create the Index_lag*/
scalar Index_lagnode[];
vector Index_lag_id[];

/*Create the dimension ratios for non-cubic domain*/
#ifndef L0_X
  #define L0_X L0
#endif
#ifndef L0_Y
  #define L0_Y L0
#endif
#ifndef L0_Z
  #define L0_Z L0
#endif

coord L0_ratio = {L0_X/L0, L0_Y/L0, L0_Z/L0};

/**
## Structure of the mesh

In the Lagrangian mesh, each node is assigned several attributes:

* ```pos```, the node coordinates
* ```lagVel```, the node velocity vector
* ```normal```, the vector normal to the membrane at the node coordinates
* ```curv```, the (mean) curvature of the membrane at the node coordinates
* ```ref_curv```, the reference (mean) curvature of the membrane at the node coordinates (default is zero)
* ```lagForce```, an elastic and/or bending force exerted by the membrane on its surrounding fluid
* ```stencil```, a [Cache](http://basilisk.fr/src/grid/tree.h#82) structure used for averaging the neighboring velocities and spreading the Lagrangian force as a body force on the neighboring Eulerian nodes
* in case of MPI simulations, ```pid```, the rank of the processor owning the Eulerian cell which contains the lagNode
* ```edge_ids```, the IDs of its connecting edges: 2 in 2D, up to 6 in 3D (because every considered shape is derived by subdividing in icosahedron, leading to 5 or 6 neighbors only).

In case of 3D simulations, other attributes are introduced:

* ```nb_neighbors```, the number of neighbors of the node, which should only be 5 or 6
* ```neighbor_ids```, the IDs of the node neighbors
* ```nb_triangles```, the number of triangles having the node as a vertex
* ```triangle_ids```, the ids of the above triangles
* ```gcurv```, the Gaussian curvature of the membrane at the node coordinates
* ```nb_fit_iterations```, the number of iterations needed to compute the membrane curvature and normal vector to the desired convergence threshold.

*/
typedef struct lagNode {
  coord pos;
  coord lagVel;
  coord normal;
  double curv;
  double ref_curv;
  coord lagForce;
  Cache stencil;
  Cache eulcell;
  #if _MPI
    int pid;
  #endif
  #if dimension < 3
    int edge_ids[2];
  #else
    int nb_neighbors;
    int neighbor_ids[6];
    int edge_ids[6];
    int nb_triangles;
    int triangle_ids[6];
    double gcurv;
    int nb_fit_iterations;
  #endif
} lagNode;

/** We specify the size of the 3x3(x3) or 5x5(x5) stencil in 2D(3D). */
#if dimension < 3
  #if STENCIL_TYPE == 3
    #define STENCIL_SIZE 9
  #else 
    #define STENCIL_SIZE 25
  #endif
#else
  #if STENCIL_TYPE == 3
    #define STENCIL_SIZE 27
  #else 
    #define STENCIL_SIZE 125
  #endif
#endif

/** Similarly, the edges of the mesh are assigned:

* ```node_ids```, the IDs of the two nodes they connect
* In case of 3D simulations, ```triangle_ids```, the IDs of the two triangles they separate
* ```l0```, the length of the edge in the initial stress-free configuration
* ```length```, the current length of the edge
* ```normal```, a vector normal to the membrane at the edge midpoint.

*/
typedef struct Edge {
  int node_ids[2];
  #if dimension > 2
    int triangle_ids[2];
  #endif
  double l0;
  double length;
  coord normal;
} Edge;

/** In case of 3D simulations, we define triangle faces. Each ```Triangle``` structure has the following attributes:

* ```node_ids```, the IDs of the triangle vertices
* ```edge_ids```, the IDs of the triangle edges
* ```area```, the current area of the triangle
* ```normal```, the normal vector to the triangle pointing outside the membrane
* ```centroid```, the coordinates of the triangle centroid
* ```refShape```, the coordinates of the vertices in the Common Plane in the stress-free configuration. By convention the first vertex is always placed at $(0, 0, 0)$, so only the coordinates of the second and third vertex are stored in ```refShape```
* ```sfc```, the shape function coefficients used in the finite elements method and computed in [store_initial_configuration](elasticity-ft.h#store_initial_configuration). There are two coefficients per vertex, hence the $3 \times 2$ array structure of ```sfc```.

 */
#if dimension > 2
  typedef struct Triangle {
    int node_ids[3];
    int edge_ids[3];
    double area;
    coord normal;
    coord centroid;
    coord refShape[2];
    double sfc[3][2]; // sfc for "shape function coefficients"
    double stretch[2];
    double tension[2];
  } Triangle;
#endif

/** The ```lagMesh``` structure contains arrays of the previously introduced nodes, edges and triangles. It defines an unstructured mesh, the membrane of our capsule. Its attributes are:

*  ```nln``` the number of Lagrangian nodes
* ```nodes```, the array of Lagrangian nodes
* ```nle```, the number of Lagrangian edges
* ```edges```, the array of Lagrangian edges
* In case of 3D simulations:
    * ```nlt```, the number of Lagrangian triangles
    * ```triangles```, the array of Lagrangian triangles
* ```updated_stretches```, a boolean used to check if the current length of the edges has been updated since the last advection of the Lagrangian nodes
* ```updated_normals```, a similar boolean telling if the nodal normal vectors should be recomputed
* ```updated_curvatures```, a last boolean telling if the nodal curvatures should be recomputed.
* ```isactive```, a boolean indicating if the capsule exists in the flow (useful
when capsules are introduced during a simulation)
*/

typedef struct lagMesh {
  int cap_id;
  int cap_type;
  double cap_es;
  double cap_radius;
  int nln;
  lagNode* nodes;
  int nle;
  Edge* edges;
  #if dimension > 2
    int nlt;
    Triangle* triangles;
  #endif
  coord centroid;
  coord ang_vel;
  double initial_volume;
  double volume;
  double circum_radius;
  double taylor_deform;
  bool updated_stretches;
  bool updated_normals;
  bool updated_curvatures;
  bool isactive;
} lagMesh;

/** We denote by ```NCAPS``` the number of Lagrangian meshes, or capsules, in
the simulation. It is one by default. */
#ifndef NCAPS
  #define NCAPS 1
#endif
#ifndef RESTART_CASE
  #define RESTART_CASE 0
#endif

/** The Lagrangian mesh is accessible in the code thanks to the structure
below, which is simply an array of Lagrangian meshes (useful when several of
them are considered). The macro $CAPS(k)$ can be used as a shortcut to access the
$k^{th}$ membrane. */
typedef struct Capsules {
  lagMesh caps[NCAPS];
  int nbcaps;
} Capsules;
Capsules allCaps;
#define CAPS(i) (allCaps.caps[i])


/**
## Initialization, memory management and useful macros.
*/
void initialize_empty_capsule(lagMesh* mesh) {
  mesh->cap_es = 1.;
  mesh->cap_radius = 1.;
  mesh->cap_id = -1;
  mesh->cap_type = -1;
  mesh->nln = 0;
  mesh->nle = 0;
  mesh->nodes = NULL;
  mesh->edges = NULL;
  #if dimension > 2
    mesh->nlt = 0;
    mesh->triangles = NULL;
  #endif
  mesh->updated_stretches = false;
  mesh->updated_normals = false;
  mesh->updated_curvatures = false;
  mesh->isactive = false;
}

void free_one_caps(lagMesh* mesh) {
  for(int i=0; i<mesh->nln; i++) free(mesh->nodes[i].stencil.p);
  for(int i=0; i<mesh->nln; i++) free(mesh->nodes[i].eulcell.p);
  free(mesh->nodes);
  free(mesh->edges);
  #if dimension > 2
    free(mesh->triangles);
  #endif
}

void free_all_caps(Capsules* caps) {
  for(int i=0; i<caps->nbcaps; i++)
    if (CAPS(i).isactive)
      free_one_caps(&(caps->caps[i]));
}

/**
## Adding capsules after restart

In case we add membranes to a simulation after restart, we need to call
the function below. */
void initialize_capsules() {
  allCaps.nbcaps = NCAPS;
  for(int i=0; i<NCAPS; i++) {
    if (CAPS(i).isactive) initialize_empty_capsule(&CAPS(i));
  }
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face() a.x[] = 0.;
  }
}

void initialize_capsule_stencils(lagMesh* mesh) {
  for(int j=0; j<mesh->nln; j++) {
    mesh->nodes[j].stencil.n = STENCIL_SIZE;
    mesh->nodes[j].stencil.nm = STENCIL_SIZE;
    mesh->nodes[j].stencil.p = (Index*) malloc(STENCIL_SIZE*sizeof(Index));
    mesh->nodes[j].eulcell.n = 1;
    mesh->nodes[j].eulcell.nm = 1;
    mesh->nodes[j].eulcell.p = (Index*) malloc(sizeof(Index));
  }
  
}

void initialize_all_capsules_stencils() {
  for(int i=0; i<NCAPS; i++)
    if (CAPS(i).isactive)
      initialize_capsule_stencils(&CAPS(i));
}

void initialize_active_capsule(lagMesh* mesh, int cap_id, int cap_type) {
  initialize_empty_capsule(mesh);
  mesh->cap_type = cap_type;
  mesh->cap_id = cap_id;
  mesh->isactive = true;
  initialize_capsule_stencils(mesh);
}

/** The next few macros are useful to compute signed distances and averages
across periodic boundaries. We assume for this purpose that the length of
the edges are less that half the domain size, which in practice should always
be the case. */
#define ACROSS_PERIODIC(a,b,L) (fabs(a - b) > L/2.)
#define PERIODIC_1DIST(a,b,L) (fabs(a - L - b) > L/2. ? a + L - b : a - L - b)
#define GENERAL_1DIST(a,b,L) (ACROSS_PERIODIC(a,b,L) ? PERIODIC_1DIST(a,b,L) : a - b)
#define PERIODIC_1DAVG(a,b,L) (fabs(a - L - b) > L/2. ? a + L + b : a - L + b)
#define GENERAL_1DAVG(a,b,L) (ACROSS_PERIODIC(a,b,L) ? PERIODIC_1DAVG(a,b,L) : a + b)
#define GENERAL_SQNORM(a,b) (sq(GENERAL_1DIST(a.x, b.x, L0)) + \
  sq(GENERAL_1DIST(a.y, b.y, L0*Dimensions.y/Dimensions.x)) + sq(GENERAL_1DIST(a.z, b.z, L0*Dimensions.z/Dimensions.x)))

#if dimension < 3
  #define cnorm(a) (sqrt(sq(a.x) + sq(a.y)))
  #define cdot(a,b) (a.x*b.x + a.y*b.y)
#else
  #define cnorm(a) (sqrt(sq(a.x) + sq(a.y) + sq(a.z)))
  #define cdot(a,b) (a.x*b.x + a.y*b.y + a.z*b.z)
#endif

/**
## Advection of the mesh

The advection have to be compatible with MPI, so we start by including the MPI
communication functions. We also include [reg-dirac.h](reg-dirac.h), which
implements the interpolation of the Eulerian velocities onto the Lagrangian
nodes.
*/

/** By default, the mesh is advected using a second-order two-step Runge Kutta
scheme. If the following macro is set to 0, a first-order forward Euler schme
is used instead. */
#ifndef ADVECT_LAG_RK2
  #define ADVECT_LAG_RK2 1
#endif

#if _MPI
  #include "capsule-ft-mpi.h"
#endif
#include "ibm-ft.h"
#ifndef CONSERVE_VOLUME
  #define CONSERVE_VOLUME 1
#endif

/**
A lot of triangulated mesh-related functions are defined in a separate header
file.
*/
#include "geometry-ft.h"

#if CONSERVE_VOLUME
  #include "volume-conservation-ft.h"
#endif

/* Utilities for different numerical simulations */
#include "plugins-ft.h"



//----------------------------------------------------------------------------
trace void synchronize (scalar * list)
//----------------------------------------------------------------------------
{
  for (scalar s in list)
    s.dirty = true;
  boundary(list);
}


/**
The function below advects each Lagrangian node by
interpolating the velocities around the node of interest. By default, a
second-order Runge Kutta scheme is used. By setting the macro
```ADVECT_LAG_RK2``` to 0, a simple forward Euler scheme is used.
*/
trace
void advect_lagMesh(lagMesh* mesh) {

  #if !(ADVECT_LAG_RK2)
    for(int i=0; i < mesh->nln; i++) {
      foreach_dimension() {
        mesh->nodes[i].pos.x += dt*mesh->nodes[i].lagVel.x;
      }
    }
  #else
    lagMesh buffer_mesh;
    buffer_mesh.isactive = true;
    buffer_mesh.nln = mesh->nln;
    buffer_mesh.nodes = malloc(mesh->nln*sizeof(lagNode));
    for(int i=0; i<mesh->nln; i++) {
      // Step 1 of RK2
      foreach_dimension()
        buffer_mesh.nodes[i].pos.x = mesh->nodes[i].pos.x +
          .5*dt*mesh->nodes[i].lagVel.x;
    }
    correct_lag_pos(&buffer_mesh);
    for(int j=0; j<buffer_mesh.nln; j++) {
      buffer_mesh.nodes[j].stencil.n = STENCIL_SIZE;
      buffer_mesh.nodes[j].stencil.nm = STENCIL_SIZE;
      buffer_mesh.nodes[j].stencil.p = malloc(STENCIL_SIZE*sizeof(Index));
      buffer_mesh.nodes[j].eulcell.n = 1;
      buffer_mesh.nodes[j].eulcell.nm = 1;
      buffer_mesh.nodes[j].eulcell.p = malloc(sizeof(Index));
    }
    
    generate_lag_stencils_one_caps(&buffer_mesh, true);
    eul2lag(&buffer_mesh);
    for(int i=0; i<mesh->nln; i++) {
      // Step 2 of RK2
      foreach_dimension()
        mesh->nodes[i].pos.x += dt*buffer_mesh.nodes[i].lagVel.x;
    }
    for(int i=0; i<buffer_mesh.nln; i++) free(buffer_mesh.nodes[i].stencil.p);
    for(int i=0; i<buffer_mesh.nln; i++) free(buffer_mesh.nodes[i].eulcell.p);
    free(buffer_mesh.nodes);
  #endif

  correct_lag_pos(mesh);
  #if CONSERVE_VOLUME
    enforce_optimal_volume_conservation(mesh);
  #endif
  comp_centroid(mesh);
  comp_volume(mesh);
  comp_capsule_geodynamics(mesh);
}





/**
## Putting the pieces together

Below, we call the above functions at the appropriate time using the Basilisk
event syntax.

We start by creating empty Lagrangian meshes, and allocating an acceleration
field in case it isn't done yet by another Basilisk solver.
*/
event defaults (i = 0) {
  allCaps.nbcaps = NCAPS;
  // #if (RESTART_CASE == 0)
    for(int i=0; i<NCAPS; i++) {
      initialize_empty_capsule(&CAPS(i));
    }
  // #endif
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face() a.x[] = 0.;
  }
}



/*Repulsive lubrication nodal force*/
void repulsive_vel() 
{
  /*Compute the cell size in the grid*/
  #if MULT_GRID == 1   
    double delta = (L0/(1 << grid->maxdepth)/Dimensions.x);
  #else
    double delta = (L0/(1 << grid->maxdepth));
  #endif

  
  /*The value of K_lub is up to the */
  // double K_lub = 0.001/(E_S);

  for(int i = 0; i < NCAPS; i++) {
    if (CAPS(i).isactive) 
    {
      lagMesh* mesh = &(CAPS(i));

      for(int j=0; j<mesh->nln; j++) 
      { 
        foreach_cache(mesh->nodes[j].eulcell)
        {
          // int lagnode_id = (int)Index_lag_id.x[];
          int lagnode_id = j;
          coord lub_vel = {0};  
          double K_lub = 0.;
          foreach_dimension()
            K_lub += sq(mesh->nodes[lagnode_id].lagVel.x);
          K_lub = sqrt(K_lub)*0.25; //half the velocity as coefficient

          if(point.level>-1)
          {        
              coord lagpt = {0};
              lagpt.x = mesh->nodes[lagnode_id].pos.x;
              lagpt.y = mesh->nodes[lagnode_id].pos.y;
              lagpt.z = mesh->nodes[lagnode_id].pos.z;
                 
              foreach_neighbor()
              {
                if(point.level >-1)
                {
                  if(((int)Index_lagnode[] > -1) && ((mesh->cap_id) != (int)Index_lagnode[])) 
                  {        
                    coord checkpt = {0};
                    checkpt.x = CAPS((int)Index_lagnode[]).nodes[(int)Index_lag_id.x[]].pos.x;
                    checkpt.y = CAPS((int)Index_lagnode[]).nodes[(int)Index_lag_id.x[]].pos.y;
                    checkpt.z = CAPS((int)Index_lagnode[]).nodes[(int)Index_lag_id.x[]].pos.z;

                    coord lub_dir = {0};
                    
                    double lub_norm = sqrt(GENERAL_SQNORM(lagpt, checkpt));
            
                    foreach_dimension() lub_dir.x = GENERAL_1DIST(lagpt.x, checkpt.x, L0*L0_ratio.x)/lub_norm;
                    if(lub_norm < 2.*delta)
                    {
                      foreach_dimension() lub_vel.x += lub_dir.x * K_lub * (sq(2.*delta/lub_norm) - 1.);
                    }
                  }
                }
              }
            /** The velocity of the node is adjusted as the repulsive effect applies. */
            foreach_dimension() mesh->nodes[lagnode_id].lagVel.x += 0.5*lub_vel.x; 
      
            
            /*A special case where the two nodes from different caps lie in the same cell, we push the node of the other capsule as well */
            if(((int)Index_lagnode[] > -1) && ((int)Index_lag_id.y[] > -1) && ((mesh->cap_id) != (int)Index_lag_id.y[]))
            {
              coord checkpt = {0};
              checkpt.x = CAPS((int)Index_lag_id.y[]).nodes[(int)Index_lag_id.z[]].pos.x;
              checkpt.y = CAPS((int)Index_lag_id.y[]).nodes[(int)Index_lag_id.z[]].pos.y;
              checkpt.z = CAPS((int)Index_lag_id.y[]).nodes[(int)Index_lag_id.z[]].pos.z;

              coord lub_dir = {0};
              double lub_norm = sqrt(GENERAL_SQNORM(lagpt, checkpt));
              foreach_dimension() lub_dir.x = GENERAL_1DIST(lagpt.x, checkpt.x, L0*L0_ratio.x)/lub_norm;
              if(lub_norm < 2.*delta)
              {                 
                CAPS((int)Index_lag_id.y[]).nodes[(int)Index_lag_id.z[]].lagVel.x -= 0.5*lub_dir.x * K_lub * (sq(2.*delta/lub_norm) - 1.);
                CAPS((int)Index_lag_id.y[]).nodes[(int)Index_lag_id.z[]].lagVel.y -= 0.5*lub_dir.y * K_lub * (sq(2.*delta/lub_norm) - 1.);
                CAPS((int)Index_lag_id.y[]).nodes[(int)Index_lag_id.z[]].lagVel.z -= 0.5*lub_dir.z * K_lub * (sq(2.*delta/lub_norm) - 1.);
              }
            }
            
          }
        }
      }
    }
  }
}


/** Below, we advect each Lagrangian node using the interpolated Eulerian
velocities. We also use this loop as an opportunity to
re-initialize the Lagrangian forces to zero. */

coord proc_max = {-HUGE, -HUGE, -HUGE};
coord proc_min = {HUGE, HUGE, HUGE};

event tracer_advection(i++) {  

  /* Distribute velocity to the lagNodes */
  for(int i=0; i<NCAPS; i++) 
  {
      if (CAPS(i).isactive) 
        eul2lag(&CAPS(i));
  }   

 
  /*We synchronize the eul field and make sure that it is updated before applying repulsive velocity */
  
#if (LUBR_VEL==1)
  synchronize({Index_lagnode, Index_lag_id}); 
  repulsive_vel();
#endif 
  /**
  In case of parallel simulations, we communicate the Lagrangian velocity
  so that all processes have the same Lagrangian velocities.
  */
  reduce_alllagVel();


  /* Advection of the lagNode */
  for(int i=0; i<NCAPS; i++) {
    if (CAPS(i).isactive) {
      advect_lagMesh(&CAPS(i));
      for(int j=0; j<CAPS(i).nln; j++)
        foreach_dimension() CAPS(i).nodes[j].lagForce.x = 0.;
    }
  }

  /* Compute borders of the curren proc */
  compute_proc_borders(&proc_max, &proc_min);

  /*Clean the index field before generating the stencils*/
  foreach()
  {
    if (cm[] > 1.e-20) 
    { Index_lagnode[] = -1;
      foreach_dimension() Index_lag_id.x[] = -1;
    }
  }

  /* Generate new stencils in corresponding procs */
  for(int i=0; i<NCAPS; i++) {
    if (CAPS(i).isactive)
      if(is_capsule_in_boundingbox(proc_max, proc_min, &CAPS(i))) 
        generate_lag_stencils_one_caps(&CAPS(i), true);
  }

}


/*Repulsive lubrication nodal force*/
void lubrication_force() 
{
  /*Compute the cell size in the grid*/
  #if MULT_GRID == 1   
    double delta = (L0/(1 << grid->maxdepth)/Dimensions.x);
    
  #else
    double delta = (L0/(1 << grid->maxdepth));
  #endif

  /*The value of K_lub is up to the */
  // double K_lub = 0.001/(E_S);

  for(int i = 0; i < NCAPS; i++) {
    if (CAPS(i).isactive) 
    {
      lagMesh* mesh = &(CAPS(i));

      for(int j=0; j<mesh->nln; j++) 
      { 
        foreach_cache(mesh->nodes[j].eulcell)
        {
          // int lagnode_id = (int)Index_lag_id.x[];
          int lagnode_id = j;
          coord lub_force = {0};  
          double K_lub = 0.;
          foreach_dimension()
            K_lub += sq(mesh->nodes[lagnode_id].lagForce.x);
          K_lub = 2.*sqrt(K_lub);

          if(point.level>-1)
          {        
              coord lagpt = {0};
              lagpt.x = mesh->nodes[lagnode_id].pos.x;
              lagpt.y = mesh->nodes[lagnode_id].pos.y;
              lagpt.z = mesh->nodes[lagnode_id].pos.z;
                 
              foreach_neighbor()
              {
                if(point.level >-1)
                {
                  if(((int)Index_lagnode[] > -1) && ((mesh->cap_id) != (int)Index_lagnode[])) 
                  {        
                    coord checkpt = {0};
                    checkpt.x = CAPS((int)Index_lagnode[]).nodes[(int)Index_lag_id.x[]].pos.x;
                    checkpt.y = CAPS((int)Index_lagnode[]).nodes[(int)Index_lag_id.x[]].pos.y;
                    checkpt.z = CAPS((int)Index_lagnode[]).nodes[(int)Index_lag_id.x[]].pos.z;

                    coord lub_dir = {0};
                    
                    double lub_norm = sqrt(GENERAL_SQNORM(lagpt, checkpt));
                    foreach_dimension() lub_dir.x = GENERAL_1DIST(lagpt.x, checkpt.x, L0*L0_ratio.x)/lub_norm;
                    if(lub_norm < 2*delta)
                    {
                      foreach_dimension() lub_force.x += lub_dir.x * K_lub * (sq(2*delta/lub_norm) - 1.);
                    }
                  }  
                }
              }
            /** The lubrication force is ready to be added to the Lagrangian force of the considered node. */
            foreach_dimension() mesh->nodes[lagnode_id].lagForce.x += lub_force.x;   
          }
        }
      }
    }
  }
}


/** In the acceleration event, we transfer the Lagrangian forces to the fluid
using a regularized Dirac function. The acceleration is stored on the cell
faces, and will be fed as a source term to the Navier-Stokes solver. */
vector forcing[];
event acceleration (i++) {

  /*We synchronize the eul field and make sure that it is updated before applying repulsive force */
  synchronize({Index_lagnode, Index_lag_id});

  /*We add the repulsive lubrication force for a better numerical stability*/
  # if LUBR_FORCE == 1  
  // lubrication_force(); 
  # endif

  face vector ae = a;
  foreach()
    if (cm[] > 1.e-20) foreach_dimension() forcing.x[] = 0.;
  for(int i=0; i<NCAPS; i++) {
    if (CAPS(i).isactive) lag2eul(forcing, &CAPS(i));
  }
  foreach_face()
    if (fm.x[] > 1.e-20) ae.x[] += .5*alpha.x[]*(forcing.x[] + forcing.x[-1]);
}

/** At the end of the simulation, we free the allocated memory.*/
event cleanup (t = end) {
  free_all_caps(&allCaps);
}


/**
## Additional functionalities
*/
#if dimension > 2
  #include "dump-ft.h"
  #include "post-processing-ft.h"
#endif


/**
## Tests
[advect_caps.c](../../tests/lagrangian_caps/advect_caps.c): Tests the
convergence of the advection scheme.


[curvature.c](../../tests/lagrangian_caps/curvature.c): Tests the computation of
the curvature at the Lagrangian nodes. Since the curvature depends on the
normals, this case also validates the computation of the normal vectors.
*/

