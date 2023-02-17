/**
# Lagrangian mesh

In this file, we implement a Lagrangian mesh for the front-tracking method,
meant to track the position and compute the stresses of an elasitc membrane.


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
  } Triangle;
#endif

/** The ```lagMesh``` structure contains arrays of the previously introduced nodes, edges and triangles. It defines an unstructured mesh, the membrane of our capsule. Its attributes are:

*  ```nlp``` the number of Lagrangian nodes (or Lagrangian points, hence the "p")
* ```nodes```, the array of Lagrangian nodes
* ```nle```, the number of Lagrangian edges
* ```edges```, the array of Lagrangian edges
* In case of 3D simulations:
    * ```nlt```, the number of Lagrangian triangles
    * ```triangles```, the array of Lagrangian triangles
* ```updated_stretches```, a boolean used to check if the current length of the edges has been updated since the last advection of the Lagrangian nodes
* ```updated_normals```, a similar boolean telling if the nodal normal vectors should be recomputed
* ```updated_curvatures```, a last boolean telling if the nodal curvatures should be recomputed.
* ```active```, a boolean indicating if the capsule exists in the flow (useful
when capsules are introduced during a simulation)
*/

typedef struct lagMesh {
  int nlp;
  lagNode* nodes;
  int nle;
  Edge* edges;
  #if dimension > 2
    int nlt;
    Triangle* triangles;
  #endif
  coord centroid;
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
them are considered). The macro $MB(k)$ can be used as a shortcut to access the
$k^{th}$ membrane. */
typedef struct Capsules {
  lagMesh mb[NCAPS];
  int nbmb;
} Capsules;
Capsules mbs;
#define MB(i) (mbs.mb[i])


/**
## Initialization, memory management and useful macros.
*/
void initialize_empty_mb(lagMesh* mesh) {
  mesh->nlp = 0;
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

void free_mesh(lagMesh* mesh) {
  for(int i=0; i<mesh->nlp; i++) free(mesh->nodes[i].stencil.p);
  free(mesh->nodes);
  free(mesh->edges);
  #if dimension > 2
    free(mesh->triangles);
  #endif
}

void free_caps(Capsules* caps) {
  for(int i=0; i<caps->nbmb; i++)
    if (mbs.mb[i].isactive)
      free_mesh(&(caps->mb[i]));
}

/** By default, the mesh is advected using a second-order two-step Runge Kutta
scheme. If the following macro is set to 0, a first-order forward Euler schme
is used instead. */
#ifndef ADVECT_LAG_RK2
  #define ADVECT_LAG_RK2 1
#endif

/** We specify the size of the 5x5(x5) stencil in 2D or 3D. */
#if dimension < 3
  #define STENCIL_SIZE 25
#else
  #define STENCIL_SIZE 125
#endif

/** The next few macros are useful to compute signed distances and averages
across periodic boundaries. We assume for this purpose that the length of
the edges are less that half the domain size, which in practice should always
be the case. */
#define ACROSS_PERIODIC(a,b) (fabs(a - b) > L0/2.)
#define PERIODIC_1DIST(a,b) (fabs(a - L0 - b) > L0/2. ? a + L0 - b : a - L0 - b)
#define GENERAL_1DIST(a,b) (ACROSS_PERIODIC(a,b) ? PERIODIC_1DIST(a,b) : a - b)
#define PERIODIC_1DAVG(a,b) (fabs(a - L0 - b) > L0/2. ? a + L0 + b : a - L0 + b)
#define GENERAL_1DAVG(a,b) (ACROSS_PERIODIC(a,b) ? PERIODIC_1DAVG(a,b) : a + b)

#if dimension < 3
  #define cnorm(a) (sqrt(sq(a.x) + sq(a.y)))
  #define cdot(a,b) (a.x*b.x + a.y*b.y)
#else
  #define cnorm(a) (sqrt(sq(a.x) + sq(a.y) + sq(a.z)))
  #define cdot(a,b) (a.x*b.x + a.y*b.y + a.z*b.z)
#endif



/**
## Geometric computations

The function below computes the length of an edge. It takes as arguments
a pointer to the mesh as well as the ID of the edge of interest.
*/
double edge_length(lagMesh* mesh, int i) {
  double length = 0.;
  int v1, v2;
  v1 = mesh->edges[i].node_ids[0];
  v2 = mesh->edges[i].node_ids[1];
  foreach_dimension() {
    length += sq(GENERAL_1DIST(mesh->nodes[v1].pos.x, mesh->nodes[v2].pos.x));
  }
  return sqrt(length);
}

/**
The function ```compute_lengths``` below computes the lengths of all edges. It
takes as an argument a pointer to the mesh. If the optional argument
```force``` is set to ```true```, the edges' lengths are computed no matter the value of ```updated_stretches```.
*/
struct _compute_lengths{
  lagMesh* mesh;
  bool force;
};

void compute_lengths(struct _compute_lengths p) {
  lagMesh* mesh = p.mesh;
  bool force = (p.force) ? p.force : false;
  if (force || !mesh->updated_stretches) {
    for(int i=0; i < mesh->nle; i++)
      mesh->edges[i].length = edge_length(mesh, i);
    mesh->updated_stretches = true;
  }
}

#if dimension < 3
/**
The two functions below compute the outward normal vector to all the edges of
a Lagrangian mesh, for 2D simulations.
*/
void comp_edge_normal(lagMesh* mesh, int i) {
    int node_id[2];
    for(int j=0; j<2; j++) node_id[j] = mesh->edges[i].node_ids[j];
    mesh->edges[i].normal.y = GENERAL_1DIST(mesh->nodes[node_id[0]].pos.x,
      mesh->nodes[node_id[1]].pos.x);
    mesh->edges[i].normal.x = GENERAL_1DIST(mesh->nodes[node_id[1]].pos.y,
      mesh->nodes[node_id[0]].pos.y);
    double normn = sqrt(sq(mesh->edges[i].normal.x)
      + sq(mesh->edges[i].normal.y));
    foreach_dimension() mesh->edges[i].normal.x /= normn;
}

void comp_edge_normals(lagMesh* mesh) {
  for(int i=0; i<mesh->nle; i++) comp_edge_normal(mesh, i);
}
#else // dimension > 2
/** In 3D simulations, the function below assumes that the Lagrangian mesh contains the origin and
is convex, and swaps the order of the nodes in order to compute an outward
normal vector. This only need to be performed at the creation of the mesh since
the outward property of the normal vectors won't change through the simulation.
*/
void comp_initial_area_normals(lagMesh* mesh) {
  for(int i=0; i<mesh->nlt; i++) {
    int nid[3]; // node ids
    coord centroid; /** Note: the centroid is only valid if the triangle is not
    across periodic boundaries, which is fine for this function since it is
    assumed the center of the membrane is at the origin. */
    foreach_dimension() centroid.x = 0.;
    for(int j=0; j<3; j++) {
      nid[j] = mesh->triangles[i].node_ids[j];
      foreach_dimension() centroid.x += mesh->nodes[nid[j]].pos.x/3;
    }
    coord normal, e[2];
    for(int j=0; j<2; j++)
      foreach_dimension()
        e[j].x = GENERAL_1DIST(mesh->nodes[nid[0]].pos.x,
          mesh->nodes[nid[j+1]].pos.x);
    foreach_dimension() normal.x = e[0].y*e[1].z - e[0].z*e[1].y;
    double norm = sqrt(sq(normal.x) + sq(normal.y) + sq(normal.z));
    double dp = 0.; // dp for "dot product"
    foreach_dimension() dp += normal.x*centroid.x;
    /** If the dot product is negative, the computed normal is inward and we
    need to swap two nodes of the triangle.*/
    if (dp < 0) {
      mesh->triangles[i].node_ids[1] = nid[2];
      mesh->triangles[i].node_ids[2] = nid[1];
      foreach_dimension() normal.x *= -1;
    }
    foreach_dimension() {
      mesh->triangles[i].centroid.x = centroid.x;
      mesh->triangles[i].normal.x = normal.x/norm;
    }
    mesh->triangles[i].area = norm/2;
  }
}

/** The two functions below compute the outward normal vector to all the
triangles of the mesh, for 3D simulations. */
void comp_triangle_area_normal(lagMesh* mesh, int i) {
  int nid[3]; // node ids
  for(int j=0; j<3; j++) nid[j] = mesh->triangles[i].node_ids[j];
  /** The next 15 lines compute the centroid of the triangle, making sure it is
  valid when the triangle lies across periodic boundaries. */
  foreach_dimension() mesh->triangles[i].centroid.x = 0.;
  for(int j=0; j<3; j++) {
    foreach_dimension() {
      mesh->triangles[i].centroid.x +=
        ACROSS_PERIODIC(mesh->nodes[nid[j]].pos.x/3,
        mesh->nodes[nid[0]].pos.x/3) ? mesh->nodes[nid[j]].pos.x/3 - L0 :
        mesh->nodes[nid[j]].pos.x/3;
    }
  }
  coord origin = {X0 + L0/2, Y0 + L0/2, Z0 + L0/2};
  foreach_dimension() {
    if (fabs(mesh->triangles[i].centroid.x - origin.x) > L0/2.) {
      if (mesh->triangles[i].centroid.x - origin.x > 0)
        mesh->triangles[i].centroid.x -= L0;
      else mesh->triangles[i].centroid.x += L0;
    }
  }
  coord normal, e[2];
  for(int j=0; j<2; j++)
    foreach_dimension()
      e[j].x = GENERAL_1DIST(mesh->nodes[nid[0]].pos.x,
        mesh->nodes[nid[j+1]].pos.x);
  foreach_dimension() normal.x = e[0].y*e[1].z - e[0].z*e[1].y;
  double norm = sqrt(sq(normal.x) + sq(normal.y) + sq(normal.z));
  foreach_dimension() mesh->triangles[i].normal.x = normal.x/norm;
  mesh->triangles[i].area = norm/2.;
}

void comp_triangle_area_normals(lagMesh* mesh) {
  for(int i=0; i<mesh->nlt; i++) comp_triangle_area_normal(mesh, i);
}
#endif

/**
The function below computes the centroid of the capsule as the average of the
coordinates of all its nodes. The centroid is stored as an attribute of the
caps structure.
*/
void comp_centroid(lagMesh* mesh) {
  foreach_dimension() mesh->centroid.x = 0.;
  for(int i=0; i<mesh->nlp; i++)
    foreach_dimension() mesh->centroid.x += mesh->nodes[i].pos.x;
  foreach_dimension() mesh->centroid.x /= mesh->nlp;
}

/** The function below updates the normal vectors on all the nodes as well as
the lengths and midpoints of all the edges (in 2D) or the area and centroids of
all the triangles (in 3D). */
void comp_normals(lagMesh* mesh) {
  if (!mesh->updated_normals) {
    #if dimension < 3
    compute_lengths(mesh);
    for(int i=0; i<mesh->nlp; i++) {
      coord n[2];
      double l[2];
      double normn;
      for(int j=0; j<2; j++) {
        int edge_id;
        edge_id = mesh->nodes[i].edge_ids[j];
        l[j] = mesh->edges[edge_id].length;
        comp_edge_normal(mesh, edge_id);
        foreach_dimension() n[j].x = mesh->edges[edge_id].normal.x;
      }
      /** the normal vector at a node is the weighted average of the normal
      vectors of its edges. The average is weighted by the distance of the node
      to each of the edges' centers. */
      double epsilon = l[1]/(l[0] + l[1]);
      normn = 0.;
      foreach_dimension() {
        mesh->nodes[i].normal.x = epsilon*n[0].x + (1. - epsilon)*n[1].x;
        normn += sq(mesh->nodes[i].normal.x);
      }
      normn = sqrt(normn);
      foreach_dimension() mesh->nodes[i].normal.x /= normn;
    }
    #else // dimension == 3
    comp_triangle_area_normals(mesh);
    for(int i=0; i<mesh->nlp; i++) {
      foreach_dimension() mesh->nodes[i].normal.x = 0.;
      double sw = 0.; // sum of the weights
      for(int j=0; j<mesh->nodes[i].nb_triangles; j++) {
        int tid = mesh->nodes[i].triangle_ids[j];
        double dist = 0.;
        foreach_dimension()
          dist += sq(mesh->nodes[i].pos.x - mesh->triangles[tid].centroid.x);
        dist = sqrt(dist);
        sw += dist;
        foreach_dimension()
          mesh->nodes[i].normal.x += mesh->triangles[tid].normal.x*dist;
      }
      foreach_dimension()
        mesh->nodes[i].normal.x /= sw;
      double normn = cnorm(mesh->nodes[i].normal);
      foreach_dimension() mesh->nodes[i].normal.x /= normn;
    }
    #endif
    mesh->updated_normals = true;
  }
}

/**
If a Lagrangian node falls exactly on an edge or a vertex of the Eulerian
mesh, some issues arise when checking for periodic boundary conditions. As a
quick fix, if this is the case we shift the point position by $10^{-10}$, as
is done in the two functions below.
*/
bool on_face(double p, int n, double l0) {
  if ((fabs(p/(l0/n)) - ((int)fabs(p/(l0/n)))) < 1.e-10) return true;
  else return false;
}

void correct_lag_pos(lagMesh* mesh) {
  for(int i=0; i < mesh->nlp; i++) {
    coord origin = {X0 + L0/2, Y0 + L0/2, Z0 + L0/2};
    foreach_dimension() {
      if (on_face(mesh->nodes[i].pos.x, N, L0))
        mesh->nodes[i].pos.x += 1.e-10;
      if (mesh->nodes[i].pos.x > origin.x + L0/2)
        mesh->nodes[i].pos.x -= L0;
      else if (mesh->nodes[i].pos.x < origin.x - L0/2)
        mesh->nodes[i].pos.x += L0;
    }
  }
  mesh->updated_stretches = false;
  mesh->updated_normals = false;
  mesh->updated_curvatures = false;
}


/**
## Advection of the mesh

The advection have to be compatible with MPI, so we start by including the MPI
communication functions. We also include [reg-dirac.h](reg-dirac.h), which
implements the interpolation of the Eulerian velocities onto the Lagrangian
nodes.
*/

#if _MPI
  #include "lag-mesh-mpi.h"
#endif
#include "reg-dirac.h"

/**
The function below advects each Lagrangian node by
interpolating the velocities around the node of interest. By default, a
second-order Runge Kutta scheme is used. By setting the macro
```ADVECT_LAG_RK2``` to 0, a simple forward Euler scheme is used.
*/
trace
void advect_lagMesh(lagMesh* mesh) {
  eul2lag(mesh);
  #if !(ADVECT_LAG_RK2)
    for(int i=0; i < mesh->nlp; i++) {
      foreach_dimension() {
        mesh->nodes[i].pos.x += dt*mesh->nodes[i].lagVel.x;
      }
    }
  #else
    lagMesh buffer_mesh;
    buffer_mesh.isactive = true;
    buffer_mesh.nlp = mesh->nlp;
    buffer_mesh.nodes = malloc(mesh->nlp*sizeof(lagNode));
    for(int i=0; i<mesh->nlp; i++) {
      // Step 1 of RK2
      foreach_dimension()
        buffer_mesh.nodes[i].pos.x = mesh->nodes[i].pos.x +
          .5*dt*mesh->nodes[i].lagVel.x;
    }
    correct_lag_pos(&buffer_mesh);
    for(int j=0; j<buffer_mesh.nlp; j++) {
      buffer_mesh.nodes[j].stencil.n = STENCIL_SIZE;
      buffer_mesh.nodes[j].stencil.nm = STENCIL_SIZE;
      buffer_mesh.nodes[j].stencil.p = malloc(STENCIL_SIZE*sizeof(Index));
    }
    generate_lag_stencils_one_caps(&buffer_mesh);
    eul2lag(&buffer_mesh);
    for(int i=0; i<mesh->nlp; i++) {
      // Step 2 of RK2
      foreach_dimension()
        mesh->nodes[i].pos.x += dt*buffer_mesh.nodes[i].lagVel.x;
    }
    for(int i=0; i<buffer_mesh.nlp; i++) free(buffer_mesh.nodes[i].stencil.p);
    free(buffer_mesh.nodes);
  #endif
  correct_lag_pos(mesh);
  generate_lag_stencils_one_caps(mesh);
  comp_centroid(mesh);
}


/**
## Putting the pieces together

Below, we call the above functions at the appropriate time using the Basilisk
event syntax.

We start by creating empty Lagrangian meshes, and allocating an acceleration
field in case it isn't done yet by another Basilisk solver.
*/
event defaults (i = 0) {
  mbs.nbmb = NCAPS;
  #if (RESTART_CASE == 0)
    for(int i=0; i<mbs.nbmb; i++) {
      initialize_empty_mb(&mbs.mb[i]);
    }
  #endif
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face() a.x[] = 0.;
  }
}

/** Below, we advect each Lagrangian node using the interpolated Eulerian
velocities. We also use this loop as an opportunity to
re-initialize the Lagrangian forces to zero. */
event tracer_advection(i++) {
  for(int i=0; i<mbs.nbmb; i++) {
    if (mbs.mb[i].isactive) {
      advect_lagMesh(&mbs.mb[i]);
      for(int j=0; j<mbs.mb[i].nlp; j++)
        foreach_dimension() mbs.mb[i].nodes[j].lagForce.x = 0.;
    }
  }
}

/** In the acceleration event, we transfer the Lagrangian forces to the fluid
using a regularized Dirac function. The acceleration is stored on the cell
faces, and will be fed as a source term to the Navier-Stokes solver. */
vector forcing[];
event acceleration (i++) {
  face vector ae = a;
  foreach()
    if (cm[] > 1.e-20) foreach_dimension() forcing.x[] = 0.;
  for(int i=0; i<mbs.nbmb; i++) {
    if (mbs.mb[i].isactive) lag2eul(forcing, &mbs.mb[i]);
  }
  foreach_face()
    if (fm.x[] > 1.e-20) ae.x[] += .5*alpha.x[]*(forcing.x[] + forcing.x[-1]);
}

/** At the end of the simulation, we free the allocated memory.*/
event cleanup (t = end) {
  free_caps(&mbs);
}

#if dimension > 2
  #include "mesh-toolbox.h"
#endif

/**
## Tests
[advect_caps.c](../../tests/lagrangian_caps/advect_caps.c): Tests the
convergence of the advection scheme.


[curvature.c](../../tests/lagrangian_caps/curvature.c): Tests the computation of
the curvature at the Lagrangian nodes. Since the curvature depends on the
normals, this case also validates the computation of the normal vectors.
*/
