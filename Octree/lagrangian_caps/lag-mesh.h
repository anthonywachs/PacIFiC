/**
# Lagrangian mesh

In this file, we implement a Lagrangian mesh meant to track the position and
compute the stresses of an elasitc membrane.


## Structure of the mesh

In the Lagrangian mesh, each node is assigned coordinates, the IDs of its
two connecting edges (in 2D), an elastic force, a velocity, a normal vector,
a current and a reference curvature (i.e. the curvature of the unstressed
resting shape) and a stencil used for averaging the neighboring velocities and
spreading the Lagrangian force as a body force on the neighboring Eulerian
nodes. In the case of MPI, we also store the rank of the processor which owns
the Eulerian cell containing the Node.
*/
typedef struct lagNode {
  coord pos;
  coord lagForce;
  coord lagVel;
  coord normal;
  double curv;
  double ref_curv;
  Cache stencil;
  #if _MPI
    int pid;
  #endif
  #if dimension < 3
    int edge_ids[2];
  #else
    int nb_neighbors;
    int neighbor_ids[6];
    int nb_edges;
    int edge_ids[6];
    int nb_triangles;
    int triangle_ids[6];
  #endif
} lagNode;

/** Similarly, the edges of the mesh are assigned the IDs of the two nodes they
connect, an undeformed and a deformed lengths $l_0$ and $l$, and a normal
vector.  */
typedef struct Edge {
  int node_ids[2];
  #if dimension > 2
    int triangle_ids[2];
  #endif
  double l0, length; // Initial edge length, current stretch
  coord normal;
} Edge;

/** In 3 dimensions, we define triangle faces */
#if dimension > 2
  typedef struct Triangle {
    int node_ids[3];
    int edge_ids[3];
    double area, refArea;
    coord normal;
    coord centroid;
    coord refShape[2];
    double sfc[3][2]; // sfc for "shape function coefficients"
  } Triangle;
#endif

/** The lagMesh struct stores an array of nodes, edges as well as their
respective sizes. The three booleans are used to check if the stretches, normal
vectors and curvatures have been re-computed since the position of the mesh was
advected. */
typedef struct lagMesh {
  int nlp;  // Number of Lagrangian points
  int nle;  // Number of Lagrangian Edges
  lagNode* nodes;  // Array of nodes
  Edge* edges;  // Array of edges
  #if dimension > 2
    int nlt; // Number of Lagrangian triangles
    Triangle* triangles;
  #endif
  bool updated_stretches;
  bool updated_normals;
  bool updated_curvatures;
} lagMesh;

/** We denote by $NCAPS$ the number of Lagrangian meshes, or capsules, in the
simulation. It is one by default. */
#ifndef NCAPS
  #define NCAPS 1
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
  for(int i=0; i<caps->nbmb; i++) free_mesh(&(caps->mb[i]));
}

/** By default, the mesh is advected using a second-order two-step Runge Kutta
scheme. If the following macro is set to $0$, a first-order forward Euler schme
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
The function below computes the membrane stretch $\lambda = \frac{l}{l_0}$,
with $l$ the current length of an edge and $l_0$ its reference length
(in the untressed resting shape).
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
a Lagrangian mesh.
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
/** The function below assumes that the Lagrangian mesh contains the origin and
is convex, and swaps the order of the nodes in order to compute an outward
normal vector. This only need to be performed at the creation of the mesh since
the outward property of the normal vectors won't change through the simulation.
*/
void comp_initial_area_normals(lagMesh* mesh) {
  for(int i=0; i<mesh->nlt; i++) {
    int nid[3]; // node ids
    coord centroid;
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

/** The two function below compute the outward normal vector to all the
triangles of the mesh. */
void comp_triangle_area_normal(lagMesh* mesh, int i) {
  int nid[3]; // node ids
  for(int j=0; j<3; j++) nid[j] = mesh->triangles[i].node_ids[j];
  foreach_dimension() mesh->triangles[i].centroid.x = 0.;
  for(int j=0; j<3; j++) {
    nid[j] = mesh->triangles[i].node_ids[j];
    foreach_dimension()
      mesh->triangles[i].centroid.x += mesh->nodes[nid[j]].pos.x/3;
  }
  coord normal, e[2];
  for(int j=0; j<2; j++)
    foreach_dimension()
      e[j].x = mesh->nodes[nid[0]].pos.x - mesh->nodes[nid[j+1]].pos.x;
  foreach_dimension() normal.x = e[0].y*e[1].z - e[0].z*e[1].y;
  double norm = sqrt(sq(normal.x) + sq(normal.y) + sq(normal.z));
  foreach_dimension() mesh->triangles[i].normal.x = normal.x/norm;
  mesh->triangles[i].area = norm/2.;
}

void comp_triangle_area_normals(lagMesh* mesh) {
  for(int i=0; i<mesh->nlt; i++) comp_triangle_area_normal(mesh, i);
}
#endif

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
    #else // dimension = 3
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
    }
    #endif
    mesh->updated_normals = true;
  }
}


/**
The function below computes the signed curvature of the Lagrangian mesh at each
node. It scales in a second-order fashion with the number of Lagrangian points.
It only works in two dimensions.
*/
void comp_curvature(lagMesh* mesh) {
  if (!mesh->updated_curvatures) {
    comp_normals(mesh);
    #if dimension < 3
    bool up; // decide if we switch the x and y axes
    lagNode* cn; // current node
    for(int i=0; i<mesh->nlp; i++) {
      cn = &(mesh->nodes[i]);
      up = (fabs(cn->normal.y) > fabs(cn->normal.x)) ? true : false;
      coord p[5]; // store the coordinates of the current node and of its
                  // neighbors'
      for(int j=0; j<5; j++) {
        int index = (mesh->nlp + i - 2 + j)%mesh->nlp;
        foreach_dimension() p[j].x = up ? mesh->nodes[index].pos.x :
          mesh->nodes[index].pos.y;
      }
      /** If one of the neighboring nodes is across a periodic boundary, we
correct its position */
      for(int j=0; j<5; j++) {
        if (j!=2) {
          foreach_dimension() {
            if (ACROSS_PERIODIC(p[j].x,p[2].x)) {
              p[j].x += (ACROSS_PERIODIC(p[j].x + L0, p[2].x)) ? -L0 : L0;
            }
          }
        }
      }

/** Since the bending force will involve taking the laplacian of the curvature,
we seek a cuvrature to fourth order accuracy, so we need to interpolate the
membrane with a fourth-degree polynomial, which we need to differentiate twice
to get the curvature:
$$P_4(x) = \sum_{j=i-2}^{i+2} y_j \prod_{k \neq j} \frac{x - x_j}{x_k - x_j} $$

$$P'_4(x) = \sum_{j=i-2}^{i+2} y_j \left( \prod_{k \neq j} \frac{1}{x_k - x_j}
\right) \sum_{l \neq j}\prod_{m \neq j, m \neq l} x - x_m $$

$$P''_4(x) = \sum_{j=i-2}^{i+2} y_j \left( \prod_{k \neq j} \frac{1}{x_k - x_j}
\right) \sum_{l \neq j}\sum_{m \neq j, m \neq l}\sum_{n \neq j, n \neq l, n
\neq m}x - x_n $$
*/
      double dy = 0.;
      double ddy = 0.;
      for(int j=0; j<5; j++) {
        double b1 = 0.; double b2 = 0.;
        for(int l=0; l<5; l++) {
          if (l!=j) {
            double c1 = 1.; double c2 = 0.;
            for(int m=0; m<5; m++) {
              if (m!=j && m!=l) {
                double d2 = 1.;
                for(int n=0; n<5; n++) {
                  if (n!=j && n!=l && n!= m) {
                    d2 *= p[2].x - p[n].x;
                  }
                }
                c1 *= p[2].x - p[m].x;
                c2 += d2;
              }
            }
            b1 += c1;
            b2 += c2;
          }
        }
        for(int k=0; k<5; k++) {
          if (k!=j) {
            b1 /= (p[k].x - p[j].x);
            b2 /= (p[k].x - p[j].x);
          }
        }
        dy += b1*p[j].y;
        ddy += b2*p[j].y;
      }

      /** The formula for the signed curvature of a function y(x) is
  $$ \kappa = \frac{y''}{(1 + y'^2)^{\frac{3}{2}}}. $$
  The sign is dertemined from a parametrization of the curve: walking
  anticlockwise along the curve, if we turn left the curvature is positive. This
  statement can be easily written as a dot product between edge i's normal
  vector and edge (i+1)'s direction vector.*/
      coord a, b;
      foreach_dimension() {
        a.x = mesh->edges[mesh->nodes[i].edge_ids[0]].normal.x;
        b.x = mesh->edges[mesh->nodes[i].edge_ids[1]].normal.x;
      }
      int s = (a.x*b.x + a.y*b.y > 0) ? 1 : -1;
      cn->curv = s*fabs(ddy)/cube(sqrt(1 + sq(dy)));
    }
    #endif
    mesh->updated_curvatures = true;
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
    foreach_dimension() {
      if (on_face(mesh->nodes[i].pos.x, N, L0))
        mesh->nodes[i].pos.x += 1.e-10;
      if(fabs(mesh->nodes[i].pos.x) > L0/2)
        mesh->nodes[i].pos.x -= L0*sign(mesh->nodes[i].pos.x);
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
second-order Runge Kutta scheme is used. By setting the macro $ADVECT\_LAG\_RK2$
to zero, a simple forward Euler scheme is used as a scheme.
*/
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
      buffer_mesh.nodes[j].stencil.nm = STENCIL_SIZE*128;
      buffer_mesh.nodes[j].stencil.p = malloc(STENCIL_SIZE*sizeof(Index));
    }
    generate_lag_stencils(&buffer_mesh);
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
  generate_lag_stencils(mesh);
}


/**
## Putting the pieces together

Below with call the above functions at the appropriate time using the Basilisk
event syntax.

We start by creating empty Lagrangian meshes, and allocating an acceleration
field in case it isn't done yet by another Basilisk module.
*/
event defaults (i = 0) {
  mbs.nbmb = NCAPS;
  for(int i=0; i<mbs.nbmb; i++) {
    initialize_empty_mb(&mbs.mb[i]);
  }
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face()
      a.x[] = 0.;
    #if OLD_QCC
    boundary ((scalar *){a});
    #endif
  }
}

/** Before the iterations start, we allocate memory for the stencils and
generate them. Note that this implementation assumes the membrane was
initialized in the init event. */
event init (i = 0) {
  for(int i=0; i<mbs.nbmb; i++) {
    for(int j=0; j<MB(i).nlp; j++) {
      MB(i).nodes[j].stencil.n = STENCIL_SIZE;
      MB(i).nodes[j].stencil.nm = STENCIL_SIZE*128;
      MB(i).nodes[j].stencil.p = (Index*) malloc(STENCIL_SIZE*sizeof(Index));
    }
    generate_lag_stencils(&MB(i));
  }
}

/** Below, we advect each Lagrangian node using the interpolated Eulerian
velocities. We also use this loop as an opportunity to
re-initialize the Lagrangian forces to zero. */
event tracer_advection(i++) {
  for(int i=0; i<mbs.nbmb; i++) {
    advect_lagMesh(&mbs.mb[i]);
    for(int j=0; j<mbs.mb[i].nlp; j++)
      foreach_dimension() mbs.mb[i].nodes[j].lagForce.x = 0.;
  }
}

/** In the acceleration event, we transfer the Lagrangian forces to the fluid
using a regularized Dirac function. The acceleration is stored on the cell
faces, and will be fed as a source term to the Navier-Stokes solver. */
vector forcing[];
event acceleration (i++) {
  face vector ae = a;
  foreach() foreach_dimension() forcing.x[] = 0.;
  for(int i=0; i<mbs.nbmb; i++) lag2eul(forcing, &mbs.mb[i]);
  foreach_face() ae.x[] += .5*alpha.x[]*(forcing.x[] + forcing.x[-1]);
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
