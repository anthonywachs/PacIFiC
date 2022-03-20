#ifndef NLP
  #define NLP 31
#endif
#ifndef RADIUS
  #define RADIUS .25
#endif
#ifndef NCAPS
  #define NCAPS 1
#endif

/** In the Lagrangian mesh, each node is assigned coordinates, the IDs of its
two connecting edges (in 2D), an elastic force, and a velocity. */
typedef struct lagNode {
  coord pos;
  int edge_ids[2];
  coord lagForce;
  coord lagVel;
} lagNode;

/** Similarly, the edges of the mesh are assigned the IDs of the two nodes they
connect, an undeformed length $l_0$ and a stretch ratio $\lamba =
\frac{l}{l_0}$, where $l$ is the current length of the edge. */
typedef struct Edge {
  int vertex_ids[2];
  double l0, st; // Initial edge length, current stretch
} Edge;

/** The lagMesh struct stores an array of nodes, edges as well as their
respective sizes. */
typedef struct lagMesh {
  int nlp;  // Number of Lagrangian points
  int nle;  // Number of Lagrangian Edges
  lagNode* nodes;  // Array of nodes
  Edge* edges;  // Array of edges
} lagMesh;

/** The function below computes the length of an edge. It takes as arguments
a pointer to the mesh as well as the ID of the edge of interest. */
double comp_length(lagMesh* mesh, int i) {
  double length = 0.;
  int v1, v2;
  v1 = mesh->edges[i].vertex_ids[0];
  v2 = mesh->edges[i].vertex_ids[1];
  foreach_dimension() {
    /** This if-statement deals with the case of periodic boundaries. */
    if (fabs(mesh->nodes[v1].pos.x - mesh->nodes[v2].pos.x) > L0/2.) {
      length += (fabs(mesh->nodes[v1].pos.x - L0
        - mesh->nodes[v2].pos.x) > L0/2.) ?
        sq(mesh->nodes[v1].pos.x + L0 - mesh->nodes[v2].pos.x) :
        sq(mesh->nodes[v1].pos.x - L0 - mesh->nodes[v2].pos.x) ;
      }
    else length += sq(mesh->nodes[v1].pos.x - mesh->nodes[v2].pos.x);
  }
  return sqrt(length);
}

void comp_mb_stretch(lagMesh* mesh) {
  for(int i=0; i < mesh->nle; i++)
    mesh->edges[i].st = comp_length(mesh, i)/mesh->edges[i].l0;
}

/** If a Lagrangian node falls exactly on an edge or a vertex of the Eulerian
mesh, some issues arise when checking for periodic boundary conditions. As a
quick fix, if this is the case we shift the point position by $10^{-10}$, as
is done in the two functions below. */
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
}

void initialize_empty_mb(lagMesh* mesh) {
  mesh->nlp = 0;
  mesh->nle = 0;
  mesh->nodes = NULL;
  mesh->edges = NULL;
}

struct _initialize_circular_mb {
  lagMesh* mesh;
  double radius;
  double nlp;
};

void initialize_circular_mb(struct _initialize_circular_mb p) {
  double radius = (p.radius) ? p.radius : RADIUS;
  int nlp = (p.nlp) ? p.nlp : NLP;
  p.mesh->nlp = nlp;
  p.mesh->nle = nlp - 1;
  p.mesh->nodes = malloc(nlp*sizeof(lagNode));
  p.mesh->edges = malloc((nlp-1)*sizeof(Edge));


  double alpha = 2*pi/(nlp-1);
  /** Fill the array of nodes */
  for(int i=0; i<nlp; i++) {
    p.mesh->nodes[i].pos.x = radius*cos(alpha*i);
    p.mesh->nodes[i].pos.y = radius*sin(alpha*i);
    p.mesh->nodes[i].edge_ids[0] = -1;
    p.mesh->nodes[i].edge_ids[1] = -1;
    foreach_dimension() {
      p.mesh->nodes[i].lagForce.x = 0.;
      p.mesh->nodes[i].lagVel.x = 0.;
    }
  }
  /** Fill the array of edges.
  For the last edge, the next vertex id is 0 */
  for(int i=0; i<nlp-1; i++) {
    p.mesh->edges[i].vertex_ids[0] = i;
    if (p.mesh->nodes[i].edge_ids[0] < 0)
      p.mesh->nodes[i].edge_ids[0] = i;
    else
      p.mesh->nodes[i].edge_ids[1] = i;
    int next_vertex_id = (i+1<nlp) ? i+1 : 0;
    p.mesh->edges[i].vertex_ids[1] = next_vertex_id;
    if (p.mesh->nodes[next_vertex_id].edge_ids[0] < 0)
      p.mesh->nodes[next_vertex_id].edge_ids[0] = i;
    else
      p.mesh->nodes[next_vertex_id].edge_ids[1] = i;
    p.mesh->edges[i].l0 = comp_length(p.mesh, i);
    p.mesh->edges[i].st = 1.;
  }
  correct_lag_pos(p.mesh);
}

#include "reg-dirac.h"

/** The function below advects each Lagrangian node of a capsule mesh by
interpolating the velocities around the node of interest. A simple forward Euler
scheme is used as a scheme. */
void advect_lagMesh(lagMesh* mesh) {
  eul2lag(mesh);
  for(int i=0; i < mesh->nlp; i++) {
    foreach_dimension() {
      mesh->nodes[i].pos.x += dt*mesh->nodes[i].lagVel.x;
    }
  }
  correct_lag_pos(mesh);
}

typedef struct Capsules {
  lagMesh mb[NCAPS];
  int nbmb;
} Capsules;


Capsules mbs;
event defaults (i = 0) {
  mbs.nbmb = NCAPS;
  for(int i=0; i<mbs.nbmb; i++) initialize_empty_mb(&mbs.mb[i]);
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face()
      a.x[] = 0.;
    boundary ((scalar *){a});
  }
}

/** Below, we advect each Lagrangian node using the interpolated Eulerian
velocities. We also take this loop onto the nodes as an opportunity to
re-initialize the Lagrangian forces to zero. */
event tracer_advection(i++) {
  for(int i=0; i<mbs.nbmb; i++) {
    advect_lagMesh(&mbs.mb[i]);
    for(int j=0; j<mbs.mb->nlp; j++)
      foreach_dimension() mbs.mb[i].nodes[j].lagForce.x = 0.;
  }
}

/** In the acceleration event, we transfer the Lagrangian forces to the fluid
using a regularized Dirac function. The acceleration is stored on the cell
faces, and will fed as a source term to the Navier-Stokes solver. */
event acceleration (i++) {
  face vector ae = a;
  for(int i=0; i<mbs.nbmb; i++) lag2eul(ae,&mbs.mb[i]);
}
