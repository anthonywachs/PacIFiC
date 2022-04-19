#ifndef NLP
  #define NLP 31
#endif
#ifndef RADIUS
  #define RADIUS .25
#endif
#ifndef NCAPS
  #define NCAPS 1
#endif
#ifndef ADVECT_LAG_RK2
  #define ADVECT_LAG_RK2 1
#endif
#define MB(i) (mbs.mb[i])

/** In the Lagrangian mesh, each node is assigned coordinates, the IDs of its
two connecting edges (in 2D), an elastic force, and a velocity. */
typedef struct lagNode {
  coord pos;
  int edge_ids[2];
  coord lagForce;
  coord lagVel;
  coord normal;
  double curv;
  double ref_curv;
  double shear_tension;
} lagNode;

/** Similarly, the edges of the mesh are assigned the IDs of the two nodes they
connect, an undeformed length $l_0$ and a stretch ratio $\lamba =
\frac{l}{l_0}$, where $l$ is the current length of the edge. */
typedef struct Edge {
  int vertex_ids[2];
  double l0, st; // Initial edge length, current stretch
  coord normal;
} Edge;

/** The lagMesh struct stores an array of nodes, edges as well as their
respective sizes. */
typedef struct lagMesh {
  int nlp;  // Number of Lagrangian points
  int nle;  // Number of Lagrangian Edges
  lagNode* nodes;  // Array of nodes
  Edge* edges;  // Array of edges
  bool updated_stretches;
  bool updated_normals;
  bool updated_curvatures;
} lagMesh;

void free_mesh(lagMesh* mesh) {
  free(mesh->nodes);
  free(mesh->edges);
}

#define ACROSS_PERIODIC(a,b) (fabs(a - b) > L0/2.)
#define PERIODIC_1DIST(a,b) (fabs(a - L0 - b) > L0/2. ? a + L0 - b : a - L0 - b)
#define GENERAL_1DIST(a,b) (ACROSS_PERIODIC(a,b) ? PERIODIC_1DIST(a,b) : a - b)
#define PERIODIC_1DAVG(a,b) (fabs(a - L0 - b) > L0/2. ? a + L0 + b : a - L0 + b)
#define GENERAL_1DAVG(a,b) (ACROSS_PERIODIC(a,b) ? PERIODIC_1DAVG(a,b) : a + b)

/** The function below computes the length of an edge. It takes as arguments
a pointer to the mesh as well as the ID of the edge of interest. */
double comp_length(lagMesh* mesh, int i) {
  double length = 0.;
  int v1, v2;
  v1 = mesh->edges[i].vertex_ids[0];
  v2 = mesh->edges[i].vertex_ids[1];
  foreach_dimension() {
    length += sq(GENERAL_1DIST(mesh->nodes[v1].pos.x, mesh->nodes[v2].pos.x));
  }
  return sqrt(length);
}

void comp_mb_stretch(lagMesh* mesh) {
  if (!mesh->updated_stretches) {
    for(int i=0; i < mesh->nle; i++)
      mesh->edges[i].st = comp_length(mesh, i)/mesh->edges[i].l0;
    mesh->updated_stretches = true;
  }
}

void comp_edge_normal(lagMesh* mesh, int i) {
    int node_id[2];
    for(int j=0; j<2; j++) node_id[j] = mesh->edges[i].vertex_ids[j];
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

/** The function below updates the normal vectors on all the nodes as well as
on the midpoints of all the edges. */
void comp_normals(lagMesh* mesh) {
  if (!mesh->updated_normals) {
    comp_mb_stretch(mesh);
    for(int i=0; i<mesh->nlp; i++) {
      coord n[2];
      double l[2];
      double normn;
      for(int j=0; j<2; j++) {
        int edge_id;
        edge_id = mesh->nodes[i].edge_ids[j];
        l[j] = mesh->edges[edge_id].st*mesh->edges[edge_id].l0;
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
    mesh->updated_normals = true;
  }
}

void comp_curvature(lagMesh* mesh) {
  if (!mesh->updated_curvatures) {
    comp_normals(mesh);
    bool up; // decide if we switch the x and y axes
    lagNode* cn; // current node
    for(int i=0; i<mesh->nlp; i++) {
      cn = &(mesh->nodes[i]);
      up = (fabs(cn->normal.y) > fabs(cn->normal.x)) ? true : false;
      coord p[3]; // store the coordinates of the current node and of its
                  // neighbors'
      for(int j=0; j<3; j++) {
        int index = (i==0 && j==0) ? mesh->nlp-1 : (i-1+j)%mesh->nlp;
        foreach_dimension() p[j].x = up ? mesh->nodes[index].pos.x :
          mesh->nodes[index].pos.y;
      }
      /* If one of the neighboring nodes is across a periodic boundary, we
correct its position */
      for(int j=0; j<3; j+=2) {
        foreach_dimension() {
          if (ACROSS_PERIODIC(p[j].x,p[1].x)) {
            p[j].x += (ACROSS_PERIODIC(p[j].x + L0, p[1].x)) ? -L0 : L0;
          }
        }
      }

      double dy, ddy; dy = 0.; ddy = 0.;
      for(int j=0; j<3; j++) {
        dy += p[j].y*(2*p[1].x - p[(j+1)%3].x - p[(j+2)%3].x)/
          ((p[j].x - p[(j+1)%3].x)*(p[j].x - p[(j+2)%3].x));
        ddy += 2*p[j].y/((p[j].x - p[(j+1)%3].x)*(p[j].x - p[(j+2)%3].x));
      }
      /** The formula for the signed curvature of a function y(x) is
  $$ \kappa = \frac{y''}{(1 + y'^2)^{\frac{3}{2}}} $$.
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
    mesh->updated_curvatures = true;
  }
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
  mesh->updated_stretches = false;
  mesh->updated_normals = false;
  mesh->updated_curvatures = false;
}

void initialize_empty_mb(lagMesh* mesh) {
  mesh->nlp = 0;
  mesh->nle = 0;
  mesh->nodes = NULL;
  mesh->edges = NULL;
  mesh->updated_stretches = false;
  mesh->updated_normals = false;
  mesh->updated_curvatures = false;
}

#include "reg-dirac.h"

/** The function below advects each Lagrangian node of a capsule mesh by
interpolating the velocities around the node of interest. A simple forward Euler
scheme is used as a scheme. */
void advect_lagMesh(lagMesh* mesh) {
  eul2lag(mesh);
  /** By default the Lagrangian mesh is advected with a second order Runge-Kutta
  scheme. Otherwise, we use a first order forward Euler advection scheme. */
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
    eul2lag(&buffer_mesh);
    for(int i=0; i<mesh->nlp; i++) {
      // Step 2 or RK2
      foreach_dimension()
        mesh->nodes[i].pos.x += dt*buffer_mesh.nodes[i].lagVel.x;
    }
    free(buffer_mesh.nodes);
  #endif
  correct_lag_pos(mesh);
}

typedef struct Capsules {
  lagMesh mb[NCAPS];
  int nbmb;
} Capsules;

void free_caps(Capsules* caps) {
  for(int i=0; i<caps->nbmb; i++) free_mesh(&(caps->mb[i]));
}

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
vector forcing[];
event acceleration (i++) {
  face vector ae = a;
  foreach() foreach_dimension() forcing.x[] = 0.;
  for(int i=0; i<mbs.nbmb; i++) lag2eul(forcing, &mbs.mb[i]);
  foreach_face() ae.x[] += .5*alpha.x[]*(forcing.x[] + forcing.x[-1]);
}

event cleanup (t = end) {
  free_caps(&mbs);
}
