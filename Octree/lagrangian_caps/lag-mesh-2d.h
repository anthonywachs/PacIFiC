#ifndef NLP
  #define NLP 31
#endif
#ifndef RADIUS
  #define RADIUS .25
#endif

typedef struct edge {
  coord* vertex[2];
  double l0, l; // Initial edge length, current edge length.
} edge;

double compute_length(edge e) {
  /**
  TODO: deal with periodic boundary conditions
  */
  return sqrt(sq(e.vertex[0]->x - e.vertex[1]->x)
            + sq(e.vertex[0]->y - e.vertex[1]->y)
  #if dimension > 2
            + sq(e.vertex[0]->z - e.vertex[1]->z));
  #else
        );
  #endif
}

typedef struct lagMesh {
  int nlp;  // Number of Lagrangian vertices
  int nle;  // Number of Lagrangian edges
  coord* vertices;  // Array of vertices
  edge* edges;  // Array of edges
  coord* lagForces;  // Array of Lagrangian forces
  coord* lagVel; // Array of Lagrangian velocities
} lagMesh;

bool on_face(double p, int n, double l0) {
  if ((fabs(p/(l0/n)) - ((int)fabs(p/(l0/n)))) < 1.e-10) {
    return true;
  }
  else {
    return false;
  }
}

void correct_lag_pos(lagMesh* mesh) {
  for(int i=0; i < mesh->nlp; i++) {
    foreach_dimension() {
      if (on_face(mesh->vertices[i].x,N,L0)) mesh->vertices[i].x += 1.e-10;
      if(mesh->vertices[i].x > L0/2) mesh->vertices[i].x -= L0;
    }
  }
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
  p.mesh->vertices = malloc(nlp*sizeof(coord));
  p.mesh->edges = malloc((nlp-1)*sizeof(edge));
  p.mesh->lagForces = malloc(nlp*sizeof(coord));
  p.mesh->lagVel = malloc(nlp*sizeof(coord));


  double alpha = 2*pi/(nlp-1);
  for(int i=0; i<nlp; i++) {
    p.mesh->vertices[i].x = radius*cos(alpha*i);
    p.mesh->vertices[i].y = radius*sin(alpha*i);
  }
  for(int i=0; i<nlp-1; i++) {
    p.mesh->edges[i].vertex[0] = &(p.mesh->vertices[i]);
    p.mesh->edges[i].vertex[1] = (i+1<nlp) ? &(p.mesh->vertices[i+1]) :
      &(p.mesh->vertices[0]);
    p.mesh->edges[i].l0 = compute_length(p.mesh->edges[i]);
    p.mesh->edges[i].l = compute_length(p.mesh->edges[i]);
    foreach_dimension() {
      p.mesh->lagForces[i].x = 0.;
      p.mesh->lagVel[i].x = 0.;
    }
  }
  correct_lag_pos(p.mesh);
}


void advect_lagMesh(lagMesh* mesh) {
  for(int i=0; i < mesh->nlp; i++) {
    foreach_dimension() {
      mesh->vertices[i].x += dt*mesh->lagVel[i].x;
    }
  }
  correct_lag_pos(mesh);
}
