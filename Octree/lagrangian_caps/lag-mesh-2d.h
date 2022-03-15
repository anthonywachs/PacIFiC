#ifndef NLP
  #define NLP 31
#endif
#ifndef RADIUS
  #define RADIUS .25
#endif

typedef struct edge {
  coord vertex[2];
  double l0, l; // Initial edge length, current edge length.
} edge;

double compute_length(edge e) {
  return sqrt(sq(e.vertex[0].x - e.vertex[1].x)
            + sq(e.vertex[0].y - e.vertex[1].y)
  #if dimension > 2
            + sq(e.vertex[0].z - e.vertex[1].z));
  #else
        );
  #endif
}

typedef struct lagMesh {
  int nlp;  // Number of Lagrangian vertices
  int nle;  // Number of Lagrangian edges
  coord* vertices;  // Array of vertices
  edge* edges;  // Array of edges
} lagMesh;

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

  double alpha = 2*pi/(nlp-1);
  for(int i=0; i<nlp; i++) {
    p.mesh->vertices[i].x = radius*cos(alpha*i);
    p.mesh->vertices[i].y = radius*sin(alpha*i);
  }
  for(int i=0; i<nlp-1; i++) {
    p.mesh->edges[i].vertex[0] = p.mesh->vertices[i];
    p.mesh->edges[i].vertex[1] = (i+1<nlp) ? p.mesh->vertices[i+1] : p.mesh->vertices[0];
    p.mesh->edges[i].l0 = compute_length(p.mesh->edges[i]);
    p.mesh->edges[i].l = compute_length(p.mesh->edges[i]);
  }
}
