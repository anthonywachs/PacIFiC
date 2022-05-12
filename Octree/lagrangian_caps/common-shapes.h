/**
Some common initial shapes for vesicles
*/

#ifndef NLP
  #define NLP 100
#endif
#ifndef RADIUS
  #define RADIUS 1.
#endif

struct _initialize_circular_mb {
  lagMesh* mesh;
  double radius;
  double nlp;
  double inclination;
  coord shift;
};

void initialize_circular_mb(struct _initialize_circular_mb p) {
  double radius = (p.radius) ? p.radius : RADIUS;
  int nlp = (p.nlp) ? p.nlp : NLP;
  coord shift;
  if (p.shift.x || p.shift.y || p.shift.z)
    {shift.x = p.shift.x; shift.y = p.shift.y; shift.z = p.shift.z;}
  else {shift.x = 0.; shift.y = 0.; shift.z = 0.;}
  p.mesh->nlp = nlp;
  p.mesh->nle = nlp;
  p.mesh->nodes = malloc(nlp*sizeof(lagNode));
  p.mesh->edges = malloc(nlp*sizeof(Edge));

  double alpha = 2*pi/(nlp);
  /** Fill the array of nodes */
  for(int i=0; i<nlp; i++) {
    p.mesh->nodes[i].pos.x = radius*cos(alpha*i) + shift.x;
    p.mesh->nodes[i].pos.y = radius*sin(alpha*i) + shift.y;
    p.mesh->nodes[i].edge_ids[0] = -1;
    p.mesh->nodes[i].edge_ids[1] = -1;
    foreach_dimension() {
      p.mesh->nodes[i].lagForce.x = 0.;
      p.mesh->nodes[i].lagVel.x = 0.;
    }
  }
  /** Fill the array of edges.
  For the last edge, the next vertex id is 0 */
  for(int i=0; i<p.mesh->nle; i++) {
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
  }
  /** The above procedure switches the two egde ids for the first node, which
  we correct below */
  p.mesh->nodes[0].edge_ids[0] = p.mesh->nle-1;
  p.mesh->nodes[0].edge_ids[1] = 0;
  correct_lag_pos(p.mesh);
  for(int i=0.; i<p.mesh->nle; i++) {
    p.mesh->edges[i].l0 = comp_length(p.mesh, i);
    p.mesh->edges[i].length = p.mesh->edges[i].l0;
  }

  #ifdef CAPS_VISCOSITY
    fraction(prevI, sq(radius) - sq(x - shift.x) - sq(y - shift.y));
  #endif
}

void initialize_biconcave_mb(struct _initialize_circular_mb p) {
  double radius = (p.radius) ? p.radius : RADIUS;
  int nlp = (p.nlp) ? p.nlp : NLP;
  double inclination = (p.inclination) ? p.inclination : 0.;
  coord shift;
  if (p.shift.x || p.shift.y || p.shift.z)
    {shift.x = p.shift.x; shift.y = p.shift.y; shift.z = p.shift.z;}
  else {shift.x = 0.; shift.y = 0.; shift.z = 0.;}
  double c = 1.3858189;
  p.mesh->nlp = nlp;
  p.mesh->nle = nlp;
  p.mesh->nodes = malloc(nlp*sizeof(lagNode));
  p.mesh->edges = malloc(nlp*sizeof(Edge));

  double alpha = 2*pi/(nlp);
  /** Fill the array of nodes */
  for(int i=0; i<nlp; i++) {
    double nrx = radius*c*cos(alpha*i);
    double nry = .5*radius*c*sin(alpha*i)*(0.207 +
      2.003*sq(cos(2*pi-alpha*i)) - 1.123*sq(sq(cos(2*pi-alpha*i))));
    p.mesh->nodes[i].pos.x = nrx*cos(inclination) - nry*sin(inclination) +
      shift.x;
    p.mesh->nodes[i].pos.y = nrx*sin(inclination) + nry*cos(inclination) +
      shift.y;
    p.mesh->nodes[i].edge_ids[0] = -1;
    p.mesh->nodes[i].edge_ids[1] = -1;
    foreach_dimension() {
      p.mesh->nodes[i].lagForce.x = 0.;
      p.mesh->nodes[i].lagVel.x = 0.;
    }
  }
  /** Fill the array of edges.
  For the last edge, the next vertex id is 0 */
  for(int i=0; i<p.mesh->nle; i++) {
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
  }
  /** The above procedure switches the two egde ids for the first node, which
  we correct below */
  p.mesh->nodes[0].edge_ids[0] = p.mesh->nle-1;
  p.mesh->nodes[0].edge_ids[1] = 0;
  correct_lag_pos(p.mesh);
  for(int i=0.; i<p.mesh->nle; i++) {
    p.mesh->edges[i].l0 = comp_length(p.mesh, i);
    p.mesh->edges[i].length = p.mesh->edges[i].l0;
  }

  #ifdef CAPS_VISCOSITY
    // We define below the local coordinates of the RBC and the parametric angle
    #define MY_X ((x - shift.x)*cos(inclination) + \
      (y - shift.y)*sin(inclination))
    #define MY_Y (-(x - shift.x)*sin(inclination) + \
      (y - shift.y)*cos(inclination))
    #define MY_Z z
    #define COSPHI2 ((sq(MY_X)+sq(MY_Z))/sq(radius*c))
    #define LAMBDA (0.207 + 2.003*COSPHI2 - 1.123*sq(COSPHI2))
    fraction(prevI, 1. - sq(MY_X/(radius*c)) -
      sq(2*MY_Y/(LAMBDA*radius*c)) -
      sq(MY_Z/(radius*c)));
  #endif
}


struct _initialize_elliptic_mb {
  lagMesh* mesh;
  double a;
  double b;
  double nlp;
  double inclination;
};

void initialize_elliptic_mb(struct _initialize_elliptic_mb p) {
  double a = (p.a) ? p.a : RADIUS;
  double b = (p.b) ? p.b : RADIUS;
  int nlp = (p.nlp) ? p.nlp : NLP;
  // double inclination = (p.inclination) ? p.inclination : 0.;
  p.mesh->nlp = nlp;
  p.mesh->nle = nlp;
  p.mesh->nodes = malloc(nlp*sizeof(lagNode));
  p.mesh->edges = malloc(nlp*sizeof(Edge));

  double alpha = 2*pi/(nlp);
  /** Fill the array of nodes */
  for(int i=0; i<nlp; i++) {
    p.mesh->nodes[i].pos.x = a*cos(alpha*i);
    p.mesh->nodes[i].pos.y = b*sin(alpha*i);
    p.mesh->nodes[i].edge_ids[0] = -1;
    p.mesh->nodes[i].edge_ids[1] = -1;
    foreach_dimension() {
      p.mesh->nodes[i].lagForce.x = 0.;
      p.mesh->nodes[i].lagVel.x = 0.;
    }
  }
  /** Fill the array of edges.
  For the last edge, the next vertex id is 0 */
  for(int i=0; i<p.mesh->nle; i++) {
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
  }
  /** The above procedure switches the two egde ids for the first node, which
  we correct below */
  p.mesh->nodes[0].edge_ids[0] = p.mesh->nle-1;
  p.mesh->nodes[0].edge_ids[1] = 0;
  correct_lag_pos(p.mesh);
  for(int i=0.; i<p.mesh->nle; i++) {
    p.mesh->edges[i].l0 = comp_length(p.mesh, i);
    p.mesh->edges[i].length = p.mesh->edges[i].l0;
  }

  #ifdef CAPS_VISCOSITY
    fraction(prevI, 1 - sq(x/a) - sq(y/b));
  #endif
}

#if dimension > 2

#define GET_LD(NODE) ((fabs(fabs(NODE.pos.x) - ll) < 1.e-8) ? 0 : \
  ((fabs(fabs(NODE.pos.y) - ll) < 1.e-8 ? 1 : 2)))
#define GET_LD_SIGN(NODE) ((GET_LD(NODE) == 0) ? sign(NODE.pos.x) : \
  ((GET_LD(NODE) == 1) ? sign(NODE.pos.y) : sign(NODE.pos.z)))
#define GET_SD(NODE) ((fabs(fabs(NODE.pos.x) - sl) < 1.e-8) ? 0 : \
  ((fabs(fabs(NODE.pos.y) - sl) < 1.e-8 ? 1 : 2)))
#define GET_SD_SIGN(NODE) ((GET_SD(NODE) == 0) ? sign(NODE.pos.x) : \
  ((GET_SD(NODE) == 1) ? sign(NODE.pos.y) : sign(NODE.pos.z)))
#define GET_ZD(NODE) ((fabs(NODE.pos.x) < 1.e-8) ? 0 : \
  ((fabs(NODE.pos.y) < 1.e-8 ? 1 : 2)))

bool edge_exists(lagMesh* mesh, int j, int k) {
  for(int i=0; i<mesh->nle; i++) {
    if ((mesh->edges[i].vertex_ids[0] == j && mesh->edges[i].vertex_ids[1] == k)
      || (mesh->edges[i].vertex_ids[0] == k
      && mesh->edges[i].vertex_ids[1] == j)) return true;
  }
  return false;
}

bool write_edge(lagMesh* mesh, int i, int j, int k) {
  if (edge_exists(mesh, j, k)) return false;
  else {
    mesh->edges[i].vertex_ids[0] = j;
    mesh->edges[i].vertex_ids[1] = k;
    for(int ii=0; ii<2; ii++) mesh->edges[i].triangle_ids[ii] = -1;
    mesh->nodes[j].edge_ids[mesh->nodes[j].nb_edges] = i;
    mesh->nodes[j].nb_edges++;
    mesh->nodes[j].neighbor_ids[mesh->nodes[j].nb_neighbors] = k;
    mesh->nodes[j].nb_neighbors++;
    mesh->nodes[k].edge_ids[mesh->nodes[k].nb_edges] = i;
    mesh->nodes[k].nb_edges++;
    mesh->nodes[k].neighbor_ids[mesh->nodes[k].nb_neighbors] = j;
    mesh->nodes[k].nb_neighbors++;
    return true;
  }
}

bool triangle_exists(lagMesh* mesh, int i, int j, int k) {
  for(int t=0; t<mesh->nlt; t++) {
    for(int a=0; a<3; a++) {
      if (mesh->triangles[t].node_ids[a] == i) {
        for(int b=0; b<3; b++) {
          if (b != a && mesh->triangles[t].node_ids[b] == j) {
            for(int c=0; c<3; c++) {
              if (c != a && c != b && mesh->triangles[t].node_ids[c] == k) {
                return true;
              }
            }
          }
        }
      }
    }
  }
  return false;
}

bool write_triangle(lagMesh* mesh, int tid, int i, int j, int k) {
  if (triangle_exists(mesh, i, j, k)) return false;
  else {
    mesh->triangles[tid].node_ids[0] = i;
    mesh->triangles[tid].node_ids[1] = j;
    mesh->triangles[tid].node_ids[2] = k;
    mesh->nodes[j].triangle_ids[mesh->nodes[i].nb_triangles] = tid;
    mesh->nodes[j].nb_triangles++;
    mesh->nodes[k].triangle_ids[mesh->nodes[k].nb_triangles] = tid;
    mesh->nodes[k].nb_triangles++;
    int c = 0;
    for(int a=0; a<3; a++) {
      int va = mesh->triangles[tid].node_ids[a];
      for(int b=0; b<3; b++) {
        if (b != a) {
          int vb = mesh->triangles[tid].node_ids[b];
          for(int m=0; m<mesh->nodes[va].nb_edges; m++) {
            for(int n=0; n<mesh->nodes[vb].nb_edges; n++) {
              if (mesh->nodes[va].edge_ids[m] == mesh->nodes[vb].edge_ids[n]) {
                mesh->triangles[tid].edge_ids[c] = mesh->nodes[va].edge_ids[m];
                c++;
                int p =
                (mesh->edges[mesh->nodes[va].edge_ids[m]].triangle_ids[0] == -1)
                ? 0 : 1;
                mesh->edges[mesh->nodes[va].edge_ids[m]].triangle_ids[p] = tid;
              }
            }
          }
        }
      }
    }
    return true;
  }
}

void initialize_icosahedron(struct _initialize_circular_mb p) {
  double radius = (p.radius) ? p.radius : RADIUS;
  coord shift;
  if (p.shift.x || p.shift.y || p.shift.z)
    {shift.x = p.shift.x; shift.y = p.shift.y; shift.z = p.shift.z;}
  else {shift.x = 0.; shift.y = 0.; shift.z = 0.;}
  p.mesh->nlp = 12;
  p.mesh->nodes = malloc(p.mesh->nlp*sizeof(lagNode));
  p.mesh->nle = 0;
  p.mesh->edges = malloc(30*sizeof(Edge));
  p.mesh->nlt = 0;
  p.mesh->triangles = malloc(20*sizeof(Triangle));

  /** Create the Lagrangian nodes */
  double sl, ll, r;
  r = sqrt(1. + sq(.5*(1. + sqrt(5))));
  sl = radius*1./r;
  ll = radius*.5*(1 + sqrt(5))/r;
  coord c[3] = {{0, sl, ll}, {sl, ll, 0}, {ll, 0, sl}};
  int s[9] = {0, 1, 1, -1, -1, -1, 1, -1, 1};
  for(int i=0; i<3; i++) {
    for(int j=0; j<4; j++) {
      p.mesh->nodes[i*4+j].nb_neighbors = 0;
      p.mesh->nodes[i*4+j].nb_edges = 0;
      p.mesh->nodes[i*4+j].nb_triangles = i;
      foreach_dimension()
        p.mesh->nodes[i*4+j].pos.x = c[i].x*
          s[(1+j+4*((fabs(c[i].x - sl) < 1.e-8) ? 0 : 1))*((fabs(c[i].x) < 1.e-8) ? 0 : 1)];
    }
  }

  /** Create the edges */
  for(int i=0; i<p.mesh->nlp; i++) {
    int my_gr = p.mesh->nodes[i].nb_triangles;
    int my_ld_sign = GET_LD_SIGN(p.mesh->nodes[i]);
    for(int j=0; j<p.mesh->nlp; j++) {
      if (p.mesh->nodes[j].nb_triangles == my_gr && j != i) {
        if (my_ld_sign == GET_LD_SIGN(p.mesh->nodes[j])) {
          if (write_edge(p.mesh, p.mesh->nle, i, j)) p.mesh->nle++;
        }
      }
      else if (GET_LD(p.mesh->nodes[j]) == GET_ZD(p.mesh->nodes[i])) {
        if (GET_SD_SIGN(p.mesh->nodes[j]) == GET_LD_SIGN(p.mesh->nodes[i])) {
          if (write_edge(p.mesh, p.mesh->nle, i, j)) p.mesh->nle++;
        }
      }
    }
  }

  /** Create the triangular faces */
  for(int i=0; i<p.mesh->nlp; i++) {
    p.mesh->nodes[i].nb_triangles = 0;
  }
  for(int i=0; i<p.mesh->nle; i++) {
    int nid[2];
    for(int j=0; j<2; j++) nid[j] = p.mesh->edges[i].vertex_ids[j];
    for(int j=0; j<p.mesh->nodes[nid[0]].nb_neighbors; j++) {
      for(int k=0; k<p.mesh->nodes[nid[1]].nb_neighbors; k++) {
        if (p.mesh->nodes[nid[0]].neighbor_ids[j] ==
          p.mesh->nodes[nid[1]].neighbor_ids[k])
          if (write_triangle(p.mesh, p.mesh->nlt, nid[0], nid[1],
            p.mesh->nodes[nid[0]].neighbor_ids[j])) {
            fprintf(stderr, "writing triangle %d %d %d\n", nid[0], nid[1],
              p.mesh->nodes[nid[0]].neighbor_ids[j]);
            p.mesh->nlt++;
          }
      }
    }
  }
  fprintf(stderr, "nlt = %d\n", p.mesh->nlt);
}

void initialize_spherical_mb(struct _initialize_circular_mb p) {
  initialize_icosahedron(p);

  double radius = (p.radius) ? p.radius : RADIUS;
  int nlp = (p.nlp) ? p.nlp : NLP;
  int nle, nlt;
  coord shift;
  if (p.shift.x || p.shift.y || p.shift.z)
    {shift.x = p.shift.x; shift.y = p.shift.y; shift.z = p.shift.z;}
  else {shift.x = 0.; shift.y = 0.; shift.z = 0.;}
  int i=1;
  int cn = 12;
  while(cn < nlp) {
    nle = 30*pow(3,i);
    cn = nle/3; //this is an upper bound
    i++;
  }
  nlp = cn;
  nlt = 20*pow(4,i);
  fprintf(stderr, "Number of Lagrangian nodes = %d\n", nlp);
  fprintf(stderr, "Number of Lagrangian edges = %d\n", nle);
  fprintf(stderr, "Number of Lagrangian triangles = %d\n", nlt);
  // p.mesh->nodes = malloc(nlp*sizeof(lagNode));
}

#endif
