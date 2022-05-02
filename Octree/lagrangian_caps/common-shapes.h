/** Some common initial shapes for vesicles */

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
    p.mesh->edges[i].st = 1.;
  }

  #ifdef MUP
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
    p.mesh->edges[i].st = 1.;
  }

  #ifdef MUP
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
    p.mesh->edges[i].st = 1.;
  }

  #ifdef MUP
    fraction(prevI, 1 - sq(x/a) - sq(y/b));
  #endif
}
