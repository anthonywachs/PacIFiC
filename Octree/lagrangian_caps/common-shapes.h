/** Some common initial shapes for vesicles */

struct _initialize_circular_mb {
  lagMesh* mesh;
  double radius;
  double nlp;
};

void initialize_circular_mb(struct _initialize_circular_mb p) {
  double radius = (p.radius) ? p.radius : RADIUS;
  int nlp = (p.nlp) ? p.nlp : NLP;
  p.mesh->nlp = nlp;
  p.mesh->nle = nlp;
  p.mesh->nodes = malloc(nlp*sizeof(lagNode));
  p.mesh->edges = malloc(nlp*sizeof(Edge));

  double alpha = 2*pi/(nlp);
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
}

void initialize_biconcave_mb(struct _initialize_circular_mb p) {
  double radius = (p.radius) ? p.radius : RADIUS;
  int nlp = (p.nlp) ? p.nlp : NLP;
  double c = 1.3858189;
  p.mesh->nlp = nlp;
  p.mesh->nle = nlp;
  p.mesh->nodes = malloc(nlp*sizeof(lagNode));
  p.mesh->edges = malloc(nlp*sizeof(Edge));

  double alpha = 2*pi/(nlp);
  /** Fill the array of nodes */
  for(int i=0; i<nlp; i++) {
    p.mesh->nodes[i].pos.y = radius*c*cos(2*pi-alpha*i) - 5*radius;
    p.mesh->nodes[i].pos.x = .5*radius*c*sin(2*pi-alpha*i)*(0.207 + 2.003*sq(cos(2*pi-alpha*i)) - 1.123*sq(sq(cos(2*pi-alpha*i))));
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
}
