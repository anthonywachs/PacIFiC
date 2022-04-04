#ifndef MUP
  #define MUP 1.;
#endif
#ifndef MUC
  #define MUC 1.;
#endif

vector G[];
void construct_divG(scalar divG, lagMesh* mesh) {
  comp_edge_normals(mesh);
  foreach() {
    foreach_dimension() G.x[] = 0.;
    coord sdist;
    coord cpos, pos;
    cpos.x = x; cpos.y = y;
    for(int i=0; i<mesh->nle; i++) {
      coord midpoint;
      int node_id[2];
      for(int j=0; j<2; j++) node_id[j] = mesh->edges[i].vertex_ids[j];
      midpoint.x = .5*(mesh->nodes[node_id[0]].pos.x
        + mesh->nodes[node_id[1]].pos.x);
      midpoint.y = .5*(mesh->nodes[node_id[0]].pos.y
        + mesh->nodes[node_id[1]].pos.y);
      // coord normal;
      // double theta = fabs(midpoint.x) > 0. ? atan2(midpoint.y, midpoint.x) :
      //   atan2(midpoint.y, midpoint.x + 1.e-6);
      // normal.x = cos(theta);
      // normal.y = sin(theta);
      if ((fabs(cpos.x - midpoint.x)<Delta/2.)
        && (fabs(cpos.y - midpoint.y)<Delta/2.)) {
        double length = mesh->edges[i].st*mesh->edges[i].l0;
        foreach_neighbor() {
          pos.x = x; pos.y = y;
          foreach_dimension() sdist.x = sq(pos.x - midpoint.x);
          if ((sdist.x <= sq(2*Delta)) && (sdist.y <= sq(2*Delta))) {
            /** We need the inward normal instead of the outward normal, hence
            the minus sign */
            foreach_dimension() G.x[] -= //why do I have to multiply by 2??
              (1 + cos(.5*pi*(pos.x - midpoint.x)/Delta))
              *(1 + cos(.5*pi*(pos.y - midpoint.y)/Delta))
              *mesh->edges[i].normal.x*length/(16.*sq(Delta));
              // *normal.x*(0.7853981633974483/mesh->nle)/(16.*sq(Delta));
          }
        }
      }
    }
  }
  foreach() {
    divG[] = 0.;
    foreach_dimension() divG[] += (G.x[1] - G.x[-1])/(2*Delta);
  }
}

double muc, mup;
scalar I[];
scalar prevI[];
scalar divG[];
event defaults (i = 0) {
  mu = new face vector;
  mup = MUP;
  muc = MUC;
}
I[top] = dirichlet(0.);
I[bottom] = dirichlet(0.);
prevI[top] = dirichlet(0.);
prevI[bottom] = dirichlet(0.);
divG[top] = dirichlet(0.);
divG[bottom] = dirichlet(0.);

event properties (i++) {
  construct_divG(divG, &(mbs.mb[0]));
  poisson (I, divG, tolerance = 1.e-10);

  // Simple clamping of I:
  foreach() {
    if (fabs(divG[]) > 1.e-10) {
      I[] = clamp(I[], 0, 1);
      prevI[] = I[];
    }
    else {
      prevI[] = round(prevI[]);
      I[] = prevI[];
    }
  }

  foreach_face() {
    face vector muv = mu;
    muv.x[] = mup + (muc - mup)*.5*(I[] + I[-1]);
    // muv.x[] = mup + (muc - mup)*.5*(I[]/maxI + I[-1]/maxI);
  }
}
