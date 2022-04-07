#ifndef MUP
  #define MUP 1.;
#endif
#ifndef MUC
  #define MUC 1.;
#endif

/** We define the "grid gradient" G, according to Tryggvason, JCP 2003. Defining
G as a face vector as opposed to a centered vector will lead to a sharper
viscosity jump - which we are happy with -, although it complicates slightly
the implementation below.*/
face vector G[];
void construct_divG(scalar divG, lagMesh* mesh) {
  comp_normals(mesh);
  comp_mb_stretch(mesh);
  foreach() {
    foreach_dimension() G.x[] = 0.;
    divG[] = 0.;
  }
  foreach() {
    coord cpos = {x, y};
    for(int i=0; i<mesh->nle; i++) {
      coord midpoint;
      int node_id[2];
      for(int j=0; j<2; j++) node_id[j] = mesh->edges[i].vertex_ids[j];
      midpoint.x = .5*(mesh->nodes[node_id[0]].pos.x
        + mesh->nodes[node_id[1]].pos.x);
      midpoint.y = .5*(mesh->nodes[node_id[0]].pos.y
        + mesh->nodes[node_id[1]].pos.y);
      if ((fabs(cpos.x - midpoint.x)<=Delta/2.)
        && (fabs(cpos.y - midpoint.y)<=Delta/2.)) {
        double length = mesh->edges[i].st*mesh->edges[i].l0;
        foreach_neighbor() {
          /** Since we made the choice to have G as a face vector, we need to
          check successively if the left and bottom edges of a neighboring cell
          are inside the stencil of the regularized Dirac function. */
          coord xpos = {x - .5*Delta, y};
          coord ypos = {x, y - .5*Delta};
          coord sxdist, sydist;
          foreach_dimension() {
            sxdist.x = sq(xpos.x - midpoint.x);
            sydist.x = sq(ypos.x - midpoint.x);
          }
          if ((sxdist.x <= sq(2*Delta)) && (sxdist.y <= sq(2*Delta))) {
            double prefactor = (1 + cos(.5*pi*(xpos.x - midpoint.x)/Delta))*
              (1 + cos(.5*pi*(xpos.y - midpoint.y)/Delta))*length/
              (16.*sq(Delta));
            G.x[] -= prefactor*mesh->edges[i].normal.x;
          }
          if ((sydist.x <= sq(2*Delta)) && (sydist.y <= sq(2*Delta))) {
            double prefactor = (1 + cos(.5*pi*(ypos.x - midpoint.x)/Delta))*
              (1 + cos(.5*pi*(ypos.y - midpoint.y)/Delta))*length/
              (16.*sq(Delta));
            G.y[] -= prefactor*mesh->edges[i].normal.y;
          }
        }
      }
    }
  }
  boundary((scalar*){G});
  foreach() {
    foreach_dimension() divG[] += (G.x[1] - G.x[0])/(Delta);
  }
  boundary({divG});
}

double muc, mup;
scalar I[];
scalar prevI[];
scalar divG[];
event defaults (i = 0) {
  mu = new face vector;
  mup = MUP;
  muc = MUC;
  I[top] = dirichlet(0.);
  I[bottom] = dirichlet(0.);
  prevI[top] = dirichlet(0.);
  prevI[bottom] = dirichlet(0.);
  divG[top] = dirichlet(0.);
  divG[bottom] = dirichlet(0.);
}

event properties (i++) {
  construct_divG(divG, &(mbs.mb[0]));
  poisson(I, divG, tolerance = 1.e-10);

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
    I[] = clamp(I[], 0, 1);
  }
  boundary({I, prevI});

  foreach_face() {
    face vector muv = mu;
    muv.x[] = mup + (muc - mup)*.5*(I[] + I[-1]);
  }
}
