/**
# Variable viscosity across a membrane
*/

#ifndef MUP
  #define MUP 1.;
#endif
#ifndef MUC
  #define MUC 1.;
#endif
#define CAPS_VISCOSITY 1.

/** We define the "grid gradient" G, according to Tryggvason, JCP 2003.*/
vector G[];
// face vector G[];
void construct_divG(scalar divG, lagMesh* mesh) {
  comp_normals(mesh);
  comp_mb_stretch(mesh);
  // foreach_face() foreach_dimension() G.x[] = 0.;
  foreach() {
    foreach_dimension() G.x[] = 0.;
    divG[] = 0.;
  }
  boundary((scalar*){divG, G});

  for(int i=0; i<mesh->nle; i++) {
    // compute the grid gradient on the midpoint of the edge
    coord gg; // grid gradient
    int en[2];
    en[0] = mesh->edges[i].vertex_ids[0];
    en[1] = mesh->edges[i].vertex_ids[1];
    foreach_dimension() gg.x = mesh->edges[i].normal.x*mesh->edges[i].length;
    for(int j=0; j<2; j++) {
      foreach_cache(mesh->nodes[en[j]].stencil) {
        //spread half the grid gradient of the edge from each of its nodes
        if (point.level >= 0) {
          // coord xpos = {x - .5*Delta, y};
          // coord ypos = {x, y - .5*Delta};
          // coord xdist, ydist;
          // foreach_dimension() {
          //   xdist.x = GENERAL_1DIST(xpos.x, mesh->nodes[en[j]].pos.x);
          //   ydist.x = GENERAL_1DIST(ypos.x, mesh->nodes[en[j]].pos.x);
          // }
          // if (sq(xdist.x) <= sq(2*Delta) && sq(xdist.y) <= sq(2*Delta)) {
          //   double weight =
          //     (1 + cos(.5*pi*xdist.x/Delta))*(1 + cos(.5*pi*xdist.y/Delta))/
          //     (16.*sq(Delta));
          //   G.x[] -= weight*.5*gg.x;
          // }
          // if (sq(ydist.x) <= sq(2*Delta) && sq(ydist.y) <= sq(2*Delta)) {
          //   double weight =
          //     (1 + cos(.5*pi*ydist.x/Delta))*(1 + cos(.5*pi*ydist.y/Delta))/
          //     (16.*sq(Delta));
          //   G.y[] -= weight*.5*gg.y;
          // }
        coord dist;
          dist.x = GENERAL_1DIST(x, mesh->nodes[en[j]].pos.x);
          dist.y = GENERAL_1DIST(y, mesh->nodes[en[j]].pos.y);
          if (sq(dist.x) <= sq(2*Delta) && sq(dist.y) <= sq(2*Delta)) {
            double weight =
              (1 + cos(.5*pi*dist.x/Delta))*(1 + cos(.5*pi*dist.y/Delta))/
              (16.*sq(Delta));
            foreach_dimension() G.x[] -= weight*.5*gg.x;
          }
        }
      }
    }
  }
  boundary((scalar*){G});
  foreach() foreach_dimension() divG[] += (G.x[1] - G.x[-1])/(2.*Delta);
  // foreach() foreach_dimension() divG[] += (G.x[1] - G.x[])/(Delta);
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
