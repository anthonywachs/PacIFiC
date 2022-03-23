#ifndef MUP
  #define MUP 1.;
#endif
#ifndef MUC
  #define MUC 1.;
#endif

void construct_divG(scalar divG, lagMesh* mesh) {
  vector G[];
  comp_normals(mesh);
  foreach() {
    coord sdist;
    coord cpos, pos;
    cpos.x = x; cpos.y = y;
    for(int i=0; i<mesh->nlp; i++) {
      if ((fabs(cpos.x - mesh->nodes[i].pos.x)<Delta/2.)
        && (fabs(cpos.y - mesh->nodes[i].pos.y)<Delta/2.)) {
        foreach_neighbor() {
          pos.x = x; pos.y = y;
          foreach_dimension() sdist.x = sq(pos.x - mesh->nodes[i].pos.x);
          if ((sdist.x <= sq(2*Delta)) && (sdist.y <= sq(2*Delta))) {
            foreach_dimension() G.x[] +=
                (1 + cos(.5*pi*(pos.x - mesh->nodes[i].pos.x)/Delta))
                *(1 + cos(.5*pi*(pos.y - mesh->nodes[i].pos.y)/Delta))
                *mesh->nodes[i].normal.x/(16.*sq(Delta));
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
event defaults (i = 0) {
  mu = new face vector;
  mup = MUP;
  muc = MUC;
}
I[left] = dirichlet(0.);
I[top] = dirichlet(0.);

event properties (i++) {
  scalar divG[];

  construct_divG(divG, &(mbs.mb[0]));
  mgstats s = poisson (I, divG, tolerance = 1.e-6);

  foreach_face() {
    face vector muv = mu;
    muv.x[] = mup + (muc - mup)*.5*(I[] + I[-1]);
  }
  view(fov = 20, bg = {1,1,1});
  clear();
  cells();
  squares("I", linear=false);
  draw_lag(&(mbs.mb[0]), lc = {1.,0.,0.}, vc = {1., 0., 0.}, vs = 3.);
  save("indicator.png");
  view(fov = 20, bg = {1,1,1});
  clear();
  cells();
  squares("mu.x", linear=false);
  draw_lag(&(mbs.mb[0]), lc = {1.,0.,0.}, vc = {1., 0., 0.}, vs = 3.);
  save("viscosity.png");
}
