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

typedef struct intcoord {
  int x, y, z;
} intcoord;

intcoord get_geoloc(coord xd, coord yd, double Delta) {
  intcoord res;
  intcoord cgeoloc;
  cgeoloc.y = sign(xd.y)*round(fabs(xd.y)/Delta);
  cgeoloc.x = sign(xd.x)*round(fabs(xd.x/Delta + .5));
  res.x = (cgeoloc.x + 2)*5 + cgeoloc.y + 2;
  cgeoloc.y = sign(yd.y)*round(fabs(yd.y)/Delta);
  cgeoloc.x = sign(yd.x)*round(fabs(yd.x/Delta + .5));
  res.y = (cgeoloc.x + 2)*5 + cgeoloc.y + 2;
  return res;
}

void construct_divG(scalar divG, lagMesh* mesh) {
  comp_normals(mesh);
  comp_mb_stretch(mesh);
  foreach() {
    foreach_dimension() G.x[] = 0.;
    divG[] = 0.;
  }
  for(int i=0; i<mesh->nle; i++) {
    coord midpoint;
    int node_id[2];
    for(int j=0; j<2; j++) node_id[j] = mesh->edges[i].vertex_ids[j];
    midpoint.x = .5*GENERAL_1DAVG(mesh->nodes[node_id[0]].pos.x,
      mesh->nodes[node_id[1]].pos.x);
    midpoint.y = .5*GENERAL_1DAVG(mesh->nodes[node_id[0]].pos.y,
      mesh->nodes[node_id[1]].pos.y);
    double length = mesh->edges[i].st*mesh->edges[i].l0;
    int ni = 0;
    int nstencil = 0;
    intcoord incr_stencil[25];
    for(int k=0; k<25; k++) {
      incr_stencil[k].x = 0;
      incr_stencil[k].y = 0;
    }
    while (nstencil < 50 && ni < 2) {
      foreach_cache(mesh->nodes[node_id[ni]].stencil) {
        /** Since we made the choice to have G as a face vector, we need to
        check successively if the left and bottom edges of a neighboring cell
        are inside the stencil of the regularized Dirac function. */
        coord xpos = {x - .5*Delta, y};
        coord ypos = {x, y - .5*Delta};
        coord xdist, ydist;
        foreach_dimension() {
          xdist.x = xpos.x - midpoint.x;
          ydist.x = ypos.x - midpoint.x;
        }
        intcoord geoloc = get_geoloc(xdist, ydist, Delta);
        if ((fabs(xdist.x) <= 2*Delta) && (fabs(xdist.y) <= 2*Delta) &&
          incr_stencil[geoloc.x].x == 0 ) {
          double prefactor = (1 + cos(.5*pi*xdist.x/Delta))*
            (1 + cos(.5*pi*xdist.y/Delta))*length/(16.*sq(Delta));
          G.x[] -= prefactor*mesh->edges[i].normal.x;
          nstencil++;
          incr_stencil[geoloc.x].x = 1;
        }
        if ((fabs(ydist.x) <= 2*Delta) && (fabs(ydist.y) <= 2*Delta) &&
          (incr_stencil[geoloc.y].y == 0)) {
          double prefactor = (1 + cos(.5*pi*ydist.x/Delta))*
            (1 + cos(.5*pi*ydist.y/Delta))*length/(16.*sq(Delta));
          G.y[] -= prefactor*mesh->edges[i].normal.y;
          nstencil++;
          incr_stencil[geoloc.y].y = 1;
        }
      }
      ni++;
    }
    fprintf(stderr, "nstencil=%d\t, incr_stencil.x=[%d", nstencil, incr_stencil[0].x);
    for(int is=1; is<25; is++) fprintf(stderr, ", %d", incr_stencil[is].x);
    fprintf(stderr, "]\n");
    if (ni == 2) fprintf(stderr, "Error: trying to access a third edge node\n");
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
