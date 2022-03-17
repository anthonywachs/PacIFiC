/**
Definition of the IBM regularized Dirac
*/

#define BGHOSTS 2

/**
The function below fills the scalar $forcing$ with the regularized source term
of magnitude $force$ located at $x0$. This assumes that $forcing$ is "clean",
i.e. that it is zero everywhere (or that if there is already some non-zero terms
the intention is to include them in the forcing).

TODO:
* implement this function for an array of source terms, i.e. an array of
coord $x0$ and an array of magnitudes $force$.
* in case of trees, make sure that a cell located far away with a big Delta is
not included in the stencil if it shouldn't.
*/
void lag2eul(scalar forcing, coord p, double force) {
  foreach() {
    coord sdist;
    coord pos;
    pos.x = x; pos.y = y;
    foreach_dimension() sdist.x = sq(pos.x - p.x);
    if ((sdist.x <= sq(2*Delta)) && (sdist.y <= sq(2*Delta))) {
      forcing[] += (1 + cos(.5*pi*(pos.x - p.x)/Delta))*(1 +
        cos(.5*pi*(pos.y - p.y)/Delta))*force/(16*sq(Delta));
      // forcing[] += (1 + cos(.5*pi*(pos.x - p.x)/Delta))*(1 +
      //   cos(.5*pi*(pos.y - p.y)/Delta))*force/16.;
    }
  }
}

/**
The function below interpolates the eulerian velocities onto the vertices of
the Lagrangian mesh.
*/
void eul2lag(lagMesh* mesh) {
  for(int i=0; i<mesh->nlp; i++){
    foreach_dimension() mesh->lagVel[i].x = 0.;
  }
  foreach() {
    coord sdist;
    coord cpos, pos;
    cpos.x = x; cpos.y = y;
    for(int i=0; i<mesh->nlp; i++) {
      if ((fabs(cpos.x - mesh->vertices[i].x)<Delta/2.)
        && (fabs(cpos.y - mesh->vertices[i].y)<Delta/2.)) {
        foreach_neighbor() {
          pos.x = x; pos.y = y;
          foreach_dimension() sdist.x = sq(pos.x - mesh->vertices[i].x);
          if ((sdist.x <= sq(2*Delta)) && (sdist.y <= sq(2*Delta))) {
            foreach_dimension() {
              mesh->lagVel[i].x += (1 + cos(.5*pi*(pos.x - mesh->vertices[i].x)
                /Delta))*(1 + cos(.5*pi*(pos.y - mesh->vertices[i].y)/Delta))
                *u.x[]/16.;
            }
          }
        }
      }
    }
  }
}
