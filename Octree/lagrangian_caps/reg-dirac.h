/**
Definition of the IBM regularized Dirac
*/

#define BGHOSTS 2

/**
The function below fills the scalar $forcing$ with the regularized source term
of magnitude $force$ located at $x0$. This assumes that $forcing$ is "clean",
i.e. that it is zero everywhere (or that if there is already some non-zero terms
the intention is to include them in the forcing).
*/
void lag2eul(vector forcing, lagMesh* mesh) {
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
            foreach_dimension() forcing.x[] +=
                (1 + cos(.5*pi*(pos.x - mesh->nodes[i].pos.x)/Delta))
                *(1 + cos(.5*pi*(pos.y - mesh->nodes[i].pos.y)/Delta))
                *mesh->nodes[i].lagForce.x/(16.*sq(Delta));
          }
        }
      }
    }
  }
}

/**
The function below interpolates the eulerian velocities onto the nodes of
the Lagrangian mesh.
*/
void eul2lag(lagMesh* mesh) {
  for(int i=0; i<mesh->nlp; i++){
    foreach_dimension() mesh->nodes[i].lagVel.x = 0.;
  }
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
            foreach_dimension() {
              mesh->nodes[i].lagVel.x += (1 + cos(.5*pi*(pos.x -
                mesh->nodes[i].pos.x)/Delta))
                *(1 + cos(.5*pi*(pos.y - mesh->nodes[i].pos.y)/Delta))
                *u.x[]/16.;
            }
          }
        }
      }
    }
  }
}
