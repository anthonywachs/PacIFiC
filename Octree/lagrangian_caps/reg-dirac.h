/**
# Definition of the IBM regularized Dirac
*/

#define BGHOSTS 2

static void change_cache_entry(Cache* s, int i, Point pt, int flag) {
  if (i > s->n) fprintf(stderr, "Error: Cache index out of range.\n");
  s->p[i].i = pt.i;
  s->p[i].j = pt.j;
  #if dimension >= 3
    s->p[i].k = pt.k;
  #endif
  s->p[i].level = pt.level;
  s->p[i].flags = flag;
}

void generate_lag_stencils(lagMesh* mesh) {
  for(int i=0; i<mesh->nlp; i++) {
    int c = 0;
    double delta = (L0/(1 << grid->maxdepth));
    for(int ni=-2; ni<=2; ni++) {
      for(int nj=-2; nj<=2; nj++) {
        Point point;
        point = locate(mesh->nodes[i].pos.x + ni*delta,
          mesh->nodes[i].pos.y + nj*delta);
        if (point.level >= 0 && point.level != grid->maxdepth)
          fprintf(stderr, "Warning: Lagrangian stencil not fully resolved.\n");
        change_cache_entry(&(mesh->nodes[i].stencil), c, point, 0);
        #if _MPI
          if (ni == 0 && nj == 0) {
            if (point.level >= 0) mesh->nodes[i].pid = cell.pid;
            else mesh->nodes[i].pid = -1;
          }
        #endif
        c++;
      }
    }
  }
}


/**
The function below fills the scalar $forcing$ with the regularized source term
of magnitude $force$ located at $x0$. This assumes that $forcing$ is "clean",
i.e. that it is zero everywhere (or that if there is already some non-zero terms
the intention is to include them in the forcing).
*/
void lag2eul(vector forcing, lagMesh* mesh) {
  for(int i=0; i<mesh->nlp; i++) {
    foreach_cache(mesh->nodes[i].stencil) {
      if (point.level >= 0) {
        coord sdist;
        sdist.x = sq(x - mesh->nodes[i].pos.x);
        sdist.y = sq(y - mesh->nodes[i].pos.y);
        if (sdist.x <= sq(2*Delta) && sdist.y <= sq(2*Delta)) {
          double weight = (1 + cos(.5*pi*(x - mesh->nodes[i].pos.x)/Delta))*
            (1 + cos(.5*pi*(y - mesh->nodes[i].pos.y)/Delta))/(16.*sq(Delta));
          foreach_dimension() forcing.x[] += weight*mesh->nodes[i].lagForce.x;
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
  for(int i=0; i<mesh->nlp; i++) {
    foreach_dimension() mesh->nodes[i].lagVel.x = 0.;
    foreach_cache(mesh->nodes[i].stencil) {
      if (point.level >= 0) {
        coord dist;
        dist.x = GENERAL_1DIST(x, mesh->nodes[i].pos.x);
        dist.y = GENERAL_1DIST(y, mesh->nodes[i].pos.y);
        if (sq(dist.x) <= sq(2*Delta) && sq(dist.y) <= sq(2*Delta)) {
          double weight = (1 + cos(.5*pi*dist.x/Delta))*
            (1 + cos(.5*pi*dist.y/Delta))/16.;
          foreach_dimension() {
            mesh->nodes[i].lagVel.x += weight*u.x[];
          }
        }
      }
    }
  }

  #if _MPI
    if (mpi_npe > 1) reduce_lagVel(mesh);
  #endif
}

scalar stencils[];
void tag_ibm_stencils(lagMesh* mesh) {
  foreach() stencils[] = 0.;
  for(int i=0; i<mesh->nlp; i++) {
    foreach_cache(mesh->nodes[i].stencil) {
      if (point.level >= 0) {
        coord dist;
        dist.x = GENERAL_1DIST(x, mesh->nodes[i].pos.x);
        dist.y = GENERAL_1DIST(y, mesh->nodes[i].pos.y);
        if (sq(dist.x) <= sq(2*Delta) && sq(dist.y) <= sq(2*Delta))
          stencils[] = sq(dist.x + dist.y)/sq(2.*Delta)*(2.+noise());
      }
    }
  }
  boundary({stencils});
}
