/**
# Definition of the IBM regularized Dirac
*/

#define BGHOSTS 2

static void change_cache_entry(Cache* s, int i, Point pt, int flag) {
  if (i > s->n) fprintf(stderr, "Error: Cache index out of range.\n");
  s->p[i].i = pt.i;
  s->p[i].j = pt.j;
  #if dimension > 2
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
        #if dimension < 3
        Point point = locate(mesh->nodes[i].pos.x + ni*delta,
          mesh->nodes[i].pos.y + nj*delta);
        #else
        for(int nk=-2; nk<=2; nk++) {
          Point point = locate(mesh->nodes[i].pos.x + ni*delta,
            mesh->nodes[i].pos.y + nj*delta, mesh->nodes[i].pos.z + nk*delta);
        #endif
        if (point.level >= 0 && point.level != grid->maxdepth)
          fprintf(stderr, "Warning: Lagrangian stencil not fully resolved.\n");
        change_cache_entry(&(mesh->nodes[i].stencil), c, point, 0);
        #if _MPI
        #if dimension < 3
        if (ni == 0 && nj == 0) {
        #else
        if (ni == 0 && nj == 0 && nk == 0) {
        #endif
          if (point.level >= 0) mesh->nodes[i].pid = cell.pid;
          else mesh->nodes[i].pid = -1;
        }
        #endif
        c++;
        #if dimension == 3
        }
        #endif
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
        coord dist;
        dist.x = GENERAL_1DIST(x, mesh->nodes[i].pos.x);
        dist.y = GENERAL_1DIST(y, mesh->nodes[i].pos.y);
        #if dimension > 2
        dist.z = GENERAL_1DIST(z, mesh->nodes[i].pos.z);
        #endif
        #if dimension < 3
        if (fabs(dist.x) <= 2*Delta && fabs(dist.y) <= 2*Delta) {
          double weight =
            (1 + cos(.5*pi*dist.x/Delta))*(1 + cos(.5*pi*dist.y/Delta))
            /(sq(4*Delta));
        #else
        if (fabs(dist.x) <= 2*Delta && fabs(dist.y) <= 2*Delta &&
          fabs(dist.z) <= 2*Delta) {
          double weight =
            (1 + cos(.5*pi*dist.x/Delta))*(1 + cos(.5*pi*dist.y/Delta))
            *(1 + cos(.5*pi*dist.z/Delta))/(cube(4*Delta));
        #endif
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
  for(int ii=0; ii<mesh->nlp; ii++) {
    foreach_dimension() mesh->nodes[ii].lagVel.x = 0.;
    foreach_cache(mesh->nodes[ii].stencil) {
      if (point.level >= 0) {
        coord dist;
        dist.x = GENERAL_1DIST(x, mesh->nodes[ii].pos.x);
        dist.y = GENERAL_1DIST(y, mesh->nodes[ii].pos.y);
        #if dimension > 2
        dist.z = GENERAL_1DIST(z, mesh->nodes[ii].pos.z);
        #endif
        #if dimension < 3
        if (fabs(dist.x) <= 2*Delta && fabs(dist.y) <= 2*Delta) {
          double weight = (1 + cos(.5*pi*dist.x/Delta))*
            (1 + cos(.5*pi*dist.y/Delta))/16.;
        #else
        if (fabs(dist.x) <= 2*Delta && fabs(dist.y) <= 2*Delta
          && fabs(dist.z) <= 2*Delta) {
          double weight = (1 + cos(.5*pi*dist.x/Delta))*
            (1 + cos(.5*pi*dist.y/Delta))*(1 + cos(.5*pi*dist.z/Delta))/64.;
        #endif
        foreach_dimension() mesh->nodes[ii].lagVel.x += weight*u.x[];
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
        #if dimension > 2
        dist.z = GENERAL_1DIST(z, mesh->nodes[i].pos.z);
        #endif
        #if dimension < 3
        if (fabs(dist.x) <= 2*Delta && fabs(dist.y) <= 2*Delta) {
          stencils[] = sq(dist.x + dist.y)/sq(2.*Delta)*(2.+noise());
        #else
        if (fabs(dist.x) <= 2*Delta && fabs(dist.y) <= 2*Delta
          && fabs(dist.z) <= 2*Delta) {
          stencils[] = sq(dist.x + dist.y + dist.z)/cube(2.*Delta)*(2.+noise());
        #endif
        }
      }
    }
  }
  #if OLD_QCC
  boundary({stencils});
  #endif
}
