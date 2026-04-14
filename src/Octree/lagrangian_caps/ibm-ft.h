/**
# Definition of the IBM regularized Dirac

The Lagrangian and Eulerian meshes communicate thanks to regularized Dirac-delta
functions: this is the core of the Immersed Boundary Method (IBM) as introduced
by [Peskin](Peskin1997).

This implementation of the IBM in Basilisk relies on the Cache structure to
efficiently loop through the Eulerian cells in the vicinity of the Lagrangian
nodes. Since the Cache structure is only defined for trees, the current
implementation is \textbf{only compatible with quadtree and octree grids, and
not with Cartesian nor multigrids}.
*/

#define BGHOSTS 2 // Having two layers of ghost cells should not be mandatory since the implementation relies on Caches

/* We re-design the locate_stencil function in the tree structure */
#if !(MULT_GRID == 1)
  Point locate_stencil (double xp = 0, double yp = 0, double zp = 0)
  {
    // We assume all stencils are at the maximal level 
    int l = grid->maxdepth;

      Point point = {0};
      point.level = l;
      int n = 1 << point.level;
      point.i = (xp - X0)/L0*L0_ratio.x*n + GHOSTS;
  #if dimension >= 2
      point.j = (yp - Y0)/L0*L0_ratio.y*n + GHOSTS;
  #endif
  #if dimension >= 3
      point.k = (zp - Z0)/L0*L0_ratio.z*n + GHOSTS;
  #endif
      if (point.i >= 0 && point.i < n + 2*GHOSTS
  #if dimension >= 2
    && point.j >= 0 && point.j < n + 2*GHOSTS
  #endif
  #if dimension >= 3
    && point.k >= 0 && point.k < n + 2*GHOSTS
  #endif
    ) {
        if (allocated(0) && is_local(cell) && is_leaf(cell))
          return point;
    }
    Point lpoint = {0};
    lpoint.level = -1;
    return lpoint;
  }
#endif

struct generate_lag_stencils_one_caps_type{
  lagMesh* mesh;
  bool no_warning;
};
/**
The function below loops through the Lagrangian nodes and "caches" the Eulerian
cells in a 5x5(x5) stencil around each node. In case of parallel simulations,
the cached cells are tagged with the process id.
*/

trace
void generate_lag_stencils_one_caps(lagMesh* mesh, bool no_warning) 
{
//  lagMesh* mesh = p.mesh;
  
//  bool no_warning = p.no_warning;
  for(int i=0; i<mesh->nln; i++) 
  {
    mesh->nodes[i].stencil.n = 0;
    mesh->nodes[i].eulcell.n = 0;
    
    /**
    The current implementation assumes that the Eulerian cells around Lagrangian
    node are all at the maximum level.
    */

    /*The definition of maxdepth is one level higher in multigrid. Dimensions denote the dimensions of the proc arrangement in each direction, 
    note that L0 is always defined along x-axis,  */
    #if MULT_GRID == 1   
    #define locate_stencil(X, Y, Z) locate(X, Y, Z)
      double delta = (L0/(1 << grid->maxdepth)/Dimensions.x);
  #else
      double delta = (L0/(1 << grid->maxdepth));
    #endif

    int n_half;
    #if STENCIL_TYPE == 3
    n_half = 1;
    #else 
    n_half = 2;
    #endif

    Cache c = {0};

    for(int ni = -n_half; ni <= n_half; ni++) {
      for(int nj = -n_half; nj <= n_half; nj++) {
        #if dimension < 3
        Point point = locate_stencil(POS_PBC_X(mesh->nodes[i].pos.x + ni*delta),
          POS_PBC_Y(mesh->nodes[i].pos.y + nj*delta));
        #else
        for(int nk = -n_half; nk <= n_half; nk++) {
          Point point = locate_stencil(POS_PBC_X(mesh->nodes[i].pos.x + ni*delta),
            POS_PBC_Y(mesh->nodes[i].pos.y + nj*delta),
            POS_PBC_Z(mesh->nodes[i].pos.z + nk*delta));
        #endif
        if (!(no_warning) && point.level >= 0 && point.level != grid->maxdepth)
          fprintf(stderr, "Warning: Lagrangian stencil not fully resolved.\n");
        cache_append(&(mesh->nodes[i].stencil), point, 0);

        #if _MPI
        #if dimension < 3
        if (ni == 0 && nj == 0) {
        #else
        if (ni == 0 && nj == 0 && nk == 0) {
        #endif
          if (point.level >= 0) cache_append(&(mesh->nodes[i].eulcell), point, 0);
          if (point.level >= 0) cache_append( &c, point, 0);
          if (point.level >= 0) mesh->nodes[i].pid = pid();
          else mesh->nodes[i].pid = -1;
        }
        #endif
        #if dimension == 3
        }
        #endif
      }
    }

    /* We color the Index_lagnode field with their cap_id + 1 , 
    the vector Index_lag_id can contain up to 3 indices of the lag nodes, 
    currently we only save one lagnode */
    foreach_cache(c)
    {
      if ((int)Index_lagnode[] == -1)
      {
        Index_lagnode[] = mesh->cap_id;
        // if((int)Index_lag_id.x[] < 0)
        Index_lag_id.x[] = i;
        Index_lag_id.y[] = -1;
        Index_lag_id.y[] = -1;
      }
      else if ((int)Index_lagnode[] > -1 && (int)Index_lagnode[] != mesh->cap_id)// if the current cell has a node from another cap
      {
        Index_lag_id.y[] = mesh->cap_id;
        Index_lag_id.z[] = i; 
      }
    }
     free(c.p); 
  }
}

struct generate_lag_stencils_type
{
  bool no_warning;
}; 


trace
void generate_lag_stencils(bool no_warnings) {
  
  /*Clear the index field before generating stencils */  
  foreach()
  {
    if (cm[] > 1.e-20) 
    { Index_lagnode[] = -1;
      foreach_dimension() Index_lag_id.x[] = -1;
    }
  }

  for(int k=0; k<NCAPS; k++)
  {
    if (CAPS(k).isactive)
      generate_lag_stencils_one_caps(mesh=&CAPS(k), no_warning=no_warnings);
  }
}


/**
## Introduction
In this file, you can find the discrete delta functions used for the force spreading and velocity interpolation in the feedback Immersed Boundary Method (IBM). Here I provide four different types of delta function:

* Regular 2-point (IBM_stencil=1) and 4-point (IBM_stencil=4) delta functions
* Smoothed 2-point (IBM_stencil=11) and 4-point delta functions  (IBM_stencil=14) proposed by Yang et, al [\[1\]](#Yang2009).

## Variable list
* stencil: delta function value
* IBM_stencil: delta function type
* dist:    distance between the cell center and the Lagrangian node on boundary.
*/
#ifndef IBM_stencil
  #if STENCIL_TYPE == 3
    #define IBM_stencil (11)
  #else 
    #define IBM_stencil (14)
  #endif
#endif


double weight_stencil(double dist)
{

    double stencil = 0.;

    if (IBM_stencil == 1)
    {
        if (fabs(dist) <= 1.)
        {
            stencil = 1. - fabs(dist);
        }
        else
        {
            stencil = 0.;
        }   
    }
    else if (IBM_stencil == 4)
    {
        if (fabs(dist) < 1.)
        {
            stencil = 0.125 * (3. - 2 * fabs(dist) + sqrt(1. + 4 * fabs(dist) - 4 * sq(dist)));
        }
        else if (fabs(dist) <= 2. && fabs(dist) > 1.)
        {
            stencil = 0.125 * (5. - 2 * fabs(dist) - sqrt(-7. + 12 * fabs(dist) - 4 * sq(dist)));
        }
    }
    else if (IBM_stencil == 11)
    {
        if (fabs(dist) <= 0.5)
        {
            stencil = 3./4. - sq(dist);
        }
        else if (fabs(dist) > 0.5 && fabs(dist) <= 1.5)
        {
            stencil = 9./8. - 3.*fabs(dist)/2.+ sq(dist)/2.;
        }
        else
        {
            stencil = 0.;
        }
    }
    else if (IBM_stencil == 14)
    {
        if (fabs(dist) <= 0.5)
        {
            stencil = 3./8. + M_PI/32. - sq(dist)/4.;
        }
        else if (fabs(dist) > 0.5 && fabs(dist) <= 1.5)
        {
            stencil = 1./4. + (1-fabs(dist))*sqrt(-2.+8.*fabs(dist)-4*sq(dist))/8. - asin(sqrt(2.)*(fabs(dist)-1.))/8.;
        }
        else if (fabs(dist) > 1.5 && fabs(dist) <= 2.5)
        {
            stencil = 17./16.-M_PI/64.-3.*fabs(dist)/4.+sq(dist)/8.+(fabs(dist)-2.)*sqrt(-14.+16.*fabs(dist)-4.*sq(dist))/16.+asin(sqrt(2.)*(fabs(dist)-2.))/16.;
        }
        else
        {
            stencil = 0.;
        }
    }

    return stencil;
}

/**
The function below fills the scalar $forcing$ with the regularized source term
of magnitude $force$ located at $x0$. This assumes that $forcing$ is "clean",
i.e. that it is zero everywhere (or that if there is already some non-zero terms
the intention is to include them in the forcing).
*/
trace
void lag2eul(vector forcing, lagMesh* mesh) {
  for(int i=0; i<mesh->nln; i++) {
    foreach_cache(mesh->nodes[i].stencil) {
      if (point.level >= 0) {
        #if EMBED
          if (cs[] > 1.e-10) {
        #endif
        coord dist;
        dist.x = GENERAL_1DIST(x, mesh->nodes[i].pos.x, L0*L0_ratio.x);
        dist.y = GENERAL_1DIST(y, mesh->nodes[i].pos.y, L0*L0_ratio.y);
        #if dimension > 2
        dist.z = GENERAL_1DIST(z, mesh->nodes[i].pos.z, L0*L0_ratio.z);
        #endif
        #if dimension < 3
        if (fabs(dist.x) <= 2*Delta && fabs(dist.y) <= 2*Delta) {
          // double weight =
          //   (1 + cos(.5*pi*dist.x/Delta))*(1 + cos(.5*pi*dist.y/Delta))
          //   /(sq(4*Delta));
        double weight =  weight_stencil(dist.x/Delta) * weight_stencil(dist.y/Delta) / sq(Delta);
        #else
        if (fabs(dist.x) <= 2*Delta && fabs(dist.y) <= 2*Delta &&
          fabs(dist.z) <= 2*Delta) {
          // double weight =
          //   (1 + cos(.5*pi*dist.x/Delta))*(1 + cos(.5*pi*dist.y/Delta))
          //   *(1 + cos(.5*pi*dist.z/Delta))/(cube(4*Delta));
        double weight =  weight_stencil(dist.x/Delta) * weight_stencil(dist.y/Delta) * weight_stencil(dist.z/Delta) /cube(Delta);
        #endif
          foreach_dimension() forcing.x[] += weight*mesh->nodes[i].lagForce.x;
        }
        #if EMBED
          }
        #endif
      }
    }
  }
}

/**
The function below interpolates the eulerian velocities onto the nodes of
the Lagrangian mesh.
*/
trace
void eul2lag(lagMesh* mesh) {
  for(int ii=0; ii<mesh->nln; ii++) {
    foreach_dimension() mesh->nodes[ii].lagVel.x = 0.;
    foreach_cache(mesh->nodes[ii].stencil) {
      if (point.level >= 0) {
        coord dist;
        dist.x = GENERAL_1DIST(x, mesh->nodes[ii].pos.x, L0*L0_ratio.x);
        dist.y = GENERAL_1DIST(y, mesh->nodes[ii].pos.y, L0*L0_ratio.y);
        #if dimension > 2
        dist.z = GENERAL_1DIST(z, mesh->nodes[ii].pos.z, L0*L0_ratio.z);
        #endif
        #if dimension < 3
        if (fabs(dist.x) <= 2*Delta && fabs(dist.y) <= 2*Delta) {
          // double weight = (1 + cos(.5*pi*dist.x/Delta))*
          //   (1 + cos(.5*pi*dist.y/Delta))/16.;
        double weight =  weight_stencil(dist.x/Delta) * weight_stencil(dist.y/Delta);            
        #else
        if (fabs(dist.x) <= 2*Delta && fabs(dist.y) <= 2*Delta
          && fabs(dist.z) <= 2*Delta) {
          // double weight = (1 + cos(.5*pi*dist.x/Delta))*
          //   (1 + cos(.5*pi*dist.y/Delta))*(1 + cos(.5*pi*dist.z/Delta))/64.;
        double weight =  weight_stencil(dist.x/Delta) * weight_stencil(dist.y/Delta) * weight_stencil(dist.z/Delta);            
        #endif
        foreach_dimension() mesh->nodes[ii].lagVel.x += weight*u.x[];
        }
      }
    }
  }
}

/**
The function below fills a scalar field "stencils" with noise in all "cached"
cells. Passing this scalar to the \textit{adapt_wavelet} function ensure all
the 5x5(x5) stencils around the Lagrangian nodes are at the same level.
*/
scalar stencils[];
trace
void tag_ibm_stencils_one_caps(lagMesh* mesh) {
  for(int i=0; i<mesh->nln; i++) {
    foreach_cache(mesh->nodes[i].stencil) {
      if (point.level >= 0) {
        coord dist;
        dist.x = GENERAL_1DIST(x, mesh->nodes[i].pos.x, L0*L0_ratio.x);
        dist.y = GENERAL_1DIST(y, mesh->nodes[i].pos.y, L0*L0_ratio.y);
        #if dimension > 2
        dist.z = GENERAL_1DIST(z, mesh->nodes[i].pos.z, L0*L0_ratio.z);
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
}

trace
void tag_ibm_stencils() {
  foreach() stencils[] = 0.;
  for(int k=0; k<NCAPS; k++)
    if (CAPS(k).isactive) tag_ibm_stencils_one_caps(&CAPS(k));
}

/**
## References

~~~bib
@Article{Peskin1977,
  author    = {Peskin, C.S.},
  title     = {{Numerical analysis of blood flow in the heart}},
  journal   = {Journal of Computational Physics},
  year      = {1977},
  volume    = {25},
  number    = {3},
  pages     = {220--252},
  file      = {:files/Peskin1977 - Numerical Analysis of Blood Flow in the Heart.pdf:PDF},
  groups    = {FluidSolid flows},
  timestamp = {2013.07.18},
}
~~~
*/
