/**
# Eulerian elasticity

In this file we define and advect the tensors that track the elongation of a
region of space $\Gamma$ occupied by an elastic solid. It is meant to be
combined with [capsule.h](capsule.h) to provide a description of biological
membranes.

It was shown by [Barthès-Biesel and Rallison](#berthes1981time) that the stress
tensor of an elastic membrane can be described in an Eulerian fashion with the
left Cauchy-Green deformation tensor $\bm{B}$ and the local projector onto the
membrane $\mathbf{P} = \mathbf{I} - \mathbf{n}\mathbf{n}$, with
$\mathbf{n}$ denoting the interface unit normal vector.

Recently, it was demonstrated by [Ii et al.](#ii2012full) that this
Eulerian expression of the stress tensor can be successfully implemented in order to
simulate elastic capsules without meshing their membrane and without requiring a
boundary-fitted fluid mesh or an immersed boundary method (IBM). They report
that it is more stable to advect the *modified* left Cauchy-Green deformation
tensor $\mathbf{G} = J^{-1}\mathbf{B}$, where
$J = \sqrt{(tr(\mathbf{B})^2 - tr(\mathbf{B}^2))/2}$ is the surface Jacobian
of the membrane, i.e. the local membrane area change [](#ii2012computational).

Both $J$ and $\mathbf{G}$ are only defined in the capsule region $\Gamma$, and
in this region they follow the advection equations:

$$ \partial_t J + \mathbf{u} \cdot \nabla J = (\nabla_s \cdot \mathbf{u}) J $$
$$ \partial_t \mathbf{G} + \mathbf{u} \cdot \nabla \mathbf{G} =
\mathbf{G} \cdot \nabla_s \mathbf{u} +  (\nabla_s \mathbf{u})^T \cdot \mathbf{G}
- (\nabla_s \cdot \mathbf{u}) \mathbf{G}
$$
$$\text{with} \quad \nabla_s = \mathbf{P} \cdot \nabla$$
denoting by $\mathbf{n}$ the interface normal, $\nabla_s$ the surface gradient and
$\nabla_s \cdot$ the surface divergence.
*/

/**
We have the option of defining a pre-inflated capsule by setting $J_0$ to a
value greater than 1. By default, there is no stress on the capsule in the
reference configuration.
*/
#ifndef J0
  #define J0 1.
#endif

/**
When the membrane region moves and a cell becomes part of the membrane, it needs
to have some initialized values for the Eulerian quantities $\mathbf{G}$ and
$J$. Therefore we create a layer of membrane cells that will be initialized to
the average values of its neighbouring cells already belonging to the membrane.
When the default BCG advection scheme is used, we need this layer of initialized
cells to have a thickness of $2 \Delta$ for a proper computation of fluxes.
*/
#ifndef NB_MB_EXTRA_LAYERS
  #define NB_MB_EXTRA_LAYERS 2
#endif

/**
We can also advect the Eulerian quantities in the sole region of the membrane
(instead of going through all grid cells and computing zero fluxes when we are
out of the membrane area). In that case, we only need one layer of initialized
cells close the the membrane
*/
#ifndef LOCAL_BCG
  #define LOCAL_BCG 0
  #ifndef NB_MB_EXTRA_LAYERS
    #define NB_MB_EXTRA_LAYERS 1
  #endif
#endif

#if (LOCAL_BCG)
  #include "local_bcg.h"
#endif

/**
We also define $G_s$, the elastic modulus of the material...
*/
#ifndef G_SHEAR
  #define G_SHEAR 2.5e-6
#endif

/**
We now declare the centered Eulerian quantities: the modified left Cauchy-Green
deformation tensor $\mathbf{G}$ and its source term in its advection equation
$\mathbf{\text{Source}_G}$, the velocity gradient and surface gradient,
a vector defining normals to the interface in the whole membrane region $\Gamma$,
the surface Jacobian $J$ and its source term in its advection equation
$\text{Source}_J$, and finally the elastic stress tensor $\mathbf{T_s}$.
*/

scalar J[], source_J[];
tensor ext_buf[];

event defaults (i = 0) {
  foreach() {
    J[] = J0;
  }
}

event init (i = 0) {
  boundary({J, source_J});
}

/**
--- OBSOLETE (keeping for backward compatibility) ---
The stability event below is a pure copy-paste from
[tension.h](http://basilisk.fr/src/tension.h), with $\sigma$
replaced by *G_SHEAR*.
*/
// #ifndef RESTRICT_DT
//   #define RESTRICT_DT 1.
// #endif
//
// event stability (i++) {
//   double amin = HUGE, amax = -HUGE, dmin = HUGE;
//   foreach_face (reduction(min:amin) reduction(max:amax) reduction(min:dmin)) {
//     if (alpha.x[]/fm.x[] > amax) amax = alpha.x[]/fm.x[];
//     if (alpha.x[]/fm.x[] < amin) amin = alpha.x[]/fm.x[];
//     if (Delta < dmin) dmin = Delta;
//   }
//   double rhom = (1./amin + 1./amax)/2.;
//   double dt = RESTRICT_DT*sqrt (rhom*cube(dmin)/(pi*G_SHEAR));
//   if (dt < dtmax)
//   dtmax = dt;
// }
// --- END OF OBSOLETE CODE --- //


event tracer_advection (i++) {
  foreach() {
    /**
    We loop through the membrane cells, and prepare construction of RHS of
    the advection equations of J and G:
    */
    source_J[] = 0.;

    if (GAMMA) {
      /**
      We construct the velocity gradient tensor at the center of the cells.
      */
      foreach_dimension() {
        my_grad_u.x.x[] = (u.x[1,0] - u.x[-1,0])/(2.*Delta);
        my_grad_u.x.y[] = (u.x[0,1] - u.x[0,-1])/(2.*Delta);
      }

      foreach_dimension() {
      /**
      We then construct velocity surface gradient and divergence
      at the center of the cells.
      */
        sgrad_u.x.x[] = (1 - sq(extended_n.x[]))*my_grad_u.x.x[] -
                        extended_n.x[]*extended_n.y[]*my_grad_u.x.y[];
        sgrad_u.x.y[] = (1 - sq(extended_n.y[]))*my_grad_u.x.y[] -
                        extended_n.x[]*extended_n.y[]*my_grad_u.x.x[];
        source_J[] += sgrad_u.x.x[];
      }
      source_J[] *= J[];
    } // end if the cell is in the membrane region
  } // foreach

  /**
  To advect J and G, we use bcg.h. The right hand side of the advection equation
  is added later.
  */
  boundary({J, source_J});
  #if (LOCAL_BCG)
    local_advection({J}, uf, dt);
  #else
    advection({J}, uf, dt);
  #endif

  /**
  Now we consider the right-hand-side of the advection equation, that we simply
  add to each advected quantity in a forward-Euler fashion.
  */
  foreach() {
    if (GAMMA) {
      J[] += dt*source_J[];
    }
    else {
      J[] = J0;
      foreach_dimension() {
        my_grad_u.x.x[] = 0.;
        my_grad_u.x.y[] = 0.;
        sgrad_u.x.x[] = 0.;
        sgrad_u.x.y[] = 0.;
      }
    }
  }
  boundary({J});
  normal_scalar_extension({J}, (scalar *) {ext_buf});
} // tracer_advection

/**
# References
~~~bib
@Article{ii2018continuum,
  author    = {Ii, Satoshi and Shimizu, Kazuya and Sugiyama, Kazuyasu and Takagi, Shu},
  journal   = {Journal of Computational Physics},
  title     = {Continuum and stochastic approach for cell adhesion process based on Eulerian fluid-capsule coupling with Lagrangian markers},
  year      = {2018},
  pages     = {769--786},
  volume    = {374},
  file      = {:ii2018continuum - Continuum and Stochastic Approach for Cell Adhesion Process Based on Eulerian Fluid Capsule Coupling with Lagrangian Markers.pdf:PDF},
  groups    = {Biological flows},
  publisher = {Elsevier},
}

@Article{ii2012computational,
  author    = {Ii, Satoshi and Sugiyama, Kazuyasu and Takagi, Shu and Matsumoto, Yoichiro},
  journal   = {Journal of Biomechanical Science and Engineering},
  title     = {A computational blood flow analysis in a capillary vessel including multiple red blood cells and platelets},
  year      = {2012},
  number    = {1},
  pages     = {72--83},
  volume    = {7},
  file      = {:ii2012computational - A Computational Blood Flow Analysis in a Capillary Vessel Including Multiple Red Blood Cells and Platelets.pdf:PDF},
  groups    = {Biological flows},
  publisher = {The Japan Society of Mechanical Engineers},
}

@Article{ii2012full,
  author  = {Ii, Satoshi and Gong, Xiaobo and Sugiyama, Kazuyasu and Wu, Jinbiao and Huang, Huaxiong and Takagi, Shu},
  journal = {Communications in Computational Physics},
  title   = {A full Eulerian fluid-membrane coupling method with a smoothed volume-of-fluid approach},
  year    = {2012},
  number  = {2},
  pages   = {544},
  volume  = {12},
  file    = {:ii2012full - A Full Eulerian Fluid Membrane Coupling Method with a Smoothed Volume of Fluid Approach.pdf:PDF},
  groups  = {Biological flows},
}

@Article{barthes1981time,
  author    = {Barthes-Biesel, Dominique and Rallison, JM},
  journal   = {Journal of Fluid Mechanics},
  title     = {The time-dependent deformation of a capsule freely suspended in a linear shear flow},
  year      = {1981},
  pages     = {251--267},
  volume    = {113},
  file      = {:barthes1981time - The Time Dependent Deformation of a Capsule Freely Suspended in a Linear Shear Flow.pdf:PDF},
  groups    = {Biological flows},
  publisher = {Cambridge University Press},
}
~~~
*/
