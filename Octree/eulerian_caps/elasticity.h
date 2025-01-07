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
tensor $\mathbf{B} = J^{-1}\mathbf{B}$, where
$J = \sqrt{(tr(\mathbf{B})^2 - tr(\mathbf{B}^2))/2}$ is the surface Jacobian
of the membrane, i.e. the local membrane area change [](#ii2012computational).

Both $J$ and $\mathbf{B}$ are only defined in the capsule region $\Gamma$, and
in this region they follow the advection equations:

$$ \partial_t J + \mathbf{u} \cdot \nabla J = (\nabla_s \cdot \mathbf{u}) J $$
$$ \partial_t \mathbf{B} + \mathbf{u} \cdot \nabla \mathbf{B} =
\mathbf{B} \cdot \nabla_s \mathbf{u} +  (\nabla_s \mathbf{u})^T \cdot \mathbf{B}
- (\nabla_s \cdot \mathbf{u}) \mathbf{B}
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

// /**
// We can also advect the Eulerian quantities in the sole region of the membrane
// (instead of going through all grid cells and computing zero fluxes when we are
// out of the membrane area). In that case, we only need one layer of initialized
// cells close the the membrane
// */
// #ifndef LOCAL_BCG
//   #define LOCAL_BCG 0
// #endif
//
// #if (LOCAL_BCG)
//   #include "local_bcg.h"
// #endif

#ifndef EXTEND_MB_ATTRIBUTES
  #define EXTEND_MB_ATTRIBUTES 1
#endif

/**
We now declare the centered Eulerian quantities: the modified left Cauchy-Green
deformation tensor $\mathbf{B}$ and its source term in its advection equation
$\mathbf{\text{Source}_B}$, the velocity gradient and surface gradient,
a vector defining normals to the interface in the whole membrane region $\Gamma$,
the surface Jacobian $J$ and its source term in its advection equation
$\text{Source}_J$, and finally the elastic stress tensor $\mathbf{tau}$.
*/

symmetric tensor B[], sB[];

event defaults (i = 0) {
  for (scalar s in {B, sB}) {
    s.v.x.i = -1;
    foreach_dimension() {
      s[left] = neumann(0.);
      s[right] = neumann(0.);
    }
  }
  foreach() {
    pseudo_t P;
    foreach_dimension() {
      P.x.x = 1 - sq(extended_n.x[]);
      P.x.y = -extended_n.x[]*extended_n.y[];
      P.x.z = -extended_n.x[]*extended_n.z[];
    }
    foreach_dimension() {
      B.x.x[] = P.x.x;
    }
    B.x.y[] = P.x.y;
    B.x.z[] = P.x.z;
    B.y.z[] = P.y.z;
  }
}

event init_again (i = 1) {
  foreach() {
    pseudo_t P;
    foreach_dimension() {
      P.x.x = 1 - sq(extended_n.x[]);
      P.x.y = -extended_n.x[]*extended_n.y[];
      P.x.z = -extended_n.x[]*extended_n.z[];
    }
    foreach_dimension() {
      B.x.x[] = P.x.x;
    }
    B.x.y[] = P.x.y;
    B.x.z[] = P.x.z;
    B.y.z[] = P.y.z;
  }
}

event init (i = 0) {
  boundary((scalar *) {B, sB});
}

event tracer_advection (i++) {
  foreach() {
    /**
    We loop through the membrane cells, and prepare construction of RHS of
    the advection equations of J and B:
    */
    foreach_dimension() {
      sB.x.x[] = 0.;
    }
    sB.x.y[] = 0.;
    sB.x.z[] = 0.;
    sB.y.z[] = 0.;

    #if EXTEND_MB_ATTRIBUTES
      if (IS_INTERFACE_CELL(point,f)) {
    #else
      if (GAMMA) {
    #endif
      /**
      We construct the velocity gradient tensor at the center of the cells.
      */
      pseudo_t grad_u;
      foreach_dimension() {
        grad_u.x.x = (u.x[1,0,0] - u.x[-1,0,0])/(2.*Delta);
        grad_u.x.y = (u.x[0,1,0] - u.x[0,-1,0])/(2.*Delta);
        grad_u.x.z = (u.x[0,0,1] - u.x[0,0,-1])/(2.*Delta);
      }
      foreach_dimension() {
        sB.x.x[] = 2*(grad_u.x.x*B.x.x[] + grad_u.x.y*B.y.x[]
          + grad_u.x.z*B.z.x[]);
      }
      sB.x.y[] = grad_u.x.x*B.x.y[] + grad_u.x.y*B.y.y[] + grad_u.x.z*B.z.y[]
        + grad_u.y.x*B.x.x[] + grad_u.y.y*B.x.y[] + grad_u.y.z*B.x.z[];
      sB.x.z[] = grad_u.x.x*B.x.z[] + grad_u.x.y*B.y.z[] + grad_u.x.z*B.z.z[]
        + B.x.x[]*grad_u.z.x + B.x.y[]*grad_u.z.y + B.x.z[]*grad_u.z.z;
      sB.y.z[] = grad_u.z.x*B.x.z[] + grad_u.z.y*B.y.x[] + grad_u.z.z*B.z.z[]
        + B.z.x[]*grad_u.z.x + B.z.y[]*grad_u.z.y + B.z.z[]*grad_u.z.z;
    } // end if the cell is in the membrane region
  } // foreach

  /** We now use the Bell-Colella-Glaz scheme to compute the material
  derivative */
  boundary((scalar *) {B, sB});
  advection((scalar *) {B}, uf, dt);

  foreach() {
    #if EXTEND_MB_ATTRIBUTES
      if (IS_INTERFACE_CELL(point,f)) {
    #else
      if (GAMMA) {
    #endif
      /** We advance in time J and \mathbf{B} using their respective
      source terms*/
      foreach_dimension() {
        B.x.x[] += dt*sB.x.x[];
      }
      B.x.y[] += dt*sB.x.y[];
      B.x.z[] += dt*sB.x.z[];
      B.y.z[] += dt*sB.y.z[];
    }
    if (!GAMMA) {
      foreach_dimension() {
        B.x.x[] = 1.;
      }
      B.x.y[] = 0.;
      B.x.z[] = 0.;
      B.y.z[] = 0.;
    }
  }
  boundary((scalar *){B});

  /**Then we extend J and \mathbf{B} in the normal direction from the interface
  in order to force these quantities to be constant on the interface. In this
  process, we re-use sJ and sB as temporary field values.*/
  #if EXTEND_MB_ATTRIBUTES
    normal_scalar_extension((scalar *) {B}, (scalar *) {sB}, nb_iter_extension);
  #endif
} // tracer_advection
