/**
# Neo-hookean law for the Eulerian elasticity framework.

In this file we implement the neo-hookean law to compute the elastic surface stresses
of biological capsules. The computation of the stresses rely on the Eulerian
quantities defined in [elasticity.h](elasticity.h). The force are transferred
to the fluid in the whole membrane region computed in [capsule.h](capsule.h).
*/

/**
We define $B_s$, the elastic modulus of the material
*/
#ifndef E_S
  #define E_S 1.
#endif

event acceleration (i++) {
  foreach() {
    if (GAMMA) {
/**
Below, we compute the stress tensor for the neo-Hookean law:
$$
\mathbf{Ts} = \frac{E_s}{3 J} (\mathbf{B} - \frac{1}{J^2}\mathbf{P})
$$
Note that the elastic modulus $E_s$ will be multiplied to the stress tensor
during the acceleration step in [elasticity.h](http://basilisk.fr/_edit/sandbox/huet/src/elasticity.h),
and that the acceleration will be computed using the continuum surface force
(CSF) formulation, for which we need to multiply the divergence of the stress
tensor $\nabla_s \cdot \mathbf{Ts}$ by the norm of the gradient of the smoothed color
function $|\nabla \phi|$ (in our case, the smoothed color function is defined in [capsule.h](http://basilisk.fr/_edit/sandbox/huet/src/capsule.h) and named *caps* $\equiv \phi$).
It was shown in [Ii et al.](ii2012full) that since
$|\nabla \phi|$ is a smoothed 1-dimensional Dirac distribution, we have
$|\nabla \phi| \nabla_s \cdot \mathbf{Ts} = \nabla_s \cdot (|\nabla \phi| \mathbf{Ts})$.
As a result, we multiply the stress tensor by $|\nabla \phi|$ at this
step.
*/
      pseudo_t P, BP, PBP;
      foreach_dimension() {
        P.x.x = 1 - sq(extended_n.x[]);
        P.x.y = -extended_n.x[]*extended_n.y[];
        P.x.z = -extended_n.x[]*extended_n.z[];
      }
      foreach_dimension() {
        BP.x.x = B.x.x[]*P.x.x + B.x.y[]*P.y.x + B.x.z[]*P.z.x;
        BP.x.y = B.x.x[]*P.x.y + B.x.y[]*P.y.y + B.x.z[]*P.z.y;
        BP.x.z = B.x.x[]*P.x.z + B.x.y[]*P.y.z + B.x.z[]*P.z.z;
      }

      foreach_dimension() {
        PBP.x.x = P.x.x*BP.x.x + P.x.y*BP.y.x + P.x.z*BP.z.x;
        PBP.x.y = P.x.x*BP.x.y + P.x.y*BP.y.y + P.x.z*BP.z.y;
        PBP.x.z = P.x.x*BP.x.z + P.x.y*BP.y.z + P.x.z*BP.z.z;
      }

      double SqTrB = sq(PBP.x.x + PBP.y.y + PBP.z.z);
      double TrSqB = 0.;
      foreach_dimension() {
        TrSqB += sq(PBP.x.x) + sq(PBP.x.y) + sq(PBP.x.z);
      }
      double J = sqrt(.5*(SqTrB - TrSqB));

      foreach_dimension() {
        Ts.x.x[] = E_S*(PBP.x.x/J - P.x.x/cube(J))/3.;
      }
      Ts.x.y[] = E_S*(PBP.x.y/J - P.x.y/cube(J))/3.;
      Ts.x.z[] = E_S*(PBP.x.z/J - P.x.z/cube(J))/3.;
      Ts.y.z[] = E_S*(PBP.y.z/J - P.y.z/cube(J))/3.;
    }
    else {
      foreach_dimension() {
        Ts.x.x[] = 0.;
      }
      Ts.x.y[] = 0.;
      Ts.x.z[] = 0.;
      Ts.y.z[] = 0.;
    }
  }
  boundary((scalar *){Ts});
}

/**
# References
~~~bib
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
~~~
*/
