/**
# Neo-hookean law for the Eulerian elasticity framework.

In this file we implement the neo-hookean law to compute the elastic surface stresses
of biological capsules. The computation of the stresses rely on the Eulerian
quantities defined in [elasticity.h](elasticity.h). The force are transferred
to the fluid in the whole membrane region computed in [capsule.h](capsule.h).
*/

event acceleration (i++)
{
  foreach()
  {
    if (GAMMA)
    {
/**
Below, we compute the stress tensor for the neo-Hookean law:
$$
\mathbf{T_s} = \frac{G_s}{3} (\mathbf{G} - \frac{1}{J^3}\mathbf{P})
$$
Note that the elastic modulus $G_s$ will be multiplied to the stress tensor
during the acceleration step in [elasticity.h](http://basilisk.fr/_edit/sandbox/huet/src/elasticity.h),
and that the acceleration will be computed using the continuum surface force
(CSF) formulation, for which we need to multiply the divergence of the stress
tensor $\nabla_s \cdot \mathbf{T_s}$ by the norm of the gradient of the smoothed color
function $|\nabla \phi|$ (in our case, the smoothed color function is defined in [capsule.h](http://basilisk.fr/_edit/sandbox/huet/src/capsule.h) and named *caps* $\equiv \phi$).
It was shown in [Ii et al.](ii2012full) that since
$|\nabla \phi|$ is a smoothed 1-dimensional Dirac distribution, we have
$|\nabla \phi| \nabla_s \cdot \mathbf{T_s} = \nabla_s \cdot (|\nabla \phi| \mathbf{T_s})$.
As a result, we multiply the stress tensor by $|\nabla \phi|$ at this
step.
*/
      foreach_dimension()
      {
        T_s.x.x[] += (G_SHEAR/3.)*(1. - 1./(cube(J[])))*
                    (sq(extended_n.y[]));
        T_s.x.y[] += -(G_SHEAR/3.)*(1. - 1./(cube(J[])))*
                    (extended_n.x[]*extended_n.y[]);
      }
    }
  }
  boundary((scalar *){T_s});
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
