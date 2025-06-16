/**
# Bending stress for the Eulerian capsule framework.

In this file we implement the bending stress introduced by Pozrikidis (2001). The computation of the stresses rely on the Eulerian quantities defined in [elasticity.h](elasticity.h). The force are transferred to the fluid in the whole membrane region computed in [capsule.h](capsule.h).
*/

#ifndef E_BEND
  #define E_BEND (7.e-19)
#endif

typedef struct { double x, y;}   pseudo_v;

/**
The stability event below is a pure copy-paste from
[tension.h](http://basilisk.fr/src/tension.h), with $\sigma$
replaced by *E_BEND*.
*/
event stability (i++) {
  double amin = HUGE, amax = -HUGE, dmin = HUGE;
  foreach_face (reduction(min:amin) reduction(max:amax) reduction(min:dmin)) {
    if (alpha.x[]/fm.x[] > amax) amax = alpha.x[]/fm.x[];
    if (alpha.x[]/fm.x[] < amin) amin = alpha.x[]/fm.x[];
    if (Delta < dmin) dmin = Delta;
  }
  double rhom = (1./amin + 1./amax)/2.;
  double dt = RESTRICT_DT*sqrt (rhom*cube(dmin)/(pi*E_BEND));
  if (dt < dtmax)
  dtmax = dt;
}


foreach_dimension() {
  double my_gradient_x( Point point, scalar s) {
    if (BIGAMMA(-1,0) && BIGAMMA(1,0)) return ((s[1,0] - s[-1,0])/(2.*Delta));
    else if (BIGAMMA(1,0)) {
      return ((-s[2,0] + 4*s[1,0] - 3*s[])/(2.*Delta));
    }
    else {
      return ((s[-2,0] - 4*s[-1,0] + 3*s[])/(2.*Delta));
    }
  }
}

tensor m_bend[];
vector q_bend[];

event acceleration (i++)
{
  foreach() {
    if (GAMMA) {
      // Step 1: compute the tensor of bending moments m_bend = grad(extended_n);
      foreach_dimension() {
        m_bend.x.x[] = E_BEND*my_gradient_x(point, extended_n.x);
        m_bend.x.y[] = E_BEND*my_gradient_y(point, extended_n.x);
      }
    }
  }
  boundary((scalar *) {m_bend});
  foreach() {
    pseudo_v sdiv_m;
    double dxmxx,dymxx, dxmyx, dymyx;
    if (GAMMA) {
    // Step 2: compute the transverse shear vector q_bend = div_s(m_bend) * P
      foreach_dimension() {
        dxmxx = my_gradient_x(point, m_bend.x.x);
        dymxx = my_gradient_y(point, m_bend.x.x);
        dxmyx = my_gradient_x(point, m_bend.y.x);
        dymyx = my_gradient_y(point, m_bend.y.x);
        sdiv_m.x = sq(extended_n.y[])*dxmxx - extended_n.x[]*extended_n.y[]*(dymxx + dxmyx) + sq(extended_n.x[])*dymyx;
      }
      foreach_dimension() {
        q_bend.x[] = sq(extended_n.y[])*sdiv_m.x - extended_n.x[]*extended_n.y[]*sdiv_m.y;
      }
    }
  }
  boundary((scalar *) {q_bend});

  foreach() {
    if (GAMMA)
    {
      foreach_dimension()
      {
        T_s.x.x[] += q_bend.x[]*extended_n.x[];
        T_s.x.y[] += q_bend.x[]*extended_n.y[];
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
