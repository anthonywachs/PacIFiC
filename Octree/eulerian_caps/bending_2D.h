/**
# Bending stress for the Eulerian capsule framework.

In this file we implement the bending stress introduced by Pozrikidis (2001). The computation of the stresses rely on the Eulerian quantities defined in [elasticity.h](elasticity.h). The force are transferred to the fluid in the whole membrane region computed in [capsule.h](capsule.h).
*/

#ifndef E_BEND
  #define E_BEND 1.
#endif
#ifndef REF_CURV
  #define REF_CURV 0
#endif


typedef struct { double x, y;}   pseudo_v;

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

#if (REF_CURV)
  tensor m_ref[];
  event init (i = 0) {
    foreach() {
      foreach_dimension() {
        m_ref.x.x[] = 0.;
        m_ref.x.y[] = 0.;
      }
    }
    boundary((scalar *) {m_ref});
  }

  event init_pp_vof (i = 1) {
    foreach() {
      if (GAMMA) {
        foreach_dimension() {
          m_ref.x.x[] = my_gradient_x(point, extended_n.x);
          m_ref.x.y[] = my_gradient_y(point, extended_n.x);
        }
      }
      else {
        foreach_dimension() {
          m_ref.x.x[] = 0.;
          m_ref.x.y[] = 0.;
        }
      }
    }
    boundary((scalar *) {m_ref});
    normal_scalar_extension((scalar *) {m_ref}, (scalar *) {ext_buf});
  }
#endif

event acceleration (i++)
{
  if (E_BEND > 1.e-22) {
    #if (REF_CURV)
    foreach() {
      if (GAMMA) {
        foreach_dimension() {
          m_bend.x.x[] = E_BEND*(my_gradient_x(point, extended_n.x) - m_ref.x.x[]);
          m_bend.x.y[] = E_BEND*(my_gradient_y(point, extended_n.x) - m_ref.x.y[]);
        }
      }
      else {
        foreach_dimension() {
          m_bend.x.x[] = 0.;
          m_bend.x.y[] = 0.;
        }
      }
    }
    #else
      foreach() {
        if (GAMMA) {
          foreach_dimension() {
            m_bend.x.x[] = E_BEND*my_gradient_x(point, extended_n.x);
            m_bend.x.y[] = E_BEND*my_gradient_y(point, extended_n.x);
          }
        }
        else {
          foreach_dimension() {
            m_bend.x.x[] = 0.;
            m_bend.x.y[] = 0.;
          }
        }
      }
    #endif
    boundary((scalar *) {m_bend});
    normal_scalar_extension((scalar *) {m_bend}, (scalar *) {ext_buf});

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
      else {
        foreach_dimension() {
          q_bend.x[] = 0.;
        }
      }
    }
    boundary((scalar *) {q_bend});
    normal_scalar_extension((scalar *) {q_bend}, (scalar *) {ext_buf});

    foreach() {
      if (GAMMA) {
        foreach_dimension()
        {
          T_s.x.x[] += q_bend.x[]*extended_n.x[];
          T_s.x.y[] += q_bend.x[]*extended_n.y[];
        }
      }
      else {
        foreach_dimension() {
          T_s.x.x[] = 0.;
          T_s.x.y[] = 0.;
        }
      }
    }
    boundary((scalar *) {T_s});
    normal_scalar_extension((scalar *) {T_s}, (scalar *) {ext_buf});
  }
}

#if (REF_CURV)
  event tracer_advection (i++) {
    advection({m_ref.x.x, m_ref.x.y, m_ref.y.x, m_ref.y.y}, uf, dt);
    boundary((scalar *) {m_ref});
    normal_scalar_extension((scalar *) {m_ref}, (scalar *) {ext_buf});
  }
#endif

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
