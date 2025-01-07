#ifndef LEVEL
  #define LEVEL 7
#endif
#ifndef SIZE
  #define SIZE 8
#endif
#ifndef RADIUS
  #define RADIUS 1.
#endif
#ifndef SHEAR_RATE
  #define SHEAR_RATE 1.
#endif
#ifndef RHO
  #define RHO 1.
#endif
#ifndef MU
  #define MU 1.
#endif
#ifndef CA
  #define CA 0.05
#endif
#ifndef G_SHEAR
  #define G_SHEAR (MU*SHEAR_RATE*RADIUS/CA)
#endif
#ifndef ND_EB
  #define ND_EB 0.
#endif
#ifndef E_BEND
  #define E_BEND (ND_EB*RADIUS*RADIUS*G_SHEAR)
#endif
#ifndef T_END
  #define T_END 4.
#endif
#ifndef START_SHEAR
  #define START_SHEAR 1
#endif
#ifndef DT_MAX
  #define DT_MAX (6.25e-5)
#endif
#ifndef MY_TOLERANCE
  #define MY_TOLERANCE 1.e-10
#endif
#ifndef U_TOL
  #define U_TOL (1.e-2)
#endif
#ifndef OUTPUT_FREQ
  #define OUTPUT_FREQ (t += 0.01)
#endif
#ifndef AMR
  #define AMR 0
#endif


#if (AMR)
  #include "grid/quadtree.h"
  #ifndef U_TOL
    #define U_TOL 1.e-2
  #endif
#else
  #include "grid/multigrid.h"
#endif
#include "eulerian_caps/navier-stokes/my_centered.h"
#include "eulerian_caps/capsule.h"
#include "eulerian_caps/elasticity_2D.h"
#include "eulerian_caps/bending_2D.h"
#include "eulerian_caps/neo_hookean_2D.h"
#include "view.h"

FILE * fp = NULL;
scalar noise_gamma[], sdivu[];

int main() {
  L0 = SIZE;
  origin(-.5*L0, -.5*L0);
  stokes = true;
  rho1 = RHO;
  rho2 = RHO;
  mu1 = MU;
  mu2 = MU;
  N = 1 << LEVEL;
  TOLERANCE = MY_TOLERANCE;
  DT = DT_MAX;
  periodic(right);
  run();
}

u.n[bottom] = dirichlet(0);
u.n[top] = dirichlet(0);
u.t[bottom] = dirichlet(SHEAR_RATE*y);
u.t[top] = dirichlet(SHEAR_RATE*y);

event init (t=0) {
  if (START_SHEAR) {
    foreach() {
      u.x[] = SHEAR_RATE*y;
      u.y[] = 0.;
    }
  }
  fraction(f, 1. - sq(x/(RADIUS)) - sq(y/(RADIUS)));
}

#if (AMR)
  event adapt (i++) {
    foreach() {
      if (GAMMA) noise_gamma[] = noise();
      else noise_gamma[] = 0.;
    }
    boundary({noise_gamma});
    adapt_wavelet ({noise_gamma,u}, (double []){1.e-6,U_TOL, U_TOL}, maxlevel = LEVEL);
  }
#endif

event logfile (i+=OUTPUT_FREQ) {
  scalar un[];
  int nb_mb_cells = 0;
  double sdivu_avg = 0., J_avg = 0., T_s_avg_xx = 0., T_s_avg_xy = 0.,
    T_s_avg_yy = 0., aen_avg = 0., abs_J_avg = 0., abs_T_s_avg_xx = 0.,
    abs_T_s_avg_xy = 0., abs_T_s_avg_yy = 0.;
  double rmax = -HUGE, rmin = HUGE ;
  foreach() {
    un[] = norm(u);
    if (GAMMA) {
      nb_mb_cells++;
      sdivu[] = sgrad_u.x.x[] + sgrad_u.y.y[];
      sdivu_avg += sdivu[];
      J_avg += J[];
      T_s_avg_xx += T_s.x.x[];
      T_s_avg_xy += T_s.x.y[];
      T_s_avg_yy += T_s.y.y[];
      aen_avg += aen[];
      abs_J_avg += fabs(J[]-1);
      abs_T_s_avg_xx += fabs(T_s.x.x[]);
      abs_T_s_avg_xy += fabs(T_s.x.y[]);
      abs_T_s_avg_yy += fabs(T_s.y.y[]);
    }
    else {
      sdivu[] = 0.;
    }
    if (f[] > 0 && f[] < 1) {
      coord p;
      coord n = mycs (point, f);
      double alpha = plane_alpha (f[], n);
      plane_area_center (n, alpha, &p);
      double rad  = sqrt(sq(x + Delta*p.x) + sq(y + Delta*p.y));
      if (rad > rmax)
      rmax = rad;
      if (rad < rmin)
      rmin = rad;
    }
  }
  double D = (rmax - rmin)/(rmax + rmin);
  if (nb_mb_cells != 0) {
    sdivu_avg /= nb_mb_cells;
    J_avg /= nb_mb_cells;
    T_s_avg_xx /= nb_mb_cells;
    T_s_avg_xy /= nb_mb_cells;
    T_s_avg_yy /= nb_mb_cells;
    aen_avg /= nb_mb_cells;
    abs_J_avg /= nb_mb_cells;
    abs_T_s_avg_xx /= nb_mb_cells;
    abs_T_s_avg_xy /= nb_mb_cells;
    abs_T_s_avg_yy /= nb_mb_cells;
  }
  fprintf (stdout, "%.17g %.17g %.17g %.17g %.17g %.17g %.17g "
    "%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g "
    "%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
    t, normf(un).max, normf(un).avg, normf(u.x).max, normf(u.x).avg,
    normf(u.y).max, normf(u.y).avg,
    normf(J).max, J_avg,
    normf(T_s.x.x).max, T_s_avg_xx, normf(T_s.x.y).max,
    T_s_avg_xy, normf(T_s.y.y).max, T_s_avg_yy,
    normf(aen).max, aen_avg, abs_J_avg, abs_T_s_avg_xx, abs_T_s_avg_xy, abs_T_s_avg_yy, sdivu_avg, D);
}

event movie (OUTPUT_FREQ) {
  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("u.x", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("ux.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("J", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("J.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("m_bend.x.x", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("m_bend_xx.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("q_bend.x", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("q_bend_x.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("T_s.x.x", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("Ts_xx.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("centered_ae.x", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("ax.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("aen", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("an.mp4");
}

event end (t = T_END) {
  return 1.;
}
