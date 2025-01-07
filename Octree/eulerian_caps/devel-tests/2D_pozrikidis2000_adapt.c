// #define JACOBI 1
#define EPS 1.e-14
#ifndef SIZE
  #define SIZE 4.
#endif
#define SHEAR_RATE 1.
#define RADIUS 1.
#define MU 1.
#ifndef CA
  #define CA 0.003125
#endif
#define G_SHEAR (MU*RADIUS*SHEAR_RATE/CA)
#define RHO 1.
#define T_BC (1./SHEAR_RATE)
#ifndef TEND
  #define TEND (10./SHEAR_RATE)
#endif
#define AREA_DILATATION_MODULUS 1.
#ifndef RESTRICT_DT
  #define RESTRICT_DT .05
#endif
#ifndef LEVEL
  #define LEVEL 6
#endif
#ifndef U_TOL
  #define UTOL 1.e-2
#endif
#ifndef START_SHEAR
  #define START_SHEAR 0
#endif

#define LOCAL_BCG 0
#ifndef NB_MB_EXTRA_LAYERS
  #define NB_MB_EXTRA_LAYERS 2
#endif
#define PROJECT_G 1
#ifndef OUTPUT_FREQ
  #define OUTPUT_FREQ 100
#endif
#ifndef NORMAL_EXTENSION
  #define NORMAL_EXTENSION 0
#endif

// #include "grid/multigrid.h"
#include "eulerian_caps/navier-stokes/my_centered.h"
#include "eulerian_caps/capsule.h"
#include "eulerian_caps/elasticity_2D.h"
#include "eulerian_caps/neo_hookean_2D.h"
#include "view.h"
#include "navier-stokes/perfs.h"

FILE * fp = NULL;
scalar sdivu[];
scalar noise_gamma[];
scalar my_gamma[];
// int LEVEL;

int main() {
  L0 = SIZE;
  origin (-0.5*L0, -0.5*L0);

  stokes = true;
  TOLERANCE = 1.e-10;
  rho1 = RHO;
  rho2 = RHO;
  mu1 = MU;
  mu2 = MU;

  N = 1 << LEVEL;

  periodic(right);

  run();
}

u.n[bottom] = dirichlet(0);
u.n[top] = dirichlet(0);
// u.t[bottom] = dirichlet((t < T_BC ? SHEAR_RATE*y : 0.));
// u.t[top] = dirichlet((t < T_BC ? SHEAR_RATE*y : 0.));
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

event adapt (i++) {
  foreach() {
    if (GAMMA) {noise_gamma[] = noise(); my_gamma[] = 1.;}
    else {noise_gamma[] = 0.; my_gamma[] = 0.;}
  }
  boundary({noise_gamma});
  adapt_wavelet ({noise_gamma,u}, (double []){1.e-6,U_TOL, U_TOL}, maxlevel = LEVEL);
}

#if TRACE > 1
event profiling (i += 20)
{
  static FILE * fp = fopen ("profiling", "a"); // In case of restart
  trace_print (fp, 1); // Display functions taking more than 1% of runtime.
}
#endif // TRACE

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

event movie (i += OUTPUT_FREQ) {
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
  squares ("noise_gamma", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("noise_gamma.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("my_gamma", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("gamma.mp4");
}

event end (t = TEND) {
  return 1.;
}
