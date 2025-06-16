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

#include "grid/multigrid.h"
#include "eulerian_caps/navier-stokes/my_centered.h"
#include "eulerian_caps/capsule.h"
#include "eulerian_caps/elasticity_2D.h"
#include "eulerian_caps/neo_hookean_2D.h"
// #include "eulerian_caps/neo_hookean.h"
#include "view.h"

FILE * fp = NULL;
scalar sdivu[];
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

event logfile (i+=OUTPUT_FREQ) {
  double rmax = -HUGE, rmin = HUGE ;
  foreach (reduction(max:rmax) reduction(min:rmin))
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
  double D = (rmax - rmin)/(rmax + rmin);
  fprintf (stderr, "%g %g %g %g\n", t, rmin, rmax, D);
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
}

event end (t = TEND) {
  return 1.;
}
