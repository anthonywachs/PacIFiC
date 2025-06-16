#ifndef LEVEL
  #define LEVEL 7
#endif
#ifndef SIZE
  #define SIZE 8
#endif
#ifndef RADIUS
  #define RADIUS 1.
#endif
#ifndef RHO
  #define RHO 1.
#endif
#ifndef MU
  #define MU 1.
#endif
#ifndef MY_TOLERANCE
  #define MY_TOLERANCE 1.e-10
#endif
#ifndef AMR
  #define AMR 1
#endif

#if AMR
  #include "grid/octree.h"
  #ifndef U_TOL
    #define U_TOL (0.01)
  #endif
#else
  #include "grid/multigrid.h"
#endif
#include "eulerian_caps/navier-stokes/my_centered.h"
#include "eulerian_caps/capsule.h"
#include "view.h"

FILE * fp = NULL;
scalar noise_gamma[], my_gamma[];

int main() {
  L0 = SIZE;
  origin(-0.5*L0, -.5*L0);
  stokes = true;
  rho1 = RHO;
  rho2 = RHO;
  mu1 = MU;
  mu2 = MU;
  N = 1 << LEVEL;
  TOLERANCE = MY_TOLERANCE;
  periodic(right);
  run();
}

event init (i=0) {
  fraction(f, 1. - sq(x/(RADIUS)) - sq(y/(RADIUS)) - sq(z/RADIUS));
}

event def_gamma (i++) {
  foreach() {
    if (GAMMA) my_gamma[] = 1.;
    else my_gamma[] = 0.;
  }
}

#if AMR
  event adapt (i++) {
    foreach() {
      if (GAMMA) noise_gamma[] = noise();
      else noise_gamma[] = 0.;
    }
    boundary({noise_gamma});
    adapt_wavelet ({noise_gamma,u}, (double []){1.e-6,U_TOL, U_TOL}, maxlevel = LEVEL);
  }
#endif

event output (i = 5) {
  fprintf(stdout, "%.17g %.17g\n", normf(f).avg, normf(f).max);
}

event movie (i = 5) {
  view (fov = 10, phi=pi/10, theta=pi/2+pi/4+pi/10, tx = 0., ty = -0., tz=-8., bg = {0.3,0.4,0.6});
  draw_vof ("f", lw=2);
  squares("my_gamma",n = {1,0,0}, alpha = 1.e-12, linear=false);
  cells (n = {1,0,0}, alpha = 1.e-12);
  squares("my_gamma",n = {0,1,0}, alpha = 1.e-12, linear=false);
  cells (n = {0,1,0}, alpha = 1.e-12);
  save ("gamma.png");

  view (fov = 10, phi=pi/10, theta=pi/2+pi/4+pi/10, tx = 0., ty = -0., tz=-8., bg = {0.3,0.4,0.6});
  draw_vof ("f", lw=2);
  squares("extended_n.x",n = {1,0,0}, alpha = 1.e-12, linear=false, min=-1, max=1);
  cells (n = {1,0,0}, alpha = 1.e-12);
  squares("extended_n.x",n = {0,1,0}, alpha = 1.e-12, linear=false, min=-1, max=1);
  cells (n = {0,1,0}, alpha = 1.e-12);
  save ("nx.png");

  view (fov = 10, phi=pi/10, theta=pi/2+pi/4+pi/10, tx = 0., ty = -0., tz=-8., bg = {0.3,0.4,0.6});
  draw_vof ("f", lw=2);
  squares("extended_n.y",n = {1,0,0}, alpha = 1.e-12, linear=false, min=-1, max=1);
  cells (n = {1,0,0}, alpha = 1.e-12);
  squares("extended_n.y",n = {0,1,0}, alpha = 1.e-12, linear=false, min=-1, max=1);
  cells (n = {0,1,0}, alpha = 1.e-12);
  save ("ny.png");
}

event end (i = 5) {
  return 1.;
}
