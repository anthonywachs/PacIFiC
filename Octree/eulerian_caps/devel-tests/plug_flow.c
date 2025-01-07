// #define JACOBI 1
#define EPS 1.e-14
#define SIZE 4.
#define UX 1.
#define RADIUS (1. + sqrt(2)*1.e-3)
#define MU 1.
#define CA 0.025
#define G_SHEAR (MU*UX/CA)
#define RHO 1.
#define TEND 2.
#define T_PLUG .5
#define AREA_DILATATION_MODULUS 1.
#define RESTRICT_DT .05

#define LOCAL_BCG 0
#define NB_MB_EXTRA_LAYERS 0
#define PROJECT_G 1
#define OUTPUT_FREQ 100
#define LEVEL 6

#include "grid/multigrid.h"
#include "eulerian_caps/navier-stokes/my_centered.h"
#include "eulerian_caps/capsule.h"
#include "eulerian_caps/elasticity_2D.h"
#include "eulerian_caps/neo_hookean_2D.h"
#include "view.h"

FILE * fp = NULL;
scalar sdivu[];

void impose_velocity() {
  foreach() {
    u.x[] = t < T_PLUG ? y : UX;
  }
  boundary((scalar *){u});
}

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
u.t[bottom] = dirichlet((t < T_PLUG ? y : UX));
u.t[top] = dirichlet((t < T_PLUG ? y : UX));

event init (t=0) {
  fraction(f, 1. - sq(x/(RADIUS)) - sq(y/(RADIUS)));
}

event acceleration (i++) {
  impose_velocity();
}

event viscous_term (i++) {
  impose_velocity();
}

event advection_term (i++) {
  impose_velocity();
}

event end_timestep (i++) {
  impose_velocity();
}

/*
We substract the acceleration terms from uf since this test assumes 0 acceleration.
**/
event projection (i++) {
  face vector af = a;
  foreach_face() {
    uf.x[] -= fm.x[]*dt*af.x[];
    af.x[] = 0.;
  }
  boundary((scalar *){af, uf});
}

scalar my_gamma[];
event output_gamma (i++) {
  foreach() {
    if (GAMMA) my_gamma[] = 1.;
    else my_gamma[] = 0.;
  }
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
  squares ("my_gamma", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("gamma.mp4");
}

event end (t = TEND) {
  return 1.;
}
