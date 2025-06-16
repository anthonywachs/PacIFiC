// #define JACOBI 1
#define EPS 1.e-14
#define SIZE 4.
#define SHEAR_RATE 1.
#define RADIUS (1. + sqrt(2)*1.e-3)
#define MU 1.
#define CA 0.00125
// #define E_BEND (MU*SIZE*SHEAR_RATE/CA)
#define E_BEND 1
#define RHO 1.
#define T_BC (1./SHEAR_RATE)
#define TEND (10./SHEAR_RATE)
#define AREA_DILATATION_MODULUS 1.
#define RESTRICT_DT .0005
#define LEVEL 6

#define NORMAL_EXTENSION 1
#define OUTPUT_FREQ 1

#include "grid/multigrid.h"
#include "eulerian_caps/navier-stokes/my_centered.h"
#include "eulerian_caps/capsule.h"
#include "eulerian_caps/bending_2D.h"
#include "view.h"

FILE * fp = NULL;
scalar sdivu[];

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
u.t[bottom] = dirichlet((t < T_BC ? SHEAR_RATE*y : 0.));
u.t[top] = dirichlet((t < T_BC ? SHEAR_RATE*y : 0.));

event init (t=0) {
  fraction(f, 1. - sq(x/(RADIUS)) - sq(y/(RADIUS)));
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
  squares ("u.y", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("uy.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("aen", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("ae.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("T_s.x.x", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("Tsxx.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("T_s.x.y", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("Tsxy.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("T_s.y.y", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("Tsyy.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("m_bend.x.x", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("mxx.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("m_bend.x.y", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("mxy.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("m_bend.y.y", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("myy.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("q_bend.x", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("qx.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("q_bend.y", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("qy.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("extended_n.x", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("nx.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("extended_n.y", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("ny.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("centered_ae.x", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("aex.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("centered_ae.y", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("aey.mp4");
}

event end (t = TEND) {
  return 1.;
}
