// #define JACOBI 1
#define EPS 1.e-14
#define SIZE 4.
#define SHEAR_RATE 1.
#define RADIUS (1. + sqrt(2)*1.e-3)
#define MU 1.
#define CA 0.025
#define G_SHEAR (MU*SIZE*SHEAR_RATE/CA)
#define RHO 1.
#define T_BC (1./SHEAR_RATE)
#define TEND (10./SHEAR_RATE)
#define AREA_DILATATION_MODULUS 1.
#define RESTRICT_DT .05
#define LEVEL 6

#define LOCAL_BCG 0
#define NB_MB_EXTRA_LAYERS 2
#define PROJECT_G 1
#define NORMAL_EXTENSION 0
#define OUTPUT_FREQ 100

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
u.t[bottom] = dirichlet((t < T_BC ? SHEAR_RATE*y : 0.));
u.t[top] = dirichlet((t < T_BC ? SHEAR_RATE*y : 0.));

event init (t=0) {
  fraction(f, 1. - sq(x/(RADIUS)) - sq(y/(RADIUS)));
}

event logfile (i+=OUTPUT_FREQ) {
  scalar un[];
  int nb_mb_cells = 0;
  double sdivu_avg = 0., J_avg = 0., T_s_avg_xx = 0., T_s_avg_xy = 0.,
    T_s_avg_yy = 0., aen_avg = 0., abs_J_avg = 0., abs_T_s_avg_xx = 0.,
    abs_T_s_avg_xy = 0., abs_T_s_avg_yy = 0.;
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
  }
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
    "%.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
    t, normf(un).max, normf(un).avg, normf(u.x).max, normf(u.x).avg,
    normf(u.y).max, normf(u.y).avg,
    normf(J).max, J_avg,
    normf(T_s.x.x).max, T_s_avg_xx, normf(T_s.x.y).max,
    T_s_avg_xy, normf(T_s.y.y).max, T_s_avg_yy,
    normf(aen).max, aen_avg, abs_J_avg, abs_T_s_avg_xx, abs_T_s_avg_xy, abs_T_s_avg_yy, sdivu_avg);
}

event movie (i += OUTPUT_FREQ) {
  // view (fov = 20, camera = "front", bg = {1,1,1});
  // cells ();
  // squares ("sdivu", linear = false);
  // draw_vof ("f", lw = 2);
  // box (notics = true);
  // save ("sdivu.mp4");

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

  // view (fov = 20, camera = "front", bg = {1,1,1});
  // cells ();
  // squares ("my_grad_u.x.x", linear = false);
  // draw_vof ("f", lw = 2);
  // box (notics = true);
  // save ("grad_xx.mp4");
  //
  // view (fov = 20, camera = "front", bg = {1,1,1});
  // cells ();
  // squares ("my_grad_u.y.y", linear = false);
  // draw_vof ("f", lw = 2);
  // box (notics = true);
  // save ("grad_yy.mp4");
  //
  // view (fov = 20, camera = "front", bg = {1,1,1});
  // cells ();
  // squares ("sgrad_u.x.x", linear = false);
  // draw_vof ("f", lw = 2);
  // box (notics = true);
  // save ("sgrad_xx.mp4");
  //
  // view (fov = 20, camera = "front", bg = {1,1,1});
  // cells ();
  // squares ("sgrad_u.y.y", linear = false);
  // draw_vof ("f", lw = 2);
  // box (notics = true);
  // save ("sgrad_yy.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("aen", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("ae.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("J", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("J.mp4");

  // view (fov = 20, camera = "front", bg = {1,1,1});
  // cells ();
  // squares ("G.x.x", linear = false);
  // draw_vof ("f", lw = 2);
  // box (notics = true);
  // save ("Gxx.mp4");
  //
  // view (fov = 20, camera = "front", bg = {1,1,1});
  // cells ();
  // squares ("G.x.y", linear = false);
  // draw_vof ("f", lw = 2);
  // box (notics = true);
  // save ("Gxy.mp4");
  //
  // view (fov = 20, camera = "front", bg = {1,1,1});
  // cells ();
  // squares ("G.y.y", linear = false);
  // draw_vof ("f", lw = 2);
  // box (notics = true);
  // save ("Gyy.mp4");
  //
  // view (fov = 20, camera = "front", bg = {1,1,1});
  // cells ();
  // squares ("T_s.x.x", linear = false);
  // draw_vof ("f", lw = 2);
  // box (notics = true);
  // save ("Tsxx.mp4");
  //
  // view (fov = 20, camera = "front", bg = {1,1,1});
  // cells ();
  // squares ("T_s.x.y", linear = false);
  // draw_vof ("f", lw = 2);
  // box (notics = true);
  // save ("Tsxy.mp4");
  //
  // view (fov = 20, camera = "front", bg = {1,1,1});
  // cells ();
  // squares ("T_s.y.y", linear = false);
  // draw_vof ("f", lw = 2);
  // box (notics = true);
  // save ("Tsyy.mp4");

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

  // view (fov = 20, camera = "front", bg = {1,1,1});
  // cells ();
  // squares ("extended_n.x", linear = false);
  // draw_vof ("f", lw = 2);
  // box (notics = true);
  // save ("nx.mp4");
  //
  // view (fov = 20, camera = "front", bg = {1,1,1});
  // cells ();
  // squares ("extended_n.y", linear = false);
  // draw_vof ("f", lw = 2);
  // box (notics = true);
  // save ("ny.mp4");
}

event end (t = TEND) {
  return 1.;
}
