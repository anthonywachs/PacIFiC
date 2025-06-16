#ifndef LEVEL
  #define LEVEL 10
#endif
#ifndef SIZE
  #define SIZE 10
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
#ifndef MU_RATIO
  #define MU_RATIO 5.
#endif
#ifndef U_MAX
  #define U_MAX 1.
#endif
#ifndef GRAD_P
  #define GRAD_P (-8*MU*U_MAX/(SIZE*SIZE))
#endif
#ifndef CA
  #define CA 0.2
#endif
#ifndef G_SHEAR
  #define G_SHEAR (MU*U_MAX/CA)
#endif
#ifndef ND_EB
  #define ND_EB 0.0005
#endif
#ifndef E_BEND
  #define E_BEND (ND_EB*RADIUS*RADIUS*G_SHEAR)
#endif
#ifndef REF_CURV
  #define REF_CURV 1
#endif
#ifndef DT_MAX
  #define DT_MAX (2.5e-4)
#endif
#ifndef T_END
  #define T_END 100.
#endif
#ifndef MY_TOLERANCE
  #define MY_TOLERANCE 1.e-10
#endif
#ifndef U_TOL
  #define U_TOL (U_MAX*1.e-2)
#endif
#ifndef OUTPUT_FREQ
  #define OUTPUT_FREQ (t += 0.1)
#endif
#ifndef AMR
  #define AMR 1
#endif
#ifndef STOKES_FLOW
  #define STOKES_FLOW true
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
scalar noise_gamma[], sdivu[], levelset[];

int main() {
  periodic(right);
  L0 = SIZE;
  origin(-.5*L0, -.5*L0);
  stokes = STOKES_FLOW;
  rho1 = RHO;
  rho2 = RHO;
  mu2 = MU;
  mu1 = MU*MU_RATIO;
  N = 1 << LEVEL;
  TOLERANCE = MY_TOLERANCE;
  DT = DT_MAX;
  run();
}

u.n[bottom] = dirichlet(0);
u.n[top] = dirichlet(0);
u.t[bottom] = dirichlet(0);
u.t[top] = dirichlet(0);

event init (t=0) {
  double a, c;
  c = 1.3858189;
  a = RADIUS/c;
  // We define below the local coordinates of the RBC and the parametric angle
  #define MY_X (y + .475*SIZE - RADIUS)
  #define MY_Y x
  #define COSPHI2 (sq(MY_X/(a*c)))
  #define Y_PREF (0.207 + 2.003*COSPHI2 - 1.123*sq(COSPHI2))
  fraction(f, sq(Y_PREF) - sq(Y_PREF*MY_X/(a*c)) - sq(2*MY_Y/(a*c)));
}

event acceleration (i++) {
  face vector ap = a;
  foreach_face(x) {
    ap.x[] -= GRAD_P;
  }
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

event logfile (OUTPUT_FREQ) {
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
    "%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
    t, normf(un).max, normf(un).avg, normf(u.x).max, normf(u.x).avg,
    normf(u.y).max, normf(u.y).avg,
    normf(J).max, J_avg,
    normf(T_s.x.x).max, T_s_avg_xx, normf(T_s.x.y).max,
    T_s_avg_xy, normf(T_s.y.y).max, T_s_avg_yy,
    normf(aen).max, aen_avg, abs_J_avg, abs_T_s_avg_xx, abs_T_s_avg_xy, abs_T_s_avg_yy, sdivu_avg, D, normf(m_bend.x.x).avg, normf(m_bend.x.x).max);
}

event movie (OUTPUT_FREQ) {
  // view (fov = 20, camera = "front", bg = {1,1,1});
  view (ty = (SIZE/2 - 2.2*RADIUS)/SIZE, fov = 1, camera = "front", bg = {1,1,1});
  cells ();
  squares ("u.x", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("ux_cells.mp4");

  // view (fov = 20, camera = "front", bg = {1,1,1});
  view (ty = (SIZE/2 - 2.2*RADIUS)/SIZE, fov = 1, camera = "front", bg = {1,1,1});
  // cells ();
  squares ("u.x", linear = true);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("ux.mp4");

  // view (fov = 20, camera = "front", bg = {1,1,1});
  view (ty = (SIZE/2 - 2.2*RADIUS)/SIZE, fov = 1, camera = "front", bg = {1,1,1});
  // cells ();
  squares ("J", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("J.mp4");

  // view (fov = 20, camera = "front", bg = {1,1,1});
  view (ty = (SIZE/2 - 2.2*RADIUS)/SIZE, fov = 1, camera = "front", bg = {1,1,1});
  // cells ();
  squares ("m_bend.x.x", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("mxx.mp4");

  view (ty = (SIZE/2 - 2.2*RADIUS)/SIZE, fov = 1, camera = "front", bg = {1,1,1});
  // view (fov = 20, camera = "front", bg = {1,1,1});
  // cells ();
  squares ("u.y", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("uy.mp4");

  // // view (fov = 20, camera = "front", bg = {1,1,1});
  // view (ty = (SIZE/2 - 2.2*RADIUS)/SIZE, fov = 1, camera = "front", bg = {1,1,1});
  // cells ();
  // squares ("u.x", linear = false);
  // draw_vof ("f", lw = 2);
  // box (notics = true);
  // save ("ux_cells.png");
  //
  // // view (fov = 20, camera = "front", bg = {1,1,1});
  // view (ty = (SIZE/2 - 2.2*RADIUS)/SIZE, fov = 1, camera = "front", bg = {1,1,1});
  // // cells ();
  // squares ("u.x", linear = true);
  // draw_vof ("f", lw = 2);
  // box (notics = true);
  // save ("ux.png");
  //
  // // view (fov = 20, camera = "front", bg = {1,1,1});
  // view (ty = (SIZE/2 - 2.2*RADIUS)/SIZE, fov = 1, camera = "front", bg = {1,1,1});
  // // cells ();
  // squares ("J", linear = false);
  // draw_vof ("f", lw = 2);
  // box (notics = true);
  // save ("J.png");
  //
  // // view (fov = 20, camera = "front", bg = {1,1,1});
  // view (ty = (SIZE/2 - 2.2*RADIUS)/SIZE, fov = 1, camera = "front", bg = {1,1,1});
  // // cells ();
  // squares ("m_bend.x.x", linear = false);
  // draw_vof ("f", lw = 2);
  // box (notics = true);
  // save ("mxx.png");
  //
  // // view (fov = 20, camera = "front", bg = {1,1,1});
  // view (ty = (SIZE/2 - 2.2*RADIUS)/SIZE, fov = 1, camera = "front", bg = {1,1,1});
  // cells ();
  // squares ("u.y", linear = false);
  // draw_vof ("f", lw = 2);
  // box (notics = true);
  // save ("uy.png");
}

event end (t = T_END) {
// event end (i = 5000) {
  return 1.;
}
