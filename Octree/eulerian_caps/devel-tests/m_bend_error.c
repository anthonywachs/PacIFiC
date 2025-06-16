// #define JACOBI 1
#define EPS 1.e-14
#define SIZE 4.
#define SHEAR_RATE 1.
#define RADIUS (1. + sqrt(2)*1.e-3)
#define MU 1.
#define CA 0.00125
// #define E_BEND (MU*SIZE*SHEAR_RATE/CA)
#define E_BEND 1.
#define RHO 1.
#define T_BC (1./SHEAR_RATE)
#define TEND (10./SHEAR_RATE)
#define AREA_DILATATION_MODULUS 1.
#define RESTRICT_DT .05
#define LEVEL 6
#define DT_MAX 1.e-7

#define NORMAL_EXTENSION 1
#define OUTPUT_FREQ 1
#define IMPOSE_CIRCULAR_NORMALS 0

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
  DT = DT_MAX;

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

tensor m_error[];
vector q_error[];
vector q_exact[];
typedef struct { pseudo_v x, y;} pseudo_t;
pseudo_t avg_error, max_error, avg_m;
pseudo_v q_avg, q_max, a_avg, a_max, q_err_avg, q_err_max;
event compute_error (i++) {
  foreach_dimension() {
    avg_m.x.x = 0;
    avg_m.x.y = 0;
    avg_error.x.x = 0;
    avg_error.x.y = 0;
    max_error.x.x = 0;
    max_error.x.y = 0;
    q_avg.x = 0;
    q_max.x = 0;
    q_err_avg.x = 0;
    q_err_max.x = 0;
    a_avg.x = 0;
    a_max.x = 0;
  }
  int nb_mb_cells = 0;
  foreach() {
    if (GAMMA) {
      nb_mb_cells++;
      double my_alpha;
      if (fabs(x)>SEPS) my_alpha = atan2(y,x);
      else my_alpha = y > 0 ? pi/2 : -pi/2;
      double r2 = sq(x) + sq(y);
      m_error.x.x[] = fabs(-(y*sin(my_alpha))*E_BEND/(r2) - m_bend.x.x[]);
      m_error.x.y[] = fabs(-(-x*sin(my_alpha))*E_BEND/(r2) - m_bend.x.y[]);
      m_error.y.x[] = fabs(-(-y*cos(my_alpha))*E_BEND/(r2) - m_bend.y.x[]);
      m_error.y.y[] = fabs(-(x*cos(my_alpha))*E_BEND/(r2) - m_bend.y.y[]);
      // q_exact.x[] = (2*E_BEND/r2)*(cos(my_alpha)*sq(sin(my_alpha))*cos(2*my_alpha));
      // q_exact.y[] = (2*E_BEND/r2)*(-sq(cos(my_alpha))*sin(my_alpha)*cos(2*my_alpha));
      q_exact.x[] = (2*x*y*E_BEND/sq(r2))*sin(my_alpha)*cos(2*my_alpha);
      q_exact.y[] = -(2*x*y*E_BEND/sq(r2))*cos(my_alpha)*cos(2*my_alpha);
      // q_exact.x[] = E_BEND/r2*(sin(my_alpha)*sin(2*my_alpha) + cos(my_alpha)*cos(2*my_alpha));
      // q_exact.y[] = E_BEND/r2*(cos(my_alpha)*sin(2*my_alpha) - sin(my_alpha)*cos(2*my_alpha));
      q_error.x[] = fabs(q_exact.x[] - q_bend.x[]);
      q_error.y[] = fabs(q_exact.y[] - q_bend.y[]);
      foreach_dimension() {
        avg_m.x.x += fabs(m_bend.x.x[]);
        avg_m.x.y += fabs(m_bend.x.y[]);
        avg_error.x.x += m_error.x.x[];
        avg_error.x.y += m_error.x.y[];
        q_avg.x += fabs(q_bend.x[]);
        q_err_avg.x += fabs(q_error.x[]);
        a_avg.x += fabs(centered_ae.x[]);
        if (m_error.x.x[]>max_error.x.x) max_error.x.x = m_error.x.x[];
        if (m_error.x.y[]>max_error.x.y) max_error.x.y = m_error.x.y[];
        if (q_bend.x[] > q_max.x) q_max.x = q_bend.x[];
        if (q_error.x[] > q_err_max.x) q_err_max.x = q_error.x[];
        if (centered_ae.x[] > a_max.x) a_max.x = centered_ae.x[];
      }
    }
    else {
      foreach_dimension() {
        m_error.x.x[] = 0;
        m_error.x.y[] = 0;
        q_error.x[] = 0;
      }
    }
  }
  if (nb_mb_cells != 0) {
    foreach_dimension() {
      avg_error.x.x /= nb_mb_cells;
      avg_error.x.y /= nb_mb_cells;
      avg_m.x.x /= nb_mb_cells;
      avg_m.x.y /= nb_mb_cells;
      q_avg.x /= nb_mb_cells;
      q_err_avg.x /= nb_mb_cells;
      a_avg.x /= nb_mb_cells;
    }
  }
  fprintf(stdout, "%g %g %g %g %g %g %g %g %g %g %g %g %g ", t, avg_m.x.x, avg_m.x.y, avg_m.y.x, avg_m.y.y, avg_error.x.x, avg_error.x.y, avg_error.y.x, avg_error.y.y, max_error.x.x, max_error.x.y,
  max_error.y.x, max_error.y.y);
  fprintf(stdout, "%g %g %g %g ", q_avg.x, q_avg.y, q_max.x, q_max.y);
  fprintf(stdout, "%g %g %g %g ", a_avg.x, a_avg.y, a_max.x, a_max.y);
  fprintf(stdout, "%g %g %g %g ", q_err_avg.x, q_err_avg.y, q_err_max.x, q_err_max.y);
  fprintf(stdout, "\n");
}

event pictures (i = 2) {
  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("u.x", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("ux.png");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("u.y", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("uy.png");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("aen", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("ae.png");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("T_s.x.x", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("Tsxx.png");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("T_s.x.y", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("Tsxy.png");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("T_s.y.y", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("Tsyy.png");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("m_bend.x.x", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("mxx.png");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("m_bend.x.y", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("mxy.png");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("m_bend.y.y", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("myy.png");

  view (fov = 20, camera = "front", bg = {1,1,1});
  // cells ();
  // squares ("q_bend.x", linear = false, min = -normf(q_exact.x).max, max = normf(q_exact.x).max);
  squares ("q_bend.x", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("qx.png");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("q_bend.y", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("qy.png");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("extended_n.x", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("nx.png");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("extended_n.y", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("ny.png");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("centered_ae.x", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("aex.png");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("centered_ae.y", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("aey.png");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("m_error.x.x", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("m_err_xx.png");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("m_error.x.y", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("m_err_xy.png");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("m_error.y.x", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("m_err_yx.png");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("m_error.y.y", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("m_err_yy.png");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("q_error.x", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("q_err_x.png");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("q_error.y", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("q_err_y.png");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("q_exact.x", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("q_exact_x.png");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("q_exact.y", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("q_exact_y.png");
}

event end (i = 2) {
  return 1.;
}
