#define EPS 1.e-14
#define SIZE 4.
#define RADIUS (1. + sqrt(2)*1.e-3)
#define MU 1.
#define RHO 1.
#define LEVEL 9
#define MAX_ITER_EXTENSION 100

scalar qs[];
vector qv[];

#include "grid/multigrid.h"
#include "eulerian_caps/navier-stokes/my_centered.h"
#include "eulerian_caps/capsule.h"
#include "eulerian_caps/normal_extension.h"
#include "view.h"

FILE * fp = NULL;

int main() {
  L0 = SIZE;
  origin (-0.5*L0, -0.5*L0);

  stokes = true;
  TOLERANCE = 1.e-10;
  rho1 = RHO;
  rho2 = RHO;
  mu1 = MU;
  mu2 = MU;

  int level = LEVEL;
  N = 1 << level;

  run();
}

u.n[bottom] = dirichlet(0);
u.n[top] = dirichlet(0);
u.n[left] = dirichlet(0);
u.n[right] = dirichlet(0);

event init (i = 0) {
  fraction(f, 1. - sq(x/(RADIUS)) - sq(y/(RADIUS)));
}

event init_q (i = 1; i+=2) {
  foreach() {
    if (GAMMA) {
      double my_alpha;
      if (fabs(x)>SEPS) my_alpha = atan2(y,x);
      else my_alpha = y > 0 ? pi/2 : -pi/2;
      qs[] = cos(my_alpha)*sin(my_alpha);
      qv.x[] = sin(my_alpha);
      qv.y[] = cos(my_alpha);
      // if (!(IS_INTERFACE_CELL(point,f)) && (i%2==0)){
      if (!(IS_INTERFACE_CELL(point,f))){
        qs[] = -1.;
        foreach_dimension() qv.x[] = -1.;
      }
    }
    else {
      qs[] = -1.;
      foreach_dimension() qv.x[] = -1.;
    }
  }
}

event extension_event (i = 2; i+=2) {
  if (i%2==0) {
    foreach() {
      if (GAMMA) {
        if (i == 2) {
          double my_alpha;
          if (fabs(x)>SEPS) my_alpha = atan2(y,x);
          else my_alpha = y > 0 ? pi/2 : -pi/2;
          extended_n.x[] = cos(my_alpha);
          extended_n.y[] = sin(my_alpha);
        }
        // if (i == 4) {
        //   extended_n.x[] = - grad_wide_caps.x[]/ng_wide_caps[];
        //   extended_n.y[] = - grad_wide_caps.y[]/ng_wide_caps[];
        // }
      }
    }
    normal_scalar_extension((scalar *){qs, qv});
  }
}

scalar error_qs[];
vector error_qv[];
double errormax_qs = -HUGE;
double errormax_qv = -HUGE;
event output (i=2) {
  double errorl2_qs = 0;
  double errorl2_qv = 0;
  int nb_mb_cells = 0;
  foreach() {
    if (GAMMA) {
      nb_mb_cells++;
      double my_alpha;
      if (fabs(x)>SEPS) my_alpha = atan2(y,x);
      else my_alpha = y > 0 ? pi/2 : -pi/2;
      error_qs[] = sqrt(sq(sin(my_alpha)*cos(my_alpha) - qs[]));
      errorl2_qs += error_qs[];
      if (error_qs[] > errormax_qs) errormax_qs = error_qs[];
      error_qv.x[] = sqrt(sq(sin(my_alpha) - qv.x[]));
      error_qv.y[] = sqrt(sq(cos(my_alpha) - qv.y[]));
      errorl2_qv += error_qv.x[] + error_qv.y[];
      if (error_qs[] > errormax_qs) errormax_qs = error_qs[];
      if (max(error_qv.x[],error_qv.y[]) > errormax_qv) errormax_qv = max(error_qv.x[],error_qv.y[]);
    }
    else {
      error_qs[] = 0;
      foreach_dimension() error_qv.x[] = 0;
    }
  }
  if (nb_mb_cells != 0) {
    errorl2_qs /= nb_mb_cells;
    errorl2_qv /= 2*nb_mb_cells;
  }
  fprintf(stdout, "%g %g %g %g %g\n", t, errorl2_qs, errormax_qs, errorl2_qv, errormax_qv);
}

event picture (i++) {
  view (fov = 20, camera = "front", bg = {1,1,1});
  squares ("qs", linear = false, min=-1, max=1);
  draw_vof ("f", lw = 2);
  box (notics = true);
  if (i==1) save ("qs_init.png");
  if (i==2) {
    save ("qs_extended.png");

    view (fov = 20, camera = "front", bg = {1,1,1});
    squares ("error_qs", linear = false, min=0, max=errormax_qs );
    draw_vof ("f", lw = 2);
    box (notics = true);
    save ("error_qs.png");
  }

  view (fov = 20, camera = "front", bg = {1,1,1});
  squares ("qv.x", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  if (i==1) save ("qv.x_init.png");
  if (i==2) {
    save ("qv_x_extended.png");

    view (fov = 20, camera = "front", bg = {1,1,1});
    squares ("error_qv.x", linear = false, min=0, max=errormax_qv );
    draw_vof ("f", lw = 2);
    box (notics = true);
    save ("error_qv_x.png");
  }

  view (fov = 20, camera = "front", bg = {1,1,1});
  squares ("qv.y", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  if (i==1) save ("qv.y_init.png");
  if (i==2) {
    save ("qv_y_extended.png");

    view (fov = 20, camera = "front", bg = {1,1,1});
    squares ("error_qv.y", linear = false, min=0, max=errormax_qv );
    draw_vof ("f", lw = 2);
    box (notics = true);
    save ("error_qv_y.png");
  }
}

event end (i = 2) {
  return 1.;
}
