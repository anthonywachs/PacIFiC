#define EPS 1.e-14
#define SIZE 4.
#define RADIUS (1. + sqrt(2)*1.e-3)
#define MU 1.
#define RHO 1.
#define LEVEL level

scalar q[];
double error_min = 0., error_max = 0.;
int level;

#include "grid/multigrid.h"
#include "eulerian_caps/navier-stokes/my_centered.h"
#include "eulerian_caps/capsule.h"
#include "eulerian_caps/normal_extension.h"
#include "view.h"

FILE * fp = NULL;

int main() {
  L0 = SIZE;
  origin (-0.5*L0, -0.5*L0);

  DT = 1.;

  stokes = true;
  TOLERANCE = 1.e-10;
  rho1 = RHO;
  rho2 = RHO;
  mu1 = MU;
  mu2 = MU;

  for (LEVEL = 4; LEVEL <= 9; LEVEL++) {
    N = 1 << LEVEL;
    run();
  }
}

u.n[bottom] = dirichlet(0);
u.n[top] = dirichlet(0);
u.n[left] = dirichlet(0);
u.n[right] = dirichlet(0);

event init (i = 0) {
  fraction(f, 1. - sq(x/(RADIUS)) - sq(y/(RADIUS)));
}

scalar interface_distance[];
event extension_event (i++) {
  foreach() {
    if (GAMMA) {
      extended_n.x[] = - grad_wide_caps.x[]/ng_wide_caps[];
      extended_n.y[] = - grad_wide_caps.y[]/ng_wide_caps[];

      if (i == 3) {
        double my_alpha;
        if (fabs(x)>SEPS) my_alpha = atan2(y,x);
        else my_alpha = y > 0 ? pi/2 : -pi/2;
        extended_n.x[] = cos(my_alpha);
        extended_n.y[] = sin(my_alpha);
      }
    }
  }
  if (i == 2) {
    initialize_normals_for_extension(extended_n);
    normal_vector_extension(extended_n);
  }
  if (i == 3) {
    foreach() {
      double r = sqrt(sq(x) + sq(y));
      interface_distance[] = r - RADIUS;
    }
    boundary({interface_distance});
    foreach() {
      if (GAMMA) {
        foreach_dimension() extended_n.x[] = (interface_distance[1,0] - interface_distance[-1,0])/(2*Delta);
        foreach_dimension() extended_n.x[] /= norm(extended_n);
      }
    }
  }
}

scalar error_normals[];
scalar error_normals_interface[];
double my_delta = 0.;
event output (i++) {
  foreach() {
    my_delta = Delta;
    if (GAMMA) {
      double my_alpha;
      if (fabs(x)>SEPS) my_alpha = atan2(y,x);
      else my_alpha = y > 0 ? pi/2 : -pi/2;
      error_normals[] = sq(extended_n.x[] - cos(my_alpha));
      error_normals[] += sq(extended_n.y[] - sin(my_alpha));
      error_normals[] = sqrt(error_normals[]);
      if (interfacial(point, f)) {
        error_normals_interface[] = sq(extended_n.x[] - cos(my_alpha));
        error_normals_interface[] += sq(extended_n.y[] - sin(my_alpha));
        error_normals_interface[] = sqrt(error_normals_interface[]);
      }
    }
    else {
      foreach_dimension() error_normals[] = 0.;
    }
  }
  if (i == 1) {
    error_min = 0.;
    error_max = normf(error_normals).max;
  }
  fprintf(stdout, "%g %g %g %g %g %g\n", t, (2*RADIUS)/my_delta, normf(error_normals).avg, normf(error_normals).max, normf(error_normals_interface).avg, normf(error_normals_interface).max);
  fflush(stdout);
}

event picture (i++) {
  if (LEVEL==6) {
    view (fov = 20, camera = "front", bg = {1,1,1});
    cells ();
    squares ("extended_n.x", linear = false);
    draw_vof ("f", lw = 2);
    box (notics = true);
    if (i == 1) {
      save ("nx_init.png");
    }
    else if (i == 2) {
      save ("nx_extended.png");
    }
    else if (i == 3) {
      save ("nx_exact.png");
    }

    view (fov = 20, camera = "front", bg = {1,1,1});
    cells ();
    squares ("extended_n.y", linear = false);
    draw_vof ("f", lw = 2);
    box (notics = true);
    if (i == 1) {
      save ("ny_init.png");
    }
    else if (i == 2) {
      save ("ny_extended.png");
    }
    else if (i == 3) {
      save ("ny_exact.png");
    }
  }
}

event end (i = 3) {
  return 1.;
}
