#define EPS 1.e-14
#define RHO 1.
#define DIAMETER 2.
#define SIZE 4.
#define MU 0.001
#define G_SHEAR 1.

#define PROJECT_G 1
#define LOCAL_BCG 0
#define NB_MB_EXTRA_LAYERS 2

vector error_acceleration[];

#include "grid/multigrid.h"
#include "eulerian_caps/navier-stokes/my_centered.h"
#include "eulerian_caps/capsule.h"
#include "eulerian_caps/elasticity.h"
#include "view.h"


FILE * acc_file[2];
int LEVEL;
scalar layers[];

int main() {
  L0 = SIZE;
  origin (-0.5*L0, -0.5*L0);

  stokes = true;
  TOLERANCE = 1.e-10;
  rho1 = RHO;
  rho2 = RHO;
  mu1 = MU;
  mu2 = MU;

  LEVEL = 6;
  N = 1 << LEVEL;

  run();
}

u.n[bottom] = dirichlet(0);
u.n[top] = dirichlet(0);
u.n[left] = dirichlet(0);
u.n[right] = dirichlet(0);
u.t[bottom] = dirichlet(0);
u.t[top] = dirichlet(0);
u.t[left] = dirichlet(0);
u.t[right] = dirichlet(0);

event init (t=0) {
  fraction(f, 1. - sq(x/(DIAMETER/2)) - sq(y/(DIAMETER/2)));
  char buf[0x100];
  for (int k=0; k<2; k++) {
    snprintf(buf, sizeof(buf), "acceleration_%d.csv", k+1);
    acc_file[k] = fopen(buf,"w");
  }
}

/** Below we impose the membrane stress $$ \bm{T_s} = \begin{bmatrix}
sin(\theta) & cos(\theta) \\
sin(\theta) & cos(\theta)
\end{bmatrix} $$*/
event acceleration (i++)
{
  foreach()
  {
    double my_alpha;
    if (fabs(x)>SEPS) my_alpha = atan2(y,x);
    else my_alpha = y > 0 ? pi/2 : -pi/2;
    T_s.x.x[] = sin(my_alpha);
    T_s.x.y[] = cos(my_alpha);
    T_s.y.x[] = sin(my_alpha);
    T_s.y.y[] = cos(my_alpha);

    if (GAMMA){
      if (i<2) {
        extended_n.x[] = cos(my_alpha);
        extended_n.y[] = sin(my_alpha);
      }
    }
  }
  boundary((scalar *){T_s});
}

event impose_u (i++) {
  fraction(f, 1. - sq(x/(DIAMETER/2)) - sq(y/(DIAMETER/2)));
  foreach() {
    foreach_dimension()
    u.x[] = 0.;
  }
}

event color_layers (i++) {
  foreach() {
    if (GAMMA) {
      double R = sqrt(sq(x) + sq(y));
      double DR = fabs(DIAMETER/2 - R);
      if (DR < Delta*sqrt(2)/2.) {
        layers[] = 4.;
      }
      else if (DR < 3*Delta*sqrt(2)/2.) {
          layers[] = -4.;
      }
      else if (DR < 5*Delta*sqrt(2)/2.) {
          layers[] = 4.;
      }
      else if (DR > 5*Delta*sqrt(2)/2.) {
          layers[] = -4.;
      }
    }
    else layers[] = 0;
  }
}

event logfile (i++) {
  foreach() {
    double my_alpha;
    if (fabs(x)>SEPS) my_alpha = atan2(y,x);
    else my_alpha = y > 0 ? pi/2 : -pi/2;
    double R = sqrt(sq(x) + sq(y));
    double ax_exact = 0., ay_exact = 0.;
    if (BIGAMMA(0,0) && BIGAMMA(1,0) && BIGAMMA(-1,0)) {
    // if (BIGAMMA(0,0) && BIGAMMA(0,-1) && BIGAMMA(0,-1)) {
    // if (BIGAMMA(0,0)) {
      ax_exact = (1 + sqrt(2)*cos(2*my_alpha + pi/4))/(2*R);
      ay_exact = (1 - sqrt(2)*sin(2*my_alpha + pi/4))/(2*R);
      error_acceleration.x[] = centered_ae.x[] - ax_exact;
      error_acceleration.y[] = centered_ae.y[] - ay_exact;
      if (i>1) fprintf(acc_file[i-2], "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
        R, my_alpha, centered_ae.x[], ax_exact, centered_ae.y[], ay_exact,
        Delta, DIAMETER);
    }
  }
  fprintf(stdout, "%g\t%g\n", normf(error_acceleration.x).max,
    normf(error_acceleration.y).max);
}

event pictures (i++) {
  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("error_acceleration.x", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("error_acceleration_x.png");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("error_acceleration.y", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("error_acceleration_y.png");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells ();
  squares ("layers", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  save ("mb_layers.png");
}

event end (i = 3) {
  return 1.;
}
