#define U_RIGHT 1.
#define LX 1.
#define RHO 1.
#define MU 1.
#define OUTPUT_FREQ 0.01
#define TEND 2.
#define LEVEL 6
#define AREA_DILATATION_MODULUS 1.

#define PROJECT_G 0
#define LOCAL_BCG 0
#define NB_MB_EXTRA_LAYERS 2
#define BGHOSTS 2

#ifndef THETA
  #define THETA (pi/6.)
#endif

#define U_X(X,Y) ((sq(cos(THETA)) - sq(sin(THETA)))*(X) + 2*sin(THETA)*cos(THETA)*(Y))
#define U_Y(X,Y) (-(sq(cos(THETA)) - sq(sin(THETA)))*(Y) + 2*sin(THETA)*cos(THETA)*(X))

#include "grid/multigrid.h"
#include "eulerian_caps/navier-stokes/my_centered.h"
#include "eulerian_caps/capsule.h"
#include "eulerian_caps/elasticity.h"
#include "eulerian_caps/skalak.h"
#include "view.h"

scalar my_gamma[];

int main() {
  size(LX);
  origin(-.5,-.5);
  rho1 = RHO;
  rho2 = RHO;
  mu1 = MU;
  mu2 = MU;

  DT = 1.e-3;

  stokes = true;
  TOLERANCE = 1.e-10;

  run();
}

/* Boundary conditions:*/
u.n[left] = dirichlet(U_X(x,y));
u.n[right] = dirichlet(U_X(x,y));
u.n[top] = dirichlet(U_Y(x,y));
u.n[bottom] = dirichlet(U_Y(x,y));
u.t[left] = dirichlet(U_Y(x,y));
u.t[right] = dirichlet(U_Y(x,y));
u.t[top] = dirichlet(U_X(x,y));
u.t[bottom] = dirichlet(U_X(x,y));
p[top] = neumann(0.);
p[bottom] = neumann(0.);
p[right] = neumann(0.);
p[left] = neumann(0.);

event init(i = 0) {
  fraction(f, y - tan(THETA)*x);
  foreach() {
    u.x[] = U_X(x,y);
    u.y[] = U_Y(x,y);
  }
}

event impose_ux (i++) {
  fraction(f, y - tan(THETA)*x);
  foreach() {
    u.x[] = U_X(x,y);
    u.y[] = U_Y(x,y);
  }
}

event end (t = TEND) {
  return (1.);
}

event output (t += OUTPUT_FREQ) {
  double J_avg = 0., Gxx_avg = 0., Gyy_avg = 0., Gxy_avg = 0., Ts_xx_avg = 0.,
        Ts_yy_avg = 0., Ts_xy_avg = 0.;
  int nb_mb_cells = 0., nb_interface_cells = 0.;
  foreach() {
    my_gamma[] = GAMMA;
    if ((GAMMA) && fabs(x)<.25) {
      nb_mb_cells += 1;
      if ((fabs(f[0,1] - f[]) > .5) || (fabs(f[0,-1] - f[]) > .5)) {
        J_avg += J[];
        Gxx_avg += G.x.x[];
        Gxy_avg += G.x.y[];
        Gyy_avg += G.y.y[];
        Ts_xx_avg += T_s.x.x[]/ngcaps[];
        Ts_xy_avg += T_s.x.y[]/ngcaps[];
        Ts_yy_avg += T_s.y.y[]/ngcaps[];
        nb_interface_cells += 1;
      }
    }
  }
  if (nb_interface_cells) {
      J_avg /= nb_interface_cells;
      Gxx_avg /= nb_interface_cells;
      Gxy_avg /= nb_interface_cells;
      Gyy_avg /= nb_interface_cells;
      Ts_xx_avg /= nb_interface_cells;
      Ts_xy_avg /= nb_interface_cells;
      Ts_yy_avg /= nb_interface_cells;
  }
  printf("%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n", t, THETA,
    J_avg, Gxx_avg, Gxy_avg, Gyy_avg, Ts_xx_avg, Ts_xy_avg, Ts_yy_avg);
  fflush(fout);
}
