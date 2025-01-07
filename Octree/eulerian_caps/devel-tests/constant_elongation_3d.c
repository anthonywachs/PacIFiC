#define U_RIGHT 1.
#define L0 1.
#define RHO 1.
#define MU 1.
#define OUTPUT_FREQ i++
#define TEND 1.
#define LEVEL 6
#define U_TOL 0.01

#ifndef THETA
  #define THETA (0.)
#endif
#ifndef PHI
  #define PHI (pi/3.)
#endif
#ifndef PSI
  #define PSI (0.)
#endif

// #define U_X(X,Y) ((sq(cos(THETA)) - sq(sin(THETA)))*(X) + 2*sin(THETA)*cos(THETA)*(Y))
// #define U_Y(X,Y) (-(sq(cos(THETA)) - sq(sin(THETA)))*(Y) + 2*sin(THETA)*cos(THETA)*(X))

#define U1_X(X,Y,Z) ((sq(cos(THETA)) - sq(sin(THETA)))*(X) + 2*sin(THETA)*cos(THETA)*(Y))
#define U1_Y(X,Y,Z) (-(sq(cos(THETA)) - sq(sin(THETA)))*(Y) + 2*sin(THETA)*cos(THETA)*(X))
#define U1_Z(X,Y,Z) (0.)

#define U2_X(X,Y,Z) ((sq(cos(PHI)) - sq(sin(PHI)))*(U1_X(X,Y,Z)) - 2*sin(PHI)*cos(PHI)*(U1_Z(X,Y,Z)))
#define U2_Y(X,Y,Z) (U1_Y(X,Y,Z))
#define U2_Z(X,Y,Z) ((sq(cos(PHI)) - sq(sin(PHI)))*(U1_Z(X,Y,Z)) + 2*sin(PHI)*cos(PHI)*(U1_X(X,Y,Z)))

// #define U3_X(X,Y,Z) ((sq(cos(PSI)) - sq(sin(PSI)))*(U2_X(X,Y,Z)) - 2*sin(PSI)*cos(PSI)*(U2_Z(X,Y,Z)))
// #define U3_Y(X,Y,Z) (U2_Y(X,Y,Z))
// #define U3_Z(X,Y,Z) ((sq(cos(PSI)) - sq(sin(PSI)))*(U2_Z(X,Y,Z)) + 2*sin(PSI)*cos(PSI)*(U2_X(X,Y,Z)))

// #define U_X(X,Y,Z) U2_X(X,Y,Z)
// #define U_Y(X,Y,Z) U2_Y(X,Y,Z)
// #define U_Z(X,Y,Z) U2_Z(X,Y,Z)

#define U_X(X,Y,Z) 0.
#define U_Y(X,Y,Z) 0.
#define U_Z(X,Y,Z) 0.

#include "grid/octree.h"
#include "eulerian_caps/navier-stokes/my_centered.h"
#include "eulerian_caps/capsule.h"
#include "eulerian_caps/elasticity.h"
#include "view.h"

scalar my_gamma[];

void impose_u() {
  foreach() {
    u.x[] = U_X(x,y,z);
    u.y[] = U_Y(x,y,z);
    u.z[] = U_Z(x,y,z);
  }
}

int main() {
  size(L0);
  origin(-.5*L0, -.5*L0, -.5*L0);
  rho1 = RHO;
  rho2 = RHO;
  mu1 = MU;
  mu2 = MU;
  DT = 1.e-2;
  stokes = true;
  TOLERANCE = 1.e-10;
  run();
}

/* Boundary conditions:*/
u.n[left] = dirichlet(U_X(x,y,z));
u.n[right] = dirichlet(U_X(x,y,z));
u.n[top] = dirichlet(U_Y(x,y,z));
u.n[bottom] = dirichlet(U_Y(x,y,z));
u.n[front] = dirichlet(U_Z(x,y,z));
u.n[back] = dirichlet(U_Z(x,y,z));
u.t[left] = dirichlet(U_Y(x,y,z));
u.t[right] = dirichlet(U_Y(x,y,z));
u.t[top] = dirichlet(U_Z(x,y,z));
u.t[bottom] = dirichlet(U_Z(x,y,z));
u.t[front] = dirichlet(U_Y(x,y,z));
u.t[back] = dirichlet(U_Y(x,y,z));
u.r[left] = dirichlet(U_Z(x,y,z));
u.r[right] = dirichlet(U_Z(x,y,z));
u.r[top] = dirichlet(U_X(x,y,z));
u.r[bottom] = dirichlet(U_X(x,y,z));
u.r[front] = dirichlet(U_Z(x,y,z));
u.r[back] = dirichlet(U_Z(x,y,z));
p[top] = neumann(0.);
p[bottom] = neumann(0.);
p[right] = neumann(0.);
p[left] = neumann(0.);
p[front] = neumann(0.);
p[back] = neumann(0.);

event init (i = 0) {
  fraction(f, cos(PSI)*(-x*sin(THETA)+y*cos(THETA)) + sin(PSI)*(sin(PHI)*(x*cos(THETA) + y*sin(THETA)) - z*cos(PHI)));
  // fraction(f, sin(THETA)*(x*cos(PHI)+sin(PHI)*(y*cos(PSI)-z*sin(PSI))) + cos(THETA)*(y*cos(PSI)-z*sin(PSI)));
}

event impose_fraction (i++) {
  fraction(f, cos(PSI)*(-x*sin(THETA)+y*cos(THETA)) + sin(PSI)*(sin(PHI)*(x*cos(THETA) + y*sin(THETA)) - z*cos(PHI)));
  // fraction(f, sin(THETA)*(x*cos(PHI)+sin(PHI)*(y*cos(PSI)-z*sin(PSI))) + cos(THETA)*(y*cos(PSI)-z*sin(PSI)));
}

event properties (i++) {
  fraction(f, cos(PSI)*(-x*sin(THETA)+y*cos(THETA)) + sin(PSI)*(sin(PHI)*(x*cos(THETA) + y*sin(THETA)) - z*cos(PHI)));
  // fraction(f, sin(THETA)*(x*cos(PHI)+sin(PHI)*(y*cos(PSI)-z*sin(PSI))) + cos(THETA)*(y*cos(PSI)-z*sin(PSI)));
}

event viscous_term (i++) {
  impose_u();
}

event acceleration (i++) {
  impose_u();
}

scalar noise_gamma[];
event adapt (i++) {
  foreach() {
    if (GAMMA) noise_gamma[] = noise();
    else noise_gamma[] = 0.;
  }
  boundary({noise_gamma});
  adapt_wavelet ({noise_gamma,u}, (double []){1.e-6,U_TOL, U_TOL}, maxlevel = LEVEL);
}

event output (OUTPUT_FREQ) {
  double J_avg = 0., J_int=0.;
  pseudo_t G_avg, G_int;
  int nb_mb_cells = 0., nb_interface_cells = 0.;
  foreach_dimension() {
    G_avg.x.x = 0.;
    G_avg.x.y = 0.;
    G_avg.x.z = 0.;
    G_int.x.x = 0.;
    G_int.x.y = 0.;
    G_int.x.z = 0.;
  }
  foreach() {
    my_gamma[] = GAMMA;
    if ((GAMMA) && fabs(x)<.25) {
      nb_mb_cells += 1;
      J_avg += J[];
      foreach_dimension() {
        G_avg.x.x += G.x.x[];
        G_avg.x.y += G.x.y[];
        G_avg.x.z += G.x.z[];
      }
      if (IS_INTERFACE_CELL(point,f)) {
        nb_interface_cells += 1;
        J_int += J[];
        foreach_dimension() {
          G_int.x.x += G.x.x[];
          G_int.x.y += G.x.y[];
          G_int.x.z += G.x.z[];
        }
      }
    }
  }
  if (nb_mb_cells) {
    J_avg /= nb_mb_cells;
    foreach_dimension() {
      G_avg.x.x /= nb_mb_cells;
      G_avg.x.y /= nb_mb_cells;
      G_avg.x.z /= nb_mb_cells;
    }
  }
  if (nb_interface_cells) {
    J_int /= nb_interface_cells;
    foreach_dimension() {
      G_int.x.x /= nb_interface_cells;
      G_int.x.y /= nb_interface_cells;
      G_int.x.z /= nb_interface_cells;
    }
  }
  fprintf(stdout, "%.17g %.17g %.17g %.17g ", t, (THETA), J_int, J_avg);
  foreach_dimension() {
    fprintf(stdout, "%.17g %.17g %.17g %.17g %.17g %.17g ", G_avg.x.x, G_int.x.x, G_avg.x.y, G_int.x.y, G_avg.x.z, G_int.x.z);
  }
  fprintf(stdout,"\n");
  fflush(fout);
}

event movie (i = 2) {
  // view (fov = 20, theta=pi/10, phi=pi/5, psi=0, tx = 0., ty = -0., tz=-8., bg = {0.3,0.4,0.6});
  view (fov = 20, theta=pi/10, phi=pi/10, psi=0, bg = {0.3,0.4,0.6});
  draw_vof ("f", lw=2);
  squares("my_gamma",n = {cos(THETA)*sin(PHI)*sin(-PSI)+sin(-THETA)*cos(PSI), sin(-THETA)*sin(PHI)*sin(-PSI)+cos(THETA)*cos(PSI), cos(PHI)*sin(-PSI)}, alpha = 1.e-12, linear=false);
  // squares("my_gamma",n = {0,0,1}, alpha = 1.e-12, linear=false);
  cells (n = {1,0,0}, alpha = 1.e-12);
  cells (n = {0,1,0}, alpha = 1.e-12);
  cells (n = {0,0,1}, alpha = 1.e-12);
  save ("ux.png");
}

// event end (t = TEND) {
event end (i = 2) {
  return (1.);
}
