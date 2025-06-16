#define U_RIGHT 1.
#define L0 1.
#define RHO 1.
#define MU 1.
#define OUTPUT_FREQ i++
#define TEND 1.
#define LEVEL 6
#define U_TOL 0.01
#define EPS 1.e-10

#ifndef ROT_X
  #define ROT_X 0
#endif
#ifndef ROT_Y
  #define ROT_Y 0
#endif
#ifndef ROT_Z
  #define ROT_Z 0
#endif
#if (!(ROT_X)&&!(ROT_Y)&&!(ROT_Z))
  #define ROT_Z 1
#endif
#ifndef THETA
  #define THETA (0.)
#endif

#define UX(X,Y,Z) ((sq(cos(THETA)) - sq(sin(THETA)))*(X) + 2*sin(THETA)*cos(THETA)*(Y))
#define UY(X,Y,Z) (-(sq(cos(THETA)) - sq(sin(THETA)))*(Y) + 2*sin(THETA)*cos(THETA)*(X))
#define UZ(X,Y,Z) 0.

#if ROT_Z
  #define MY_X x
  #define MY_Y y
  #define MY_Z z
  #define U_X(X,Y,Z) UX(X,Y,Z)
  #define U_Y(X,Y,Z) UY(X,Y,Z)
  #define U_Z(X,Y,Z) UZ(X,Y,Z)
#endif
#if ROT_Y
  #define MY_X y
  #define MY_Y z
  #define MY_Z x
  #define U_X(X,Y,Z) UZ(X,Y,Z)
  #define U_Y(X,Y,Z) UX(X,Y,Z)
  #define U_Z(X,Y,Z) UY(X,Y,Z)
#endif
#if ROT_X
  #define MY_X z
  #define MY_Y x
  #define MY_Z y
  #define U_X(X,Y,Z) UY(X,Y,Z)
  #define U_Y(X,Y,Z) UZ(X,Y,Z)
  #define U_Z(X,Y,Z) UX(X,Y,Z)
#endif



#include "eulerian_caps/grid/my_octree.h"
#include "eulerian_caps/navier-stokes/my_centered.h"
#include "eulerian_caps/capsule.h"
#include "eulerian_caps/elasticity.h"
#include "eulerian_caps/neo_hookean.h"
#include "view.h"

scalar my_gamma[];

void impose_u() {
  foreach() {
    u.x[] = U_X(MY_X, MY_Y, MY_Z);
    u.y[] = U_Y(MY_X, MY_Y, MY_Z);
    u.z[] = U_Z(MY_X, MY_Y, MY_Z);
  }
  boundary((scalar*){u});
}

int main() {
  size(L0);
  origin(-.5*L0, -.5*L0, -.5*L0);
  rho1 = RHO;
  rho2 = RHO;
  mu1 = MU;
  mu2 = MU;
  DT = 5.e-3;
  stokes = true;
  TOLERANCE = HUGE;
  run();
}

/* Boundary conditions:*/
u.n[left] = dirichlet(U_X(MY_X, MY_Y, MY_Z));
u.n[right] = dirichlet(U_X(MY_X, MY_Y, MY_Z));
u.n[top] = dirichlet(U_Y(MY_X, MY_Y, MY_Z));
u.n[bottom] = dirichlet(U_Y(MY_X, MY_Y, MY_Z));
u.n[front] = dirichlet(U_Z(MY_X, MY_Y, MY_Z));
u.n[back] = dirichlet(U_Z(MY_X, MY_Y, MY_Z));
u.t[left] = dirichlet(U_Y(MY_X, MY_Y, MY_Z));
u.t[right] = dirichlet(U_Y(MY_X, MY_Y, MY_Z));
u.t[top] = dirichlet(U_Z(MY_X, MY_Y, MY_Z));
u.t[bottom] = dirichlet(U_Z(MY_X, MY_Y, MY_Z));
u.t[front] = dirichlet(U_X(MY_X, MY_Y, MY_Z));
u.t[back] = dirichlet(U_X(MY_X, MY_Y, MY_Z));
u.r[left] = dirichlet(U_Z(MY_X, MY_Y, MY_Z));
u.r[right] = dirichlet(U_Z(MY_X, MY_Y, MY_Z));
u.r[top] = dirichlet(U_X(MY_X, MY_Y, MY_Z));
u.r[bottom] = dirichlet(U_X(MY_X, MY_Y, MY_Z));
u.r[front] = dirichlet(U_Y(MY_X, MY_Y, MY_Z));
u.r[back] = dirichlet(U_Y(MY_X, MY_Y, MY_Z));
p[top] = neumann(0.);
p[bottom] = neumann(0.);
p[right] = neumann(0.);
p[left] = neumann(0.);
p[front] = neumann(0.);
p[back] = neumann(0.);


event init (i = 0) {
  fraction(f, -MY_X*sin(THETA) + MY_Y*cos(THETA) + EPS);
  impose_u();
}

event impose_fraction (i++) {
  fraction(f, -MY_X*sin(THETA) + MY_Y*cos(THETA) + EPS);
}

event properties (i++) {
  fraction(f,  -MY_X*sin(THETA) + MY_Y*cos(THETA) + EPS);
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
  pseudo_t G_avg, G_int, Ts_avg, Ts_int;
  int nb_mb_cells = 0., nb_interface_cells = 0.;
  foreach_dimension() {
    G_avg.x.x = 0.;
    G_avg.x.y = 0.;
    G_avg.x.z = 0.;
    G_int.x.x = 0.;
    G_int.x.y = 0.;
    G_int.x.z = 0.;
    Ts_avg.x.x = 0.;
    Ts_avg.x.y = 0.;
    Ts_avg.x.z = 0.;
    Ts_int.x.x = 0.;
    Ts_int.x.y = 0.;
    Ts_int.x.z = 0.;
  }
  foreach() {
    my_gamma[] = GAMMA;
    if ((GAMMA)) {
      nb_mb_cells += 1;
      J_avg += J[];
      foreach_dimension() {
        G_avg.x.x += G.x.x[];
        G_avg.x.y += G.x.y[];
        G_avg.x.z += G.x.z[];
        Ts_avg.x.x += Ts.x.x[];
        Ts_avg.x.y += Ts.x.y[];
        Ts_avg.x.z += Ts.x.z[];
      }
      if (IS_INTERFACE_CELL(point,f)) {
        nb_interface_cells += 1;
        J_int += J[];
        foreach_dimension() {
          G_int.x.x += G.x.x[];
          G_int.x.y += G.x.y[];
          G_int.x.z += G.x.z[];
          Ts_int.x.x += Ts.x.x[];
          Ts_int.x.y += Ts.x.y[];
          Ts_int.x.z += Ts.x.z[];
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
      Ts_avg.x.x /= nb_mb_cells;
      Ts_avg.x.y /= nb_mb_cells;
      Ts_avg.x.z /= nb_mb_cells;
    }
  }
  if (nb_interface_cells) {
    J_int /= nb_interface_cells;
    foreach_dimension() {
      G_int.x.x /= nb_interface_cells;
      G_int.x.y /= nb_interface_cells;
      G_int.x.z /= nb_interface_cells;
      Ts_int.x.x /= nb_interface_cells;
      Ts_int.x.y /= nb_interface_cells;
      Ts_int.x.z /= nb_interface_cells;
    }
  }
  fprintf(stdout, "%.17g %.17g %d %d %d %.17g %.17g ", t, (THETA), (ROT_X), (ROT_Y), (ROT_Z), J_avg, J_int);
  foreach_dimension() {
    fprintf(stdout, "%.17g %.17g %.17g %.17g %.17g %.17g ", G_avg.x.x, G_int.x.x, G_avg.x.y, G_int.x.y, G_avg.x.z, G_int.x.z);
    fprintf(stdout, "%.17g %.17g %.17g %.17g %.17g %.17g ", Ts_avg.x.x, Ts_int.x.x, Ts_avg.x.y, Ts_int.x.y, Ts_avg.x.z, Ts_int.x.z);
  }
  fprintf(stdout,"\n");
  fflush(fout);
}

// event movie (i++) {
//   // view (fov = 20, theta=pi/10, phi=pi/5, psi=0, tx = 0., ty = -0., tz=-8., bg = {0.3,0.4,0.6});
//   char name[20];
//   sprintf(name, "J%d.png",i);
//   view (fov = 20, theta=pi/10, phi=pi/10, psi=0, bg = {0.3,0.4,0.6});
//   draw_vof ("f", lw=2);
//   squares("J",n = {0,0,1}, alpha = 1.e-12, linear=false);
//   squares("J",n = {0,1,0}, alpha = 1.e-12, linear=false);
//   squares("J",n = {1,0,0}, alpha = 1.e-12, linear=false);
//   cells (n = {1,0,0}, alpha = 1.e-12);
//   cells (n = {0,1,0}, alpha = 1.e-12);
//   cells (n = {0,0,1}, alpha = 1.e-12);
//   save (name);
//
//   sprintf(name, "uz%d.png",i);
//   view (fov = 20, theta=pi/10, phi=pi/10, psi=0, bg = {0.3,0.4,0.6});
//   // draw_vof ("f", lw=2);
//   squares("u.z",n = {0,0,1}, alpha = 1.e-12, linear=false);
//   squares("u.z",n = {0,1,0}, alpha = 1.e-12, linear=false);
//   squares("u.z",n = {1,0,0}, alpha = 1.e-12, linear=false);
//   // cells (n = {1,0,0}, alpha = 1.e-12);
//   cells (n = {0,1,0}, alpha = 1.e-12);
//   cells (n = {0,0,1}, alpha = 1.e-12);
//   save (name);
// }

event end (t = TEND) {
  return (1.);
}
