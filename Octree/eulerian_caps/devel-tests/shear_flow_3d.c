/* Definition of numerical parameters and global variables **/
#define DT_MAX 1.e-3
#define LEVEL 7
#define U_TOL 0.01
#define OUTPUT_FREQ 10
#define RENORMALIZE_INTERFACIAL_NORMALS 1
#define RENORMALIZE_EXTENDED_NORMALS 1
#define USE_HEIGHT_FUNCTIONS 1
#define EXTEND_MB_ATTRIBUTES 1
#define SHARP_DIRAC 1
#define MOVIE_FIELDS u.x,u.y,u.z,B.x.x,B.y.y,B.z.z,Ts.x.x,Ts.y.y,Ts.z.z,\
centered_ae.x,centered_ae.y,aen,extended_n.x,extended_n.y,extended_n.z
// #define MOVIE_FIELDS u.x,u.y,u.z,Bf.x.x,Bf.x.y, Bf.x.z, Bf.y.z,Bf.y.y,Bf.z.z,Ts.x.x,Ts.y.y,Ts.z.z,centered_ae.x,centered_ae.y,aen,extended_n.x,extended_n.y,extended_n.z,BfP.x.x,BfP.x.y,BfP.x.z, BfP.y.z,BfP.y.y,BfP.z.z,PBfP.x.x,PBfP.x.y,PBfP.x.z, PBfP.y.z,PBfP.y.y,PBfP.z.z, Proj.x.x,Proj.x.y,Proj.x.z,Proj.y.y,Proj.y.z,Proj.z.z,Jf,PBfP.y.x,PBfP.y.z,PBfP.z.x,PBfP.z.y
FILE* output_file = NULL;

/* Definition of physical parameters of the simulation **/
#define L0 4.
#define RHO 1.
#define MU 1.
#define CAPILLARY 0.2
#define SHEAR_RATE 1.
#define RADIUS (1. + 1.e-8)
#define E_S ((MU*SHEAR_RATE*RADIUS)/(CAPILLARY))
// #define TEND (2./SHEAR_RATE)
#define TEND (20.*DT_MAX)


/* We use an octree grid, the centered Navier-Stokes solver, and the
Eulerian capsule solver, with a Neo-Hookean law for the elastic membrane (but
with no bending stresses). Additionally, we output movies using bview.**/
#include "eulerian_caps/grid/my_octree.h"
#include "eulerian_caps/navier-stokes/my_centered.h"
#include "eulerian_caps/capsule.h"
#include "eulerian_caps/elasticity.h"
#include "eulerian_caps/neo_hookean.h"
#include "view.h"

int main(int argc, char* argv[]) {
  init_grid(1 << LEVEL);
  origin(-.5*L0, -.5*L0, -.5*L0);
  periodic(left);
  periodic(front);
  stokes = true;
  TOLERANCE = 1.e-10;
  DT = DT_MAX;
  rho1 = RHO, rho2 = RHO;
  mu1 = MU, mu2 = MU;
  run();
}

/* Definition of boundary conditions **/
u.n[top] = dirichlet(0.);
u.n[bottom] = dirichlet(0.);
u.t[top] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);
u.r[top] = dirichlet(SHEAR_RATE*y);
u.r[bottom] = dirichlet(SHEAR_RATE*y);


/** Define the spherical coordinates for the analytical normals at t = 0 */
#define RAD (sqrt(sq(x) + sq(y) + sq(z)) > 1.e-30 ? \
  sqrt(sq(x) + sq(y) + sq(z)) : 0.)
#if dimension > 2
  #define THETA (fabs(z) > 1.e-30 ? atan2(sqrt(sq(x) + sq(y)),z) : pi/2.)
#else
 #define THETA (pi/2.)
#endif
#define PHI (fabs(x) > 1.e-30 ? atan2(y,x) : ((y > 0 ? (pi/2.) : (3*pi/2.))))

event init (i = 0) {
  foreach() {
    u.x[] = SHEAR_RATE*y;
    u.y[] = 0.;
    u.z[] = 0.;
    extended_n.x[] = cos(PHI)*sin(THETA);
    extended_n.y[] = sin(PHI)*sin(THETA);
    extended_n.z[] = cos(THETA);
  }
  fraction(f, sq(RADIUS) - sq(x) - sq(y) - sq(z));
  output_file = fopen("output.txt","w");
}

/* We adapt according to both the velocity field and around the interface **/
event adapt (i++) {
  scalar noise_gamma[];
  foreach() {
    if (GAMMA) noise_gamma[] = noise();
    else noise_gamma[] = 0.;
  }
  boundary({noise_gamma});
  adapt_wavelet({noise_gamma, u}, (double[]){1.e-6, U_TOL, U_TOL, U_TOL},  maxlevel = LEVEL);
}

/* We output the Taylor deformation parameter to quantitatively compare our
method to the literature, in particular to the Boundary Element Method
(Pozrikidis, JFM, 1995), the Immersed Boundary Method (Eggleton and Popel, PoF,
1998) and to the first implementation of this fully Eulerian method (Ii et al.,
Comm. Comp. Phys., 2012). We follow the post-processing method of the first two
references, i.e. we compute the minimum and maximum radii in the shear plane
z = 0 only. **/
event output( i+= OUTPUT_FREQ) {
  double rmin, rmax, D;
  rmin = HUGE, rmax = -HUGE;
  foreach(reduction(min:rmin) reduction(max:rmax)) {
    if ((fabs(z) < Delta) && (f[] > 0) && (f[] < 1)) {
      coord p;
      coord n = mycs (point, f);
      double alpha = plane_alpha (f[], n);
      plane_area_center (n, alpha, &p);
      double radius;
      radius = sqrt(sq(x + p.x*Delta) + sq(y + p.y*Delta) + sq(z + p.z*Delta));
      if (radius < rmin) rmin = radius;
      if (radius > rmax) rmax = radius;
    }
  }
  D = (rmax - rmin)/(rmax + rmin);
  fprintf(output_file, "%.17g %.17g %.17g %.17g\n", t, rmin, rmax, D);
  fflush(output_file);
}

event movie ( i+= OUTPUT_FREQ) {
  scalar * movie_fields = {MOVIE_FIELDS};
  char movie_name[64];
  for (scalar s in movie_fields) {
    view(fov = 20, camera = "front", bg = {1,1,1});
    clear();
    squares(s.name, linear = false);
    cells(n = {0,0,1}, alpha = 1.e-12);
    sprintf(movie_name, "%s_%d.png", s.name, i/OUTPUT_FREQ);
    save(movie_name);
  }
  view(fov = 20, camera = "front", bg = {1,1,1});
  clear();
  draw_vof("f", lw = 2);
  squares("u.x", linear = false);
  cells(n = {0,0,1}, alpha = 1.e-12);
  sprintf(movie_name, "ux_cells_%d.png", i/OUTPUT_FREQ);
  save(movie_name);
}

event end (t = TEND) {
  fclose(output_file);
  return 1.;
}
