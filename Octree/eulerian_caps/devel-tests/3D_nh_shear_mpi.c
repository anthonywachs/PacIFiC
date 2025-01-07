#define EPS 1.e-14
#ifndef SIZE
  #define SIZE 4.
#endif
#define SHEAR_RATE 1.
#define RADIUS 1.
#define MU 1.
#ifndef CA
  #define CA 0.2
#endif
#define E_S (MU*RADIUS*SHEAR_RATE/CA)
#define RHO 1.
#ifndef TEND
  #define TEND (2./SHEAR_RATE)
#endif
#ifndef LEVEL
  #define LEVEL 7
#endif
#ifndef START_SHEAR
  #define START_SHEAR 1
#endif
#ifndef DT_MAX
  #define DT_MAX 1.e-4
#endif
#ifndef U_TOL
  #define U_TOL 0.01
#endif
#ifndef OUTPUT_FREQ
  #define OUTPUT_FREQ 10
#endif

#include "eulerian_caps/grid/my_octree.h"
#include "eulerian_caps/navier-stokes/my_centered.h"
#include "eulerian_caps/capsule.h"
#include "eulerian_caps/elasticity.h"
#include "eulerian_caps/neo_hookean.h"
#include "view.h"

FILE* foutput = NULL;

int main() {
  L0 = SIZE;
  origin (-0.5*L0, -0.5*L0, -0.5*L0);
  stokes = true;
  TOLERANCE = 1.e-10;
  rho1 = RHO;
  rho2 = RHO;
  mu1 = MU;
  mu2 = MU;
  DT = DT_MAX;
  periodic(right);
  periodic(front);
  N = 1 << LEVEL;
  run();
}

u.n[bottom] = dirichlet(0);
u.n[top] = dirichlet(0);
u.t[bottom] = dirichlet(0.);
u.t[top] = dirichlet(0.);
u.r[top] = dirichlet(SHEAR_RATE*y);
u.r[bottom] =  dirichlet(SHEAR_RATE*y);

event init (t=0) {
  if (START_SHEAR) {
    foreach() {
      u.x[] = SHEAR_RATE*y;
      u.y[] = 0.;
      u.z[] = 0.;
    }
  }
  fraction(f, 1. - sq(x/(RADIUS)) - sq(y/(RADIUS)) - sq(z/(RADIUS)));
  foutput  = fopen("output.txt","w");
}

scalar noise_gamma[];
event adapt (i++) {
  foreach() {
    if (GAMMA) noise_gamma[] = noise();
    else noise_gamma[] = 0.;
  }
  boundary({noise_gamma});
  adapt_wavelet ({noise_gamma,u}, (double[]){1.e-6, U_TOL, U_TOL, U_TOL}, maxlevel = LEVEL);
}

event logfile (i+=OUTPUT_FREQ) {
  double rmax = -HUGE, rmin = HUGE ;
  foreach (reduction(max:rmax) reduction(min:rmin))
    if (IS_INTERFACE_CELL(point,f) && (fabs(z) < Delta)) {
      coord p;
      coord n = mycs (point, f);
      double alpha = plane_alpha (f[], n);
      plane_area_center (n, alpha, &p);
      // double rad  = sqrt(sq(x + Delta*p.x) + sq(y + Delta*p.y) + sq(z + Delta*p.z));
      double rad  = sqrt(sq(x + Delta*p.x) + sq(y + Delta*p.y));
      if (rad > rmax)
	     rmax = rad;
      if (rad < rmin)
	     rmin = rad;
    }
  double D = (rmax - rmin)/(rmax + rmin);
  fprintf (foutput, "%g %g %g %g\n", t, rmin, rmax, D);
}

event pictures (i+=OUTPUT_FREQ) {
  char name[100];

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells (n = {0,0,1}, alpha = 1.e-12);
  squares ("u.x", n = {0,0,1}, alpha = 1.e-12, linear = false);
  draw_vof ("f");
  box (notics = true);
  sprintf(name,"ux_vof%d.png",i/OUTPUT_FREQ);
  // save (name);
  save("ux_vof.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells (n = {0,0,1}, alpha = 1.e-12);
  squares ("u.x", n = {0,0,1}, alpha = 1.e-12, linear = false);
  box (notics = true);
  sprintf(name,"ux%d.png",i/OUTPUT_FREQ);
  // save (name);
  save("ux.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells (n = {0,0,1}, alpha = 1.e-12);
  squares ("u.y", n = {0,0,1}, alpha = 1.e-12, linear = false);
  box (notics = true);
  sprintf(name,"uy%d.png",i/OUTPUT_FREQ);
  // save (name);
  save("uy.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells (n = {0,0,1}, alpha = 1.e-12);
  squares ("u.z", n = {0,0,1}, alpha = 1.e-12, linear = false);
  box (notics = true);
  sprintf(name,"uz%d.png",i/OUTPUT_FREQ);
  // save (name);
  save("uz.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells (n = {0,0,1}, alpha = 1.e-12);
  squares ("centered_ae.x", n = {0,0,1}, alpha = 1.e-12, linear = false);
  box (notics = true);
  sprintf(name,"ax%d.png",i/OUTPUT_FREQ);
  // save (name);
  save("ax.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells (n = {0,0,1}, alpha = 1.e-12);
  squares ("centered_ae.y", n = {0,0,1}, alpha = 1.e-12, linear = false);
  box (notics = true);
  sprintf(name,"ay%d.png",i/OUTPUT_FREQ);
  // save (name);
  save("ay.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells (n = {0,0,1}, alpha = 1.e-12);
  squares ("centered_ae.z", n = {0,0,1}, alpha = 1.e-12, linear = false);
  box (notics = true);
  sprintf(name,"az%d.png",i/OUTPUT_FREQ);
  // save (name);
  save("az.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells (n = {0,0,1}, alpha = 1.e-12);
  squares ("J", n = {0,0,1}, alpha = 1.e-12, linear = false);
  box (notics = true);
  sprintf(name,"J%d.png",i/OUTPUT_FREQ);
  // save (name);
  save("J.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells (n = {0,0,1}, alpha = 1.e-12);
  squares ("G.x.x", n = {0,0,1}, alpha = 1.e-12, linear = false);
  box (notics = true);
  sprintf(name,"Gxx%d.png",i/OUTPUT_FREQ);
  // save (name);
  save("Gxx.mp4");

  // view (fov = 20, camera = "front", bg = {1,1,1});
  // cells (n = {0,0,1}, alpha = 1.e-12);
  // squares ("G.x.y", n = {0,0,1}, alpha = 1.e-12, linear = false);
  // box (notics = true);
  // sprintf(name,"Gxy%d.png",i/OUTPUT_FREQ);
  // save (name);

  // view (fov = 20, camera = "front", bg = {1,1,1});
  // cells (n = {0,0,1}, alpha = 1.e-12);
  // squares ("G.x.z", n = {0,0,1}, alpha = 1.e-12, linear = false);
  // box (notics = true);
  // sprintf(name,"Gxz%d.png",i/OUTPUT_FREQ);
  // save (name);

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells (n = {0,0,1}, alpha = 1.e-12);
  squares ("G.y.y", n = {0,0,1}, alpha = 1.e-12, linear = false);
  box (notics = true);
  sprintf(name,"Gyy%d.png",i/OUTPUT_FREQ);
  // save (name);
  save("Gyy.mp4");

  // view (fov = 20, camera = "front", bg = {1,1,1});
  // cells (n = {0,0,1}, alpha = 1.e-12);
  // squares ("G.y.z", n = {0,0,1}, alpha = 1.e-12, linear = false);
  // box (notics = true);
  // sprintf(name,"Gyz%d.png",i/OUTPUT_FREQ);
  // save (name);

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells (n = {0,0,1}, alpha = 1.e-12);
  squares ("G.z.z", n = {0,0,1}, alpha = 1.e-12, linear = false);
  box (notics = true);
  sprintf(name,"Gzz%d.png",i/OUTPUT_FREQ);
  // save (name);
  save("Gzz.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells (n = {0,0,1}, alpha = 1.e-12);
  squares ("Ts.x.x", n = {0,0,1}, alpha = 1.e-12, linear = false);
  box (notics = true);
  sprintf(name,"Txx%d.png",i/OUTPUT_FREQ);
  // save (name);
  save("Txx.mp4");

  // view (fov = 20, camera = "front", bg = {1,1,1});
  // cells (n = {0,0,1}, alpha = 1.e-12);
  // squares ("Ts.x.y", n = {0,0,1}, alpha = 1.e-12, linear = false);
  // box (notics = true);
  // sprintf(name,"Txy%d.png",i/OUTPUT_FREQ);
  // save (name);

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells (n = {0,0,1}, alpha = 1.e-12);
  squares ("extended_n.x", n = {0,0,1}, alpha = 1.e-12, linear = false);
  sprintf(name,"nx%d.png",i/OUTPUT_FREQ);
  // save (name);
  save("nx.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells (n = {0,0,1}, alpha = 1.e-12);
  squares ("extended_n.y", n = {0,0,1}, alpha = 1.e-12, linear = false);
  sprintf(name,"ny%d.png",i/OUTPUT_FREQ);
  // save (name);
  save("ny.mp4");

  view (fov = 20, camera = "front", bg = {1,1,1});
  cells (n = {0,0,1}, alpha = 1.e-12);
  squares ("extended_n.z", n = {0,0,1}, alpha = 1.e-12, linear = false);
  sprintf(name,"nz%d.png",i/OUTPUT_FREQ);
  // save (name);
  save("nz.mp4");
}

event end (t = TEND) {
  fclose(foutput);
  return 1.;
}
