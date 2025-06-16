#define EPS 1.e-14
#define SIZE 4.
#define RADIUS (1. + sqrt(2)*1.e-3)
#define MU 1.
#define RHO 1.
#define LEVEL 9
#define EXTENSION_OUTPUT 1

vector extended_n[];
scalar q[];

static FILE * output_file[2];
int file_nb;
#define OUTPUT_FILE output_file[file_nb]

#include "grid/multigrid.h"
#include "eulerian_caps/navier-stokes/my_centered.h"
#include "eulerian_caps/capsule.h"
#include "view.h"
#include "eulerian_caps/normal_extension.h"

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
  char buf[0x100];
  if (mpi_rank == 0)
    for (int k=0; k<2; k++) {
      snprintf(buf, sizeof(buf), "q_extension_%d.csv", k+1);
      output_file[k] = fopen(buf,"w");
    }
}

event init_q (i = 1; i+=2) {
  foreach() {
    if (GAMMA) {
      if ((fabs(f[1,0] - f[]) > .5) || (fabs(f[-1,0] - f[]) > .5) || (fabs(f[0,1] - f[]) > .5) || (fabs(f[0,-1] - f[]) > .5)) {
        double my_alpha;
        if (fabs(x)>SEPS) my_alpha = atan2(y,x);
        else my_alpha = y > 0 ? pi/2 : -pi/2;
        q[] = cos(my_alpha);
      }
      else q[] = 0;
    }
    else q[] = -1.;
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
        if (i == 4) {
          extended_n.x[] = - grad_wide_caps.x[]/ng_wide_caps[];
          extended_n.y[] = - grad_wide_caps.y[]/ng_wide_caps[];
        }
      }
    }
    file_nb = i/2-1;
    normal_scalar_extension(q);
  }
}

event picture (i++) {
  view (fov = 20, camera = "front", bg = {1,1,1});
  squares ("q", linear = false);
  draw_vof ("f", lw = 2);
  box (notics = true);
  if (i==1) save ("q_init.png");
  if (i==2) save ("q_extended_1.png");
  if (i==4) save ("q_extended_2.png");
}

event end (i = 4) {
  if (mpi_rank == 0) for (int k=0; k<2; k++) fclose(output_file[k]);
  return 1.;
}
