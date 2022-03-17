// #define LEVEL 5
#define MIN_LEVEL 3
#define MAX_LEVEL 7
#define NLP 31
#define RADIUS .125
#define L0 1.
#define MY_DT (.25*L0/(1 << MAX_LEVEL))
// #define MY_DT (.25*L0/(1 << level)) // Uncomment this line to see convergence in time
#define T_MAX 1.

double dt;
vector u[];

#include "grid/quadtree.h"
#include "lagrangian_caps/lag-mesh-2d.h"
#include "lagrangian_caps/reg-dirac.h"
#include "lagrangian_caps/view-ft.h"

int main(int argc, char* argv[]) {
  fprintf(stdout, "level\tavg_err\tmax_err\n");

  double t;
  for(int level=MIN_LEVEL; level <= MAX_LEVEL; level++) {
    dt = MY_DT;

    fprintf(stdout, "level=%d\n", level);
    N = 1 << level;
    origin(-.5*L0, -.5*L0);
    init_grid(N);

    lagMesh mb;
    initialize_circular_mb(&mb);
    coord ref_data[NLP];
    for(int i=0; i < NLP; i++){
      foreach_dimension() ref_data[i].x = mb.vertices[i].x;
    }

    t = 0.;
    int c = 0;
    while (t <= T_MAX) {
      foreach() {
        double x0 = x+L0/2.;
        double y0 = y+L0/2.;
        u.x[] = -2*sq(sin(pi*x0))*sin(pi*y0)*cos(pi*y0)*cos(pi*t/T_MAX);
        u.y[] = -2*sin(pi*x0)*cos(pi*x0)*sq(cos(pi*y0))*cos(pi*t/T_MAX);
      }
      boundary((scalar*){u});

      eul2lag(&mb);
      advect_lagMesh(&mb);

      if ((level == MAX_LEVEL) && (c%10 == 0)) {
        view(fov = 20, bg = {1,1,1});
        clear();
        if (level < 7) cells();
        squares("u.x", min = -2., max = 2.);
        draw_lag(&mb, lc = {1.,0.,0.}, vc = {1.,0.,0.});
        save("advect_caps.mp4");
      }

      t += dt;
      c++;
    }

    /**
    Compute the error
    */
    double avg_err, max_err;
    avg_err = 0.; max_err = -HUGE;
    for(int i=0; i < NLP; i++){
      double err = 0.;
      foreach_dimension() err += sq(ref_data[i].x - mb.vertices[i].x);
      err = sqrt(err);
      avg_err += err;
      if (err > max_err) max_err = err;
    }
    avg_err /= NLP;
    fprintf(stdout, "%d\t%g\t%g\n", level, avg_err, max_err);
  }

  return 0.;
}
