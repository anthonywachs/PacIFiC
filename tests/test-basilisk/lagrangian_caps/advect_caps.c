// #define LEVEL 5
#define MIN_LEVEL 5
#define MAX_LEVEL 7
#define NLP 47
#define RADIUS .125
#define L0 1.
#define MY_DT (.25*L0/(1 << MAX_LEVEL))
// #define MY_DT (.25*L0/(1 << level)) // Uncomment this line to see convergence in time
#define T_MAX 1.
#define PERIODIC 1

double dt;
vector u[];
FILE* foutput = NULL;

#include "grid/quadtree.h"
#include "lagrangian_caps/lag-mesh-2d.h"
#include "lagrangian_caps/reg-dirac.h"
#include "lagrangian_caps/view-ft.h"

int main(int argc, char* argv[]) {
  char name[64];
  char bc_type[16];
  #if PERIODIC
    sprintf(bc_type,"periodic");
    foutput = fopen("output_periodic.txt","w");
  #else
    sprintf(bc_type,"non_periodic");
    foutput = fopen("output_non_periodic.txt","w");
  #endif
  fprintf(foutput, "level\tavg_err\tmax_err\n");

  double t;
  for(int level=MIN_LEVEL; level <= MAX_LEVEL; level++) {
    dt = MY_DT;

    fprintf(foutput, "level=%d\n", level);
    N = 1 << level;
    origin(-.5*L0, -.5*L0);
    init_grid(N);
    #if PERIODIC
      periodic(left);
      periodic(top);
    #endif

    lagMesh mb;
    initialize_circular_mb(&mb);
    lagMesh mb2;
    initialize_circular_mb(&mb2);
    #if PERIODIC
      for(int i=0; i < mb.nlp; i++)
        foreach_dimension() mb.vertices[i].x += .501*L0;
      for(int i=0; i < mb2.nlp; i++)
        foreach_dimension() mb2.vertices[i].x += .501*L0;
    #else
      for(int i=0; i < mb.nlp; i++)
        foreach_dimension() mb.vertices[i].x += .001*L0;
      for(int i=0; i < mb2.nlp; i++)
        foreach_dimension() mb2.vertices[i].x += .001*L0;
    #endif
    correct_lag_pos(&mb);
    correct_lag_pos(&mb2);
    coord ref_data[NLP];
    for(int i=0; i < NLP; i++){
      foreach_dimension() ref_data[i].x = mb.vertices[i].x;
    }

    if (level == 6) {
      view(fov = 20, bg = {1,1,1});
      clear();
      vectors("u",scale = L0/(1 << level)/20.);
      draw_lag(&mb, lc = {1.,0.,0.}, vc = {1.,0.,0.});
      draw_lag(&mb2);
      sprintf(name, "advect_caps_start_%s.png", bc_type);
      save(name);
    }

    t = 0.;
    int c = 0;
    while (t <= T_MAX) {
      foreach() {
        #if PERIODIC
          double x0 = x + 0.001*L0;
          double y0 = y + 0.001*L0;
        #else
          double x0 = x+1.002*L0/2.;
          double y0 = y+1.002*L0/2.;
        #endif
        u.x[] = -2*sq(sin(pi*x0))*sin(pi*y0)*cos(pi*y0)*cos(pi*t/T_MAX);
        u.y[] = -2*sin(pi*x0)*cos(pi*x0)*sq(cos(pi*y0))*cos(pi*t/T_MAX);
      }
      boundary((scalar*){u});

      eul2lag(&mb);
      advect_lagMesh(&mb);

      if ((level == 6)) {
        view(fov = 20, bg = {1,1,1});
        clear();
        // if (level < 7) cells();
        // squares("u.x", min = -2., max = 2.);
        vectors("u",scale = L0/(1 << level)/20.);
        draw_lag(&mb, lc = {1.,0.,0.}, vc = {1.,0.,0.});
        draw_lag(&mb2);
        sprintf(name, "advect_caps_%s.mp4", bc_type);
        save(name);
      }

      t += dt;
      c++;
    }

    if(level == 6) {
      for(int i=0; i<5; i++) {
        view(fov = 20, bg = {1,1,1});
        clear();
        vectors("u",scale = L0/(1 << level)/20.);
        draw_lag(&mb, lc = {1.,0.,0.}, vc = {1.,0.,0.});
        draw_lag(&mb2);
        sprintf(name, "advect_caps_%s.mp4", bc_type);
        save(name);
      }
      view(fov = 20, bg = {1,1,1});
      clear();
      vectors("u",scale = L0/(1 << level)/20.);
      draw_lag(&mb, lc = {1.,0.,0.}, vc = {1.,0.,0.});
      draw_lag(&mb2);
      sprintf(name, "advect_caps_end_%s.mp4", bc_type);
      save(name);
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
    fprintf(foutput, "%d\t%g\t%g\n", level, avg_err, max_err);
  }

  return 0.;
}
