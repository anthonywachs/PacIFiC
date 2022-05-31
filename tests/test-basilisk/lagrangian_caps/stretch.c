// #define LEVEL 5
#define MIN_LEVEL 6
#define MAX_LEVEL 6
#define NLP 31
#define RADIUS .125
#define L0 1.
#define MY_DT (1.e-2*.25*L0/(1 << MAX_LEVEL))
// #define MY_DT (.25*L0/(1 << level)) // Uncomment this line to see convergence in time
#define T_MAX 1.
#define SD (.5*T_MAX)

#define MT ((t > .5*T_MAX) ? T_MAX-t : t)
#define X1 (fabs(x-L0*MT) <= .5*L0 ? (x-L0*MT) : \
  (x-L0*MT > .5*L0) ? (x-L0*MT-L0) : (x-L0*MT+L0))
#define Y1 (fabs(y-L0*MT) <= .5*L0 ? (y-L0*MT) : \
  (y-L0*MT > .5*L0) ? (y-L0*MT-L0) : (y-L0*MT+L0))
#define ALPHA ((fabs(X1) < 1.e-14) ? (Y1>0. ? pi/2. : 3*pi/2.) : atan2(Y1,X1))
#define X2 (fabs(x-.5*L0) <= .5*L0 ? (x-.5*L0) : \
  (x-.5*L0 > .5*L0) ? (x-.5*L0-L0) : (x-.5*L0+L0))
#define Y2 (fabs(y-.5*L0) <= .5*L0 ? (y-.5*L0) : \
  (y-.5*L0 > .5*L0) ? (y-.5*L0-L0) : (y-.5*L0+L0))
#define THETA ((fabs(X2) < 1.e-14) ? (Y2>0. ? pi/2. : 3*pi/2.) : atan2(Y2,X2))
#define RMAX (L0/3.)

double dt;
vector u[];
FILE* stf = NULL;
face vector a[];

#include "grid/quadtree.h"
#include "lagrangian_caps/lag-mesh-2d.h"
#include "lagrangian_caps/view-ft.h"

void output_stretch(lagMesh* mesh, FILE* file) {
  comp_mb_stretch(mesh);
  double avg_st, max_st, min_st;
  avg_st = 0.; max_st = -HUGE; min_st = HUGE;
  for(int i=0; i<mesh->nle; i++){
    avg_st += mesh->edges[i].st;
    if (mesh->edges[i].st > max_st) max_st = mesh->edges[i].st;
    else if (mesh->edges[i].st < min_st) min_st = mesh->edges[i].st;
  }
  avg_st /= mesh->nle;
  fprintf(file, "%g %g %g\n", avg_st, min_st, max_st);
}

int main(int argc, char* argv[]) {
  fprintf(stdout, "level\tavg_err\tmax_err\n");
  stf = fopen("stretch.txt","w");
  fprintf(stf, "avg_stretch min_stretch max_stretch\n");

  double t;
  for(int level=MIN_LEVEL; level <= MAX_LEVEL; level++) {
    dt = MY_DT;

    fprintf(stdout, "level=%d\n", level);
    N = 1 << level;
    origin(-.5*L0, -.5*L0);
    init_grid(N);
    periodic(left);
    periodic(top);

    lagMesh mb;
    initialize_circular_mb(&mb);
    coord ref_data[NLP];
    for(int i=0; i < NLP; i++){
      foreach_dimension() ref_data[i].x = mb.nodes[i].pos.x;
    }
    lagMesh mb2;
    initialize_circular_mb(&mb2);

    output_stretch(&mb, stf);

    t = 0.;
    int c = 0;
    while (t <= T_MAX) {
      foreach() {
        if (t < .5*T_MAX) {
          u.x[] = (L0 + cos(ALPHA)/(3));
          u.y[] = (L0 + sin(ALPHA)/(3));
        }
        else {
          u.x[] = -(L0 + cos(ALPHA)/(3));
          u.y[] = -(L0 + sin(ALPHA)/(3));
        }
      }
      boundary((scalar*){u});

      eul2lag(&mb);
      advect_lagMesh(&mb);
      output_stretch(&mb, stf);

      if (level == MAX_LEVEL && ((c%100) == 0)) {
        view(fov = 20, bg = {1,1,1});
        clear();
        if (level < 7) cells();
        draw_lag(&mb, lc = {1.,0.,0.}, vc = {1.,0.,0.});
        draw_lag(&mb2);
        save("stretch.mp4");
      }

      t += dt;
      c++;
    }


    if (level == MAX_LEVEL) {
      for(int i=0; i<20; i++) {
        view(fov = 20, bg = {1,1,1});
        clear();
        if (level < 7) cells();
        draw_lag(&mb, lc = {1.,0.,0.}, vc = {1.,0.,0.});
        draw_lag(&mb2);
        save("stretch.mp4");
      }
    }

    view(fov = 20, bg = {1,1,1});
    clear();
    if (level < 7) cells();
    draw_lag(&mb, lc = {1.,0.,0.}, vc = {1.,0.,0.});
    draw_lag(&mb2);
    char name[64];
    sprintf(name,"err_adv_lvl%d.png",level);
    save(name);

    output_stretch(&mb, stf);

    /**
    Compute the error
    */
    double avg_err, max_err;
    avg_err = 0.; max_err = -HUGE;
    for(int i=0; i < NLP; i++){
      double err = 0.;
      foreach_dimension() err += sq(ref_data[i].x - mb.nodes[i].pos.x);
      err = sqrt(err);
      avg_err += err;
      if (err > max_err) max_err = err;
    }
    avg_err /= NLP;
    fprintf(stdout, "%d\t%g\t%g\n", level, avg_err, max_err);
  }

  return 0.;
}
