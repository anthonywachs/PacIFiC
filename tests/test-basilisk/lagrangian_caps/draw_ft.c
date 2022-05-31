#define LEVEL 4
#define NLP 31
#define RADIUS .25

#include "grid/quadtree.h"
#include "lagrangian_caps/lag-mesh-2d.h"
#include "lagrangian_caps/view-ft.h"


int main(int argc, char* argv[]) {
  lagMesh mb_mesh;
  initialize_circular_mb(&mb_mesh);

  origin(-.5*L0, -.5*L0);
  N = 1 << LEVEL;
  init_grid(N);

  view(fov = 20, bg = {1,1,1});
  clear();
  cells();
  draw_lag(&mb_mesh, lc = {1.,0.,0.}, vc = {1.,0.,0.});
  save("draw_ft.png");

  return 0.;
}
