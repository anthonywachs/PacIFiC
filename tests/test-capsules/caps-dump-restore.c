/**
# Test case for the dump and restore of a capsule

*/
#define L0 1.
#define RADIUS .2
#define LEVEL 5
#define LAG_LEVEL 4
#define NCAPS 2

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "lagrangian_caps/capsule-ft.h"
#include "lagrangian_caps/neo-hookean-ft.h"
#include "lagrangian_caps/bending-ft.h"
#include "lagrangian_caps/common-shapes-ft.h"

int main(int argc, char* argv[]) {
  origin(-.5*L0, -.5*L0, -.5*L0);
  N = 1 << LEVEL;
  run();
}

event init (i = 0) {
  activate_spherical_capsule(&CAPS(0), radius = RADIUS, 
    level = LAG_LEVEL);

  comp_volume(&CAPS(0));
  CAPS(0).initial_volume = CAPS(0).volume;

  double c0, c1, c2;
  c0 = 0.2072; c1 = 2.0026; c2 = -1.1228;
  double radius = RADIUS;
  for(int i=0; i<CAPS(0).nln; i++) {
    double rho = sqrt(sq(CAPS(0).nodes[i].pos.x) +
    sq(CAPS(0).nodes[i].pos.z))/radius;
    rho = (rho > 1) ? 1 : rho;
    int sign = (CAPS(0).nodes[i].pos.y > 0.) ? 1 : -1;
    CAPS(0).nodes[i].pos.y = sign*.5*radius*sqrt(1 - sq(rho))*
    (c0 + c1*sq(rho) + c2*sq(sq(rho)));
  }
  correct_lag_pos(&CAPS(0));
  comp_normals(&CAPS(0));
  comp_centroid(&CAPS(0));
  comp_volume(&CAPS(0));

  FILE* file = fopen("caps0.dump", "w");
  dump_lagmesh(file, &CAPS(0));
  fclose(file);
  file = fopen("caps0.dump", "r");
  restore_lagmesh(file, &CAPS(1));
  fclose(file);
  file = fopen("caps1.dump", "w");
  dump_lagmesh(file, &CAPS(1));
  fclose(file);
  dump_capsules(fp=stderr);
}

