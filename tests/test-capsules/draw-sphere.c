/**
# Capsule creation test

In this case we simply create a spherical capsule of radius 1 and ensure that
all the membrane nodes are at the right location. We don't test the connectivity
as it should be a necessary condition for the nodes position to be correct.

This file essentially tests most functions in mesh-toolbox.h and 
common-shapes.h.
*/

#define L0 1.
#define RADIUS .25
#define NCAPS 1
#define LAG_LEVEL 4

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "lagrangian_caps/capsule-ft.h"
#include "lagrangian_caps/common-shapes-ft.h"
#include "lagrangian_caps/view-ft.h"

int main() {
  origin(-.5*L0, -.5*L0, -.5*L0);
  run();
}

event init (i = 0) {
  activate_spherical_capsule(&CAPS(0), radius=RADIUS, level=LAG_LEVEL);
  foreach() u.x[] = sin(2*pi*x*y/L0);

  view(fov = 18.9, bg = {1,1,1});
  clear();
  draw_lags(lw = .5, ns = 3, facets = true);
  squares("u.x", n = {0,0,1});
  cells(n = {0,0,1});
  save(fp=stderr);

  exit(0);
}