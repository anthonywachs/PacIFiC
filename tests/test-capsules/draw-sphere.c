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
#include "lagrangian_caps/lag-mesh.h"
#include "lagrangian_caps/common-shapes.h"
#include "lagrangian_caps/view-ft.h"

int main() {
    origin(-.5*L0, -.5*L0, -.5*L0);
    run();
}

event init (i = 0) {
    activate_spherical_capsule(&MB(0), radius=RADIUS, level=LAG_LEVEL);

    view(fov = 18.9, bg = {1,1,1});
    clear();
    draw_lag(&MB(0), lw = .5, ns = 3, facets = true);
    save(fp=stderr);

    exit(0);
}