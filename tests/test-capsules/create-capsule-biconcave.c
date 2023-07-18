/**
# Capsule creation test: biconcave shape

In this case we simply create a biconcave capsule of radius 1 and ensure that
all the membrane nodes are at the right location. We don't test the connectivity
as it should be a necessary condition for the nodes position to be correct.

This file essentially tests most functions in mesh-toolbox.h and 
common-shapes.h.
*/

#define RADIUS 1.4
#define NCAPS 1
#define LAG_LEVEL 4

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "lagrangian_caps/lag-mesh.h"
#include "lagrangian_caps/common-shapes.h"

int main() {
    run();
}

event init (i = 0) {
    activate_spherical_capsule(&MB(0), radius=RADIUS, level=LAG_LEVEL, 
        shift={1., 2. -4.});
    for(int i=0; i<MB(0).nlp; i++) {
        foreach_dimension() fprintf(stderr, "%g ", MB(0).nodes[i].pos.x);
        fprintf(stderr, "\n");
    }
    exit(0);
}