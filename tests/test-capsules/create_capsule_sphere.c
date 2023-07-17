#define RADIUS 1.
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