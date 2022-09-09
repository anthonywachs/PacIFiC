# Spheres settling

This test case checks that the translational and rotational contact models with history of Grains3D leads to a static equilibrium, with velocities approaching the machine epsilon (i.e. less than 1.e-10).

## Set up
20 glued-spheres are randomly placed in a box and settle under gravity. Time integration is performed by the velocity verlet algorithm. This test runs on one processor.

## Success assessment
We check that the maximum translational velociy is lower than 1.e-10. In practice it approaches 1.e-15.
