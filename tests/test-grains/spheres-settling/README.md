# Spheres settling

This test case checks that the translational and rotational contact models with history of Grains3D leads to a static equilibrium, with velocities approaching the machine epsilon (i.e. less than 1.e-10).

## Set up
20 spheres which initial positions are defined in Grains/Init/initial_position.result are placed in a box and settle under gravity. Time integration is performed by the velocity verlet algorithm. No paraview output is generated. This test runs on one processor.

## Success assessment
1) The first test outputs the maximum translational velociy, which must be lower than 1.e-10.
2) Then, all the files located in Reference_results are compared, which include the binary restart files.

## Remarks
Any change in the output files would result in a fail of this test case, regardless of the ability of the contact model with history to lead to a static equilibrium.
Hence the first success assessement: if it succeeds but the overall test case fails, it is probably due to an updated formatting of the output files.
