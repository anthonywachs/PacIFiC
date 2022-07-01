/** Bending force for front-tracking */

event acceleration (i++) {
  for(int i=0; i<mbs.nbmb; i++) {
    lagMesh* mesh = &(MB(i));
    comp_curvature(&MB(0));
    double curv = mesh->nodes[i].curv;
    double gcurv = mesh->nodes[i].gcurv;
    double lbcurv = laplace_beltrami(mesh, i, true);
    double bending_force = 2*(2*curv*(sq(curv) - sq(gcurv)) + lbcurv);
    foreach_dimension()
      mesh->nodes[i].lagForce.x -= mesh->nodes[i].normal.x*bending_force;
  }
}
