/**
# Bending force for front-tracking membranes


*/

event acceleration (i++) {
  for(int i=0; i<mbs.nbmb; i++) {
    lagMesh* mesh = &(MB(i));
    comp_curvature(&MB(0));
    double curv = mesh->nodes[i].curv;
    double gcurv = mesh->nodes[i].gcurv;
    double lbcurv = laplace_beltrami(mesh, i, true);
    double bending_surface_force = 2*(2*curv*(sq(curv) - sq(gcurv)) + lbcurv);
    /** We now have to compute the area associated to each node */
    double area = 0;
    for(int j=0; j<MB(0).nodes[i].nb_triangles; j++) {
      int tid = &(MB(0).nodes[i].triangles[j]);
      if (is_obtuse_triangle(&MB(0), tid) {
        area += (is_obtuse_node(&MB(0), tid, i)) ? MB(0).triangles[tid].area/2 :
          MB(0).triangles[tid].area/4;
      }
      else {
        voronoi_area = 0;
        for(int k=0; k<3; k++) {
          int nid = &MB(0).triangles[tid].node_ids[k];
          if (nid != i) {
            int onid[2]; // onid for "opposite node ids"
            // find onid[0], onid[1]
            // compute their angle facing the relevant triangle
            // compute the cotangent of that angle
            // compute the squared length of [i:nid]
            // compute the voronoi area
          }
        }
        area += voronoi_area/16
      }
    }
    /** The bending force is ready to be added to the Lagrangian force of the
    considered node. */
    foreach_dimension()
      mesh->nodes[i].lagForce.x -=
        mesh->nodes[i].normal.x*bending_surface_force*area;
  }
}
