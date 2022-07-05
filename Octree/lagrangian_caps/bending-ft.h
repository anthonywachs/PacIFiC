/**
# Bending force for front-tracking membranes


*/
#ifndef E_B
  #define E_B 1.
#endif
#ifndef REF_CURV
  #define REF_CURV 1
#endif
#define cot(x) (cos(x)/sin(x))

#include "curvature-ft.h"

/**
The function below computes the nodal area of node $i$.
*/
double compute_node_area(lagMesh* mesh, int i) {
  double area = 0.;
  for(int j=0; j<mesh->nodes[i].nb_triangles; j++) {
    int tid = mesh->nodes[i].triangle_ids[j];
    if (is_obtuse_triangle(mesh, tid)) {
      area += (is_obtuse_node(mesh, tid, i)) ? mesh->triangles[tid].area/2 :
        mesh->triangles[tid].area/4;
    }
    else {
      double voronoi_area = 0.;
      for(int k=0; k<3; k++) {
        int nid = mesh->triangles[tid].node_ids[k];
        if (nid != i) {
          int eid = -1; // eid for "edge id", connecting i and nid
          for (int l=0; l<3; l++) {
            int teid = mesh->triangles[tid].edge_ids[l]; // temporary edge id
            if ((mesh->edges[teid].node_ids[0] == i ||
              mesh->edges[teid].node_ids[1] == i) &&
              (mesh->edges[teid].node_ids[0] == nid ||
              mesh->edges[teid].node_ids[1] == nid)) eid = teid;
          }
          int onid[2]; // onid for "opposite node ids"
          // find onid[0], onid[1]
          for(int l=0; l<2; l++) {
            // tneid for "triangle neighboring eid"
            int tneid = mesh->edges[eid].triangle_ids[l];
            for(int m=0; m<3; m++)
              if (mesh->triangles[tneid].node_ids[m] != i &&
                mesh->triangles[tneid].node_ids[m] != nid)
                onid[l] = mesh->triangles[tneid].node_ids[m];
          }
          // compute their angle facing the relevant triangle
          double theta[2];
          for(int l=0; l<2; l++) {
            theta[l] = comp_angle(&(mesh->nodes[onid[l]]), &(mesh->nodes[i]),
              &(mesh->nodes[nid]));
          }
          // compute the squared length of [i:nid]
          double edge_length = 0.;
          foreach_dimension()
            edge_length += sq(mesh->nodes[i].pos.x - mesh->nodes[nid].pos.x);
          voronoi_area += (cot(theta[0]) + cot(theta[1]))*edge_length;
        }
      }
      area += voronoi_area/16.;
    }
  }
  return area;
}

event acceleration (i++) {
  for(int i=0; i<mbs.nbmb; i++) {
    lagMesh* mesh = &(MB(i));
    comp_curvature(mesh);
    for(int j=0; j<mesh->nlp; j++) {
      double curv = mesh->nodes[j].curv;
      double rcurv = mesh->nodes[i].ref_curv;
      double gcurv = mesh->nodes[j].gcurv;
      double lbcurv = laplace_beltrami(mesh, j, true);
      double bending_surface_force =
        2*E_B*(2*(curv - rcurv)*(sq(curv) - gcurv + rcurv*curv) + lbcurv);
      /** We now have to compute the area associated to each node */
      double area = compute_node_area(mesh, j);
      /** The bending force is ready to be added to the Lagrangian force of the
      considered node. */
      foreach_dimension()
        mesh->nodes[j].lagForce.x -=
          mesh->nodes[j].normal.x*bending_surface_force*area;
    }
  }
}

event init (i = 0) {
  for(int i=0; i<mbs.nbmb; i++) {
    lagMesh* mesh = &(MB(i));
    #if REF_CURV
    comp_curvature(mesh);
    #endif
    for(int j=0; j<mesh->nlp; j++) {
      #if REF_CURV
      mesh->nodes[j].ref_curv = mesh->nodes[j].curv;
      #else
      mesh->nodes[j].ref_curv = 0.;
      #endif
    }
  }
}
