/**
In this file we compute the bending force on each Lagrangian node of the
membrane(s). The expression for the bending force is the following:
$$ \bm{f_b} = q\bm{n} \quad \text{with} \, q = \frac{dm}{dl} \quad \text{and}
\, m = E_b (\kappa - \kappa_R) $$
where $\kappa$ is the curvature, $\kappa_R$ is the reference curvature in the
unstressed configuration, and $l$ is a line element along the membrane.
*/

#ifndef E_B
  #define E_B 1.
#endif
#ifndef REF_CURV
  #define REF_CURV 0
#endif

  event init (i=0) {
    for(int k=0; k<mbs.nbmb; k++) {
      comp_curvature(&mbs.mb[k]);
      for(int i=0; i<mbs.mb[k].nlp; i++) {
        #if REF_CURV
          mbs.mb[k].nodes[i].ref_curv = mbs.mb[k].nodes[i].curv;
        #else
          mbs.mb[k].nodes[i].ref_curv = 0.;
        #endif
      }
    }
  }

void bending(lagMesh* mesh) {
  comp_curvature(mesh);
  coord* raw_bending = malloc(mesh->nlp*sizeof(coord));
  for(int i=0; i<mesh->nlp; i++) {
    double l[2];
    for(int j=0; j<2; j++)
      l[j] = mesh->edges[mesh->nodes[i].edge_ids[j]].length;
    double lavg = .5*(l[0] + l[1]);
    double curv[3];
    for(int j=0; j<3; j++) {
      int neighbor_index = (i-1+j) > 0 ? (i-1+j)%mesh->nlp : mesh->nlp - 1;
      curv[j] = mesh->nodes[neighbor_index].curv - mesh->nodes[neighbor_index].ref_curv;
    }
    double ddcurv = (curv[0]/(l[0]*lavg) - 2.*curv[1]/(l[0]*l[1])
      + curv[2]/(l[1]*lavg));
    foreach_dimension() raw_bending[i].x = -E_B*(ddcurv
      + .5*cube(curv[1]))*mesh->nodes[i].normal.x*lavg;
  }
  /** We apply a smoothing to the bending force, from Longuet-Higgins &
  Cokelet, Proc. R. Soc. Long., 1976; as advised by Pozrikidis in his book
  Computational hydrodynamics of capsules and biological cells (2010). */
  for(int i=0; i<mesh->nlp; i++) {
    // double coeff[5] = {-1./16., 1./4., 5./8., 1./4., -1./16.};
    double coeff[7] = {-1./32., 0., 9./32., 1./2., 9./32., 0., -1./32.};
    for(int j=0; j<7; j++)
      foreach_dimension()
        mesh->nodes[i].lagForce.x +=
          coeff[j]*raw_bending[(mesh->nlp + i - 3 + j)%mesh->nlp].x;
  }
  free(raw_bending);
}

// void bending(lagMesh* mesh) {
//   compute_lengths(mesh);
//   comp_curvature(mesh);
//   for(int i=0; i<mesh->nlp; i++) {
//     /** For each node, we compute the bending moment at the midpoint of each
//     connecting edge. The current and reference curvature at the midpoint of the
//     edges are simply the average of those at their nodes. */
//     double m[2]; // the bending moment
//     double l[2]; // the length of the edges
//     for(int j=0; j<2; j++) {
//       int edge_id = mesh->nodes[i].edge_ids[j];
//       int edge_nodes[2];
//       edge_nodes[0] = mesh->edges[edge_id].node_ids[0];
//       edge_nodes[1] = mesh->edges[edge_id].node_ids[1];
//       l[j] = mesh->edges[edge_id].length;
//       m[j] = .5*E_B*(mesh->nodes[edge_nodes[0]].curv +
//         mesh->nodes[edge_nodes[1]].curv - mesh->nodes[edge_nodes[0]].ref_curv -
//           mesh->nodes[edge_nodes[1]].ref_curv);
//     }
//     /** We then differentiate the bending moment to obtain q=dm/dl at the
//     considered node. */
//     double q = (m[1] - m[0])/(.5*(l[0] + l[1]));
//     foreach_dimension() mesh->nodes[i].lagForce.x += q*mesh->nodes[i].normal.x;
//   }
// }



event acceleration (i++) {
  for(int i=0; i<mbs.nbmb; i++) bending(&mbs.mb[i]);
}
