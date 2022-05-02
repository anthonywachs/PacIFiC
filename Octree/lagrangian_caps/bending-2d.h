/** In this file we compute the bending force on each Lagrangian node of the
membrane(s). The expression for the bending force is the following:
$$ \bm{f_b} = q\bm{n} \quad \text{with} \, q = \frac{dm}{dl} \quad \text{and}
\, m = E_b (\kappa - \kappa_R) $$
where $\kappa$ is the curvature, $\kappa_R$ is the reference curvature in the
unstressed configuration, and $l$ is a line element along the membrane.
*/

#ifndef E_B
  #define E_B 1.
#endif

event init (i=0) {
  for(int k=0; k<mbs.nbmb; k++) {
    comp_curvature(&mbs.mb[k]);
    for(int i=0; i<mbs.mb[k].nlp; i++) {
      mbs.mb[k].nodes[i].ref_curv = mbs.mb[k].nodes[i].curv;
    }
  }
}

void bending(lagMesh* mesh) {
  comp_curvature(mesh);
  for(int i=0; i<mesh->nlp; i++) {
    /** For each node, we compute the bending moment at the midpoint of each
    connecting edge. The current and reference curvature at the midpoint of the
    edges are simply the average of those at their nodes. */
    double m[3]; // the bending moment
    double q[2]; // the transverse shear tension
    double l[2]; // the length of the edges
    for(int j=0; j<2; j++) {
      int edge_id = mesh->nodes[i].edge_ids[j];
      int edge_nodes[2];
      edge_nodes[0] = mesh->edges[edge_id].vertex_ids[0];
      edge_nodes[1] = mesh->edges[edge_id].vertex_ids[1];
      l[j] = mesh->edges[edge_id].length;
      m[0] = E_B*(mesh->nodes[edge_nodes[0]].curv -
        mesh->nodes[edge_nodes[0]].ref_curv);
      m[1] = E_B*(mesh->nodes[edge_nodes[1]].curv -
        mesh->nodes[edge_nodes[1]].ref_curv);
      q[j] = (m[1] - m[0])/l[j];
    }
    /** We then differentiate the transverse shear
    tension q=dm/dl at the considered node. */
    coord fb;
    foreach_dimension()
      fb.x = (q[1]*mesh->edges[mesh->nodes[i].edge_ids[1]].normal.x -
         q[0]*mesh->edges[mesh->nodes[i].edge_ids[0]].normal.x)/
         (.5*(l[0] + l[1]));
    foreach_dimension() mesh->nodes[i].lagForce.x += fb.x;
  }
}

event acceleration (i++) {
  for(int i=0; i<mbs.nbmb; i++) bending(&mbs.mb[i]);
}
