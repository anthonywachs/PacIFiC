/** In this file we compute the neo-Hookean force on each Lagrangian node of
the membrane(s).*/

#ifndef E_S
  #define E_S 1.
#endif

void neo_hookean(lagMesh* mesh) {
  compute_lengths(mesh);
  for(int i=0; i<mesh->nlp; i++) {
    coord T[2];
    for(int j=0; j<2; j++) {
      int edge_id, edge_node1, edge_node2;
      edge_id = mesh->nodes[i].edge_ids[j];
      double stretch_cube =
        cube(mesh->edges[edge_id].length/mesh->edges[edge_id].l0);
      double tension_norm = (fabs(stretch_cube) > 1.e-10) ?
        E_S*(stretch_cube - 1.)/sqrt(stretch_cube) : 0.;
      /** We compute the direction vector $e$ for the tension */
      edge_node1 = mesh->edges[edge_id].node_ids[0];
      edge_node2 = mesh->edges[edge_id].node_ids[1];
      coord e;
      double ne = 0.;
      foreach_dimension() {
        double x1 = mesh->nodes[edge_node1].pos.x;
        double x2 = mesh->nodes[edge_node2].pos.x;
        e.x = (fabs(x1 - x2) < L0/2.) ? x1 - x2 : ((fabs(x1 - L0 - x2) > L0/2.)
          ? x1 + L0 - x2 : x1 - L0 - x2) ;
        ne += sq(e.x);
      }
      ne = sqrt(ne);
      /** $\bm{T_i} = \frac{E_s}{\lamba_i^{\frac{3}{2}}} (\lambda^3 - 1)
      \bm{e_i}$ */
      foreach_dimension()
        T[j].x = (fabs(ne) > 1.e-10) ? tension_norm*e.x/ne : 0.;
    }
    foreach_dimension() mesh->nodes[i].lagForce.x += T[0].x - T[1].x;
  }
}

event acceleration (i++) {
  for(int i=0; i<mbs.nbmb; i++) neo_hookean(&mbs.mb[i]);
}
