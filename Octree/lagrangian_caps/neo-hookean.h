/** In this file we compute the neo-Hookean force on each Lagrangian node of
the membrane(s).*/

#ifndef E_S
  #define E_S 1.
#endif

void neo_hookean(lagMesh* mesh) {
  for(int i=0; i<mesh->nlp; i++) {
    coord t[2];
    for(int j=0; j<2; j++) {
      int edge_id, edge_node1, edge_node2;
      edge_id = mesh->nodes[i].edge_ids[j];
      double stretch_cube = cube(mesh->edges[edge_id].st);
      double tension_norm = E_S*(stretch_cube - 1.)/sqrt(stretch_cube);
      /** We compute the direction vector $e$ for the tension */
      edge_node1 = mesh->edges[edge_id].vertex_ids[0];
      edge_node2 = mesh->edges[edge_id].vertex_ids[1];
      coord e;
      double ne = 0.;
      foreach_dimension() {
        e.x = mesh->nodes[edge_node1].pos.x - mesh->nodes[edge_node2].pos.x;
        ne += sq(e.x);
      }
      ne = sqrt(ne);
      /** $\bm{T_i} = \frac{E_s}{\lamba_i^{\frac{3}{2}}} (\lambda^3 - 1)
      \bm{e_i}$ */
      foreach_dimension()
        t[j].x = tension_norm*e.x/ne;
    }
    foreach_dimension() mesh->lagForces[i].x = t[0].x - t[1].x;
  }
}

event acceleration (i++) {
  for(int i=0; i<mbs.nbmb; i++) neo_hookean(&mbs.mb[i]);
}
