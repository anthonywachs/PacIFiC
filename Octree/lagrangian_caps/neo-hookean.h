/** In this file we compute the neo-Hookean force on each Lagrangian node of
the membrane(s).*/

#ifndef E_S
  #define E_S 1.
#endif

void neo_hookean(lagMesh* mesh) {
  for(int node=0; node<mesh->nlp; node++) {
    coord t[2];
    for(int edge=0; edge<2; edge++){
      // check if t[0] and t[1] need to be swapped
      // fill t[0], t[1]
      double stretch = mesh->vertices[node]->edges[edge]
      double tension = E_S*(
    }
    foreach_dimension() mesh->lagForces[node].x = t[0].x - t[1].x;
  }
}

event acceleration (i++) {
  for(int i=0; i<mbs.nbmb; i++) neo_hookean(&mbs.mb[i]);
}
