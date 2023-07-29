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
      if (mbs.mb[k].isactive) {
        comp_curvature(&mbs.mb[k]);
        for(int i=0; i<mbs.mb[k].nln; i++) {
          #if REF_CURV
            mbs.mb[k].nodes[i].ref_curv = mbs.mb[k].nodes[i].curv;
          #else
            mbs.mb[k].nodes[i].ref_curv = 0.;
          #endif
        }
      }
    }
  }

void bending(lagMesh* mesh) {
  comp_curvature(mesh);
  coord* raw_bending = malloc(mesh->nln*sizeof(coord));
  for(int i=0; i<mesh->nln; i++) {
    double l[2];
    for(int j=0; j<2; j++)
      l[j] = mesh->edges[mesh->nodes[i].edge_ids[j]].length;
    double lavg = .5*(l[0] + l[1]);
    double curv[3];
    for(int j=0; j<3; j++) {
      int neighbor_index = (i-1+j) > 0 ? (i-1+j)%mesh->nln : mesh->nln - 1;
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
  for(int i=0; i<mesh->nln; i++) {
    // double coeff[5] = {-1./16., 1./4., 5./8., 1./4., -1./16.};
    double coeff[7] = {-1./32., 0., 9./32., 1./2., 9./32., 0., -1./32.};
    for(int j=0; j<7; j++)
      foreach_dimension()
        mesh->nodes[i].lagForce.x +=
          coeff[j]*raw_bending[(mesh->nln + i - 3 + j)%mesh->nln].x;
  }
  free(raw_bending);
}

event acceleration (i++) {
  for(int i=0; i<mbs.nbmb; i++) if (mbs.mb[i].isactive) bending(&mbs.mb[i]);
}
