/**
In this file, we use the Lagrangian mesh to compute the stretches, and the
stresses associated to a specific elastic law. Default is the Neo-Hookean model.
*/

#ifndef DWDL1
  #ifndef E_S
    #define E_S 1.
  #endif
  #define DWDL1(L1, L2) (E_S*L1*(1. - 1./(sq(L1*L2)))/3.)
  #define DWDL2(L1, L2) (E_S*L2*(1. - 1./(sq(L1*L2)))/3.)
#endif

#if dimension < 3
  #define cnorm(a) (sqrt(sq(a.x) + sq(a.y)))
  #define cdot(a,b) (a.x*b.x + a.y*b.y)
#else
  #define cnorm(a) (sqrt(sq(a.x) + sq(a.y) + sq(a.z)))
  #define cdot(a,b) (a.x*b.x + a.y*b.y + a.z*b.z)
#endif

#if dimension > 2
/** The function below returns the vertices of the triangle $tid$, rotated to
the reference plane. Since node 0 is always located at $(0,0,0)$, this
function only returns the coordinates of nodes 1 and 2.
*/
void rotate_to_reference_plane(lagMesh* mesh, int tid, coord* rn) {
  if (!mesh->updated_normals) comp_triangle_area_normals(mesh);
  int nodes[3];
  for(int i=0; i<3; i++) nodes[i] = mesh->triangles[tid].node_ids[i];

  /** Step 1. compute the rotation matrix */
  coord er[3], ec[3]; // reference and current coordinate systems
  er[0].x = 1.; er[0].y = 0.; er[0].z = 0.;
  er[1].x = 0.; er[1].y = 1.; er[1].z = 0.;
  er[2].x = 0.; er[2].y = 0.; er[2].z = 1.;
  /** Following Doddi and Bagchi (2008), $ec[0]$ is in the direction of
  [node0, node2]; $ec[2]$ is the vector normal to the triangle and
  $ec[1] = ec[2] x ec[0]$*/
  foreach_dimension() {
    ec[0].x = mesh->nodes[nodes[2]].pos.x - mesh->nodes[nodes[0]].pos.x;
    ec[2].x = mesh->triangles[tid].normal.x;
  }
  double enorm = cnorm(ec[0]);
  foreach_dimension() ec[0].x /= enorm;
  foreach_dimension() ec[1].x = ec[2].y*ec[0].z - ec[2].z*ec[0].y;
  enorm = cnorm(ec[1]);
  foreach_dimension() ec[1].x /= enorm;
  double M[9];
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      M[3*i+j] = cdot(ec[i], er[j]);

  /** Step 2. rotate the triangle to the reference basis */
  double refNode[3];
  for(int k=0; k<2; k++) {
    double cv[3]; // cv for "current vector"
    cv[0] = mesh->nodes[nodes[k+1]].pos.x - mesh->nodes[nodes[0]].pos.x;
    cv[1] = mesh->nodes[nodes[k+1]].pos.y - mesh->nodes[nodes[0]].pos.y;
    cv[2] = mesh->nodes[nodes[k+1]].pos.z - mesh->nodes[nodes[0]].pos.z;
    for(int i=0; i<3; i++) {
      refNode[i] = 0.;
      for(int j=0; j<3; j++)
        refNode[i] += M[3*i+j]*cv[j];
    }
    rn[k].x = refNode[0];
    rn[k].y = refNode[1];
    rn[k].z = refNode[2];
  }
}

void store_initial_configuration(lagMesh* mesh) {
  for(int i=0; i<mesh->nlt; i++)
    rotate_to_reference_plane(mesh, i, mesh->triangles[i].refShape);
}

event init (i = 0) {
  for(int j=0; j<NCAPS; j++) store_initial_configuration(&MB(j));
}
#endif

void comp_elastic_stress(lagMesh* mesh) {
  #if dimension < 3
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
#else
  /** Loop through each triangle */
  for(int i=0; i<mesh->nlt; i++) {
    /** 1. Rotate triangle to common plane using the rotation matrix $\bm{R}$ */
    coord cn[2];
    rotate_to_reference_plane(mesh, i, cn);

    /** 2. Compute the shape functions of the triangle */

    /** 3. Compute the displacement $v_k$ of each node */

    /** 4. Compute the right Cauchy-Green deformation tensor from $v_k$:
    $$ \bm{C} = \bm{F^t}\bm{F}\;, \quad \bm{F} =  frac{\partial v_k}{\partial
    \bm{x^P}} $$
    with $\bm{x^P} = [x_p, y_p]$ the two-dimensional coordinates of the common-
    plane
    */

    /** 5. Compute the two principal stretches $\epsilon_1$, $\epsilon_2$ from
    \bm{C} */

    // /** 6. For each node of the triangle, compute the force in the common plane,
    // then rotate it and add it to the Lagrangian force of the node */
    // for(int j=0; j<3; j++) {
    //   /** 6.1 Compute $\frac{\partial \epsilon_1}{\partial \bm{v_j}}$ and
    //   $\frac{\partial \epsilon_2}{\partial \bm{v_j}}$ */
    //
    //   /** 6.2 Compute $f_j^P = \frac{\partial W}{\partial \epsilon_1}
    //   \frac{\partial \epsilon_1}{\partial \bm{v_j}} +
    //   \frac{\partial W}{\partial \epsilon_2}
    //   \frac{\partial \epsilon_2}{\partial \bm{v_j}} $*/
    //
    //   /** 6.3 Rotate the force in the common plane to the current plane:
    //   $f_j = \bm{R^T} f_j^P$ */
    // }
  }
#endif
}
