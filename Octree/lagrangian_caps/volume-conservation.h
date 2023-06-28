/**
# Enforcing volume conservation

In this file, we follow the method of Sigüenza et al. to enforce exact volume 
conservation of the capsule. The details can be found in appendix A of [1](#siguenza2016validation).
*/

#include "smallest_root_cubic.h"

coord* correct_periodic_nodes_pos(coord* result, coord a, coord b, coord ref) {
    foreach_dimension() {
        result[0].x = (fabs(a.x - ref.x) < L0/2) ? a.x : 
            (a.x > ref.x) ? a.x - L0 : a.x + L0;
        result[1].x = (fabs(b.x - ref.x) < L0/2) ? b.x : 
            (b.x > ref.x) ? b.x - L0 : b.x + L0;
    }
}

foreach_dimension()
double periodic_friendly_cross_product_x(coord a, coord b, coord ref) {
    coord nodes[2];
    correct_periodic_nodes_pos(nodes, a, b, ref);
    return nodes[0].y*nodes[1].z - nodes[0].z*nodes[1].y;
}

double periodic_friendly_dot_product(coord a, coord b, coord ref) {
    coord nodes[2];
    correct_periodic_nodes_pos(nodes, a, b, ref);
    return cdot(nodes[0], nodes[1]);
}

trace
coord compute_alpha_m(lagMesh* mesh, int m) {
    coord alpha = {0., 0., 0.};

    for(int i=0; i<mesh->nodes[m].nb_triangles; i++) {
        int tid = mesh->nodes[m].triangle_ids[i];
        int local_id = 0;
        while (mesh->triangles[tid].node_ids[local_id] != m && local_id < 3) 
            local_id++;
        int nids[2];    // `nids` for "neighbor ids"
        nids[0] = mesh->triangles[tid].node_ids[(local_id + 1)%3];
        nids[1] = mesh->triangles[tid].node_ids[(local_id + 2)%3];

        foreach_dimension()
            alpha.x += -periodic_friendly_cross_product_x(
                mesh->nodes[nids[0]].pos, mesh->nodes[nids[1]].pos, 
                mesh->centroid)/12;
    }

    return alpha;
}

trace
void enforce_optimal_volume_conservation(lagMesh* mesh) {
    coord* alpha = malloc(mesh->nlp*sizeof(coord));
    for(int m=0; m<mesh->nlp; m++)
        alpha[m] = compute_alpha_m(mesh, m);

    double coeff_polynomial[4];
    for(int j=0; j<mesh->nlt; j++) {
        int tn[3]; // `tn` for "triangle nodes"
        for(int k=0; k<3; k++) 
            tn[k] = mesh->triangles[j].node_ids[k];
        for(int k=0; k<3; k++) {
            coord cp0; coord cp2; // `cp` for "cross-product"
            foreach_dimension() {
                cp0.x = alpha[tn[(k+1)%3]].y*alpha[tn[(k+2)%3]].z
                    - alpha[tn[(k+1)%3]].z*alpha[tn[(k+2)%3]].y;
                cp2.x = periodic_friendly_cross_product_x(
                    mesh->nodes[tn[(k+1)%3]].pos, mesh->nodes[tn[(k+2)%3]].pos,
                    mesh->centroid);
            }
            coeff_polynomial[3] += cdot(alpha[tn[k]], cp0);
            coeff_polynomial[2] += periodic_friendly_dot_product(
                mesh->nodes[tn[k]].pos, cp0, mesh->centroid);
            coeff_polynomial[1] += cdot(alpha[tn[k]], cp2);
        }
    }
    coeff_polynomial[3] /= 18;
    coeff_polynomial[2] /= 6;
    coeff_polynomial[1] /= 6;
    coeff_polynomial[0] = mesh->volume - mesh->initial_volume;
    double normalize_factor = 
        max(max(fabs(coeff_polynomial[3]), coeff_polynomial[2]), 
        max(fabs(coeff_polynomial[1]), coeff_polynomial[0]));
    coeff_polynomial[3] /= normalize_factor;
    coeff_polynomial[2] /= normalize_factor;
    coeff_polynomial[1] /= normalize_factor;
    coeff_polynomial[0] /= normalize_factor;
    double lambda = find_smallest_real_root(coeff_polynomial);

    for(int i=0; i<mesh->nlp; i++) {
        foreach_dimension() 
            mesh->nodes[i].pos.x += lambda*alpha[i].x;
    }
    free(alpha);
    correct_lag_pos(mesh);
}

/**
## References

~~~bib
@Article{siguenza2016validation,
  author    = {Sig{\"u}enza, Julien and Mendez, Simon and Ambard, Dominique and Dubois, Fr{\'e}d{\'e}ric and Jourdan, Franck and Mozul, R{\'e}my and Nicoud, Franck},
  journal   = {Journal of Computational Physics},
  title     = {Validation of an immersed thick boundary method for simulating fluid--structure interactions of deformable membranes},
  year      = {2016},
  pages     = {723--746},
  volume    = {322},
  file      = {:1-s2.0-S0021999116302662-main.pdf:PDF},
  publisher = {Elsevier},
}

~~~
*/