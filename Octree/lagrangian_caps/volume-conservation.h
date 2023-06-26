/**
# Enforcing volume conservation

In this file, we follow the method of Sigüenza et al. to enforce exact volume 
conservation of the capsule. The details can be found in appendix A of [1](#siguenza2016validation).
*/

double third_order_polynomial(double[4] c, double x) {
    return (c[0]*cube(x) + c[1]*sq(x) + c[2]*x + c[3]);
}

double derivative_third_order_polynomial(double[4] c, double x) {
    return (3*c[0]*sq(x) + 2*c[1]*x + c[2]);
}

double find_smallest_real_root(double[4] c) {
    /** Solve 3rd-order polynomial:
    We use Newton's method with initial guess 0 */
    double tolerance = 1e-5;
    int i = 0;
    int max_iterations 25;
    double xn = 0;
    double fxn = third_order_polynomial(c, xn);
    while ((fabs(fxn) > tolerance) && (i < max_iterations)) {
        double dfxn = derivative_third_order_polynomial(c, xn);
        if (dfxn < 1.e-10) {
            fprintf(stderr, "Error: zero derivative in Newton's method.");
            return xn;
        }
        xn = xn - fxn/dfxn; 
        fxn = third_order_polynomial(xn);
        i++;
    }
    return xn;
}

void compute_alpha_m(lagMesh* mesh, int m) {
    double alpha = 0;

    

    return -alpha/12;
}

void enforce_optimal_volume_conservation(lagMesh* mesh) {
    coord* alpha = malloc(mesh->nlp*sizeof(coord));
    for(int m=0; m<mesh->nlp; m++)
        alpha[m] = compute_alpha_m(mesh, m);

    for(int i=0; i<mesh->nlp; i++) {
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
                    cp2.x = 
                        mesh->nodes[tn[(k+1)%3]].y*mesh->nodes[tn[(k+2)%3]].z
                        - mesh->nodes[tn[(k+1)%3]].z*mesh->nodes[tn[(k+1)%3]].y;
                }
                coeff_polynomial[0] += cdot(alpha[tn[k]], cp0);
                coeff_polynomial[1] += cdot(mesh->nodes[tn[k]], cp0);
                coeff_polynomial[2] += cdot(alpha[tn[k]], cp2);
            }
        }
        coeff_polynomial[0] /= 18;
        coeff_polynomial[1] /= 6;
        coeff_polynomial[2] /= 6;
        coeff_polynomial[4] = mesh->volume - mesh->initial_volume;

        double lambda = find_smallest_real_root(coeff_polynomial);

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