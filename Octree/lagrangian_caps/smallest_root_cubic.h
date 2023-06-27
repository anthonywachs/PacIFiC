typedef struct realQuadraticRoots {
    bool real;
    double roots[2];
} realQuadraticRoots;

realQuadraticRoots real_quadratic_roots(double* p) {
    realQuadraticRoots result;
    double a, b, c;
    a = p[2]; b = p[1]; c = p[0];
    if (fabs(a) < 1.e-10) {
        result.real = true;
        result.roots[0] = (fabs(b) > 1.e-10) ? -c/b : HUGE;
        result.roots[1] = HUGE;
        return result;
    }

    double delta = sq(b) - 4*a*c;
    if (fabs(delta) > 1.e-10 && delta < 0) {
        result.real = false;
        result.roots[0] = HUGE;
        result.roots[1] = HUGE;
    }
    if (fabs(delta) > 1.e-10 && delta > 0) {
        result.real = true;
        double r1, r2;
        r1 = (-b - sqrt(delta))/(2*a);
        r2 = (-b + sqrt(delta))/(2*a);
        result.roots[0] = fabs(r1) > fabs(r2) ? r2 : r1;
        result.roots[1] = fabs(r1) > fabs(r2) ? r1 : r2;
    }
    if (fabs(delta) < 1.e-10) {
        result.real = true;
        result.roots[0] = -b/(2*a);
        result.roots[1] = -b/(2*a);
    }
    return result;
}

double compute_b0(double* a, double u, double v) {
    return (a[1] - v*a[3] - u*a[2] + sq(u)*a[3]);
}

double compute_b1(double* a, double u, double v) {
    return (a[2] - u*a[3]);
}

double compute_b2(double* a, double u, double v) {
    return (a[3]);
}

double compute_c(double* a, double u, double v) {
    // return (a[0] - v*(a[2] - u*a[3]) - u*(a[1] - v*a[3] - u*a[2] + sq(u)*a[3]));
    double b1 = compute_b1(a, u, v);
    double b0 = compute_b0(a, u, v);
    return (a[0] - v*b1 - u*b0);
}

double compute_d(double* a, double u, double v) {
    // return (-v*(a[1] - v*a[3] - u*a[2] + sq(u)*a[3]));
    double b0 = compute_b0(a, u, v);
    return (-v*b0); // for us a_0 = 0
}

trace
double find_smallest_real_root(double* a) {
    double tolerance = 1e-10;
    double epsilon = 1.e-6;
        
    if (fabs(a[3]) > epsilon) {
        /** Solve 3rd-order polynomial:

        We use Berstow's method to factor out a 2nd-order polynomial. We then 
        compute its roots as well as the root of the linear function, and 
        return the real root of smallest modulus. */
        int i = 0;
        int max_iterations = 200;
        
        /** Step 1: factorization */
        double ui, vi, uj, vj;
        ui = a[2]/a[3];
        vi = a[1]/a[3];
        double c = compute_c(a, ui, vi);
        double d = compute_d(a, ui, vi);

        double b[3];
        while ((fabs(c) > tolerance || fabs(d) > tolerance)
            && i < max_iterations) {
            double g, h;
            b[0] = compute_b0(a, ui, vi);
            b[1] = compute_b1(a, ui, vi);
            b[2] = a[3];
            h = b[0] - vi*b[2];
            g = b[1] - ui*b[2];
            double det = vi*sq(g) + h*(h - ui*g);
            if (fabs(det) < epsilon) {
                fprintf(stderr, "Error: zero determinant in Berstow's method.");
                return HUGE;
            }
            uj = ui - (-h*c + g*d)/det;
            vj = vi - (-g*vi*c + (g*ui - h)*d)/det;
            ui = uj;
            vi = vj;
            c = compute_c(a, ui, vi);
            d = compute_d(a, ui, vi);
            i++;
        }
        fprintf(stderr, "Number of iterations in Bairstow's method: %d\n", i);
        b[0] = compute_b0(a, ui, vi);
        b[1] = compute_b1(a, ui, vi);
        b[2] = a[3];

        double quadratic_coeff[3];
        quadratic_coeff[0] = vi;
        quadratic_coeff[1] = ui;
        quadratic_coeff[2] = 1.;
        realQuadraticRoots my_roots = real_quadratic_roots(quadratic_coeff);
        // printf("vi=%g\n",fabs(vi));
        double root1 = fabs(vi) > epsilon ? my_roots.roots[0] : 
            my_roots.roots[1];
        // printf("root1=%g, root2=%g, ", my_roots.roots[0], my_roots.roots[1]);
        quadratic_coeff[0] = b[0];
        quadratic_coeff[1] = b[1];
        quadratic_coeff[2] = b[2];
        my_roots = real_quadratic_roots(quadratic_coeff);
        // printf("root3=%g, root4=%g\n", my_roots.roots[0], my_roots.roots[1]);
        double root2 = fabs(vi) > epsilon ? my_roots.roots[1] : 
            my_roots.roots[0];
        double smallest_real_root = (fabs(root1) > fabs(root2)) ? root2 : root1;
        return smallest_real_root;
    }

    if (fabs(a[2]) > epsilon) {
        /** Solve 2nd-order polynomial */
        realQuadraticRoots my_roots = real_quadratic_roots(a);
        if (my_roots.real == false) {
            fprintf(stderr, "Warning: no real roots found in \
                find_smallest_real_root");
            return 0.;
        }
        double root1, root2;
        root1 = my_roots.roots[0];
        root2 = my_roots.roots[1];
        double smallest_real_root = (fabs(root1) > fabs(root2)) ? root2 : root1;
        return smallest_real_root;
    }

    if (fabs(a[1]) > epsilon) {
        /** Solve linear function */
        return -a[0]/a[1];
    }

    fprintf(stderr, "Warning: giving constant function instead of third-order \
        polynomial in find_smallest_real_root.");
    return 0.;
}