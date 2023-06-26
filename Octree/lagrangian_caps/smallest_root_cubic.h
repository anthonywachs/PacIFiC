typedef struct realQuadraticRoots {
    bool real;
    double roots[2];
} realQuadraticRoots;

realQuadraticRoots real_quadratic_roots(double* p) {
    realQuadraticRoots result;
    double a, b, c;
    a = p[2]; b = p[1]; c = p[0];
    double delta = sq(b) - 4*a*c;
    if (fabs(delta) > 1.e-10 && delta < 0) {
        result.real = false;
        result.roots[0] = HUGE;
        result.roots[1] = HUGE;
    }
    if (fabs(delta) > 1.e-10 && delta > 0) {
        result.real = true;
        result.roots[0] = (-b - sqrt(delta))/(2*a);
        result.roots[1] = (-b + sqrt(delta))/(2*a);
    }
    if (fabs(delta) < 1.e-10) {
        result.real = true;
        result.roots[0] = -b/(2*a);
        result.roots[1] = -b/(2*a);
    }
    return result;
}

double find_smallest_real_root(double* a) {
    double tolerance = 1e-10;
    double epsilon = 1.e-10;
        
    if (fabs(a[3]) > epsilon) {
        /** Solve 3rd-order polynomial:

        We use Berstow's method to factor out a 2nd-order polynomial. We then 
        compute its roots as well as the root of the linear function, and 
        return the real root of smallest modulus. */
        int i = 0;
        int max_iterations = 100;
        
        /** Step 1: factorization */
        double ui, vi, uj, vj;
        ui = a[2]/a[3];
        vi = a[1]/a[3];
        double c = a[3]*sq(ui) - a[2]*ui - a[3]*vi + a[1];
        double d = a[3]*ui*vi -a[2]*vi + a[0];
        printf("c = %g\n", fabs(c));
        printf("d = %g\n", fabs(d));

        while ((fabs(c) > tolerance || fabs(d) > tolerance)
            && i < max_iterations) {
            printf("hi! c = %g and d = %g\n", fabs(c), fabs(d));
            double det = (a[3]*ui - a[2])*(2*a[3]*ui - a[2]) + sq(a[3])*vi;
            if (fabs(det) < epsilon) {
                fprintf(stderr, "Error: zero determinant in Berstow's method.");
                return HUGE;
            }
            uj = ui - ((a[3]*ui - a[2])*c + a[3]*d)/det;
            vj = vi - (-a[3]*vi*c + (2*a[3]*ui - a[2])*d)/det;
            ui = uj;
            vi = vj;
            c = a[3]*sq(ui) - a[2]*ui - a[3]*vi + a[1];
            d = a[3]*ui*vi -a[2]*vi + a[0];
            i++;
        }

        double root0 = ui - a[2]/a[3];
        printf("root0 = %g\n", root0);
        double quadratic_coeff[3];
        quadratic_coeff[0] = a[2] - ui*a[3];
        quadratic_coeff[1] = a[3];
        quadratic_coeff[2] = 1;
        realQuadraticRoots my_roots = real_quadratic_roots(quadratic_coeff);
        double root1 = my_roots.roots[0];
        double smallest_real_root = (fabs(root1) > fabs(root0)) ? root0 : root1;
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
        return my_roots.roots[0];
    }

    if (fabs(a[1]) > epsilon) {
        /** Solve linear function */
        return -a[0]/a[1];
    }

    fprintf(stderr, "Warning: giving constant function instead of third-order \
        polynomial in find_smallest_real_root.");
    return 0.;
}