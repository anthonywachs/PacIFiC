//IBM interpolation stencil type
#ifndef IBM_stencil
#define IBM_stencil (1)
#endif

double stencil(double dist)
{

    double stencil = 0.;

    if (IBM_stencil == 1)
    {
        if (fabs(dist) <= 1.)
        {
            stencil = 1. - fabs(dist);
        }
        else
        {
            stencil = 0.;
        }   
    }
    else if (IBM_stencil == 4)
    {
        if (fabs(dist) < 1.)
        {
            stencil = 0.125 * (3. - 2 * fabs(dist) + sqrt(1. + 4 * fabs(dist) - 4 * sq(dist)));
        }
        else if (fabs(dist) <= 2. && fabs(dist) > 1.)
        {
            stencil = 0.125 * (5. - 2 * fabs(dist) - sqrt(-7. + 12 * fabs(dist) - 4 * sq(dist)));
        }
    }
    else if (IBM_stencil == 11)
    {
        if (fabs(dist) <= 0.5)
        {
            stencil = 3./4. - sq(dist);
        }
        else if (fabs(dist) > 0.5 && fabs(dist) <= 1.5)
        {
            stencil = 9./8. - 3.*fabs(dist)/2.+ sq(dist)/2.;
        }
        else
        {
            stencil = 0.;
        }
    }
    else if (IBM_stencil == 14)
    {
        if (fabs(dist) <= 0.5)
        {
            stencil = 3./8. + M_PI/32. - sq(dist)/4.;
        }
        else if (fabs(dist) > 0.5 && fabs(dist) <= 1.5)
        {
            stencil = 1./4. + (1-fabs(dist))*sqrt(-2.+8.*fabs(dist)-4*sq(dist))/8. - asin(sqrt(2.)*(fabs(dist)-1.))/8.;
        }
        else if (fabs(dist) > 1.5 && fabs(dist) <= 2.5)
        {
            stencil = 17./16.-M_PI/64.-3.*fabs(dist)/4.+sq(dist)/8.+(fabs(dist)-2.)*sqrt(-14.+16.*fabs(dist)-4.*sq(dist))/16.+asin(sqrt(2.)*(fabs(dist)-2.))/16.;
        }
        else
        {
            stencil = 0.;
        }
    }

    return stencil;
}
