
#include "GJK_SV.hh"


double rel_error1 = 1.e-8;   // relative error in the computed distance
double abs_error1 = 1.e-12;  // absolute error if the distance is almost zero


/* ========================================================================== */
/*                             Low-Level Methods                              */
/* ========================================================================== */
// TODO: COMMENT ???
static inline double determinant( Vector3 const p,
                             Vector3 const q,
                             Vector3 const r )
{
  return ( p[0] * ( (q[1] * r[2] ) - ( r[1] * q[2] ) ) - 
           p[1] * ( q[0] * r[2] - r[0] * q[2] ) +
           p[2] * ( q[0] * r[1] - r[0] * q[1] ) );
}




// -----------------------------------------------------------------------------
// TODO: COMMENT ???
static inline void projectOnLine( Vector3 const p,
                                  Vector3 const q,
                                  Vector3& v )
{
    Vector3 pq = p - q;
    double const tmp = ( p * pq ) / ( pq * pq);
    v = p - tmp * pq;
}




// -----------------------------------------------------------------------------
// TODO: COMMENT ???
static inline void projectOnPlane( Vector3 const p,
                                   Vector3 const q,
                                   Vector3 const r, 
                                   Vector3& v )
{
    Vector3 n = ( Vector3( p - q ) ) ^ ( Vector3( p - r ) );
    double const tmp = ( n * p ) / ( n * n );
    v = tmp * n;
}




// -----------------------------------------------------------------------------
// TODO: COMMENT ???
static inline int hff1( Vector3 const p, 
                        Vector3 const q )
{
    if ( Norm2( p ) - p * q > 0 )
        return 1; // keep q
    return 0;
}




// -----------------------------------------------------------------------------
// TODO: COMMENT ???
static inline int hff2( Vector3 const p,
                        Vector3 const q,
                        Vector3 const r )
{
    Vector3 pq = q - p;
    Vector3 nTemp = pq ^ ( r - p );
    Vector3 n = pq ^ nTemp;
    return ( p * n < 0 ); // Discard r if true
}




// -----------------------------------------------------------------------------
// TODO: COMMENT ???
static inline int hff3( Vector3 const p,
                        Vector3 const q,
                        Vector3 const r )
{
    Vector3 n = ( Vector3( q - p ) ) ^ ( Vector3( r - p ) );
    return ( p * n <= 0 ); // discard s if true
}




// -----------------------------------------------------------------------------
// Handling the case where the simplex is of kind 1-simplex ( line )
static inline void S1D( Simplex& s, 
                        Vector3& v )
{
    Vector3 const s1p = s.vrtx[1];
    Vector3 const s2p = s.vrtx[0];

    if ( hff1( s1p, s2p ) ) 
    {
        // update v, no need to update s, return V{1,2}
        projectOnLine( s1p, s2p, v );
        return;
    } 
    else 
    {
        // Update v and s, return V{1}
        v = s.vrtx[1];
        s.nvrtx = 1;
        s.vrtx[0] = s.vrtx[1];
        return;
    }
}




// -----------------------------------------------------------------------------
// Handling the case where the simplex is of kind 2-simplex ( triangle )
static inline void S2D( Simplex& s,
                        Vector3& v ) 
{
    Vector3 const s1p = s.vrtx[2];
    Vector3 const s2p = s.vrtx[1];
    Vector3 const s3p = s.vrtx[0];
    int const hff1f_s12 = hff1( s1p, s2p );
    int const hff1f_s13 = hff1( s1p, s3p );

    if ( hff1f_s12 ) 
    {
        int const hff2f_23 = !hff2( s1p, s2p, s3p );
        if ( hff2f_23 )
        {
            if ( hff1f_s13 )
            {
                int const hff2f_32 = !hff2( s1p, s3p, s2p );
                if ( hff2f_32 )
                {
                    projectOnPlane( s1p, s2p, s3p, v ); // update s, 
                                                        // no need to update c
                    return;                             // return V{1,2,3}
                } 
                else
                {
                    // update v, update s, return V{1,3}
                    projectOnLine( s1p, s3p, v );
                    s.nvrtx = 2;
                    s.vrtx[1] = s.vrtx[2];
                    return;
                }
            } 
            else
            {
                projectOnPlane( s1p, s2p, s3p, v ); // update s, 
                                                    // no need to update c
                return;                             // return V{1,2,3}
            }
        } 
        else
        {
            // update v, update s, return V{1,2}
            projectOnLine( s1p, s2p, v );
            s.nvrtx = 2;
            s.vrtx[0] = s.vrtx[2];
            return;
        }
    } 
    else if ( hff1f_s13 ) 
    {
        int const hff2f_32 = !hff2( s1p, s3p, s2p );
        if ( hff2f_32 ) 
        {
            projectOnPlane( s1p, s2p, s3p, v ); // update s, no need to update v
            return;                             // return V{1,2,3}
        }
        else
        {
            // update v, update s, return V{1,3}
            projectOnLine( s1p, s3p, v );
            s.nvrtx = 2;
            s.vrtx[1] = s.vrtx[2];
            return;
        }
    }
    else
    {
        // update s and v, return V{1}
        v = s.vrtx[2];
        s.nvrtx = 1;
        s.vrtx[0] = s.vrtx[2];
        return;
    }
}





// -----------------------------------------------------------------------------
// Handling the case where the simplex is of kind 3-simplex ( tetrahedron )
static inline void S3D( Simplex& s, 
                        Vector3& v )
{
    Vector3 s1, s2, s3, s4, s1s2, s1s3, s1s4;
    Vector3 si, sj, sk;
    int testLineThree, testLineFour, testPlaneTwo, 
        testPlaneThree, testPlaneFour, dotTotal;
    int i, j, k;

    s1 = s.vrtx[3];
    s2 = s.vrtx[2];
    s3 = s.vrtx[1];
    s4 = s.vrtx[0];
    s1s2 = s2 - s1;
    s1s3 = s3 - s1;
    s1s4 = s4 - s1;

    int hff1_tests[3];
    hff1_tests[2] = hff1( s1, s2 );
    hff1_tests[1] = hff1( s1, s3 );
    hff1_tests[0] = hff1( s1, s4 );
    testLineThree = hff1_tests[1];
    testLineFour  = hff1_tests[0];

    dotTotal = hff1_tests[2] + testLineThree + testLineFour;
    if ( dotTotal == 0 )
    {
        v = s1;
        s.nvrtx = 1;
        s.vrtx[0] = s1;
        return;
    }

    double const det134 = determinant( s1s3, s1s4, s1s2 );
    int const sss = ( det134 <= 0 );

    testPlaneTwo = hff3( s1, s3, s4 ) - sss;
    testPlaneTwo = testPlaneTwo * testPlaneTwo;
    testPlaneThree = hff3( s1, s4, s2 ) - sss;
    testPlaneThree = testPlaneThree * testPlaneThree;
    testPlaneFour = hff3( s1, s2, s3 ) - sss;
    testPlaneFour = testPlaneFour * testPlaneFour;

    switch ( testPlaneTwo + testPlaneThree + testPlaneFour )
    {
        case 3:
            v = Vector3( 0., 0., 0. );
            s.nvrtx = 4;
            break;
        case 2:
            // Only one facing the origin
            // 1,i,j, are the indices of the points on the triangle and remove k
            // from simplex
            s.nvrtx = 3;
            if ( !testPlaneTwo ) // removes s2
                s.vrtx[2] = s.vrtx[3];
            else if ( !testPlaneThree ) // removes s3
            { 
                s.vrtx[1] = s2;
                s.vrtx[2] = s.vrtx[3];
            }
            else if ( !testPlaneFour ) // removes s4 - no need to reorder
            {
                s.vrtx[0] = s3;
                s.vrtx[1] = s2;
                s.vrtx[2] = s.vrtx[3];
            }
            // Call S2D
            S2D( s, v );
            break;
        case 1:
            // Two triangles face the origins:
            // The only positive hff3 is for triangle 1,i,j, therefore k must be
            // in the solution as it supports the the point of minimum norm.
            // 1,i,j, are the indices of the points on the triangle and remove k
            // from simplex
            s.nvrtx = 3;
            if ( testPlaneTwo )
            {
                k = 2; // s2
                i = 1;
                j = 0;
            }
            else if ( testPlaneThree )
            {
                k = 1; // s3
                i = 0;
                j = 2;
            }
            else
            {
                k = 0; // s4
                i = 2;
                j = 1;
            }

            si = s.vrtx[i];
            sj = s.vrtx[j];
            sk = s.vrtx[k];

            if ( dotTotal == 1 )
            {
                if ( hff1_tests[k] )
                {
                    if ( !hff2( s1, sk, si ) )
                    {
                        s.nvrtx = 3;
                        s.vrtx[2] = s.vrtx[3];
                        s.vrtx[1] = si;
                        s.vrtx[0] = sk;
                        projectOnPlane( s1, si, sk, v );
                    }
                    else if ( !hff2( s1, sk, sj ) )
                    {
                        s.nvrtx = 3;
                        s.vrtx[2] = s.vrtx[3];
                        s.vrtx[1] = sj;
                        s.vrtx[0] = sk;
                        projectOnPlane( s1, sj, sk, v );
                    } 
                    else 
                    {
                        s.nvrtx = 2;
                        s.vrtx[1] = s.vrtx[3];
                        s.vrtx[0] = sk;
                        projectOnLine(s1, sk, v);
                    }
                } 
                else if ( hff1_tests[i] )
                {
                    if ( !hff2( s1, si, sk ) )
                    {
                        s.nvrtx = 3;
                        s.vrtx[2] = s.vrtx[3];
                        s.vrtx[1] = si;
                        s.vrtx[0] = sk;
                        projectOnPlane( s1, si, sk, v );
                    } 
                    else
                    {
                        s.nvrtx = 2;
                        s.vrtx[1] = s.vrtx[3];
                        s.vrtx[0] = si;
                        projectOnLine( s1, si, v );
                    }
                } 
                else
                {
                    if ( !hff2( s1, sj, sk ) )
                    {
                        s.nvrtx = 3;
                        s.vrtx[2] = s.vrtx[3];
                        s.vrtx[1] = sj;
                        s.vrtx[0] = sk;
                        projectOnPlane( s1, sj, sk, v );
                    }
                    else
                    {
                        s.nvrtx = 2;
                        s.vrtx[1] = s.vrtx[3];
                        s.vrtx[0] = sj;
                        projectOnLine( s1, sj, v );
                    }
                }
            }
            else if ( dotTotal == 2 )
            {
                // Two edges have positive hff1, meaning that for two edges the
                // origin's project fall on the segement.
                // Certainly the edge 1,k supports the the point of minimum norm,
                // and so hff1_1k is positive.
                if ( hff1_tests[i] )
                {
                    if ( !hff2( s1, sk, si ) )
                    {
                        if ( !hff2( s1, si, sk ) ) 
                        {
                            s.nvrtx = 3;
                            s.vrtx[2] = s.vrtx[3];
                            s.vrtx[1] = si;
                            s.vrtx[0] = sk;
                            projectOnPlane( s1, si, sk, v );
                        }
                        else
                        {
                            s.nvrtx = 2;
                            s.vrtx[1] = s.vrtx[3];
                            s.vrtx[0] = sk;
                            projectOnLine( s1, sk, v );
                        }
                    }
                    else
                    {
                        if ( !hff2( s1, sk, sj ) )
                        {
                            s.nvrtx = 3;
                            s.vrtx[2] = s.vrtx[3];
                            s.vrtx[1] = sj;
                            s.vrtx[0] = sk;
                            projectOnPlane( s1, sj, sk, v );
                        } 
                        else
                        {
                            s.nvrtx = 2;
                            s.vrtx[1] = s.vrtx[3];
                            s.vrtx[0] = sk;
                            projectOnLine( s1, sk, v );
                        }
                    }
                }
                else if ( hff1_tests[j] ) // there is no other choice
                { 
                    if ( !hff2( s1, sk, sj ) )
                    {
                        if (!hff2(s1, sj, sk))
                        {
                            s.nvrtx = 3;
                            s.vrtx[2] = s.vrtx[3];
                            s.vrtx[1] = sj;
                            s.vrtx[0] = sk;
                            projectOnPlane( s1, sj, sk, v );
                        }
                        else
                        {
                            s.nvrtx = 2;
                            s.vrtx[1] = s.vrtx[3];
                            s.vrtx[0] = sj;
                            projectOnLine( s1, sj, v );
                        }
                    } 
                    else
                    {
                        if (!hff2(s1, sk, si)) 
                        {
                            s.nvrtx = 3;
                            s.vrtx[2] = s.vrtx[3];
                            s.vrtx[1] = si;
                            s.vrtx[0] = sk;
                            projectOnPlane( s1, si, sk, v );
                        }
                        else
                        {
                            s.nvrtx = 2;
                            s.vrtx[1] = s.vrtx[3];
                            s.vrtx[0] = sk;
                            projectOnLine( s1, sk, v );
                        }
                    }
                } 
                else
                {
                    // ERROR;
                }

            } 
            else if ( dotTotal == 3 )
            {
                // MM : ALL THIS HYPHOTESIS IS FALSE
                // sk is s.t. hff3 for sk < 0. So, sk must support the origin 
                // because there are 2 triangles facing the origin.
                int hff2_ik = hff2( s1, si, sk );
                int hff2_jk = hff2( s1, sj, sk );
                int hff2_ki = hff2( s1, sk, si );
                int hff2_kj = hff2( s1, sk, sj );

                // if ( hff2_ki == 0 && hff2_kj == 0 )
                // {
                //     mexPrintf("\n\n UNEXPECTED VALUES!!! \n\n");
                // }
                if ( hff2_ki == 1 && hff2_kj == 1 )
                {
                    s.nvrtx = 2;
                    s.vrtx[1] = s.vrtx[3];
                    s.vrtx[0] = sk;
                    projectOnLine( s1, sk, v );
                } 
                else if ( hff2_ki ) // discard i
                {
                    if ( hff2_jk ) // discard k
                    {       
                        s.nvrtx = 2;
                        s.vrtx[1] = s.vrtx[3];
                        s.vrtx[0] = sj;
                        projectOnLine( s1, sj, v );
                    }
                    else
                    {
                        s.nvrtx = 3;
                        s.vrtx[2] = s.vrtx[3];
                        s.vrtx[1] = sj;
                        s.vrtx[0] = sk;
                        projectOnPlane( s1, sk, sj, v );
                    }
                } 
                else // discard j
                {
                    if ( hff2_ik ) // discard k
                    {
                        s.nvrtx = 2;
                        s.vrtx[1] = s.vrtx[3];
                        s.vrtx[0] = si;
                        projectOnLine( s1, si, v );
                    }
                    else
                    {
                        s.nvrtx = 3;
                        s.vrtx[2] = s.vrtx[3];
                        s.vrtx[1] = si;
                        s.vrtx[0] = sk;
                        projectOnPlane( s1, sk, si, v );
                    }
                }
            }
            break;
        case 0:
            // The origin is outside all 3 triangles
            if ( dotTotal == 1 )
            {
                if (testLineThree) // Here si is set such that hff(s1,si) > 0
                {
                    k = 2;
                    i = 1; // s3
                    j = 0;
                }
                else if ( testLineFour ) 
                {
                    k = 1; // s3
                    i = 0;
                    j = 2;
                }
                else
                {
                    k = 0;
                    i = 2; // s2
                    j = 1;
                }
                si = s.vrtx[i];
                sj = s.vrtx[j];
                sk = s.vrtx[k];

                if ( !hff2( s1, si, sj ) )
                {
                    s.nvrtx = 3;
                    s.vrtx[2] = s.vrtx[3];
                    s.vrtx[1] = si;
                    s.vrtx[0] = sj;
                    projectOnPlane( s1, si, sj, v );
                }
                else if ( !hff2( s1, si, sk ) )
                {
                    s.nvrtx = 3;
                    s.vrtx[2] = s.vrtx[3];
                    s.vrtx[1] = si;
                    s.vrtx[0] = sk;
                    projectOnPlane( s1, si, sk, v );
                }
                else
                {
                    s.nvrtx = 2;
                    s.vrtx[1] = s.vrtx[3];
                    s.vrtx[0] = si;
                    projectOnLine( s1, si, v );
                }
            }
            else if ( dotTotal == 2 )
            {
                // Here si is set such that hff(s1,si) < 0
                s.nvrtx = 3;
                if ( !testLineThree )
                {
                    k = 2;
                    i = 1; // s3
                    j = 0;
                } 
                else if ( !testLineFour ) 
                {
                    k = 1;
                    i = 0; // s4
                    j = 2;
                } 
                else
                {
                    k = 0;
                    i = 2; // s2
                    j = 1;
                }
                si = s.vrtx[i];
                sj = s.vrtx[j];
                sk = s.vrtx[k];

                if ( !hff2( s1, sj, sk ) )
                {
                    if ( !hff2( s1, sk, sj ) ) 
                    {
                        s.nvrtx = 3;
                        s.vrtx[2] = s.vrtx[3];
                        s.vrtx[1] = sj;
                        s.vrtx[0] = sk;
                        projectOnPlane( s1, sj, sk, v );
                    } 
                    else if ( !hff2( s1, sk, si ) )
                    {
                        s.nvrtx = 3;
                        s.vrtx[2] = s.vrtx[3];
                        s.vrtx[1] = si;
                        s.vrtx[0] = sk;
                        projectOnPlane( s1, sk, si, v );
                    } 
                    else 
                    {
                        s.nvrtx = 2;
                        s.vrtx[1] = s.vrtx[3];
                        s.vrtx[0] = sk;
                        projectOnLine( s1, sk, v );
                    }
                } 
                else if ( !hff2( s1, sj, si ) )
                {
                    s.nvrtx = 3;
                    s.vrtx[2] = s.vrtx[3];
                    s.vrtx[1] = si;
                    s.vrtx[0] = sj;
                    projectOnPlane( s1, si, sj, v );
                }
                else
                {
                    s.nvrtx = 2;
                    s.vrtx[1] = s.vrtx[3];
                    s.vrtx[0] = sj;
                    projectOnLine( s1, sj, v );
                }
            }
            break;
        default:
            printf("\nERROR:\tunhandled");
    }
}




// -----------------------------------------------------------------------------
// Distance subalgorithm using the signed volume method
static inline void subalgorithm( Simplex& s,
                                 Vector3& v)
{
    switch ( s.nvrtx )
    {
        case 4:
            S3D( s, v );
            break;
        case 3:
            S2D( s, v );
            break;
        case 2:
            S1D( s, v );
            break;
        default:
            printf("\nERROR:\t invalid simplex\n");
    }
}




/* ========================================================================== */
/*                            High-Level Methods                              */
/* ========================================================================== */
// Returns whether 2 convex shapes intersect using the GJK_SV algorithm
double closest_points_GJK_SV( Convex const& a, 
                              Convex const& b, 
                              Transform const& a2w,
                              Transform const& b2w, 
                              Point3& pa,
                              Point3& pb,
                              int& nbIter ) 
{
    unsigned int k = 0;                  // iteration counter
    int const maxk = 30;                 // maximum number of GJK iterations
    unsigned int i;
    
    /* Initialise search direction */
    // Vector3 c2c = b2w.getOrigin() - a2w.getOrigin();
    // Vector3 w = a2w( a.support( ( -c2c ) * a2w.getBasis() ) ) - 
    //             b2w( b.support( (  c2c ) * b2w.getBasis() ) );
    Vector3 w = a2w( a.support( Vector3Null ) ) - 
                b2w( b.support( Vector3Null ) );
    Vector3 v = w;
    Vector3 d = w;
    double dist2 = Norm2( v );
    double momentum = 0., oneMinusMomentum = 1.;
    bool acceleration = true;

    /* Initialise simplex */
    Simplex s = { 1, { Vector3( 0., 0., 0. ) } };
    s.vrtx[0] = w;

    /* Begin GJK iteration */
    do {
        k++;
        if ( acceleration )
        {
            momentum = k / ( k + 2. );
            oneMinusMomentum = 1. - momentum;
            d = momentum * d + 
                momentum * oneMinusMomentum * w +
                oneMinusMomentum * oneMinusMomentum * v;
            w = a2w( a.support( ( -d ) * a2w.getBasis() ) ) - 
                b2w( b.support( (  d ) * b2w.getBasis() ) );
            
            if ( dist2 - v * w <= sqrt( dist2 ) * rel_error1 / 2. )
            {
                if ( sqrt( dist2 ) * Norm( d ) - v * d <= abs_error1 )
                    break;
                w = a2w( a.support( ( -v ) * a2w.getBasis() ) ) - 
                    b2w( b.support( (  v ) * b2w.getBasis() ) );
                acceleration = false;
            }
        }
        else
        {
            w = a2w( a.support( ( -v ) * a2w.getBasis() ) ) - 
                b2w( b.support( (  v ) * b2w.getBasis() ) );

            double exeedtol_rel = dist2 - ( v * w );
            if ( exeedtol_rel <= sqrt( dist2 ) * rel_error1 / 2. || 
                 exeedtol_rel <= abs_error1 )
                break;
        }

        // Add new vertex to simplex
        i = s.nvrtx;
        s.vrtx[i] = w;
        s.nvrtx++;

        /* Invoke distance sub-algorithm */
        subalgorithm( s, v );
        dist2 = Norm2( v );

    } while ( ( s.nvrtx != 4 ) && ( k != maxk ) );
    pa = a.support( ( -v ) * a2w.getBasis() );
    pb = a.support( (  v ) * b2w.getBasis() );
    nbIter = k;
    return ( sqrt( dist2 ) );
}