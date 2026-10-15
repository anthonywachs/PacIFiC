#include "GJK.hh"
#include "MatrixMath.hh"
#include "MiscMath.hh"
#include "QuaternionMath.hh"

/* ============================================================================================== */
/* Johnson Low-Level Methods                                                                      */
/* ============================================================================================== */
template <typename T>
__HOSTDEVICE__ static INLINE void computeDet(const uint bits,
                                             const uint last,
                                             const uint last_bit,
                                             const uint all_bits,
                                             const Vector3<T> (&y)[4],
                                             T dp[4][4],
                                             T det[16][4])
{
    for(uint i = 0, bit = 1; i < 4; ++i, bit <<= 1)
        if(bits & bit)
            dp[i][last] = dp[last][i] = y[i] * y[last];
    dp[last][last] = y[last] * y[last];

    det[last_bit][last] = T(1);
    for(uint j = 0, sj = 1; j < 4; ++j, sj <<= 1)
    {
        if(bits & sj)
        {
            uint s2       = sj | last_bit;
            det[s2][j]    = dp[last][last] - dp[last][j];
            det[s2][last] = dp[j][j] - dp[j][last];
            for(uint k = 0, sk = 1; k < j; ++k, sk <<= 1)
            {
                if(bits & sk)
                {
                    int s3     = sk | s2;
                    det[s3][k] = det[s2][j] * (dp[j][j] - dp[j][k])
                                 + det[s2][last] * (dp[last][j] - dp[last][k]);
                    det[s3][j] = det[sk | last_bit][k] * (dp[k][k] - dp[k][j])
                                 + det[sk | last_bit][last] * (dp[last][k] - dp[last][j]);
                    det[s3][last] = det[sk | sj][k] * (dp[k][k] - dp[k][last])
                                    + det[sk | sj][j] * (dp[j][k] - dp[j][last]);
                }
            }
        }
    }

    if(all_bits == 15)
    {
        det[15][0] = det[14][1] * (dp[1][1] - dp[1][0]) + det[14][2] * (dp[2][1] - dp[2][0])
                     + det[14][3] * (dp[3][1] - dp[3][0]);
        det[15][1] = det[13][0] * (dp[0][0] - dp[0][1]) + det[13][2] * (dp[2][0] - dp[2][1])
                     + det[13][3] * (dp[3][0] - dp[3][1]);
        det[15][2] = det[11][0] * (dp[0][0] - dp[0][2]) + det[11][1] * (dp[1][0] - dp[1][2])
                     + det[11][3] * (dp[3][0] - dp[3][2]);
        det[15][3] = det[7][0] * (dp[0][0] - dp[0][3]) + det[7][1] * (dp[1][0] - dp[1][3])
                     + det[7][2] * (dp[2][0] - dp[2][3]);
    }
}

// -------------------------------------------------------------------------------------------------
template <typename T>
__HOSTDEVICE__ static INLINE bool valid(const uint s, const uint all_bits, const T det[16][4])
{
    for(uint i = 0, bit = 1; i < 4; ++i, bit <<= 1)
    {
        if(all_bits & bit)
        {
            if(s & bit)
            {
                if(det[s][i] < HIGHEPS<T>)
                    return (false);
            }
            else if(det[s | bit][i] > T(0))
                return (false);
        }
    }
    return (true);
}

// -------------------------------------------------------------------------------------------------
// Unified computeVector implementation that works for both algorithms
template <typename T, typename WeightType>
__HOSTDEVICE__ static INLINE void computeVector(const uint bits,
                                                const Vector3<T> (&y)[4],
                                                const WeightType& weights,
                                                Vector3<T>&       v)
{
    v.setValue(T(0), T(0), T(0));

    // Process based on weight type (array dimensions indicate algorithm)
    if constexpr(std::is_pointer<WeightType>::value
                 || (std::is_array<WeightType>::value && std::extent<WeightType, 0>::value > 4))
    {
        // Johnson algorithm (det[16][4])
        T sum = T(0);
        for(uint i = 0, bit = 1; i < 4; ++i, bit <<= 1)
        {
            if(bits & bit)
            {
                sum += weights[bits][i];
                v += weights[bits][i] * y[i];
            }
        }
        v *= T(1) / sum;
    }
    else
    {
        // SignedVolume algorithm (lambdas[4])
        for(uint i = 0; i < 4; ++i)
        {
            if(bits & (1 << i))
            {
                v += static_cast<T>(weights[i]) * y[i];
            }
        }
    }
}

// -------------------------------------------------------------------------------------------------
// Unified computePoints implementation that handles both Johnson and
// SignedVolume algorithms
template <typename T, typename WeightType>
__HOSTDEVICE__ static INLINE void computePoints(const uint bits,
                                                const Vector3<T> (&p)[4],
                                                const Vector3<T> (&q)[4],
                                                const WeightType& weights,
                                                Vector3<T>&       p1,
                                                Vector3<T>&       p2)
{
    T sum = T(0);
    p1.setValue(T(0), T(0), T(0));
    p2.setValue(T(0), T(0), T(0));

    // Process based on weight type (array dimensions indicate algorithm)
    if constexpr(std::is_pointer<WeightType>::value
                 || (std::is_array<WeightType>::value && std::extent<WeightType, 0>::value > 4))
    {
        // Johnson algorithm (det[16][4])
        for(uint i = 0, bit = 1; i < 4; ++i, bit <<= 1)
        {
            if(bits & bit)
            {
                sum += weights[bits][i];
                p1 += weights[bits][i] * p[i];
                p2 += weights[bits][i] * q[i];
            }
        }
        T s = T(1) / sum;
        p1 *= s;
        p2 *= s;
    }
    else  // SignedVolume algorithm (lambdas[4])
    {
        for(uint i = 0; i < 4; ++i)
        {
            if(bits & (1 << i))
            {
                p1 += static_cast<T>(weights[i]) * p[i];
                p2 += static_cast<T>(weights[i]) * q[i];
            }
        }
    }
}

// -------------------------------------------------------------------------------------------------
template <typename T>
__HOSTDEVICE__ static INLINE bool proper(const uint s, const T det[16][4])
{
    for(uint i = 0, bit = 1; i < 4; ++i, bit <<= 1)
        if((s & bit) && det[s][i] < HIGHEPS<T>)
            return (false);
    return (true);
}

// -------------------------------------------------------------------------------------------------
template <typename T>
__HOSTDEVICE__ static INLINE bool closest(uint&      bits,
                                          const uint last,
                                          const uint last_bit,
                                          const uint all_bits,
                                          const Vector3<T> (&y)[4],
                                          T           dp[4][4],
                                          T           det[16][4],
                                          Vector3<T>& v)
{
    uint s;
    computeDet(bits, last, last_bit, all_bits, y, dp, det);
    for(s = bits; s; --s)
    {
        if((s & bits) == s)
        {
            if(valid<T>(s | last_bit, all_bits, det))
            {
                bits = s | last_bit;
                computeVector<T>(bits, y, det, v);
                return (true);
            }
        }
    }
    if(valid<T>(last_bit, all_bits, det))
    {
        bits = last_bit;
        v    = y[last];
        return (true);
    }
    // Original GJK calls the backup procedure at this point.
    T min_dist2 = INFINITY;
    for(s = all_bits; s; --s)
    {
        if((s & all_bits) == s)
        {
            if(proper<T>(s, det))
            {
                Vector3<T> u;
                computeVector<T>(s, y, det, u);
                T dist2 = norm2(u);
                if(dist2 < min_dist2)
                {
                    min_dist2 = dist2;
                    bits      = s;
                    v         = u;
                }
            }
        }
    }
    return (false);
}

// -------------------------------------------------------------------------------------------------
// The next function is used for detecting degenerate cases that cause
// termination problems due to rounding errors.
template <typename T>
__HOSTDEVICE__ static INLINE bool
    degenerate(const uint bits, const Vector3<T> (&y)[4], const Vector3<T>& w)
{
    constexpr T err = HIGHEPS<T>;
    for(uint i = 0, bit = 1; i < 4; ++i, bit <<= 1)
    {
        if(bits & bit)
        {
            // Use ::fabs instead of fabs to avoid the template error
            const T yi2 = y[i] * y[i];
            if(::fabs(w * y[i] - yi2) < err * std::max(yi2, err))
                return (true);
        }
    }
    return (false);
}

// -------------------------------------------------------------------------------------------------
// For num_iterations > 1000
__HOSTDEVICE__
void catch_me()
{
    printf("closestPointsGJK: Exceeding 1000 iterations.\n");
}

// -------------------------------------------------------------------------------------------------
/** @brief Fast reciprocal square root: single hardware instruction on GPU
    (rsqrtf / rsqrt), standard 1/sqrt on CPU.
    @param x  value to invert-sqrt; undefined for x < 0 */
template <typename T>
__HOSTDEVICE__ static INLINE T inverseSqrt(T x)
{
#ifdef __CUDA_ARCH__
    if constexpr(std::is_same_v<T, float>)
        return rsqrtf(x);
    else
        return rsqrt(x);
#else
    return T(1) / sqrt(x);
#endif
}

/* ============================================================================================== */
/* Signed Volume Low-Level Methods                                                                */
/* ============================================================================================== */
template <typename T>
__HOSTDEVICE__ static INLINE uint compareSigns(T a, T b)
{
    // Maybe there's a faster way to deal with this set of operations?
    return static_cast<uint>(!((a > 0) ^ (b > 0)));
}

// -------------------------------------------------------------------------------------------------
template <typename T>
__HOSTDEVICE__ static INLINE void s1d(const Vector3<T> (&y)[4], uint& bits, T (&lambdas)[4])
{
    // Identify the appropriate indices
    bool s1_set = false;
    uint i1 = 0xffffffff, i2 = 0xffffffff;
    for(uint i = 0; i < 4; ++i)
    {
        if(bits & (1 << i))
        {
            if(s1_set)
            {
                i2 = i;
                break;
            }
            else
            {
                i1     = i;
                s1_set = true;
            }
        }
    }

    // Calculate the signed volume of the simplex.
    Vector3<T> t      = y[i2] - y[i1];
    uint       I      = 0;
    T          neg_tI = -t[0];

    if(fabs(t[1]) > fabs(neg_tI))
    {
        I      = 1;
        neg_tI = -t[1];
    }

    if(fabs(t[2]) > fabs(neg_tI))
    {
        I      = 2;
        neg_tI = -t[2];
    }

    T pI = (y[i2] * t) / norm2(t) * neg_tI + y[i2][I];

    // Identify the signed volume resulting from replacing each point by the
    // origin.
    T    C[2]                = {-y[i2][I] + pI, y[i1][I] - pI};
    uint sign_comparisons[2] = {compareSigns(neg_tI, C[0]), compareSigns(neg_tI, C[1])};

    // If all signed volumes are identical, the origin lies inside the simplex.
    if(sign_comparisons[0] + sign_comparisons[1] == 2)
    {
        lambdas[i1] = C[0] / neg_tI;
        lambdas[i2] = C[1] / neg_tI;
    }
    else
    {
        // The point to retain is the one whose sign matches. In the
        // first case, the origin lies past the first point.
        if(sign_comparisons[0])
        {
            bits &= ~(1 << i2);
            lambdas[i1] = T(1);
        }
        else
        {
            bits &= ~(1 << i1);
            lambdas[i2] = T(1);
        }
    }
}

// -------------------------------------------------------------------------------------------------
template <typename T>
__HOSTDEVICE__ static INLINE void s2d(const Vector3<T> (&y)[4], uint& bits, T (&lambdas)[4])
{
    uint counter = 0, point0_idx = 0, point1_idx = 0, point2_idx = 0;
    for(uint i = 0; i < 4; ++i)
    {
        if(bits & (1 << i))
        {
            if(counter == 0)
                point0_idx = i;
            else if(counter == 1)
                point1_idx = i;
            else
                point2_idx = i;
            counter += 1;
        }
    }

    Vector3<T> n  = (y[point1_idx] - y[point0_idx]) ^ (y[point2_idx] - y[point0_idx]);
    Vector3<T> p0 = (y[point0_idx] * n / norm2(n)) * n;

    // Choose maximum area plane to project onto.
    // Make sure to store the *signed* area of the plane.
    // This loop is unrolled to save a few extra ops (assigning
    // an initial area of zero, an extra abs, etc)
    uint idx_x  = 1;
    uint idx_y  = 2;
    T    mu_max = (y[point1_idx][1] * y[point2_idx][2] + y[point0_idx][1] * y[point1_idx][2]
                + y[point2_idx][1] * y[point0_idx][2] - y[point1_idx][1] * y[point0_idx][2]
                - y[point2_idx][1] * y[point1_idx][2] - y[point0_idx][1] * y[point2_idx][2]);

    // This term is multiplied by -1.
    T mu = (y[point1_idx][2] * y[point0_idx][0] + y[point2_idx][2] * y[point1_idx][0]
            + y[point0_idx][2] * y[point2_idx][0] - y[point1_idx][2] * y[point2_idx][0]
            - y[point0_idx][2] * y[point1_idx][0] - y[point2_idx][2] * y[point0_idx][0]);
    if(fabs(mu) > fabs(mu_max))
    {
        mu_max = mu;
        idx_x  = 0;
    }

    mu = (y[point1_idx][0] * y[point2_idx][1] + y[point0_idx][0] * y[point1_idx][1]
          + y[point2_idx][0] * y[point0_idx][1] - y[point1_idx][0] * y[point0_idx][1]
          - y[point2_idx][0] * y[point1_idx][1] - y[point0_idx][0] * y[point2_idx][1]);
    if(fabs(mu) > fabs(mu_max))
    {
        mu_max = mu;
        idx_x  = 0;
        idx_y  = 1;
    }

    // Compute the signed areas of each of the simplices formed by replacing an
    // index with a projection of the origin onto the area in this plane
    T    C[3]                = {T(0)};
    bool sign_comparisons[3] = {false};

    C[0]                = (p0[idx_x] * y[point1_idx][idx_y] + p0[idx_y] * y[point2_idx][idx_x]
            + y[point1_idx][idx_x] * y[point2_idx][idx_y] - p0[idx_x] * y[point2_idx][idx_y]
            - p0[idx_y] * y[point1_idx][idx_x] - y[point2_idx][idx_x] * y[point1_idx][idx_y]);
    sign_comparisons[0] = compareSigns(mu_max, C[0]);

    C[1]                = (p0[idx_x] * y[point2_idx][idx_y] + p0[idx_y] * y[point0_idx][idx_x]
            + y[point2_idx][idx_x] * y[point0_idx][idx_y] - p0[idx_x] * y[point0_idx][idx_y]
            - p0[idx_y] * y[point2_idx][idx_x] - y[point0_idx][idx_x] * y[point2_idx][idx_y]);
    sign_comparisons[1] = compareSigns(mu_max, C[1]);

    C[2]                = (p0[idx_x] * y[point0_idx][idx_y] + p0[idx_y] * y[point1_idx][idx_x]
            + y[point0_idx][idx_x] * y[point1_idx][idx_y] - p0[idx_x] * y[point1_idx][idx_y]
            - p0[idx_y] * y[point0_idx][idx_x] - y[point1_idx][idx_x] * y[point0_idx][idx_y]);
    sign_comparisons[2] = compareSigns(mu_max, C[2]);

    if(sign_comparisons[0] + sign_comparisons[1] + sign_comparisons[2] == 3)
    {
        lambdas[point0_idx] = C[0] / mu_max;
        lambdas[point1_idx] = C[1] / mu_max;
        lambdas[point2_idx] = C[2] / mu_max;
    }
    else
    {
        T          d = std::numeric_limits<T>::max();
        Vector3<T> new_point;
        uint       new_bits = 0;
        for(uint j = 0; j < 3; ++j)
        {
            if(!sign_comparisons[j])
            {
                uint new_used = bits;
                // Test removal of the current point.
                if(j == 0)
                    new_used &= ~(1 << point0_idx);
                else if(j == 1)
                    new_used &= ~(1 << point1_idx);
                else
                    new_used &= ~(1 << point2_idx);

                T new_lambdas[4] = {T(0)};

                s1d(y, new_used, new_lambdas);
                // Consider resetting in place if possible.
                new_point[0] = 0;
                new_point[1] = 0;
                new_point[2] = 0;
                for(uint i = 0; i < 4; ++i)
                {
                    if(new_used & (1 << i))
                        new_point += new_lambdas[i] * y[i];
                }
                T d_star = new_point * new_point;
                if(d_star < d)
                {
                    new_bits = new_used;
                    d        = d_star;
                    for(uint i = 0; i < 4; ++i)
                        lambdas[i] = new_lambdas[i];
                }
            }
        }
        bits = new_bits;
    }
}

// -------------------------------------------------------------------------------------------------
template <typename T>
__HOSTDEVICE__ static INLINE void s3d(const Vector3<T> (&y)[4], uint& bits, T (&lambdas)[4])
{
    T C[4] = {0.};

    // Compute all minors and the total determinant of the matrix M,
    // which is the transpose of the y matrix with an extra row of
    // ones at the bottom. Since the indexing is nontrivial and the
    // array is small (and we can save on some negation), all the
    // computations are done directly rather than with a loop.
    // C[0] and C[2] are negated due to the (-1)^(i+j+1) prefactor,
    // where i is always 4 because we're expanding about the 4th row.
    C[0] = y[3][0] * y[2][1] * y[1][2] + y[2][0] * y[1][1] * y[3][2] + y[1][0] * y[3][1] * y[2][2]
           - y[1][0] * y[2][1] * y[3][2] - y[2][0] * y[3][1] * y[1][2]
           - y[3][0] * y[1][1] * y[2][2];
    C[1] = y[0][0] * y[2][1] * y[3][2] + y[2][0] * y[3][1] * y[0][2] + y[3][0] * y[0][1] * y[2][2]
           - y[3][0] * y[2][1] * y[0][2] - y[2][0] * y[0][1] * y[3][2]
           - y[0][0] * y[3][1] * y[2][2];
    C[2] = y[3][0] * y[1][1] * y[0][2] + y[1][0] * y[0][1] * y[3][2] + y[0][0] * y[3][1] * y[1][2]
           - y[0][0] * y[1][1] * y[3][2] - y[1][0] * y[3][1] * y[0][2]
           - y[3][0] * y[0][1] * y[1][2];
    C[3] = y[0][0] * y[1][1] * y[2][2] + y[1][0] * y[2][1] * y[0][2] + y[2][0] * y[0][1] * y[1][2]
           - y[2][0] * y[1][1] * y[0][2] - y[1][0] * y[0][1] * y[2][2]
           - y[0][0] * y[2][1] * y[1][2];
    T dM = C[0] + C[1] + C[2] + C[3];

    uint sign_comparisons[4] = {0};
    sign_comparisons[0]      = compareSigns(dM, C[0]);
    sign_comparisons[1]      = compareSigns(dM, C[1]);
    sign_comparisons[2]      = compareSigns(dM, C[2]);
    sign_comparisons[3]      = compareSigns(dM, C[3]);

    if((sign_comparisons[0] + sign_comparisons[1] + sign_comparisons[2] + sign_comparisons[3]) == 4)
    {
        for(uint i = 0; i < 4; ++i)
            lambdas[i] = C[i] / dM;
    }
    else
    {
        T          d = std::numeric_limits<T>::max(), d_star = T(0);
        Vector3<T> new_point;
        uint       new_bits = 0;
        for(uint j = 0; j < 4; ++j)
        {
            if(!sign_comparisons[j])
            {
                // Test removal of the current point.
                uint new_used = bits;
                new_used &= ~(1 << j);
                T new_lambdas[4] = {T(0)};

                s2d(y, new_used, new_lambdas);

                new_point = Vector3<T>();
                for(uint i = 0; i < 4; ++i)
                {
                    if(new_used & (1 << i))
                        new_point += new_lambdas[i] * y[i];
                }
                d_star = new_point * new_point;
                if(d_star < d)
                {
                    new_bits = new_used;
                    d        = d_star;
                    for(uint i = 0; i < 4; ++i)
                        lambdas[i] = new_lambdas[i];
                }
            }
        }
        bits = new_bits;
    }
}

// -------------------------------------------------------------------------------------------------
template <typename T>
__HOSTDEVICE__ static INLINE void
    sv_subalgorithm(const Vector3<T> (&y)[4], uint& bits, T (&lambdas)[4], Vector3<T>& v)
{
    // The y array is never modified by this function.  The bits may be
    // modified if necessary, and the lambdas will be updated.  All the other
    // functions (if they need to make deeper calls e.g. s3d->s2d) will have to
    // make copies of bits to avoid overwriting that data incorrectly.
#ifdef __CUDA_ARCH__
    const uint num_used = __popc(bits);
#else
    const uint num_used = static_cast<uint>(__builtin_popcount(bits));
#endif

    // Start with the most common cases.
    if(num_used == 1)
    {
        for(uint i = 0; i < 4; ++i)
        {
            if(bits & (1 << i))
                lambdas[i] = T(1);
        }
    }
    else if(num_used == 2)
        s1d(y, bits, lambdas);
    else if(num_used == 3)
        s2d(y, bits, lambdas);
    else
        s3d(y, bits, lambdas);

    computeVector<T>(bits, y, lambdas, v);
}

/* ============================================================================================== */
/* High-Level Methods                                                                             */
/* ============================================================================================== */
// Returns whether 2 convex shapes intersect using the GJK algorithm - relative
// transformation
template <typename T>
__HOSTDEVICE__ bool intersectGJK(const Convex<T>& a, const Convex<T>& b, const Transform3<T>& b2a)
{
    uint       bits     = 0;         // identifies current simplex
    uint       last     = 0;         // identifies last found support point
    uint       last_bit = 0;         // last_bit = 1<<last
    uint       all_bits = 0;         // all_bits = bits|last_bit
    Vector3<T> y[4];                 // support points of A-B in world
    T          det[16][4] = {T(0)};  // cached sub-determinants
    T          dp[4][4]   = {T(0)};

    Vector3<T> v(b2a.getOrigin());
    Vector3<T> w;
    T          prod;

    do
    {
        last     = 0;
        last_bit = 1;
        while(bits & last_bit)
        {
            ++last;
            last_bit <<= 1;
        }
        w    = a.support(-v) - b2a(b.support(v * b2a.getBasis()));
        prod = v * w;
        if(prod > T(0) || fabs(prod) < HIGHEPS<T>)
            return (false);
        if(degenerate(all_bits, y, w))
            return (false);
        y[last]  = w;
        all_bits = bits | last_bit;
        if(!closest(bits, last, last_bit, all_bits, y, dp, det, v))
            return (false);
    } while(bits < 15 && !isApproxZero(v));
    return (true);
}

// -------------------------------------------------------------------------------------------------
// Returns whether 2 convex shapes intersect using the GJK algorithm
template <typename T>
__HOSTDEVICE__ bool intersectGJK(const Convex<T>&     a,
                                 const Convex<T>&     b,
                                 const Transform3<T>& a2w,
                                 const Transform3<T>& b2w)
{
    uint       bits     = 0;         // identifies current simplex
    uint       last     = 0;         // identifies last found support point
    uint       last_bit = 0;         // last_bit = 1<<last
    uint       all_bits = 0;         // all_bits = bits|last_bit
    Vector3<T> y[4];                 // support points of A-B in world
    T          det[16][4] = {T(0)};  // cached sub-determinants
    T          dp[4][4]   = {T(0)};

    Vector3<T> v(b2w.getOrigin() - a2w.getOrigin());
    Vector3<T> w;
    T          prod;

    do
    {
        last     = 0;
        last_bit = 1;
        while(bits & last_bit)
        {
            ++last;
            last_bit <<= 1;
        }
        w    = a2w(a.support((-v) * a2w.getBasis())) - b2w(b.support(v * b2w.getBasis()));
        prod = v * w;
        if(prod > T(0) || fabs(prod) < HIGHEPS<T>)
            return (false);
        if(degenerate(all_bits, y, w))
            return (false);
        y[last]  = w;
        all_bits = bits | last_bit;
        if(!closest(bits, last, last_bit, all_bits, y, dp, det, v))
            return (false);
    } while(bits < 15 && !isApproxZero(v));
    return (true);
}

// -------------------------------------------------------------------------------------------------
// Returns whether 2 convex shapes intersect using the GJK algorithm - relative
// transformation
template <typename T>
__HOSTDEVICE__ bool intersectGJK(const Convex<T>&     a,
                                 const Convex<T>&     b,
                                 const Vector3<T>&    v_b2a,
                                 const Quaternion<T>& q_b2a)
{
    uint       bits     = 0;         // identifies current simplex
    uint       last     = 0;         // identifies last found support point
    uint       last_bit = 0;         // last_bit = 1<<last
    uint       all_bits = 0;         // all_bits = bits|last_bit
    Vector3<T> y[4];                 // support points of A-B in world
    T          det[16][4] = {T(0)};  // cached sub-determinants
    T          dp[4][4]   = {T(0)};

    Vector3<T> v(v_b2a);
    Vector3<T> w;
    T          prod;

    Vector3<T> p, q;
    do
    {
        last     = 0;
        last_bit = 1;
        while(bits & last_bit)
        {
            ++last;
            last_bit <<= 1;
        }
        // w = a.support(-v) - q_b2a(b.support(v ^ q_b2a));
        p = a.support(-v);
        q = b.support(q_b2a << v);
        FusedMinkowskiDifference(p, q, v_b2a, q_b2a, w);
        prod = v * w;
        if(prod > T(0) || fabs(prod) < HIGHEPS<T>)
            return (false);
        if(degenerate(all_bits, y, w))
            return (false);
        y[last]  = w;
        all_bits = bits | last_bit;
        if(!closest(bits, last, last_bit, all_bits, y, dp, det, v))
            return (false);
    } while(bits < 15 && !isApproxZero(v));
    return (true);
}

// -------------------------------------------------------------------------------------------------
// Returns whether 2 convex shapes intersect using the GJK algorithm
template <typename T>
__HOSTDEVICE__ bool intersectGJK(const Convex<T>&     a,
                                 const Convex<T>&     b,
                                 const Vector3<T>&    v_a2w,
                                 const Vector3<T>&    v_b2w,
                                 const Quaternion<T>& q_a2w,
                                 const Quaternion<T>& q_b2w)
{
    uint       bits     = 0;         // identifies current simplex
    uint       last     = 0;         // identifies last found support point
    uint       last_bit = 0;         // last_bit = 1<<last
    uint       all_bits = 0;         // all_bits = bits|last_bit
    Vector3<T> y[4];                 // support points of A-B in world
    T          det[16][4] = {T(0)};  // cached sub-determinants
    T          dp[4][4]   = {T(0)};

    Vector3<T> v(v_b2w - v_a2w);
    Vector3<T> w;
    T          prod;

    Vector3<T> p, q;
    do
    {
        last     = 0;
        last_bit = 1;
        while(bits & last_bit)
        {
            ++last;
            last_bit <<= 1;
        }
        // w = q_a2w(a.support((-v) ^ q_a2w)) - q_b2w(b.support(v ^ q_b2w));
        p = a.support(q_a2w << (-v));
        q = b.support(q_b2w << v);
        FusedMinkowskiDifference(p, q, v_a2w, v_b2w, q_a2w, q_b2w, w);
        prod = v * w;
        if(prod > T(0) || fabs(prod) < HIGHEPS<T>)
            return (false);
        if(degenerate(all_bits, y, w))
            return (false);
        y[last]  = w;
        all_bits = bits | last_bit;
        if(!closest(bits, last, last_bit, all_bits, y, dp, det, v))
            return (false);
    } while(bits < 15 && !isApproxZero(v));
    return (true);
}

// -------------------------------------------------------------------------------------------------
// Johnson implementation of closest points algorithm
template <typename T, bool Acceleration>
__HOSTDEVICE__ T computeClosestPoints_GJK_Johnson(const Convex<T>&     a,
                                                  const Convex<T>&     b,
                                                  const Transform3<T>& b2a,
                                                  const T              crustA,
                                                  const T              crustB,
                                                  Vector3<T>&          pa,
                                                  Vector3<T>&          pb,
                                                  uint&                nbIter)
{
    // Constants
    constexpr T    relError = EPS<T>;               // relative tolerance
    constexpr T    absError = T(1.e-4 * relError);  // absolute tolerance
    constexpr uint MAXITERS = 1000;

    // Johnson-specific variables
    uint       bits     = 0;            // identifies current simplex
    uint       last     = 0;            // identifies last found support point
    uint       last_bit = 0;            // last_bit = 1<<last
    uint       all_bits = 0;            // all_bits = bits|last_bit
    Vector3<T> p[4];                    // support points of A in local
    Vector3<T> q[4];                    // support points of B in local
    Vector3<T> y[4];                    // support points of A-B in world
    T          det[16][4]    = {T(0)};  // cached sub-determinants
    T          dp[4][4]      = {T(0)};  // cached dot products
    T          mu            = T(0);    // optimality gap
    uint       numIterations = 0;       // No. iterations

    // Acceleration-specific variables
    T momentum = T(0);  // Only used if Acceleration is true

    // Initializing vectors
    Vector3<T> v(b2a.getOrigin());
    Vector3<T> w;
    T          dist2   = norm2(v);
    T          invDist = (dist2 > T(0)) ? inverseSqrt(dist2) : T(0);
    T          dist    = dist2 * invDist;

    while(bits < 15 && dist > HIGHEPS<T> && numIterations < MAXITERS)
    {
        // Updating the bits
        last     = 0;
        last_bit = 1;
        while(bits & last_bit)
        {
            ++last;
            last_bit <<= 1;
        }

        // Support points (eroded by crust in local frame, reusing local dir)
        const Vector3<T> vB = v * b2a.getBasis();
        p[last]             = a.support(-v, crustA, invDist);
        q[last]             = b.support(vB, crustB, invDist);
        w                   = p[last] - b2a(q[last]);

        // termination criteria -- optimality gap
        mu = dist - (v * w) * invDist;
        if(mu <= dist * relError || mu < absError)
            break;

        // termination criteria -- degenerate case
        if(degenerate(all_bits, y, w))
            break;

        // if not terminated, get ready for the next iteration
        y[last]  = w;
        all_bits = bits | last_bit;
        if(!closest(bits, last, last_bit, all_bits, y, dp, det, v))
            break;

        ++numIterations;
        dist2   = norm2(v);
        invDist = (dist2 > T(0)) ? inverseSqrt(dist2) : T(0);
        dist    = dist2 * invDist;
    }

    computePoints<T>(bits, p, q, det, pa, pb);

    if(numIterations >= 1000)
        catch_me();
    else
        nbIter = numIterations;

    return (dist);
}

// -------------------------------------------------------------------------------------------------
// SignedVolume implementation of closest points algorithm
template <typename T, bool Acceleration>
__HOSTDEVICE__ T computeClosestPoints_GJK_SignedVolume(const Convex<T>&     a,
                                                       const Convex<T>&     b,
                                                       const Transform3<T>& b2a,
                                                       const T              crustA,
                                                       const T              crustB,
                                                       Vector3<T>&          pa,
                                                       Vector3<T>&          pb,
                                                       uint&                nbIter)
{
    // Constants
    constexpr T    relError = EPS<T>;               // relative tolerance
    constexpr T    absError = T(1.e-4 * relError);  // absolute tolerance
    constexpr uint MAXITERS = 1000;                 // Maximum iterations

    // SignedVolume-specific variables
    uint bits = 0;  // identifies current simplex
    uint last = 0;  // identifies last found support point

    Vector3<T> p[4];                 // support points of A in local
    Vector3<T> q[4];                 // support points of B in local
    Vector3<T> y[4];                 // support points of A-B in world
    T          lambdas[4] = {T(0)};  // Weights

    T   mu            = T(0);  // optimality gap
    int numIterations = 0;     // No. iterations

    // Acceleration-specific variables
    T momentum = T(0);  // Only used if Acceleration is true

    // Initializing vectors
    Vector3<T> v(b2a.getOrigin());
    Vector3<T> w;
    T          dist2   = norm2(v);
    T          invDist = (dist2 > T(0)) ? inverseSqrt(dist2) : T(0);
    T          dist    = dist2 * invDist;

    while(bits < 15 && dist > HIGHEPS<T> && numIterations < MAXITERS)
    {
        // Updating the bits
        for(uint new_index = 0; new_index < 4; ++new_index)
        {
            // At least one of these must be empty, otherwise overlap.
            if(!(bits & (1 << new_index)))
            {
                last = new_index;
                break;
            }
        }

        // Support points (eroded by crust in local frame, reusing local dir)
        const Vector3<T> vB = v * b2a.getBasis();
        p[last]             = a.support(-v, crustA, invDist);
        q[last]             = b.support(vB, crustB, invDist);
        w                   = p[last] - b2a(q[last]);

        // termination criteria -- optimality gap
        mu = dist - (v * w) * invDist;
        if(mu <= dist * relError || mu < absError)
            break;

        // termination criteria -- degenerate case
        if(degenerate(bits, y, w))
            break;

        // if not terminated, get ready for the next iteration
        y[last] = w;
        bits |= (1 << last);
        sv_subalgorithm(y, bits, lambdas, v);

        ++numIterations;
        dist2   = norm2(v);
        invDist = (dist2 > T(0)) ? inverseSqrt(dist2) : T(0);
        dist    = dist2 * invDist;
    }

    computePoints<T>(bits, p, q, lambdas, pa, pb);

    if(numIterations >= 1000)
        catch_me();
    else
        nbIter = numIterations;

    return (dist);
}

// -------------------------------------------------------------------------------------------------
// Dispatcher function that routes to the appropriate implementation
template <typename T, GJKType GJKType, bool Acceleration>
__HOSTDEVICE__ T computeClosestPoints_GJK(const Convex<T>&     a,
                                          const Convex<T>&     b,
                                          const Transform3<T>& b2a,
                                          const T              crustA,
                                          const T              crustB,
                                          Vector3<T>&          pa,
                                          Vector3<T>&          pb,
                                          uint&                nbIter)
{
    static_assert(GJKType == GJKType::JOHNSON || GJKType == GJKType::SIGNEDVOLUME,
                  "GJKType must be either Johnson or SignedVolume");

    if constexpr(GJKType == GJKType::JOHNSON)
    {
        return computeClosestPoints_GJK_Johnson<T, Acceleration>(a,
                                                                 b,
                                                                 b2a,
                                                                 crustA,
                                                                 crustB,
                                                                 pa,
                                                                 pb,
                                                                 nbIter);
    }
    else
    {
        return computeClosestPoints_GJK_SignedVolume<T, Acceleration>(a,
                                                                      b,
                                                                      b2a,
                                                                      crustA,
                                                                      crustB,
                                                                      pa,
                                                                      pb,
                                                                      nbIter);
    }
}

// -------------------------------------------------------------------------------------------------
// Johnson implementation for two transforms
template <typename T, bool Acceleration>
__HOSTDEVICE__ T computeClosestPoints_GJK_Johnson(const Convex<T>&     a,
                                                  const Convex<T>&     b,
                                                  const Transform3<T>& a2w,
                                                  const Transform3<T>& b2w,
                                                  const T              crustA,
                                                  const T              crustB,
                                                  Vector3<T>&          pa,
                                                  Vector3<T>&          pb,
                                                  uint&                nbIter)
{
    // Constants
    constexpr T relError = EPS<T>;               // relative tolerance
    constexpr T absError = T(1.e-4 * relError);  // absolute tolerance

    // Johnson-specific variables
    uint bits     = 0;  // identifies current simplex
    uint last     = 0;  // identifies last found support point
    uint last_bit = 0;  // last_bit = 1<<last
    uint all_bits = 0;  // all_bits = bits|last_bit

    Vector3<T> p[4];                 // support points of A in local
    Vector3<T> q[4];                 // support points of B in local
    Vector3<T> y[4];                 // support points of A-B in world
    T          det[16][4] = {T(0)};  // cached sub-determinants
    T          dp[4][4]   = {T(0)};  // cached dot products

    T   mu            = T(0);  // optimality gap
    int numIterations = 0;     // No. iterations

    // Acceleration-specific variables
    T momentum = T(0);  // Only used if Acceleration is true

    // Initializing vectors
    Vector3<T> v(a2w.getOrigin() - b2w.getOrigin());
    Vector3<T> w;
    T          dist2   = norm2(v);
    T          invDist = (dist2 > T(0)) ? inverseSqrt(dist2) : T(0);
    T          dist    = dist2 * invDist;

    while(bits < 15 && dist > HIGHEPS<T> && numIterations < 1000)
    {
        // Updating the bits
        last     = 0;
        last_bit = 1;
        while(bits & last_bit)
        {
            ++last;
            last_bit <<= 1;
        }

        // Support points (eroded by crust in local frame, reusing local dirs)
        const Vector3<T> vA = (-v) * a2w.getBasis();
        const Vector3<T> vB = v * b2w.getBasis();
        p[last]             = a.support(vA, crustA, invDist);
        q[last]             = b.support(vB, crustB, invDist);
        w                   = a2w(p[last]) - b2w(q[last]);

        // termination criteria -- optimality gap
        mu = dist - (v * w) * invDist;
        if(mu <= dist * relError || mu < absError)
            break;

        // termination criteria -- degenerate case
        if(degenerate(all_bits, y, w))
            break;

        // if not terminated, get ready for the next iteration
        y[last]  = w;
        all_bits = bits | last_bit;
        if(!closest(bits, last, last_bit, all_bits, y, dp, det, v))
            break;

        ++numIterations;
        dist2   = norm2(v);
        invDist = (dist2 > T(0)) ? inverseSqrt(dist2) : T(0);
        dist    = dist2 * invDist;
    }

    computePoints<T>(bits, p, q, det, pa, pb);

    if(numIterations >= 1000)
        catch_me();
    else
        nbIter = numIterations;

    return (dist);
}

// -------------------------------------------------------------------------------------------------
// SignedVolume implementation for two transforms
template <typename T, bool Acceleration>
__HOSTDEVICE__ T computeClosestPoints_GJK_SignedVolume(const Convex<T>&     a,
                                                       const Convex<T>&     b,
                                                       const Transform3<T>& a2w,
                                                       const Transform3<T>& b2w,
                                                       const T              crustA,
                                                       const T              crustB,
                                                       Vector3<T>&          pa,
                                                       Vector3<T>&          pb,
                                                       uint&                nbIter)
{
    // Constants
    constexpr T relError = EPS<T>;               // relative tolerance
    constexpr T absError = T(1.e-4 * relError);  // absolute tolerance

    // SignedVolume-specific variables
    uint bits = 0;  // identifies current simplex
    uint last = 0;  // identifies last found support point

    Vector3<T> p[4];                 // support points of A in local
    Vector3<T> q[4];                 // support points of B in local
    Vector3<T> y[4];                 // support points of A-B in world
    T          lambdas[4] = {T(0)};  // Weights

    T   mu            = T(0);  // optimality gap
    int numIterations = 0;     // No. iterations

    // Acceleration-specific variables
    T momentum = T(0);  // Only used if Acceleration is true

    // Initializing vectors
    Vector3<T> v(a2w.getOrigin() - b2w.getOrigin());
    Vector3<T> w;
    T          dist2   = norm2(v);
    T          invDist = (dist2 > T(0)) ? inverseSqrt(dist2) : T(0);
    T          dist    = dist2 * invDist;

    while(bits < 15 && dist > HIGHEPS<T> && numIterations < 1000)
    {
        // Updating the bits
        for(uint new_index = 0; new_index < 4; ++new_index)
        {
            // At least one of these must be empty, otherwise overlap.
            if(!(bits & (1 << new_index)))
            {
                last = new_index;
                break;
            }
        }

        // Support points (eroded by crust in local frame, reusing local dirs)
        const Vector3<T> vA = (-v) * a2w.getBasis();
        const Vector3<T> vB = v * b2w.getBasis();
        p[last]             = a.support(vA, crustA, invDist);
        q[last]             = b.support(vB, crustB, invDist);
        w                   = a2w(p[last]) - b2w(q[last]);

        // termination criteria -- optimality gap
        mu = dist - (v * w) * invDist;
        if(mu <= dist * relError || mu < absError)
            break;

        // termination criteria -- degenerate case
        if(degenerate(bits, y, w))
            break;

        // if not terminated, get ready for the next iteration
        y[last] = w;
        bits |= (1 << last);
        sv_subalgorithm(y, bits, lambdas, v);

        ++numIterations;
        dist2   = norm2(v);
        invDist = (dist2 > T(0)) ? inverseSqrt(dist2) : T(0);
        dist    = dist2 * invDist;
    }

    computePoints<T>(bits, p, q, lambdas, pa, pb);

    if(numIterations >= 1000)
        catch_me();
    else
        nbIter = numIterations;

    return (dist);
}

// -------------------------------------------------------------------------------------------------
// Dispatcher function for two transforms
template <typename T, GJKType GJKType, bool Acceleration>
__HOSTDEVICE__ T computeClosestPoints_GJK(const Convex<T>&     a,
                                          const Convex<T>&     b,
                                          const Transform3<T>& a2w,
                                          const Transform3<T>& b2w,
                                          const T              crustA,
                                          const T              crustB,
                                          Vector3<T>&          pa,
                                          Vector3<T>&          pb,
                                          uint&                nbIter)
{
    static_assert(GJKType == GJKType::JOHNSON || GJKType == GJKType::SIGNEDVOLUME,
                  "GJKType must be either Johnson or SignedVolume");

    if constexpr(GJKType == GJKType::JOHNSON)
    {
        return computeClosestPoints_GJK_Johnson<T, Acceleration>(a,
                                                                 b,
                                                                 a2w,
                                                                 b2w,
                                                                 crustA,
                                                                 crustB,
                                                                 pa,
                                                                 pb,
                                                                 nbIter);
    }
    else
    {
        return computeClosestPoints_GJK_SignedVolume<T, Acceleration>(a,
                                                                      b,
                                                                      a2w,
                                                                      b2w,
                                                                      crustA,
                                                                      crustB,
                                                                      pa,
                                                                      pb,
                                                                      nbIter);
    }
}

// -------------------------------------------------------------------------------------------------
// Johnson implementation for Vector/Quaternion
template <typename T, bool Acceleration>
__HOSTDEVICE__ T computeClosestPoints_GJK_Johnson(const Convex<T>&     a,
                                                  const Convex<T>&     b,
                                                  const Vector3<T>&    v_b2a,
                                                  const Quaternion<T>& q_b2a,
                                                  const T              crustA,
                                                  const T              crustB,
                                                  Vector3<T>&          pa,
                                                  Vector3<T>&          pb,
                                                  uint&                nbIter)
{
    // Constants
    constexpr T relError = EPS<T>;               // relative tolerance
    constexpr T absError = T(1.e-4 * relError);  // absolute tolerance

    // Johnson-specific variables
    uint bits     = 0;  // identifies current simplex
    uint last     = 0;  // identifies last found support point
    uint last_bit = 0;  // last_bit = 1<<last
    uint all_bits = 0;  // all_bits = bits|last_bit

    Vector3<T> p[4];                 // support points of A in local
    Vector3<T> q[4];                 // support points of B in local
    Vector3<T> y[4];                 // support points of A-B in world
    T          det[16][4] = {T(0)};  // cached sub-determinants
    T          dp[4][4]   = {T(0)};  // cached dot products

    T   mu            = T(0);  // optimality gap
    int numIterations = 0;     // No. iterations

    // Acceleration-specific variables
    T momentum = T(0);  // Only used if Acceleration is true

    // Initializing vectors
    Vector3<T> v(v_b2a);
    Vector3<T> w;
    T          dist2   = norm2(v);
    T          invDist = (dist2 > T(0)) ? inverseSqrt(dist2) : T(0);
    T          dist    = dist2 * invDist;

    while(bits < 15 && dist > HIGHEPS<T> && numIterations < 1000)
    {
        // Updating the bits
        last     = 0;
        last_bit = 1;
        while(bits & last_bit)
        {
            ++last;
            last_bit <<= 1;
        }

        // Support points (eroded by crust in local frame, reusing local dirs)
        const Vector3<T> vB = q_b2a << v;
        p[last]             = a.support(-v, crustA, invDist);
        q[last]             = b.support(vB, crustB, invDist);
        FusedMinkowskiDifference(p[last], q[last], v_b2a, q_b2a, w);

        // termination criteria -- optimality gap
        mu = dist - (v * w) * invDist;
        if(mu <= dist * relError || mu < absError)
            break;

        // termination criteria -- degenerate case
        if(degenerate(all_bits, y, w))
            break;

        // if not terminated, get ready for the next iteration
        y[last]  = w;
        all_bits = bits | last_bit;
        if(!closest(bits, last, last_bit, all_bits, y, dp, det, v))
            break;

        ++numIterations;
        dist2   = norm2(v);
        invDist = (dist2 > T(0)) ? inverseSqrt(dist2) : T(0);
        dist    = dist2 * invDist;
    }

    computePoints<T>(bits, p, q, det, pa, pb);

    if(numIterations >= 1000)
        catch_me();
    else
        nbIter = numIterations;

    return (dist);
}

// -------------------------------------------------------------------------------------------------
// SignedVolume implementation for Vector/Quaternion
template <typename T, bool Acceleration>
__HOSTDEVICE__ T computeClosestPoints_GJK_SignedVolume(const Convex<T>&     a,
                                                       const Convex<T>&     b,
                                                       const Vector3<T>&    v_b2a,
                                                       const Quaternion<T>& q_b2a,
                                                       const T              crustA,
                                                       const T              crustB,
                                                       Vector3<T>&          pa,
                                                       Vector3<T>&          pb,
                                                       uint&                nbIter)
{
    // Constants
    constexpr T relError = EPS<T>;               // relative tolerance
    constexpr T absError = T(1.e-4 * relError);  // absolute tolerance

    // SignedVolume-specific variables
    uint bits = 0;  // identifies current simplex
    uint last = 0;  // identifies last found support point

    Vector3<T> p[4];                 // support points of A in local
    Vector3<T> q[4];                 // support points of B in local
    Vector3<T> y[4];                 // support points of A-B in world
    T          lambdas[4] = {T(0)};  // Weights

    T   mu            = T(0);  // optimality gap
    int numIterations = 0;     // No. iterations

    // Acceleration-specific variables
    T momentum = T(0);  // Only used if Acceleration is true

    // Initializing vectors
    Vector3<T> v(v_b2a);
    Vector3<T> w;
    T          dist2   = norm2(v);
    T          invDist = (dist2 > T(0)) ? inverseSqrt(dist2) : T(0);
    T          dist    = dist2 * invDist;

    while(bits < 15 && dist > HIGHEPS<T> && numIterations < 1000)
    {
        // Updating the bits
        for(uint new_index = 0; new_index < 4; ++new_index)
        {
            // At least one of these must be empty, otherwise overlap.
            if(!(bits & (1 << new_index)))
            {
                last = new_index;
                break;
            }
        }

        // Support points (eroded by crust in local frame, reusing local dirs)
        const Vector3<T> vB = q_b2a << v;
        p[last]             = a.support(-v, crustA, invDist);
        q[last]             = b.support(vB, crustB, invDist);
        FusedMinkowskiDifference(p[last], q[last], v_b2a, q_b2a, w);

        // termination criteria -- optimality gap
        mu = dist - (v * w) * invDist;
        if(mu <= dist * relError || mu < absError)
            break;

        // termination criteria -- degenerate case
        if(degenerate(bits, y, w))
            break;

        // if not terminated, get ready for the next iteration
        y[last] = w;
        bits |= (1 << last);
        sv_subalgorithm(y, bits, lambdas, v);

        ++numIterations;
        dist2   = norm2(v);
        invDist = (dist2 > T(0)) ? inverseSqrt(dist2) : T(0);
        dist    = dist2 * invDist;
    }

    computePoints<T>(bits, p, q, lambdas, pa, pb);

    if(numIterations >= 1000)
        catch_me();
    else
        nbIter = numIterations;

    return (dist);
}

// -------------------------------------------------------------------------------------------------
// Dispatcher function for Vector/Quaternion
template <typename T, GJKType GJKType, bool Acceleration>
__HOSTDEVICE__ T computeClosestPoints_GJK(const Convex<T>&     a,
                                          const Convex<T>&     b,
                                          const Vector3<T>&    v_b2a,
                                          const Quaternion<T>& q_b2a,
                                          const T              crustA,
                                          const T              crustB,
                                          Vector3<T>&          pa,
                                          Vector3<T>&          pb,
                                          uint&                nbIter)
{
    static_assert(GJKType == GJKType::JOHNSON || GJKType == GJKType::SIGNEDVOLUME,
                  "GJKType must be either Johnson or SignedVolume");

    if constexpr(GJKType == GJKType::JOHNSON)
    {
        return computeClosestPoints_GJK_Johnson<T, Acceleration>(a,
                                                                 b,
                                                                 v_b2a,
                                                                 q_b2a,
                                                                 crustA,
                                                                 crustB,
                                                                 pa,
                                                                 pb,
                                                                 nbIter);
    }
    else
    {
        return computeClosestPoints_GJK_SignedVolume<T, Acceleration>(a,
                                                                      b,
                                                                      v_b2a,
                                                                      q_b2a,
                                                                      crustA,
                                                                      crustB,
                                                                      pa,
                                                                      pb,
                                                                      nbIter);
    }
}

// -------------------------------------------------------------------------------------------------
// Johnson implementation for Vector/Vector/Quaternion/Quaternion
template <typename T, bool Acceleration>
__HOSTDEVICE__ T computeClosestPoints_GJK_Johnson(const Convex<T>&     a,
                                                  const Convex<T>&     b,
                                                  const Vector3<T>&    v_a2w,
                                                  const Vector3<T>&    v_b2w,
                                                  const Quaternion<T>& q_a2w,
                                                  const Quaternion<T>& q_b2w,
                                                  const T              crustA,
                                                  const T              crustB,
                                                  Vector3<T>&          pa,
                                                  Vector3<T>&          pb,
                                                  uint&                nbIter)
{
    // Constants
    constexpr T relError = EPS<T>;               // relative tolerance
    constexpr T absError = T(1.e-4 * relError);  // absolute tolerance

    // Johnson-specific variables
    uint bits     = 0;  // identifies current simplex
    uint last     = 0;  // identifies last found support point
    uint last_bit = 0;  // last_bit = 1<<last
    uint all_bits = 0;  // all_bits = bits|last_bit

    Vector3<T> p[4];                 // support points of A in local
    Vector3<T> q[4];                 // support points of B in local
    Vector3<T> y[4];                 // support points of A-B in world
    T          det[16][4] = {T(0)};  // cached sub-determinants
    T          dp[4][4]   = {T(0)};  // cached dot products

    T   mu            = T(0);  // optimality gap
    int numIterations = 0;     // No. iterations

    // Acceleration-specific variables
    T momentum = T(0);  // Only used if Acceleration is true

    // Initializing vectors
    Vector3<T> v(v_a2w - v_b2w);
    Vector3<T> w;
    T          dist2   = norm2(v);
    T          invDist = (dist2 > T(0)) ? inverseSqrt(dist2) : T(0);
    T          dist    = dist2 * invDist;

    while(bits < 15 && dist > HIGHEPS<T> && numIterations < 1000)
    {
        // Updating the bits
        last     = 0;
        last_bit = 1;
        while(bits & last_bit)
        {
            ++last;
            last_bit <<= 1;
        }

        // Support points (eroded by crust in local frame, reusing local dirs)
        const Vector3<T> vA = q_a2w << (-v);
        const Vector3<T> vB = q_b2w << v;
        p[last]             = a.support(vA, crustA, invDist);
        q[last]             = b.support(vB, crustB, invDist);
        FusedMinkowskiDifference(p[last], q[last], v_a2w, v_b2w, q_a2w, q_b2w, w);

        // termination criteria -- optimality gap
        mu = dist - (v * w) * invDist;
        if(mu <= dist * relError || mu < absError)
            break;

        // termination criteria -- degenerate case
        if(degenerate(all_bits, y, w))
            break;

        // if not terminated, get ready for the next iteration
        y[last]  = w;
        all_bits = bits | last_bit;
        if(!closest(bits, last, last_bit, all_bits, y, dp, det, v))
            break;

        ++numIterations;
        dist2   = norm2(v);
        invDist = (dist2 > T(0)) ? inverseSqrt(dist2) : T(0);
        dist    = dist2 * invDist;
    }

    computePoints<T>(bits, p, q, det, pa, pb);

    if(numIterations >= 1000)
        catch_me();
    else
        nbIter = numIterations;

    return (dist);
}

// -------------------------------------------------------------------------------------------------
// SignedVolume implementation for Vector/Vector/Quaternion/Quaternion
template <typename T, bool Acceleration>
__HOSTDEVICE__ T computeClosestPoints_GJK_SignedVolume(const Convex<T>&     a,
                                                       const Convex<T>&     b,
                                                       const Vector3<T>&    v_a2w,
                                                       const Vector3<T>&    v_b2w,
                                                       const Quaternion<T>& q_a2w,
                                                       const Quaternion<T>& q_b2w,
                                                       const T              crustA,
                                                       const T              crustB,
                                                       Vector3<T>&          pa,
                                                       Vector3<T>&          pb,
                                                       uint&                nbIter)
{
    // Constants
    constexpr T relError = EPS<T>;               // relative tolerance
    constexpr T absError = T(1.e-4 * relError);  // absolute tolerance

    // SignedVolume-specific variables
    uint bits = 0;  // identifies current simplex
    uint last = 0;  // identifies last found support point

    Vector3<T> p[4];                 // support points of A in local
    Vector3<T> q[4];                 // support points of B in local
    Vector3<T> y[4];                 // support points of A-B in world
    T          lambdas[4] = {T(0)};  // Weights

    T   mu            = T(0);  // optimality gap
    int numIterations = 0;     // No. iterations

    // Acceleration-specific variables
    T momentum = T(0);  // Only used if Acceleration is true

    // Initializing vectors
    Vector3<T> v(v_a2w - v_b2w);
    Vector3<T> w;
    T          dist2   = norm2(v);
    T          invDist = (dist2 > T(0)) ? inverseSqrt(dist2) : T(0);
    T          dist    = dist2 * invDist;

    while(bits < 15 && dist > HIGHEPS<T> && numIterations < 1000)
    {
        // Updating the bits
        for(uint new_index = 0; new_index < 4; ++new_index)
        {
            // At least one of these must be empty, otherwise overlap.
            if(!(bits & (1 << new_index)))
            {
                last = new_index;
                break;
            }
        }

        // Support points (eroded by crust in local frame, reusing local dirs)
        const Vector3<T> vA = q_a2w << (-v);
        const Vector3<T> vB = q_b2w << v;
        p[last]             = a.support(vA, crustA, invDist);
        q[last]             = b.support(vB, crustB, invDist);
        FusedMinkowskiDifference(p[last], q[last], v_a2w, v_b2w, q_a2w, q_b2w, w);

        // termination criteria -- optimality gap
        mu = dist - (v * w) * invDist;
        if(mu <= dist * relError || mu < absError)
            break;

        // termination criteria -- degenerate case
        if(degenerate(bits, y, w))
            break;

        // if not terminated, get ready for the next iteration
        y[last] = w;
        bits |= (1 << last);
        sv_subalgorithm(y, bits, lambdas, v);

        ++numIterations;
        dist2   = norm2(v);
        invDist = (dist2 > T(0)) ? inverseSqrt(dist2) : T(0);
        dist    = dist2 * invDist;
    }

    computePoints<T>(bits, p, q, lambdas, pa, pb);

    if(numIterations >= 1000)
        catch_me();
    else
        nbIter = numIterations;

    return (dist);
}

// -------------------------------------------------------------------------------------------------
// Dispatcher function for Vector/Vector/Quaternion/Quaternion
template <typename T, GJKType GJKType, bool Acceleration>
__HOSTDEVICE__ T computeClosestPoints_GJK(const Convex<T>&     a,
                                          const Convex<T>&     b,
                                          const Vector3<T>&    v_a2w,
                                          const Vector3<T>&    v_b2w,
                                          const Quaternion<T>& q_a2w,
                                          const Quaternion<T>& q_b2w,
                                          const T              crustA,
                                          const T              crustB,
                                          Vector3<T>&          pa,
                                          Vector3<T>&          pb,
                                          uint&                nbIter)
{
    static_assert(GJKType == GJKType::JOHNSON || GJKType == GJKType::SIGNEDVOLUME,
                  "GJKType must be either Johnson or SignedVolume");

    if constexpr(GJKType == GJKType::JOHNSON)
    {
        return computeClosestPoints_GJK_Johnson<T, Acceleration>(a,
                                                                 b,
                                                                 v_a2w,
                                                                 v_b2w,
                                                                 q_a2w,
                                                                 q_b2w,
                                                                 crustA,
                                                                 crustB,
                                                                 pa,
                                                                 pb,
                                                                 nbIter);
    }
    else
    {
        return computeClosestPoints_GJK_SignedVolume<T, Acceleration>(a,
                                                                      b,
                                                                      v_a2w,
                                                                      v_b2w,
                                                                      q_a2w,
                                                                      q_b2w,
                                                                      crustA,
                                                                      crustB,
                                                                      pa,
                                                                      pb,
                                                                      nbIter);
    }
}

// -------------------------------------------------------------------------------------------------
// Johnson implementation for vec/quat
template <typename T, bool Acceleration>
__HOSTDEVICE__ T computeClosestPoints_GJK_Johnson(const ShapeData<T>&  sdA,
                                                  const ShapeData<T>&  sdB,
                                                  const Vector3<T>&    v_b2a,
                                                  const Quaternion<T>& q_b2a,
                                                  const T              crustA,
                                                  const T              crustB,
                                                  Vector3<T>&          pa,
                                                  Vector3<T>&          pb,
                                                  uint&                nbIter)
{
    // Constants
    constexpr T relError = EPS<T>;               // relative tolerance
    constexpr T absError = T(1.e-4 * relError);  // absolute tolerance

    // Johnson-specific variables
    uint bits     = 0;  // identifies current simplex
    uint last     = 0;  // identifies last found support point
    uint last_bit = 0;  // last_bit = 1<<last
    uint all_bits = 0;  // all_bits = bits|last_bit

    Vector3<T> p[4];                 // support points of A in local
    Vector3<T> q[4];                 // support points of B in local
    Vector3<T> y[4];                 // support points of A-B in world
    T          det[16][4] = {T(0)};  // cached sub-determinants
    T          dp[4][4]   = {T(0)};  // cached dot products

    T   mu            = T(0);  // optimality gap
    int numIterations = 0;     // No. iterations

    // Acceleration-specific variables
    T momentum = T(0);  // Only used if Acceleration is true

    // Initializing vectors
    Vector3<T> v(v_b2a);
    Vector3<T> w;
    T          dist2   = norm2(v);
    T          invDist = (dist2 > T(0)) ? inverseSqrt(dist2) : T(0);
    T          dist    = dist2 * invDist;

    while(bits < 15 && dist > HIGHEPS<T> && numIterations < 1000)
    {
        // Updating the bits
        last     = 0;
        last_bit = 1;
        while(bits & last_bit)
        {
            ++last;
            last_bit <<= 1;
        }

        // Support points (vtable-free via ShapeData)
        const Vector3<T> vB = q_b2a << v;
        p[last]             = device_support(sdA, -v, crustA, invDist);
        q[last]             = device_support(sdB, vB, crustB, invDist);
        FusedMinkowskiDifference(p[last], q[last], v_b2a, q_b2a, w);

        // termination criteria -- optimality gap
        mu = dist - (v * w) * invDist;
        if(mu <= dist * relError || mu < absError)
            break;

        // termination criteria -- degenerate case
        if(degenerate(all_bits, y, w))
            break;

        // if not terminated, get ready for the next iteration
        y[last]  = w;
        all_bits = bits | last_bit;
        if(!closest(bits, last, last_bit, all_bits, y, dp, det, v))
            break;

        ++numIterations;
        dist2   = norm2(v);
        invDist = (dist2 > T(0)) ? inverseSqrt(dist2) : T(0);
        dist    = dist2 * invDist;
    }

    computePoints<T>(bits, p, q, det, pa, pb);

    if(numIterations >= 1000)
        catch_me();
    else
        nbIter = numIterations;

    return (dist);
}

// -------------------------------------------------------------------------------------------------
// SignedVolume implementation for vec/quat
template <typename T, bool Acceleration>
__HOSTDEVICE__ T computeClosestPoints_GJK_SignedVolume(const ShapeData<T>&  sdA,
                                                       const ShapeData<T>&  sdB,
                                                       const Vector3<T>&    v_b2a,
                                                       const Quaternion<T>& q_b2a,
                                                       const T              crustA,
                                                       const T              crustB,
                                                       Vector3<T>&          pa,
                                                       Vector3<T>&          pb,
                                                       uint&                nbIter)
{
    // Constants
    constexpr T relError = EPS<T>;               // relative tolerance
    constexpr T absError = T(1.e-4 * relError);  // absolute tolerance

    // SignedVolume-specific variables
    uint bits = 0;  // identifies current simplex
    uint last = 0;  // identifies last found support point

    Vector3<T> p[4];                 // support points of A in local
    Vector3<T> q[4];                 // support points of B in local
    Vector3<T> y[4];                 // support points of A-B in world
    T          lambdas[4] = {T(0)};  // Weights

    T   mu            = T(0);  // optimality gap
    int numIterations = 0;     // No. iterations

    // Acceleration-specific variables
    T momentum = T(0);  // Only used if Acceleration is true

    // Initializing vectors
    Vector3<T> v(v_b2a);
    Vector3<T> w;
    T          dist2   = norm2(v);
    T          invDist = (dist2 > T(0)) ? inverseSqrt(dist2) : T(0);
    T          dist    = dist2 * invDist;

    while(bits < 15 && dist > HIGHEPS<T> && numIterations < 1000)
    {
        // Updating the bits
        for(uint new_index = 0; new_index < 4; ++new_index)
        {
            // At least one of these must be empty, otherwise overlap.
            if(!(bits & (1 << new_index)))
            {
                last = new_index;
                break;
            }
        }

        // Support points (vtable-free via ShapeData)
        const Vector3<T> vB = q_b2a << v;
        p[last]             = device_support(sdA, -v, crustA, invDist);
        q[last]             = device_support(sdB, vB, crustB, invDist);
        FusedMinkowskiDifference(p[last], q[last], v_b2a, q_b2a, w);

        // termination criteria -- optimality gap
        mu = dist - (v * w) * invDist;
        if(mu <= dist * relError || mu < absError)
            break;

        // termination criteria -- degenerate case
        if(degenerate(bits, y, w))
            break;

        // if not terminated, get ready for the next iteration
        y[last] = w;
        bits |= (1 << last);
        sv_subalgorithm(y, bits, lambdas, v);

        ++numIterations;
        dist2   = norm2(v);
        invDist = (dist2 > T(0)) ? inverseSqrt(dist2) : T(0);
        dist    = dist2 * invDist;
    }

    computePoints<T>(bits, p, q, lambdas, pa, pb);

    if(numIterations >= 1000)
        catch_me();
    else
        nbIter = numIterations;

    return (dist);
}

// -------------------------------------------------------------------------------------------------
// Dispatcher function for vec/quat
template <typename T, GJKType GJKType, bool Acceleration>
__HOSTDEVICE__ T computeClosestPoints_GJK(const ShapeData<T>&  sdA,
                                          const ShapeData<T>&  sdB,
                                          const Vector3<T>&    v_b2a,
                                          const Quaternion<T>& q_b2a,
                                          const T              crustA,
                                          const T              crustB,
                                          Vector3<T>&          pa,
                                          Vector3<T>&          pb,
                                          uint&                nbIter)
{
    static_assert(GJKType == GJKType::JOHNSON || GJKType == GJKType::SIGNEDVOLUME,
                  "GJKType must be either Johnson or SignedVolume");

    if constexpr(GJKType == GJKType::JOHNSON)
    {
        return computeClosestPoints_GJK_Johnson<T, Acceleration>(sdA,
                                                                 sdB,
                                                                 v_b2a,
                                                                 q_b2a,
                                                                 crustA,
                                                                 crustB,
                                                                 pa,
                                                                 pb,
                                                                 nbIter);
    }
    else
    {
        return computeClosestPoints_GJK_SignedVolume<T, Acceleration>(sdA,
                                                                      sdB,
                                                                      v_b2a,
                                                                      q_b2a,
                                                                      crustA,
                                                                      crustB,
                                                                      pa,
                                                                      pb,
                                                                      nbIter);
    }
}

// -------------------------------------------------------------------------------------------------
// Johnson sub-algorithm, world frame, vtable-free via ShapeData
template <typename T, bool Acceleration>
__HOSTDEVICE__ T computeClosestPoints_GJK_Johnson(const ShapeData<T>&  sdA,
                                                  const ShapeData<T>&  sdB,
                                                  const Vector3<T>&    v_a2w,
                                                  const Vector3<T>&    v_b2w,
                                                  const Quaternion<T>& q_a2w,
                                                  const Quaternion<T>& q_b2w,
                                                  const T              crustA,
                                                  const T              crustB,
                                                  Vector3<T>&          pa,
                                                  Vector3<T>&          pb,
                                                  uint&                nbIter)
{
    // Constants
    constexpr T relError = EPS<T>;               // relative tolerance
    constexpr T absError = T(1.e-4 * relError);  // absolute tolerance

    // Johnson-specific variables
    uint bits     = 0;  // identifies current simplex
    uint last     = 0;  // identifies last found support point
    uint last_bit = 0;  // last_bit = 1<<last
    uint all_bits = 0;  // all_bits = bits|last_bit

    Vector3<T> p[4];                 // support points of A in local
    Vector3<T> q[4];                 // support points of B in local
    Vector3<T> y[4];                 // support points of A-B in world
    T          det[16][4] = {T(0)};  // cached sub-determinants
    T          dp[4][4]   = {T(0)};  // cached dot products

    T   mu            = T(0);  // optimality gap
    int numIterations = 0;     // No. iterations

    // Acceleration-specific variables
    T momentum = T(0);  // Only used if Acceleration is true

    // Initializing vectors
    Vector3<T> v(v_a2w - v_b2w);
    Vector3<T> w;
    T          dist2   = norm2(v);
    T          invDist = (dist2 > T(0)) ? inverseSqrt(dist2) : T(0);
    T          dist    = dist2 * invDist;

    while(bits < 15 && dist > HIGHEPS<T> && numIterations < 1000)
    {
        // Updating the bits
        last     = 0;
        last_bit = 1;
        while(bits & last_bit)
        {
            ++last;
            last_bit <<= 1;
        }

        // Support points (vtable-free via ShapeData)
        const Vector3<T> vA = q_a2w << (-v);
        const Vector3<T> vB = q_b2w << v;
        p[last]             = device_support(sdA, vA, crustA, invDist);
        q[last]             = device_support(sdB, vB, crustB, invDist);
        FusedMinkowskiDifference(p[last], q[last], v_a2w, v_b2w, q_a2w, q_b2w, w);

        // termination criteria -- optimality gap
        mu = dist - (v * w) * invDist;
        if(mu <= dist * relError || mu < absError)
            break;

        // termination criteria -- degenerate case
        if(degenerate(all_bits, y, w))
            break;

        // if not terminated, get ready for the next iteration
        y[last]  = w;
        all_bits = bits | last_bit;
        if(!closest(bits, last, last_bit, all_bits, y, dp, det, v))
            break;

        ++numIterations;
        dist2   = norm2(v);
        invDist = (dist2 > T(0)) ? inverseSqrt(dist2) : T(0);
        dist    = dist2 * invDist;
    }

    computePoints<T>(bits, p, q, det, pa, pb);

    if(numIterations >= 1000)
        catch_me();
    else
        nbIter = numIterations;

    return (dist);
}

// -------------------------------------------------------------------------------------------------
// SignedVolume sub-algorithm, world frame, vtable-free via ShapeData
template <typename T, bool Acceleration>
__HOSTDEVICE__ T computeClosestPoints_GJK_SignedVolume(const ShapeData<T>&  sdA,
                                                       const ShapeData<T>&  sdB,
                                                       const Vector3<T>&    v_a2w,
                                                       const Vector3<T>&    v_b2w,
                                                       const Quaternion<T>& q_a2w,
                                                       const Quaternion<T>& q_b2w,
                                                       const T              crustA,
                                                       const T              crustB,
                                                       Vector3<T>&          pa,
                                                       Vector3<T>&          pb,
                                                       uint&                nbIter)
{
    // Constants
    constexpr T relError = EPS<T>;               // relative tolerance
    constexpr T absError = T(1.e-4 * relError);  // absolute tolerance

    // SignedVolume-specific variables
    uint bits = 0;  // identifies current simplex
    uint last = 0;  // identifies last found support point

    Vector3<T> p[4];                 // support points of A in local
    Vector3<T> q[4];                 // support points of B in local
    Vector3<T> y[4];                 // support points of A-B in world
    T          lambdas[4] = {T(0)};  // Weights

    T   mu            = T(0);  // optimality gap
    int numIterations = 0;     // No. iterations

    // Acceleration-specific variables
    T momentum = T(0);  // Only used if Acceleration is true

    // Initializing vectors
    Vector3<T> v(v_a2w - v_b2w);
    Vector3<T> w;
    T          dist2   = norm2(v);
    T          invDist = (dist2 > T(0)) ? inverseSqrt(dist2) : T(0);
    T          dist    = dist2 * invDist;

    while(bits < 15 && dist > HIGHEPS<T> && numIterations < 1000)
    {
        // Updating the bits
        for(uint new_index = 0; new_index < 4; ++new_index)
        {
            // At least one of these must be empty, otherwise overlap.
            if(!(bits & (1 << new_index)))
            {
                last = new_index;
                break;
            }
        }

        // Support points (vtable-free via ShapeData)
        const Vector3<T> vA = q_a2w << (-v);
        const Vector3<T> vB = q_b2w << v;
        p[last]             = device_support(sdA, vA, crustA, invDist);
        q[last]             = device_support(sdB, vB, crustB, invDist);
        FusedMinkowskiDifference(p[last], q[last], v_a2w, v_b2w, q_a2w, q_b2w, w);

        // termination criteria -- optimality gap
        mu = dist - (v * w) * invDist;
        if(mu <= dist * relError || mu < absError)
            break;

        // termination criteria -- degenerate case
        if(degenerate(bits, y, w))
            break;

        // if not terminated, get ready for the next iteration
        y[last] = w;
        bits |= (1 << last);
        sv_subalgorithm(y, bits, lambdas, v);

        ++numIterations;
        dist2   = norm2(v);
        invDist = (dist2 > T(0)) ? inverseSqrt(dist2) : T(0);
        dist    = dist2 * invDist;
    }

    computePoints<T>(bits, p, q, lambdas, pa, pb);

    if(numIterations >= 1000)
        catch_me();
    else
        nbIter = numIterations;

    return (dist);
}

// -------------------------------------------------------------------------------------------------
// Dispatcher world frame
template <typename T, GJKType GJKType, bool Acceleration>
__HOSTDEVICE__ T computeClosestPoints_GJK(const ShapeData<T>&  sdA,
                                          const ShapeData<T>&  sdB,
                                          const Vector3<T>&    v_a2w,
                                          const Vector3<T>&    v_b2w,
                                          const Quaternion<T>& q_a2w,
                                          const Quaternion<T>& q_b2w,
                                          const T              crustA,
                                          const T              crustB,
                                          Vector3<T>&          pa,
                                          Vector3<T>&          pb,
                                          uint&                nbIter)
{
    static_assert(GJKType == GJKType::JOHNSON || GJKType == GJKType::SIGNEDVOLUME,
                  "GJKType must be either Johnson or SignedVolume");

    if constexpr(GJKType == GJKType::JOHNSON)
    {
        return computeClosestPoints_GJK_Johnson<T, Acceleration>(sdA,
                                                                 sdB,
                                                                 v_a2w,
                                                                 v_b2w,
                                                                 q_a2w,
                                                                 q_b2w,
                                                                 crustA,
                                                                 crustB,
                                                                 pa,
                                                                 pb,
                                                                 nbIter);
    }
    else
    {
        return computeClosestPoints_GJK_SignedVolume<T, Acceleration>(sdA,
                                                                      sdB,
                                                                      v_a2w,
                                                                      v_b2w,
                                                                      q_a2w,
                                                                      q_b2w,
                                                                      crustA,
                                                                      crustB,
                                                                      pa,
                                                                      pb,
                                                                      nbIter);
    }
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
#define X(T)                                                               \
    template __HOSTDEVICE__ bool intersectGJK(const Convex<T>&     a,      \
                                              const Convex<T>&     b,      \
                                              const Transform3<T>& b2a);   \
    template __HOSTDEVICE__ bool intersectGJK(const Convex<T>&     a,      \
                                              const Convex<T>&     b,      \
                                              const Transform3<T>& a2w,    \
                                              const Transform3<T>& b2w);   \
    template __HOSTDEVICE__ bool intersectGJK(const Convex<T>&     a,      \
                                              const Convex<T>&     b,      \
                                              const Vector3<T>&    v_b2a,  \
                                              const Quaternion<T>& q_b2a); \
    template __HOSTDEVICE__ bool intersectGJK(const Convex<T>&     a,      \
                                              const Convex<T>&     b,      \
                                              const Vector3<T>&    v_a2w,  \
                                              const Vector3<T>&    v_b2w,  \
                                              const Quaternion<T>& q_a2w,  \
                                              const Quaternion<T>& q_b2w);
X(float)
X(double)
#undef X

#define X(T, GJK, ACC)                                                                           \
    template __HOSTDEVICE__ T computeClosestPoints_GJK<T, GJK, ACC>(const Convex<T>&     a,      \
                                                                    const Convex<T>&     b,      \
                                                                    const Transform3<T>& b2a,    \
                                                                    const T              crustA, \
                                                                    const T              crustB, \
                                                                    Vector3<T>&          pa,     \
                                                                    Vector3<T>&          pb,     \
                                                                    uint&                nbIter);               \
    template __HOSTDEVICE__ T computeClosestPoints_GJK<T, GJK, ACC>(const Convex<T>&     a,      \
                                                                    const Convex<T>&     b,      \
                                                                    const Transform3<T>& a2w,    \
                                                                    const Transform3<T>& b2w,    \
                                                                    const T              crustA, \
                                                                    const T              crustB, \
                                                                    Vector3<T>&          pa,     \
                                                                    Vector3<T>&          pb,     \
                                                                    uint&                nbIter);               \
    template __HOSTDEVICE__ T computeClosestPoints_GJK<T, GJK, ACC>(const Convex<T>&     a,      \
                                                                    const Convex<T>&     b,      \
                                                                    const Vector3<T>&    v_b2a,  \
                                                                    const Quaternion<T>& q_b2a,  \
                                                                    const T              crustA, \
                                                                    const T              crustB, \
                                                                    Vector3<T>&          pa,     \
                                                                    Vector3<T>&          pb,     \
                                                                    uint&                nbIter);               \
    template __HOSTDEVICE__ T computeClosestPoints_GJK<T, GJK, ACC>(const Convex<T>&     a,      \
                                                                    const Convex<T>&     b,      \
                                                                    const Vector3<T>&    v_a2w,  \
                                                                    const Vector3<T>&    v_b2w,  \
                                                                    const Quaternion<T>& q_a2w,  \
                                                                    const Quaternion<T>& q_b2w,  \
                                                                    const T              crustA, \
                                                                    const T              crustB, \
                                                                    Vector3<T>&          pa,     \
                                                                    Vector3<T>&          pb,     \
                                                                    uint&                nbIter);               \
    template __HOSTDEVICE__ T computeClosestPoints_GJK<T, GJK, ACC>(const ShapeData<T>&  sdA,    \
                                                                    const ShapeData<T>&  sdB,    \
                                                                    const Vector3<T>&    v_b2a,  \
                                                                    const Quaternion<T>& q_b2a,  \
                                                                    const T              crustA, \
                                                                    const T              crustB, \
                                                                    Vector3<T>&          pa,     \
                                                                    Vector3<T>&          pb,     \
                                                                    uint&                nbIter);               \
    template __HOSTDEVICE__ T computeClosestPoints_GJK<T, GJK, ACC>(const ShapeData<T>&  sdA,    \
                                                                    const ShapeData<T>&  sdB,    \
                                                                    const Vector3<T>&    v_a2w,  \
                                                                    const Vector3<T>&    v_b2w,  \
                                                                    const Quaternion<T>& q_a2w,  \
                                                                    const Quaternion<T>& q_b2w,  \
                                                                    const T              crustA, \
                                                                    const T              crustB, \
                                                                    Vector3<T>&          pa,     \
                                                                    Vector3<T>&          pb,     \
                                                                    uint&                nbIter);
X(float, GJKType::JOHNSON, false)
X(double, GJKType::JOHNSON, false)
X(float, GJKType::JOHNSON, true)
X(double, GJKType::JOHNSON, true)
X(float, GJKType::SIGNEDVOLUME, false)
X(double, GJKType::SIGNEDVOLUME, false)
X(float, GJKType::SIGNEDVOLUME, true)
X(double, GJKType::SIGNEDVOLUME, true)
#undef X