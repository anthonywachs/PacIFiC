#ifndef _OBC_HH_
#define _OBC_HH_

#include "MatrixMath.hh"
#include "Quaternion.hh"
#include "Transform3.hh"

// =================================================================================================
/** @brief The header for the oriented bounding cylinder (OBC) collision detection.

    Oriented Bounding Cylinders (OBC) routines to find whether bounding cylinders are in contact
    or not.

    @author A.Yazdani - 2026 - Construction */
// =================================================================================================
/** @name OBC: Internal helpers */
//@{
// -------------------------------------------------------------------------------------------------
/** @brief Returns the sign of a value: -1, 0, or +1
    @param val The input value */
template <typename T>
__HOSTDEVICE__ static INLINE int obc_sgn(T val)
{
    return ((T(0) < val) - (val < T(0)));
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns the real roots of the quartic x^4 + b*x^3 + c*x^2 + d*x + e = 0.
    @param b  Coefficient of x^3
    @param c  Coefficient of x^2
    @param d  Coefficient of x
    @param e  Constant term
    @param sol Output array of up to 4 roots (caller must provide space for 4 elements)
    @param nbRoots Number of real roots found */
template <typename T>
__HOSTDEVICE__ static INLINE void solveQuarticOBC(T b, T c, T d, T e, T sol[4], int& nbRoots)
{
    nbRoots = 0;

    // Depressed quartic: y^4 + p*y^2 + q*y + r = 0
    const T b2 = b * b;
    const T p  = c - T(3) * b2 / T(8);
    const T q  = b2 * b / T(8) - b * c / T(2) + d;
    const T r  = -T(3) * b2 * b2 / T(256) + e - b * d / T(4) + b2 * c / T(16);
    const T p2 = p * p;

    if(fabs(q) < EPS<T>)
    {
        // Reduce to biquadratic: solve y^2 = (-p +/- sqrt(p^2/4 - r))
        const T del = p2 / T(4) - r;
        if(del < T(0))
            return;

        const T m1 = -p / T(2) + sqrt(del);
        const T m2 = -p / T(2) - sqrt(del);
        if(m1 > T(0))
        {
            sol[nbRoots++] = sqrt(m1) - b / T(4);
            sol[nbRoots++] = -sqrt(m1) - b / T(4);
        }
        if(m2 > T(0))
        {
            sol[nbRoots++] = sqrt(m2) - b / T(4);
            sol[nbRoots++] = -sqrt(m2) - b / T(4);
        }
    }
    else
    {
        // Resolvent cubic: find a real root m of
        // x^3 + (p/2)*x^2 + (p^2/16 - r/4)*x - q^2/64 ... (recast form used below)
        const T u   = -p2 / T(36) - r / T(3);
        const T v   = -p2 * p / T(216) + r * p / T(6) - q * q / T(16);
        const T del = u * u * u + v * v;

        T m = T(0);
        if(del < T(0))
            m = T(2) * sqrt(-u) * cos(acos(v / sqrt(-u) / u) / T(3)) - p / T(3);
        else
        {
            m = cbrt(-v + sqrt(del));
            m = m - u / m - p / T(3);
        }

        if(m < T(0))
            return;

        const T sqrt_mhalf = sqrt(m / T(2));
        const T first_var  = -p / T(2) - m / T(2) - q / sqrt_mhalf / T(4);
        const T second_var = first_var + q / sqrt_mhalf / T(2);

        if(first_var > T(0))
        {
            sol[nbRoots++] = sqrt_mhalf + sqrt(first_var) - b / T(4);
            sol[nbRoots++] = sqrt_mhalf - sqrt(first_var) - b / T(4);
        }
        if(second_var > T(0))
        {
            sol[nbRoots++] = -sqrt_mhalf + sqrt(second_var) - b / T(4);
            sol[nbRoots++] = -sqrt_mhalf - sqrt(second_var) - b / T(4);
        }
    }
}
//@}

/** @name OBC: External methods */
//@{
// -------------------------------------------------------------------------------------------------
/** @brief Returns whether two oriented bounding cylinders are in contact - absolute transformation.
    @param r1   Radius of the first cylinder
    @param h1   Half-height of the first cylinder
    @param ori1 Axis direction of the first cylinder in its local frame
    @param r2   Radius of the second cylinder
    @param h2   Half-height of the second cylinder
    @param ori2 Axis direction of the second cylinder in its local frame
    @param trA2W Transformation from A's local space to world space
    @param trB2W Transformation from B's local space to world space */
template <typename T>
__HOSTDEVICE__ bool intersectOrientedBoundingCylinder(T                    r1,
                                                      T                    h1,
                                                      const Vector3<T>&    ori1,
                                                      T                    r2,
                                                      T                    h2,
                                                      const Vector3<T>&    ori2,
                                                      const Transform3<T>& trA2W,
                                                      const Transform3<T>& trB2W)
{
    // Cylinder axis of A in world space
    const Vector3<T> eZ = trA2W.getBasis() * ori1;
    Vector3<T>       e2 = trB2W.getBasis() * ori2;

    // Vector from center of A to center of B
    const Vector3<T> xPt = trB2W.getOrigin() - trA2W.getOrigin();

    // Local frame for A: eX is perpendicular to both axes, eY completes the frame
    const Vector3<T> eX = (e2 ^ eZ).normalized();
    const Vector3<T> eY = (eZ ^ eX).normalized();

    const T x  = eX * xPt;
    const T y  = eY * xPt;
    const T z  = eZ * xPt;
    const T ey = eY * e2;
    const T ez = eZ * e2;

    /* Step 1: Shortest distance -- Projection onto XZ */
    if(fabs(x) >= r1 + r2)
        return false;
    else
    {
        const T s2 = y / ey;
        const T s1 = z + s2 * ez;
        if(fabs(s1) < h1 && fabs(s2) < h2)
            return true;
    }

    /* Step 2: Projection onto YZ -- check rectangle overlap */
    {
        const T fey = fabs(ey) + LOWEPS<T>;
        const T fez = fabs(ez) + LOWEPS<T>;
        // Projection onto y-axis
        if(fabs(y) > r1 + h2 * fey + r2 * fez)
            return false;
        // Projection onto z-axis
        if(fabs(z) > h1 + h2 * fez + r2 * fey)
            return false;
        // Projection onto e_normal
        if(fabs(y * ez - z * ey) > r2 + h1 * fey + r1 * fez)
            return false;
        // Projection onto e
        if(fabs(y * ey + z * ez) > h2 + h1 * fez + r1 * fey)
            return false;
    }

    /* Step 3: Projection onto XY -- check ellipse/circle overlap */
    if(fabs(ez) < T(0.99) && fabs(ez) > T(0.01))
    {
        // Secondary cylinder projects as an ellipse onto A's XY plane
        {
            const T topy    = y + h2 * ey;
            const T bottomy = y - h2 * ey;
            if((topy >= T(0)) == (bottomy >= T(0)))
            {
                const T x0 = x;
                const T y0 = fabs(topy) > fabs(bottomy) ? bottomy : topy;
                if(x0 * x0 + y0 * y0 / ez / ez > r2 * r2)
                {
                    const T ey_sq = ey * ey;
                    const T A     = r2 * r2 * ey_sq * ey_sq;
                    const T B     = T(2) * y0 * ez / r2 / ey_sq;
                    const T C     = x0 * x0 / A;
                    const T D     = y0 * y0 * ez * ez / A;

                    T   sint[4];
                    int nbRoots;
                    solveQuarticOBC(-B, C + D - T(1), B, -D, sint, nbRoots);

                    for(int i = 0; i < nbRoots; i++)
                    {
                        if(((sint[i] >= T(0)) != (B >= T(0))))
                        {
                            const T cost = -obc_sgn(x0) * sqrt(T(1) - sint[i] * sint[i]);
                            const T ptX  = x0 + r2 * cost;
                            const T ptY  = y0 + r2 * ez * sint[i];
                            if(ptX * ptX + ptY * ptY > r1 * r1)
                                return false;
                        }
                    }
                }
            }
        }

        // Primary cylinder projects as an ellipse onto B's XY plane
        {
            const Vector3<T> eY2     = (eX ^ e2).normalized();
            const T          y2      = -(eY2 * xPt);
            const T          ey2     = eY2 * eZ;
            const T          topy    = y2 + h1 * ey2;
            const T          bottomy = y2 - h1 * ey2;
            if((topy >= T(0)) == (bottomy >= T(0)))
            {
                const T x0 = -x;
                const T y0 = fabs(topy) > fabs(bottomy) ? bottomy : topy;
                if(x0 * x0 + y0 * y0 / ez / ez > r1 * r1)
                {
                    const T ey_sq = ey2 * ey2;
                    const T A     = r1 * r1 * ey_sq * ey_sq;
                    const T B     = T(2) * y0 * ez / r1 / ey_sq;
                    const T C     = x0 * x0 / A;
                    const T D     = y0 * y0 * ez * ez / A;

                    T   sint[4];
                    int nbRoots;
                    solveQuarticOBC(-B, C + D - T(1), B, -D, sint, nbRoots);

                    for(int i = 0; i < nbRoots; i++)
                    {
                        if((sint[i] >= T(0)) != (B >= T(0)))
                        {
                            const T cost = -obc_sgn(x0) * sqrt(T(1) - sint[i] * sint[i]);
                            const T ptX  = x0 + r1 * cost;
                            const T ptY  = y0 + r1 * ez * sint[i];
                            if(ptX * ptX + ptY * ptY > r2 * r2)
                                return false;
                        }
                    }
                }
            }
        }
    }

    return true;
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns whether two oriented bounding cylinders are in contact - relative
    transformation.
    @param r1    Radius of the first cylinder
    @param h1    Half-height of the first cylinder
    @param ori1  Axis direction of the first cylinder in its local frame
    @param r2    Radius of the second cylinder
    @param h2    Half-height of the second cylinder
    @param ori2  Axis direction of the second cylinder in its local frame
    @param trB2A Transformation from B's local space to A's local space */
template <typename T>
__HOSTDEVICE__ bool intersectOrientedBoundingCylinder(T                    r1,
                                                      T                    h1,
                                                      const Vector3<T>&    ori1,
                                                      T                    r2,
                                                      T                    h2,
                                                      const Vector3<T>&    ori2,
                                                      const Transform3<T>& trB2A)
{
    // A's axis is already expressed in A's own frame
    const Vector3<T>& eZ = ori1;
    // B's axis and center-to-center vector expressed in A's frame
    const Vector3<T> e2  = trB2A.getBasis() * ori2;
    const Vector3<T> xPt = trB2A.getOrigin();

    // Local frame for A: eX is perpendicular to both axes, eY completes the frame
    const Vector3<T> eX = (e2 ^ eZ).normalized();
    const Vector3<T> eY = (eZ ^ eX).normalized();

    const T x  = eX * xPt;
    const T y  = eY * xPt;
    const T z  = eZ * xPt;
    const T ey = eY * e2;
    const T ez = eZ * e2;

    /* Step 1: Shortest distance -- Projection onto XZ */
    if(fabs(x) >= r1 + r2)
        return false;
    else
    {
        const T s2 = y / ey;
        const T s1 = z + s2 * ez;
        if(fabs(s1) < h1 && fabs(s2) < h2)
            return true;
    }

    /* Step 2: Projection onto YZ -- check rectangle overlap */
    {
        const T fey = fabs(ey) + LOWEPS<T>;
        const T fez = fabs(ez) + LOWEPS<T>;
        // Projection onto y-axis
        if(fabs(y) > r1 + h2 * fey + r2 * fez)
            return false;
        // Projection onto z-axis
        if(fabs(z) > h1 + h2 * fez + r2 * fey)
            return false;
        // Projection onto e_normal
        if(fabs(y * ez - z * ey) > r2 + h1 * fey + r1 * fez)
            return false;
        // Projection onto e
        if(fabs(y * ey + z * ez) > h2 + h1 * fez + r1 * fey)
            return false;
    }

    /* Step 3: Projection onto XY -- check ellipse/circle overlap */
    if(fabs(ez) < T(0.99) && fabs(ez) > T(0.01))
    {
        // Secondary cylinder projects as an ellipse onto A's XY plane
        {
            const T topy    = y + h2 * ey;
            const T bottomy = y - h2 * ey;
            if((topy >= T(0)) == (bottomy >= T(0)))
            {
                const T x0 = x;
                const T y0 = fabs(topy) > fabs(bottomy) ? bottomy : topy;
                if(x0 * x0 + y0 * y0 / ez / ez > r2 * r2)
                {
                    const T ey_sq = ey * ey;
                    const T A     = r2 * r2 * ey_sq * ey_sq;
                    const T B     = T(2) * y0 * ez / r2 / ey_sq;
                    const T C     = x0 * x0 / A;
                    const T D     = y0 * y0 * ez * ez / A;

                    T   sint[4];
                    int nbRoots;
                    solveQuarticOBC(-B, C + D - T(1), B, -D, sint, nbRoots);

                    for(int i = 0; i < nbRoots; i++)
                    {
                        if((sint[i] >= T(0)) != (B >= T(0)))
                        {
                            const T cost = -obc_sgn(x0) * sqrt(T(1) - sint[i] * sint[i]);
                            const T ptX  = x0 + r2 * cost;
                            const T ptY  = y0 + r2 * ez * sint[i];
                            if(ptX * ptX + ptY * ptY > r1 * r1)
                                return false;
                        }
                    }
                }
            }
        }

        // Primary cylinder projects as an ellipse onto B's XY plane
        {
            const Vector3<T> eY2     = (eX ^ e2).normalized();
            const T          y2      = -(eY2 * xPt);
            const T          ey2     = eY2 * eZ;
            const T          topy    = y2 + h1 * ey2;
            const T          bottomy = y2 - h1 * ey2;
            if((topy >= T(0)) == (bottomy >= T(0)))
            {
                const T x0 = -x;
                const T y0 = fabs(topy) > fabs(bottomy) ? bottomy : topy;
                if(x0 * x0 + y0 * y0 / ez / ez > r1 * r1)
                {
                    const T ey_sq = ey2 * ey2;
                    const T A     = r1 * r1 * ey_sq * ey_sq;
                    const T B     = T(2) * y0 * ez / r1 / ey_sq;
                    const T C     = x0 * x0 / A;
                    const T D     = y0 * y0 * ez * ez / A;

                    T   sint[4];
                    int nbRoots;
                    solveQuarticOBC(-B, C + D - T(1), B, -D, sint, nbRoots);

                    for(int i = 0; i < nbRoots; i++)
                    {
                        if((sint[i] >= T(0)) != (B >= T(0)))
                        {
                            const T cost = -obc_sgn(x0) * sqrt(T(1) - sint[i] * sint[i]);
                            const T ptX  = x0 + r1 * cost;
                            const T ptY  = y0 + r1 * ez * sint[i];
                            if(ptX * ptX + ptY * ptY > r2 * r2)
                                return false;
                        }
                    }
                }
            }
        }
    }

    return true;
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns whether two oriented bounding cylinders are in contact - quaternion version.
    @param r1    Radius of the first cylinder
    @param h1    Half-height of the first cylinder
    @param ori1  Axis direction of the first cylinder in its local frame
    @param r2    Radius of the second cylinder
    @param h2    Half-height of the second cylinder
    @param ori2  Axis direction of the second cylinder in its local frame
    @param v_a2w Translation from A's local space to world space
    @param v_b2w Translation from B's local space to world space
    @param q_a2w Rotation from A's local space to world space
    @param q_b2w Rotation from B's local space to world space */
template <typename T>
__HOSTDEVICE__ bool intersectOrientedBoundingCylinder(T                    r1,
                                                      T                    h1,
                                                      const Vector3<T>&    ori1,
                                                      T                    r2,
                                                      T                    h2,
                                                      const Vector3<T>&    ori2,
                                                      const Vector3<T>&    v_a2w,
                                                      const Vector3<T>&    v_b2w,
                                                      const Quaternion<T>& q_a2w,
                                                      const Quaternion<T>& q_b2w)
{
    // Cylinder axes in world space (direct quaternion rotation, no matrix construction)
    const Vector3<T> eZ = q_a2w >> ori1;
    const Vector3<T> e2 = q_b2w >> ori2;

    // Vector from center of A to center of B
    const Vector3<T> xPt = v_b2w - v_a2w;

    // Local frame for A
    const Vector3<T> eX = (e2 ^ eZ).normalized();
    const Vector3<T> eY = (eZ ^ eX).normalized();

    const T x  = eX * xPt;
    const T y  = eY * xPt;
    const T z  = eZ * xPt;
    const T ey = eY * e2;
    const T ez = eZ * e2;

    /* Step 1: Shortest distance -- Projection onto XZ */
    if(fabs(x) >= r1 + r2)
        return false;
    else
    {
        const T s2 = y / ey;
        const T s1 = z + s2 * ez;
        if(fabs(s1) < h1 && fabs(s2) < h2)
            return true;
    }

    /* Step 2: Projection onto YZ -- check rectangle overlap */
    {
        const T fey = fabs(ey) + LOWEPS<T>;
        const T fez = fabs(ez) + LOWEPS<T>;
        if(fabs(y) > r1 + h2 * fey + r2 * fez)
            return false;
        if(fabs(z) > h1 + h2 * fez + r2 * fey)
            return false;
        if(fabs(y * ez - z * ey) > r2 + h1 * fey + r1 * fez)
            return false;
        if(fabs(y * ey + z * ez) > h2 + h1 * fez + r1 * fey)
            return false;
    }

    /* Step 3: Projection onto XY -- check ellipse/circle overlap */
    if(fabs(ez) < T(0.99) && fabs(ez) > T(0.01))
    {
        // Secondary cylinder projects as an ellipse
        {
            const T topy    = y + h2 * ey;
            const T bottomy = y - h2 * ey;
            if((topy >= T(0)) == (bottomy >= T(0)))
            {
                const T x0 = x;
                const T y0 = fabs(topy) > fabs(bottomy) ? bottomy : topy;
                if(x0 * x0 + y0 * y0 / ez / ez > r2 * r2)
                {
                    const T ey_sq = ey * ey;
                    const T A     = r2 * r2 * ey_sq * ey_sq;
                    const T B     = T(2) * y0 * ez / r2 / ey_sq;
                    const T C     = x0 * x0 / A;
                    const T D     = y0 * y0 * ez * ez / A;

                    T   sint[4];
                    int nbRoots;
                    solveQuarticOBC(-B, C + D - T(1), B, -D, sint, nbRoots);

                    for(int i = 0; i < nbRoots; i++)
                    {
                        if((sint[i] >= T(0)) != (B >= T(0)))
                        {
                            const T cost = -obc_sgn(x0) * sqrt(T(1) - sint[i] * sint[i]);
                            const T ptX  = x0 + r2 * cost;
                            const T ptY  = y0 + r2 * ez * sint[i];
                            if(ptX * ptX + ptY * ptY > r1 * r1)
                                return false;
                        }
                    }
                }
            }
        }

        // Primary cylinder projects as an ellipse
        {
            const Vector3<T> eY2     = (eX ^ e2).normalized();
            const T          y2      = -(eY2 * xPt);
            const T          ey2     = eY2 * eZ;
            const T          topy    = y2 + h1 * ey2;
            const T          bottomy = y2 - h1 * ey2;
            if((topy >= T(0)) == (bottomy >= T(0)))
            {
                const T x0 = -x;
                const T y0 = fabs(topy) > fabs(bottomy) ? bottomy : topy;
                if(x0 * x0 + y0 * y0 / ez / ez > r1 * r1)
                {
                    const T ey_sq = ey2 * ey2;
                    const T A     = r1 * r1 * ey_sq * ey_sq;
                    const T B     = T(2) * y0 * ez / r1 / ey_sq;
                    const T C     = x0 * x0 / A;
                    const T D     = y0 * y0 * ez * ez / A;

                    T   sint[4];
                    int nbRoots;
                    solveQuarticOBC(-B, C + D - T(1), B, -D, sint, nbRoots);

                    for(int i = 0; i < nbRoots; i++)
                    {
                        if((sint[i] >= T(0)) != (B >= T(0)))
                        {
                            const T cost = -obc_sgn(x0) * sqrt(T(1) - sint[i] * sint[i]);
                            const T ptX  = x0 + r1 * cost;
                            const T ptY  = y0 + r1 * ez * sint[i];
                            if(ptX * ptX + ptY * ptY > r2 * r2)
                                return false;
                        }
                    }
                }
            }
        }
    }

    return true;
}

// -------------------------------------------------------------------------------------------------
/** @brief Returns whether two oriented bounding cylinders are in contact - quaternion relative
    version.
    @param r1    Radius of the first cylinder
    @param h1    Half-height of the first cylinder
    @param ori1  Axis direction of the first cylinder in its local frame
    @param r2    Radius of the second cylinder
    @param h2    Half-height of the second cylinder
    @param ori2  Axis direction of the second cylinder in its local frame
    @param v_b2a Translation from B's local space to A's local space
    @param q_b2a Rotation from B's local space to A's local space */
template <typename T>
__HOSTDEVICE__ bool intersectOrientedBoundingCylinder(T                    r1,
                                                      T                    h1,
                                                      const Vector3<T>&    ori1,
                                                      T                    r2,
                                                      T                    h2,
                                                      const Vector3<T>&    ori2,
                                                      const Vector3<T>&    v_b2a,
                                                      const Quaternion<T>& q_b2a)
{
    // A's axis is already expressed in A's own frame
    const Vector3<T>& eZ = ori1;
    // B's axis in A's frame (direct quaternion rotation, no matrix construction)
    const Vector3<T> e2  = q_b2a >> ori2;
    const Vector3<T> xPt = v_b2a;

    // Local frame for A: eX is perpendicular to both axes, eY completes the frame
    const Vector3<T> eX = (e2 ^ eZ).normalized();
    const Vector3<T> eY = (eZ ^ eX).normalized();

    const T x  = eX * xPt;
    const T y  = eY * xPt;
    const T z  = eZ * xPt;
    const T ey = eY * e2;
    const T ez = eZ * e2;

    /* Step 1: Shortest distance -- Projection onto XZ */
    if(fabs(x) >= r1 + r2)
        return false;
    else
    {
        const T s2 = y / ey;
        const T s1 = z + s2 * ez;
        if(fabs(s1) < h1 && fabs(s2) < h2)
            return true;
    }

    /* Step 2: Projection onto YZ -- check rectangle overlap */
    {
        const T fey = fabs(ey) + LOWEPS<T>;
        const T fez = fabs(ez) + LOWEPS<T>;
        // Projection onto y-axis
        if(fabs(y) > r1 + h2 * fey + r2 * fez)
            return false;
        // Projection onto z-axis
        if(fabs(z) > h1 + h2 * fez + r2 * fey)
            return false;
        // Projection onto e_normal
        if(fabs(y * ez - z * ey) > r2 + h1 * fey + r1 * fez)
            return false;
        // Projection onto e
        if(fabs(y * ey + z * ez) > h2 + h1 * fez + r1 * fey)
            return false;
    }

    /* Step 3: Projection onto XY -- check ellipse/circle overlap */
    if(fabs(ez) < T(0.99) && fabs(ez) > T(0.01))
    {
        // Secondary cylinder projects as an ellipse onto A's XY plane
        {
            const T topy    = y + h2 * ey;
            const T bottomy = y - h2 * ey;
            if((topy >= T(0)) == (bottomy >= T(0)))
            {
                const T x0 = x;
                const T y0 = fabs(topy) > fabs(bottomy) ? bottomy : topy;
                if(x0 * x0 + y0 * y0 / ez / ez > r2 * r2)
                {
                    const T ey_sq = ey * ey;
                    const T A     = r2 * r2 * ey_sq * ey_sq;
                    const T B     = T(2) * y0 * ez / r2 / ey_sq;
                    const T C     = x0 * x0 / A;
                    const T D     = y0 * y0 * ez * ez / A;

                    T   sint[4];
                    int nbRoots;
                    solveQuarticOBC(-B, C + D - T(1), B, -D, sint, nbRoots);

                    for(int i = 0; i < nbRoots; i++)
                    {
                        if((sint[i] >= T(0)) != (B >= T(0)))
                        {
                            const T cost = -obc_sgn(x0) * sqrt(T(1) - sint[i] * sint[i]);
                            const T ptX  = x0 + r2 * cost;
                            const T ptY  = y0 + r2 * ez * sint[i];
                            if(ptX * ptX + ptY * ptY > r1 * r1)
                                return false;
                        }
                    }
                }
            }
        }

        // Primary cylinder projects as an ellipse onto B's XY plane
        {
            const Vector3<T> eY2     = (eX ^ e2).normalized();
            const T          y2      = -(eY2 * xPt);
            const T          ey2     = eY2 * eZ;
            const T          topy    = y2 + h1 * ey2;
            const T          bottomy = y2 - h1 * ey2;
            if((topy >= T(0)) == (bottomy >= T(0)))
            {
                const T x0 = -x;
                const T y0 = fabs(topy) > fabs(bottomy) ? bottomy : topy;
                if(x0 * x0 + y0 * y0 / ez / ez > r1 * r1)
                {
                    const T ey_sq = ey2 * ey2;
                    const T A     = r1 * r1 * ey_sq * ey_sq;
                    const T B     = T(2) * y0 * ez / r1 / ey_sq;
                    const T C     = x0 * x0 / A;
                    const T D     = y0 * y0 * ez * ez / A;

                    T   sint[4];
                    int nbRoots;
                    solveQuarticOBC(-B, C + D - T(1), B, -D, sint, nbRoots);

                    for(int i = 0; i < nbRoots; i++)
                    {
                        if((sint[i] >= T(0)) != (B >= T(0)))
                        {
                            const T cost = -obc_sgn(x0) * sqrt(T(1) - sint[i] * sint[i]);
                            const T ptX  = x0 + r1 * cost;
                            const T ptY  = y0 + r1 * ez * sint[i];
                            if(ptX * ptX + ptY * ptY > r2 * r2)
                                return false;
                        }
                    }
                }
            }
        }
    }

    return true;
}
//@}

#endif
