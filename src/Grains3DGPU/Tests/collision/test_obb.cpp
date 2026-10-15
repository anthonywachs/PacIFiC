#include "BoundingBox.hh"
#include "OBB.hh"
#include "Quaternion.hh"
#include "QuaternionMath.hh"
#include "Transform3.hh"
#include "Vector3.hh"
#include <cmath>
#include <gtest/gtest.h>
#include <random>

class OBBSimpleTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // Create bounding boxes for testing using extents (half-lengths)
        // Box A: centered at origin, dimensions 2x2x2 (extent 1,1,1)
        Vector3<double> extentA(1.0, 1.0, 1.0);
        boundingBoxA = BoundingBox<double>(extentA);

        // Box B: smaller box, dimensions 1x1x1 (extent 0.5,0.5,0.5)
        Vector3<double> extentB(0.5, 0.5, 0.5);
        boundingBoxB = BoundingBox<double>(extentB);

        // Identity transform
        identity_transform = Transform3<double>(Quaternion<double>(0.0, 0.0, 0.0, 1.0),
                                                Vector3<double>(0.0, 0.0, 0.0));
    }

    BoundingBox<double> boundingBoxA;
    BoundingBox<double> boundingBoxB;
    Transform3<double>  identity_transform;
    const double        EPSILON = 1e-10;
};

// Test basic bounding box creation and properties
TEST_F(OBBSimpleTest, BasicBoundingBoxProperties)
{
    // Test extent values
    Vector3<double> extentA = boundingBoxA.getExtent();
    EXPECT_NEAR(extentA[0], 1.0, EPSILON);
    EXPECT_NEAR(extentA[1], 1.0, EPSILON);
    EXPECT_NEAR(extentA[2], 1.0, EPSILON);

    Vector3<double> extentB = boundingBoxB.getExtent();
    EXPECT_NEAR(extentB[0], 0.5, EPSILON);
    EXPECT_NEAR(extentB[1], 0.5, EPSILON);
    EXPECT_NEAR(extentB[2], 0.5, EPSILON);
}

// Test transform creation and properties
TEST_F(OBBSimpleTest, BasicTransformProperties)
{
    // Test identity transform
    Vector3<double> origin = identity_transform.getOrigin();
    EXPECT_NEAR(origin[0], 0.0, EPSILON);
    EXPECT_NEAR(origin[1], 0.0, EPSILON);
    EXPECT_NEAR(origin[2], 0.0, EPSILON);
}

// Test various bounding box sizes
TEST_F(OBBSimpleTest, VariousBoundingBoxSizes)
{
    // Very small box (extent 0.001)
    Vector3<double>     tinyExtent(0.001, 0.001, 0.001);
    BoundingBox<double> tinyBox(tinyExtent);

    Vector3<double> tinyResult = tinyBox.getExtent();
    EXPECT_NEAR(tinyResult[0], 0.001, EPSILON);
    EXPECT_NEAR(tinyResult[1], 0.001, EPSILON);
    EXPECT_NEAR(tinyResult[2], 0.001, EPSILON);

    // Very large box (extent 1000.0)
    Vector3<double>     hugeExtent(1000.0, 1000.0, 1000.0);
    BoundingBox<double> hugeBox(hugeExtent);

    Vector3<double> hugeResult = hugeBox.getExtent();
    EXPECT_NEAR(hugeResult[0], 1000.0, EPSILON);
    EXPECT_NEAR(hugeResult[1], 1000.0, EPSILON);
    EXPECT_NEAR(hugeResult[2], 1000.0, EPSILON);
}

// Test different constructors
TEST_F(OBBSimpleTest, DifferentConstructors)
{
    // Constructor with individual components
    BoundingBox<double> box1(2.0, 3.0, 4.0);
    Vector3<double>     extent1 = box1.getExtent();
    EXPECT_NEAR(extent1[0], 2.0, EPSILON);
    EXPECT_NEAR(extent1[1], 3.0, EPSILON);
    EXPECT_NEAR(extent1[2], 4.0, EPSILON);

    // Constructor with vector
    Vector3<double>     extent_vector(5.0, 6.0, 7.0);
    BoundingBox<double> box2(extent_vector);
    Vector3<double>     extent2 = box2.getExtent();
    EXPECT_NEAR(extent2[0], 5.0, EPSILON);
    EXPECT_NEAR(extent2[1], 6.0, EPSILON);
    EXPECT_NEAR(extent2[2], 7.0, EPSILON);

    // Default constructor
    BoundingBox<double> box3;
    Vector3<double>     extent3 = box3.getExtent();
    // Default extents should be defined values
    EXPECT_TRUE(extent3[0] >= 0.0);
    EXPECT_TRUE(extent3[1] >= 0.0);
    EXPECT_TRUE(extent3[2] >= 0.0);
}

// =================================================================================================
// OBB intersection tests
//
// intersectOrientedBoundingBox has 4 overloads:
//   1. Absolute Transform3:  (a, b, trA2W, trB2W)
//   2. Relative Transform3:  (a, b, trB2A)
//   3. Absolute quat/vec:    (a, b, v_a2w, v_b2w, q_a2w, q_b2w)
//   4. Relative quat/vec:    (a, b, v_b2a, q_b2a)
//
// Coordinate conventions:
//   q_b2a = inverse(q_a2w) * q_b2w
//   v_b2a = q_a2w << (v_b2w - v_a2w)     i.e. R_A^T (p_B - p_A)
//   trA2W = Transform3(q_a2w, v_a2w)
//   trB2A = Transform3(trA2W, trB2W)      (convenience constructor for relative transform)
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/** @brief Builds a unit quaternion from an axis and angle (radians). */
static Quaternion<double> obbMakeQuat(double ax, double ay, double az, double angle)
{
    const double s   = std::sin(angle * 0.5);
    const double c   = std::cos(angle * 0.5);
    const double len = std::sqrt(ax * ax + ay * ay + az * az);
    return Quaternion<double>((ax / len) * s, (ay / len) * s, (az / len) * s, c);
}

// -------------------------------------------------------------------------------------------------
/** @brief Calls all four OBB overloads from identical world-space inputs and asserts they all
    return `expected`.  Also checks that the four results are mutually consistent. */
static void checkAllOBBOverloads(const Vector3<double>&    a,
                                 const Vector3<double>&    b,
                                 const Vector3<double>&    v_a2w,
                                 const Quaternion<double>& q_a2w,
                                 const Vector3<double>&    v_b2w,
                                 const Quaternion<double>& q_b2w,
                                 bool                      expected,
                                 const std::string&        label)
{
    // ---- Overload 1: two absolute Transform3 ----
    const Transform3<double> trA2W(q_a2w, v_a2w);
    const Transform3<double> trB2W(q_b2w, v_b2w);
    const bool               res1 = intersectOrientedBoundingBox(a, b, trA2W, trB2W);

    // ---- Overload 2: relative Transform3 (b expressed in a's frame) ----
    const Transform3<double> trB2A(trA2W, trB2W);
    const bool               res2 = intersectOrientedBoundingBox(a, b, trB2A);

    // ---- Overload 3: absolute quat/vec ----
    const bool res3 = intersectOrientedBoundingBox(a, b, v_a2w, v_b2w, q_a2w, q_b2w);

    // ---- Overload 4: relative quat/vec ----
    const Quaternion<double> q_b2a = inverse(q_a2w) * q_b2w;
    const Vector3<double>    v_b2a = q_a2w << (v_b2w - v_a2w);
    const bool               res4  = intersectOrientedBoundingBox(a, b, v_b2a, q_b2a);

    EXPECT_EQ(res1, expected) << "[" << label << "] overload 1 (abs Transform3) wrong";
    EXPECT_EQ(res2, expected) << "[" << label << "] overload 2 (rel Transform3) wrong";
    EXPECT_EQ(res3, expected) << "[" << label << "] overload 3 (abs quat/vec) wrong";
    EXPECT_EQ(res4, expected) << "[" << label << "] overload 4 (rel quat/vec) wrong";

    EXPECT_EQ(res1, res2) << "[" << label << "] overloads 1 vs 2 disagree";
    EXPECT_EQ(res1, res3) << "[" << label << "] overloads 1 vs 3 disagree";
    EXPECT_EQ(res1, res4) << "[" << label << "] overloads 1 vs 4 disagree";
}

// =================================================================================================
// Fixture for OBB intersection tests
// =================================================================================================
class OBBIntersectionTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        q_id      = Quaternion<double>(0.0, 0.0, 0.0, 1.0);
        q_rot90Y  = obbMakeQuat(0.0, 1.0, 0.0, M_PI / 2.0);
        q_rot90X  = obbMakeQuat(1.0, 0.0, 0.0, M_PI / 2.0);
        q_rot45Y  = obbMakeQuat(0.0, 1.0, 0.0, M_PI / 4.0);
        q_rot45Z  = obbMakeQuat(0.0, 0.0, 1.0, M_PI / 4.0);
        q_oblique = obbMakeQuat(1.0, 1.0, 0.0, M_PI / 6.0);
    }

    Quaternion<double> q_id;
    Quaternion<double> q_rot90Y;
    Quaternion<double> q_rot90X;
    Quaternion<double> q_rot45Y;
    Quaternion<double> q_rot45Z;
    Quaternion<double> q_oblique;
};

// -------------------------------------------------------------------------------------------------
// Both boxes at the same origin with identity orientation -> always intersect.
// Expected: true
// -------------------------------------------------------------------------------------------------
TEST_F(OBBIntersectionTest, SamePosition_IdentityOrientation)
{
    const Vector3<double> a(1.0, 1.0, 1.0);
    const Vector3<double> b(0.5, 0.5, 0.5);
    const Vector3<double> origin(0.0, 0.0, 0.0);
    checkAllOBBOverloads(a,
                         b,
                         origin,
                         q_id,
                         origin,
                         q_id,
                         true,
                         "SamePosition_IdentityOrientation");
}

// -------------------------------------------------------------------------------------------------
// B separated beyond extent sum along X.
// Expected: false
// -------------------------------------------------------------------------------------------------
TEST_F(OBBIntersectionTest, SeparatedAlongX)
{
    const Vector3<double> a(1.0, 1.0, 1.0);
    const Vector3<double> b(1.0, 1.0, 1.0);
    checkAllOBBOverloads(a,
                         b,
                         Vector3<double>(0.0, 0.0, 0.0),
                         q_id,
                         Vector3<double>(5.0, 0.0, 0.0),
                         q_id,
                         false,
                         "SeparatedAlongX");
}

// -------------------------------------------------------------------------------------------------
// B just touching A along X (gap = a[0] + b[0] exactly): treated as intersecting
// by the SAT using LOWEPS noise.
// Expected: true
// -------------------------------------------------------------------------------------------------
TEST_F(OBBIntersectionTest, TouchingAlongX)
{
    const Vector3<double> a(1.0, 1.0, 1.0);
    const Vector3<double> b(1.0, 1.0, 1.0);
    checkAllOBBOverloads(a,
                         b,
                         Vector3<double>(0.0, 0.0, 0.0),
                         q_id,
                         Vector3<double>(2.0, 0.0, 0.0),
                         q_id,
                         true,
                         "TouchingAlongX");
}

// -------------------------------------------------------------------------------------------------
// B separated exactly by one extent sum + 0.01 margin: clearly separated.
// Expected: false
// -------------------------------------------------------------------------------------------------
TEST_F(OBBIntersectionTest, JustSeparatedAlongX)
{
    const Vector3<double> a(1.0, 1.0, 1.0);
    const Vector3<double> b(1.0, 1.0, 1.0);
    checkAllOBBOverloads(a,
                         b,
                         Vector3<double>(0.0, 0.0, 0.0),
                         q_id,
                         Vector3<double>(2.01, 0.0, 0.0),
                         q_id,
                         false,
                         "JustSeparatedAlongX");
}

// -------------------------------------------------------------------------------------------------
// B rotated 90 deg around Y but placed at same center -> still intersecting.
// Expected: true
// -------------------------------------------------------------------------------------------------
TEST_F(OBBIntersectionTest, Rotated90Y_SameCenter)
{
    const Vector3<double> a(1.0, 1.0, 1.0);
    const Vector3<double> b(1.0, 1.0, 1.0);
    const Vector3<double> center(0.0, 0.0, 0.0);
    checkAllOBBOverloads(a, b, center, q_id, center, q_rot90Y, true, "Rotated90Y_SameCenter");
}

// -------------------------------------------------------------------------------------------------
// B rotated 45 deg around Z and displaced far in X.
// Expected: false
// -------------------------------------------------------------------------------------------------
TEST_F(OBBIntersectionTest, Rotated45Z_FarAlong_X)
{
    const Vector3<double> a(1.0, 1.0, 1.0);
    const Vector3<double> b(1.0, 1.0, 1.0);
    checkAllOBBOverloads(a,
                         b,
                         Vector3<double>(0.0, 0.0, 0.0),
                         q_id,
                         Vector3<double>(10.0, 0.0, 0.0),
                         q_rot45Z,
                         false,
                         "Rotated45Z_FarAlongX");
}

// -------------------------------------------------------------------------------------------------
// Asymmetric extents: thin slab A, large box B overlapping.
// Expected: true
// -------------------------------------------------------------------------------------------------
TEST_F(OBBIntersectionTest, AsymmetricExtents_Overlapping)
{
    const Vector3<double> a(0.1, 5.0, 5.0);  // thin slab in X
    const Vector3<double> b(3.0, 3.0, 3.0);
    checkAllOBBOverloads(a,
                         b,
                         Vector3<double>(0.0, 0.0, 0.0),
                         q_id,
                         Vector3<double>(0.05, 0.0, 0.0),
                         q_rot45Y,
                         true,
                         "AsymmetricExtents_Overlapping");
}

// -------------------------------------------------------------------------------------------------
// Asymmetric extents: B displaced past A's thin dimension.
// Expected: false
// -------------------------------------------------------------------------------------------------
TEST_F(OBBIntersectionTest, AsymmetricExtents_Separated)
{
    const Vector3<double> a(0.1, 5.0, 5.0);
    const Vector3<double> b(3.0, 3.0, 3.0);
    checkAllOBBOverloads(a,
                         b,
                         Vector3<double>(0.0, 0.0, 0.0),
                         q_id,
                         Vector3<double>(5.0, 0.0, 0.0),
                         q_rot45Y,
                         false,
                         "AsymmetricExtents_Separated");
}

// -------------------------------------------------------------------------------------------------
// Non-identity A orientation, B close: both boxes at same world position, A rotated.
// Expected: true
// -------------------------------------------------------------------------------------------------
TEST_F(OBBIntersectionTest, NonIdentityA_SameCenter)
{
    const Vector3<double> a(1.0, 2.0, 3.0);
    const Vector3<double> b(1.0, 2.0, 3.0);
    const Vector3<double> center(3.0, -1.0, 2.0);
    checkAllOBBOverloads(a, b, center, q_rot90X, center, q_rot45Z, true, "NonIdentityA_SameCenter");
}

// -------------------------------------------------------------------------------------------------
// Non-identity A orientation, B far away.
// Expected: false
// -------------------------------------------------------------------------------------------------
TEST_F(OBBIntersectionTest, NonIdentityA_BFarAway)
{
    const Vector3<double> a(1.0, 1.0, 1.0);
    const Vector3<double> b(1.0, 1.0, 1.0);
    checkAllOBBOverloads(a,
                         b,
                         Vector3<double>(0.0, 0.0, 0.0),
                         q_rot90X,
                         Vector3<double>(20.0, 0.0, 0.0),
                         q_oblique,
                         false,
                         "NonIdentityA_BFarAway");
}

// -------------------------------------------------------------------------------------------------
// Symmetry: swapping A and B must give the same result (intersecting case).
// -------------------------------------------------------------------------------------------------
TEST_F(OBBIntersectionTest, SymmetryAB_Intersecting)
{
    const Vector3<double> a(1.0, 2.0, 0.5);
    const Vector3<double> b(0.5, 1.0, 2.0);

    const bool ab = intersectOrientedBoundingBox(
        a,
        b,
        Transform3<double>(q_id, Vector3<double>(0.0, 0.0, 0.0)),
        Transform3<double>(q_rot45Y, Vector3<double>(1.0, 0.0, 0.0)));

    const bool ba
        = intersectOrientedBoundingBox(b,
                                       a,
                                       Transform3<double>(q_rot45Y, Vector3<double>(1.0, 0.0, 0.0)),
                                       Transform3<double>(q_id, Vector3<double>(0.0, 0.0, 0.0)));

    EXPECT_EQ(ab, ba) << "OBB result must be symmetric under A<->B swap (intersecting case)";
    EXPECT_TRUE(ab);
}

// -------------------------------------------------------------------------------------------------
// Symmetry: swapping A and B must give the same result (separated case).
// -------------------------------------------------------------------------------------------------
TEST_F(OBBIntersectionTest, SymmetryAB_Separated)
{
    const Vector3<double> a(1.0, 1.0, 1.0);
    const Vector3<double> b(1.0, 1.0, 1.0);

    const bool ab = intersectOrientedBoundingBox(
        a,
        b,
        Transform3<double>(q_id, Vector3<double>(0.0, 0.0, 0.0)),
        Transform3<double>(q_rot45Z, Vector3<double>(10.0, 0.0, 0.0)));

    const bool ba = intersectOrientedBoundingBox(
        b,
        a,
        Transform3<double>(q_rot45Z, Vector3<double>(10.0, 0.0, 0.0)),
        Transform3<double>(q_id, Vector3<double>(0.0, 0.0, 0.0)));

    EXPECT_EQ(ab, ba) << "OBB result must be symmetric under A<->B swap (separated case)";
    EXPECT_FALSE(ab);
}

// =================================================================================================
// Random consistency test: 100 random configs, all four overloads must agree.
// =================================================================================================
TEST_F(OBBIntersectionTest, RandomConfigs_AllOverloadsConsistent)
{
    std::mt19937_64                        rng(20260315ULL);
    std::uniform_real_distribution<double> distPos(-5.0, 5.0);
    std::uniform_real_distribution<double> distDim(0.1, 3.0);
    std::uniform_real_distribution<double> distAxis(-1.0, 1.0);
    std::uniform_real_distribution<double> distAngle(0.0, 2.0 * M_PI);

    auto randUnitAxis = [&]() -> std::tuple<double, double, double> {
        double x, y, z, len;
        do
        {
            x   = distAxis(rng);
            y   = distAxis(rng);
            z   = distAxis(rng);
            len = std::sqrt(x * x + y * y + z * z);
        } while(len < 1e-6);
        return {x / len, y / len, z / len};
    };

    auto randQuat = [&]() -> Quaternion<double> {
        auto [ax, ay, az]  = randUnitAxis();
        const double angle = distAngle(rng);
        const double s     = std::sin(angle * 0.5);
        const double c     = std::cos(angle * 0.5);
        return Quaternion<double>(ax * s, ay * s, az * s, c);
    };

    for(int trial = 0; trial < 100; ++trial)
    {
        // Random half-extents
        const Vector3<double> a(distDim(rng), distDim(rng), distDim(rng));
        const Vector3<double> b(distDim(rng), distDim(rng), distDim(rng));

        // Random world poses
        const Vector3<double>    v_a2w(distPos(rng), distPos(rng), distPos(rng));
        const Vector3<double>    v_b2w(distPos(rng), distPos(rng), distPos(rng));
        const Quaternion<double> q_a2w = randQuat();
        const Quaternion<double> q_b2w = randQuat();

        // Overload 1: absolute Transform3
        const Transform3<double> trA2W(q_a2w, v_a2w);
        const Transform3<double> trB2W(q_b2w, v_b2w);
        const bool               res1 = intersectOrientedBoundingBox(a, b, trA2W, trB2W);

        // Overload 2: relative Transform3
        const Transform3<double> trB2A(trA2W, trB2W);
        const bool               res2 = intersectOrientedBoundingBox(a, b, trB2A);

        // Overload 3: absolute quat/vec
        const bool res3 = intersectOrientedBoundingBox(a, b, v_a2w, v_b2w, q_a2w, q_b2w);

        // Overload 4: relative quat/vec
        const Quaternion<double> q_b2a = inverse(q_a2w) * q_b2w;
        const Vector3<double>    v_b2a = q_a2w << (v_b2w - v_a2w);
        const bool               res4  = intersectOrientedBoundingBox(a, b, v_b2a, q_b2a);

        EXPECT_EQ(res1, res2) << "trial " << trial << ": overloads 1 vs 2 disagree";
        EXPECT_EQ(res1, res3) << "trial " << trial << ": overloads 1 vs 3 disagree";
        EXPECT_EQ(res1, res4) << "trial " << trial << ": overloads 1 vs 4 disagree";
    }
}
