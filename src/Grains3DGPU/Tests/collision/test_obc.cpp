// =================================================================================================
// test_obc.cpp -- Unit tests for intersectOrientedBoundingCylinder (OBC)
//
// The function has 4 overloads that differ in how the transforms are supplied:
//   1. Absolute Transform3:  (trA2W, trB2W)      - two absolute Transform3 objects
//   2. Relative Transform3:  (trB2A)              - one relative Transform3 object
//   3. Absolute quat/vec:    (v_a2w, v_b2w, q_a2w, q_b2w)
//   4. Relative quat/vec:    (v_b2a, q_b2a)
//
// Because the overloads use different coordinate systems, for every test we:
//   (a) fix world-space parameters,
//   (b) derive each overload's inputs from the world parameters,
//   (c) assert all four overloads agree AND match the known expected result.
//
// Coordinate conventions used below:
//   q_a2w  : quaternion rotating A's LOCAL frame to WORLD
//   v_a2w  : A's center in WORLD
//   q_b2a  = inverse(q_a2w) * q_b2w
//   v_b2a  = q_a2w << (v_b2w - v_a2w)          i.e. R_A^T (p_B - p_A)
//   trA2W  = Transform3(q_a2w, v_a2w)
//   trB2W  = Transform3(q_b2w, v_b2w)
//   trB2A  = Transform3(trA2W, trB2W)           convenience constructor for relative transform
//
// =================================================================================================

#include "OBC.hh"
#include "Quaternion.hh"
#include "QuaternionMath.hh"
#include "Transform3.hh"
#include "Vector3.hh"

#include <cmath>
#include <gtest/gtest.h>
#include <random>

// =================================================================================================
/** @brief Builds a unit quaternion from an axis (ax, ay, az) and angle (radians). */
// =================================================================================================
static Quaternion<double> makeQuat(double ax, double ay, double az, double angle)
{
    const double s = std::sin(angle * 0.5);
    const double c = std::cos(angle * 0.5);
    // Normalize the axis first
    const double len = std::sqrt(ax * ax + ay * ay + az * az);
    return Quaternion<double>((ax / len) * s, (ay / len) * s, (az / len) * s, c);
}

// =================================================================================================
/** @brief Helper that calls all four OBC overloads with equivalent inputs derived from the
    world-space description and asserts that they all return `expected`.
    @param r1     Radius of cylinder A
    @param h1     Half-height of cylinder A
    @param ori1   Axis direction of A in A's LOCAL frame
    @param r2     Radius of cylinder B
    @param h2     Half-height of cylinder B
    @param ori2   Axis direction of B in B's LOCAL frame
    @param v_a2w  A's center in world space
    @param q_a2w  Quaternion rotating A's local frame to world
    @param v_b2w  B's center in world space
    @param q_b2w  Quaternion rotating B's local frame to world
    @param expected Expected intersection result
    @param label  Human-readable test name for error messages */
// =================================================================================================
static void checkAllOBCOverloads(double                    r1,
                                 double                    h1,
                                 const Vector3<double>&    ori1,
                                 double                    r2,
                                 double                    h2,
                                 const Vector3<double>&    ori2,
                                 const Vector3<double>&    v_a2w,
                                 const Quaternion<double>& q_a2w,
                                 const Vector3<double>&    v_b2w,
                                 const Quaternion<double>& q_b2w,
                                 bool                      expected,
                                 const std::string&        label)
{
    // ---- Overload 1: two absolute Transform3 objects ----
    const Transform3<double> trA2W(q_a2w, v_a2w);
    const Transform3<double> trB2W(q_b2w, v_b2w);
    const bool res1 = intersectOrientedBoundingCylinder(r1, h1, ori1, r2, h2, ori2, trA2W, trB2W);

    // ---- Overload 2: one relative Transform3 (b2a) ----
    // Convenience constructor: Transform3(t1, t2) = t2 expressed in t1's local frame
    const Transform3<double> trB2A(trA2W, trB2W);
    const bool res2 = intersectOrientedBoundingCylinder(r1, h1, ori1, r2, h2, ori2, trB2A);

    // ---- Overload 3: absolute quaternion + translation ----
    const bool res3
        = intersectOrientedBoundingCylinder(r1, h1, ori1, r2, h2, ori2, v_a2w, v_b2w, q_a2w, q_b2w);

    // ---- Overload 4: relative quaternion + translation ----
    // q_b2a = inv(q_a2w) * q_b2w
    // v_b2a = R_A^T * (p_B - p_A), computed via inverse quaternion rotation
    const Quaternion<double> q_b2a = inverse(q_a2w) * q_b2w;
    const Vector3<double>    v_b2a = q_a2w << (v_b2w - v_a2w);
    const bool res4 = intersectOrientedBoundingCylinder(r1, h1, ori1, r2, h2, ori2, v_b2a, q_b2a);

    // --- Assertions: all overloads must match the expected result ---
    EXPECT_EQ(res1, expected) << "[" << label << "] overload 1 (abs Transform3) wrong";
    EXPECT_EQ(res2, expected) << "[" << label << "] overload 2 (rel Transform3) wrong";
    EXPECT_EQ(res3, expected) << "[" << label << "] overload 3 (abs quat/vec) wrong";
    EXPECT_EQ(res4, expected) << "[" << label << "] overload 4 (rel quat/vec) wrong";

    // --- Cross-consistency: all overloads must agree with each other ---
    EXPECT_EQ(res1, res2) << "[" << label << "] overloads 1 vs 2 disagree";
    EXPECT_EQ(res1, res3) << "[" << label << "] overloads 1 vs 3 disagree";
    EXPECT_EQ(res1, res4) << "[" << label << "] overloads 1 vs 4 disagree";
}

// =================================================================================================
// Test fixture
// =================================================================================================
class OBCSimpleTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // Both cylinders use their local Z-axis as the cylinder axis
        ori_z = Vector3<double>(0.0, 0.0, 1.0);

        // Identity quaternion: no rotation
        q_id = Quaternion<double>(0.0, 0.0, 0.0, 1.0);

        // 90 deg rotation around Y: maps local +Z to world +X
        // q = (0, sin45 deg, 0, cos45 deg)
        q_rot90Y = makeQuat(0.0, 1.0, 0.0, M_PI / 2.0);

        // 90 deg rotation around X: maps local +Z to world -Y
        // Rx(90 deg) * (0,0,1) = (0, -1, 0)
        q_rot90X = makeQuat(1.0, 0.0, 0.0, M_PI / 2.0);

        // 45 deg rotation around Y: maps local +Z to world (sqrt(2)/2, 0, sqrt(2)/2)
        q_rot45Y = makeQuat(0.0, 1.0, 0.0, M_PI / 4.0);

        // 30 deg rotation around an oblique axis (1/sqrt(2), 1/sqrt(2), 0) -- arbitrary orientation
        q_rot30_oblique = makeQuat(1.0, 1.0, 0.0, M_PI / 6.0);
    }

    Vector3<double>    ori_z;            // local cylinder axis = +Z
    Quaternion<double> q_id;             // identity
    Quaternion<double> q_rot90Y;         // 90 deg around world Y
    Quaternion<double> q_rot90X;         // 90 deg around world X
    Quaternion<double> q_rot45Y;         // 45 deg around world Y
    Quaternion<double> q_rot30_oblique;  // 30 deg around (1,1,0)/sqrt(2)
};

// =================================================================================================
// Tests -- Cylinder A is Z-aligned at origin (identity), cylinder B varies.
// Both axes are non-parallel in all tests so the algorithm is well-conditioned.
// =================================================================================================

// -------------------------------------------------------------------------------------------------
// Clearly separated: B is X-aligned (A perp B) and displaced far to the side (+Y direction).
// The lateral offset (Y) exceeds r1+r2, so the cylindrical hulls cannot overlap.
// Expected: false
// -------------------------------------------------------------------------------------------------
TEST_F(OBCSimpleTest, PerpendicularAxes_RadiallySeparated)
{
    // A: Z-aligned at origin   B: X-aligned (rot90Y) at (0, 5, 0)
    checkAllOBCOverloads(1.0,
                         1.0,
                         ori_z,  // cylinder A
                         1.0,
                         1.0,
                         ori_z,  // cylinder B (local axis also Z, rotated to X in world)
                         Vector3<double>(0.0, 0.0, 0.0),
                         q_id,  // A: world center & orientation
                         Vector3<double>(0.0, 5.0, 0.0),
                         q_rot90Y,  // B: world center & orientation
                         false,
                         "PerpendicularAxes_RadiallySeparated");
}

// -------------------------------------------------------------------------------------------------
// Clearly overlapping: B is X-aligned with the same center as A.
// Axes are perpendicular and centers coincide -> deep intersection.
// Expected: true
// -------------------------------------------------------------------------------------------------
TEST_F(OBCSimpleTest, PerpendicularAxes_SameCenter)
{
    // A: Z-aligned at origin   B: X-aligned at same origin
    checkAllOBCOverloads(1.0,
                         2.0,
                         ori_z,
                         1.0,
                         2.0,
                         ori_z,
                         Vector3<double>(0.0, 0.0, 0.0),
                         q_id,
                         Vector3<double>(0.0, 0.0, 0.0),
                         q_rot90Y,
                         true,
                         "PerpendicularAxes_SameCenter");
}

// -------------------------------------------------------------------------------------------------
// Separated axially: B is X-aligned but displaced far from A along Z.
// A extends +/-1 in Z; B is far above at z=5, outside A's extent.
// Expected: false
// -------------------------------------------------------------------------------------------------
TEST_F(OBCSimpleTest, PerpendicularAxes_AxiallySeparated)
{
    // A: Z-aligned at origin  B: X-aligned, displaced 5 units along Z
    checkAllOBCOverloads(1.0,
                         1.0,
                         ori_z,
                         1.0,
                         1.0,
                         ori_z,
                         Vector3<double>(0.0, 0.0, 0.0),
                         q_id,
                         Vector3<double>(0.0, 0.0, 5.0),
                         q_rot90Y,
                         false,
                         "PerpendicularAxes_AxiallySeparated");
}

// -------------------------------------------------------------------------------------------------
// Partial overlap: B is X-aligned and slightly offset in Y (within reach).
// B's center is only 1.5 units from A's axis, which is less than r1+r2=2.
// Expected: true
// -------------------------------------------------------------------------------------------------
TEST_F(OBCSimpleTest, PerpendicularAxes_PartialOverlap)
{
    // A: Z-aligned at origin  B: X-aligned, center at (0, 1.5, 0)
    checkAllOBCOverloads(1.0,
                         2.0,
                         ori_z,
                         1.0,
                         2.0,
                         ori_z,
                         Vector3<double>(0.0, 0.0, 0.0),
                         q_id,
                         Vector3<double>(0.0, 1.5, 0.0),
                         q_rot90Y,
                         true,
                         "PerpendicularAxes_PartialOverlap");
}

// -------------------------------------------------------------------------------------------------
// 45-degree axes, same center: largest possible overlap scenario.
// Expected: true
// -------------------------------------------------------------------------------------------------
TEST_F(OBCSimpleTest, FortyFiveDegreeAxes_SameCenter)
{
    // A: Z-aligned at origin  B: tilted 45 deg around Y (world axis = (sqrt(2)/2, 0, sqrt(2)/2))
    checkAllOBCOverloads(1.0,
                         2.0,
                         ori_z,
                         1.0,
                         2.0,
                         ori_z,
                         Vector3<double>(0.0, 0.0, 0.0),
                         q_id,
                         Vector3<double>(0.0, 0.0, 0.0),
                         q_rot45Y,
                         true,
                         "FortyFiveDegreeAxes_SameCenter");
}

// -------------------------------------------------------------------------------------------------
// 45-degree axes, B far away: clearly no intersection.
// Expected: false
// -------------------------------------------------------------------------------------------------
TEST_F(OBCSimpleTest, FortyFiveDegreeAxes_ClearSeparation)
{
    // A: Z-aligned at origin  B: tilted 45 deg around Y, displaced 10 units in X
    checkAllOBCOverloads(1.0,
                         1.0,
                         ori_z,
                         1.0,
                         1.0,
                         ori_z,
                         Vector3<double>(0.0, 0.0, 0.0),
                         q_id,
                         Vector3<double>(10.0, 0.0, 0.0),
                         q_rot45Y,
                         false,
                         "FortyFiveDegreeAxes_ClearSeparation");
}

// -------------------------------------------------------------------------------------------------
// B axis along -Y (rot90X), same center as A.
// Axes are perpendicular in a different plane than the rot90Y case.
// Expected: true
// -------------------------------------------------------------------------------------------------
TEST_F(OBCSimpleTest, PerpendicularAxes_Rotated90X_SameCenter)
{
    // A: Z-aligned at origin  B: rotated 90 deg around X -> world axis = (0, -1, 0)
    checkAllOBCOverloads(1.0,
                         2.0,
                         ori_z,
                         1.0,
                         2.0,
                         ori_z,
                         Vector3<double>(0.0, 0.0, 0.0),
                         q_id,
                         Vector3<double>(0.0, 0.0, 0.0),
                         q_rot90X,
                         true,
                         "PerpendicularAxes_Rotated90X_SameCenter");
}

// -------------------------------------------------------------------------------------------------
// B axis along -Y (rot90X), B is displaced far in Z.
// Expected: false
// -------------------------------------------------------------------------------------------------
TEST_F(OBCSimpleTest, PerpendicularAxes_Rotated90X_AxiallySeparated)
{
    // A: Z-aligned at origin  B: -Y-aligned, displaced 10 units in Z
    checkAllOBCOverloads(1.0,
                         1.0,
                         ori_z,
                         1.0,
                         1.0,
                         ori_z,
                         Vector3<double>(0.0, 0.0, 0.0),
                         q_id,
                         Vector3<double>(0.0, 0.0, 10.0),
                         q_rot90X,
                         false,
                         "PerpendicularAxes_Rotated90X_AxiallySeparated");
}

// =================================================================================================
// Tests with non-identity A orientation -- stress-tests the coordinate transformation.
// A is -Y-aligned at (1, 2, 3) (rotated 90 deg around X).  B varies.
// =================================================================================================

// -------------------------------------------------------------------------------------------------
// A is -Y-aligned at (1,2,3); B is Z-aligned at (1,2,3.5).
// Physically: B's center is close to A's center, well inside both hulls -> intersect.
// Expected: true
// -------------------------------------------------------------------------------------------------
TEST_F(OBCSimpleTest, NonIdentityA_BClose_Intersecting)
{
    //  A: center (1,2,3), -Y-aligned (rot90X)
    //  B: center (1,2,3.5), Z-aligned (identity)
    checkAllOBCOverloads(1.0,
                         1.0,
                         ori_z,
                         1.0,
                         1.0,
                         ori_z,
                         Vector3<double>(1.0, 2.0, 3.0),
                         q_rot90X,
                         Vector3<double>(1.0, 2.0, 3.5),
                         q_id,
                         true,
                         "NonIdentityA_BClose_Intersecting");
}

// -------------------------------------------------------------------------------------------------
// A is -Y-aligned at (1,2,3); B is Z-aligned at (1,2,7).
// B is 4 units above A's center in Z. A's half-height is 1 -> axes are far apart.
// Expected: false
// -------------------------------------------------------------------------------------------------
TEST_F(OBCSimpleTest, NonIdentityA_BFar_Separated)
{
    //  A: center (1,2,3), -Y-aligned (rot90X)
    //  B: center (1,2,7), Z-aligned (identity)
    checkAllOBCOverloads(1.0,
                         1.0,
                         ori_z,
                         1.0,
                         1.0,
                         ori_z,
                         Vector3<double>(1.0, 2.0, 3.0),
                         q_rot90X,
                         Vector3<double>(1.0, 2.0, 7.0),
                         q_id,
                         false,
                         "NonIdentityA_BFar_Separated");
}

// =================================================================================================
// Tests with non-identity orientation on BOTH cylinders -- maximum transformation stress.
// =================================================================================================

// -------------------------------------------------------------------------------------------------
// Both cylinders rotated and translated to non-trivial world poses, but placed
// at the same center.  Any two cylinders sharing a center must intersect as long
// as their dimensions are positive.
// Expected: true
// -------------------------------------------------------------------------------------------------
TEST_F(OBCSimpleTest, BothRotated_SameCenter_MustIntersect)
{
    //  A: center (3,-1,2), rotated 30 deg around oblique axis (1,1,0)/sqrt(2)
    //  B: center (3,-1,2), rotated 90 deg around Y
    const Vector3<double> center(3.0, -1.0, 2.0);
    checkAllOBCOverloads(1.5,
                         2.0,
                         ori_z,
                         1.5,
                         2.0,
                         ori_z,
                         center,
                         q_rot30_oblique,
                         center,
                         q_rot90Y,
                         true,
                         "BothRotated_SameCenter_MustIntersect");
}

// -------------------------------------------------------------------------------------------------
// Both cylinders rotated, but B is far from A along the world X axis.
// The separation (20 units) dwarfs any cylinder dimension -> always false.
// Expected: false
// -------------------------------------------------------------------------------------------------
TEST_F(OBCSimpleTest, BothRotated_FarApart_MustNotIntersect)
{
    //  A: center (-5,0,0), rotated 90 deg around X
    //  B: center (15,0,0), rotated 45 deg around Y
    checkAllOBCOverloads(2.0,
                         3.0,
                         ori_z,
                         2.0,
                         3.0,
                         ori_z,
                         Vector3<double>(-5.0, 0.0, 0.0),
                         q_rot90X,
                         Vector3<double>(15.0, 0.0, 0.0),
                         q_rot45Y,
                         false,
                         "BothRotated_FarApart_MustNotIntersect");
}

// -------------------------------------------------------------------------------------------------
// A rotated 30 deg (oblique) at (-2,4,-1), B rotated 90 degX at (-2,4,-1.3).
// Centers are only 0.3 units apart while both radii are 1.0 and half-heights 1.5.
// Expected: true
// -------------------------------------------------------------------------------------------------
TEST_F(OBCSimpleTest, BothRotated_CloseCenters_Intersecting)
{
    checkAllOBCOverloads(1.0,
                         1.5,
                         ori_z,
                         1.0,
                         1.5,
                         ori_z,
                         Vector3<double>(-2.0, 4.0, -1.0),
                         q_rot30_oblique,
                         Vector3<double>(-2.0, 4.0, -1.3),
                         q_rot90X,
                         true,
                         "BothRotated_CloseCenters_Intersecting");
}

// -------------------------------------------------------------------------------------------------
// A rotated 30 deg (oblique) at (0,0,0), B rotated 90 degX at (0,0,10).
// Clearly separated: 10-unit gap far exceeds both half-heights.
// Expected: false
// -------------------------------------------------------------------------------------------------
TEST_F(OBCSimpleTest, BothRotated_LargeAxialGap_Separated)
{
    checkAllOBCOverloads(1.0,
                         1.0,
                         ori_z,
                         1.0,
                         1.0,
                         ori_z,
                         Vector3<double>(0.0, 0.0, 0.0),
                         q_rot30_oblique,
                         Vector3<double>(0.0, 0.0, 10.0),
                         q_rot90X,
                         false,
                         "BothRotated_LargeAxialGap_Separated");
}

// =================================================================================================
// Tests with different cylinder sizes
// =================================================================================================

// -------------------------------------------------------------------------------------------------
// Large A, tiny B at A's center -- B is completely inside A.
// Expected: true
// -------------------------------------------------------------------------------------------------
TEST_F(OBCSimpleTest, LargeCylinderContainsSmall_SameCenter)
{
    checkAllOBCOverloads(5.0,
                         5.0,
                         ori_z,  // large cylinder A
                         0.1,
                         0.1,
                         ori_z,  // tiny cylinder B, rotated to be X-aligned
                         Vector3<double>(0.0, 0.0, 0.0),
                         q_id,
                         Vector3<double>(0.0, 0.0, 0.0),
                         q_rot90Y,
                         true,
                         "LargeCylinderContainsSmall_SameCenter");
}

// -------------------------------------------------------------------------------------------------
// Tiny B clearly outside large A.
// Expected: false
// -------------------------------------------------------------------------------------------------
TEST_F(OBCSimpleTest, LargeCylinder_TinyBOutside)
{
    checkAllOBCOverloads(5.0,
                         5.0,
                         ori_z,
                         0.1,
                         0.1,
                         ori_z,
                         Vector3<double>(0.0, 0.0, 0.0),
                         q_id,
                         Vector3<double>(50.0, 0.0, 0.0),
                         q_rot90Y,
                         false,
                         "LargeCylinder_TinyBOutside");
}

// =================================================================================================
// Test that swapping A and B produces the same collision result (symmetry).
// =================================================================================================

// -------------------------------------------------------------------------------------------------
// Symmetry check: result must be the same when A and B are swapped.
// -------------------------------------------------------------------------------------------------
TEST_F(OBCSimpleTest, SymmetryAB_Intersecting)
{
    // Configuration AB: A Z-aligned at origin, B X-aligned at (0,0,0)
    const bool ab = intersectOrientedBoundingCylinder(
        1.0,
        2.0,
        ori_z,
        1.0,
        2.0,
        ori_z,
        Transform3<double>(q_id, Vector3<double>(0.0, 0.0, 0.0)),
        Transform3<double>(q_rot90Y, Vector3<double>(0.0, 0.0, 0.0)));

    // Configuration BA: B first (X-aligned at origin), A second (Z-aligned at origin)
    const bool ba = intersectOrientedBoundingCylinder(
        1.0,
        2.0,
        ori_z,
        1.0,
        2.0,
        ori_z,
        Transform3<double>(q_rot90Y, Vector3<double>(0.0, 0.0, 0.0)),
        Transform3<double>(q_id, Vector3<double>(0.0, 0.0, 0.0)));

    EXPECT_EQ(ab, ba) << "OBC result must be symmetric under A<->B swap (intersecting case)";
    EXPECT_TRUE(ab);
}

// -------------------------------------------------------------------------------------------------
// Symmetry check for a separated case.
// -------------------------------------------------------------------------------------------------
TEST_F(OBCSimpleTest, SymmetryAB_Separated)
{
    // Configuration AB: A Z-aligned at origin, B X-aligned far away
    const bool ab = intersectOrientedBoundingCylinder(
        1.0,
        1.0,
        ori_z,
        1.0,
        1.0,
        ori_z,
        Transform3<double>(q_id, Vector3<double>(0.0, 0.0, 0.0)),
        Transform3<double>(q_rot90Y, Vector3<double>(10.0, 0.0, 0.0)));

    // Configuration BA: swap A and B
    const bool ba = intersectOrientedBoundingCylinder(
        1.0,
        1.0,
        ori_z,
        1.0,
        1.0,
        ori_z,
        Transform3<double>(q_rot90Y, Vector3<double>(10.0, 0.0, 0.0)),
        Transform3<double>(q_id, Vector3<double>(0.0, 0.0, 0.0)));

    EXPECT_EQ(ab, ba) << "OBC result must be symmetric under A<->B swap (separated case)";
    EXPECT_FALSE(ab);
}

// =================================================================================================
// Test that all four overloads are mutually consistent even when ori1/ori2 are not the
// canonical Z-axis, i.e. the cylinder is not "upright" in its own local frame.
// =================================================================================================

// -------------------------------------------------------------------------------------------------
// Non-Z local axis: both cylinders lie along their local X-axis (ori = (1,0,0)).
// A stays at identity (so world axis = X), B is rotated 90 deg around Z
// (world axis = Y).  They share the same center -> must intersect.
// Expected: true
// -------------------------------------------------------------------------------------------------
TEST_F(OBCSimpleTest, NonZLocalAxis_SameCenter)
{
    const Vector3<double> ori_x(1.0, 0.0, 0.0);

    // 90 deg rotation around Z maps local +X to world +... let's use the identity
    // quaternion so A's world axis stays +X.
    // B is rotated 90 deg around Z: Rz(90 deg) maps +X -> +Y.
    const Quaternion<double> q_rot90Z = makeQuat(0.0, 0.0, 1.0, M_PI / 2.0);

    checkAllOBCOverloads(1.0,
                         2.0,
                         ori_x,
                         1.0,
                         2.0,
                         ori_x,
                         Vector3<double>(0.0, 0.0, 0.0),
                         q_id,
                         Vector3<double>(0.0, 0.0, 0.0),
                         q_rot90Z,
                         true,
                         "NonZLocalAxis_SameCenter");
}

// -------------------------------------------------------------------------------------------------
// Non-Z local axis, B separated far in X.
// Expected: false
// -------------------------------------------------------------------------------------------------
TEST_F(OBCSimpleTest, NonZLocalAxis_Separated)
{
    const Vector3<double>    ori_x(1.0, 0.0, 0.0);
    const Quaternion<double> q_rot90Z = makeQuat(0.0, 0.0, 1.0, M_PI / 2.0);

    checkAllOBCOverloads(1.0,
                         1.0,
                         ori_x,
                         1.0,
                         1.0,
                         ori_x,
                         Vector3<double>(0.0, 0.0, 0.0),
                         q_id,
                         Vector3<double>(20.0, 0.0, 0.0),
                         q_rot90Z,
                         false,
                         "NonZLocalAxis_Separated");
}

// =================================================================================================
// Random consistency test: generate 100 random configurations and verify all four overloads
// agree on the result.  No expected boolean is asserted -- the test only checks mutual
// consistency between overloads.  A fixed seed makes failures reproducible.
// =================================================================================================
TEST_F(OBCSimpleTest, RandomConfigs_AllOverloadsConsistent)
{
    // Fixed seed for reproducibility
    std::mt19937_64                        rng(20260315ULL);
    std::uniform_real_distribution<double> distPos(-5.0, 5.0);          // world position components
    std::uniform_real_distribution<double> distDim(0.1, 3.0);           // radii / half-heights
    std::uniform_real_distribution<double> distAxis(-1.0, 1.0);         // raw axis components
    std::uniform_real_distribution<double> distAngle(0.0, 2.0 * M_PI);  // rotation angle

    // Helper: random unit vector (resamples if near-zero)
    auto randUnitVec = [&]() -> Vector3<double> {
        double x, y, z, len;
        do
        {
            x   = distAxis(rng);
            y   = distAxis(rng);
            z   = distAxis(rng);
            len = std::sqrt(x * x + y * y + z * z);
        } while(len < 1e-6);
        return Vector3<double>(x / len, y / len, z / len);
    };

    // Helper: random unit quaternion from axis-angle
    auto randQuat = [&]() -> Quaternion<double> {
        const Vector3<double> ax    = randUnitVec();
        const double          angle = distAngle(rng);
        const double          s     = std::sin(angle * 0.5);
        const double          c     = std::cos(angle * 0.5);
        return Quaternion<double>(ax[0] * s, ax[1] * s, ax[2] * s, c);
    };

    for(int trial = 0; trial < 100; ++trial)
    {
        // Random cylinder dimensions
        const double r1 = distDim(rng);
        const double h1 = distDim(rng);
        const double r2 = distDim(rng);
        const double h2 = distDim(rng);

        // Random local cylinder axes (unit vectors)
        const Vector3<double> ori1 = randUnitVec();
        const Vector3<double> ori2 = randUnitVec();

        // Random world poses
        const Vector3<double>    v_a2w(distPos(rng), distPos(rng), distPos(rng));
        const Vector3<double>    v_b2w(distPos(rng), distPos(rng), distPos(rng));
        const Quaternion<double> q_a2w = randQuat();
        const Quaternion<double> q_b2w = randQuat();

        // ----- Overload 1: absolute Transform3 -----
        const Transform3<double> trA2W(q_a2w, v_a2w);
        const Transform3<double> trB2W(q_b2w, v_b2w);
        const bool               res1
            = intersectOrientedBoundingCylinder(r1, h1, ori1, r2, h2, ori2, trA2W, trB2W);

        // ----- Overload 2: relative Transform3 -----
        const Transform3<double> trB2A(trA2W, trB2W);
        const bool res2 = intersectOrientedBoundingCylinder(r1, h1, ori1, r2, h2, ori2, trB2A);

        // ----- Overload 3: absolute quat/vec -----
        const bool res3 = intersectOrientedBoundingCylinder(r1,
                                                            h1,
                                                            ori1,
                                                            r2,
                                                            h2,
                                                            ori2,
                                                            v_a2w,
                                                            v_b2w,
                                                            q_a2w,
                                                            q_b2w);

        // ----- Overload 4: relative quat/vec -----
        const Quaternion<double> q_b2a = inverse(q_a2w) * q_b2w;
        const Vector3<double>    v_b2a = q_a2w << (v_b2w - v_a2w);
        const bool               res4
            = intersectOrientedBoundingCylinder(r1, h1, ori1, r2, h2, ori2, v_b2a, q_b2a);

        // All four must agree
        EXPECT_EQ(res1, res2) << "trial " << trial << ": overloads 1 vs 2 disagree";
        EXPECT_EQ(res1, res3) << "trial " << trial << ": overloads 1 vs 3 disagree";
        EXPECT_EQ(res1, res4) << "trial " << trial << ": overloads 1 vs 4 disagree";
    }
}
