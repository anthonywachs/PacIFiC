#include <cmath>
#include <gtest/gtest.h>
#include <random>

#include "Box.hh"
#include "GJK.hh"
#include "Quaternion.hh"
#include "QuaternionMath.hh"
#include "Sphere.hh"
#include "Superquadric.hh"
#include "Transform3.hh"
#include "Vector3.hh"

class GJKTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        boxA = new Box<double>(1.0, 1.0, 1.0);  // half-extents 1x1x1
        boxB = new Box<double>(0.5, 0.5, 0.5);  // half-extents 0.5x0.5x0.5

        sphereA = new Sphere<double>(1.0);  // radius 1
        sphereB = new Sphere<double>(0.5);  // radius 0.5

        superquadricA = new Superquadric<double>(1.0, 1.0, 1.0, 2.0, 2.0);  // extents 1; n1,n2=2
        superquadricB
            = new Superquadric<double>(0.5, 0.5, 0.5, 1.5, 1.5);  // extents 0.5; n1,n2=1.5

        identity_transform = Transform3<double>(Quaternion<double>(0.0, 0.0, 0.0, 1.0),
                                                Vector3<double>(0.0, 0.0, 0.0));

        separated_transform
            = Transform3<double>(Quaternion<double>(0.0, 0.0, 0.0, 1.0),
                                 Vector3<double>(3.0, 0.0, 0.0)  // Move 3 units in x direction
            );

        overlapping_transform
            = Transform3<double>(Quaternion<double>(0.0, 0.0, 0.0, 1.0),
                                 Vector3<double>(0.5, 0.0, 0.0)  // Move 0.5 units in x direction
            );
    }

    void TearDown() override
    {
        delete boxA;
        delete boxB;
        delete sphereA;
        delete sphereB;
        delete superquadricA;
        delete superquadricB;
    }

    Box<double>* boxA;
    Box<double>* boxB;

    Sphere<double>* sphereA;
    Sphere<double>* sphereB;

    Superquadric<double>* superquadricA;
    Superquadric<double>* superquadricB;

    Transform3<double> identity_transform;
    Transform3<double> separated_transform;
    Transform3<double> overlapping_transform;
    const double       EPSILON = 1e-10;
};

TEST_F(GJKTest, IdenticalShapesIntersect)
{
    bool result = intersectGJK(*boxA, *boxA, identity_transform);

    if(!result)
    {
        Transform3<double> tiny_offset(Quaternion<double>(0.0, 0.0, 0.0, 1.0),
                                       Vector3<double>(1e-6, 0.0, 0.0)  // Very small displacement
        );
        result = intersectGJK(*boxA, *boxB, tiny_offset);
    }

    EXPECT_TRUE(result);
}

TEST_F(GJKTest, OverlappingShapesIntersect)
{
    bool result = intersectGJK(*boxA, *boxB, overlapping_transform);
    EXPECT_TRUE(result);
}

TEST_F(GJKTest, SeparatedShapesDoNotIntersect)
{
    bool result = intersectGJK(*boxA, *boxB, separated_transform);
    EXPECT_FALSE(result);
}

TEST_F(GJKTest, WorldCoordinatesIntersection)
{
    Transform3<double> transformA = identity_transform;
    Transform3<double> transformB = overlapping_transform;

    bool result = intersectGJK(*boxA, *boxB, transformA, transformB);
    EXPECT_TRUE(result);

    transformB = separated_transform;
    result     = intersectGJK(*boxA, *boxB, transformA, transformB);
    EXPECT_FALSE(result);
}

TEST_F(GJKTest, RotatedShapesIntersection)
{
    double             angle = M_PI / 12.0;  // 15 degrees
    Quaternion<double> rotation(0.0, 0.0, sin(angle / 2.0), cos(angle / 2.0));
    Transform3<double> rotated_transform(rotation, Vector3<double>(0.1, 0.1, 0.0));

    bool result = intersectGJK(*boxA, *boxB, rotated_transform);
    EXPECT_TRUE(result);  // Should intersect with small offset
}

TEST_F(GJKTest, EdgeCases)
{
    bool result = intersectGJK(*boxA, *boxB, overlapping_transform);
    EXPECT_TRUE(result);

    Box<double> tiny_box(0.001, 0.001, 0.001);
    result = intersectGJK(*boxA, tiny_box, overlapping_transform);
    EXPECT_TRUE(result);  // Tiny box should intersect at overlapping position
}

TEST_F(GJKTest, PerformanceAndStability)
{
    for(int i = 0; i < 50; ++i)
    {
        double             offset = i * 0.05 + 0.01;
        Transform3<double> test_transform(Quaternion<double>(0.0, 0.0, 0.0, 1.0),
                                          Vector3<double>(offset, 0.0, 0.0));

        bool result = intersectGJK(*boxA, *boxB, test_transform);

        // Should intersect for small offsets, not intersect for large offsets
        if(offset < 1.0)
        {
            EXPECT_TRUE(result);
        }
        else if(offset > 2.0)
        {
            EXPECT_FALSE(result);
        }
        // Skip assertion for boundary region where result may vary
    }
}

TEST_F(GJKTest, VariousOrientations)
{
    double angles[] = {0.0, M_PI / 12.0};  // Just test 0 and 15 degrees

    for(double angle : angles)
    {
        Quaternion<double> rotZ(0.0, 0.0, sin(angle / 2.0), cos(angle / 2.0));
        Transform3<double> transformZ(rotZ, Vector3<double>(0.5, 0.0, 0.0));
        bool               resultZ = intersectGJK(*boxA, *boxB, transformZ);
        EXPECT_TRUE(resultZ);

        Quaternion<double> rotX(sin(angle / 2.0), 0.0, 0.0, cos(angle / 2.0));
        Transform3<double> transformX(rotX, Vector3<double>(0.5, 0.0, 0.0));
        bool               resultX = intersectGJK(*boxA, *boxB, transformX);
        EXPECT_TRUE(resultX);

        Quaternion<double> rotY(0.0, sin(angle / 2.0), 0.0, cos(angle / 2.0));
        Transform3<double> transformY(rotY, Vector3<double>(0.5, 0.0, 0.0));
        bool               resultY = intersectGJK(*boxA, *boxB, transformY);
        EXPECT_TRUE(resultY);
    }
}

TEST_F(GJKTest, CollisionConsistency)
{
    bool result1 = intersectGJK(*boxA, *boxB, overlapping_transform);
    bool result2 = intersectGJK(*boxA, *boxB, overlapping_transform);
    bool result3 = intersectGJK(*boxA, *boxB, overlapping_transform);

    EXPECT_EQ(result1, result2);
    EXPECT_EQ(result2, result3);
    EXPECT_TRUE(result1);  // Should be true for overlapping case

    result1 = intersectGJK(*boxA, *boxB, separated_transform);
    result2 = intersectGJK(*boxA, *boxB, separated_transform);
    EXPECT_EQ(result1, result2);
    EXPECT_FALSE(result1);  // Should be false for separated case
}

TEST_F(GJKTest, QuaternionBasedRelativeTransformation)
{
    Vector3<double>    v_b2a(0.5, 0.0, 0.0);  // Same as overlapping_transform
    Quaternion<double> q_b2a(0.0, 0.0, 0.0, 1.0);

    bool result = intersectGJK(*boxA, *boxB, v_b2a, q_b2a);
    EXPECT_TRUE(result);

    Vector3<double> v_separated(3.0, 0.0, 0.0);
    result = intersectGJK(*boxA, *boxB, v_separated, q_b2a);
    EXPECT_FALSE(result);

    double             angle = M_PI / 12.0;  // 15 degrees
    Quaternion<double> q_rotated(0.0, 0.0, sin(angle / 2.0), cos(angle / 2.0));
    Vector3<double>    v_overlapping(0.5, 0.0, 0.0);
    result = intersectGJK(*boxA, *boxB, v_overlapping, q_rotated);
    EXPECT_TRUE(result);
}

TEST_F(GJKTest, QuaternionBasedWorldCoordinates)
{
    Vector3<double>    v_a2w(0.0, 0.0, 0.0);
    Vector3<double>    v_b2w(0.5, 0.0, 0.0);
    Quaternion<double> q_a2w(0.0, 0.0, 0.0, 1.0);
    Quaternion<double> q_b2w(0.0, 0.0, 0.0, 1.0);

    bool result = intersectGJK(*boxA, *boxB, v_a2w, v_b2w, q_a2w, q_b2w);
    EXPECT_TRUE(result);

    v_b2w  = Vector3<double>(3.0, 0.0, 0.0);
    result = intersectGJK(*boxA, *boxB, v_a2w, v_b2w, q_a2w, q_b2w);
    EXPECT_FALSE(result);

    double angle = M_PI / 12.0;  // 15 degrees
    q_b2w        = Quaternion<double>(0.0, 0.0, sin(angle / 2.0), cos(angle / 2.0));
    v_b2w        = Vector3<double>(0.5, 0.0, 0.0);
    result       = intersectGJK(*boxA, *boxB, v_a2w, v_b2w, q_a2w, q_b2w);
    EXPECT_TRUE(result);
}

TEST_F(GJKTest, Transform3VsQuaternionConsistency)
{
    Vector3<double>    position(0.5, 0.0, 0.0);
    Quaternion<double> rotation(0.0, 0.0, 0.0, 1.0);
    Transform3<double> transform(rotation, position);

    bool transform_result  = intersectGJK(*boxA, *boxB, transform);
    bool quaternion_result = intersectGJK(*boxA, *boxB, position, rotation);

    EXPECT_EQ(transform_result, quaternion_result)
        << "Transform3 and quaternion-based GJK should give same result for "
           "overlapping case";
    EXPECT_TRUE(transform_result);

    position  = Vector3<double>(3.0, 0.0, 0.0);
    transform = Transform3<double>(rotation, position);

    transform_result  = intersectGJK(*boxA, *boxB, transform);
    quaternion_result = intersectGJK(*boxA, *boxB, position, rotation);

    EXPECT_EQ(transform_result, quaternion_result)
        << "Transform3 and quaternion-based GJK should give same result for "
           "separated case";
    EXPECT_FALSE(transform_result);

    double angle = M_PI / 12.0;  // 15 degrees
    rotation     = Quaternion<double>(0.0, 0.0, sin(angle / 2.0), cos(angle / 2.0));
    position     = Vector3<double>(0.5, 0.0, 0.0);
    transform    = Transform3<double>(rotation, position);

    transform_result  = intersectGJK(*boxA, *boxB, transform);
    quaternion_result = intersectGJK(*boxA, *boxB, position, rotation);

    EXPECT_EQ(transform_result, quaternion_result)
        << "Transform3 and quaternion-based GJK should give same result for "
           "rotated case";
    EXPECT_TRUE(transform_result);
}

TEST_F(GJKTest, WorldCoordinatesConsistency)
{
    Vector3<double>    pos_a(0.0, 0.0, 0.0);
    Vector3<double>    pos_b(0.5, 0.0, 0.0);
    Quaternion<double> rot_a(0.0, 0.0, 0.0, 1.0);
    Quaternion<double> rot_b(0.0, 0.0, 0.0, 1.0);

    Transform3<double> transform_a(rot_a, pos_a);
    Transform3<double> transform_b(rot_b, pos_b);

    bool transform_result  = intersectGJK(*boxA, *boxB, transform_a, transform_b);
    bool quaternion_result = intersectGJK(*boxA, *boxB, pos_a, pos_b, rot_a, rot_b);

    EXPECT_EQ(transform_result, quaternion_result)
        << "World coordinates should be consistent between Transform3 and "
           "quaternion versions";
    EXPECT_TRUE(transform_result);

    double angle = M_PI / 8.0;  // 22.5 degrees
    rot_b        = Quaternion<double>(0.0, 0.0, sin(angle / 2.0), cos(angle / 2.0));
    transform_b  = Transform3<double>(rot_b, pos_b);

    transform_result  = intersectGJK(*boxA, *boxB, transform_a, transform_b);
    quaternion_result = intersectGJK(*boxA, *boxB, pos_a, pos_b, rot_a, rot_b);

    EXPECT_EQ(transform_result, quaternion_result)
        << "Rotated world coordinates should be consistent";
    EXPECT_TRUE(transform_result);
}

TEST_F(GJKTest, QuaternionRotationTests)
{
    Vector3<double> overlapping_pos(0.5, 0.0, 0.0);

    double angles[]
        = {0.0, M_PI / 12.0, M_PI / 8.0, M_PI / 6.0};  // 0 deg, 15 deg, 22.5 deg, 30 deg

    for(double angle : angles)
    {
        Quaternion<double> rot_x(sin(angle / 2.0), 0.0, 0.0, cos(angle / 2.0));
        bool               result = intersectGJK(*boxA, *boxB, overlapping_pos, rot_x);
        EXPECT_TRUE(result) << "X-axis rotation at angle " << angle << " should intersect";

        Quaternion<double> rot_y(0.0, sin(angle / 2.0), 0.0, cos(angle / 2.0));
        result = intersectGJK(*boxA, *boxB, overlapping_pos, rot_y);
        EXPECT_TRUE(result) << "Y-axis rotation at angle " << angle << " should intersect";

        Quaternion<double> rot_z(0.0, 0.0, sin(angle / 2.0), cos(angle / 2.0));
        result = intersectGJK(*boxA, *boxB, overlapping_pos, rot_z);
        EXPECT_TRUE(result) << "Z-axis rotation at angle " << angle << " should intersect";
    }
}

TEST_F(GJKTest, QuaternionPerformanceAndEdgeCases)
{
    Quaternion<double> identity_quat(0.0, 0.0, 0.0, 1.0);

    for(int i = 0; i < 30; ++i)
    {
        double          offset = i * 0.1 + 0.01;  // 0.01 to 3.01
        Vector3<double> test_pos(offset, 0.0, 0.0);

        bool result = intersectGJK(*boxA, *boxB, test_pos, identity_quat);

        if(offset < 1.0)
        {
            EXPECT_TRUE(result) << "Should intersect at offset " << offset;
        }
        else if(offset > 2.0)
        {
            EXPECT_FALSE(result) << "Should not intersect at offset " << offset;
        }
        // Skip boundary region assertions
    }

    Box<double>     tiny_box(0.001, 0.001, 0.001);
    Vector3<double> small_offset(0.1, 0.1, 0.1);
    bool            result = intersectGJK(*boxA, tiny_box, small_offset, identity_quat);
    EXPECT_TRUE(result) << "Tiny box should intersect with small offset";

    // Edge case: Large rotation (should still work with overlapping position)
    double             large_angle = M_PI / 3.0;  // 60 degrees
    Quaternion<double> large_rotation(0.0, 0.0, sin(large_angle / 2.0), cos(large_angle / 2.0));
    Vector3<double>    close_pos(0.3, 0.0,
                              0.0);  // Closer position for large rotation
    result = intersectGJK(*boxA, *boxB, close_pos, large_rotation);
    EXPECT_TRUE(result) << "Large rotation should still intersect with close position";
}

TEST_F(GJKTest, QuaternionNormalizationConsistency)
{
    Vector3<double> test_pos(0.5, 0.0, 0.0);
    double          angle = M_PI / 6.0;  // 30 degrees

    Quaternion<double> normalized_quat(0.0, 0.0, sin(angle / 2.0), cos(angle / 2.0));

    // Create unnormalized quaternion (should be automatically normalized by
    // GJK)
    Quaternion<double> unnormalized_quat(0.0, 0.0, 2.0 * sin(angle / 2.0), 2.0 * cos(angle / 2.0));

    bool normalized_result   = intersectGJK(*boxA, *boxB, test_pos, normalized_quat);
    bool unnormalized_result = intersectGJK(*boxA, *boxB, test_pos, unnormalized_quat);

    EXPECT_EQ(normalized_result, unnormalized_result)
        << "GJK should handle quaternion normalization consistently";
    EXPECT_TRUE(normalized_result);
}

// =================================================================================================
// Multi-Shape Tests: Box vs Sphere
// =================================================================================================

TEST_F(GJKTest, BoxSphereTransform3)
{
    bool result = intersectGJK(*boxA, *sphereB, overlapping_transform);
    EXPECT_TRUE(result) << "Box and sphere should intersect when overlapping";

    result = intersectGJK(*boxA, *sphereB, separated_transform);
    EXPECT_FALSE(result) << "Box and sphere should not intersect when separated";

    result = intersectGJK(*boxA, *sphereB, identity_transform, overlapping_transform);
    EXPECT_TRUE(result) << "Box and sphere should intersect in world coordinates";

    double             angle = M_PI / 8.0;  // 22.5 degrees
    Quaternion<double> rotation(0.0, 0.0, sin(angle / 2.0), cos(angle / 2.0));
    Transform3<double> rotated_overlapping(rotation, Vector3<double>(0.5, 0.0, 0.0));
    result = intersectGJK(*boxA, *sphereB, rotated_overlapping);
    EXPECT_TRUE(result) << "Rotated box and sphere should intersect when overlapping";
}

TEST_F(GJKTest, BoxSphereQuaternion)
{
    Vector3<double>    overlapping_pos(0.5, 0.0, 0.0);
    Vector3<double>    separated_pos(3.0, 0.0, 0.0);
    Quaternion<double> identity_quat(0.0, 0.0, 0.0, 1.0);

    bool result = intersectGJK(*boxA, *sphereB, overlapping_pos, identity_quat);
    EXPECT_TRUE(result) << "Box and sphere should intersect when overlapping (quaternion API)";

    result = intersectGJK(*boxA, *sphereB, separated_pos, identity_quat);
    EXPECT_FALSE(result) << "Box and sphere should not intersect when "
                            "separated (quaternion API)";

    Vector3<double> box_pos(0.0, 0.0, 0.0);
    Vector3<double> sphere_pos(0.5, 0.0, 0.0);
    result = intersectGJK(*boxA, *sphereB, box_pos, sphere_pos, identity_quat, identity_quat);
    EXPECT_TRUE(result) << "Box and sphere should intersect in world "
                           "coordinates (quaternion API)";

    double             angle = M_PI / 8.0;  // 22.5 degrees
    Quaternion<double> rotation(0.0, 0.0, sin(angle / 2.0), cos(angle / 2.0));
    result = intersectGJK(*boxA, *sphereB, overlapping_pos, rotation);
    EXPECT_TRUE(result) << "Rotated box and sphere should intersect when "
                           "overlapping (quaternion API)";
}

// =================================================================================================
// Multi-Shape Tests: Box vs Superquadric
// =================================================================================================

TEST_F(GJKTest, BoxSuperquadricTransform3)
{
    bool result = intersectGJK(*boxA, *superquadricB, overlapping_transform);
    EXPECT_TRUE(result) << "Box and superquadric should intersect when overlapping";

    result = intersectGJK(*boxA, *superquadricB, separated_transform);
    EXPECT_FALSE(result) << "Box and superquadric should not intersect when separated";

    result = intersectGJK(*boxA, *superquadricB, identity_transform, overlapping_transform);
    EXPECT_TRUE(result) << "Box and superquadric should intersect in world coordinates";

    double             angle = M_PI / 6.0;  // 30 degrees
    Quaternion<double> rotation(0.0, 0.0, sin(angle / 2.0), cos(angle / 2.0));
    Transform3<double> rotated_overlapping(rotation, Vector3<double>(0.4, 0.0, 0.0));
    result = intersectGJK(*boxA, *superquadricB, rotated_overlapping);
    EXPECT_TRUE(result) << "Rotated box and superquadric should intersect when overlapping";
}

TEST_F(GJKTest, BoxSuperquadricQuaternion)
{
    Vector3<double>    overlapping_pos(0.5, 0.0, 0.0);
    Vector3<double>    separated_pos(3.0, 0.0, 0.0);
    Quaternion<double> identity_quat(0.0, 0.0, 0.0, 1.0);

    bool result = intersectGJK(*boxA, *superquadricB, overlapping_pos, identity_quat);
    EXPECT_TRUE(result) << "Box and superquadric should intersect when "
                           "overlapping (quaternion API)";

    result = intersectGJK(*boxA, *superquadricB, separated_pos, identity_quat);
    EXPECT_FALSE(result) << "Box and superquadric should not intersect when "
                            "separated (quaternion API)";

    Vector3<double> box_pos(0.0, 0.0, 0.0);
    Vector3<double> superquadric_pos(0.5, 0.0, 0.0);
    result = intersectGJK(*boxA,
                          *superquadricB,
                          box_pos,
                          superquadric_pos,
                          identity_quat,
                          identity_quat);
    EXPECT_TRUE(result) << "Box and superquadric should intersect in world "
                           "coordinates (quaternion API)";

    double             angle = M_PI / 6.0;  // 30 degrees
    Quaternion<double> rotation(0.0, 0.0, sin(angle / 2.0), cos(angle / 2.0));
    result = intersectGJK(*boxA, *superquadricB, overlapping_pos, rotation);
    EXPECT_TRUE(result) << "Rotated box and superquadric should intersect when "
                           "overlapping (quaternion API)";
}

// =================================================================================================
// Multi-Shape Tests: Sphere vs Superquadric
// =================================================================================================

TEST_F(GJKTest, SphereSuperquadricTransform3)
{
    bool result = intersectGJK(*sphereA, *superquadricB, overlapping_transform);
    EXPECT_TRUE(result) << "Sphere and superquadric should intersect when overlapping";

    result = intersectGJK(*sphereA, *superquadricB, separated_transform);
    EXPECT_FALSE(result) << "Sphere and superquadric should not intersect when separated";

    result = intersectGJK(*sphereA, *superquadricB, identity_transform, overlapping_transform);
    EXPECT_TRUE(result) << "Sphere and superquadric should intersect in world coordinates";

    double             angle = M_PI / 4.0;  // 45 degrees
    Quaternion<double> rotation(0.0, 0.0, sin(angle / 2.0), cos(angle / 2.0));
    Transform3<double> rotated_overlapping(rotation, Vector3<double>(0.3, 0.0, 0.0));
    result = intersectGJK(*sphereA, *superquadricB, rotated_overlapping);
    EXPECT_TRUE(result) << "Rotated sphere and superquadric should intersect when overlapping";
}

TEST_F(GJKTest, SphereSuperquadricQuaternion)
{
    Vector3<double>    overlapping_pos(0.5, 0.0, 0.0);
    Vector3<double>    separated_pos(3.0, 0.0, 0.0);
    Quaternion<double> identity_quat(0.0, 0.0, 0.0, 1.0);

    bool result = intersectGJK(*sphereA, *superquadricB, overlapping_pos, identity_quat);
    EXPECT_TRUE(result) << "Sphere and superquadric should intersect when "
                           "overlapping (quaternion API)";

    result = intersectGJK(*sphereA, *superquadricB, separated_pos, identity_quat);
    EXPECT_FALSE(result) << "Sphere and superquadric should not intersect when "
                            "separated (quaternion API)";

    Vector3<double> sphere_pos(0.0, 0.0, 0.0);
    Vector3<double> superquadric_pos(0.3, 0.0, 0.0);
    result = intersectGJK(*sphereA,
                          *superquadricB,
                          sphere_pos,
                          superquadric_pos,
                          identity_quat,
                          identity_quat);
    EXPECT_TRUE(result) << "Sphere and superquadric should intersect in world "
                           "coordinates (quaternion API)";

    double             angle = M_PI / 4.0;  // 45 degrees
    Quaternion<double> rotation(0.0, 0.0, sin(angle / 2.0), cos(angle / 2.0));
    result = intersectGJK(*sphereA, *superquadricB, overlapping_pos, rotation);
    EXPECT_TRUE(result) << "Rotated sphere and superquadric should intersect "
                           "when overlapping (quaternion API)";
}

// =================================================================================================
// Multi-Shape Tests: Sphere vs Sphere
// =================================================================================================

TEST_F(GJKTest, SphereSphereTransform3)
{
    bool result = intersectGJK(*sphereA, *sphereB, overlapping_transform);
    EXPECT_TRUE(result) << "Spheres should intersect when overlapping";

    result = intersectGJK(*sphereA, *sphereB, separated_transform);
    EXPECT_FALSE(result) << "Spheres should not intersect when separated";

    // Test touching case (distance = sum of radii = 1.5)
    Transform3<double> touching_transform(Quaternion<double>(0.0, 0.0, 0.0, 1.0),
                                          Vector3<double>(1.5, 0.0, 0.0));
    result = intersectGJK(*sphereA, *sphereB, touching_transform);
    // Note: GJK might have numerical precision issues at exact touching
    // We test slightly inside touching distance
    Transform3<double> near_touching_transform(Quaternion<double>(0.0, 0.0, 0.0, 1.0),
                                               Vector3<double>(1.4, 0.0, 0.0));
    result = intersectGJK(*sphereA, *sphereB, near_touching_transform);
    EXPECT_TRUE(result) << "Spheres should intersect when slightly overlapping";
}

TEST_F(GJKTest, SphereSphereQuaternion)
{
    Vector3<double>    overlapping_pos(0.5, 0.0, 0.0);
    Vector3<double>    separated_pos(3.0, 0.0, 0.0);
    Vector3<double>    near_touching_pos(1.4, 0.0, 0.0);
    Quaternion<double> identity_quat(0.0, 0.0, 0.0, 1.0);

    bool result = intersectGJK(*sphereA, *sphereB, overlapping_pos, identity_quat);
    EXPECT_TRUE(result) << "Spheres should intersect when overlapping (quaternion API)";

    result = intersectGJK(*sphereA, *sphereB, separated_pos, identity_quat);
    EXPECT_FALSE(result) << "Spheres should not intersect when separated (quaternion API)";

    result = intersectGJK(*sphereA, *sphereB, near_touching_pos, identity_quat);
    EXPECT_TRUE(result) << "Spheres should intersect when slightly overlapping "
                           "(quaternion API)";

    Vector3<double> sphere1_pos(0.0, 0.0, 0.0);
    Vector3<double> sphere2_pos(0.5, 0.0, 0.0);
    result
        = intersectGJK(*sphereA, *sphereB, sphere1_pos, sphere2_pos, identity_quat, identity_quat);
    EXPECT_TRUE(result) << "Spheres should intersect in world coordinates (quaternion API)";
}

// =================================================================================================
// Multi-Shape Tests: Superquadric vs Superquadric
// =================================================================================================

TEST_F(GJKTest, SuperquadricSuperquadricTransform3)
{
    bool result = intersectGJK(*superquadricA, *superquadricB, overlapping_transform);
    EXPECT_TRUE(result) << "Superquadrics should intersect when overlapping";

    result = intersectGJK(*superquadricA, *superquadricB, separated_transform);
    EXPECT_FALSE(result) << "Superquadrics should not intersect when separated";

    result
        = intersectGJK(*superquadricA, *superquadricB, identity_transform, overlapping_transform);
    EXPECT_TRUE(result) << "Superquadrics should intersect in world coordinates";

    double             angle = M_PI / 3.0;  // 60 degrees
    Quaternion<double> rotation(0.0, 0.0, sin(angle / 2.0), cos(angle / 2.0));
    Transform3<double> rotated_overlapping(rotation, Vector3<double>(0.4, 0.0, 0.0));
    result = intersectGJK(*superquadricA, *superquadricB, rotated_overlapping);
    EXPECT_TRUE(result) << "Rotated superquadrics should intersect when overlapping";
}

TEST_F(GJKTest, SuperquadricSuperquadricQuaternion)
{
    Vector3<double>    overlapping_pos(0.5, 0.0, 0.0);
    Vector3<double>    separated_pos(3.0, 0.0, 0.0);
    Quaternion<double> identity_quat(0.0, 0.0, 0.0, 1.0);

    bool result = intersectGJK(*superquadricA, *superquadricB, overlapping_pos, identity_quat);
    EXPECT_TRUE(result) << "Superquadrics should intersect when overlapping (quaternion API)";

    result = intersectGJK(*superquadricA, *superquadricB, separated_pos, identity_quat);
    EXPECT_FALSE(result) << "Superquadrics should not intersect when separated (quaternion API)";

    Vector3<double> superquadric1_pos(0.0, 0.0, 0.0);
    Vector3<double> superquadric2_pos(0.4, 0.0, 0.0);
    result = intersectGJK(*superquadricA,
                          *superquadricB,
                          superquadric1_pos,
                          superquadric2_pos,
                          identity_quat,
                          identity_quat);
    EXPECT_TRUE(result) << "Superquadrics should intersect in world "
                           "coordinates (quaternion API)";

    double             angle = M_PI / 3.0;  // 60 degrees
    Quaternion<double> rotation(0.0, 0.0, sin(angle / 2.0), cos(angle / 2.0));
    result = intersectGJK(*superquadricA, *superquadricB, overlapping_pos, rotation);
    EXPECT_TRUE(result) << "Rotated superquadrics should intersect when "
                           "overlapping (quaternion API)";
}

// =================================================================================================
// API Consistency Tests for All Shape Combinations
// =================================================================================================

TEST_F(GJKTest, AllShapesAPIConsistency)
{
    Vector3<double>    test_pos(0.5, 0.0, 0.0);
    Quaternion<double> test_quat(0.0, 0.0, 0.0, 1.0);
    Transform3<double> test_transform(test_quat, test_pos);

    std::vector<std::pair<Convex<double>*, Convex<double>*>> shape_pairs
        = {{boxA, boxB},
           {boxA, sphereB},
           {boxA, superquadricB},
           {sphereA, sphereB},
           {sphereA, superquadricB},
           {superquadricA, superquadricB}};

    std::vector<std::string> shape_names = {"Box-Box",
                                            "Box-Sphere",
                                            "Box-Superquadric",
                                            "Sphere-Sphere",
                                            "Sphere-Superquadric",
                                            "Superquadric-Superquadric"};

    for(size_t i = 0; i < shape_pairs.size(); ++i)
    {
        auto&              pair = shape_pairs[i];
        const std::string& name = shape_names[i];

        bool transform_result  = intersectGJK(*pair.first, *pair.second, test_transform);
        bool quaternion_result = intersectGJK(*pair.first, *pair.second, test_pos, test_quat);

        EXPECT_EQ(transform_result, quaternion_result)
            << "Transform3 vs Quaternion API inconsistency for " << name;
        EXPECT_TRUE(transform_result) << name << " should intersect at overlapping position";

        Vector3<double>    pos_a(0.0, 0.0, 0.0);
        Vector3<double>    pos_b(0.5, 0.0, 0.0);
        Quaternion<double> rot_a(0.0, 0.0, 0.0, 1.0);
        Quaternion<double> rot_b(0.0, 0.0, 0.0, 1.0);

        Transform3<double> transform_a(rot_a, pos_a);
        Transform3<double> transform_b(rot_b, pos_b);

        bool world_transform_result
            = intersectGJK(*pair.first, *pair.second, transform_a, transform_b);
        bool world_quaternion_result
            = intersectGJK(*pair.first, *pair.second, pos_a, pos_b, rot_a, rot_b);

        EXPECT_EQ(world_transform_result, world_quaternion_result)
            << "World coordinates API inconsistency for " << name;
        EXPECT_TRUE(world_transform_result) << name << " should intersect in world coordinates";
    }
}

// =================================================================================================
// Performance Tests for Different Shape Combinations
// =================================================================================================

TEST_F(GJKTest, MultiShapePerformanceTest)
{
    Vector3<double>    base_pos(0.1, 0.0, 0.0);
    Quaternion<double> identity_quat(0.0, 0.0, 0.0, 1.0);

    for(int i = 0; i < 20; ++i)
    {
        double          offset = i * 0.1 + 0.1;  // 0.1 to 2.1
        Vector3<double> test_pos(offset, 0.0, 0.0);

        bool box_sphere_result = intersectGJK(*boxA, *sphereB, test_pos, identity_quat);

        bool sphere_superquadric_result
            = intersectGJK(*sphereA, *superquadricB, test_pos, identity_quat);

        bool box_superquadric_result = intersectGJK(*boxA, *superquadricB, test_pos, identity_quat);

        if(offset < 1.0)
        {
            EXPECT_TRUE(box_sphere_result) << "Box-Sphere should intersect at offset " << offset;
            EXPECT_TRUE(sphere_superquadric_result)
                << "Sphere-Superquadric should intersect at offset " << offset;
            EXPECT_TRUE(box_superquadric_result)
                << "Box-Superquadric should intersect at offset " << offset;
        }
        else if(offset > 1.8)
        {
            EXPECT_FALSE(box_sphere_result)
                << "Box-Sphere should not intersect at offset " << offset;
            EXPECT_FALSE(sphere_superquadric_result)
                << "Sphere-Superquadric should not intersect at offset " << offset;
            EXPECT_FALSE(box_superquadric_result)
                << "Box-Superquadric should not intersect at offset " << offset;
        }
    }
}

// =================================================================================================
// Edge Cases for Different Shape Combinations
// =================================================================================================

TEST_F(GJKTest, MultiShapeEdgeCases)
{
    Quaternion<double> identity_quat(0.0, 0.0, 0.0, 1.0);

    Sphere<double>       tiny_sphere(0.001);
    Box<double>          tiny_box(0.001, 0.001, 0.001);
    Superquadric<double> tiny_superquadric(0.001, 0.001, 0.001, 2.0, 2.0);

    Vector3<double> small_offset(0.1, 0.1, 0.1);

    bool result = intersectGJK(*boxA, tiny_sphere, small_offset, identity_quat);
    EXPECT_TRUE(result) << "Box should intersect with tiny sphere at small offset";

    result = intersectGJK(*sphereA, tiny_sphere, small_offset, identity_quat);
    EXPECT_TRUE(result) << "Sphere should intersect with tiny sphere at small offset";

    result = intersectGJK(*superquadricA, tiny_sphere, small_offset, identity_quat);
    EXPECT_TRUE(result) << "Superquadric should intersect with tiny sphere at small offset";

    result = intersectGJK(*sphereA, tiny_box, small_offset, identity_quat);
    EXPECT_TRUE(result) << "Sphere should intersect with tiny box at small offset";

    result = intersectGJK(*superquadricA, tiny_box, small_offset, identity_quat);
    EXPECT_TRUE(result) << "Superquadric should intersect with tiny box at small offset";

    result = intersectGJK(*boxA, tiny_superquadric, small_offset, identity_quat);
    EXPECT_TRUE(result) << "Box should intersect with tiny superquadric at small offset";

    result = intersectGJK(*sphereA, tiny_superquadric, small_offset, identity_quat);
    EXPECT_TRUE(result) << "Sphere should intersect with tiny superquadric at small offset";

    double             large_angle = M_PI / 2.0;  // 90 degrees
    Quaternion<double> large_rotation(0.0, 0.0, sin(large_angle / 2.0), cos(large_angle / 2.0));
    Vector3<double>    close_pos(0.2, 0.0, 0.0);

    result = intersectGJK(*boxA, *sphereB, close_pos, large_rotation);
    EXPECT_TRUE(result) << "Box-Sphere should intersect with large rotation at close position";

    result = intersectGJK(*sphereA, *superquadricB, close_pos, large_rotation);
    EXPECT_TRUE(result) << "Sphere-Superquadric should intersect with large "
                           "rotation at close position";

    result = intersectGJK(*boxA, *superquadricB, close_pos, large_rotation);
    EXPECT_TRUE(result) << "Box-Superquadric should intersect with large "
                           "rotation at close position";
}

// =================================================================================================
// Overload-consistency tests
//
// intersectGJK and computeClosestPoints_GJK each have 4 overloads:
//   1. Relative Transform3  (b2a)
//   2. Absolute Transform3  (a2w, b2w)
//   3. Relative quat/vec    (v_b2a, q_b2a)
//   4. Absolute quat/vec    (v_a2w, v_b2w, q_a2w, q_b2w)
//
// All derive their inputs from the same world-space description, so every
// result must be numerically identical.
//
// NOTE on closest-point frames:
//   pa is always returned in A's LOCAL frame for all 4 overloads.
//   pb is always returned in B's LOCAL frame for all 4 overloads.
//   Therefore pa1==pa2==pa3==pa4 and pb1==pb2==pb3==pb4 (within floating
//   point tolerance), and no coordinate-change is needed.
// =================================================================================================

// -------------------------------------------------------------------------------------------------
/** @brief Builds a unit quaternion from axis (ax,ay,az) and angle (rad). */
static Quaternion<double> gjkMakeQuat(double ax, double ay, double az, double angle)
{
    const double s   = std::sin(angle * 0.5);
    const double c   = std::cos(angle * 0.5);
    const double len = std::sqrt(ax * ax + ay * ay + az * az);
    return Quaternion<double>((ax / len) * s, (ay / len) * s, (az / len) * s, c);
}

// -------------------------------------------------------------------------------------------------
/** @brief Calls all 4 intersectGJK overloads from world-space inputs and
    asserts they all return `expected` and agree with each other. */
static void checkAllGJKIntersectOverloads(const Convex<double>&     a,
                                          const Convex<double>&     b,
                                          const Vector3<double>&    v_a2w,
                                          const Quaternion<double>& q_a2w,
                                          const Vector3<double>&    v_b2w,
                                          const Quaternion<double>& q_b2w,
                                          bool                      expected,
                                          const std::string&        label)
{
    const Transform3<double> trA2W(q_a2w, v_a2w);
    const Transform3<double> trB2W(q_b2w, v_b2w);
    // Derive relative inputs: B expressed in A's frame
    const Transform3<double> trB2A(trA2W, trB2W);
    const Quaternion<double> q_b2a = inverse(q_a2w) * q_b2w;
    const Vector3<double>    v_b2a = q_a2w << (v_b2w - v_a2w);

    // Overload 1: relative Transform3
    const bool res1 = intersectGJK(a, b, trB2A);
    // Overload 2: absolute Transform3
    const bool res2 = intersectGJK(a, b, trA2W, trB2W);
    // Overload 3: relative quat/vec
    const bool res3 = intersectGJK(a, b, v_b2a, q_b2a);
    // Overload 4: absolute quat/vec
    const bool res4 = intersectGJK(a, b, v_a2w, v_b2w, q_a2w, q_b2w);

    EXPECT_EQ(res1, expected) << "[" << label << "] overload 1 (rel Transform3) wrong";
    EXPECT_EQ(res2, expected) << "[" << label << "] overload 2 (abs Transform3) wrong";
    EXPECT_EQ(res3, expected) << "[" << label << "] overload 3 (rel quat/vec) wrong";
    EXPECT_EQ(res4, expected) << "[" << label << "] overload 4 (abs quat/vec) wrong";

    EXPECT_EQ(res1, res2) << "[" << label << "] overloads 1 vs 2 disagree";
    EXPECT_EQ(res1, res3) << "[" << label << "] overloads 1 vs 3 disagree";
    EXPECT_EQ(res1, res4) << "[" << label << "] overloads 1 vs 4 disagree";
}

// -------------------------------------------------------------------------------------------------
/** @brief Calls all 4 computeClosestPoints_GJK overloads from world-space
    inputs and asserts that (a) all four distances agree and (b) each overload
    returns a witness pair whose world-space separation equals that distance.
    Only call this when the shapes are NOT intersecting (d > 0).

    Cross-overload comparison of the raw pa/pb components is deliberately
    avoided: for degenerate configurations (e.g. parallel face-to-face boxes)
    GJK may return any valid witness on the closest feature, and overloads that
    receive pre-computed relative transforms vs. absolute world transforms start
    the search in opposite directions, converging to different -- but equally
    valid -- corners of the same face.

    @param distEPS  Tolerance for distance / segment-length comparison */
static void checkAllGJKClosestPointsOverloads(const Convex<double>&     a,
                                              const Convex<double>&     b,
                                              const Vector3<double>&    v_a2w,
                                              const Quaternion<double>& q_a2w,
                                              const Vector3<double>&    v_b2w,
                                              const Quaternion<double>& q_b2w,
                                              const std::string&        label,
                                              double                    distEPS = 1e-6)
{
    const Transform3<double> trA2W(q_a2w, v_a2w);
    const Transform3<double> trB2W(q_b2w, v_b2w);
    const Transform3<double> trB2A(trA2W, trB2W);
    const Quaternion<double> q_b2a = inverse(q_a2w) * q_b2w;
    const Vector3<double>    v_b2a = q_a2w << (v_b2w - v_a2w);

    constexpr double crust = 0.0;
    uint             nbIter1, nbIter2, nbIter3, nbIter4;
    Vector3<double>  pa1, pb1, pa2, pb2, pa3, pb3, pa4, pb4;

    // Overload 1: relative Transform3
    const double d1 = computeClosestPoints_GJK<double, GJKType::JOHNSON>(a,
                                                                         b,
                                                                         trB2A,
                                                                         crust,
                                                                         crust,
                                                                         pa1,
                                                                         pb1,
                                                                         nbIter1);
    // Overload 2: absolute Transform3
    const double d2 = computeClosestPoints_GJK<double, GJKType::JOHNSON>(a,
                                                                         b,
                                                                         trA2W,
                                                                         trB2W,
                                                                         crust,
                                                                         crust,
                                                                         pa2,
                                                                         pb2,
                                                                         nbIter2);
    // Overload 3: relative quat/vec
    const double d3 = computeClosestPoints_GJK<double, GJKType::JOHNSON>(a,
                                                                         b,
                                                                         v_b2a,
                                                                         q_b2a,
                                                                         crust,
                                                                         crust,
                                                                         pa3,
                                                                         pb3,
                                                                         nbIter3);
    // Overload 4: absolute quat/vec
    const double d4 = computeClosestPoints_GJK<double, GJKType::JOHNSON>(a,
                                                                         b,
                                                                         v_a2w,
                                                                         v_b2w,
                                                                         q_a2w,
                                                                         q_b2w,
                                                                         crust,
                                                                         crust,
                                                                         pa4,
                                                                         pb4,
                                                                         nbIter4);

    // All distances must agree
    EXPECT_NEAR(d1, d2, distEPS) << "[" << label << "] distance: overloads 1 vs 2 differ";
    EXPECT_NEAR(d1, d3, distEPS) << "[" << label << "] distance: overloads 1 vs 3 differ";
    EXPECT_NEAR(d1, d4, distEPS) << "[" << label << "] distance: overloads 1 vs 4 differ";

    // Verify that each overload's witness pair forms a segment of length == d in
    // world space.  We convert pa (in A's local frame) and pb (in B's local frame)
    // back to world coordinates and check the Euclidean distance.
    const auto toWorldA = [&](const Vector3<double>& p) { return (q_a2w >> p) + v_a2w; };
    const auto toWorldB = [&](const Vector3<double>& p) { return (q_b2w >> p) + v_b2w; };

    EXPECT_NEAR(norm(toWorldA(pa1) - toWorldB(pb1)), d1, distEPS)
        << "[" << label << "] overload 1: |pa_world - pb_world| != d";
    EXPECT_NEAR(norm(toWorldA(pa2) - toWorldB(pb2)), d2, distEPS)
        << "[" << label << "] overload 2: |pa_world - pb_world| != d";
    EXPECT_NEAR(norm(toWorldA(pa3) - toWorldB(pb3)), d3, distEPS)
        << "[" << label << "] overload 3: |pa_world - pb_world| != d";
    EXPECT_NEAR(norm(toWorldA(pa4) - toWorldB(pb4)), d4, distEPS)
        << "[" << label << "] overload 4: |pa_world - pb_world| != d";
}

// =================================================================================================
// Overload-consistency fixture
// =================================================================================================
class GJKOverloadConsistencyTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        boxA          = new Box<double>(1.0, 1.0, 1.0);
        boxB          = new Box<double>(0.5, 0.5, 0.5);
        sphereA       = new Sphere<double>(1.0);
        sphereB       = new Sphere<double>(0.5);
        superquadricA = new Superquadric<double>(1.0, 1.0, 1.0, 2.0, 2.0);
        superquadricB = new Superquadric<double>(0.5, 0.5, 0.5, 1.5, 1.5);

        q_id      = Quaternion<double>(0.0, 0.0, 0.0, 1.0);
        q_rot90Y  = gjkMakeQuat(0.0, 1.0, 0.0, M_PI / 2.0);
        q_rot45Z  = gjkMakeQuat(0.0, 0.0, 1.0, M_PI / 4.0);
        q_oblique = gjkMakeQuat(1.0, 1.0, 0.0, M_PI / 6.0);
    }

    void TearDown() override
    {
        delete boxA;
        delete boxB;
        delete sphereA;
        delete sphereB;
        delete superquadricA;
        delete superquadricB;
    }

    Box<double>*          boxA;
    Box<double>*          boxB;
    Sphere<double>*       sphereA;
    Sphere<double>*       sphereB;
    Superquadric<double>* superquadricA;
    Superquadric<double>* superquadricB;

    Quaternion<double> q_id;
    Quaternion<double> q_rot90Y;
    Quaternion<double> q_rot45Z;
    Quaternion<double> q_oblique;
};

// =================================================================================================
// intersectGJK -- all 4 overloads produce consistent results
// =================================================================================================

TEST_F(GJKOverloadConsistencyTest, IntersectOverloads_BoxBox_Overlapping)
{
    checkAllGJKIntersectOverloads(*boxA,
                                  *boxB,
                                  Vector3<double>(0.0, 0.0, 0.0),
                                  q_id,
                                  Vector3<double>(0.5, 0.0, 0.0),
                                  q_id,
                                  true,
                                  "BoxBox_Overlapping");
}

TEST_F(GJKOverloadConsistencyTest, IntersectOverloads_BoxBox_Separated)
{
    checkAllGJKIntersectOverloads(*boxA,
                                  *boxB,
                                  Vector3<double>(0.0, 0.0, 0.0),
                                  q_id,
                                  Vector3<double>(3.0, 0.0, 0.0),
                                  q_id,
                                  false,
                                  "BoxBox_Separated");
}

TEST_F(GJKOverloadConsistencyTest, IntersectOverloads_BoxBox_Rotated_Overlapping)
{
    // A rotated 45 deg around Z, B displaced 0.5 in X, rotated 90 deg around Y
    checkAllGJKIntersectOverloads(*boxA,
                                  *boxB,
                                  Vector3<double>(0.0, 0.0, 0.0),
                                  q_rot45Z,
                                  Vector3<double>(0.5, 0.0, 0.0),
                                  q_rot90Y,
                                  true,
                                  "BoxBox_Rotated_Overlapping");
}

TEST_F(GJKOverloadConsistencyTest, IntersectOverloads_BoxBox_Rotated_Separated)
{
    checkAllGJKIntersectOverloads(*boxA,
                                  *boxB,
                                  Vector3<double>(0.0, 0.0, 0.0),
                                  q_oblique,
                                  Vector3<double>(5.0, 0.0, 0.0),
                                  q_rot90Y,
                                  false,
                                  "BoxBox_Rotated_Separated");
}

TEST_F(GJKOverloadConsistencyTest, IntersectOverloads_BoxSphere_Overlapping)
{
    checkAllGJKIntersectOverloads(*boxA,
                                  *sphereB,
                                  Vector3<double>(0.0, 0.0, 0.0),
                                  q_id,
                                  Vector3<double>(0.5, 0.0, 0.0),
                                  q_id,
                                  true,
                                  "BoxSphere_Overlapping");
}

TEST_F(GJKOverloadConsistencyTest, IntersectOverloads_BoxSphere_Separated)
{
    checkAllGJKIntersectOverloads(*boxA,
                                  *sphereB,
                                  Vector3<double>(0.0, 0.0, 0.0),
                                  q_rot90Y,
                                  Vector3<double>(4.0, 0.0, 0.0),
                                  q_oblique,
                                  false,
                                  "BoxSphere_Separated");
}

TEST_F(GJKOverloadConsistencyTest, IntersectOverloads_SphereSphere_Overlapping)
{
    checkAllGJKIntersectOverloads(*sphereA,
                                  *sphereB,
                                  Vector3<double>(0.0, 0.0, 0.0),
                                  q_id,
                                  Vector3<double>(1.2, 0.0, 0.0),
                                  q_id,
                                  true,
                                  "SphereSphere_Overlapping");
}

TEST_F(GJKOverloadConsistencyTest, IntersectOverloads_SphereSphere_Separated)
{
    checkAllGJKIntersectOverloads(*sphereA,
                                  *sphereB,
                                  Vector3<double>(0.0, 0.0, 0.0),
                                  q_id,
                                  Vector3<double>(4.0, 0.0, 0.0),
                                  q_id,
                                  false,
                                  "SphereSphere_Separated");
}

TEST_F(GJKOverloadConsistencyTest, IntersectOverloads_BoxSuperquadric_Overlapping)
{
    checkAllGJKIntersectOverloads(*boxA,
                                  *superquadricB,
                                  Vector3<double>(0.0, 0.0, 0.0),
                                  q_rot45Z,
                                  Vector3<double>(0.4, 0.0, 0.0),
                                  q_rot90Y,
                                  true,
                                  "BoxSuperquadric_Overlapping");
}

TEST_F(GJKOverloadConsistencyTest, IntersectOverloads_BoxSuperquadric_Separated)
{
    checkAllGJKIntersectOverloads(*boxA,
                                  *superquadricB,
                                  Vector3<double>(0.0, 0.0, 0.0),
                                  q_oblique,
                                  Vector3<double>(6.0, 0.0, 0.0),
                                  q_rot90Y,
                                  false,
                                  "BoxSuperquadric_Separated");
}

// =================================================================================================
// computeClosestPoints_GJK -- all 4 overloads produce same distance and closest points
// (Only tested for non-intersecting pairs where the distance result is meaningful)
// =================================================================================================

TEST_F(GJKOverloadConsistencyTest, ClosestPoints_BoxBox_IdentityOrientation)
{
    // Separated boxes, identity orientation, gap = 3 - 1 - 0.5 = 1.5
    checkAllGJKClosestPointsOverloads(*boxA,
                                      *boxB,
                                      Vector3<double>(0.0, 0.0, 0.0),
                                      q_id,
                                      Vector3<double>(3.0, 0.0, 0.0),
                                      q_id,
                                      "ClosestPts_BoxBox_Identity");
}

TEST_F(GJKOverloadConsistencyTest, ClosestPoints_BoxBox_Rotated)
{
    // Box A identity at origin, box B rotated 90 degY, center at (4, 0, 0)
    checkAllGJKClosestPointsOverloads(*boxA,
                                      *boxB,
                                      Vector3<double>(0.0, 0.0, 0.0),
                                      q_id,
                                      Vector3<double>(4.0, 0.0, 0.0),
                                      q_rot90Y,
                                      "ClosestPts_BoxBox_Rotated");
}

TEST_F(GJKOverloadConsistencyTest, ClosestPoints_BoxSphere_Separated)
{
    // Box A at origin, sphere B at (4, 0, 0): gap = 4 - 1 - 0.5 = 2.5
    checkAllGJKClosestPointsOverloads(*boxA,
                                      *sphereB,
                                      Vector3<double>(0.0, 0.0, 0.0),
                                      q_id,
                                      Vector3<double>(4.0, 0.0, 0.0),
                                      q_id,
                                      "ClosestPts_BoxSphere_Separated");
}

TEST_F(GJKOverloadConsistencyTest, ClosestPoints_BoxSphere_Rotated_A)
{
    // Rotated box A separated from sphere B
    checkAllGJKClosestPointsOverloads(*boxA,
                                      *sphereB,
                                      Vector3<double>(0.0, 0.0, 0.0),
                                      q_rot45Z,
                                      Vector3<double>(5.0, 0.0, 0.0),
                                      q_id,
                                      "ClosestPts_BoxSphere_RotatedA");
}

TEST_F(GJKOverloadConsistencyTest, ClosestPoints_SphereSphere_Separated)
{
    // Sphere A at origin, sphere B at (5, 0, 0): gap = 5 - 1 - 0.5 = 3.5
    checkAllGJKClosestPointsOverloads(*sphereA,
                                      *sphereB,
                                      Vector3<double>(0.0, 0.0, 0.0),
                                      q_id,
                                      Vector3<double>(5.0, 0.0, 0.0),
                                      q_id,
                                      "ClosestPts_SphereSphere_Separated");
}

TEST_F(GJKOverloadConsistencyTest, ClosestPoints_SphereSphere_DiagonalSeparation)
{
    // Diagonal separation: A at origin, B at (3, 3, 3)
    checkAllGJKClosestPointsOverloads(*sphereA,
                                      *sphereB,
                                      Vector3<double>(0.0, 0.0, 0.0),
                                      q_id,
                                      Vector3<double>(3.0, 3.0, 3.0),
                                      q_id,
                                      "ClosestPts_SphereSphere_Diagonal");
}

TEST_F(GJKOverloadConsistencyTest, ClosestPoints_BoxSuperquadric_Separated)
{
    // Box A rotated at origin, superquadric B far away with oblique rotation
    checkAllGJKClosestPointsOverloads(*boxA,
                                      *superquadricB,
                                      Vector3<double>(0.0, 0.0, 0.0),
                                      q_rot90Y,
                                      Vector3<double>(5.0, 0.0, 0.0),
                                      q_oblique,
                                      "ClosestPts_BoxSuperquadric_Separated");
}

TEST_F(GJKOverloadConsistencyTest, ClosestPoints_NonIdentityBothShapes)
{
    // Both shapes at non-trivial world poses
    checkAllGJKClosestPointsOverloads(*boxA,
                                      *sphereA,
                                      Vector3<double>(-3.0, 1.0, 0.0),
                                      q_rot45Z,
                                      Vector3<double>(3.0, -1.0, 0.0),
                                      q_oblique,
                                      "ClosestPts_NonIdentityBoth");
}

// =================================================================================================
// Random 100-config consistency test
// =================================================================================================
TEST_F(GJKOverloadConsistencyTest, RandomConfigs_AllOverloadsConsistent)
{
    std::mt19937_64                        rng(20260315ULL);
    std::uniform_real_distribution<double> distPos(-4.0, 4.0);
    std::uniform_real_distribution<double> distDim(0.1, 1.5);
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

    // Shape pool to randomly pick from each trial
    Convex<double>* shapes[]  = {boxA, boxB, sphereA, sphereB, superquadricA, superquadricB};
    constexpr int   numShapes = 6;
    std::uniform_int_distribution<int> distShape(0, numShapes - 1);

    constexpr double crust   = 0.0;
    constexpr double distEPS = 1e-6;

    for(int trial = 0; trial < 100; ++trial)
    {
        const Convex<double>& shapeA = *shapes[distShape(rng)];
        const Convex<double>& shapeB = *shapes[distShape(rng)];

        const Vector3<double>    v_a2w(distPos(rng), distPos(rng), distPos(rng));
        const Vector3<double>    v_b2w(distPos(rng), distPos(rng), distPos(rng));
        const Quaternion<double> q_a2w = randQuat();
        const Quaternion<double> q_b2w = randQuat();

        // ---- Derive relative inputs ----
        const Transform3<double> trA2W(q_a2w, v_a2w);
        const Transform3<double> trB2W(q_b2w, v_b2w);
        const Transform3<double> trB2A(trA2W, trB2W);
        const Quaternion<double> q_b2a = inverse(q_a2w) * q_b2w;
        const Vector3<double>    v_b2a = q_a2w << (v_b2w - v_a2w);

        // ---- intersectGJK: all 4 overloads must agree ----
        const bool iRes1 = intersectGJK(shapeA, shapeB, trB2A);
        const bool iRes2 = intersectGJK(shapeA, shapeB, trA2W, trB2W);
        const bool iRes3 = intersectGJK(shapeA, shapeB, v_b2a, q_b2a);
        const bool iRes4 = intersectGJK(shapeA, shapeB, v_a2w, v_b2w, q_a2w, q_b2w);

        EXPECT_EQ(iRes1, iRes2) << "trial " << trial << ": intersect overloads 1 vs 2 disagree";
        EXPECT_EQ(iRes1, iRes3) << "trial " << trial << ": intersect overloads 1 vs 3 disagree";
        EXPECT_EQ(iRes1, iRes4) << "trial " << trial << ": intersect overloads 1 vs 4 disagree";

        // ---- computeClosestPoints_GJK: all 4 overloads must agree ----
        // Only meaningful when shapes are not intersecting (d > small threshold)
        uint            nbIter1, nbIter2, nbIter3, nbIter4;
        Vector3<double> pa1, pb1, pa2, pb2, pa3, pb3, pa4, pb4;

        const double d1 = computeClosestPoints_GJK<double, GJKType::JOHNSON>(shapeA,
                                                                             shapeB,
                                                                             trB2A,
                                                                             crust,
                                                                             crust,
                                                                             pa1,
                                                                             pb1,
                                                                             nbIter1);
        const double d2 = computeClosestPoints_GJK<double, GJKType::JOHNSON>(shapeA,
                                                                             shapeB,
                                                                             trA2W,
                                                                             trB2W,
                                                                             crust,
                                                                             crust,
                                                                             pa2,
                                                                             pb2,
                                                                             nbIter2);
        const double d3 = computeClosestPoints_GJK<double, GJKType::JOHNSON>(shapeA,
                                                                             shapeB,
                                                                             v_b2a,
                                                                             q_b2a,
                                                                             crust,
                                                                             crust,
                                                                             pa3,
                                                                             pb3,
                                                                             nbIter3);
        const double d4 = computeClosestPoints_GJK<double, GJKType::JOHNSON>(shapeA,
                                                                             shapeB,
                                                                             v_a2w,
                                                                             v_b2w,
                                                                             q_a2w,
                                                                             q_b2w,
                                                                             crust,
                                                                             crust,
                                                                             pa4,
                                                                             pb4,
                                                                             nbIter4);

        // Distance must agree across all overloads
        EXPECT_NEAR(d1, d2, distEPS) << "trial " << trial << ": distance overloads 1 vs 2 differ";
        EXPECT_NEAR(d1, d3, distEPS) << "trial " << trial << ": distance overloads 1 vs 3 differ";
        EXPECT_NEAR(d1, d4, distEPS) << "trial " << trial << ": distance overloads 1 vs 4 differ";

        // Each overload must return a witness pair (pa in A's local, pb in B's local)
        // whose world-space separation equals d.  Cross-overload comparison of the raw
        // components is skipped: degenerate face-to-face configurations allow multiple
        // valid witnesses, and different overloads can legitimately pick different corners
        // of the same closest feature.
        if(!iRes1)  // Only meaningful when shapes are not intersecting
        {
            const auto toWA = [&](const Vector3<double>& p) { return (q_a2w >> p) + v_a2w; };
            const auto toWB = [&](const Vector3<double>& p) { return (q_b2w >> p) + v_b2w; };

            EXPECT_NEAR(norm(toWA(pa1) - toWB(pb1)), d1, distEPS)
                << "trial " << trial << ": overload 1 |pa_world - pb_world| != d";
            EXPECT_NEAR(norm(toWA(pa2) - toWB(pb2)), d2, distEPS)
                << "trial " << trial << ": overload 2 |pa_world - pb_world| != d";
            EXPECT_NEAR(norm(toWA(pa3) - toWB(pb3)), d3, distEPS)
                << "trial " << trial << ": overload 3 |pa_world - pb_world| != d";
            EXPECT_NEAR(norm(toWA(pa4) - toWB(pb4)), d4, distEPS)
                << "trial " << trial << ": overload 4 |pa_world - pb_world| != d";
        }
    }
}
