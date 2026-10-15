#include <cmath>
#include <gtest/gtest.h>

#include "Box.hh"
#include "CollisionDetection.hh"
#include "Quaternion.hh"
#include "RigidBody.hh"
#include "Vector3.hh"

class CollisionDetectionTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        boxA = new Box<double>(1.0, 1.0, 1.0);
        boxB = new Box<double>(0.5, 0.5, 0.5);

        rigidBodyA = new RigidBody<double>(boxA, 0.1, 1000.0, 1);
        rigidBodyB = new RigidBody<double>(boxB, 0.1, 1000.0, 1);

        origin               = Vector3<double>(0.0, 0.0, 0.0);
        separated_position   = Vector3<double>(3.0, 0.0, 0.0);
        overlapping_position = Vector3<double>(0.5, 0.0, 0.0);
        touching_position    = Vector3<double>(1.5, 0.0, 0.0);

        identity_quaternion = Quaternion<double>(0.0, 0.0, 0.0, 1.0);

        // 45 degree rotation around Z axis
        double angle       = M_PI / 4.0;
        rotated_quaternion = Quaternion<double>(0.0, 0.0, sin(angle / 2.0), cos(angle / 2.0));
    }

    void TearDown() override
    {
        delete rigidBodyA;
        delete rigidBodyB;
    }

    Box<double>*       boxA;
    Box<double>*       boxB;
    RigidBody<double>* rigidBodyA;
    RigidBody<double>* rigidBodyB;
    Vector3<double>    origin;
    Vector3<double>    separated_position;
    Vector3<double>    overlapping_position;
    Vector3<double>    touching_position;
    Quaternion<double> identity_quaternion;
    Quaternion<double> rotated_quaternion;
    const double       EPSILON = EPS<double>;
};

TEST_F(CollisionDetectionTest, RelativeTransformationIntersection)
{
    bool result
        = intersectRigidBodies(*rigidBodyA, *rigidBodyB, overlapping_position, identity_quaternion);
    EXPECT_TRUE(result);

    result
        = intersectRigidBodies(*rigidBodyA, *rigidBodyB, separated_position, identity_quaternion);
    EXPECT_FALSE(result);
}

TEST_F(CollisionDetectionTest, WorldCoordinatesIntersection)
{
    bool result = intersectRigidBodies(*rigidBodyA,
                                       *rigidBodyB,
                                       origin,
                                       overlapping_position,
                                       identity_quaternion,
                                       identity_quaternion);
    EXPECT_TRUE(result);

    result = intersectRigidBodies(*rigidBodyA,
                                  *rigidBodyB,
                                  origin,
                                  separated_position,
                                  identity_quaternion,
                                  identity_quaternion);
    EXPECT_FALSE(result);
}

TEST_F(CollisionDetectionTest, RotatedRigidBodiesIntersection)
{
    bool result
        = intersectRigidBodies(*rigidBodyA, *rigidBodyB, separated_position, rotated_quaternion);
    EXPECT_FALSE(result);

    result
        = intersectRigidBodies(*rigidBodyA, *rigidBodyB, overlapping_position, rotated_quaternion);
    EXPECT_TRUE(result);
}

TEST_F(CollisionDetectionTest, CollisionConsistency)
{
    bool result1
        = intersectRigidBodies(*rigidBodyA, *rigidBodyB, overlapping_position, identity_quaternion);
    bool result2
        = intersectRigidBodies(*rigidBodyA, *rigidBodyB, overlapping_position, identity_quaternion);
    bool result3
        = intersectRigidBodies(*rigidBodyA, *rigidBodyB, overlapping_position, identity_quaternion);

    EXPECT_EQ(result1, result2);
    EXPECT_EQ(result2, result3);
    EXPECT_TRUE(result1);  // Should be true for overlapping case

    result1
        = intersectRigidBodies(*rigidBodyA, *rigidBodyB, separated_position, identity_quaternion);
    result2
        = intersectRigidBodies(*rigidBodyA, *rigidBodyB, separated_position, identity_quaternion);
    EXPECT_EQ(result1, result2);
    EXPECT_FALSE(result1);  // Should be false for separated case
}

TEST_F(CollisionDetectionTest, StressTest)
{
    int numTests     = 1000;
    int successCount = 0;

    for(int i = 0; i < numTests; ++i)
    {
        double          x = (i % 100) * 0.03 + 0.5;  // Range from 0.5 to ~3.5
        Vector3<double> testPos(x, 0.0, 0.0);

        bool result = intersectRigidBodies(*rigidBodyA, *rigidBodyB, testPos, identity_quaternion);
        if(result)
            successCount++;
    }

    EXPECT_GT(successCount, 0);
    EXPECT_LT(successCount, numTests);
}
