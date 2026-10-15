#include <cmath>
#include <gtest/gtest.h>

#include "Vector3.hh"
#include "VectorMath.hh"

class Vector3Test : public ::testing::Test
{
protected:
    void SetUp() override
    {
        v1   = Vector3<double>(1.0, 2.0, 3.0);
        v2   = Vector3<double>(4.0, 5.0, 6.0);
        zero = Vector3<double>(0.0, 0.0, 0.0);

        zero_vec     = Vector3<double>(0.0, 0.0, 0.0);
        unit_x       = Vector3<double>(1.0, 0.0, 0.0);
        unit_y       = Vector3<double>(0.0, 1.0, 0.0);
        unit_z       = Vector3<double>(0.0, 0.0, 1.0);
        test_vec1    = Vector3<double>(3.0, 4.0, 5.0);
        test_vec2    = Vector3<double>(1.0, 2.0, 3.0);
        negative_vec = Vector3<double>(-1.5, -2.5, -3.5);
        large_vec    = Vector3<double>(1e6, 2e6, 3e6);
        small_vec    = Vector3<double>(1e-6, 2e-6, 3e-6);
    }

    Vector3<double> v1, v2, zero;
    Vector3<double> zero_vec, unit_x, unit_y, unit_z;
    Vector3<double> test_vec1, test_vec2, negative_vec;
    Vector3<double> large_vec, small_vec;
    const double    EPSILON = 1e-10;
};

TEST_F(Vector3Test, ConstructorAndAccessors)
{
    EXPECT_DOUBLE_EQ(v1[0], 1.0);
    EXPECT_DOUBLE_EQ(v1[1], 2.0);
    EXPECT_DOUBLE_EQ(v1[2], 3.0);
}

TEST_F(Vector3Test, CopyConstructor)
{
    Vector3<double> copy(v1);
    EXPECT_DOUBLE_EQ(copy[0], v1[0]);
    EXPECT_DOUBLE_EQ(copy[1], v1[1]);
    EXPECT_DOUBLE_EQ(copy[2], v1[2]);
}

TEST_F(Vector3Test, AssignmentOperator)
{
    Vector3<double> assigned = v2;
    EXPECT_DOUBLE_EQ(assigned[0], v2[0]);
    EXPECT_DOUBLE_EQ(assigned[1], v2[1]);
    EXPECT_DOUBLE_EQ(assigned[2], v2[2]);
}

TEST_F(Vector3Test, SetValueMethods)
{
    Vector3<double> test;
    test.setValue(7.0, 8.0, 9.0);
    EXPECT_DOUBLE_EQ(test[0], 7.0);
    EXPECT_DOUBLE_EQ(test[1], 8.0);
    EXPECT_DOUBLE_EQ(test[2], 9.0);
}

TEST_F(Vector3Test, NormalizationMethod)
{
    Vector3<double> test(3.0, 4.0, 0.0);
    test.normalize();

    double mag_squared = test[0] * test[0] + test[1] * test[1] + test[2] * test[2];
    EXPECT_NEAR(mag_squared, 1.0, EPSILON);

    EXPECT_NEAR(test[0], 0.6, EPSILON);
    EXPECT_NEAR(test[1], 0.8, EPSILON);
    EXPECT_NEAR(test[2], 0.0, EPSILON);
}

TEST_F(Vector3Test, ResetMethod)
{
    Vector3<double> test(1.0, 2.0, 3.0);
    test.reset();
    EXPECT_DOUBLE_EQ(test[0], 0.0);
    EXPECT_DOUBLE_EQ(test[1], 0.0);
    EXPECT_DOUBLE_EQ(test[2], 0.0);
}

TEST_F(Vector3Test, NormFunction)
{
    EXPECT_NEAR(norm(unit_x), 1.0, EPSILON);
    EXPECT_NEAR(norm(unit_y), 1.0, EPSILON);
    EXPECT_NEAR(norm(unit_z), 1.0, EPSILON);

    EXPECT_NEAR(norm(zero_vec), 0.0, EPSILON);

    EXPECT_NEAR(norm(test_vec1), sqrt(50.0), EPSILON);

    EXPECT_NEAR(norm(negative_vec), sqrt(1.5 * 1.5 + 2.5 * 2.5 + 3.5 * 3.5), EPSILON);
}

TEST_F(Vector3Test, Norm2Function)
{
    EXPECT_NEAR(norm2(unit_x), 1.0, EPSILON);
    EXPECT_NEAR(norm2(unit_y), 1.0, EPSILON);
    EXPECT_NEAR(norm2(unit_z), 1.0, EPSILON);

    EXPECT_NEAR(norm2(zero_vec), 0.0, EPSILON);

    EXPECT_NEAR(norm2(test_vec1), 50.0, EPSILON);

    EXPECT_NEAR(norm2(test_vec2), 14.0, EPSILON);
}

TEST_F(Vector3Test, IsApproxZeroFunction)
{
    EXPECT_TRUE(isApproxZero(zero_vec, EPSILON));

    EXPECT_TRUE(isApproxZero(small_vec, 1e-5));
    EXPECT_FALSE(isApproxZero(small_vec, 1e-7));

    EXPECT_FALSE(isApproxZero(unit_x, EPSILON));
    EXPECT_FALSE(isApproxZero(unit_y, EPSILON));
    EXPECT_FALSE(isApproxZero(unit_z, EPSILON));

    EXPECT_FALSE(isApproxZero(test_vec1, EPSILON));
    EXPECT_FALSE(isApproxZero(large_vec, EPSILON));
}

TEST_F(Vector3Test, RoundFunction)
{
    Vector3<double> decimal_vec(1e-12, 2e-12, -3e-12);
    round(decimal_vec, 1e-10);

    EXPECT_DOUBLE_EQ(decimal_vec[0], 0.0);
    EXPECT_DOUBLE_EQ(decimal_vec[1], 0.0);
    EXPECT_DOUBLE_EQ(decimal_vec[2], 0.0);

    Vector3<double> large_vec(0.1, -0.2, 0.3);
    round(large_vec, 1e-10);
    EXPECT_DOUBLE_EQ(large_vec[0], 0.1);
    EXPECT_DOUBLE_EQ(large_vec[1], -0.2);
    EXPECT_DOUBLE_EQ(large_vec[2], 0.3);
}

TEST_F(Vector3Test, ArithmeticOperators)
{
    Vector3<double> sum = test_vec1 + test_vec2;
    EXPECT_DOUBLE_EQ(sum[0], 4.0);
    EXPECT_DOUBLE_EQ(sum[1], 6.0);
    EXPECT_DOUBLE_EQ(sum[2], 8.0);

    Vector3<double> diff = test_vec1 - test_vec2;
    EXPECT_DOUBLE_EQ(diff[0], 2.0);
    EXPECT_DOUBLE_EQ(diff[1], 2.0);
    EXPECT_DOUBLE_EQ(diff[2], 2.0);

    Vector3<double> scaled = 2.0 * test_vec1;
    EXPECT_DOUBLE_EQ(scaled[0], 6.0);
    EXPECT_DOUBLE_EQ(scaled[1], 8.0);
    EXPECT_DOUBLE_EQ(scaled[2], 10.0);

    Vector3<double> divided = test_vec1 / 2.0;
    EXPECT_DOUBLE_EQ(divided[0], 1.5);
    EXPECT_DOUBLE_EQ(divided[1], 2.0);
    EXPECT_DOUBLE_EQ(divided[2], 2.5);

    Vector3<double> negated = -test_vec1;
    EXPECT_DOUBLE_EQ(negated[0], -3.0);
    EXPECT_DOUBLE_EQ(negated[1], -4.0);
    EXPECT_DOUBLE_EQ(negated[2], -5.0);
}

TEST_F(Vector3Test, DotProduct)
{
    EXPECT_NEAR(unit_x * unit_y, 0.0, EPSILON);
    EXPECT_NEAR(unit_x * unit_z, 0.0, EPSILON);
    EXPECT_NEAR(unit_y * unit_z, 0.0, EPSILON);

    EXPECT_NEAR(unit_x * unit_x, 1.0, EPSILON);
    EXPECT_NEAR(test_vec1 * test_vec1, norm2(test_vec1), EPSILON);

    EXPECT_NEAR(test_vec1 * test_vec2, 26.0, EPSILON);

    EXPECT_NEAR(test_vec1 * zero_vec, 0.0, EPSILON);
}

TEST_F(Vector3Test, CrossProduct)
{
    Vector3<double> cross_xy = unit_x ^ unit_y;
    EXPECT_NEAR(cross_xy[0], 0.0, EPSILON);
    EXPECT_NEAR(cross_xy[1], 0.0, EPSILON);
    EXPECT_NEAR(cross_xy[2], 1.0, EPSILON);

    Vector3<double> cross_yz = unit_y ^ unit_z;
    EXPECT_NEAR(cross_yz[0], 1.0, EPSILON);
    EXPECT_NEAR(cross_yz[1], 0.0, EPSILON);
    EXPECT_NEAR(cross_yz[2], 0.0, EPSILON);

    Vector3<double> cross_zx = unit_z ^ unit_x;
    EXPECT_NEAR(cross_zx[0], 0.0, EPSILON);
    EXPECT_NEAR(cross_zx[1], 1.0, EPSILON);
    EXPECT_NEAR(cross_zx[2], 0.0, EPSILON);

    Vector3<double> cross1 = test_vec1 ^ test_vec2;
    Vector3<double> cross2 = test_vec2 ^ test_vec1;
    EXPECT_NEAR(cross1[0], -cross2[0], EPSILON);
    EXPECT_NEAR(cross1[1], -cross2[1], EPSILON);
    EXPECT_NEAR(cross1[2], -cross2[2], EPSILON);

    Vector3<double> self_cross = test_vec1 ^ test_vec1;
    EXPECT_NEAR(self_cross[0], 0.0, EPSILON);
    EXPECT_NEAR(self_cross[1], 0.0, EPSILON);
    EXPECT_NEAR(self_cross[2], 0.0, EPSILON);

    Vector3<double> result = test_vec1 ^ test_vec2;
    EXPECT_NEAR(result[0], 2.0, EPSILON);
    EXPECT_NEAR(result[1], -4.0, EPSILON);
    EXPECT_NEAR(result[2], 2.0, EPSILON);
}

TEST_F(Vector3Test, ComparisonOperators)
{
    Vector3<double> v_copy = test_vec1;
    Vector3<double> v_different(3.0, 4.0, 5.1);

    EXPECT_TRUE(test_vec1 == v_copy);
    EXPECT_FALSE(test_vec1 == test_vec2);
    EXPECT_FALSE(test_vec1 == v_different);

    EXPECT_FALSE(test_vec1 != v_copy);
    EXPECT_TRUE(test_vec1 != test_vec2);
    EXPECT_TRUE(test_vec1 != v_different);

    EXPECT_TRUE(zero_vec == zero_vec);
    EXPECT_FALSE(zero_vec == unit_x);
}

TEST_F(Vector3Test, CompoundAssignmentOperators)
{
    Vector3<double> test = test_vec1;

    test += test_vec2;
    EXPECT_DOUBLE_EQ(test[0], 4.0);
    EXPECT_DOUBLE_EQ(test[1], 6.0);
    EXPECT_DOUBLE_EQ(test[2], 8.0);

    test -= test_vec2;
    EXPECT_DOUBLE_EQ(test[0], 3.0);
    EXPECT_DOUBLE_EQ(test[1], 4.0);
    EXPECT_DOUBLE_EQ(test[2], 5.0);

    test *= 2.0;
    EXPECT_DOUBLE_EQ(test[0], 6.0);
    EXPECT_DOUBLE_EQ(test[1], 8.0);
    EXPECT_DOUBLE_EQ(test[2], 10.0);

    test /= 2.0;
    EXPECT_DOUBLE_EQ(test[0], 3.0);
    EXPECT_DOUBLE_EQ(test[1], 4.0);
    EXPECT_DOUBLE_EQ(test[2], 5.0);
}

TEST_F(Vector3Test, EdgeCases)
{
    double large_norm = norm(large_vec);
    EXPECT_TRUE(std::isfinite(large_norm));
    EXPECT_GT(large_norm, 0.0);

    double small_norm = norm(small_vec);
    EXPECT_TRUE(std::isfinite(small_norm));
    EXPECT_GT(small_norm, 0.0);

    Vector3<double> sum = large_vec + small_vec;
    EXPECT_TRUE(std::isfinite(sum[0]));
    EXPECT_TRUE(std::isfinite(sum[1]));
    EXPECT_TRUE(std::isfinite(sum[2]));

    Vector3<double> cross = unit_x ^ unit_y;
    EXPECT_NEAR(norm(cross), 1.0, EPSILON);
}

TEST_F(Vector3Test, BufferAccess)
{
    const double* buffer = v1.getBuffer();
    EXPECT_DOUBLE_EQ(buffer[0], 1.0);
    EXPECT_DOUBLE_EQ(buffer[1], 2.0);
    EXPECT_DOUBLE_EQ(buffer[2], 3.0);
}
