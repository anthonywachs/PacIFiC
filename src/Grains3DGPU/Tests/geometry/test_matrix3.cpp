#include <cmath>
#include <gtest/gtest.h>

#include "Matrix3.hh"
#include "MatrixMath.hh"
#include "Vector3.hh"

class Matrix3Test : public ::testing::Test
{
protected:
    void SetUp() override
    {
        testMatrix = Matrix3<double>(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
        identity   = Matrix3<double>();
        invertible = Matrix3<double>(2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 4.0);
        rotation_z = Matrix3<double>(0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        v1         = Vector3<double>(1.0, 2.0, 3.0);
        v2         = Vector3<double>(4.0, 5.0, 6.0);
    }

    Matrix3<double> identity, testMatrix, invertible, rotation_z;
    Vector3<double> v1, v2;
    const double    EPSILON = 1e-10;
};

TEST_F(Matrix3Test, ElementAccess)
{
    EXPECT_DOUBLE_EQ(testMatrix(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(testMatrix(0, 1), 2.0);
    EXPECT_DOUBLE_EQ(testMatrix(0, 2), 3.0);
    EXPECT_DOUBLE_EQ(testMatrix(1, 0), 4.0);
    EXPECT_DOUBLE_EQ(testMatrix(1, 1), 5.0);
    EXPECT_DOUBLE_EQ(testMatrix(1, 2), 6.0);
    EXPECT_DOUBLE_EQ(testMatrix(2, 0), 7.0);
    EXPECT_DOUBLE_EQ(testMatrix(2, 1), 8.0);
    EXPECT_DOUBLE_EQ(testMatrix(2, 2), 9.0);
}

TEST_F(Matrix3Test, RowAccess)
{
    Vector3<double> row0 = testMatrix[0];
    Vector3<double> row1 = testMatrix[1];
    Vector3<double> row2 = testMatrix[2];

    EXPECT_DOUBLE_EQ(row0[0], 1.0);
    EXPECT_DOUBLE_EQ(row0[1], 2.0);
    EXPECT_DOUBLE_EQ(row0[2], 3.0);

    EXPECT_DOUBLE_EQ(row1[0], 4.0);
    EXPECT_DOUBLE_EQ(row1[1], 5.0);
    EXPECT_DOUBLE_EQ(row1[2], 6.0);

    EXPECT_DOUBLE_EQ(row2[0], 7.0);
    EXPECT_DOUBLE_EQ(row2[1], 8.0);
    EXPECT_DOUBLE_EQ(row2[2], 9.0);
}

TEST_F(Matrix3Test, LinearIndexing)
{
    EXPECT_DOUBLE_EQ(testMatrix(0), 1.0);
    EXPECT_DOUBLE_EQ(testMatrix(1), 2.0);
    EXPECT_DOUBLE_EQ(testMatrix(2), 3.0);
    EXPECT_DOUBLE_EQ(testMatrix(3), 4.0);
    EXPECT_DOUBLE_EQ(testMatrix(8), 9.0);
}

TEST_F(Matrix3Test, IdentityMatrix)
{
    EXPECT_DOUBLE_EQ(identity(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(identity(1, 1), 1.0);
    EXPECT_DOUBLE_EQ(identity(2, 2), 1.0);
    EXPECT_DOUBLE_EQ(identity(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(identity(0, 2), 0.0);
    EXPECT_DOUBLE_EQ(identity(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(identity(1, 2), 0.0);
    EXPECT_DOUBLE_EQ(identity(2, 0), 0.0);
    EXPECT_DOUBLE_EQ(identity(2, 1), 0.0);
}

TEST_F(Matrix3Test, CopyConstructor)
{
    Matrix3<double> copy(testMatrix);
    EXPECT_DOUBLE_EQ(copy(1, 1), testMatrix(1, 1));
    EXPECT_DOUBLE_EQ(copy(2, 0), testMatrix(2, 0));
    EXPECT_DOUBLE_EQ(copy(0, 2), testMatrix(0, 2));
}

TEST_F(Matrix3Test, AssignmentOperator)
{
    Matrix3<double> assigned = testMatrix;
    EXPECT_DOUBLE_EQ(assigned(0, 0), testMatrix(0, 0));
    EXPECT_DOUBLE_EQ(assigned(2, 2), testMatrix(2, 2));
    EXPECT_DOUBLE_EQ(assigned(1, 2), testMatrix(1, 2));
}

TEST_F(Matrix3Test, SetValueMethods)
{
    Matrix3<double> test;
    test.setValue(10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0);
    EXPECT_DOUBLE_EQ(test(0, 0), 10.0);
    EXPECT_DOUBLE_EQ(test(1, 1), 14.0);
    EXPECT_DOUBLE_EQ(test(2, 2), 18.0);
    EXPECT_DOUBLE_EQ(test(0, 1), 11.0);
    EXPECT_DOUBLE_EQ(test(2, 1), 17.0);
}

TEST_F(Matrix3Test, ConstAccess)
{
    const Matrix3<double> const_matrix = testMatrix;

    EXPECT_DOUBLE_EQ(const_matrix[0][0], 1.0);
    EXPECT_DOUBLE_EQ(const_matrix[0][1], 2.0);
    EXPECT_DOUBLE_EQ(const_matrix[0][2], 3.0);
    EXPECT_DOUBLE_EQ(const_matrix[1][0], 4.0);
    EXPECT_DOUBLE_EQ(const_matrix[1][1], 5.0);
    EXPECT_DOUBLE_EQ(const_matrix[1][2], 6.0);
    EXPECT_DOUBLE_EQ(const_matrix[2][0], 7.0);
    EXPECT_DOUBLE_EQ(const_matrix[2][1], 8.0);
    EXPECT_DOUBLE_EQ(const_matrix[2][2], 9.0);
}

TEST_F(Matrix3Test, FabsFunction)
{
    const Matrix3<double> negative(-1.0, -2.0, 3.0, 4.0, -5.0, -6.0, -7.0, 8.0, -9.0);

    Matrix3<double> result = fabs<double>(negative);

    EXPECT_DOUBLE_EQ(result(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(result(0, 1), 2.0);
    EXPECT_DOUBLE_EQ(result(0, 2), 3.0);
    EXPECT_DOUBLE_EQ(result(1, 0), 4.0);
    EXPECT_DOUBLE_EQ(result(1, 1), 5.0);
    EXPECT_DOUBLE_EQ(result(1, 2), 6.0);
    EXPECT_DOUBLE_EQ(result(2, 0), 7.0);
    EXPECT_DOUBLE_EQ(result(2, 1), 8.0);
    EXPECT_DOUBLE_EQ(result(2, 2), 9.0);
}

TEST_F(Matrix3Test, FabsInPlace)
{
    Matrix3<double> negative(-1.0, -2.0, 3.0, 4.0, -5.0, -6.0, -7.0, 8.0, -9.0);

    fabs(negative);

    EXPECT_DOUBLE_EQ(negative(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(negative(0, 1), 2.0);
    EXPECT_DOUBLE_EQ(negative(0, 2), 3.0);
    EXPECT_DOUBLE_EQ(negative(1, 0), 4.0);
    EXPECT_DOUBLE_EQ(negative(1, 1), 5.0);
    EXPECT_DOUBLE_EQ(negative(1, 2), 6.0);
    EXPECT_DOUBLE_EQ(negative(2, 0), 7.0);
    EXPECT_DOUBLE_EQ(negative(2, 1), 8.0);
    EXPECT_DOUBLE_EQ(negative(2, 2), 9.0);
}

TEST_F(Matrix3Test, Determinant)
{
    EXPECT_NEAR(determinant(identity), 1.0, EPSILON);

    EXPECT_NEAR(determinant(testMatrix), 0.0, EPSILON);

    double det = determinant(invertible);
    EXPECT_GT(std::abs(det), EPSILON);

    EXPECT_NEAR(determinant(rotation_z), 1.0, EPSILON);
}

TEST_F(Matrix3Test, Transpose)
{
    const Matrix3<double> test_const = testMatrix;
    Matrix3<double>       result     = transpose<double>(test_const);

    EXPECT_DOUBLE_EQ(result(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(result(0, 1), 4.0);
    EXPECT_DOUBLE_EQ(result(0, 2), 7.0);
    EXPECT_DOUBLE_EQ(result(1, 0), 2.0);
    EXPECT_DOUBLE_EQ(result(1, 1), 5.0);
    EXPECT_DOUBLE_EQ(result(1, 2), 8.0);
    EXPECT_DOUBLE_EQ(result(2, 0), 3.0);
    EXPECT_DOUBLE_EQ(result(2, 1), 6.0);
    EXPECT_DOUBLE_EQ(result(2, 2), 9.0);
}

TEST_F(Matrix3Test, TransposeInPlace)
{
    Matrix3<double> matrix = testMatrix;
    transpose(matrix);

    EXPECT_DOUBLE_EQ(matrix(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(matrix(0, 1), 4.0);
    EXPECT_DOUBLE_EQ(matrix(0, 2), 7.0);
    EXPECT_DOUBLE_EQ(matrix(1, 0), 2.0);
    EXPECT_DOUBLE_EQ(matrix(1, 1), 5.0);
    EXPECT_DOUBLE_EQ(matrix(1, 2), 8.0);
    EXPECT_DOUBLE_EQ(matrix(2, 0), 3.0);
    EXPECT_DOUBLE_EQ(matrix(2, 1), 6.0);
    EXPECT_DOUBLE_EQ(matrix(2, 2), 9.0);
}

TEST_F(Matrix3Test, Inverse)
{
    const Matrix3<double> invertible_const = invertible;
    Matrix3<double>       inv              = inverse<double>(invertible_const);

    EXPECT_NEAR(inv(0, 0), 0.5, EPSILON);        // 1/2
    EXPECT_NEAR(inv(1, 1), 1.0 / 3.0, EPSILON);  // 1/3
    EXPECT_NEAR(inv(2, 2), 0.25, EPSILON);       // 1/4

    EXPECT_NEAR(inv(0, 1), 0.0, EPSILON);
    EXPECT_NEAR(inv(0, 2), 0.0, EPSILON);
    EXPECT_NEAR(inv(1, 0), 0.0, EPSILON);
    EXPECT_NEAR(inv(1, 2), 0.0, EPSILON);
    EXPECT_NEAR(inv(2, 0), 0.0, EPSILON);
    EXPECT_NEAR(inv(2, 1), 0.0, EPSILON);

    Matrix3<double> product = invertible * inv;

    for(int i = 0; i < 3; ++i)
    {
        for(int j = 0; j < 3; ++j)
        {
            if(i == j)
                EXPECT_NEAR(product(i, j), 1.0, EPSILON);
            else
                EXPECT_NEAR(product(i, j), 0.0, EPSILON);
        }
    }
}

TEST_F(Matrix3Test, InverseInPlace)
{
    Matrix3<double> original = invertible;
    Matrix3<double> matrix   = invertible;
    inverse<double>(matrix);

    EXPECT_NEAR(matrix(0, 0), 0.5, EPSILON);        // 1/2
    EXPECT_NEAR(matrix(1, 1), 1.0 / 3.0, EPSILON);  // 1/3
    EXPECT_NEAR(matrix(2, 2), 0.25, EPSILON);       // 1/4

    Matrix3<double> product = original * matrix;

    for(int i = 0; i < 3; ++i)
    {
        for(int j = 0; j < 3; ++j)
        {
            if(i == j)
                EXPECT_NEAR(product(i, j), 1.0, EPSILON);
            else
                EXPECT_NEAR(product(i, j), 0.0, EPSILON);
        }
    }
}

TEST_F(Matrix3Test, Scale)
{
    Vector3<double>       scale_factors(2.0, 3.0, 4.0);
    const Matrix3<double> identity_const = identity;
    Matrix3<double>       result         = scale<double>(identity_const, scale_factors);

    EXPECT_DOUBLE_EQ(result(0, 0), 2.0);
    EXPECT_DOUBLE_EQ(result(1, 1), 3.0);
    EXPECT_DOUBLE_EQ(result(2, 2), 4.0);

    EXPECT_DOUBLE_EQ(result(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(result(0, 2), 0.0);
    EXPECT_DOUBLE_EQ(result(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(result(1, 2), 0.0);
    EXPECT_DOUBLE_EQ(result(2, 0), 0.0);
    EXPECT_DOUBLE_EQ(result(2, 1), 0.0);
}

TEST_F(Matrix3Test, ScaleInPlace)
{
    Matrix3<double> matrix = identity;
    Vector3<double> scale_factors(2.0, 3.0, 4.0);

    scale(matrix, scale_factors);

    EXPECT_DOUBLE_EQ(matrix(0, 0), 2.0);
    EXPECT_DOUBLE_EQ(matrix(1, 1), 3.0);
    EXPECT_DOUBLE_EQ(matrix(2, 2), 4.0);
}

TEST_F(Matrix3Test, ArithmeticOperations)
{
    Matrix3<double> sum = identity + testMatrix;
    EXPECT_DOUBLE_EQ(sum(0, 0), 2.0);  // 1 + 1
    EXPECT_DOUBLE_EQ(sum(0, 1), 2.0);  // 0 + 2
    EXPECT_DOUBLE_EQ(sum(1, 1), 6.0);  // 1 + 5

    Matrix3<double> diff = testMatrix - identity;
    EXPECT_DOUBLE_EQ(diff(0, 0), 0.0);  // 1 - 1
    EXPECT_DOUBLE_EQ(diff(0, 1), 2.0);  // 2 - 0
    EXPECT_DOUBLE_EQ(diff(1, 1), 4.0);  // 5 - 1

    Matrix3<double> neg = -testMatrix;
    EXPECT_DOUBLE_EQ(neg(0, 0), -1.0);
    EXPECT_DOUBLE_EQ(neg(1, 1), -5.0);
    EXPECT_DOUBLE_EQ(neg(2, 2), -9.0);

    double          scalar = 2.5;
    Matrix3<double> scaled = scalar * identity;
    EXPECT_DOUBLE_EQ(scaled(0, 0), 2.5);
    EXPECT_DOUBLE_EQ(scaled(1, 1), 2.5);
    EXPECT_DOUBLE_EQ(scaled(2, 2), 2.5);
}

TEST_F(Matrix3Test, MatrixVectorOperations)
{
    Vector3<double> result = identity * v1;
    EXPECT_DOUBLE_EQ(result[0], 1.0);
    EXPECT_DOUBLE_EQ(result[1], 2.0);
    EXPECT_DOUBLE_EQ(result[2], 3.0);

    Vector3<double> unit_x(1.0, 0.0, 0.0);
    Vector3<double> rotated = rotation_z * unit_x;
    EXPECT_NEAR(rotated[0], 0.0, EPSILON);
    EXPECT_NEAR(rotated[1], 1.0, EPSILON);
    EXPECT_NEAR(rotated[2], 0.0, EPSILON);

    Vector3<double> result2 = v1 * identity;
    EXPECT_DOUBLE_EQ(result2[0], 1.0);
    EXPECT_DOUBLE_EQ(result2[1], 2.0);
    EXPECT_DOUBLE_EQ(result2[2], 3.0);
}

TEST_F(Matrix3Test, MatrixMultiplication)
{
    Matrix3<double> result = identity * testMatrix;

    for(int i = 0; i < 3; ++i)
    {
        for(int j = 0; j < 3; ++j)
            EXPECT_DOUBLE_EQ(result(i, j), testMatrix(i, j));
    }

    Matrix3<double> double_rotation = rotation_z * rotation_z;
    Vector3<double> unit_x(1.0, 0.0, 0.0);
    Vector3<double> rotated = double_rotation * unit_x;

    EXPECT_NEAR(rotated[0], -1.0, EPSILON);
    EXPECT_NEAR(rotated[1], 0.0, EPSILON);
    EXPECT_NEAR(rotated[2], 0.0, EPSILON);
}

TEST_F(Matrix3Test, IsRotation)
{
    EXPECT_TRUE(isRotation(identity));

    EXPECT_TRUE(isRotation(rotation_z));

    EXPECT_FALSE(isRotation(testMatrix));

    Matrix3<double> scaled = 2.0 * identity;
    EXPECT_FALSE(isRotation(scaled));
}

TEST_F(Matrix3Test, EdgeCases)
{
    Matrix3<double> zero(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    EXPECT_DOUBLE_EQ(determinant(zero), 0.0);
    EXPECT_FALSE(isRotation(zero));

    Matrix3<double> small(1e-12, 0.0, 0.0, 0.0, 1e-12, 0.0, 0.0, 0.0, 1e-12);

    EXPECT_NEAR(determinant(small), 1e-36, 1e-38);
}
