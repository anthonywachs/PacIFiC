#include "Matrix3.hh"
#include "MatrixMath.hh"
#include "Quaternion.hh"
#include "QuaternionMath.hh"
#include "Vector3.hh"
#include "VectorMath.hh"
#include <cmath>
#include <gtest/gtest.h>

class RotationMathTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        test_vec1 = Vector3<double>(1.12, 0.00, 7.12);
        test_vec2 = Vector3<double>(-2.4, 5.10, 0.08);
        test_vec3 = Vector3<double>(-0.4, 0.00, 0.00);
        test_vec4 = Vector3<double>(-4.0, 0.08, 1.27);

        quat1 = Quaternion<double>(M_PI / 4., M_PI / 4., M_PI / 4.);
        mat1  = Matrix3<double>(M_PI / 4., M_PI / 4., M_PI / 4.);
        quat2 = Quaternion<double>(M_PI / 2., M_PI / 2., M_PI / 2.);
        mat2  = Matrix3<double>(M_PI / 2., M_PI / 2., M_PI / 2.);
        quat3 = Quaternion<double>(-M_PI / 3., M_PI / 2., M_PI / 6.);
        mat3  = Matrix3<double>(-M_PI / 3., M_PI / 2., M_PI / 6.);
    }

    Vector3<double>    test_vec1, test_vec2, test_vec3, test_vec4;
    Quaternion<double> quat1, quat2, quat3;
    Matrix3<double>    mat1, mat2, mat3;
    const double       EPSILON = HIGHEPS<double>;
};

TEST_F(RotationMathTest, QuatAndMatrixConversions)
{
    Matrix3<double> mat_from_quat1 = quat1.toMatrix();
    Matrix3<double> mat_from_quat2 = quat2.toMatrix();
    Matrix3<double> mat_from_quat3 = quat3.toMatrix();

    EXPECT_TRUE(mat1 == mat_from_quat1);
    EXPECT_TRUE(mat2 == mat_from_quat2);
    EXPECT_TRUE(mat3 == mat_from_quat3);

    Quaternion<double> quat_from_mat1;
    quat_from_mat1.setQuaternion(mat1);
    Quaternion<double> quat_from_mat2;
    quat_from_mat2.setQuaternion(mat2);
    Quaternion<double> quat_from_mat3;
    quat_from_mat3.setQuaternion(mat3);

    EXPECT_TRUE(quat1 == quat_from_mat1);
    EXPECT_TRUE(quat2 == quat_from_mat2);
    EXPECT_TRUE(quat3 == quat_from_mat3);
}

TEST_F(RotationMathTest, PrincipalAxesRotations)
{
    std::vector<std::pair<Quaternion<double>, Matrix3<double>>> test_rotations
        = {{quat1, mat1}, {quat2, mat2}, {quat3, mat3}};

    std::vector<Vector3<double>> test_vectors = {Vector3<double>(1.0, 0.0, 0.0),
                                                 Vector3<double>(0.0, 1.0, 0.0),
                                                 Vector3<double>(0.0, 0.0, 1.0)};

    for(const auto& [quat, matrix] : test_rotations)
    {
        for(const auto& vec : test_vectors)
        {
            Vector3<double> matrix_result = matrix * vec;
            Vector3<double> quat_result   = quat >> vec;

            EXPECT_TRUE(matrix_result == quat_result);
        }
    }
}

TEST_F(RotationMathTest, PrincipalAxesInverseRotations)
{
    std::vector<std::pair<Quaternion<double>, Matrix3<double>>> test_rotations
        = {{quat1, mat1}, {quat2, mat2}, {quat3, mat3}};

    std::vector<Vector3<double>> test_vectors = {Vector3<double>(1.0, 0.0, 0.0),
                                                 Vector3<double>(0.0, 1.0, 0.0),
                                                 Vector3<double>(0.0, 0.0, 1.0)};

    for(const auto& [quat, matrix] : test_rotations)
    {
        for(const auto& vec : test_vectors)
        {
            Vector3<double> matrix_result = transpose(matrix) * vec;
            Vector3<double> quat_result   = quat << vec;

            EXPECT_TRUE(matrix_result == quat_result);
        }
    }
}

TEST_F(RotationMathTest, ArbitraryAxesRotations)
{
    std::vector<std::pair<Quaternion<double>, Matrix3<double>>> test_rotations
        = {{quat1, mat1}, {quat2, mat2}, {quat3, mat3}};

    std::vector<Vector3<double>> test_vectors = {test_vec1, test_vec2, test_vec3, test_vec4};

    for(const auto& [quat, matrix] : test_rotations)
    {
        for(const auto& vec : test_vectors)
        {
            Vector3<double> matrix_result = matrix * vec;
            Vector3<double> quat_result   = quat >> vec;

            EXPECT_TRUE(matrix_result == quat_result);
        }
    }
}

TEST_F(RotationMathTest, ArbitraryAxesInverseRotations)
{
    std::vector<std::pair<Quaternion<double>, Matrix3<double>>> test_rotations
        = {{quat1, mat1}, {quat2, mat2}, {quat3, mat3}};

    std::vector<Vector3<double>> test_vectors = {test_vec1, test_vec2, test_vec3, test_vec4};

    for(const auto& [quat, matrix] : test_rotations)
    {
        for(const auto& vec : test_vectors)
        {
            Vector3<double> matrix_result = transpose(matrix) * vec;
            Vector3<double> quat_result   = quat << vec;

            EXPECT_TRUE(matrix_result == quat_result);
        }
    }
}

TEST_F(RotationMathTest, RotationComposition)
{
    Matrix3<double> m1 = mat1 * mat2;
    Matrix3<double> m2 = mat2 * mat3;
    Matrix3<double> m3 = mat3 * mat1;

    Quaternion<double> q1 = quat1 * quat2;
    Quaternion<double> q2 = quat2 * quat3;
    Quaternion<double> q3 = quat3 * quat1;

    Matrix3<double> mat_from_q1 = q1.toMatrix();
    Matrix3<double> mat_from_q2 = q2.toMatrix();
    Matrix3<double> mat_from_q3 = q3.toMatrix();

    EXPECT_TRUE(m1 == mat_from_q1);
    EXPECT_TRUE(m2 == mat_from_q2);
    EXPECT_TRUE(m3 == mat_from_q3);

    Quaternion<double> quat_from_m1;
    quat_from_m1.setQuaternion(m1);
    Quaternion<double> quat_from_m2;
    quat_from_m2.setQuaternion(m2);
    Quaternion<double> quat_from_m3;
    quat_from_m3.setQuaternion(m3);

    EXPECT_TRUE(q1 == quat_from_m1);
    EXPECT_TRUE(q2 == quat_from_m2);
    EXPECT_TRUE(q3 == quat_from_m3);
}