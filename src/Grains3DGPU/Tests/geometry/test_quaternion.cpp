#include <cmath>
#include <gtest/gtest.h>

#include "Quaternion.hh"
#include "QuaternionMath.hh"

class QuaternionTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        identity        = Quaternion<double>(0.0, 0.0, 0.0, 1.0);
        unit_i          = Quaternion<double>(1.0, 0.0, 0.0, 0.0);
        unit_j          = Quaternion<double>(0.0, 1.0, 0.0, 0.0);
        unit_k          = Quaternion<double>(0.0, 0.0, 1.0, 0.0);
        test_q1         = Quaternion<double>(1.0, 2.0, 3.0, 4.0);
        test_q2         = Quaternion<double>(2.0, 1.0, 0.5, 2.5);
        zero_q          = Quaternion<double>(0.0, 0.0, 0.0, 0.0);
        double norm_val = std::sqrt(1.0 + 4.0 + 9.0 + 16.0);  // norm of test_q1
        normalized_q
            = Quaternion<double>(1.0 / norm_val, 2.0 / norm_val, 3.0 / norm_val, 4.0 / norm_val);
        test_vec1 = Vector3<double>(1.0, 0.0, 0.0);
        test_vec2 = Vector3<double>(0.0, 1.0, 0.0);
        test_vec3 = Vector3<double>(1.0, 2.0, 3.0);
    }

    Quaternion<double> identity, unit_i, unit_j, unit_k;
    Quaternion<double> test_q1, test_q2, zero_q, normalized_q;
    Vector3<double>    test_vec1, test_vec2, test_vec3;
    const double       EPSILON = 1e-10;
};

TEST_F(QuaternionTest, BasicQuaternionCreation)
{
    Quaternion<double> q(1.0, 2.0, 3.0, 4.0);

    EXPECT_NEAR(q[0], 1.0, 1e-10);
    EXPECT_NEAR(q[1], 2.0, 1e-10);
    EXPECT_NEAR(q[2], 3.0, 1e-10);
    EXPECT_NEAR(q[3], 4.0, 1e-10);

    EXPECT_NEAR(q.getScalar(), 4.0, 1e-10);
    Vector3<double> vec = q.getVector();
    EXPECT_NEAR(vec[0], 1.0, 1e-10);
    EXPECT_NEAR(vec[1], 2.0, 1e-10);
    EXPECT_NEAR(vec[2], 3.0, 1e-10);
}

TEST_F(QuaternionTest, QuaternionNorm)
{
    Quaternion<double> q(1.0, 2.0, 3.0, 4.0);

    double n        = norm(q);
    double expected = std::sqrt(1.0 * 1.0 + 2.0 * 2.0 + 3.0 * 3.0 + 4.0 * 4.0);

    EXPECT_NEAR(n, expected, 1e-10);

    double n2 = norm2(q);
    EXPECT_NEAR(n2, expected * expected, 1e-10);
}

TEST_F(QuaternionTest, QuaternionMultiplication)
{
    Quaternion<double> q1(1.0, 0.0, 0.0, 0.0);
    Quaternion<double> q2(0.0, 1.0, 0.0, 0.0);
    Quaternion<double> result = q1 * q2;

    double n = norm(result);
    EXPECT_GT(n, 0.0);

    EXPECT_NO_THROW(q1 * q2);
}

TEST_F(QuaternionTest, QuaternionVectorOperations)
{
    Quaternion<double> q(0.0, 0.0, 0.0, 1.0);
    Vector3<double>    v(1.0, 2.0, 3.0);

    Quaternion<double> result = q * v;

    double n = norm(result);
    EXPECT_GT(n, 0.0);
}

TEST_F(QuaternionTest, QuaternionSetters)
{
    Quaternion<double> q;
    q.setScalar(5.0);
    Vector3<double> v(1.0, 2.0, 3.0);
    q.setVector(v);

    EXPECT_NEAR(q.getScalar(), 5.0, 1e-10);
    Vector3<double> vec_result = q.getVector();
    EXPECT_NEAR(vec_result[0], 1.0, 1e-10);
    EXPECT_NEAR(vec_result[1], 2.0, 1e-10);
    EXPECT_NEAR(vec_result[2], 3.0, 1e-10);

    q.setQuaternion(2.0, 3.0, 4.0, 6.0);
    EXPECT_NEAR(q[0], 2.0, 1e-10);
    EXPECT_NEAR(q[1], 3.0, 1e-10);
    EXPECT_NEAR(q[2], 4.0, 1e-10);
    EXPECT_NEAR(q[3], 6.0, 1e-10);
}

TEST_F(QuaternionTest, NormFunctions)
{
    EXPECT_NEAR(norm(identity), 1.0, EPSILON);
    EXPECT_NEAR(norm2(identity), 1.0, EPSILON);

    EXPECT_NEAR(norm(zero_q), 0.0, EPSILON);
    EXPECT_NEAR(norm2(zero_q), 0.0, EPSILON);

    EXPECT_NEAR(norm(unit_i), 1.0, EPSILON);
    EXPECT_NEAR(norm(unit_j), 1.0, EPSILON);
    EXPECT_NEAR(norm(unit_k), 1.0, EPSILON);

    double expected_norm = std::sqrt(30.0);
    EXPECT_NEAR(norm(test_q1), expected_norm, EPSILON);
    EXPECT_NEAR(norm2(test_q1), 30.0, EPSILON);

    EXPECT_NEAR(norm(normalized_q), 1.0, EPSILON);
}

TEST_F(QuaternionTest, ConjugateFunction)
{
    const Quaternion<double> const_identity = identity;
    Quaternion<double>       conj_identity  = conjugate(const_identity);
    EXPECT_NEAR(conj_identity[0], 0.0, EPSILON);
    EXPECT_NEAR(conj_identity[1], 0.0, EPSILON);
    EXPECT_NEAR(conj_identity[2], 0.0, EPSILON);
    EXPECT_NEAR(conj_identity[3], 1.0, EPSILON);

    const Quaternion<double> const_test_q1 = test_q1;
    Quaternion<double>       conj_test     = conjugate(const_test_q1);
    EXPECT_NEAR(conj_test[0], -1.0, EPSILON);
    EXPECT_NEAR(conj_test[1], -2.0, EPSILON);
    EXPECT_NEAR(conj_test[2], -3.0, EPSILON);
    EXPECT_NEAR(conj_test[3], 4.0, EPSILON);

    Quaternion<double> test_copy = test_q1;
    conjugate(test_copy);
    EXPECT_NEAR(test_copy[0], -1.0, EPSILON);
    EXPECT_NEAR(test_copy[1], -2.0, EPSILON);
    EXPECT_NEAR(test_copy[2], -3.0, EPSILON);
    EXPECT_NEAR(test_copy[3], 4.0, EPSILON);

    Quaternion<double>       temp_q       = conjugate(const_test_q1);
    const Quaternion<double> const_temp_q = temp_q;
    Quaternion<double>       double_conj  = conjugate(const_temp_q);
    EXPECT_NEAR(double_conj[0], test_q1[0], EPSILON);
    EXPECT_NEAR(double_conj[1], test_q1[1], EPSILON);
    EXPECT_NEAR(double_conj[2], test_q1[2], EPSILON);
    EXPECT_NEAR(double_conj[3], test_q1[3], EPSILON);
}

TEST_F(QuaternionTest, InverseFunction)
{
    const Quaternion<double> const_identity = identity;
    Quaternion<double>       inv_identity   = inverse(const_identity);
    EXPECT_NEAR(norm(inv_identity - identity), 0.0, EPSILON);

    const Quaternion<double> const_normalized_q = normalized_q;
    Quaternion<double>       inv_normalized     = inverse(const_normalized_q);

    Quaternion<double> conj_normalized = conjugate(const_normalized_q);
    EXPECT_NEAR(inv_normalized[0], conj_normalized[0], EPSILON);
    EXPECT_NEAR(inv_normalized[1], conj_normalized[1], EPSILON);
    EXPECT_NEAR(inv_normalized[2], conj_normalized[2], EPSILON);
    EXPECT_NEAR(inv_normalized[3], conj_normalized[3], EPSILON);

    Quaternion<double> product = normalized_q * inv_normalized;
    EXPECT_NEAR(product[0], 0.0, EPSILON);
    EXPECT_NEAR(product[1], 0.0, EPSILON);
    EXPECT_NEAR(product[2], 0.0, EPSILON);
    EXPECT_NEAR(product[3], 1.0, EPSILON);

    Quaternion<double> test_copy = normalized_q;
    inverse(test_copy);
    EXPECT_NEAR(test_copy[0], conj_normalized[0], EPSILON);
    EXPECT_NEAR(test_copy[1], conj_normalized[1], EPSILON);
    EXPECT_NEAR(test_copy[2], conj_normalized[2], EPSILON);
    EXPECT_NEAR(test_copy[3], conj_normalized[3], EPSILON);
}

TEST_F(QuaternionTest, AdditionOperators)
{
    Quaternion<double> sum = test_q1 + test_q2;
    EXPECT_NEAR(sum[0], 3.0, EPSILON);  // 1+2
    EXPECT_NEAR(sum[1], 3.0, EPSILON);  // 2+1
    EXPECT_NEAR(sum[2], 3.5, EPSILON);  // 3+0.5
    EXPECT_NEAR(sum[3], 6.5, EPSILON);  // 4+2.5

    Quaternion<double> sum_identity = test_q1 + identity;
    EXPECT_NEAR(sum_identity[0], 1.0, EPSILON);
    EXPECT_NEAR(sum_identity[1], 2.0, EPSILON);
    EXPECT_NEAR(sum_identity[2], 3.0, EPSILON);
    EXPECT_NEAR(sum_identity[3], 5.0, EPSILON);  // 4+1

    Quaternion<double> test_copy = test_q1;
    test_copy += test_q2;
    EXPECT_NEAR(test_copy[0], 3.0, EPSILON);
    EXPECT_NEAR(test_copy[1], 3.0, EPSILON);
    EXPECT_NEAR(test_copy[2], 3.5, EPSILON);
    EXPECT_NEAR(test_copy[3], 6.5, EPSILON);
}

TEST_F(QuaternionTest, SubtractionOperators)
{
    Quaternion<double> diff = test_q1 - test_q2;
    EXPECT_NEAR(diff[0], -1.0, EPSILON);  // 1-2
    EXPECT_NEAR(diff[1], 1.0, EPSILON);   // 2-1
    EXPECT_NEAR(diff[2], 2.5, EPSILON);   // 3-0.5
    EXPECT_NEAR(diff[3], 1.5, EPSILON);   // 4-2.5

    Quaternion<double> diff_identity = test_q1 - identity;
    EXPECT_NEAR(diff_identity[0], 1.0, EPSILON);
    EXPECT_NEAR(diff_identity[1], 2.0, EPSILON);
    EXPECT_NEAR(diff_identity[2], 3.0, EPSILON);
    EXPECT_NEAR(diff_identity[3], 3.0, EPSILON);  // 4-1

    Quaternion<double> test_copy = test_q1;
    test_copy -= test_q2;
    EXPECT_NEAR(test_copy[0], -1.0, EPSILON);
    EXPECT_NEAR(test_copy[1], 1.0, EPSILON);
    EXPECT_NEAR(test_copy[2], 2.5, EPSILON);
    EXPECT_NEAR(test_copy[3], 1.5, EPSILON);
}

TEST_F(QuaternionTest, ScalarMultiplication)
{
    Quaternion<double> scaled = 2.0 * test_q1;
    EXPECT_NEAR(scaled[0], 2.0, EPSILON);
    EXPECT_NEAR(scaled[1], 4.0, EPSILON);
    EXPECT_NEAR(scaled[2], 6.0, EPSILON);
    EXPECT_NEAR(scaled[3], 8.0, EPSILON);

    Quaternion<double> zero_scaled = 0.0 * test_q1;
    EXPECT_NEAR(zero_scaled[0], 0.0, EPSILON);
    EXPECT_NEAR(zero_scaled[1], 0.0, EPSILON);
    EXPECT_NEAR(zero_scaled[2], 0.0, EPSILON);
    EXPECT_NEAR(zero_scaled[3], 0.0, EPSILON);

    Quaternion<double> test_copy = test_q1;
    test_copy *= 3.0;
    EXPECT_NEAR(test_copy[0], 3.0, EPSILON);
    EXPECT_NEAR(test_copy[1], 6.0, EPSILON);
    EXPECT_NEAR(test_copy[2], 9.0, EPSILON);
    EXPECT_NEAR(test_copy[3], 12.0, EPSILON);
}

TEST_F(QuaternionTest, QuaternionMultiplicationComprehensive)
{
    Quaternion<double> result_left  = identity * test_q1;
    Quaternion<double> result_right = test_q1 * identity;

    EXPECT_NEAR(result_left[0], test_q1[0], EPSILON);
    EXPECT_NEAR(result_left[1], test_q1[1], EPSILON);
    EXPECT_NEAR(result_left[2], test_q1[2], EPSILON);
    EXPECT_NEAR(result_left[3], test_q1[3], EPSILON);

    EXPECT_NEAR(result_right[0], test_q1[0], EPSILON);
    EXPECT_NEAR(result_right[1], test_q1[1], EPSILON);
    EXPECT_NEAR(result_right[2], test_q1[2], EPSILON);
    EXPECT_NEAR(result_right[3], test_q1[3], EPSILON);

    // Test unit quaternion multiplication (i*j = k, j*k = i, k*i = j)
    Quaternion<double> ij = unit_i * unit_j;
    EXPECT_NEAR(ij[0], 0.0, EPSILON);
    EXPECT_NEAR(ij[1], 0.0, EPSILON);
    EXPECT_NEAR(ij[2], 1.0, EPSILON);  // k
    EXPECT_NEAR(ij[3], 0.0, EPSILON);

    Quaternion<double> jk = unit_j * unit_k;
    EXPECT_NEAR(jk[0], 1.0, EPSILON);  // i
    EXPECT_NEAR(jk[1], 0.0, EPSILON);
    EXPECT_NEAR(jk[2], 0.0, EPSILON);
    EXPECT_NEAR(jk[3], 0.0, EPSILON);

    Quaternion<double> ki = unit_k * unit_i;
    EXPECT_NEAR(ki[0], 0.0, EPSILON);
    EXPECT_NEAR(ki[1], 1.0, EPSILON);  // j
    EXPECT_NEAR(ki[2], 0.0, EPSILON);
    EXPECT_NEAR(ki[3], 0.0, EPSILON);

    Quaternion<double> test_copy = test_q1;
    test_copy *= identity;
    EXPECT_NEAR(test_copy[0], test_q1[0], EPSILON);
    EXPECT_NEAR(test_copy[1], test_q1[1], EPSILON);
    EXPECT_NEAR(test_copy[2], test_q1[2], EPSILON);
    EXPECT_NEAR(test_copy[3], test_q1[3], EPSILON);
}

TEST_F(QuaternionTest, QuaternionVectorMultiplication)
{
    Quaternion<double> result1 = identity * test_vec1;
    EXPECT_NEAR(result1[0], 1.0, EPSILON);
    EXPECT_NEAR(result1[1], 0.0, EPSILON);
    EXPECT_NEAR(result1[2], 0.0, EPSILON);
    EXPECT_NEAR(result1[3], 0.0, EPSILON);

    Quaternion<double> result2 = identity * test_vec3;
    EXPECT_NEAR(result2[0], 1.0, EPSILON);
    EXPECT_NEAR(result2[1], 2.0, EPSILON);
    EXPECT_NEAR(result2[2], 3.0, EPSILON);
    EXPECT_NEAR(result2[3], 0.0, EPSILON);

    Vector3<double> vec_copy = test_vec1;
    vec_copy *= identity;
}

TEST_F(QuaternionTest, VectorRotationOperators)
{
    Vector3<double> rotated_left  = identity << test_vec1;
    Vector3<double> rotated_right = identity >> test_vec1;

    EXPECT_NEAR(rotated_left[0], test_vec1[0], EPSILON);
    EXPECT_NEAR(rotated_left[1], test_vec1[1], EPSILON);
    EXPECT_NEAR(rotated_left[2], test_vec1[2], EPSILON);

    EXPECT_NEAR(rotated_right[0], test_vec1[0], EPSILON);
    EXPECT_NEAR(rotated_right[1], test_vec1[1], EPSILON);
    EXPECT_NEAR(rotated_right[2], test_vec1[2], EPSILON);

    Vector3<double> rotated_norm = normalized_q << test_vec1;

    EXPECT_NEAR(norm(rotated_norm), norm(test_vec1), EPSILON);

    Vector3<double> vec_copy = test_vec1;
    vec_copy <<= identity;
    EXPECT_NEAR(vec_copy[0], test_vec1[0], EPSILON);
    EXPECT_NEAR(vec_copy[1], test_vec1[1], EPSILON);
    EXPECT_NEAR(vec_copy[2], test_vec1[2], EPSILON);
}

TEST_F(QuaternionTest, ComparisonOperators)
{
    Quaternion<double> q_copy = test_q1;
    Quaternion<double> q_different(1.0, 2.0, 3.0, 4.1);

    EXPECT_TRUE(test_q1 == q_copy);
    EXPECT_FALSE(test_q1 == test_q2);
    EXPECT_FALSE(test_q1 == q_different);

    EXPECT_FALSE(test_q1 != q_copy);
    EXPECT_TRUE(test_q1 != test_q2);
    EXPECT_TRUE(test_q1 != q_different);

    EXPECT_TRUE(identity == identity);
    EXPECT_FALSE(identity == test_q1);

    EXPECT_TRUE(zero_q == zero_q);
    EXPECT_FALSE(zero_q == identity);
}

TEST_F(QuaternionTest, EdgeCasesAndProperties)
{
    Quaternion<double> product = normalized_q * normalized_q;
    EXPECT_NEAR(norm(product), 1.0, EPSILON);

    Quaternion<double> left_assoc  = (unit_i * unit_j) * unit_k;
    Quaternion<double> right_assoc = unit_i * (unit_j * unit_k);
    EXPECT_NEAR(left_assoc[0], right_assoc[0], EPSILON);
    EXPECT_NEAR(left_assoc[1], right_assoc[1], EPSILON);
    EXPECT_NEAR(left_assoc[2], right_assoc[2], EPSILON);
    EXPECT_NEAR(left_assoc[3], right_assoc[3], EPSILON);

    Quaternion<double>       prod       = test_q1 * test_q2;
    const Quaternion<double> const_prod = prod;
    Quaternion<double>       conj_prod  = conjugate(const_prod);

    const Quaternion<double> const_test_q1 = test_q1;
    const Quaternion<double> const_test_q2 = test_q2;
    Quaternion<double>       conj_q1       = conjugate(const_test_q1);
    Quaternion<double>       conj_q2       = conjugate(const_test_q2);
    Quaternion<double>       conj_rev      = conj_q2 * conj_q1;

    EXPECT_NEAR(conj_prod[0], conj_rev[0], EPSILON);
    EXPECT_NEAR(conj_prod[1], conj_rev[1], EPSILON);
    EXPECT_NEAR(conj_prod[2], conj_rev[2], EPSILON);
    EXPECT_NEAR(conj_prod[3], conj_rev[3], EPSILON);

    double norm_prod  = norm(test_q1 * test_q2);
    double prod_norms = norm(test_q1) * norm(test_q2);
    EXPECT_NEAR(norm_prod, prod_norms, EPSILON);
}
