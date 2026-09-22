#include "test_Vector2D.h"
#include "Vector2D.h"
#include <cmath>

static bool testAdditionAndSubtraction() {
    Vector2D v1(2.0, 3.0);
    Vector2D v2(4.0, -1.0);

    Vector2D sum = v1 + v2;
    TEST_ASSERT_NEAR(sum.x, 6.0, 1e-9, "Vector addition X failed");
    TEST_ASSERT_NEAR(sum.y, 2.0, 1e-9, "Vector addition Y failed");

    Vector2D diff = v1 - v2;
    TEST_ASSERT_NEAR(diff.x, -2.0, 1e-9, "Vector subtraction X failed");
    TEST_ASSERT_NEAR(diff.y, 4.0, 1e-9, "Vector subtraction Y failed");

    return true;
}

static bool testScalarMultiplicationAndDivision() {
    Vector2D v(3.0, -4.0);

    Vector2D scaled = v * 2.5;
    TEST_ASSERT_NEAR(scaled.x, 7.5, 1e-9, "Scalar multiplication X failed");
    TEST_ASSERT_NEAR(scaled.y, -10.0, 1e-9, "Scalar multiplication Y failed");

    Vector2D divided = v / 2.0;
    TEST_ASSERT_NEAR(divided.x, 1.5, 1e-9, "Scalar division X failed");
    TEST_ASSERT_NEAR(divided.y, -2.0, 1e-9, "Scalar division Y failed");

    return true;
}

static bool testDotProductAndNorm() {
    Vector2D v1(3.0, 4.0);
    Vector2D v2(1.0, 2.0);

    TEST_ASSERT_NEAR(v1.dot(v2), 11.0, 1e-9, "Dot product failed");
    TEST_ASSERT_NEAR(v1.normSq(), 25.0, 1e-9, "Squared norm failed");
    TEST_ASSERT_NEAR(v1.norm(), 5.0, 1e-9, "Norm failed");

    Vector2D normV1 = v1.normalized();
    TEST_ASSERT_NEAR(normV1.norm(), 1.0, 1e-9, "Normalized vector length failed");
    TEST_ASSERT_NEAR(normV1.x, 0.6, 1e-9, "Normalized X failed");
    TEST_ASSERT_NEAR(normV1.y, 0.8, 1e-9, "Normalized Y failed");

    return true;
}

static bool testZeroVectorNormalization() {
    Vector2D zero(0.0, 0.0);
    Vector2D normZero = zero.normalized();

    TEST_ASSERT_NEAR(normZero.x, 0.0, 1e-9, "Zero vector normalized X failed");
    TEST_ASSERT_NEAR(normZero.y, 0.0, 1e-9, "Zero vector normalized Y failed");

    return true;
}

void registerVector2DTests(TestSuite& suite) {
    suite.startSection("Vector2D Module");
    suite.runTest("Addition and Subtraction", testAdditionAndSubtraction);
    suite.runTest("Scalar Multiplication and Division", testScalarMultiplicationAndDivision);
    suite.runTest("Dot Product and Norm", testDotProductAndNorm);
    suite.runTest("Zero Vector Normalization", testZeroVectorNormalization);
}