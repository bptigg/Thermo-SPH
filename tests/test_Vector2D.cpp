#include "Vector2D.h"
#include <cassert>
#include <cmath>
#include <iostream>

constexpr double EPSILON = 1e-9;

bool approxEqual(double a, double b) {
    return std::abs(a - b) < EPSILON;
}

void testAdditionAndSubtraction() {
    Vector2D v1(2.0, 3.0);
    Vector2D v2(4.0, -1.0);

    Vector2D sum = v1 + v2;
    assert(approxEqual(sum.x, 6.0));
    assert(approxEqual(sum.y, 2.0));

    Vector2D diff = v1 - v2;
    assert(approxEqual(diff.x, -2.0));
    assert(approxEqual(diff.y, 4.0));

    std::cout << " [PASS] testAdditionAndSubtraction" << std::endl;
}

void testScalarMultiplicationAndDivision() {
    Vector2D v(3.0, -4.0);

    Vector2D scaled = v * 2.5;
    assert(approxEqual(scaled.x, 7.5));
    assert(approxEqual(scaled.y, -10.0));

    Vector2D divided = v / 2.0;
    assert(approxEqual(divided.x, 1.5));
    assert(approxEqual(divided.y, -2.0));

    std::cout << " [PASS] testScalarMultiplicationAndDivision" << std::endl;
}

void testDotProductAndNorm() {
    Vector2D v1(3.0, 4.0);
    Vector2D v2(1.0, 2.0);

    assert(approxEqual(v1.dot(v2), 11.0));
    assert(approxEqual(v1.normSq(), 25.0));
    assert(approxEqual(v1.norm(), 5.0));

    Vector2D normV1 = v1.normalized();
    assert(approxEqual(normV1.norm(), 1.0));
    assert(approxEqual(normV1.x, 0.6));
    assert(approxEqual(normV1.y, 0.8));

    std::cout << " [PASS] testDotProductAndNorm" << std::endl;
}

void testZeroVectorNormalization() {
    Vector2D zero(0.0, 0.0);
    Vector2D normZero = zero.normalized();

    assert(approxEqual(normZero.x, 0.0));
    assert(approxEqual(normZero.y, 0.0));

    std::cout << " [PASS] testZeroVectorNormalization" << std::endl;
}

int main() {
    std::cout << "=== Running Vector2D Unit Tests ===" << std::endl;
    testAdditionAndSubtraction();
    testScalarMultiplicationAndDivision();
    testDotProductAndNorm();
    testZeroVectorNormalization();
    std::cout << "All Vector2D tests passed successfully!\n" << std::endl;
    return 0;
}