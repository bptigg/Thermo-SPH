#include "kernel.h"
#include <cmath>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

double CubicSplineKernel::value(const Vector2D& r, double h) const {
    double dist = r.norm();
    double q = dist / h;
    if (q >= 2.0) return 0.0;

    double sigma = 10.0 / (7.0 * M_PI * h * h);
    if (q < 1.0) {
        return sigma * (1.0 - 1.5 * q * q + 0.75 * q * q * q);
    }
    double term = 2.0 - q;
    return sigma * 0.25 * term * term * term;
}

Vector2D CubicSplineKernel::gradient(const Vector2D& r, double h) const {
    double dist = r.norm();
    if (dist <= 1e-12) return Vector2D(0.0, 0.0);

    double q = dist / h;
    if (q >= 2.0) return Vector2D(0.0, 0.0);

    double sigma = 10.0 / (7.0 * M_PI * h * h);
    double dWdq = (q < 1.0) ? sigma * (-3.0 * q + 2.25 * q * q)
                            : sigma * (-0.75 * (2.0 - q) * (2.0 - q));

    return r * (dWdq / (h * dist));
}

double WendlandC2Kernel::value(const Vector2D& r, double h) const {
    double dist = r.norm();
    double q = dist / h;
    if (q >= 2.0) return 0.0;

    double sigma = 7.0 / (4.0 * M_PI * h * h);
    double term = 1.0 - 0.5 * q;
    return sigma * std::pow(term, 4) * (1.0 + 2.0 * q);
}

Vector2D WendlandC2Kernel::gradient(const Vector2D& r, double h) const {
    double dist = r.norm();
    if (dist <= 1e-12) return Vector2D(0.0, 0.0);

    double q = dist / h;
    if (q >= 2.0) return Vector2D(0.0, 0.0);

    double sigma = 7.0 / (4.0 * M_PI * h * h);
    double term = 1.0 - 0.5 * q;
    double dWdq = sigma * (-5.0 * q * std::pow(term, 3));

    return r * (dWdq / (h * dist));
}