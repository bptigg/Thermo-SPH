#pragma once
#include "Vector2D.h"

enum class KernelType {
    CubicSpline,
    WendlandC2
};

class Kernel {
public:
    virtual ~Kernel() = default;

    // Kernel evaluation W(r, h)
    virtual double value(const Vector2D& r, double h) const = 0;

    // Gradient evaluation grad W(r, h)
    virtual Vector2D gradient(const Vector2D& r, double h) const = 0;

    // Cutoff radius factor (e.g., 2.0 * h)
    virtual double cutoffFactor() const = 0;
    virtual KernelType type() const = 0; // GPU dispatcher query
};

class CubicSplineKernel : public Kernel {
public:
    double value(const Vector2D& r, double h) const override;
    Vector2D gradient(const Vector2D& r, double h) const override;
    double cutoffFactor() const override { return 2.0; }
    KernelType type() const override { return KernelType::CubicSpline; }
};

class WendlandC2Kernel : public Kernel {
public:
    double value(const Vector2D& r, double h) const override;
    Vector2D gradient(const Vector2D& r, double h) const override;
    double cutoffFactor() const override { return 2.0; }
    KernelType type() const override { return KernelType::CubicSpline; }
};

