#pragma once

#include "Particle.h"
#include "SpatialGrid.h"
#include "SystemVolatility.h"
#include "ThreadPool.h"
#include <vector>
#include <memory>

class VelocityVerletIntegrator {
private:
    double h_;
    double gamma_;
    double cflSafety_;
    double maxDt_;
    double minDt_;
    size_t sampleSize_;
    SpatialGrid spatialGrid_;

public:
    VelocityVerletIntegrator(
        double h,
        double gamma = 1.4,
        double cflSafety = 0.25,
        double maxDt = 1e-3,
        double minDt = 1e-7,
        size_t sampleSize = 128)
        : h_(h),
          gamma_(gamma),
          cflSafety_(cflSafety),
          maxDt_(maxDt),
          minDt_(minDt),
          sampleSize_(sampleSize),
          spatialGrid_(2.0 * h) {}

    // Parallelized Integration Steps
    void kickFirstHalf(std::vector<std::unique_ptr<Particle>>& particles, double dt, ThreadPool& pool) const;
    void drift(std::vector<std::unique_ptr<Particle>>& particles, double dt, ThreadPool& pool) const;
    void kickSecondHalf(std::vector<std::unique_ptr<Particle>>& particles, double dt, ThreadPool& pool) const;

    // Single-threaded O(1) Volatility Calculation (128-particle strided sample)
    double computeAdaptiveTimestep(
        const std::vector<std::unique_ptr<Particle>>& particles,
        double currentDt) const;

    SpatialGrid& getSpatialGrid() { return spatialGrid_; }
    const SpatialGrid& getSpatialGrid() const { return spatialGrid_; }
};