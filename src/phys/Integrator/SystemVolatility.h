#pragma once
#include "Particle.h"
#include <vector>
#include <memory>
#include <cstddef>

struct VolatilityMetrics {
    double maxSignalSpeed = 0.0;
    double maxAccel = 0.0;
    double maxRelativeEnergyRate = 0.0;
    double compositeVolatility = 0.0;
};

class SystemVolatility {
public:
    static VolatilityMetrics evaluate(
        const std::vector<std::unique_ptr<Particle>>& particles,
        double h,
        double gamma = 1.4,
        size_t sampleSize = 128); // Target sample count
};

