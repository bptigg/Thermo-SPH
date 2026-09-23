#include "SystemVolatility.h"
#include <cmath>
#include <algorithm>

VolatilityMetrics SystemVolatility::evaluate(
    const std::vector<std::unique_ptr<Particle>>& particles,
    double h,
    double gamma,
    size_t sampleSize) 
{
    VolatilityMetrics metrics;
    size_t totalParticles = particles.size();
    if (totalParticles == 0) return metrics;

    // Determine stride step: evaluate all if below sample threshold
    size_t stride = (totalParticles > sampleSize) ? (totalParticles / sampleSize) : 1;

    for (size_t i = 0; i < totalParticles; i += stride) {
        const auto& p = particles[i];
        if (!p->isDynamic()) continue;

        // 1. Sound speed c_s = sqrt(gamma * (gamma - 1) * u)
        double c_sound = std::sqrt(gamma * (gamma - 1.0) * std::max(p->u, 1e-6));
        double v_norm  = p->vel.norm();
        double signal  = c_sound + v_norm;

        if (signal > metrics.maxSignalSpeed) {
            metrics.maxSignalSpeed = signal;
        }

        // 2. Acceleration magnitude
        double a_norm = p->accel.norm();
        if (a_norm > metrics.maxAccel) {
            metrics.maxAccel = a_norm;
        }

        // 3. Thermal volatility |du/dt| / (u + eps)
        double thermalRate = std::abs(p->dudt) / (p->u + 1e-4);
        if (thermalRate > metrics.maxRelativeEnergyRate) {
            metrics.maxRelativeEnergyRate = thermalRate;
        }
    }

    // Composite Volatility metric (inverse time scale [1/s])
    double omegaCFL     = metrics.maxSignalSpeed / h;
    double omegaForce   = std::sqrt(metrics.maxAccel / h);
    double omegaThermal = metrics.maxRelativeEnergyRate;

    metrics.compositeVolatility = std::max({omegaCFL, omegaForce, omegaThermal});

    return metrics;
}