#include "thermal_equilibrium.h"

ThermalEquilibriumIC::ThermalEquilibriumIC() : ThermalEquilibriumIC(Parameters{}) {}
ThermalEquilibriumIC::ThermalEquilibriumIC(Parameters params) : params_(params) {}

std::vector<std::shared_ptr<Particle>> ThermalEquilibriumIC::generateParticles() {
    std::vector<std::shared_ptr<Particle>> particles;
    int id = 0;

    double particleMass = params_.density * params_.spacing * params_.spacing;

    for (double x = params_.domainMin.x; x <= params_.domainMax.x; x += params_.spacing) {
        for (double y = params_.domainMin.y; y <= params_.domainMax.y; y += params_.spacing) {
            double u = (x < 0.0) ? params_.uHot : params_.uCold;

            auto p = std::make_shared<FluidParticle>(
                id++,
                Vector2D(x, y),
                Vector2D(0.0, 0.0),
                particleMass,
                params_.density,
                u
            );
            particles.push_back(p);
        }
    }
    return particles;
}