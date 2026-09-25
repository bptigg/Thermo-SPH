#include "sod_shock_tube.h"
#include <cmath>

SodShockTubeIC::SodShockTubeIC() : SodShockTubeIC(Parameters{}) {}
SodShockTubeIC::SodShockTubeIC(Parameters params) : params_(params) {}

std::vector<std::shared_ptr<Particle>> SodShockTubeIC::generateParticles() {
    std::vector<std::shared_ptr<Particle>> particles;
    int id = 0;

    double dxLeft = params_.spacing;
    double dxRight = params_.spacing * std::pow(params_.rhoLeft / params_.rhoRight, 0.5);
    double particleMass = params_.rhoLeft * dxLeft * dxLeft;

    for (double x = params_.xMin; x <= params_.xMax; x += (x < 0 ? dxLeft : dxRight)) {
        for (double y = params_.yMin; y <= params_.yMax; y += (x < 0 ? dxLeft : dxRight)) {
            bool isLeft = (x < 0.0);
            double rho = isLeft ? params_.rhoLeft : params_.rhoRight;
            double u   = isLeft ? params_.uLeft   : params_.uRight;

            auto p = std::make_shared<FluidParticle>(
                id++,
                Vector2D(x, y),
                Vector2D(0.0, 0.0),
                particleMass,
                rho,
                u
            );
            particles.push_back(p);
        }
    }
    return particles;
}


