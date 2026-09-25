#pragma once
#include "IBoundary.h"

class PlanarBoundary : public IBoundary {
private:
    Vector2D normal_; 
    bool isNoSlip_;

public:
    PlanarBoundary(Vector2D normal, bool isNoSlip = true);

    std::vector<std::shared_ptr<Particle>> generateGhosts(
        const Vector2D& hostPos,
        const Vector2D& hostVel,
        double hostEnergy,
        const std::vector<std::shared_ptr<Particle>>& fluidParticles,
        double supportRadius) const override;
};