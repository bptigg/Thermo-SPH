#pragma once
#include "IBoundary.h"

class PlanarBoundary : public IBoundary {
private:
    Vector2D normal_; 
    bool isNoSlip_;
    double spacing_;
    double surfaceOffset_;

public:
    // spacing > 0: this host owns a wall segment of that width, so each fluid particle is mirrored once.
    // surfaceOffset: distance from the host to the wall surface along the normal.
    PlanarBoundary(Vector2D normal, bool isNoSlip = true, double spacing = 0.0, double surfaceOffset = 0.0);

    std::vector<std::shared_ptr<Particle>> generateGhosts(
        const Vector2D& hostPos,
        const Vector2D& hostVel,
        double hostEnergy,
        const std::vector<std::shared_ptr<Particle>>& fluidParticles,
        double supportRadius,
        const SpatialGrid& grid) const override;
};