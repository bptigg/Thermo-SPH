#pragma once
#include "IBoundary.h"

// Concave corner formed by two perpendicular walls; mirrors fluid through the corner point
// (the diagonal image that the two planar boundaries do not produce).
class CornerBoundary : public IBoundary {
private:
    Vector2D normalA_;
    Vector2D normalB_;
    Vector2D corner_;
    bool isNoSlip_;

public:
    // normalA/normalB point into the fluid; corner is the point where the two wall surfaces meet.
    CornerBoundary(Vector2D normalA, Vector2D normalB, Vector2D corner, bool isNoSlip = true);

    std::vector<std::shared_ptr<Particle>> generateGhosts(
        const Vector2D& hostPos,
        const Vector2D& hostVel,
        double hostEnergy,
        const std::vector<std::shared_ptr<Particle>>& fluidParticles,
        double supportRadius,
        const SpatialGrid& grid) const override;
};
