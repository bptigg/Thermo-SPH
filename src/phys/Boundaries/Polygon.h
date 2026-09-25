#pragma once
#include "IBoundary.h"
#include <vector>

class PolygonBoundary : public IBoundary {
private:
    std::vector<Vector2D> localVertices_; 
    bool isNoSlip_;

public:
    PolygonBoundary(const std::vector<Vector2D>& localVertices, bool isNoSlip = true);

    std::vector<std::shared_ptr<Particle>> generateGhosts(
        const Vector2D& hostPos,
        const Vector2D& hostVel,
        double hostEnergy,
        const std::vector<std::shared_ptr<Particle>>& fluidParticles,
        double supportRadius) const override;
};