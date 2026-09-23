#pragma once
#include "IBoundary.h"

class CircularBoundary : public IBoundary {
private:
    double radius_;
    bool isNoSlip_;

public:
    CircularBoundary(double radius, bool isNoSlip = true)
        : radius_(radius), isNoSlip_(isNoSlip) {}

    std::vector<std::unique_ptr<Particle>> generateGhosts(
        const Vector2D& hostPos,
        const Vector2D& hostVel,
        double hostEnergy,
        const std::vector<std::unique_ptr<Particle>>& fluidParticles,
        double supportRadius) const override;
};