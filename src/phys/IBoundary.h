#pragma once
#include <vector>
#include <memory>
#include "Vector2D.h"

class Particle;

class IBoundary {
public:
    virtual ~IBoundary() = default;

    // Evaluates proximity to fluid particles and returns boundary ghost particles
    virtual std::vector<std::shared_ptr<Particle>> generateGhosts(
        const Vector2D& hostPos,
        const Vector2D& hostVel,
        double hostEnergy,
        const std::vector<std::shared_ptr<Particle>>& fluidParticles,
        double supportRadius) const = 0;
};