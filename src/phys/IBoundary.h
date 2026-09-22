#pragma once
#include <vector>
#include <memory>

class Particle;

class IBoundary {
public:
    virtual ~IBoundary() = default;

    // Evaluates proximity to fluid particles and returns boundary ghost particles
    virtual std::vector<std::unique_ptr<Particle>> generateGhosts(
        const Vector2D& hostPos,
        const Vector2D& hostVel,
        double hostEnergy,
        const std::vector<std::unique_ptr<Particle>>& fluidParticles,
        double supportRadius) const = 0;
};