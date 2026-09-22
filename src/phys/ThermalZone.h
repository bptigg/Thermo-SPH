#pragma once
#include "Vector2D.h"
#include "Particle.h"
#include <vector>
#include <memory>

class ThermalZone {
public:
    Vector2D minBounds;
    Vector2D maxBounds;
    double energyRate; // Positive = Heating, Negative = Siphoning/Cooling

    ThermalZone(Vector2D minB, Vector2D maxB, double rate)
        : minBounds(minB), maxBounds(maxB), energyRate(rate) {}

    bool contains(const Vector2D& pos) const {
        return (pos.x >= minBounds.x && pos.x <= maxBounds.x &&
                pos.y >= minBounds.y && pos.y <= maxBounds.y);
    }

    void apply(std::vector<std::unique_ptr<Particle>>& particles) const {
        for (auto& p : particles) {
            if (p->isFluid() && contains(p->pos)) {
                p->dudt += energyRate / p->mass;
            }
        }
    }
};