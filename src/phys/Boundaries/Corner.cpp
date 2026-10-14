#include "Corner.h"
#include "Particle.h"
#include "SpatialGrid.h"
#include <cmath>

CornerBoundary::CornerBoundary(Vector2D normalA, Vector2D normalB, Vector2D corner, bool isNoSlip)
    : normalA_(normalA / normalA.norm()), normalB_(normalB / normalB.norm()),
      corner_(corner), isNoSlip_(isNoSlip) {}

std::vector<std::shared_ptr<Particle>> CornerBoundary::generateGhosts(
    const Vector2D& hostPos,
    const Vector2D& hostVel,
    double /*hostEnergy*/,
    const std::vector<std::shared_ptr<Particle>>& fluidParticles,
    double supportRadius,
    const SpatialGrid& grid) const
{
    std::vector<std::shared_ptr<Particle>> ghosts;
    int ghostID = -1000000;

    for (size_t idx : grid.getNeighborIndices(hostPos)) {
        if (idx >= fluidParticles.size()) continue;
        const auto& fluid = fluidParticles[idx];
        if (!fluid->isFluid()) continue;

        Vector2D rel = fluid->pos - corner_;
        double sA = rel.dot(normalA_);
        double sB = rel.dot(normalB_);
        if (sA >= supportRadius || sB >= supportRadius) continue;
        if (sA < 1e-5 || sB < 1e-5) continue;

        Vector2D ghostPos = fluid->pos - normalA_ * (2.0 * sA) - normalB_ * (2.0 * sB);
        Vector2D ghostVel = isNoSlip_ ? hostVel * 2.0 - fluid->vel : fluid->vel * -1.0;

        auto ghost = std::make_shared<FluidParticle>(
            ghostID--, ghostPos, ghostVel, fluid->mass, fluid->density, fluid->u);
        ghost->h = fluid->h;
        ghosts.push_back(std::move(ghost));
    }
    return ghosts;
}
