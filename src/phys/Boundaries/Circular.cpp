#include "Circular.h"
#include "Particle.h"

std::vector<std::shared_ptr<Particle>> CircularBoundary::generateGhosts(
    const Vector2D& hostPos,
    const Vector2D& hostVel,
    double hostEnergy,
    const std::vector<std::shared_ptr<Particle>>& fluidParticles,
    double supportRadius) const 
{
    std::vector<std::shared_ptr<Particle>> ghosts;
    int ghostID = -1;

    for (const auto& fluid : fluidParticles) {
        if (!fluid->isFluid()) continue;

        Vector2D relPos = fluid->pos - hostPos;
        double distFromCenter = relPos.norm();
        double distToSurface = distFromCenter - radius_;

        if (distToSurface > 0.0 && distToSurface < supportRadius && distFromCenter > 1e-12) {
            Vector2D normal = relPos / distFromCenter;
            Vector2D ghostPos = fluid->pos - normal * (2.0 * distToSurface);

            Vector2D ghostVel = isNoSlip_ 
                ? (hostVel * 2.0 - fluid->vel)
                : (fluid->vel - normal * (2.0 * (fluid->vel - hostVel).dot(normal)));

            ghosts.push_back(std::make_unique<FluidParticle>(
                ghostID--, ghostPos, ghostVel, fluid->mass, fluid->density, hostEnergy
            ));
        }
    }

    return ghosts;
}