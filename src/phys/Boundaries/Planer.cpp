#include "Planer.h"
#include "Particle.h"

PlanarBoundary::PlanarBoundary(Vector2D normal, bool isNoSlip)
    : isNoSlip_(isNoSlip) 
{
    double len = normal.norm();
    normal_ = (len > 0.0) ? normal / len : Vector2D(0.0, 1.0);
}

std::vector<std::unique_ptr<Particle>> PlanarBoundary::generateGhosts(
    const Vector2D& hostPos,
    const Vector2D& hostVel,
    double hostEnergy,
    const std::vector<std::unique_ptr<Particle>>& fluidParticles,
    double supportRadius) const 
{
    std::vector<std::unique_ptr<Particle>> ghosts;
    int ghostID = -1;

    for (const auto& fluid : fluidParticles) {
        if (!fluid->isFluid()) continue;

        Vector2D relPos = fluid->pos - hostPos;
        double dist = relPos.dot(normal_);

        if (dist > 0.0 && dist < supportRadius) {
            Vector2D ghostPos = fluid->pos - normal_ * (2.0 * dist);
            Vector2D ghostVel;

            if (isNoSlip_) {
                ghostVel = hostVel * 2.0 - fluid->vel;
            } else {
                Vector2D vRel = fluid->vel - hostVel;
                double vRelNorm = vRel.dot(normal_);
                ghostVel = fluid->vel - normal_ * (2.0 * vRelNorm);
            }

            ghosts.push_back(std::make_unique<FluidParticle>(
                ghostID--, ghostPos, ghostVel, fluid->mass, fluid->density, hostEnergy
            ));
        }
    }

    return ghosts;
}