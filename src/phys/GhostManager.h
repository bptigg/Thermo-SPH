#pragma once
#include "Particle.h"
#include <vector>
#include <memory>

class GhostManager {
public:
    static std::vector<std::shared_ptr<Particle>> generateAllGhosts(
        const std::vector<std::shared_ptr<Particle>>& worldParticles,
        const std::vector<std::shared_ptr<Particle>>& FluidParticles,
        double supportRadius) 
    {
        std::vector<std::shared_ptr<Particle>> allGhosts;

        for (const auto& p : worldParticles) {
            if (p->hasBoundary()) {
                auto ghosts = p->generateGhosts(FluidParticles, supportRadius);
                for (auto& g : ghosts) {
                    allGhosts.push_back(std::move(g));
                }
            }
        }

        return allGhosts;
    }
};