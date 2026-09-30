#pragma once
#include "Particle.h"
#include <vector>
#include <memory>

class GhostManager {

public:
    static std::vector<std::shared_ptr<Particle>> generateAllGhosts(
        const std::vector<std::shared_ptr<Particle>>& worldParticles,
        const std::vector<std::shared_ptr<Particle>>& FluidParticles,
        double supportRadius, const SpatialGrid& grid, ThreadPool& pool) 
    {

        size_t n = worldParticles.size();
        if (n == 0) return {};
        size_t chunkSize = 128;

        size_t numChunks = (n + chunkSize - 1) / chunkSize;
        std::vector<std::vector<std::shared_ptr<Particle>>> chunkGhosts(numChunks);

        pool.parallel_for(0, n, [&worldParticles, &FluidParticles, &grid, &chunkGhosts, supportRadius, chunkSize](size_t start, size_t end) {
            std::vector<std::shared_ptr<Particle>> localGhosts;
            size_t chunkIdx = start / chunkSize; // Unique to each parallel job

            for (size_t i = start; i < end; ++i) {
                const auto& p = worldParticles[i];
                if (p->hasBoundary()) {
                    auto ghosts = p->generateGhosts(FluidParticles, supportRadius, grid);
                    for (auto& g : ghosts) {
                        //chunkGhosts[chunkIdx].push_back(std::move(g));
                        localGhosts.push_back(std::move(g));
                    }
                }
            }

            chunkGhosts[chunkIdx] = std::move(localGhosts);
        }, chunkSize);

        std::vector<std::shared_ptr<Particle>> allGhosts;
        for (auto& local : chunkGhosts) {
            allGhosts.insert(allGhosts.end(), 
                               std::make_move_iterator(local.begin()), 
                               std::make_move_iterator(local.end()));
        }

        return allGhosts;

        //for (const auto& p : worldParticles) {
        //    if (p->hasBoundary()) {
        //        auto ghosts = p->generateGhosts(FluidParticles, supportRadius, grid);
        //        for (auto& g : ghosts) {
        //            allGhosts.push_back(std::move(g));
        //        }
        //    }
        //}
//
        //return allGhosts;
    }
};