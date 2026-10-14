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
        double cellSize = grid.getCellSize();

        pool.parallel_for(0, n, [&worldParticles, &FluidParticles, &grid, &chunkGhosts, supportRadius, chunkSize, cellSize](size_t start, size_t end) {
            std::vector<std::shared_ptr<Particle>> localGhosts;
            size_t chunkIdx = start / chunkSize; // Unique to each parallel job

            for (size_t i = start; i < end; ++i) {
                const auto& p = worldParticles[i];
                if(!p || !p->hasBoundary()) continue;

                int cellX = static_cast<int>(std::floor(p->pos.x / cellSize));
                int cellY = static_cast<int>(std::floor(p->pos.y / cellSize));

                bool nearFluid = false;
                for (int dx = -1; dx <= 1; ++dx) {
                    for (int dy = -1; dy <= 1; ++dy) {
                        if (grid.isCellActive(cellX + dx, cellY + dy)) {
                            nearFluid = true;
                            break;
                        }
                    }
                    if (nearFluid) break;
                }

                if (!nearFluid) continue;

                auto ghosts = p->generateGhosts(FluidParticles, supportRadius, grid);
                for (auto& g : ghosts) {
                    localGhosts.push_back(std::move(g));
                }

                //if (p->hasBoundary()) {
                //    auto ghosts = p->generateGhosts(FluidParticles, supportRadius, grid);
                //    for (auto& g : ghosts) {
                //        //chunkGhosts[chunkIdx].push_back(std::move(g));
                //        localGhosts.push_back(std::move(g));
                //    }
                //}
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