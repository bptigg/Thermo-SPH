#include "SpatialGrid.h"
#include "Particle.h"

void SpatialGrid::build(const std::vector<std::shared_ptr<Particle>> &particles)
{
    grid_.clear();
    occupiedCells_.clear();
    for (size_t i = 0; i < particles.size(); ++i) {
        int cellX = static_cast<int>(std::floor(particles[i]->pos.x / cellSize_));
        int cellY = static_cast<int>(std::floor(particles[i]->pos.y / cellSize_));
        int64_t key = hashCell(cellX, cellY);

        grid_[key].push_back(i);
        occupiedCells_.insert(key);
    }
}

void SpatialGrid::build(const std::vector<Particle*> &particles)
{
    grid_.clear();
    occupiedCells_.clear();

    for (size_t i = 0; i < particles.size(); ++i) {
        if (!particles[i]) continue;
        int cellX = static_cast<int>(std::floor(particles[i]->pos.x / cellSize_));
        int cellY = static_cast<int>(std::floor(particles[i]->pos.y / cellSize_));
        int64_t key = hashCell(cellX, cellY);

        grid_[key].push_back(i);
        occupiedCells_.insert(key);
    }
}
