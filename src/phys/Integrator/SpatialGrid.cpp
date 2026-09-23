#include "SpatialGrid.h"
#include "Particle.h"

void SpatialGrid::build(const std::vector<std::unique_ptr<Particle>> &particles)
{
    grid_.clear();
    for (size_t i = 0; i < particles.size(); ++i) {
        int cellX = static_cast<int>(std::floor(particles[i]->pos.x / cellSize_));
        int cellY = static_cast<int>(std::floor(particles[i]->pos.y / cellSize_));
        grid_[hashCell(cellX, cellY)].push_back(i);
    }
}
