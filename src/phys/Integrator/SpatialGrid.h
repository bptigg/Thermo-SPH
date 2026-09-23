#pragma once
#include <vector>
#include <memory>
#include <unordered_map>
#include <cstdint>
#include "Vector2D.h"

class Particle;

class SpatialGrid {
private:
    double cellSize_;
    
    int64_t hashCell(int cellX, int cellY) const {
        uint64_t ux = static_cast<uint32_t>(cellX);
        uint64_t uy = static_cast<uint32_t>(cellY);
        return static_cast<int64_t>((ux << 32) | uy);
    }

    std::unordered_map<int64_t, std::vector<size_t>> grid_;

public:
    explicit SpatialGrid(double cellSize) : cellSize_(cellSize) {}

    void build(const std::vector<std::unique_ptr<Particle>>& particles);

    std::vector<size_t> getNeighborIndices(const Vector2D& pos) const {
        std::vector<size_t> neighbors;
        int centerCellX = static_cast<int>(std::floor(pos.x / cellSize_));
        int centerCellY = static_cast<int>(std::floor(pos.y / cellSize_));

        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                int64_t key = hashCell(centerCellX + dx, centerCellY + dy);
                auto it = grid_.find(key);
                if (it != grid_.end()) {
                    neighbors.insert(neighbors.end(), it->second.begin(), it->second.end());
                }
            }
        }
        return neighbors;
    }
};