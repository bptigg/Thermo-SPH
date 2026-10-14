#pragma once
#include <vector>
#include <memory>
#include <unordered_map>
#include <unordered_set>
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
    std::unordered_set<int64_t> occupiedCells_;

public:
    explicit SpatialGrid(double cellSize = 0.2) : cellSize_(cellSize) {}

    void build(const std::vector<std::shared_ptr<Particle>>& particles);
    void build(const std::vector<Particle*>& particles);

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

    const std::unordered_set<int64_t>& getActiveCells() const {
        return occupiedCells_;
    }

    bool isCellActive(int cellX, int cellY) const {
        return occupiedCells_.find(hashCell(cellX, cellY)) != occupiedCells_.end();
    }

    static void decodeCell(int64_t key, int& cellX, int& cellY) {
        uint64_t ukey = static_cast<uint64_t>(key);
        cellX = static_cast<int32_t>(ukey >> 32);
        cellY = static_cast<int32_t>(ukey & 0xFFFFFFFF);
    }

    double getCellSize() const { return cellSize_; }

};