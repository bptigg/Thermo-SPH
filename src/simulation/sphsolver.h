#pragma once

#include <vector>
#include <memory>
#include "Particle.h"
#include "Vector2D.h"
#include "kernel.h"
#include "SpatialGrid.h"

class SPHSolver {
public:
    struct Parameters {
        double gamma = 1.4;
        double alpha = 1.0; 
        double beta  = 2.0; 
        double eta   = 1.3; // Kernel scale factor: h_i = eta * sqrt(m_i / rho_i)
    };

    explicit SPHSolver(Parameters params);

    static double computeSoundSpeed(double u, double gamma);

    void updateSmoothingLengths(std::vector<std::shared_ptr<Particle>>& particles) const;

    void computeDensityAndPressure(
        std::vector<std::shared_ptr<Particle>>& particles,
        const std::vector<std::shared_ptr<Particle>>& ghosts,
        const SpatialGrid& grid,
        const Kernel& kernel) const;

    void computeDerivatives(
        std::vector<std::shared_ptr<Particle>>& particles,
        const std::vector<std::shared_ptr<Particle>>& ghosts,
        const SpatialGrid& grid,
        const Kernel& kernel,
        Vector2D& netBoundaryForce) const;

private:
    Parameters params_;
};