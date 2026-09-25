#pragma once
#include "Particle.h"
#include "SpatialGrid.h"
#include "kernel.h"
#include <vector>
#include <memory>

struct Parameters {
    double gamma = 1.4;
    double alpha = 1.0;
    double beta = 1.0;
    double eta = 1.2;
    double fluidDensity = 1.0;

    bool enableGravity = false;
    Vector2D gravity = Vector2D(0.0, -9.81);
};

class SPHSolver {
public:
    explicit SPHSolver(Parameters params);

    void updateSmoothingLengths(std::vector<std::shared_ptr<Particle>>& fluidParticles) const;

    void computeDensityAndPressure(
        std::vector<std::shared_ptr<Particle>>& fluidParticles,
        const std::vector<std::shared_ptr<Particle>>& internalRigidParticles,
        const std::vector<std::shared_ptr<Particle>>& ghostParticles,
        const SpatialGrid& grid,
        const Kernel& kernel) const;

    void computeDerivatives(
        std::vector<std::shared_ptr<Particle>>& fluidParticles,
        const std::vector<std::shared_ptr<Particle>>& internalRigidParticles,
        const std::vector<std::shared_ptr<Particle>>& ghostParticles,
        const SpatialGrid& grid,
        const Kernel& kernel) const;

    void setParameters(const Parameters& params) { params_ = params; }
    const Parameters& getParameters() const { return params_; }

private:
    Parameters params_;
    double computeSoundSpeed(double u, double gamma) const;
};