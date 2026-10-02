#pragma once
#include "Particle.h"
#include "SpatialGrid.h"
#include "kernel.h"
#include "ThreadPool.h"
#include "EquationOfState.h"
#include <vector>
#include <memory>
#include <cstdint>

struct Parameters {
    double gamma = 1.4;
    double alpha = 1.0;
    double beta = 1.0;
    double eta = 1.2;
    double thermal = 0.5;
    double boltzmannConstant = 1.0;
    double thermalViscosity = 1.0;
    //double fluidDensity = 1.0;

    bool enableGravity = false;
    Vector2D gravity = Vector2D(0.0, -9.81);
};

class SPHSolver {
public:
    SPHSolver(Parameters params, std::shared_ptr<EquationOfState> eos);
    void addThreadPool(std::shared_ptr<ThreadPool> pool);

    void updateSmoothingLengths(std::vector<std::shared_ptr<Particle>>& fluidParticles) const;

    void computeDensityAndPressure(
        std::vector<std::shared_ptr<Particle>>& fluidParticles,
        const std::vector<std::shared_ptr<Particle>>& internalRigidParticles,
        const std::vector<std::shared_ptr<Particle>>& ghostParticles,
        const SpatialGrid& grid,
        const Kernel& kernel) const;

    void setThermalNoiseEnabled(bool enabled) { enableThermalNoise_ = enabled; }

    void computeDerivatives(
        std::vector<std::shared_ptr<Particle>>& fluidParticles,
        const std::vector<std::shared_ptr<Particle>>& internalRigidParticles,
        const std::vector<std::shared_ptr<Particle>>& ghostParticles,
        const SpatialGrid& grid,
        const Kernel& kernel,
        double dt,
        std::uint64_t noiseSample) const;

    void setParameters(const Parameters& params) { params_ = params; }
    const Parameters& getParameters() const { return params_; }

    void setEOS(std::shared_ptr<EquationOfState> eos) { eos_ = std::move(eos); }
    const EquationOfState& getEOS() const { return *eos_; }

private:
    std::shared_ptr<ThreadPool> threadPool_ = nullptr;
    Parameters params_;
    std::shared_ptr<EquationOfState> eos_;
    bool enableThermalNoise_ = false;
    //double computeSoundSpeed(double u, double gamma) const;
};