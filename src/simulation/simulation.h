#pragma once

#include <vector>
#include <memory>
#include "Particle.h"
#include "RigidBody.h"
#include "sphsolver.h"
#include "kernel.h"
#include "SpatialGrid.h"
#include "VelocityVerlet.h"
#include "ThreadPool.h"
#include "initialconditions.h"

class SimulationEngine {
public:
    struct Configuration {
        double tFinal = 10.0;
        double initialDt = 1e-4;
        double supportRadius = 0.2;
    };

    SimulationEngine(
        Configuration config,
        SPHSolver solver,
        std::shared_ptr<Kernel> kernel,
        std::unique_ptr<InitialConditions> initialConditions,
        size_t threadCount = 4);

    ~SimulationEngine() = default;

    void run();
    void step();

    double getCurrentTime() const { return currentTime_; }
    size_t getCurrentStep() const { return currentStep_; }
    
    const std::vector<std::shared_ptr<Particle>>& getParticles() const { return particles_; }
    const std::vector<std::shared_ptr<RigidObject>>& getRigidObjects() const { return rigidObjects_; }

    ThreadPool& getThreadPool() {return threadPool_;} 

private:
    Configuration config_;
    SPHSolver solver_;
    std::shared_ptr<Kernel> kernel_;
    SpatialGrid grid_;
    VelocityVerletIntegrator integrator_;
    ThreadPool threadPool_;

    std::vector<std::shared_ptr<Particle>> particles_;
    std::vector<std::shared_ptr<RigidObject>> rigidObjects_;
    std::vector<std::shared_ptr<Particle>> ghosts_;

    double currentTime_ = 0.0;
    double dt_ = 1e-4;
    size_t currentStep_ = 0;
};