#pragma once

#include <vector>
#include <memory>
#include "Particle.h"
#include "RigidBody.h"
#include "rigidbodysolver.h"
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
        double supportRadius = 0.1;

        // Bounding box for domain classification
        Vector2D domainMin = Vector2D(0.0, 0.0);
        Vector2D domainMax = Vector2D(1.0, 1.0);

        // Optional Gravity Parameters
        bool enableGravity = false;
        Vector2D gravity = Vector2D(0.0, -9.81);
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
    
    const std::vector<std::shared_ptr<Particle>>& getParticles() const { return fluidParticles_; }
    const std::vector<std::shared_ptr<Particle>>& getFluidParticles() const { return fluidParticles_; }
    const std::vector<std::shared_ptr<RigidObject>>& getRigidObjects() const { return rigidObjects_; }

    ThreadPool& getThreadPool() { return threadPool_; } 

private:
    void classifyRigidBodies();
    void collectInternalRigidParticles(std::vector<std::shared_ptr<Particle>>& outParticles);

    Configuration config_;
    SPHSolver solver_;
    RigidBodySolver rigidSolver_;
    std::shared_ptr<Kernel> kernel_;
    SpatialGrid grid_;
    VelocityVerletIntegrator integrator_;
    ThreadPool threadPool_;

    std::vector<std::shared_ptr<Particle>> fluidParticles_;
    std::vector<std::shared_ptr<RigidObject>> rigidObjects_;
    std::vector<std::shared_ptr<Particle>> ghosts_;

    double currentTime_ = 0.0;
    double dt_ = 1e-4;
    size_t currentStep_ = 0;
};