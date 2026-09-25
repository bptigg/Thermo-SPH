#include "simulation.h"
#include "GhostManager.h"
#include <cmath>
#include <algorithm>

SimulationEngine::SimulationEngine(
    Configuration config,
    SPHSolver solver,
    std::shared_ptr<Kernel> kernel,
    std::unique_ptr<InitialConditions> initialConditions,
    size_t threadCount)
    : config_(config),
      solver_(solver),
      kernel_(std::move(kernel)),
      integrator_(grid_),
      threadPool_(threadCount),
      dt_(config.initialDt)
{
    threadPool_.start();
    
    // Filter out fluid particles from initial conditions
    auto allParticles = initialConditions->generateParticles();
    for (const auto& p : allParticles) {
        if (p->isFluid()) {
            fluidParticles_.push_back(p);
        }
    }

    rigidObjects_ = initialConditions->getRigidObjects();

    // 1. Configure Rigid Body Solver parameters (including optional gravity)
    RigidBodySolverParameters rbParams;
    rbParams.enableGravity = config_.enableGravity;
    rbParams.gravity = config_.gravity;
    rigidSolver_.setParameters(rbParams);

    // 2. Classify rigid objects as EXTERNAL_BOUNDARY vs INTERNAL_OBJECT
    classifyRigidBodies();

    // 3. Gather internal rigid body constituent particles
    std::vector<std::shared_ptr<Particle>> internalRigidParticles;
    collectInternalRigidParticles(internalRigidParticles);

    // 4. Initial SPH evaluation using 5-argument API
    solver_.updateSmoothingLengths(fluidParticles_);
    grid_.build(fluidParticles_);
    ghosts_ = GhostManager::generateAllGhosts(fluidParticles_, config_.supportRadius);
    
    solver_.computeDensityAndPressure(fluidParticles_, internalRigidParticles, ghosts_, grid_, *kernel_);
}

void SimulationEngine::classifyRigidBodies() {
    double tolerance = 1e-3;

    for (auto& obj : rigidObjects_) {
        obj->finalizeInitialization();

        if (obj->isStatic()) {
            AABB box = obj->getAABB();

            bool touchesMinX = std::abs(box.min.x - config_.domainMin.x) < tolerance;
            bool touchesMaxX = std::abs(box.max.x - config_.domainMax.x) < tolerance;
            bool touchesMinY = std::abs(box.min.y - config_.domainMin.y) < tolerance;
            bool touchesMaxY = std::abs(box.max.y - config_.domainMax.y) < tolerance;

            if (touchesMinX || touchesMaxX || touchesMinY || touchesMaxY) {
                obj->setType(RigidBodyType::EXTERNAL_BOUNDARY);
            } else {
                obj->setType(RigidBodyType::INTERNAL_OBJECT);
            }
        } else {
            obj->setType(RigidBodyType::INTERNAL_OBJECT);
        }
    }
}

void SimulationEngine::collectInternalRigidParticles(std::vector<std::shared_ptr<Particle>>& outParticles) {
    outParticles.clear();
    for (auto& obj : rigidObjects_) {
        obj->clearForces();
        if (obj->getType() == RigidBodyType::INTERNAL_OBJECT) {
            for (const auto& p : obj->getParticles()) {
                outParticles.push_back(p);
            }
        }
    }
}

void SimulationEngine::step() {
    // 1. First half kick & drift for fluid particles
    integrator_.kickFirstHalf(fluidParticles_, dt_, threadPool_);
    integrator_.drift(fluidParticles_, dt_, threadPool_);

    // 2. Re-evaluate h_i & Grid
    solver_.updateSmoothingLengths(fluidParticles_);
    grid_.build(fluidParticles_);

    // 3. Generate Ghost Particles for External Boundaries
    ghosts_ = GhostManager::generateAllGhosts(fluidParticles_, config_.supportRadius);

    // 4. Gather internal rigid body particles & clear object accumulators
    std::vector<std::shared_ptr<Particle>> internalRigidParticles;
    collectInternalRigidParticles(internalRigidParticles);

    // 5. SPH Thermodynamics & Fluid Forces (5 arguments matching SPHSolver)
    solver_.computeDensityAndPressure(fluidParticles_, internalRigidParticles, ghosts_, grid_, *kernel_);
    solver_.computeDerivatives(fluidParticles_, internalRigidParticles, ghosts_, grid_, *kernel_);

    // 6. Solid-Solid Collisions & Newton-Euler Integration
    rigidSolver_.solveCollisions(rigidObjects_, dt_);
    rigidSolver_.integrate(rigidObjects_, dt_);

    // 7. Second half kick & Adaptive Timestep for Fluid Particles
    integrator_.kickSecondHalf(fluidParticles_, dt_, threadPool_);

    dt_ = integrator_.computeAdaptiveTimestep(fluidParticles_, dt_);
    currentTime_ += dt_;
    currentStep_++;
}

void SimulationEngine::run() {
    threadPool_.start();
    while (currentTime_ < config_.tFinal) {
        step();
    }
    threadPool_.Stop();
}