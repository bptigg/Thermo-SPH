#include "simulation.h"
#include "GhostManager.h"
#include "Piston.h"
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
    solver_.setThermalNoiseEnabled(config_.enableThermalNoise);
    
    // Filter out fluid particles from initial conditions
    auto allParticles = initialConditions->generateParticles();
    for (const auto& p : allParticles) {
        if (p->isFluid()) {
            fluidParticles_.push_back(p);
        }
        else if (p->hasBoundary()) {
            boundaryParticles_.push_back(p);
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

    // 3. Initial SPH evaluation using 5-argument API
    solver_.updateSmoothingLengths(fluidParticles_);
    grid_.build(fluidParticles_);
    ghosts_ = GhostManager::generateAllGhosts(boundaryParticles_, fluidParticles_, config_.supportRadius, grid_, threadPool_);

    // 4. Gather internal rigid body constituent particles
    std::vector<std::shared_ptr<Particle>> internalRigidParticles;
    collectInternalRigidParticles(internalRigidParticles, grid_);
    config_.bath.stateProvider = initialConditions->getThermalStateProvider();
    
    solver_.computeDensityAndPressure(fluidParticles_, internalRigidParticles, ghosts_, grid_, *kernel_);
    solver_.computeDerivatives(fluidParticles_, internalRigidParticles, ghosts_, grid_, *kernel_, dt_, 0);
}

void SimulationEngine::classifyRigidBodies() {
    double tolerance = 1e-3;

    for (auto& obj : rigidObjects_) {
        obj->finalizeInitialization();

        //if (obj->isStatic()) {
        //    AABB box = obj->getAABB();
//
        //    bool touchesMinX = std::abs(box.min.x - config_.domainMin.x) < tolerance;
        //    bool touchesMaxX = std::abs(box.max.x - config_.domainMax.x) < tolerance;
        //    bool touchesMinY = std::abs(box.min.y - config_.domainMin.y) < tolerance;
        //    bool touchesMaxY = std::abs(box.max.y - config_.domainMax.y) < tolerance;
//
        //    if (touchesMinX || touchesMaxX || touchesMinY || touchesMaxY) {
        //        obj->setType(RigidBodyType::EXTERNAL_BOUNDARY);
        //    } else {
        //        obj->setType(RigidBodyType::INTERNAL_OBJECT);
        //    }
        //} else {
        //    obj->setType(RigidBodyType::INTERNAL_OBJECT);
        //}
    }
}

void SimulationEngine::collectInternalRigidParticles(std::vector<std::shared_ptr<Particle>>& outParticles, const SpatialGrid& grid) {
    outParticles.clear();
    double cellSize = grid.getCellSize();
    for (auto& obj : rigidObjects_) {
        obj->clearForces();
        if (obj->getType() == RigidBodyType::INTERNAL_OBJECT) {
            for (const auto& p : obj->getParticles()) {
                if (!p) continue;

                int cellX = static_cast<int>(std::floor(p->pos.x / cellSize));
                int cellY = static_cast<int>(std::floor(p->pos.y / cellSize));
                bool nearFluid = false;
                for (int dx = -1; dx <= 1; ++dx) {
                    for (int dy = -1; dy <= 1; ++dy) {
                        if (grid.isCellActive(cellX + dx, cellY + dy)) {
                            nearFluid = true;
                            break;
                        }
                    }
                    if (nearFluid) break;
                }
                if (nearFluid) {
                    outParticles.push_back(p);
                }
            }
        }
    }
}

void SimulationEngine::step() {
    // 1. First half kick & drift for fluid particles
    integrator_.kickFirstHalf(fluidParticles_, dt_, threadPool_);
    integrator_.drift(fluidParticles_, dt_, threadPool_);
    threadPool_.waitFinished();
    for (auto& obj : rigidObjects_) {
        obj->constrainFluid(fluidParticles_);
    }

    // 2. Re-evaluate h_i & Grid
    solver_.updateSmoothingLengths(fluidParticles_);
    grid_.build(fluidParticles_);

    // 3. Generate Ghost Particles for External Boundaries
    ghosts_ = GhostManager::generateAllGhosts(boundaryParticles_, fluidParticles_, config_.supportRadius, grid_, threadPool_);

    // 4. Gather internal rigid body particles & clear object accumulators
    std::vector<std::shared_ptr<Particle>> internalRigidParticles;
    collectInternalRigidParticles(internalRigidParticles, grid_);

    // 5. SPH Thermodynamics & Fluid Forces (5 arguments matching SPHSolver)
    solver_.computeDensityAndPressure(fluidParticles_, internalRigidParticles, ghosts_, grid_, *kernel_);
    solver_.computeDerivatives(fluidParticles_, internalRigidParticles, ghosts_, grid_, *kernel_, dt_, currentStep_ + 1);
    ThermalBath bath = config_.bath;
    //if (bath.hotU > 0.0) {
    //    for (const auto& obj : rigidObjects_) {
    //        if (auto piston = std::dynamic_pointer_cast<Piston>(obj);
    //            piston && piston->getCenterOfMass().y < bath.hotBelowY) {
    //            bath.targetU = bath.hotU;
    //        }
    //    }
    //}
    
    bathHeat_ += bath.apply(fluidParticles_, config_.domainMin, config_.domainMax) * dt_;

    // 6. Solid-Solid Collisions & Newton-Euler Integration
    for (auto& obj : rigidObjects_) {
        obj->preStep(currentTime_, dt_);
    }
    //rigidSolver_.solveCollisions(rigidObjects_, dt_, threadPool_);
    rigidSolver_.integrate(rigidObjects_, dt_);

    // 7. Second half kick & Adaptive Timestep for Fluid Particles
    integrator_.kickSecondHalf(fluidParticles_, dt_, threadPool_);

    currentTime_ += dt_;
    dt_ = integrator_.computeAdaptiveTimestep(fluidParticles_, dt_);
    currentStep_++;
}

void SimulationEngine::run() {
    threadPool_.start();
    while (currentTime_ < config_.tFinal) {
        step();
    }
    threadPool_.Stop();
}