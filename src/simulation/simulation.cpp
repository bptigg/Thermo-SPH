#include "simulation.h"
#include "GhostManager.h"

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
    particles_ = initialConditions->generateParticles();
    rigidObjects_ = initialConditions->getRigidObjects();

    for (auto& obj : rigidObjects_) {
        obj->finalizeInitialization();
    }

    solver_.updateSmoothingLengths(particles_);
    grid_.build(particles_);
    ghosts_ = GhostManager::generateAllGhosts(particles_, config_.supportRadius);
    
    // Dereference pointer (*kernel_) to pass const Kernel& to SPHSolver methods
    solver_.computeDensityAndPressure(particles_, ghosts_, grid_, *kernel_);
}

void SimulationEngine::step() {
    Vector2D netBoundaryForce(0.0, 0.0);

    // 1. Kick & Drift
    integrator_.kickFirstHalf(particles_, dt_, threadPool_);
    integrator_.drift(particles_, dt_, threadPool_);

    // 2. Re-evaluate h_i & Grid
    solver_.updateSmoothingLengths(particles_);
    grid_.build(particles_);
    ghosts_ = GhostManager::generateAllGhosts(particles_, config_.supportRadius);

    // 3. SPH Field Derivatives (Pass *kernel_)
    solver_.computeDensityAndPressure(particles_, ghosts_, grid_, *kernel_);
    solver_.computeDerivatives(particles_, ghosts_, grid_, *kernel_, netBoundaryForce);

    // 4. Rigid Object Kinematics
    for (auto& obj : rigidObjects_) {
        obj->accumulateForces();
        obj->updateKinematics(dt_);
    }

    // 5. Final Kick
    integrator_.kickSecondHalf(particles_, dt_, threadPool_);

    dt_ = integrator_.computeAdaptiveTimestep(particles_, dt_);
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