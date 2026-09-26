#pragma once
#include "RigidBody.h"
#include <vector>
#include <memory>

struct RigidBodySolverParameters {
    bool enableGravity = false;
    Vector2D gravity = Vector2D(0.0, -9.81);
    double restitution = 0.1;
    double particleRadius = 0.05;
};

class RigidBodySolver {
public:
    explicit RigidBodySolver(RigidBodySolverParameters params = {});

    // Phase 1: Integrate external forces & gravity into linear and angular velocities
    void integrateVelocities(std::vector<std::shared_ptr<RigidObject>>& bodies, double dt);

    // Phase 2: Detect particle contact overlaps and resolve via direct impulse & position unwrapping
    bool solveCollisions(std::vector<std::shared_ptr<RigidObject>>& bodies, double dt);

    // Phase 3: Advance positions and rotations using post-collision velocities
    void integratePositions(std::vector<std::shared_ptr<RigidObject>>& bodies, double dt);

    // Full pipeline convenience step: integrateVelocities -> solveCollisions -> integratePositions
    bool integrate(std::vector<std::shared_ptr<RigidObject>>& bodies, double dt);

    void setParameters(const RigidBodySolverParameters& params) { params_ = params; }
    const RigidBodySolverParameters& getParameters() const { return params_; }

private:
    RigidBodySolverParameters params_;
};