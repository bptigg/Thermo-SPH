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

    // Solves collisions between pairs of rigid objects using AABB and particle contacts
    void solveCollisions(std::vector<std::shared_ptr<RigidObject>>& bodies, double dt);

    // Newton-Euler integration of rigid objects
    void integrate(std::vector<std::shared_ptr<RigidObject>>& bodies, double dt);

    void setParameters(const RigidBodySolverParameters& params) { params_ = params; }
    const RigidBodySolverParameters& getParameters() const { return params_; }

private:
    RigidBodySolverParameters params_;
};