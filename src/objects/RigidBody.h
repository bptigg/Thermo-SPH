#pragma once
#include "Particle.h"
#include <vector>
#include <memory>

#include <cmath>

struct AABB {
    Vector2D min;
    Vector2D max;

    bool overlaps(const AABB& other) const;
};

class RigidObject {
protected:
    std::vector<std::shared_ptr<Particle>> particles_;
    std::vector<Vector2D> localOffsets_;

    double totalMass_ = 0.0;
    double inertia_ = 0.0;

    Vector2D centerOfMass_{0.0, 0.0};
    Vector2D linearVel_{0.0, 0.0};
    double angle_ = 0.0;
    double angularVel_ = 0.0;

    Vector2D forceAccumulator_{0.0, 0.0};
    double torqueAccumulator_ = 0.0;

    bool lockRotation_ = false;
    bool lockX_ = false;
    bool lockY_ = false;

public:
    virtual ~RigidObject() = default;

    void addParticle(std::shared_ptr<Particle> p);
    const std::vector<std::shared_ptr<Particle>>& getParticles() const { return particles_; }

    void finalizeInitialization();
    void accumulateForces();
    virtual void updateKinematics(double dt);

    void applyImpulse(const Vector2D& impulse, const Vector2D& r);

    AABB getAABB() const;
    void setConstraints(bool lockX, bool lockY, bool lockRotation);

    // Getters
    const Vector2D& getCenterOfMass() const { return centerOfMass_; }
    const Vector2D& getLinearVel() const { return linearVel_; }
    double getAngularVel() const { return angularVel_; }
    double getTotalMass() const { return totalMass_; }
    double getInertia() const { return inertia_; }
    bool isRotationLocked() const { return lockRotation_; }
    virtual bool isStatic() const { return lockX_ && lockY_ && lockRotation_; }
};