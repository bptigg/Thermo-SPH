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

enum class RigidBodyType {
    EXTERNAL_BOUNDARY, // Outer walls/containers (generates ghost particles)
    INTERNAL_OBJECT    // Pistons, debris, obstacles (interacts via constituent particles)
};

class RigidObject {
protected:
    std::vector<std::shared_ptr<Particle>> particles_;
    std::vector<Vector2D> localOffsets_;

    RigidBodyType type_ = RigidBodyType::INTERNAL_OBJECT;

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
    void clearForces();
    
    // Accumulates force at a specific constituent particle position
    void addForceAtPosition(const Vector2D& force, const Vector2D& worldPos);
    
    // Synchronizes all constituent particles to the current center of mass and rotation angle
    void updateParticlePositions();

    void applyImpulse(const Vector2D& impulse, const Vector2D& r);

    AABB getAABB() const;
    void setConstraints(bool lockX, bool lockY, bool lockRotation);

    // Getters
    void setType(RigidBodyType t) { type_ = t; }
    RigidBodyType getType() const { return type_; }

    const Vector2D& getCenterOfMass() const { return centerOfMass_; }
    void setCenterOfMass(const Vector2D& pos) { centerOfMass_ = pos; }

    const Vector2D& getLinearVel() const { return linearVel_; }
    void setLinearVel(const Vector2D& vel) { linearVel_ = vel; }

    double getAngle() const { return angle_; }
    void setAngle(double a) { angle_ = a; }

    double getAngularVel() const { return angularVel_; }
    void setAngularVel(double omega) { angularVel_ = omega; }

    double getTotalMass() const { return totalMass_; }
    double getInertia() const { return inertia_; }

    Vector2D getForceAccumulator() const { return forceAccumulator_; }
    double getTorqueAccumulator() const { return torqueAccumulator_; }

    bool isRotationLocked() const { return lockRotation_; }
    bool isXLocked() const { return lockX_; }
    bool isYLocked() const { return lockY_; }
    virtual bool isStatic() const { return lockX_ && lockY_ && lockRotation_; }
};