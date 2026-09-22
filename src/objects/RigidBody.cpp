#include "RigidBody.h"
#include <cmath>
#include <algorithm>

bool AABB::overlaps(const AABB& other) const {
    return (min.x <= other.max.x && max.x >= other.min.x &&
            min.y <= other.max.y && max.y >= other.min.y);
}

void RigidObject::addParticle(std::shared_ptr<Particle> p) {
    particles_.push_back(p);
}

void RigidObject::finalizeInitialization() {
    totalMass_ = 0.0;
    centerOfMass_ = Vector2D(0.0, 0.0);

    for (const auto& p : particles_) {
        totalMass_ += p->mass;
        centerOfMass_ += p->pos * p->mass;
    }

    if (totalMass_ > 0.0) {
        centerOfMass_ /= totalMass_;
    }

    // Calculate initial local offsets and Moment of Inertia (I)
    inertia_ = 0.0;
    localOffsets_.clear();
    for (const auto& p : particles_) {
        Vector2D r = p->pos - centerOfMass_;
        localOffsets_.push_back(r);
        inertia_ += p->mass * (r.x * r.x + r.y * r.y);
    }
}

void RigidObject::accumulateForces() {
    forceAccumulator_ = Vector2D(0.0, 0.0);
    torqueAccumulator_ = 0.0;

    for (size_t i = 0; i < particles_.size(); ++i) {
        const auto& p = particles_[i];
        Vector2D force = p->accel * p->mass; // Hydrodynamic force from SPH
        Vector2D r = p->pos - centerOfMass_;

        forceAccumulator_ += force;
        // 2D Cross product (r x F) = r.x * F.y - r.y * F.x
        torqueAccumulator_ += (r.x * force.y - r.y * force.x);
    }
}

void RigidObject::updateKinematics(double dt) {
    if (totalMass_ <= 0.0) return;

    // 1. Angular Integration
    if (!lockRotation_ && inertia_ > 0.0) {
        double angularAccel = torqueAccumulator_ / inertia_;
        angularVel_ += angularAccel * dt;
        angle_ += angularVel_ * dt;
    } else {
        angularVel_ = 0.0;
    }

    // 2. Linear Integration
    Vector2D linearAccel = forceAccumulator_ / totalMass_;
    if (lockX_) linearAccel.x = 0.0;
    if (lockY_) linearAccel.y = 0.0;

    linearVel_ += linearAccel * dt;
    if (lockX_) linearVel_.x = 0.0;
    if (lockY_) linearVel_.y = 0.0;

    centerOfMass_ += linearVel_ * dt;

    // 3. Synchronize All Constituent Particles
    double cosA = std::cos(angle_);
    double sinA = std::sin(angle_);

    for (size_t i = 0; i < particles_.size(); ++i) {
        // Rotate local offset vector
        Vector2D rLocal = localOffsets_[i];
        Vector2D rRotated(
            cosA * rLocal.x - sinA * rLocal.y,
            sinA * rLocal.x + cosA * rLocal.y
        );

        // Update particle position
        particles_[i]->pos = centerOfMass_ + rRotated;

        // Update particle velocity: V_cm + (Omega x r)
        Vector2D rotVel(-angularVel_ * rRotated.y, angularVel_ * rRotated.x);
        particles_[i]->vel = linearVel_ + rotVel;
    }
}

void RigidObject::applyImpulse(const Vector2D& impulse, const Vector2D& r) {
    if (!lockX_ || !lockY_) {
        Vector2D linImpulse = impulse;
        if (lockX_) linImpulse.x = 0.0;
        if (lockY_) linImpulse.y = 0.0;
        
        if (totalMass_ > 0.0) {
            linearVel_ += linImpulse / totalMass_;
        }
    }

    if (!lockRotation_ && inertia_ > 0.0) {
        double angularImpulse = r.x * impulse.y - r.y * impulse.x;
        angularVel_ += angularImpulse / inertia_;
    }
}

AABB RigidObject::getAABB() const {
    AABB box;
    if (particles_.empty()) return box;

    box.min = particles_[0]->pos;
    box.max = particles_[0]->pos;

    for (const auto& p : particles_) {
        box.min.x = std::min(box.min.x, p->pos.x);
        box.min.y = std::min(box.min.y, p->pos.y);
        box.max.x = std::max(box.max.x, p->pos.x);
        box.max.y = std::max(box.max.y, p->pos.y);
    }
    return box;
}

void RigidObject::setConstraints(bool lockX, bool lockY, bool lockRotation) {
    lockX_ = lockX;
    lockY_ = lockY;
    lockRotation_ = lockRotation;
}