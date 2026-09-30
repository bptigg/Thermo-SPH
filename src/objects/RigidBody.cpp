#include "RigidBody.h"
#include <cmath>
#include <algorithm>
#include <atomic>

bool AABB::overlaps(const AABB& other) const {
    return !(max.x < other.min.x || min.x > other.max.x ||
             max.y < other.min.y || min.y > other.max.y);
}

void RigidObject::addParticle(std::shared_ptr<Particle> p) {
    particles_.push_back(p);
}

void RigidObject::finalizeInitialization() {
    totalMass_ = 0.0;
    centerOfMass_ = Vector2D(0.0, 0.0);

    for (const auto& p : particles_) {
        p->parentBody = this; // Link constituent particle back to this object
        totalMass_ += p->mass;
        centerOfMass_ += p->pos * p->mass;
    }

    if (totalMass_ > 0.0) {
        centerOfMass_ /= totalMass_;
    }

    inertia_ = 0.0;
    localOffsets_.clear();
    for (const auto& p : particles_) {
        Vector2D r = p->pos - centerOfMass_;
        localOffsets_.push_back(r);
        inertia_ += p->mass * (r.x * r.x + r.y * r.y);
    }
    if (inertia_ == 0.0) inertia_ = 1.0;
}

void RigidObject::clearForces() {
    forceAccumulator_ = Vector2D(0.0, 0.0);
    torqueAccumulator_ = 0.0;
}

void RigidObject::addForceAtPosition(const Vector2D& force, const Vector2D& worldPos) {
    if (isStatic()) return;

    Vector2D f = force;
    if (lockX_) f.x = 0.0;
    if (lockY_) f.y = 0.0;

    std::atomic_ref<double> forceX(forceAccumulator_.x);
    std::atomic_ref<double> forceY(forceAccumulator_.y);
    forceX.fetch_add(f.x, std::memory_order_relaxed);
    forceY.fetch_add(f.y, std::memory_order_relaxed);

    //forceAccumulator_ += f;

    if (!lockRotation_) {
        Vector2D r = worldPos - centerOfMass_;
        double torque = (r.x * f.y - r.y * f.x);
        
        std::atomic_ref<double> torqueRef(torqueAccumulator_);
        torqueRef.fetch_add(torque, std::memory_order_relaxed);
    }
}

void RigidObject::updateParticlePositions() {
    double cosA = std::cos(angle_);
    double sinA = std::sin(angle_);

    if (pivotEnabled_) {
        Vector2D rotatedPivotOffset(
            cosA * pivotLocalOffset_.x - sinA * pivotLocalOffset_.y,
            sinA * pivotLocalOffset_.x + cosA * pivotLocalOffset_.y
        );
        centerOfMass_ = pivotWorld_ - rotatedPivotOffset;
        linearVel_ = Vector2D(
            -angularVel_ * (centerOfMass_.y - pivotWorld_.y),
            angularVel_ * (centerOfMass_.x - pivotWorld_.x)
        );
    }

    for (size_t i = 0; i < particles_.size(); ++i) {
        Vector2D rLocal = localOffsets_[i];
        Vector2D rRotated(
            cosA * rLocal.x - sinA * rLocal.y,
            sinA * rLocal.x + cosA * rLocal.y
        );

        particles_[i]->pos = centerOfMass_ + rRotated;
        Vector2D rotVel(-angularVel_ * rRotated.y, angularVel_ * rRotated.x);
        particles_[i]->vel = linearVel_ + rotVel;
    }
}

void RigidObject::applyImpulse(const Vector2D& impulse, const Vector2D& r) {
    if (pivotEnabled_) {
        double inertia = getInertiaAboutPivot();
        if (!lockRotation_ && inertia > 0.0) {
            Vector2D impulseArm = r + centerOfMass_ - pivotWorld_;
            double angularImpulse = impulseArm.x * impulse.y - impulseArm.y * impulse.x;
            angularVel_ += angularImpulse / inertia;
            Vector2D pivotToCenter = centerOfMass_ - pivotWorld_;
            linearVel_ = Vector2D(-angularVel_ * pivotToCenter.y, angularVel_ * pivotToCenter.x);
        }
        return;
    }

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

void RigidObject::setPivot(const Vector2D& worldPosition) {
    pivotEnabled_ = true;
    pivotWorld_ = worldPosition;

    double cosA = std::cos(angle_);
    double sinA = std::sin(angle_);
    Vector2D pivotOffset = pivotWorld_ - centerOfMass_;
    pivotLocalOffset_ = Vector2D(
        cosA * pivotOffset.x + sinA * pivotOffset.y,
        -sinA * pivotOffset.x + cosA * pivotOffset.y
    );

    Vector2D pivotToCenter = centerOfMass_ - pivotWorld_;
    linearVel_ = Vector2D(-angularVel_ * pivotToCenter.y, angularVel_ * pivotToCenter.x);
}

double RigidObject::getInertiaAboutPivot() const {
    return inertia_ + totalMass_ * (centerOfMass_ - pivotWorld_).normSq();
}

void RigidObject::setAngularSpring(double stiffness, double damping, double restAngle) {
    angularSpringStiffness_ = std::max(0.0, stiffness);
    angularSpringDamping_ = std::max(0.0, damping);
    angularSpringRestAngle_ = restAngle;
}

double RigidObject::getAngularSpringTorque() const {
    double angleError = std::remainder(angle_ - angularSpringRestAngle_, 2.0 * std::acos(-1.0));
    return -angularSpringStiffness_ * angleError - angularSpringDamping_ * angularVel_;
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