#include "Piston.h"
#include <atomic>
#include <cmath>

Piston::Piston(DirectionLocks locks) : locks_(locks) {
    setType(RigidBodyType::INTERNAL_OBJECT);
    setConstraints(false, false, true);
}

std::shared_ptr<Piston> Piston::makeKinematic(const Vector2D& velocity) {
    auto piston = std::make_shared<Piston>(DirectionLocks::all());
    piston->setLinearVel(velocity);
    return piston;
}

void Piston::scheduleRelease(double time, const DirectionLocks& locksAfter) {
    releaseTime_ = time;
    locksAfterRelease_ = locksAfter;
    released_ = false;
}

void Piston::enableFlywheel(double crankRadius, double inertia, double initialOmega, double shaftDamping) {
    flywheel_ = true;
    crankRadius_ = crankRadius;
    inertia_ = inertia;
    omega_ = initialOmega;
    shaftDamping_ = shaftDamping;
    theta_ = 0.0;
    crankCentreY_ = centerOfMass_.y - crankRadius;
    locks_ = DirectionLocks::all();
}

void Piston::addForceAtPosition(const Vector2D& force, const Vector2D& worldPos) {
    if (flywheel_) {
        std::atomic_ref<double> gasForce(gasForceY_);
        gasForce.fetch_add(force.y, std::memory_order_relaxed);
    }
    Vector2D f = force;
    if ((f.x > 0.0 && locks_.posX) || (f.x < 0.0 && locks_.negX)) f.x = 0.0;
    if ((f.y > 0.0 && locks_.posY) || (f.y < 0.0 && locks_.negY)) f.y = 0.0;
    RigidObject::addForceAtPosition(f, worldPos);
}

void Piston::constrainFluid(const std::vector<std::shared_ptr<Particle>>& fluid) {
    if (!sealLid_ || particles_.empty()) return;

    const AABB box = getAABB();
    const double surface = box.min.y - lidClearance_;
    const double band = (box.max.y - box.min.y) + lidClearance_;

    for (const auto& p : fluid) {
        if (p->pos.x < box.min.x || p->pos.x > box.max.x) continue;
        if (p->pos.y > surface && p->pos.y < surface + band) {
            p->pos.y = surface;
            if (p->vel.y > linearVel_.y) p->vel.y = linearVel_.y;
        }
    }
}

void Piston::preStep(double time, double dt) {
    if (flywheel_) {
        // Virtual work: the force on the piston times dy/dtheta is the torque on the crank
        const double dyDtheta = -crankRadius_ * std::sin(theta_);
        const double force = gasForceY_ + load_.y;
        gasForceY_ = 0.0;
        const double torque = force * dyDtheta - shaftDamping_ * omega_;
        harvested_ += shaftDamping_ * omega_ * omega_ * dt;
        omega_ += torque / inertia_ * dt;
        theta_ += omega_ * dt;
        const double targetY = crankCentreY_ + crankRadius_ * std::cos(theta_);
        linearVel_ = Vector2D(0.0, (targetY - centerOfMass_.y) / dt);
        return;
    }

    if (!released_ && releaseTime_ >= 0.0 && time >= releaseTime_) {
        locks_ = locksAfterRelease_;
        setLinearVel(Vector2D(0.0, 0.0));
        released_ = true;
    }

    if (!isKinematic() && load_.normSq() > 0.0) {
        RigidObject::addForceAtPosition(load_, centerOfMass_);
    }

    if (!isKinematic() && ptoDamping_ > 0.0 && linearVel_.y > 0.0) {
        const double vy = linearVel_.y;
        RigidObject::addForceAtPosition(Vector2D(0.0, -ptoDamping_ * vy), centerOfMass_);
        harvested_ += ptoDamping_ * vy * vy * dt;
    }

    bool atLimit = false;
    if (centerOfMass_.y <= minY_) {
        centerOfMass_.y = minY_;
        if (linearVel_.y < 0.0) linearVel_.y = 0.0;
        atLimit = true;
    } else if (centerOfMass_.y >= maxY_) {
        centerOfMass_.y = maxY_;
        if (linearVel_.y > 0.0) linearVel_.y = 0.0;
        atLimit = true;
    }
    if (atLimit) updateParticlePositions();
}
