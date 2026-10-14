#pragma once
#include "RigidBody.h"
#include <limits>

// A locked direction discards force components pointing that way (e.g. negY: pushes downward are ignored).
struct DirectionLocks {
    bool posX = false;
    bool negX = false;
    bool posY = false;
    bool negY = false;

    static DirectionLocks all() { return {true, true, true, true}; }
    // Free to be pushed up only: the usual gas-driven piston
    static DirectionLocks upwardOnly() { return {true, true, false, true}; }

    bool isAll() const { return posX && negX && posY && negY; }
};

// Rigid body that exchanges force with the fluid through its constituent particles.
// A kinematic piston is a Piston whose locks are all on: gas forces are ignored and it keeps its set velocity.
class Piston : public RigidObject {
public:
    explicit Piston(DirectionLocks locks = DirectionLocks::upwardOnly());

    static std::shared_ptr<Piston> makeKinematic(const Vector2D& velocity);

    void setLocks(const DirectionLocks& locks) { locks_ = locks; }
    const DirectionLocks& getLocks() const { return locks_; }
    bool isKinematic() const { return locks_.isAll(); }

    // Constant force (e.g. weight / external load) applied to dynamic phases, ignoring the locks
    void setConstantLoad(const Vector2D& load) { load_ = load; }
    const Vector2D& getConstantLoad() const { return load_; }

    // Energy siphon: a drag force -c*v on the piston while it rises (the gas's power stroke), like a generator.
    // Inactive while kinematic.
    void setPowerTakeOff(double damping) { ptoDamping_ = damping; }
    double getPowerTakeOff() const { return ptoDamping_; }
    // Cumulative work removed by the power take-off
    double getHarvestedEnergy() const { return harvested_; }

    // Limits on the centre-of-mass y position; the piston stops when it reaches either
    void setTravelLimits(double minY, double maxY) { minY_ = minY; maxY_ = maxY; }

    // Slider-crank flywheel: the piston height follows y = yMid + r*cos(theta) (theta = 0 at top dead centre).
    // Gas pressure torques the flywheel through the crank, and the flywheel drives the piston back down.
    // `shaftDamping` is a generator on the shaft that removes c*omega^2 of power. The piston's own mass is neglected.
    void enableFlywheel(double crankRadius, double inertia, double initialOmega, double shaftDamping = 0.0);
    bool hasFlywheel() const { return flywheel_; }
    double getCrankAngle() const { return theta_; }
    double getFlywheelOmega() const { return omega_; }
    double getFlywheelEnergy() const { return 0.5 * inertia_ * omega_ * omega_; }

    // At simulation time t the piston switches to `locksAfter` and its velocity is zeroed
    void scheduleRelease(double time, const DirectionLocks& locksAfter);
    bool isReleased() const { return released_; }

    void addForceAtPosition(const Vector2D& force, const Vector2D& worldPos) override;
    void preStep(double time, double dt) override;

    // Fluid under the piston may not cross its lower face; `clearance` is face-to-first-particle distance.
    void setLidSeal(double clearance) { sealLid_ = true; lidClearance_ = clearance; }
    void constrainFluid(const std::vector<std::shared_ptr<Particle>>& fluid) override;
    bool participatesInCollisions() const override { return false; }
    // Never static: the solver must keep integrating it even when fully locked
    bool isStatic() const override { return false; }

private:
    DirectionLocks locks_;
    Vector2D load_{0.0, 0.0};
    double ptoDamping_ = 0.0;
    double harvested_ = 0.0;
    double minY_ = -std::numeric_limits<double>::infinity();
    double maxY_ = std::numeric_limits<double>::infinity();

    double releaseTime_ = -1.0;
    DirectionLocks locksAfterRelease_;
    bool released_ = false;

    bool sealLid_ = false;
    double lidClearance_ = 0.0;

    bool flywheel_ = false;
    double crankRadius_ = 0.0;
    double crankCentreY_ = 0.0;
    double inertia_ = 1.0;
    double theta_ = 0.0;
    double omega_ = 0.0;
    double shaftDamping_ = 0.0;
    double gasForceY_ = 0.0;
};
