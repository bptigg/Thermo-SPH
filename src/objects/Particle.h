#pragma once
#include "Vector2D.h"
#include "IBoundary.h"
#include <memory>
#include <vector>


enum class MotionType
{
    DEFAULT = 0,
    STATIC, //fixed walls
    NONSTATIC, //kinematics - has a prescribed motion
    DYNAMIC //particle, breakable wall, pistons driven by force..
};

class Particle 
{
public:
    int id;
    Vector2D pos;
    Vector2D vel;
    Vector2D accel;

    double mass = 0.0;
    double density = 0.0;
    double u = 0.0; //specific internal energy
    double dudt = 0.0; // rate of change of internal energy
    double h = 0.0; //smoothing length 

protected:
    std::shared_ptr<IBoundary> boundaryComponent_ = nullptr;
public:

    Particle(int id, Vector2D pos, Vector2D vel, double m, double r, double u)
        : id(id), pos(pos), vel(vel), accel(0.0,0.0),
          mass(m), density(r), u(u) {}

    virtual ~Particle() = default;

    virtual void kickHalf(double) = 0;
    virtual void drift(double) = 0;

    virtual MotionType getMotionType() const = 0;

    virtual bool isFluid() const { return false; }
    
    bool isStatic() const { return getMotionType() == MotionType::STATIC; }
    bool isNonStatic() const { return getMotionType() == MotionType::NONSTATIC; }
    bool isDynamic() const { return getMotionType() == MotionType::DYNAMIC; }

    void setBoundaryComponent(std::shared_ptr<IBoundary> boundary) {
        boundaryComponent_ = std::move(boundary);
    }

    bool hasBoundary() const { 
        return boundaryComponent_ != nullptr; 
    }

    std::vector<std::shared_ptr<Particle>> generateGhosts(
        const std::vector<std::shared_ptr<Particle>>& fluidParticles, 
        double supportRadius) const 
    {
        if (!boundaryComponent_) return {};
        return boundaryComponent_->generateGhosts(pos, vel, u, fluidParticles, supportRadius);
    }

    Vector2D getMomentum() const { 
        return vel * mass; 
    }

    double getKineticEnergy() const { 
        return 0.5 * mass * vel.normSq(); 
    }

    double getInternalEnergy() const { 
        return mass * u; 
    }

    double getTotalEnergy() const { 
        return getKineticEnergy() + getInternalEnergy(); 
    }
};

class FluidParticle : public Particle 
{
public:
    using Particle::Particle;

    void kickHalf(double dt) override {
        if (!isDynamic()) return;
        vel += accel * (0.5 * dt);
        u   += dudt  * (0.5 * dt);
        if (u < 1e-5) u = 1e-5;
    }

    void drift(double dt) override {
        if (!isDynamic()) return;
        pos += vel * dt;
    }

    MotionType getMotionType() const override { return MotionType::DYNAMIC; }
    bool isFluid() const override { return true; }
};

class SolidParticle : public Particle 
{
private:
    MotionType motionType_;

public:
    SolidParticle(int id, Vector2D pos, Vector2D vel, double m, double r, double u, MotionType type = MotionType::STATIC)
        : Particle(id, pos, vel, m, r, u), motionType_(type) {}

    bool isFluid() const override { return false; }
    MotionType getMotionType() const override { return motionType_; }
    void setMotionType(MotionType type) { motionType_ = type; }

    void kickHalf(double) override {}
    void drift(double) override {}
};


