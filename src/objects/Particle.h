#pragma once
#include "Vector2D.h"


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
};

class StaticParticle : public Particle 
{
public:
    using Particle::Particle;

    void kickHalf(double) override {}
    void drift(double) override {}

    MotionType getMotionType() const override { return MotionType::STATIC; }
};

class NonStaticParticle : public Particle {
protected:
    Vector2D prescribedVel_;

public:
    NonStaticParticle(int id, Vector2D pos, Vector2D prescribedVel, double mass, double density, double u)
        : Particle(id, pos, prescribedVel, mass, density, u), prescribedVel_(prescribedVel) {}

    void setPrescribedVelocity(Vector2D v) { 
        prescribedVel_ = v; 
        vel = v; 
    }

    void kickHalf(double) override {}

    void drift(double dt) override {
        pos += prescribedVel_ * dt;
    }

    MotionType getMotionType() const override { return MotionType::NONSTATIC; }
};

class DynamicParticle : public Particle {
public:
    using Particle::Particle;

    void kickHalf(double dt) override {
        vel += accel * (0.5 * dt);
        u += dudt * (0.5 * dt);
        if (u < 1e-4) u = 1e-4; 
    }

    void drift(double dt) override {
        pos += vel * dt;
    }

    MotionType getMotionType() const override { return MotionType::DYNAMIC; }
};


