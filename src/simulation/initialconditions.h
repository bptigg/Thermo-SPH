#pragma once

#include <vector>
#include <memory>
#include <functional>
#include "Particle.h"
#include "RigidBody.h"
#include "Vector2D.h"

class InitialConditions {
public:
    virtual ~InitialConditions() = default;

    virtual std::vector<std::shared_ptr<Particle>> generateParticles() = 0;
    virtual std::vector<std::shared_ptr<RigidObject>> getRigidObjects() { return {}; }
};

class CustomIC : public InitialConditions {
public:
    CustomIC() = default;

    void addParticle(std::shared_ptr<Particle> p);
    void addRigidObject(std::shared_ptr<RigidObject> obj);

    std::vector<std::shared_ptr<Particle>> generateParticles() override;
    std::vector<std::shared_ptr<RigidObject>> getRigidObjects() override;

private:
    std::vector<std::shared_ptr<Particle>> particles_;
    std::vector<std::shared_ptr<RigidObject>> rigidObjects_;
};

