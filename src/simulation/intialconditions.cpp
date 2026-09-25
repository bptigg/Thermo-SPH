#include "initialconditions.h"

void CustomIC::addParticle(std::shared_ptr<Particle> p) {
    particles_.push_back(p);
}

void CustomIC::addRigidObject(std::shared_ptr<RigidObject> obj) {
    rigidObjects_.push_back(obj);
    for (const auto& p : obj->getParticles()) {
        particles_.push_back(p);
    }
}

std::vector<std::shared_ptr<Particle>> CustomIC::generateParticles() {
    return particles_;
}

std::vector<std::shared_ptr<RigidObject>> CustomIC::getRigidObjects() {
    return rigidObjects_;
}

