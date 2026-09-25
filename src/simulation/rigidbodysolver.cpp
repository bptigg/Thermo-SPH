#include "rigidbodysolver.h"
#include <cmath>

RigidBodySolver::RigidBodySolver(RigidBodySolverParameters params) : params_(params) {}

void RigidBodySolver::solveCollisions(std::vector<std::shared_ptr<RigidObject>>& bodies, double dt) {
    size_t numBodies = bodies.size();
    double minDist = 2.0 * params_.particleRadius;
    double minDistSq = minDist * minDist;

    for (size_t a = 0; a < numBodies; ++a) {
        for (size_t b = a + 1; b < numBodies; ++b) {
            auto& bodyA = bodies[a];
            auto& bodyB = bodies[b];

            if (bodyA->isStatic() && bodyB->isStatic()) continue;

            // Broadphase AABB overlap test
            if (!bodyA->getAABB().overlaps(bodyB->getAABB())) continue;

            // Narrowphase constituent particle contact detection
            for (const auto& pA : bodyA->getParticles()) {
                for (const auto& pB : bodyB->getParticles()) {
                    Vector2D r = pA->pos - pB->pos;
                    double distSq = r.normSq();

                    if (distSq < minDistSq && distSq > 1e-12) {
                        double dist = std::sqrt(distSq);
                        Vector2D normal = r / dist;

                        Vector2D relVel = pA->vel - pB->vel;
                        double velAlongNormal = relVel.dot(normal);

                        if (velAlongNormal < 0.0) {
                            double impulseMag = -(1.0 + params_.restitution) * velAlongNormal;
                            double invMassA = bodyA->isStatic() ? 0.0 : (1.0 / bodyA->getTotalMass());
                            double invMassB = bodyB->isStatic() ? 0.0 : (1.0 / bodyB->getTotalMass());

                            if (invMassA + invMassB > 0.0) {
                                impulseMag /= (invMassA + invMassB);
                                Vector2D impulse = normal * impulseMag;

                                bodyA->addForceAtPosition(impulse / dt, pA->pos);
                                bodyB->addForceAtPosition(-impulse / dt, pB->pos);
                            }
                        }
                    }
                }
            }
        }
    }
}

void RigidBodySolver::integrate(std::vector<std::shared_ptr<RigidObject>>& bodies, double dt) {
    for (auto& body : bodies) {
        if (body->isStatic() || body->getTotalMass() <= 0.0) continue;

        Vector2D totalForce = body->getForceAccumulator();
        if (params_.enableGravity) {
            totalForce += params_.gravity * body->getTotalMass();
        }

        // Translation
        Vector2D linearAccel = totalForce / body->getTotalMass();
        if (body->isXLocked()) linearAccel.x = 0.0;
        if (body->isYLocked()) linearAccel.y = 0.0;

        Vector2D vel = body->getLinearVel() + linearAccel * dt;
        if (body->isXLocked()) vel.x = 0.0;
        if (body->isYLocked()) vel.y = 0.0;

        body->setLinearVel(vel);
        body->setCenterOfMass(body->getCenterOfMass() + vel * dt);

        // Rotation
        if (!body->isRotationLocked() && body->getInertia() > 0.0) {
            double angularAccel = body->getTorqueAccumulator() / body->getInertia();
            double omega = body->getAngularVel() + angularAccel * dt;
            body->setAngularVel(omega);
            body->setAngle(body->getAngle() + omega * dt);
        } else {
            body->setAngularVel(0.0);
        }

        // Sync constituent particles
        body->updateParticlePositions();
    }
}