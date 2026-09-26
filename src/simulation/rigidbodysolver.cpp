#include "rigidbodysolver.h"
#include <cmath>
#include <vector>
#include <algorithm>

RigidBodySolver::RigidBodySolver(RigidBodySolverParameters params) : params_(params) {}

void RigidBodySolver::integrateVelocities(std::vector<std::shared_ptr<RigidObject>>& bodies, double dt) {
    for (auto& body : bodies) {
        if (body->isStatic() || body->getTotalMass() <= 0.0) continue;

        // Apply external gravity if enabled
        if (params_.enableGravity) {
            body->addForceAtPosition(params_.gravity * body->getTotalMass(), body->getCenterOfMass());
        }

        // Advance linear velocity from forceAccumulator_
        Vector2D accel = body->getForceAccumulator() / body->getTotalMass();
        if (body->isXLocked()) accel.x = 0.0;
        if (body->isYLocked()) accel.y = 0.0;

        Vector2D vel = body->getLinearVel() + accel * dt;

        if (body->isXLocked()) vel.x = 0.0;
        if (body->isYLocked()) vel.y = 0.0;
        body->setLinearVel(vel);

        // Advance angular velocity from torqueAccumulator_
        if (!body->isRotationLocked() && body->getInertia() > 0.0) {
            double alpha = body->getTorqueAccumulator() / body->getInertia();
            body->setAngularVel(body->getAngularVel() + alpha * dt);
        } else {
            body->setAngularVel(0.0);
        }
    }
}

bool RigidBodySolver::solveCollisions(std::vector<std::shared_ptr<RigidObject>>& bodies, double dt) {
    size_t numBodies = bodies.size();
    if (numBodies == 0) return false;

    double minDist = 2.0 * params_.particleRadius;
    double minDistSq = minDist * minDist;

    struct BodyCache {
        bool isStatic;
        double invMass;
        double invInertia;
        AABB aabb;
        AABB sweptAABB;
    };

    std::vector<BodyCache> bodyCache(numBodies);

    for (size_t i = 0; i < numBodies; ++i) {
        const auto& body = bodies[i];
        bodyCache[i].isStatic = body->isStatic();
        bodyCache[i].invMass = bodyCache[i].isStatic ? 0.0 : (1.0 / body->getTotalMass());
        bodyCache[i].invInertia = (bodyCache[i].isStatic || body->isRotationLocked() || body->getInertia() <= 0.0)
                                  ? 0.0 : (1.0 / body->getInertia());
        bodyCache[i].aabb = body->getAABB();

        AABB swept = bodyCache[i].aabb;
        if (!bodyCache[i].isStatic) {
            Vector2D velDt = body->getLinearVel() * dt;
            swept.min.x = std::min(swept.min.x, swept.min.x + velDt.x);
            swept.max.x = std::max(swept.max.x, swept.max.x + velDt.x);
            swept.min.y = std::min(swept.min.y, swept.min.y + velDt.y);
            swept.max.y = std::max(swept.max.y, swept.max.y + velDt.y);
        }
        bodyCache[i].sweptAABB = swept;
    }

    struct Contact {
        size_t bodyAIdx;
        size_t bodyBIdx;
        size_t particleAIdx;
        size_t particleBIdx;
        Vector2D normal;
        double currentDistSq;
        double velAlongNormal;
    };

    std::vector<Contact> bodyPairContacts;

    // Pass 1: Contact Detection & Pair Reduction (1 primary contact per body pair)
    for (size_t a = 0; a < numBodies; ++a) {
        for (size_t b = a + 1; b < numBodies; ++b) {
            if (bodyCache[a].isStatic && bodyCache[b].isStatic) continue;

            const AABB& aabbA = bodyCache[a].aabb;
            const AABB& aabbB = bodyCache[b].aabb;

            if (!aabbA.overlaps(aabbB)) {
                const AABB& sweptA = bodyCache[a].sweptAABB;
                const AABB& sweptB = bodyCache[b].sweptAABB;

                bool aabbMayIntersect = !(sweptA.max.x < sweptB.min.x || sweptB.max.x < sweptA.min.x ||
                                          sweptA.max.y < sweptB.min.y || sweptB.max.y < sweptA.min.y);
                if (!aabbMayIntersect) continue;
            }

            const auto& bodyA = bodies[a];
            const auto& bodyB = bodies[b];
            const auto& particlesA = bodyA->getParticles();
            const auto& particlesB = bodyB->getParticles();

            bool foundContact = false;
            Contact bestContact;
            double minFoundDistSq = 1e18; // Keep deepest penetration point

            for (size_t pIdxA = 0; pIdxA < particlesA.size(); ++pIdxA) {
                const auto& pA = particlesA[pIdxA];
                for (size_t pIdxB = 0; pIdxB < particlesB.size(); ++pIdxB) {
                    const auto& pB = particlesB[pIdxB];

                    Vector2D relStart = pA->pos - pB->pos;
                    double currentDistSq = relStart.normSq();

                    if (currentDistSq < minDistSq && currentDistSq < 1e-12) {
                        continue;
                    }

                    // Compute current particle velocities based on center of mass and angular velocity
                    Vector2D rA = pA->pos - bodyA->getCenterOfMass();
                    Vector2D rB = pB->pos - bodyB->getCenterOfMass();
                    Vector2D vA = bodyA->getLinearVel() + Vector2D(-bodyA->getAngularVel() * rA.y, bodyA->getAngularVel() * rA.x);
                    Vector2D vB = bodyB->getLinearVel() + Vector2D(-bodyB->getAngularVel() * rB.y, bodyB->getAngularVel() * rB.x);

                    Vector2D relVelAtImpact = vA - vB;
                    Vector2D relVel = relVelAtImpact * dt;

                    double aVal = relVel.dot(relVel);
                    bool intersectsDuringStep = false;
                    Vector2D normal(0.0, 1.0);

                    if (aVal > 1e-12) {
                        double h = relStart.dot(relVel);
                        double cTerm = currentDistSq - minDistSq;
                        double disc = h * h - aVal * cTerm;

                        if (disc >= 0.0) {
                            double sqrtDisc = std::sqrt(disc);
                            double t1 = (-h - sqrtDisc) / aVal;
                            double t2 = (-h + sqrtDisc) / aVal;

                            double firstValid = -1.0;
                            if (t1 >= 0.0 && t1 <= 1.0) firstValid = t1;
                            if (t2 >= 0.0 && t2 <= 1.0 && (firstValid < 0.0 || t2 < firstValid)) firstValid = t2;

                            if (firstValid >= 0.0) {
                                Vector2D relAtImpact = relStart + relVel * firstValid;
                                double relLen = relAtImpact.norm();
                                if (relLen > 1e-12) {
                                    normal = relAtImpact / relLen;
                                    intersectsDuringStep = true;
                                }
                            }
                        }
                    }

                    if (!intersectsDuringStep && currentDistSq <= minDistSq) {
                        double relLen = relStart.norm();
                        if (relLen > 1e-12) {
                            normal = relStart / relLen;
                            intersectsDuringStep = true;
                        }
                    }

                    if (!intersectsDuringStep) continue;

                    double velAlongNormal = relVelAtImpact.dot(normal);
                    if (velAlongNormal >= 0.0 && currentDistSq > minDistSq) {
                        continue;
                    }

                    if (currentDistSq < minFoundDistSq) {
                        minFoundDistSq = currentDistSq;
                        bestContact = {a, b, pIdxA, pIdxB, normal, currentDistSq, velAlongNormal};
                        foundContact = true;
                    }
                }
            }

            if (foundContact) {
                bodyPairContacts.push_back(bestContact);
            }
        }
    }
    if (bodyPairContacts.empty()) return false;

    // Pass 2: Calculate contact forces and record into forceAccumulator_
    int collisions = 0;
    for (const auto& contact : bodyPairContacts) {
        size_t a = contact.bodyAIdx;
        size_t b = contact.bodyBIdx;
        auto& bodyA = bodies[a];
        auto& bodyB = bodies[b];

        double invMassA = bodyCache[a].invMass;
        double invMassB = bodyCache[b].invMass;
        double invInertiaA = bodyCache[a].invInertia;
        double invInertiaB = bodyCache[b].invInertia;

        const auto& pA = bodyA->getParticles()[contact.particleAIdx];
        const auto& pB = bodyB->getParticles()[contact.particleBIdx];

        Vector2D rA = pA->pos - bodyA->getCenterOfMass();
        Vector2D rB = pB->pos - bodyB->getCenterOfMass();

        double rACrossN = rA.x * contact.normal.y - rA.y * contact.normal.x;
        double rBCrossN = rB.x * contact.normal.y - rB.y * contact.normal.x;

        double effectiveMass = invMassA + invMassB + 
                              (rACrossN * rACrossN) * invInertiaA + 
                              (rBCrossN * rBCrossN) * invInertiaB;

        if (effectiveMass > 0.0) {
            double impulseMag = -(1.0 + params_.restitution) * contact.velAlongNormal / effectiveMass;
            Vector2D impulse = contact.normal * impulseMag;

            // Convert impulse to contact force over timestep dt and accumulate
            Vector2D contactForce = impulse / dt;
            bodyA->addForceAtPosition(contactForce, pA->pos);
            bodyB->addForceAtPosition(-contactForce, pB->pos);

            // Penetration position correction
            double invMassSum = invMassA + invMassB;
            if (invMassSum > 0.0) {
                double dist = std::sqrt(contact.currentDistSq);
                double penetration = std::max(0.0, minDist - dist);
                if (penetration > 0.0) {
                    Vector2D correction = contact.normal * (penetration / invMassSum);
                    if (!bodyCache[a].isStatic) {
                        bodyA->setCenterOfMass(bodyA->getCenterOfMass() + correction * invMassA);
                        bodyA->updateParticlePositions();
                    }
                    if (!bodyCache[b].isStatic) {
                        bodyB->setCenterOfMass(bodyB->getCenterOfMass() - correction * invMassB);
                        bodyB->updateParticlePositions();
                    }
                }
            }
            collisions++;
        }
    }
    return collisions > 0;
}

void RigidBodySolver::integratePositions(std::vector<std::shared_ptr<RigidObject>>& bodies, double dt) {
    for (auto& body : bodies) {
        if (body->isStatic()) continue;

        // Advance center of mass position
        body->setCenterOfMass(body->getCenterOfMass() + body->getLinearVel() * dt);

        // Advance rotation angle
        if (!body->isRotationLocked()) {
            body->setAngle(body->getAngle() + body->getAngularVel() * dt);
        }

        // Synchronize constituent particle positions and velocities
        body->updateParticlePositions();
    }
}

bool RigidBodySolver::integrate(std::vector<std::shared_ptr<RigidObject>>& bodies, double dt) {
    // 1. Reset forces from the previous frame
    for (auto& body : bodies) {
        body->clearForces();
    }

    bool collided = solveCollisions(bodies, dt);
    integrateVelocities(bodies, dt);
    integratePositions(bodies, dt);

    return collided;
}