#include "rigidbodysolver.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <map>

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

    //double minDist = 2.0 * params_.particleRadius;
    //double minDistSq = minDist * minDist;

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
        const Particle* pA; 
        const Particle* pB;
        Vector2D normal;
        double currentDistSq;
        double velAlongNormal;
    };

    struct SAPEndElement {
        size_t bodyIdx;
        double value;
        bool isMin;
    };

    std::vector<SAPEndElement> endpoints;
    endpoints.reserve(numBodies * 2);

    for (size_t i = 0; i < numBodies; ++i) {
        endpoints.push_back({i, bodyCache[i].sweptAABB.min.x, true});
        endpoints.push_back({i, bodyCache[i].sweptAABB.max.x, false});
    }

    std::sort(endpoints.begin(), endpoints.end(), [](const SAPEndElement& a, const SAPEndElement& b) {
        if (a.value == b.value) {
            return a.isMin && !b.isMin; // Process min endpoints before max on ties
        }
        return a.value < b.value;
    });

    std::vector<std::pair<size_t, size_t>> broadphasePairs;
    std::vector<size_t> activeBodies;

    for (const auto& endpoint : endpoints) {
        if (endpoint.isMin) {
            for (size_t activeIdx : activeBodies) {
                if (bodyCache[endpoint.bodyIdx].isStatic && bodyCache[activeIdx].isStatic) {
                    continue;
                }

                // Check 2D AABB overlap
                if (bodyCache[endpoint.bodyIdx].sweptAABB.overlaps(bodyCache[activeIdx].sweptAABB)) {
                    broadphasePairs.emplace_back(endpoint.bodyIdx, activeIdx);
                }
            }
            activeBodies.push_back(endpoint.bodyIdx);
        } 
        else {
            auto it = std::find(activeBodies.begin(), activeBodies.end(), endpoint.bodyIdx);
            if (it != activeBodies.end()) {
                activeBodies.erase(it);
            }
        }
    }

    if (broadphasePairs.size() ==0) {
        return false;
    }


    std::vector<Contact> bodyPairContacts;

    auto narrowPhaseDection = [](const std::vector<std::pair<size_t, size_t>>& broadphasePairs, 
        const std::vector<std::shared_ptr<RigidObject>>& bodies,
        double dt,
        double particleRadius) 
    {
        std::vector<Contact> contacts;
        double minDist = 2.0 * particleRadius;
        double minDistSq = minDist * minDist;
        for (const auto& [a, b] : broadphasePairs) 
        {
            const auto& bodyA = bodies[a];
            const auto& bodyB = bodies[b];
            const auto& particlesA = bodyA->getParticles();
            const auto& particlesB = bodyB->getParticles();

            for (const auto& pA : particlesA) 
            {
                for (const auto& pB : particlesB) 
                {
                    Vector2D relPos = pA->pos - pB->pos;
                    double currentDistSq = relPos.normSq();

                    // Ignore virtually identical positions to avoid division by zero
                    if (currentDistSq < 1e-12) continue;

                    Vector2D relVelAtImpact = pA->vel - pB->vel;
                    Vector2D relVel = relVelAtImpact * dt;
                    double aVal = relVel.normSq();

                    bool isColliding = false;
                    Vector2D normal(0.0, 1.0);

                    // Continuous Collision Test (Swept Sphere-Sphere)
                    if (aVal > 1e-12) 
                    {
                        double h = relPos.dot(relVel);
                        double cTerm = currentDistSq - minDistSq;
                        double disc = h * h - aVal * cTerm;

                        if (disc >= 0.0) 
                        {
                            double sqrtDisc = std::sqrt(disc);
                            double t1 = (-h - sqrtDisc) / aVal;
                            double t2 = (-h + sqrtDisc) / aVal;

                            double tValid = -1.0;
                            if (t1 >= 0.0 && t1 <= 1.0) tValid = t1;
                            if (t2 >= 0.0 && t2 <= 1.0 && (tValid < 0.0 || t2 < tValid)) tValid = t2;

                            if (tValid >= 0.0) {
                                Vector2D relAtImpact = relPos + relVel * tValid;
                                double relLen = relAtImpact.norm();
                                if (relLen > 1e-12) {
                                    normal = relAtImpact / relLen;
                                    isColliding = true;
                                }
                            }
                        }
                    }

                    // Fallback Static Overlap Test
                    if (!isColliding && currentDistSq <= minDistSq) {
                        double relLen = std::sqrt(currentDistSq);
                        if (relLen > 1e-12) {
                            normal = relPos / relLen;
                            isColliding = true;
                        }
                    }

                    if (!isColliding) continue;

                    // Relative velocity along normal direction
                    double velAlongNormal = relVelAtImpact.dot(normal);

                    // Skip if objects are moving apart and not overlapping
                    if (velAlongNormal >= 0.0 && currentDistSq > minDistSq) {
                        continue;
                    }
                    contacts.push_back({a, b, pA.get(), pB.get(), normal, currentDistSq, velAlongNormal});
                }
            }
        }
        return contacts;
    };

    auto contacts = narrowPhaseDection(broadphasePairs, bodies, dt, params_.particleRadius);
    if (contacts.empty()) {
        return false;
    }

    auto applyPenaltyForces = [&](const std::vector<Contact>& contacts) {
        double minOverlapDist = 2.0 * params_.particleRadius;

        std::map<std::pair<size_t, size_t>, size_t> pairContactCounts;
        for (const auto& contact : contacts) {
            size_t minIdx = std::min(contact.bodyAIdx, contact.bodyBIdx);
            size_t maxIdx = std::max(contact.bodyAIdx, contact.bodyBIdx);
            pairContactCounts[{minIdx, maxIdx}]++;
        }

        // Relaxation factor: 0.01 to 0.05 resolves overlap smoothly over 20-100 steps
        const double beta = 0.2; 
        const double maxAccel = 1000.0; // Maximum allowed contact acceleration (m/s^2)

        for (const auto& contact : contacts) {
            size_t idxA = contact.bodyAIdx;
            size_t idxB = contact.bodyBIdx;

            const auto& cacheA = bodyCache[idxA];
            const auto& cacheB = bodyCache[idxB];

            // Lever arms relative to body COM
            Vector2D rA = contact.pA->pos - bodies[idxA]->getCenterOfMass();
            Vector2D rB = contact.pB->pos - bodies[idxB]->getCenterOfMass();

            // 2D Cross product (r x n)
            double rA_cross_n = rA.x * contact.normal.y - rA.y * contact.normal.x;
            double rB_cross_n = rB.x * contact.normal.y - rB.y * contact.normal.x;

            // Effective inverse mass along contact normal
            double kNormal = cacheA.invMass + cacheB.invMass +
                             (rA_cross_n * rA_cross_n) * cacheA.invInertia +
                             (rB_cross_n * rB_cross_n) * cacheB.invInertia;

            if (kNormal < 1e-12) continue; // Both bodies static

            double effectiveMass = 1.0 / kNormal;

            // Penetration depth
            double currentDist = std::sqrt(contact.currentDistSq);
            double penetration = minOverlapDist - currentDist;
            if (penetration <= 0.0 && contact.velAlongNormal >= 0.0) continue;

            // Get number of contacts sharing this body pair
            size_t minIdx = std::min(idxA, idxB);
            size_t maxIdx = std::max(idxA, idxB);
            size_t numContacts = pairContactCounts[{minIdx, maxIdx}];

            // 2. Scale stiffness by relaxation factor beta and divide by active contact count
            double kSpring = params_.penaltyStiffness > 0.0 
                             ? (params_.penaltyStiffness / numContacts)
                             : (beta * effectiveMass / (dt * dt * numContacts));

            double e = std::clamp(params_.restitution, 0.0, 1.0);
            double logE = (e < 1e-4) ? -9.21 : std::log(e);
            double dampingRatio = -logE / std::sqrt(M_PI * M_PI + logE * logE);
            double cDamping = 2.0 * dampingRatio * std::sqrt(kSpring * (effectiveMass / numContacts));

            // Spring-Damper Normal Force
            double springForce = kSpring * std::max(0.0, penetration);
            double dampingForce = -cDamping * contact.velAlongNormal;
            double totalNormalForceMag = std::max(0.0, springForce + dampingForce);

            // 3. Clamp force so max acceleration per contact does not explode
            double maxForce = effectiveMass * maxAccel / numContacts;
            totalNormalForceMag = std::min(totalNormalForceMag, maxForce);

            Vector2D forceA = contact.normal * totalNormalForceMag;

            // Apply equal and opposite reaction forces
            if (!cacheA.isStatic) {
                bodies[idxA]->addForceAtPosition(forceA, contact.pA->pos);
            }
            if (!cacheB.isStatic) {
                bodies[idxB]->addForceAtPosition(-forceA, contact.pB->pos);
            }
        }
    };

    applyPenaltyForces(contacts);
    return true;
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