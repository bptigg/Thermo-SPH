#include "dam_break.h"
#include <memory>

DamBreakIC::DamBreakIC() : DamBreakIC(Parameters{}) {}
DamBreakIC::DamBreakIC(Parameters params) : params_(params) {}

void DamBreakIC::buildScenario() {
    if (built_) return;

    generatedParticles_.clear();
    generatedRigidObjects_.clear();
    int id = 0;

    const double particleMass = params_.fluidDensity * params_.spacing * params_.spacing;
    const double h = params_.spacing;
    const int layers = params_.wallLayers;

    // -------------------------------------------------------------
    // 1. Static Bottom Floor Wall (placed in negative Y layers)
    // -------------------------------------------------------------
    for (double x = params_.domainMin.x - layers * h; x <= params_.domainMax.x + layers * h; x += h) {
        for (int layer = 1; layer <= layers; ++layer) {
            double y = params_.domainMin.y - layer * h;
            auto wallP = std::make_shared<SolidParticle>(
                id++,
                Vector2D(x, y),
                Vector2D(0.0, 0.0),
                particleMass,
                params_.fluidDensity,
                params_.initialInternalEnergy,
                MotionType::STATIC
            );
            generatedParticles_.push_back(wallP);
        }
    }

    // -------------------------------------------------------------
    // 2. Static Left Containment Wall (placed in negative X layers)
    // -------------------------------------------------------------
    for (double y = params_.domainMin.y; y <= params_.domainMax.y; y += h) {
        for (int layer = 1; layer <= layers; ++layer) {
            double x = params_.domainMin.x - layer * h;
            auto wallP = std::make_shared<SolidParticle>(
                id++,
                Vector2D(x, y),
                Vector2D(0.0, 0.0),
                particleMass,
                params_.fluidDensity,
                params_.initialInternalEnergy,
                MotionType::STATIC
            );
            generatedParticles_.push_back(wallP);
        }
    }

    // -------------------------------------------------------------
    // 3. Static Right Containment Wall (retains gravity surge)
    // -------------------------------------------------------------
    if (params_.enableRightWall) {
        for (double y = params_.domainMin.y; y <= params_.domainMax.y; y += h) {
            for (int layer = 1; layer <= layers; ++layer) {
                double x = params_.domainMax.x + (layer - 1) * h;
                auto wallP = std::make_shared<SolidParticle>(
                    id++,
                    Vector2D(x, y),
                    Vector2D(0.0, 0.0),
                    particleMass,
                    params_.fluidDensity,
                    params_.initialInternalEnergy,
                    MotionType::STATIC
                );
                generatedParticles_.push_back(wallP);
            }
        }
    }

    // -------------------------------------------------------------
    // 4. Static Dam Barrier Wall (obstacle over which fluid spills)
    // -------------------------------------------------------------
    double damXMax = params_.damPos.x + params_.damSize.x;
    double damYMax = params_.damPos.y + params_.damSize.y;

    for (double x = params_.damPos.x; x <= damXMax; x += h) {
        for (double y = params_.damPos.y; y <= damYMax; y += h) {
            auto damP = std::make_shared<SolidParticle>(
                id++,
                Vector2D(x, y),
                Vector2D(0.0, 0.0),
                particleMass,
                params_.fluidDensity,
                params_.initialInternalEnergy,
                MotionType::STATIC
            );
            generatedParticles_.push_back(damP);
        }
    }

    // -------------------------------------------------------------
    // 5. Dynamic Fluid Column (collapses under gravity)
    // -------------------------------------------------------------
    for (double x = params_.fluidBoxMin.x; x <= params_.fluidBoxMax.x; x += h) {
        for (double y = params_.fluidBoxMin.y; y <= params_.fluidBoxMax.y; y += h) {
            auto fluidP = std::make_shared<FluidParticle>(
                id++,
                Vector2D(x, y),
                Vector2D(0.0, 0.0), // Fluid starts at rest; gravity accelerates it downwards
                particleMass,
                params_.fluidDensity,
                params_.initialInternalEnergy
            );
            generatedParticles_.push_back(fluidP);
        }
    }

    // -------------------------------------------------------------
    // 6. Dynamic Rigid Obstacle Downstream
    // -------------------------------------------------------------
    if (params_.enableRigidObstacle) {
        auto rigidBody = std::make_shared<RigidObject>();

        double oxMin = params_.obstaclePos.x;
        double oyMin = params_.obstaclePos.y;
        double oxMax = oxMin + params_.obstacleSize.x;
        double oyMax = oyMin + params_.obstacleSize.y;

        double obstacleParticleMass = (params_.fluidDensity * params_.obstacleDensityRatio) * h * h;

        for (double x = oxMin; x <= oxMax; x += h) {
            for (double y = oyMin; y <= oyMax; y += h) {
                auto solidP = std::make_shared<SolidParticle>(
                    id++,
                    Vector2D(x, y),
                    Vector2D(0.0, 0.0),
                    obstacleParticleMass,
                    params_.fluidDensity,
                    params_.initialInternalEnergy,
                    MotionType::DYNAMIC
                );

                // Added only to rigidBody to avoid double-counting in engine/metrics
                rigidBody->addParticle(solidP);
            }
        }

        generatedRigidObjects_.push_back(rigidBody);
    }

    built_ = true;
}

std::vector<std::shared_ptr<Particle>> DamBreakIC::generateParticles() {
    buildScenario();
    return generatedParticles_;
}

std::vector<std::shared_ptr<RigidObject>> DamBreakIC::getRigidObjects() {
    buildScenario();
    return generatedRigidObjects_;
}