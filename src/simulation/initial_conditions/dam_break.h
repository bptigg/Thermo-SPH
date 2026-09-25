#pragma once
#include "initialconditions.h"

class DamBreakIC : public InitialConditions {
public:
    struct Parameters {
        // Domain & Boundary parameters
        Vector2D domainMin = Vector2D(0.0, 0.0);
        Vector2D domainMax = Vector2D(1.2, 1.0);
        int wallLayers = 3;                       // 3 particle layers thick to prevent leakage
        double spacing = 0.02;
        double fluidDensity = 1.0;
        double initialInternalEnergy = 1.0;

        // Fluid Column (placed behind the dam, height = 0.70m)
        Vector2D fluidBoxMin = Vector2D(0.04, 0.04);
        Vector2D fluidBoxMax = Vector2D(0.38, 0.70);

        // Static Dam Barrier Wall (crest height = 0.35m < fluid height 0.70m)
        Vector2D damPos = Vector2D(0.40, 0.0);
        Vector2D damSize = Vector2D(0.04, 0.35);

        // Optional Dynamic Rigid Block downstream in the spillway
        bool enableRigidObstacle = true;
        Vector2D obstaclePos = Vector2D(0.80, 0.04);
        Vector2D obstacleSize = Vector2D(0.10, 0.20);
        double obstacleParticleMass = 0.05;
    };

    DamBreakIC();
    explicit DamBreakIC(Parameters params);

    std::vector<std::shared_ptr<Particle>> generateParticles() override;
    std::vector<std::shared_ptr<RigidObject>> getRigidObjects() override;

private:
    Parameters params_;
    std::vector<std::shared_ptr<Particle>> generatedParticles_;
    std::vector<std::shared_ptr<RigidObject>> generatedRigidObjects_;
    bool built_ = false;

    void buildScenario();
};