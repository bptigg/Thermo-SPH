#pragma once
#include "initialconditions.h"

// Closed box of ideal gas with a small, extremely over-pressured disc at its centre (2D Sedov-Taylor blast).
// The blast drives a strong cylindrical shock outward that reflects off the walls.
class ExplosionIC : public InitialConditions {
public:
    struct Parameters {
        Vector2D domainMin = Vector2D(0.0, 0.0);
        Vector2D domainMax = Vector2D(1.0, 1.0);
        int wallLayers = 3;
        double spacing = 0.02;

        double gasDensity = 1.0;
        double ambientInternalEnergy = 1.0;

        Vector2D blastCentre = Vector2D(0.5, 0.5);
        double blastRadius = 0.08;
        // Internal energy inside the blast relative to the ambient gas (the pressure ratio at equal density)
        double blastEnergyRatio = 2000.0;

        double gamma = 1.4;
        bool enableGravity = false;
        Vector2D gravity = Vector2D(0.0, 0.0);
    };

    ExplosionIC();
    explicit ExplosionIC(Parameters params);

    std::vector<std::shared_ptr<Particle>> generateParticles() override;
    std::vector<std::shared_ptr<RigidObject>> getRigidObjects() override;

    const Parameters& getParameters() const { return params_; }
    // Total internal energy released by the blast above ambient
    double getBlastEnergy();

private:
    Parameters params_;
    std::vector<std::shared_ptr<Particle>> generatedParticles_;
    std::vector<std::shared_ptr<RigidObject>> generatedRigidObjects_;
    double blastEnergy_ = 0.0;
    bool built_ = false;

    void buildScenario();
};
