#pragma once
#include "initialconditions.h"

class ThermalEquilibriumIC : public InitialConditions {
public:
    struct Parameters {
        Vector2D domainMin = Vector2D(-0.5, -0.5);
        Vector2D domainMax = Vector2D(0.5, 0.5);
        double spacing = 0.02;
        double density = 0.1;
        double uHot  = 1.0; // High thermal energy region
        double uCold = 0.1;  // Low thermal energy region

        bool lw = true;
        bool rw = true;
        bool tw = true;
        bool bw = true;

        int wallLayers = 3;
    };

    ThermalEquilibriumIC();
    explicit ThermalEquilibriumIC(Parameters params);

    std::vector<std::shared_ptr<Particle>> generateParticles() override;
    std::vector<std::shared_ptr<RigidObject>> getRigidObjects() override;

private:
    Parameters params_;
    std::vector<std::shared_ptr<Particle>> generatedParticles_;
    std::vector<std::shared_ptr<RigidObject>> generatedRigidObjects_;
    bool built_ = false;

    void buildScenario();
};