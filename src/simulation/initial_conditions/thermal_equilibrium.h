#pragma once
#include "initialconditions.h"

class ThermalEquilibriumIC : public InitialConditions {
public:
    struct Parameters {
        Vector2D domainMin = Vector2D(-0.5, -0.5);
        Vector2D domainMax = Vector2D(0.5, 0.5);
        double spacing = 0.02;
        double density = 1.0;
        double uHot  = 10.0; // High thermal energy region
        double uCold = 1.0;  // Low thermal energy region
    };

    ThermalEquilibriumIC();
    explicit ThermalEquilibriumIC(Parameters params);

    std::vector<std::shared_ptr<Particle>> generateParticles() override;

private:
    Parameters params_;
};