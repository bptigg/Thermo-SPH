#pragma once
#include "initialconditions.h"

class SodShockTubeIC : public InitialConditions {
public:
    struct Parameters {
        double xMin = -1.0, xMax = 1.0;
        double yMin = -0.2, yMax = 0.2;
        double spacing = 0.02;
        
        // High density / pressure / energy left state
        double rhoLeft = 1.0;
        double uLeft   = 2.5; 
        
        // Low density / pressure / energy right state
        double rhoRight = 0.125;
        double uRight   = 2.0;
    };

    SodShockTubeIC();
    explicit SodShockTubeIC(Parameters params);

    std::vector<std::shared_ptr<Particle>> generateParticles() override;

private:
    Parameters params_;
};