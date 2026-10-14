#pragma once
#include "initialconditions.h"
#include "Piston.h"

// Gas in an open-topped box closed by a piston.
// ConstantVelocity: the piston is driven down kinematically, then released so the gas pushes it back up.
// ConstantForce: a fixed downward force acts on a free piston that cannot rise above its starting height.
class PistonCompressionIC : public InitialConditions {
public:
    enum class DriveMode { ConstantVelocity, ConstantForce, Flywheel };

    struct Parameters {
        DriveMode driveMode = DriveMode::ConstantVelocity;
        Vector2D domainMin = Vector2D(0.0, 0.0);
        Vector2D domainMax = Vector2D(0.5, 1.0);
        int wallLayers = 3;
        double spacing = 0.02;

        double gasDensity = 1.0;
        double gasInternalEnergy = 0.5;
        double gasHeight = 0.6;

        int pistonLayers = 3;
        double pistonDensityRatio = 5.0;
        double compressionSpeed = 0.05;
        double compressionTime = 6.0;

        // Downward load on the released piston as a fraction of the initial gas force on it
        double loadFactor = 0.6;

        // ConstantForce mode: downward force as a multiple of the initial gas force on the piston
        double constantForceFactor = 0.8;
        // Lowest allowed piston centre-of-mass height
        double pistonMinY = 0.1;
        // ConstantForce mode: how far above its starting height the piston may rise
        double pistonCeilingAboveStart = 0.0;
        // Power take-off drag coefficient on the piston (0 = none)
        double powerTakeOff = 0.0;

        // Flywheel mode: slider-crank geometry and shaft
        double crankRadius = 0.12;
        double flywheelInertia = 0.02;
        double flywheelInitialOmega = 3.0;
        double shaftDamping = 0.0;
        double gamma = 1.4;

        bool enableGravity = false;
        Vector2D gravity = Vector2D(0.0, 0.0);
    };

    PistonCompressionIC();
    explicit PistonCompressionIC(Parameters params);

    std::vector<std::shared_ptr<Particle>> generateParticles() override;
    std::vector<std::shared_ptr<RigidObject>> getRigidObjects() override;

    std::function<double()> getThermalStateProvider() const override {
        if (piston_) {
            return [weakPiston = std::weak_ptr<RigidObject>(piston_)]() -> double {
                if (auto obj = weakPiston.lock()) {
                    return obj->getThermalTriggerValue(); // Calls polymorphic hook on RigidObject[cite: 10]
                }
                return 0.0;
            };
        }
        return nullptr;
    }

    std::shared_ptr<Piston> getPiston();
    const Parameters& getParameters() const { return params_; }

private:
    Parameters params_;
    std::vector<std::shared_ptr<Particle>> generatedParticles_;
    std::vector<std::shared_ptr<RigidObject>> generatedRigidObjects_;
    std::shared_ptr<Piston> piston_;
    bool built_ = false;

    void buildScenario();
};
