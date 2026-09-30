#pragma once
#include <algorithm>
#include <cmath>

// Abstract Equation of State interface
class EquationOfState {
public:
    virtual ~EquationOfState() = default;

    // Calculates pressure P given density (rho) and internal energy (u)
    [[nodiscard]] virtual double computePressure(double density, double internalEnergy) const = 0;

    // Calculates local speed of sound c_s given density (rho) and internal energy (u)
    [[nodiscard]] virtual double computeSoundSpeed(double density, double internalEnergy) const = 0;
};

class IdealGasEOS : public EquationOfState {
private:
    double gamma_;

public:
    explicit IdealGasEOS(double gamma = 1.4) : gamma_(gamma) {}

    double computePressure(double density, double internalEnergy) const override {
        double p = (gamma_ - 1.0) * density * std::max(internalEnergy, 0.0);
        return std::max(0.0, p); // Prevent non-physical negative pressures
    }

    double computeSoundSpeed(double density, double internalEnergy) const override {
        return std::sqrt(gamma_ * (gamma_ - 1.0) * std::max(internalEnergy, 1e-6));
    }

    [[nodiscard]] double getGamma() const { return gamma_; }
};

class TaitEOS : public EquationOfState {
private:
    double gamma_;        // Stiffness exponent (typically 7.0 for water)
    double rho0_;         // Reference rest density
    double soundSpeed0_;  // Reference speed of sound
    double B_;            // Stiffness coefficient: B = (rho0 * c_s0^2) / gamma

public:
    TaitEOS(double rho0 = 1.0, double soundSpeed0 = 20.0, double gamma = 7.0)
        : gamma_(gamma), rho0_(rho0), soundSpeed0_(soundSpeed0) 
    {
        B_ = (rho0_ * soundSpeed0_ * soundSpeed0_) / gamma_;
    }

    double computePressure(double density, double internalEnergy) const override {
        double relativeDensity = density / rho0_;
        double p = B_ * (std::pow(relativeDensity, gamma_) - 1.0);
        return std::max(0.0, p);
    }

    double computeSoundSpeed(double density, double internalEnergy) const override {
        double relativeDensity = std::max(density / rho0_, 1e-5);
        return soundSpeed0_ * std::pow(relativeDensity, 0.5 * (gamma_ - 1.0));
    }
};