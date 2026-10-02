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

    [[nodiscard]] virtual double computeTemperature(double internalEnergy) const = 0;
    [[nodiscard]] virtual double computeInternalEnergy(double temperature) const = 0;
    [[nodiscard]] virtual double getCv() const = 0;
};

class IdealGasEOS : public EquationOfState {
private:
    double gamma_;
    double Cv_; // Specific heat at constant volume

public:
    explicit IdealGasEOS(double gamma = 1.4) : gamma_(gamma), Cv_(1.0/(gamma-1.0)) {}

    double computePressure(double density, double internalEnergy) const override {
        double p = (gamma_ - 1.0) * density * std::max(internalEnergy, 0.0);
        return std::max(0.0, p); // Prevent non-physical negative pressures
    }

    double computeSoundSpeed(double density, double internalEnergy) const override {
        return std::sqrt(gamma_ * (gamma_ - 1.0) * std::max(internalEnergy, 1e-6));
    }

    double computeTemperature(double internalEnergy) const override {
        // T = u / c_v = (gamma - 1) * u
        return std::max(internalEnergy, 0.0) / Cv_;
    }

    double computeInternalEnergy(double temperature) const override {
        // u = c_v * T = T / (gamma - 1)
        return Cv_ * std::max(temperature, 0.0);
    }

    [[nodiscard]] double getCv() const override { return Cv_; }
    [[nodiscard]] double getGamma() const { return gamma_; }
};

class TaitEOS : public EquationOfState {
private:
    double gamma_;        // Stiffness exponent (typically 7.0 for water)
    double rho0_;         // Reference rest density
    double soundSpeed0_;  // Reference speed of sound
    double B_;            // Stiffness coefficient: B = (rho0 * c_s0^2) / gamma
    double Cv_;

public:
    TaitEOS(double rho0 = 1.0, double soundSpeed0 = 20.0, double gamma = 7.0)
        : gamma_(gamma), rho0_(rho0), soundSpeed0_(soundSpeed0)
    {
        B_ = (rho0_ * soundSpeed0_ * soundSpeed0_) / gamma_;
        Cv_ = (1.0 / (gamma - 1.0));
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

    double computeTemperature(double internalEnergy) const override {
        return std::max(internalEnergy, 0.0) / Cv_;
    }

    double computeInternalEnergy(double temperature) const override {
        return Cv_ * std::max(temperature, 0.0);
    }

    [[nodiscard]] double getCv() const override { return Cv_; }
    [[nodiscard]] double getGamma() const { return gamma_; }
};