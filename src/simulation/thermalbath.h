#pragma once
#include <memory>
#include <vector>
#include "Particle.h"
#include "Vector2D.h"

// Berendsen-style heat bath: relaxes the internal energy of nearby gas toward targetU over time tau.

enum class BathMode { DEFAULT, HOT_ONLY, COLD_ONLY };
enum class BathState { OFF, ON };

struct ThermalBath {
    bool enabled = false;
    double targetU = 0.5;
    double tau = 0.1;
    // Bath acts within this distance of the domain walls; <= 0 couples the whole gas.
    double wallThickness = 0.0;
    // Optional heat source: while a piston sits below hotBelowY the bath target becomes hotU (hotU <= 0 disables)
    double hotU = 0.0;
    double hotBelowY = 0.0;

    mutable BathMode mode = BathMode::DEFAULT;
    mutable BathState state = BathState::OFF;

    mutable double currentActiveU = 0.5;
    mutable bool isActive = false;
    
    std::function<double()> stateProvider = nullptr;

    void updateState() const {
        //if default the bath should fluctuate between hot and cold depending on the piston position
        //if cold only the bath should turn off when the piston is below hotBelowY
        //if hot only the bath should turn off when the piston is above hotBelowY
        state = BathState::ON;
        isActive = true;
        //currentActiveU = targetU;

        double triggerVal = (stateProvider ? stateProvider() : 0.0);

        if (mode == BathMode::DEFAULT) {
            if (hotU > 0.0 && triggerVal< hotBelowY) {
                currentActiveU = hotU;
            } else {
                currentActiveU = targetU;
            }
        } 
        else if (mode == BathMode::HOT_ONLY) {
            if (hotU > 0.0 && triggerVal< hotBelowY) {
                currentActiveU = hotU;
            } else {
                state = BathState::OFF;
                currentActiveU = 0.0; // bath off
            }
        } 
        else if (mode == BathMode::COLD_ONLY) {
            if (triggerVal< hotBelowY) {
                currentActiveU = targetU;
            } else {
                state = BathState::OFF;
                currentActiveU = 0.0; // bath off
            }
        }
    }

    double getCurrentActiveU() const {
        updateState();
        return currentActiveU;
    }

    double isActiveBath() const {
        updateState();
        return isActive;
    }

    bool contains(const Vector2D& pos, const Vector2D& domainMin, const Vector2D& domainMax) const {
        if (wallThickness <= 0.0) return true;
        return pos.x < domainMin.x + wallThickness || pos.x > domainMax.x - wallThickness ||
               pos.y < domainMin.y + wallThickness;
    }

    // Adds the relaxation rate to dudt and returns the heat flow into the gas (energy per time).
    double apply(const std::vector<std::shared_ptr<Particle>>& fluid,
                 const Vector2D& domainMin, const Vector2D& domainMax) {
        if (!enabled || tau <= 0.0) return 0.0;
        updateState();
        double power = 0.0;
        if (state == BathState::OFF) return 0.0;
        for (const auto& p : fluid) {
            if (!contains(p->pos, domainMin, domainMax)) continue;
            double rate = (currentActiveU - p->u) / tau;
            p->dudt += rate;
            power += p->mass * rate;
        }
        return power;
    }
};
