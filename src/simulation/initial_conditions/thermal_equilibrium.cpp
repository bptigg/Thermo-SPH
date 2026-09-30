#include "thermal_equilibrium.h"
#include "Planer.h"
#include <cmath>

ThermalEquilibriumIC::ThermalEquilibriumIC() : ThermalEquilibriumIC(Parameters{}) {}
ThermalEquilibriumIC::ThermalEquilibriumIC(Parameters params) : params_(params) {}

void ThermalEquilibriumIC::buildScenario() {
    if (built_) return;

    generatedParticles_.clear();
    generatedRigidObjects_.clear();
    int id = 0;

    const double h_cold = params_.spacing;
    const double rho_cold = params_.density;
    const double u_cold = params_.uCold;
    const double u_hot = params_.uHot;

    const double particleMass = rho_cold * h_cold * h_cold;

    const double rho_hot = rho_cold * (u_cold / u_hot);
    const double h_hot = h_cold * std::sqrt(u_hot / u_cold);
    const int layers = params_.wallLayers;

    auto registerStaticWallBody = [&](std::vector<std::shared_ptr<Particle>>& wallParticles, const Vector2D& normal) {
        if (wallParticles.empty()) return;

        auto wallBody = std::make_shared<RigidObject>();
        wallBody->setType(RigidBodyType::EXTERNAL_BOUNDARY);
        for (const auto& p : wallParticles) {
            p->setBoundaryComponent(std::make_shared<PlanarBoundary>(normal, true));
            wallBody->addParticle(p);
        }
        wallBody->finalizeInitialization();
        wallBody->setConstraints(true, true, true);
        generatedRigidObjects_.push_back(wallBody);
    };

    if (params_.bw) {
        std::vector<std::shared_ptr<Particle>> bottomWall;
        // Hot left side
        for (double x = params_.domainMin.x; x < 0.0; x += h_hot) {
            auto wallP = std::make_shared<SolidParticle>(
                id++, Vector2D(x, params_.domainMin.y), Vector2D(0.0, 0.0),
                particleMass, rho_hot, u_hot, MotionType::STATIC
            );
            generatedParticles_.push_back(wallP);
            bottomWall.push_back(wallP);
        }
        // Cold right side
        for (double x = 0.0; x <= params_.domainMax.x; x += h_cold) {
            auto wallP = std::make_shared<SolidParticle>(
                id++, Vector2D(x, params_.domainMin.y), Vector2D(0.0, 0.0),
                particleMass, rho_cold, u_cold, MotionType::STATIC
            );
            generatedParticles_.push_back(wallP);
            bottomWall.push_back(wallP);
        }
        registerStaticWallBody(bottomWall, Vector2D(0.0, 1.0));
    }
    if (params_.tw) {
        std::vector<std::shared_ptr<Particle>> topWall;
        // Hot left side
        for (double x = params_.domainMin.x; x < 0.0; x += h_hot) {
            auto wallP = std::make_shared<SolidParticle>(
                id++, Vector2D(x, params_.domainMax.y), Vector2D(0.0, 0.0),
                particleMass, rho_hot, u_hot, MotionType::STATIC
            );
            generatedParticles_.push_back(wallP);
            topWall.push_back(wallP);
        }
        // Cold right side
        for (double x = 0.0; x <= params_.domainMax.x; x += h_cold) {
            auto wallP = std::make_shared<SolidParticle>(
                id++, Vector2D(x, params_.domainMax.y), Vector2D(0.0, 0.0),
                particleMass, rho_cold, u_cold, MotionType::STATIC
            );
            generatedParticles_.push_back(wallP);
            topWall.push_back(wallP);
        }
        registerStaticWallBody(topWall, Vector2D(0.0, -1.0));
    }
    if (params_.lw) {
        std::vector<std::shared_ptr<Particle>> leftWall;
        for (double y = params_.domainMin.y; y <= params_.domainMax.y; y += h_hot) {
            auto wallP = std::make_shared<SolidParticle>(
                id++, Vector2D(params_.domainMin.x, y), Vector2D(0.0, 0.0),
                particleMass, rho_hot, u_hot, MotionType::STATIC
            );
            generatedParticles_.push_back(wallP);
            leftWall.push_back(wallP);
        }
        registerStaticWallBody(leftWall, Vector2D(1.0, 0.0));
    }
    if (params_.rw) {
        std::vector<std::shared_ptr<Particle>> rightWall;
        for (double y = params_.domainMin.y; y <= params_.domainMax.y; y += h_cold) {
            auto wallP = std::make_shared<SolidParticle>(
                id++, Vector2D(params_.domainMax.x, y), Vector2D(0.0, 0.0),
                particleMass, rho_cold, u_cold, MotionType::STATIC
            );
            generatedParticles_.push_back(wallP);
            rightWall.push_back(wallP);
        }
        registerStaticWallBody(rightWall, Vector2D(-1.0, 0.0));
    }

    for (double x = params_.domainMin.x + 0.5 * h_hot; x < 0.0; x += h_hot) {
        for (double y = params_.domainMin.y + 0.5 * h_hot; y <= params_.domainMax.y; y += h_hot) {
            auto p = std::make_shared<FluidParticle>(
                id++, Vector2D(x, y), Vector2D(0.0, 0.0),
                particleMass, rho_hot, u_hot
            );
            generatedParticles_.push_back(p);
        }
    }

    for (double x = 0.5 * h_cold; x <= params_.domainMax.x; x += h_cold) {
        for (double y = params_.domainMin.y + 0.5 * h_cold; y <= params_.domainMax.y; y += h_cold) {
            auto p = std::make_shared<FluidParticle>(
                id++, Vector2D(x, y), Vector2D(0.0, 0.0),
                particleMass, rho_cold, u_cold
            );
            generatedParticles_.push_back(p);
        }
    }
    built_ = true;
}

std::vector<std::shared_ptr<Particle>> ThermalEquilibriumIC::generateParticles() {
    buildScenario();
    return generatedParticles_;
}

std::vector<std::shared_ptr<RigidObject>> ThermalEquilibriumIC::getRigidObjects() {
    buildScenario();
    return generatedRigidObjects_;
}