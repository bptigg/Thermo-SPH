#include "explosion.h"
#include "Planer.h"
#include "Corner.h"
#include <cmath>

ExplosionIC::ExplosionIC() : ExplosionIC(Parameters{}) {}
ExplosionIC::ExplosionIC(Parameters params) : params_(params) {}

void ExplosionIC::buildScenario() {
    if (built_) return;

    generatedParticles_.clear();
    generatedRigidObjects_.clear();
    blastEnergy_ = 0.0;
    int id = 0;

    const double h = params_.spacing;
    const double particleMass = params_.gasDensity * h * h;
    const int layers = params_.wallLayers;
    const Vector2D& dmin = params_.domainMin;
    const Vector2D& dmax = params_.domainMax;

    // Only layer-1 particles mirror fluid; the surface sits half a spacing outside the first gas row.
    auto makeWall = [&](const std::vector<Vector2D>& positions, const std::vector<bool>& isHost, const Vector2D& normal) {
        auto body = std::make_shared<RigidObject>();
        body->setType(RigidBodyType::EXTERNAL_BOUNDARY);
        for (size_t i = 0; i < positions.size(); ++i) {
            auto p = std::make_shared<SolidParticle>(
                id++, positions[i], Vector2D(0.0, 0.0), 2.0 * particleMass,
                params_.gasDensity, params_.ambientInternalEnergy, MotionType::STATIC);
            if (isHost[i]) p->setBoundaryComponent(std::make_shared<PlanarBoundary>(normal, true, h, 1.5 * h));
            generatedParticles_.push_back(p);
            body->addParticle(p);
        }
        body->finalizeInitialization();
        body->setConstraints(true, true, true);
        generatedRigidObjects_.push_back(body);
    };

    const int nx = static_cast<int>(std::floor((dmax.x - dmin.x) / h + 1e-9));
    const int ny = static_cast<int>(std::floor((dmax.y - dmin.y) / h + 1e-9));

    std::vector<Vector2D> floorPos, ceilPos, leftPos, rightPos;
    std::vector<bool> floorHost, ceilHost, leftHost, rightHost;
    for (int i = -layers; i <= nx + layers; ++i) {
        for (int layer = 1; layer <= layers; ++layer) {
            floorPos.emplace_back(dmin.x + i * h, dmin.y - layer * h);
            floorHost.push_back(layer == 1);
            ceilPos.emplace_back(dmin.x + i * h, dmax.y + layer * h);
            ceilHost.push_back(layer == 1);
        }
    }
    for (int j = 0; j <= ny; ++j) {
        for (int layer = 1; layer <= layers; ++layer) {
            leftPos.emplace_back(dmin.x - layer * h, dmin.y + j * h);
            leftHost.push_back(layer == 1);
            rightPos.emplace_back(dmax.x + layer * h, dmin.y + j * h);
            rightHost.push_back(layer == 1);
        }
    }
    makeWall(floorPos, floorHost, Vector2D(0.0, 1.0));
    makeWall(ceilPos, ceilHost, Vector2D(0.0, -1.0));
    makeWall(leftPos, leftHost, Vector2D(1.0, 0.0));
    makeWall(rightPos, rightHost, Vector2D(-1.0, 0.0));

    // Diagonal images at the four corners
    const double sx0 = dmin.x + 0.5 * h, sx1 = dmax.x - 0.5 * h;
    const double sy0 = dmin.y + 0.5 * h, sy1 = dmax.y - 0.5 * h;
    struct Corner { Vector2D host, nA, nB, point; };
    const Corner corners[] = {
        {Vector2D(dmin.x - 2 * h, dmin.y - 2 * h), Vector2D(1, 0), Vector2D(0, 1), Vector2D(sx0, sy0)},
        {Vector2D(dmax.x + 2 * h, dmin.y - 2 * h), Vector2D(-1, 0), Vector2D(0, 1), Vector2D(sx1, sy0)},
        {Vector2D(dmin.x - 2 * h, dmax.y + 2 * h), Vector2D(1, 0), Vector2D(0, -1), Vector2D(sx0, sy1)},
        {Vector2D(dmax.x + 2 * h, dmax.y + 2 * h), Vector2D(-1, 0), Vector2D(0, -1), Vector2D(sx1, sy1)},
    };
    for (const auto& c : corners) {
        auto p = std::make_shared<SolidParticle>(
            id++, c.host, Vector2D(0.0, 0.0), 2.0 * particleMass,
            params_.gasDensity, params_.ambientInternalEnergy, MotionType::STATIC);
        p->setBoundaryComponent(std::make_shared<CornerBoundary>(c.nA, c.nB, c.point));
        generatedParticles_.push_back(p);
    }

    // Gas, with the blast disc at the centre
    const double u0 = params_.ambientInternalEnergy;
    const double uBlast = u0 * params_.blastEnergyRatio;
    for (double x = dmin.x + h; x <= dmax.x - h + 1e-9; x += h) {
        for (double y = dmin.y + h; y <= dmax.y - h + 1e-9; y += h) {
            const double r = std::hypot(x - params_.blastCentre.x, y - params_.blastCentre.y);
            const double u = r <= params_.blastRadius ? uBlast : u0;
            blastEnergy_ += particleMass * (u - u0);
            generatedParticles_.push_back(std::make_shared<FluidParticle>(
                id++, Vector2D(x, y), Vector2D(0.0, 0.0), particleMass, params_.gasDensity, u));
        }
    }

    built_ = true;
}

std::vector<std::shared_ptr<Particle>> ExplosionIC::generateParticles() {
    buildScenario();
    return generatedParticles_;
}

std::vector<std::shared_ptr<RigidObject>> ExplosionIC::getRigidObjects() {
    buildScenario();
    return generatedRigidObjects_;
}

double ExplosionIC::getBlastEnergy() {
    buildScenario();
    return blastEnergy_;
}
