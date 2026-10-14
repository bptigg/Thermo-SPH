#include "piston_compression.h"
#include "Planer.h"
#include "Corner.h"
#include <cmath>

PistonCompressionIC::PistonCompressionIC() : PistonCompressionIC(Parameters{}) {}
PistonCompressionIC::PistonCompressionIC(Parameters params) : params_(params) {}

void PistonCompressionIC::buildScenario() {
    if (built_) return;

    generatedParticles_.clear();
    generatedRigidObjects_.clear();
    int id = 0;

    const double h = params_.spacing;
    const double particleMass = params_.gasDensity * h * h;
    const int layers = params_.wallLayers;

    // Only layer-1 particles mirror fluid; the physical surface sits half a spacing outside the first gas row.
    auto registerStaticWallBody = [&](std::vector<std::shared_ptr<Particle>>& wallParticles,
                                      const std::vector<std::shared_ptr<Particle>>& hosts, const Vector2D& normal) {
        auto wallBody = std::make_shared<RigidObject>();
        wallBody->setType(RigidBodyType::EXTERNAL_BOUNDARY);
        for (const auto& p : hosts) {
            p->setBoundaryComponent(std::make_shared<PlanarBoundary>(normal, true, h, 1.5 * h));
        }
        for (const auto& p : wallParticles) {
            wallBody->addParticle(p);
        }
        wallBody->finalizeInitialization();
        wallBody->setConstraints(true, true, true);
        generatedRigidObjects_.push_back(wallBody);
    };

    auto makeWallParticle = [&](const Vector2D& pos, std::vector<std::shared_ptr<Particle>>& group,
                                std::vector<std::shared_ptr<Particle>>* hosts = nullptr) {
        auto p = std::make_shared<SolidParticle>(
            id++, pos, Vector2D(0.0, 0.0), 2.0 * particleMass,
            params_.gasDensity, params_.gasInternalEnergy, MotionType::STATIC);
        generatedParticles_.push_back(p);
        group.push_back(p);
        if (hosts) hosts->push_back(p);
    };

    const int nx = static_cast<int>(std::floor((params_.domainMax.x - params_.domainMin.x) / h + 1e-9));
    const int ny = static_cast<int>(std::floor((params_.domainMax.y - params_.domainMin.y) / h + 1e-9));

    // Floor
    std::vector<std::shared_ptr<Particle>> floorParticles, floorHosts;
    for (int i = -layers; i <= nx + layers; ++i) {
        for (int layer = 1; layer <= layers; ++layer) {
            makeWallParticle(Vector2D(params_.domainMin.x + i * h, params_.domainMin.y - layer * h),
                             floorParticles, layer == 1 ? &floorHosts : nullptr);
        }
    }
    registerStaticWallBody(floorParticles, floorHosts, Vector2D(0.0, 1.0));

    // Left and right walls
    std::vector<std::shared_ptr<Particle>> leftParticles, rightParticles, leftHosts, rightHosts;
    for (int j = 0; j <= ny; ++j) {
        for (int layer = 1; layer <= layers; ++layer) {
            double y = params_.domainMin.y + j * h;
            makeWallParticle(Vector2D(params_.domainMin.x - layer * h, y), leftParticles,
                             layer == 1 ? &leftHosts : nullptr);
            makeWallParticle(Vector2D(params_.domainMax.x + layer * h, y), rightParticles,
                             layer == 1 ? &rightHosts : nullptr);
        }
    }
    registerStaticWallBody(leftParticles, leftHosts, Vector2D(1.0, 0.0));
    registerStaticWallBody(rightParticles, rightHosts, Vector2D(-1.0, 0.0));

    // Diagonal images at the two floor corners
    {
        const double surfMinX = params_.domainMin.x + 0.5 * h;
        const double surfMaxX = params_.domainMax.x - 0.5 * h;
        const double surfMinY = params_.domainMin.y + 0.5 * h;
        std::vector<std::shared_ptr<Particle>> cornerHosts;
        makeWallParticle(Vector2D(params_.domainMin.x - 2.0 * h, params_.domainMin.y - 2.0 * h), cornerHosts);
        makeWallParticle(Vector2D(params_.domainMax.x + 2.0 * h, params_.domainMin.y - 2.0 * h), cornerHosts);
        cornerHosts[0]->setBoundaryComponent(std::make_shared<CornerBoundary>(
            Vector2D(1.0, 0.0), Vector2D(0.0, 1.0), Vector2D(surfMinX, surfMinY)));
        cornerHosts[1]->setBoundaryComponent(std::make_shared<CornerBoundary>(
            Vector2D(-1.0, 0.0), Vector2D(0.0, 1.0), Vector2D(surfMaxX, surfMinY)));
    }

    // Gas
    const double xMin = params_.domainMin.x + h;
    const double xMax = params_.domainMax.x - h;
    const double yMin = params_.domainMin.y + h;
    const double yMax = params_.domainMin.y + params_.gasHeight;
    for (double x = xMin; x <= xMax + 1e-9; x += h) {
        for (double y = yMin; y <= yMax + 1e-9; y += h) {
            generatedParticles_.push_back(std::make_shared<FluidParticle>(
                id++, Vector2D(x, y), Vector2D(0.0, 0.0),
                particleMass, params_.gasDensity, params_.gasInternalEnergy));
        }
    }

    // Piston sits on top of the gas
    piston_ = std::make_shared<Piston>(DirectionLocks::all());
    const double pistonMass = params_.pistonDensityRatio * particleMass;
    const double pistonBottom = yMax + 1.5 * h;
    // Edge lid particles sit on the wall surfaces so gas cannot squeeze up the side gap
    const double pistonXMin = params_.domainMin.x + 0.5 * h;
    const double pistonXMax = params_.domainMax.x - 0.5 * h;
    // Extra columns reach into the walls so the gas sees a complete lid at the top corners
    const int wallExtension = 3;
    for (double x = pistonXMin - wallExtension * h; x <= pistonXMax + wallExtension * h + 1e-9; x += h) {
        for (int layer = 0; layer < params_.pistonLayers; ++layer) {
            piston_->addParticle(std::make_shared<SolidParticle>(
                id++, Vector2D(x, pistonBottom + layer * h), Vector2D(0.0, 0.0),
                pistonMass, params_.gasDensity, params_.gasInternalEnergy, MotionType::NONSTATIC));
        }
    }
    piston_->finalizeInitialization();

    // Initial gas pressure P = (gamma - 1) rho u acting over the piston width
    const double width = pistonXMax - pistonXMin + h;
    const double gasForce = (params_.gamma - 1.0) * params_.gasDensity * params_.gasInternalEnergy * width;
    if (params_.driveMode == DriveMode::Flywheel) {
        piston_->enableFlywheel(params_.crankRadius, params_.flywheelInertia,
                                params_.flywheelInitialOmega, params_.shaftDamping);
    } else if (params_.driveMode == DriveMode::ConstantForce) {
        piston_->setLocks(DirectionLocks::upwardOnly());
        piston_->setConstantLoad(Vector2D(0.0, -params_.constantForceFactor * gasForce));
        // Starting height (plus any allowance) is the ceiling
        piston_->setTravelLimits(params_.pistonMinY,
                                 piston_->getCenterOfMass().y + params_.pistonCeilingAboveStart);
    } else {
        piston_->setConstantLoad(Vector2D(0.0, -params_.loadFactor * gasForce));
        piston_->setTravelLimits(params_.pistonMinY, params_.domainMax.y - 0.1);
        piston_->setLinearVel(Vector2D(0.0, -params_.compressionSpeed));
        piston_->scheduleRelease(params_.compressionTime, DirectionLocks::upwardOnly());
    }
    piston_->setLidSeal(0.5 * h);
    piston_->setPowerTakeOff(params_.powerTakeOff);
    generatedRigidObjects_.push_back(piston_);

    built_ = true;
}

std::vector<std::shared_ptr<Particle>> PistonCompressionIC::generateParticles() {
    buildScenario();
    return generatedParticles_;
}

std::vector<std::shared_ptr<RigidObject>> PistonCompressionIC::getRigidObjects() {
    buildScenario();
    return generatedRigidObjects_;
}

std::shared_ptr<Piston> PistonCompressionIC::getPiston() {
    buildScenario();
    return piston_;
}
