#include "test_Piston.h"
#include "Piston.h"
#include "rigidbodysolver.h"
#include "piston_compression.h"
#include "thermalbath.h"
#include "SpatialGrid.h"
#include "ThreadPool.h"
#include "GhostManager.h"
#include <map>
#include <algorithm>
#include <memory>

static std::shared_ptr<Piston> makePiston(const DirectionLocks& locks) {
    auto piston = std::make_shared<Piston>(locks);
    piston->addParticle(std::make_shared<SolidParticle>(
        1, Vector2D(0.0, 0.0), Vector2D(0.0, 0.0), 1.0, 1.0, 1.0, MotionType::NONSTATIC));
    piston->addParticle(std::make_shared<SolidParticle>(
        2, Vector2D(1.0, 0.0), Vector2D(0.0, 0.0), 1.0, 1.0, 1.0, MotionType::NONSTATIC));
    piston->finalizeInitialization();
    return piston;
}

static bool testDirectionLocksFilterForce() {
    auto piston = makePiston(DirectionLocks::upwardOnly());
    piston->addForceAtPosition(Vector2D(0.0, -5.0), piston->getCenterOfMass());
    TEST_ASSERT_NEAR(piston->getForceAccumulator().y, 0.0, 1e-12, "Downward force should be ignored");

    piston->addForceAtPosition(Vector2D(0.0, 3.0), piston->getCenterOfMass());
    TEST_ASSERT_NEAR(piston->getForceAccumulator().y, 3.0, 1e-12, "Upward force should be accepted");

    piston->addForceAtPosition(Vector2D(2.0, 0.0), piston->getCenterOfMass());
    TEST_ASSERT_NEAR(piston->getForceAccumulator().x, 0.0, 1e-12, "Sideways force should be ignored");
    return true;
}

static bool testAllLocksIgnoreEveryForce() {
    auto piston = makePiston(DirectionLocks::all());
    TEST_ASSERT(piston->isKinematic(), "All locks should make the piston kinematic");
    piston->addForceAtPosition(Vector2D(4.0, 5.0), piston->getCenterOfMass());
    piston->addForceAtPosition(Vector2D(-4.0, -5.0), piston->getCenterOfMass());
    TEST_ASSERT_NEAR(piston->getForceAccumulator().normSq(), 0.0, 1e-12, "Kinematic piston must ignore all forces");
    return true;
}

static bool testKinematicPistonKeepsPrescribedVelocity() {
    auto piston = Piston::makeKinematic(Vector2D(0.0, -0.5));
    piston->addParticle(std::make_shared<SolidParticle>(
        1, Vector2D(0.0, 1.0), Vector2D(0.0, 0.0), 1.0, 1.0, 1.0, MotionType::NONSTATIC));
    piston->finalizeInitialization();

    RigidBodySolver solver;
    std::vector<std::shared_ptr<RigidObject>> bodies{piston};
    for (int i = 0; i < 10; ++i) {
        piston->clearForces();
        piston->addForceAtPosition(Vector2D(0.0, 100.0), piston->getCenterOfMass());
        piston->preStep(0.1 * i, 0.1);
        solver.integrate(bodies, 0.1);
    }

    TEST_ASSERT_NEAR(piston->getLinearVel().y, -0.5, 1e-9, "Kinematic velocity should be unchanged by forces");
    TEST_ASSERT_NEAR(piston->getCenterOfMass().y, 0.5, 1e-9, "Kinematic piston should travel at the set speed");
    TEST_ASSERT_NEAR(piston->getParticles()[0]->pos.y, 0.5, 1e-9, "Particles should follow the piston");
    TEST_ASSERT_NEAR(piston->getParticles()[0]->vel.y, -0.5, 1e-9, "Particles should carry the piston velocity");
    return true;
}

static bool testDynamicPistonOnlyRises() {
    auto piston = makePiston(DirectionLocks::upwardOnly());
    RigidBodySolver solver;
    std::vector<std::shared_ptr<RigidObject>> bodies{piston};

    piston->clearForces();
    piston->addForceAtPosition(Vector2D(0.0, -10.0), piston->getCenterOfMass());
    solver.integrate(bodies, 0.1);
    TEST_ASSERT_NEAR(piston->getLinearVel().y, 0.0, 1e-12, "Downward push must not move the piston");

    piston->clearForces();
    piston->addForceAtPosition(Vector2D(0.0, 10.0), piston->getCenterOfMass());
    solver.integrate(bodies, 0.1);
    TEST_ASSERT_NEAR(piston->getLinearVel().y, 0.5, 1e-9, "Upward push should accelerate the piston (F/m*dt)");
    TEST_ASSERT(piston->getCenterOfMass().y > 0.0, "Piston should have moved up");
    TEST_ASSERT_NEAR(piston->getLinearVel().x, 0.0, 1e-12, "Piston should not drift sideways");
    return true;
}

static bool testScheduledRelease() {
    auto piston = makePiston(DirectionLocks::all());
    piston->setLinearVel(Vector2D(0.0, -1.0));
    piston->scheduleRelease(1.0, DirectionLocks::upwardOnly());

    piston->preStep(0.5, 0.1);
    TEST_ASSERT(piston->isKinematic() && !piston->isReleased(), "Piston should stay kinematic before release time");
    TEST_ASSERT_NEAR(piston->getLinearVel().y, -1.0, 1e-12, "Velocity should be kept before release");

    piston->preStep(1.0, 0.1);
    TEST_ASSERT(!piston->isKinematic() && piston->isReleased(), "Piston should become dynamic at release time");
    TEST_ASSERT_NEAR(piston->getLinearVel().y, 0.0, 1e-12, "Velocity should be zeroed on release");
    return true;
}

static bool testConstantLoadOnlyAffectsDynamicPiston() {
    auto dynamic = makePiston(DirectionLocks::upwardOnly());
    dynamic->setConstantLoad(Vector2D(0.0, -2.0));
    dynamic->preStep(0.0, 0.1);
    TEST_ASSERT_NEAR(dynamic->getForceAccumulator().y, -2.0, 1e-12, "Load should bypass the direction locks");

    auto kinematic = makePiston(DirectionLocks::all());
    kinematic->setConstantLoad(Vector2D(0.0, -2.0));
    kinematic->preStep(0.0, 0.1);
    TEST_ASSERT_NEAR(kinematic->getForceAccumulator().y, 0.0, 1e-12, "Kinematic piston should ignore the load");
    return true;
}

static bool testTravelLimits() {
    auto piston = makePiston(DirectionLocks::all());
    piston->setTravelLimits(-1.0, 1.0);
    piston->setLinearVel(Vector2D(0.0, 1.0));
    piston->setCenterOfMass(Vector2D(0.5, 1.2));
    piston->preStep(0.0, 0.1);

    TEST_ASSERT_NEAR(piston->getCenterOfMass().y, 1.0, 1e-12, "Piston should be clamped to the upper limit");
    TEST_ASSERT_NEAR(piston->getLinearVel().y, 0.0, 1e-12, "Piston should stop at the upper limit");
    TEST_ASSERT_NEAR(piston->getParticles()[0]->pos.y, 1.0, 1e-12, "Particles should follow the clamped position");
    return true;
}

static bool testPistonSkippedByCollisionSolver() {
    auto piston = makePiston(DirectionLocks::all());
    auto wall = std::make_shared<RigidObject>();
    wall->addParticle(std::make_shared<SolidParticle>(
        3, Vector2D(0.0, 0.02), Vector2D(0.0, 0.0), 1.0, 1.0, 1.0, MotionType::STATIC));
    wall->finalizeInitialization();
    wall->setConstraints(true, true, true);

    RigidBodySolver solver;
    std::vector<std::shared_ptr<RigidObject>> bodies{piston, wall};
    TEST_ASSERT(!solver.solveCollisions(bodies, 0.01), "Piston should not generate rigid-rigid contacts");
    return true;
}

static bool testPistonCompressionScenario() {
    PistonCompressionIC::Parameters params;
    params.compressionTime = 0.5;
    params.compressionSpeed = 0.1;
    PistonCompressionIC ic(params);
    auto particles = ic.generateParticles();
    auto bodies = ic.getRigidObjects();
    auto piston = ic.getPiston();

    TEST_ASSERT(piston != nullptr, "Scenario should provide a piston");
    TEST_ASSERT(piston->isKinematic(), "Piston should start kinematic");
    TEST_ASSERT(piston->getType() == RigidBodyType::INTERNAL_OBJECT, "Piston should interact via rigid-body particles");
    TEST_ASSERT(piston->getLinearVel().y < 0.0, "Piston should start moving down");

    bool inBodies = false;
    for (const auto& b : bodies) if (b == piston) inBodies = true;
    TEST_ASSERT(inBodies, "Piston should be registered as a rigid object");

    double maxGasY = 0.0;
    size_t gasCount = 0;
    for (const auto& p : particles) {
        TEST_ASSERT(p->parentBody == nullptr || p->parentBody != piston.get(), "Piston particles must not be in the particle list");
        if (p->isFluid()) {
            maxGasY = std::max(maxGasY, p->pos.y);
            ++gasCount;
        }
    }
    TEST_ASSERT(gasCount > 0, "Scenario should contain gas");
    double pistonBottom = piston->getAABB().min.y;
    TEST_ASSERT(pistonBottom > maxGasY, "Piston should start above the gas");

    // Compress with no gas feedback, then release: gas push up moves it, push down does not
    RigidBodySolver solver;
    std::vector<std::shared_ptr<RigidObject>> pistonOnly{piston};
    double y0 = piston->getCenterOfMass().y;
    for (double t = 0.0; t < 0.5; t += 0.01) {
        piston->clearForces();
        piston->preStep(t, 0.01);
        solver.integrate(pistonOnly, 0.01);
    }
    TEST_ASSERT(piston->getCenterOfMass().y < y0 - 0.04, "Piston should have compressed the gas");

    piston->clearForces();
    piston->preStep(0.5, 0.01);
    TEST_ASSERT(piston->isReleased() && !piston->isKinematic(), "Piston should be dynamic after compression");
    TEST_ASSERT_NEAR(piston->getLinearVel().y, 0.0, 1e-12, "Piston should be at rest after release");
    piston->addForceAtPosition(Vector2D(0.0, 5.0), piston->getCenterOfMass());
    solver.integrate(pistonOnly, 0.01);
    TEST_ASSERT(piston->getLinearVel().y > 0.0, "Gas pressure should push released piston up");
    return true;
}

static bool testWallGhostsAreNotDuplicated() {
    PistonCompressionIC ic;
    auto all = ic.generateParticles();
    std::vector<std::shared_ptr<Particle>> fluid, walls;
    for (const auto& p : all) {
        if (p->isFluid()) fluid.push_back(p);
        else if (p->hasBoundary()) walls.push_back(p);
    }

    SpatialGrid grid(0.2);
    grid.build(fluid);
    ThreadPool pool(2);
    pool.start();
    auto ghosts = GhostManager::generateAllGhosts(walls, fluid, 0.2, grid, pool);
    pool.Stop();

    std::map<std::pair<long, long>, int> multiplicity;
    for (const auto& g : ghosts) {
        multiplicity[{std::lround(g->pos.x * 1e6), std::lround(g->pos.y * 1e6)}]++;
    }
    int maxMultiplicity = 0;
    for (const auto& [pos, count] : multiplicity) maxMultiplicity = std::max(maxMultiplicity, count);

    TEST_ASSERT(!ghosts.empty(), "Walls should generate ghosts");
    // Corner particles may be mirrored once per adjacent wall, so allow at most one per wall plane.
    TEST_ASSERT(maxMultiplicity <= 2, "Each fluid particle must be mirrored once per wall plane");
    TEST_ASSERT(ghosts.size() < 4 * fluid.size(), "Ghost count should be of the order of the near-wall fluid count");
    return true;
}

static bool testLidSealStopsPenetration() {
    auto piston = makePiston(DirectionLocks::all());
    piston->setLinearVel(Vector2D(0.0, -0.5));
    piston->setLidSeal(0.1);

    auto below = std::make_shared<FluidParticle>(10, Vector2D(0.5, -0.3), Vector2D(0.0, 2.0), 1.0, 1.0, 1.0);
    auto inside = std::make_shared<FluidParticle>(11, Vector2D(0.5, -0.05), Vector2D(0.0, 2.0), 1.0, 1.0, 1.0);
    auto far = std::make_shared<FluidParticle>(12, Vector2D(0.5, -1.0), Vector2D(0.0, 2.0), 1.0, 1.0, 1.0);
    piston->constrainFluid({below, inside, far});

    TEST_ASSERT_NEAR(below->pos.y, -0.3, 1e-12, "Fluid well under the face should not move");
    TEST_ASSERT_NEAR(inside->pos.y, -0.1, 1e-12, "Fluid inside the lid should be pushed to the face");
    TEST_ASSERT_NEAR(inside->vel.y, -0.5, 1e-12, "Pushed fluid should take the piston velocity");
    TEST_ASSERT_NEAR(far->vel.y, 2.0, 1e-12, "Distant fluid should be untouched");
    return true;
}

static bool testConstantForceScenario() {
    PistonCompressionIC::Parameters params;
    params.driveMode = PistonCompressionIC::DriveMode::ConstantForce;
    PistonCompressionIC ic(params);
    auto piston = ic.getPiston();
    const double startY = piston->getCenterOfMass().y;

    TEST_ASSERT(!piston->isKinematic(), "Constant-force piston should respond to gas forces");
    TEST_ASSERT_NEAR(piston->getLinearVel().y, 0.0, 1e-12, "Piston should start at rest");
    TEST_ASSERT(piston->getConstantLoad().y < 0.0, "A downward force should act on the piston");

    RigidBodySolver solver;
    std::vector<std::shared_ptr<RigidObject>> bodies{piston};
    for (int i = 0; i < 20; ++i) {
        piston->clearForces();
        piston->preStep(0.01 * i, 0.01);
        piston->addForceAtPosition(Vector2D(0.0, 100.0), piston->getCenterOfMass());
        solver.integrate(bodies, 0.01);
    }
    piston->preStep(0.2, 0.01);
    TEST_ASSERT(piston->getCenterOfMass().y <= startY + 1e-9, "Piston must not rise above its starting height");

    for (int i = 0; i < 400; ++i) {
        piston->clearForces();
        piston->preStep(0.2 + 0.01 * i, 0.01);
        solver.integrate(bodies, 0.01);
    }
    piston->preStep(10.0, 0.01);
    TEST_ASSERT(piston->getCenterOfMass().y < startY - 0.1, "Load should push the piston well below its start");
    TEST_ASSERT(piston->getCenterOfMass().y >= params.pistonMinY - 1e-9, "Piston must stop at the lower limit");
    return true;
}

static bool testThermalBathRelaxesNearWalls() {
    ThermalBath bath;
    bath.enabled = true;
    bath.targetU = 0.5;
    bath.tau = 0.1;
    bath.wallThickness = 0.1;

    const Vector2D dmin(0.0, 0.0), dmax(1.0, 1.0);
    auto hot = std::make_shared<FluidParticle>(1, Vector2D(0.05, 0.5), Vector2D(0.0, 0.0), 2.0, 1.0, 1.0);
    auto cold = std::make_shared<FluidParticle>(2, Vector2D(0.5, 0.05), Vector2D(0.0, 0.0), 1.0, 1.0, 0.25);
    auto middle = std::make_shared<FluidParticle>(3, Vector2D(0.5, 0.5), Vector2D(0.0, 0.0), 1.0, 1.0, 2.0);

    double power = bath.apply({hot, cold, middle}, dmin, dmax);
    TEST_ASSERT_NEAR(hot->dudt, -5.0, 1e-12, "Hot gas at a wall should cool toward the bath");
    TEST_ASSERT_NEAR(cold->dudt, 2.5, 1e-12, "Cold gas at a wall should warm toward the bath");
    TEST_ASSERT_NEAR(middle->dudt, 0.0, 1e-12, "Gas away from the walls should be unaffected");
    TEST_ASSERT_NEAR(power, 2.0 * -5.0 + 1.0 * 2.5, 1e-12, "Returned power should be the mass-weighted heat flow");

    bath.wallThickness = 0.0;
    middle->dudt = 0.0;
    bath.apply({middle}, dmin, dmax);
    TEST_ASSERT(middle->dudt < 0.0, "Zero thickness should couple the whole gas");

    bath.enabled = false;
    middle->dudt = 0.0;
    TEST_ASSERT_NEAR(bath.apply({middle}, dmin, dmax), 0.0, 1e-12, "Disabled bath should do nothing");
    TEST_ASSERT_NEAR(middle->dudt, 0.0, 1e-12, "Disabled bath must not touch dudt");
    return true;
}

static bool testPowerTakeOffExtractsWork() {
    auto piston = makePiston(DirectionLocks::upwardOnly());
    piston->setPowerTakeOff(2.0);
    piston->setLinearVel(Vector2D(0.0, 3.0));

    piston->preStep(0.0, 0.1);
    TEST_ASSERT_NEAR(piston->getForceAccumulator().y, -6.0, 1e-12, "Take-off should oppose the motion with -c*v");
    TEST_ASSERT_NEAR(piston->getHarvestedEnergy(), 2.0 * 9.0 * 0.1, 1e-12, "Harvested energy should be c*v^2*dt");

    piston->clearForces();
    piston->setLinearVel(Vector2D(0.0, -3.0));
    piston->preStep(0.1, 0.1);
    TEST_ASSERT_NEAR(piston->getForceAccumulator().y, 0.0, 1e-12, "A falling piston (compression stroke) is not braked");

    auto kinematic = makePiston(DirectionLocks::all());
    kinematic->setPowerTakeOff(2.0);
    kinematic->setLinearVel(Vector2D(0.0, 3.0));
    kinematic->preStep(0.0, 0.1);
    TEST_ASSERT_NEAR(kinematic->getHarvestedEnergy(), 0.0, 1e-12, "A kinematic piston must not be braked");
    return true;
}

static bool testFlywheelCrankKinematics() {
    auto piston = makePiston(DirectionLocks::upwardOnly());
    const double startY = piston->getCenterOfMass().y;
    piston->enableFlywheel(0.5, 1.0, 2.0);
    TEST_ASSERT(piston->isKinematic(), "A flywheel piston follows its crank, not forces");

    RigidBodySolver solver;
    std::vector<std::shared_ptr<RigidObject>> bodies{piston};
    const double dt = 1e-3;
    for (int i = 0; i < 1000; ++i) {
        piston->clearForces();
        piston->preStep(i * dt, dt);
        solver.integrate(bodies, dt);
    }

    const double expectedY = startY - 0.5 + 0.5 * std::cos(piston->getCrankAngle());
    TEST_ASSERT_NEAR(piston->getFlywheelOmega(), 2.0, 1e-9, "A free flywheel keeps its speed");
    TEST_ASSERT_NEAR(piston->getCrankAngle(), 2.0, 1e-6, "theta should advance by omega*t");
    TEST_ASSERT_NEAR(piston->getCenterOfMass().y, expectedY, 2e-3, "Piston height should follow the crank");
    return true;
}

static bool testFlywheelGasTorqueAndShaftHarvest() {
    auto piston = makePiston(DirectionLocks::upwardOnly());
    piston->enableFlywheel(0.5, 1.0, 1.0);
    const double dt = 1e-3;

    // Piston moving down (omega > 0, theta > 0); an upward gas push resists it and slows the flywheel
    piston->preStep(0.0, dt);
    for (int i = 0; i < 200; ++i) piston->preStep(i * dt, dt);
    const double omegaBefore = piston->getFlywheelOmega();
    piston->addForceAtPosition(Vector2D(0.0, 10.0), piston->getCenterOfMass());
    piston->preStep(1.0, dt);
    TEST_ASSERT(piston->getFlywheelOmega() < omegaBefore, "Gas pushing against the stroke should slow the flywheel");

    auto generator = makePiston(DirectionLocks::upwardOnly());
    generator->enableFlywheel(0.5, 1.0, 2.0, 0.4);
    generator->preStep(0.0, 0.1);
    TEST_ASSERT_NEAR(generator->getHarvestedEnergy(), 0.4 * 4.0 * 0.1, 1e-12, "Shaft harvest should be c*omega^2*dt");
    TEST_ASSERT_NEAR(generator->getFlywheelOmega(), 2.0 - 0.4 * 2.0 / 1.0 * 0.1, 1e-12, "Shaft drag should slow the flywheel");
    return true;
}

void registerPistonTests(TestSuite& suite) {
    suite.startSection("Piston Module");
    suite.runTest("Flywheel Crank Kinematics", testFlywheelCrankKinematics);
    suite.runTest("Flywheel Gas Torque and Shaft", testFlywheelGasTorqueAndShaftHarvest);
    suite.runTest("Power Take-Off", testPowerTakeOffExtractsWork);
    suite.runTest("Thermal Bath Relaxation", testThermalBathRelaxesNearWalls);
    suite.runTest("Constant Force Scenario", testConstantForceScenario);
    suite.runTest("Lid Seal Stops Penetration", testLidSealStopsPenetration);
    suite.runTest("Wall Ghosts Not Duplicated", testWallGhostsAreNotDuplicated);
    suite.runTest("Direction Locks Filter Force", testDirectionLocksFilterForce);
    suite.runTest("All Locks Ignore Every Force", testAllLocksIgnoreEveryForce);
    suite.runTest("Kinematic Piston Keeps Velocity", testKinematicPistonKeepsPrescribedVelocity);
    suite.runTest("Dynamic Piston Only Rises", testDynamicPistonOnlyRises);
    suite.runTest("Scheduled Release", testScheduledRelease);
    suite.runTest("Constant Load", testConstantLoadOnlyAffectsDynamicPiston);
    suite.runTest("Travel Limits", testTravelLimits);
    suite.runTest("Piston Skipped By Collision Solver", testPistonSkippedByCollisionSolver);
    suite.runTest("Piston Compression Scenario", testPistonCompressionScenario);
}
