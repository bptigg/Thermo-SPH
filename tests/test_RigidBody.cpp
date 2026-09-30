#include "test_RigidBody.h"
#include "RigidBody.h"
#include "Particle.h"
#include "dam_break.h"
#include "rigidbodysolver.h"
#include "ThreadPool.h"
#include <memory>

static std::shared_ptr<SolidParticle> createParticle(int id, Vector2D pos, double mass = 1.0) {
    return std::make_shared<SolidParticle>(
        id, pos, Vector2D(0.0, 0.0), mass, 1000.0, 300.0, MotionType::STATIC
    );
}

static bool testFinalizeInitialization() {
    RigidObject obj;
    auto p1 = createParticle(1, Vector2D(0.0, 0.0), 2.0);
    auto p2 = createParticle(2, Vector2D(3.0, 0.0), 2.0);

    obj.addParticle(p1);
    obj.addParticle(p2);
    obj.finalizeInitialization();

    TEST_ASSERT_NEAR(obj.getTotalMass(), 4.0, 1e-6, "Total mass calculation failed");
    TEST_ASSERT_NEAR(obj.getCenterOfMass().x, 1.5, 1e-6, "Center of Mass X failed");
    TEST_ASSERT_NEAR(obj.getCenterOfMass().y, 0.0, 1e-6, "Center of Mass Y failed");
    TEST_ASSERT_NEAR(obj.getInertia(), 9.0, 1e-6, "Moment of Inertia calculation failed");

    return true;
}

static bool testForceAndImpulseState() {
    RigidObject obj;
    auto p1 = createParticle(1, Vector2D(-1.0, 0.0), 1.0);
    auto p2 = createParticle(2, Vector2D(1.0, 0.0), 1.0);

    obj.addParticle(p1);
    obj.addParticle(p2);
    obj.finalizeInitialization();
    obj.setConstraints(false, false, false);

    obj.addForceAtPosition(Vector2D(0.0, 4.0), Vector2D(0.0, 0.0));
    TEST_ASSERT_NEAR(obj.getForceAccumulator().y, 4.0, 1e-6, "Force accumulator failed");

    obj.applyImpulse(Vector2D(0.0, 2.0), Vector2D(1.0, 0.0));
    TEST_ASSERT_NEAR(obj.getLinearVel().y, 1.0, 1e-6, "Linear impulse resolution failed");
    TEST_ASSERT_NEAR(obj.getAngularVel(), 1.0, 1e-6, "Angular impulse resolution failed");

    return true;
}

static bool testRigidBodyRotatesAroundPivot() {
    auto obj = std::make_shared<RigidObject>();
    obj->addParticle(createParticle(1, Vector2D(1.0, 0.0)));
    obj->addParticle(createParticle(2, Vector2D(2.0, 0.0)));
    obj->finalizeInitialization();
    obj->setConstraints(false, false, false);
    obj->setPivot(Vector2D(0.0, 0.0));
    obj->setAngularVel(1.0);

    RigidBodySolver solver;
    std::vector<std::shared_ptr<RigidObject>> bodies{obj};
    solver.integratePositions(bodies, std::acos(-1.0) / 2.0);

    TEST_ASSERT_NEAR(obj->getPivotPosition().x, 0.0, 1e-6, "Pivot X moved");
    TEST_ASSERT_NEAR(obj->getPivotPosition().y, 0.0, 1e-6, "Pivot Y moved");
    TEST_ASSERT_NEAR(obj->getCenterOfMass().x, 0.0, 1e-6, "COM should orbit the pivot");
    TEST_ASSERT_NEAR(obj->getCenterOfMass().y, 1.5, 1e-6, "COM should orbit the pivot");
    TEST_ASSERT_NEAR(obj->getParticles()[0]->pos.x, 0.0, 1e-6, "First particle X rotation failed");
    TEST_ASSERT_NEAR(obj->getParticles()[0]->pos.y, 1.0, 1e-6, "First particle Y rotation failed");

    return true;
}

static bool testPivotForceUsesPivotInertia() {
    auto obj = std::make_shared<RigidObject>();
    obj->addParticle(createParticle(1, Vector2D(1.0, 0.0)));
    obj->addParticle(createParticle(2, Vector2D(2.0, 0.0)));
    obj->finalizeInitialization();
    obj->setConstraints(false, false, false);
    obj->setPivot(Vector2D(0.0, 0.0));
    obj->addForceAtPosition(Vector2D(0.0, 2.0), obj->getCenterOfMass());

    RigidBodySolver solver;
    std::vector<std::shared_ptr<RigidObject>> bodies{obj};
    solver.integrateVelocities(bodies, 1.0);

    TEST_ASSERT_NEAR(obj->getInertiaAboutPivot(), 5.0, 1e-6, "Parallel-axis inertia calculation failed");
    TEST_ASSERT_NEAR(obj->getAngularVel(), 0.6, 1e-6, "Pivot torque integration failed");
    TEST_ASSERT_NEAR(obj->getLinearVel().x, 0.0, 1e-6, "Pinned COM velocity X failed");
    TEST_ASSERT_NEAR(obj->getLinearVel().y, 0.9, 1e-6, "Pinned COM velocity Y failed");

    return true;
}

static bool testPivotRestoringTorque() {
    auto obj = std::make_shared<RigidObject>();
    obj->addParticle(createParticle(1, Vector2D(1.0, 0.0)));
    obj->addParticle(createParticle(2, Vector2D(2.0, 0.0)));
    obj->finalizeInitialization();
    obj->setPivot(Vector2D(0.0, 0.0));
    obj->setAngularSpring(10.0, 2.0, 0.0);
    obj->setAngle(0.2);

    TEST_ASSERT_NEAR(obj->getAngularSpringTorque(), -2.0, 1e-6,
                     "Restoring spring should apply torque toward its rest angle");

    RigidBodySolver solver;
    std::vector<std::shared_ptr<RigidObject>> bodies{obj};
    solver.integrateVelocities(bodies, 0.1);
    TEST_ASSERT(obj->getAngularVel() < 0.0,
                "Restoring spring should accelerate the displaced body toward rest");
    return true;
}

static bool testPivotCollisionProducesTorqueAboutAnchor() {
    RigidBodySolver solver(RigidBodySolverParameters{
        .enableGravity = false,
        .gravity = Vector2D(0.0, -9.81),
        .restitution = 0.1,
        .particleRadius = 0.05
    });

    auto pivoted = std::make_shared<RigidObject>();
    pivoted->addParticle(createParticle(1, Vector2D(1.0, 0.0)));
    pivoted->addParticle(createParticle(2, Vector2D(2.0, 0.1)));
    pivoted->finalizeInitialization();
    pivoted->setPivot(Vector2D(0.0, 0.0));
    pivoted->setAngularVel(1.0);
    pivoted->updateParticlePositions();

    auto obstacle = std::make_shared<RigidObject>();
    obstacle->addParticle(createParticle(3, Vector2D(1.0, 0.08)));
    obstacle->addParticle(createParticle(4, Vector2D(1.0, 0.0)));
    obstacle->finalizeInitialization();
    obstacle->setConstraints(true, true, true);

    std::vector<std::shared_ptr<RigidObject>> bodies{pivoted, obstacle};
    TEST_ASSERT(solver.solveCollisions(bodies, 0.01), "Pinned body contact should be detected");
    solver.integrateVelocities(bodies, 0.01);

    TEST_ASSERT(pivoted->getAngularVel() < 1.0, "Contact should change angular velocity about the pivot");
    TEST_ASSERT_NEAR(pivoted->getPivotPosition().x, 0.0, 1e-6, "Collision moved pivot X");
    TEST_ASSERT_NEAR(pivoted->getPivotPosition().y, 0.0, 1e-6, "Collision moved pivot Y");
    return true;
}

static bool testStaticBoundaryClassification() {
    RigidObject wall;
    auto p1 = createParticle(1, Vector2D(0.0, 0.0), 1.0);
    wall.addParticle(p1);
    wall.finalizeInitialization();
    wall.setConstraints(true, true, true);

    TEST_ASSERT(wall.isStatic(), "Static wall should be classified as static");
    TEST_ASSERT(wall.getType() == RigidBodyType::INTERNAL_OBJECT, "Default rigid body type should be internal object");
    return true;
}

static bool testAABBOverlap() {
    AABB box1{Vector2D(0.0, 0.0), Vector2D(2.0, 2.0)};
    AABB box2{Vector2D(1.5, 1.5), Vector2D(3.0, 3.0)};
    AABB box3{Vector2D(5.0, 5.0), Vector2D(6.0, 6.0)};

    TEST_ASSERT(box1.overlaps(box2), "Box 1 and Box 2 should overlap");
    TEST_ASSERT(!box1.overlaps(box3), "Box 1 and Box 3 should NOT overlap");

    return true;
}

static bool testDamBreakWallPivotParameter() {
    auto findDamWall = [](const std::vector<std::shared_ptr<RigidObject>>& bodies) {
        for (const auto& body : bodies) {
            if (body->getType() == RigidBodyType::INTERNAL_OBJECT) return body;
        }
        return std::shared_ptr<RigidObject>{};
    };

    DamBreakIC::Parameters pivotParams;
    pivotParams.enableDamWallPivot = true;
    DamBreakIC pivotIC(pivotParams);
    auto pivotWall = findDamWall(pivotIC.getRigidObjects());
    TEST_ASSERT(pivotWall != nullptr, "Dam barrier rigid body should be created");
    TEST_ASSERT(pivotWall->hasPivot(), "Enabled option should pivot the dam barrier");
    TEST_ASSERT(!pivotWall->isStatic(), "Pivoted dam barrier should be dynamic");
    TEST_ASSERT_NEAR(pivotWall->getPivotPosition().x,
                     pivotParams.damPos.x + 0.5 * pivotParams.damSize.x,
                     1e-6, "Dam pivot should be at the bottom midpoint");
    TEST_ASSERT_NEAR(pivotWall->getPivotPosition().y, pivotParams.damPos.y,
                     1e-6, "Dam pivot should be at the bottom edge");

    DamBreakIC fixedIC;
    auto fixedWall = findDamWall(fixedIC.getRigidObjects());
    TEST_ASSERT(fixedWall != nullptr, "Fixed dam barrier rigid body should be created");
    TEST_ASSERT(!fixedWall->hasPivot(), "Disabled option should not set a pivot");
    TEST_ASSERT(fixedWall->isStatic(), "Disabled option should preserve the static barrier");

    return true;
}

static bool testDamBreakWallRigidBodies() {
    DamBreakIC ic;
    auto rigidBodies = ic.getRigidObjects();

    TEST_ASSERT(!rigidBodies.empty(), "Dam-break scenario should create at least one rigid body");

    bool hasStaticBoundary = false;
    for (const auto& body : rigidBodies) {
        if (body->isStatic() && body->getType() == RigidBodyType::EXTERNAL_BOUNDARY) {
            hasStaticBoundary = true;
            break;
        }
    }

    TEST_ASSERT(hasStaticBoundary, "Dam-break retaining walls should be registered as static external boundary rigid bodies");
    return true;
}

static bool testDamDebrisFalling() {
    DamBreakIC ic;
    auto rigidBodies = ic.getRigidObjects();

    TEST_ASSERT(!rigidBodies.empty(), "Dam-break scenario should create at least one rigid body");

    bool hasStaticBoundary = false;
    for (const auto& body : rigidBodies) {
        if (body->isStatic() && body->getType() == RigidBodyType::EXTERNAL_BOUNDARY) {
            hasStaticBoundary = true;
            break;
        }
    }

    bool hitFloor = false;
    bool rebounded = false;

    RigidBodySolver solver(RigidBodySolverParameters{
        .enableGravity = true,
        .gravity = Vector2D(0.0, -9.81),
        .restitution = 0.1,
        .particleRadius = 0.05
    });

    ThreadPool pool(4);
    pool.start();

    auto falling = rigidBodies[4];

    for (int step = 0; step < 2000; ++step) {
        //solver.solveCollisions(bodies, 0.01, pool);

        double y = falling->getCenterOfMass().y;
        double x = falling->getCenterOfMass().x;
        std::cout << x << ", " << y << std::endl;

        solver.integrate(rigidBodies, 0.01);

        if (!hitFloor && y <= 0.05 && falling->getLinearVel().y < 0.0) {
            hitFloor = true;
        }
        if (hitFloor && falling->getLinearVel().y > 0.0) {
            rebounded = true;
            break;
        }
    }

    pool.Stop();
}

static bool testRigidBodySolverTwoBodyCollision() {
    RigidBodySolver solver(RigidBodySolverParameters{
        .enableGravity = false,
        .gravity = Vector2D(0.0, -9.81),
        .restitution = 0.1,
        .particleRadius = 0.05
    });

    auto bodyA = std::make_shared<RigidObject>();
    bodyA->addParticle(std::make_shared<SolidParticle>(1, Vector2D(0.00, 0.0), Vector2D(0.0, 0.0), 1.0, 1000.0, 300.0, MotionType::DYNAMIC));
    bodyA->addParticle(std::make_shared<SolidParticle>(2, Vector2D(0.10, 0.0), Vector2D(0.0, 0.0), 1.0, 1000.0, 300.0, MotionType::DYNAMIC));
    bodyA->finalizeInitialization();
    bodyA->setConstraints(false, false, false);
    bodyA->setCenterOfMass(Vector2D(0.05, 0.0));
    bodyA->setLinearVel(Vector2D(1.0, 0.0));

    auto bodyB = std::make_shared<RigidObject>();
    bodyB->addParticle(std::make_shared<SolidParticle>(3, Vector2D(0.12, 0.0), Vector2D(0.0, 0.0), 1.0, 1000.0, 300.0, MotionType::DYNAMIC));
    bodyB->addParticle(std::make_shared<SolidParticle>(4, Vector2D(0.22, 0.0), Vector2D(0.0, 0.0), 1.0, 1000.0, 300.0, MotionType::DYNAMIC));
    bodyB->finalizeInitialization();
    bodyB->setConstraints(false, false, false);
    bodyB->setCenterOfMass(Vector2D(0.17, 0.0));
    bodyB->setLinearVel(Vector2D(-0.5, 0.0));

    ThreadPool pool(4);
    pool.start();

    std::vector<std::shared_ptr<RigidObject>> bodies{bodyA, bodyB};
    bool collided = false;
    while (!collided) {
        //collided = solver.solveCollisions(bodies, 0.01, pool);
        collided = solver.integrate(bodies, 0.01);
        for (auto b : bodies)
        {
            double x = b->getCenterOfMass().x;
            std::cout << x << std::endl;
        }
    }

    pool.Stop();

    TEST_ASSERT(bodyA->getForceAccumulator().x < 0.0, "Body A should receive a negative-x impulse from the solver contact normal");
    TEST_ASSERT(bodyB->getForceAccumulator().x > 0.0, "Body B should receive the opposite positive-x impulse");
    return true;
}

static bool testRigidBodyFallsToStaticFloor() {
    RigidBodySolver solver(RigidBodySolverParameters{
        .enableGravity = true,
        .gravity = Vector2D(0.0, -9.81),
        .restitution = 0.1,
        .particleRadius = 0.05
    });

    auto floor = std::make_shared<RigidObject>();
    floor->addParticle(std::make_shared<SolidParticle>(1, Vector2D(0.0, 0.0), Vector2D(0.0, 0.0), 1.0, 1000.0, 300.0, MotionType::STATIC));
    floor->finalizeInitialization();
    floor->setConstraints(true, true, true);
    floor->setType(RigidBodyType::EXTERNAL_BOUNDARY);

    auto falling = std::make_shared<RigidObject>();
    falling->addParticle(std::make_shared<SolidParticle>(2, Vector2D(0.0, 0.25), Vector2D(0.0, 0.0), 1.0, 1000.0, 300.0, MotionType::DYNAMIC));
    falling->finalizeInitialization();
    falling->setConstraints(false, false, false);
    falling->setCenterOfMass(Vector2D(0.0, 0.25));
    falling->setLinearVel(Vector2D(0.0, 0.0));

    std::vector<std::shared_ptr<RigidObject>> bodies{floor, falling};



    bool hitFloor = false;
    bool rebounded = false;

    ThreadPool pool(4);
    pool.start();

    for (int step = 0; step < 2000; ++step) {
        //solver.solveCollisions(bodies, 0.01, pool);
        solver.integrate(bodies, 0.01);

        double y = falling->getCenterOfMass().y;
        if (!hitFloor && y <= 0.05 && falling->getLinearVel().y < 0.0) {
            hitFloor = true;
        }
        if (hitFloor && falling->getLinearVel().y > 0.0) {
            rebounded = true;
            break;
        }
    }

    pool.Stop();

    TEST_ASSERT(hitFloor, "A dynamic rigid body under gravity should fall until it contacts the static floor");
    TEST_ASSERT(rebounded, "A dynamic rigid body should bounce upward when it contacts the static floor");
    return true;
}

void registerRigidBodyTests(TestSuite& suite) {
    suite.startSection("RigidBody Module");
    suite.runTest("Finalize Initialization", testFinalizeInitialization);
    suite.runTest("Force and Impulse State", testForceAndImpulseState);
    suite.runTest("Rigid Body Rotates Around Pivot", testRigidBodyRotatesAroundPivot);
    suite.runTest("Pivot Force Uses Pivot Inertia", testPivotForceUsesPivotInertia);
    suite.runTest("Pivot Restoring Torque", testPivotRestoringTorque);
    suite.runTest("Pivot Collision Torque", testPivotCollisionProducesTorqueAboutAnchor);
    suite.runTest("Static Boundary Classification", testStaticBoundaryClassification);
    suite.runTest("AABB Overlap", testAABBOverlap);
    suite.runTest("Dam Break Wall Pivot Parameter", testDamBreakWallPivotParameter);
    //suite.runTest("DamBreak Wall Rigid Bodies", testDamBreakWallRigidBodies);
    //suite.runTest("DamBreak falling debris", testDamDebrisFalling);
    //suite.runTest("RigidBodySolver Two-Body Collision", testRigidBodySolverTwoBodyCollision);
    //suite.runTest("RigidBody Falls to Static Floor", testRigidBodyFallsToStaticFloor);
}