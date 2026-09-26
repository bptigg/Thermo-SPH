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
    suite.runTest("Static Boundary Classification", testStaticBoundaryClassification);
    suite.runTest("AABB Overlap", testAABBOverlap);
    suite.runTest("DamBreak Wall Rigid Bodies", testDamBreakWallRigidBodies);
    suite.runTest("RigidBodySolver Two-Body Collision", testRigidBodySolverTwoBodyCollision);
    suite.runTest("RigidBody Falls to Static Floor", testRigidBodyFallsToStaticFloor);
}