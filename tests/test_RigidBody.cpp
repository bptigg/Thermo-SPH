#include "test_RigidBody.h"
#include "RigidBody.h"
#include "Particle.h"
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

static bool testAccumulateForcesAndTorque() {
    RigidObject obj;
    auto p1 = createParticle(1, Vector2D(-1.0, 0.0), 1.0);
    auto p2 = createParticle(2, Vector2D(1.0, 0.0), 1.0);

    obj.addParticle(p1);
    obj.addParticle(p2);
    obj.finalizeInitialization();

    p1->accel = Vector2D(0.0, -10.0);
    p2->accel = Vector2D(0.0, 10.0);

    obj.accumulateForces();
    obj.updateKinematics(0.1);

    TEST_ASSERT_NEAR(obj.getLinearVel().x, 0.0, 1e-6, "Linear Vel X should be 0");
    TEST_ASSERT_NEAR(obj.getLinearVel().y, 0.0, 1e-6, "Linear Vel Y should be 0");
    TEST_ASSERT_NEAR(obj.getAngularVel(), 1.0, 1e-6, "Angular velocity calculation failed");

    return true;
}

static bool testPistonConstraints() {
    RigidObject piston;
    auto p1 = createParticle(1, Vector2D(0.0, 0.0), 1.0);
    piston.addParticle(p1);
    piston.finalizeInitialization();

    piston.setConstraints(true, false, true);

    p1->accel = Vector2D(100.0, -9.81);
    piston.accumulateForces();
    piston.updateKinematics(0.1);

    TEST_ASSERT_NEAR(piston.getLinearVel().x, 0.0, 1e-6, "Piston X motion should be locked");
    TEST_ASSERT_NEAR(piston.getLinearVel().y, -0.981, 1e-6, "Piston Y motion integration failed");
    TEST_ASSERT_NEAR(piston.getAngularVel(), 0.0, 1e-6, "Piston rotation should be locked");

    return true;
}

static bool testParticleSynchronization() {
    RigidObject obj;
    auto p1 = createParticle(1, Vector2D(0.0, 0.0), 1.0);
    auto p2 = createParticle(2, Vector2D(2.0, 0.0), 1.0);

    obj.addParticle(p1);
    obj.addParticle(p2);
    obj.finalizeInitialization();

    p1->accel = Vector2D(0.0, 10.0);
    p2->accel = Vector2D(0.0, 10.0);

    obj.accumulateForces();
    obj.updateKinematics(0.1);

    TEST_ASSERT_NEAR(p1->pos.y, 0.1, 1e-6, "Particle 1 position sync failed");
    TEST_ASSERT_NEAR(p2->pos.y, 0.1, 1e-6, "Particle 2 position sync failed");
    TEST_ASSERT_NEAR(p1->vel.y, 1.0, 1e-6, "Particle 1 velocity sync failed");
    TEST_ASSERT_NEAR(p2->vel.y, 1.0, 1e-6, "Particle 2 velocity sync failed");

    return true;
}

static bool testApplyImpulse() {
    RigidObject obj;
    auto p1 = createParticle(1, Vector2D(-1.0, 0.0), 1.0);
    auto p2 = createParticle(2, Vector2D(1.0, 0.0), 1.0);

    obj.addParticle(p1);
    obj.addParticle(p2);
    obj.finalizeInitialization();

    Vector2D impulse(0.0, 2.0);
    Vector2D r(1.0, 0.0);

    obj.applyImpulse(impulse, r);

    TEST_ASSERT_NEAR(obj.getLinearVel().y, 1.0, 1e-6, "Linear impulse resolution failed");
    TEST_ASSERT_NEAR(obj.getAngularVel(), 1.0, 1e-6, "Angular impulse resolution failed");

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

static bool test1DTwoBoxWallCollisions() {
    // Box 1 (Light, m = 1.0 kg)
    RigidObject box1;
    auto p1 = createParticle(1, Vector2D(1.0, 0.0), 1.0);
    box1.addParticle(p1);
    box1.finalizeInitialization();
    box1.setConstraints(false /*lockX*/, true /*lockY*/, true /*lockRotation*/);

    // Box 2 (Heavy, m = 100.0 kg, 100x heavier)
    RigidObject box2;
    auto p2 = createParticle(2, Vector2D(5.0, 0.0), 100.0);
    box2.addParticle(p2);
    box2.finalizeInitialization();
    box2.setConstraints(false /*lockX*/, true /*lockY*/, true /*lockRotation*/);

    // Initial state: Box 1 at rest, Box 2 moving left towards Box 1
    // Apply initial impulse to set Box 2 velocity to -1.0 m/s
    box2.applyImpulse(Vector2D(-100.0, 0.0), Vector2D(0.0, 0.0));

    int collisionCount = 0;
    const double dt = 0.001; // Small time step for smooth kinematic update
    const int maxSteps = 500000; // Safety guard against infinite loops
    int step = 0;

    double m1 = box1.getTotalMass(); // 1.0
    double m2 = box2.getTotalMass(); // 100.0

    while (step++ < maxSteps) {
        // Integrate positions
        box1.updateKinematics(dt);
        box2.updateKinematics(dt);

        double x1 = box1.getCenterOfMass().x;
        double x2 = box2.getCenterOfMass().x;
        double v1 = box1.getLinearVel().x;
        double v2 = box2.getLinearVel().x;

        // 1. Collision between Box 1 and Wall at x = 0
        if (x1 <= 0.0 && v1 < 0.0) {
            // Elastic impulse from impenetrable wall: J = -2 * m1 * v1
            box1.applyImpulse(Vector2D(-2.0 * m1 * v1, 0.0), Vector2D(0.0, 0.0));
            collisionCount++;
            continue;
        }

        // 2. Collision between Box 1 and Box 2 (x1 >= x2 and converging: v1 > v2)
        if (x1 >= x2 && (v1 - v2) > 0.0) {
            // 1D Elastic Impulse: J = 2 * m1 * m2 * (v2 - v1) / (m1 + m2)
            double J_x = (2.0 * m1 * m2 * (v2 - v1)) / (m1 + m2);

            box1.applyImpulse(Vector2D(J_x, 0.0), Vector2D(0.0, 0.0));
            box2.applyImpulse(Vector2D(-J_x, 0.0), Vector2D(0.0, 0.0));
            collisionCount++;
            continue;
        }

        // TERMINATION CONDITION:
        // Both boxes are moving away from the wall (v1 > 0, v2 > 0) AND 
        // Box 2 is moving at least as fast as Box 1 (v2 >= v1),
        // meaning Box 2 will forever pull away from Box 1 and no further collisions can occur.
        if (v1 > 0.0 && v2 > 0.0 && v2 >= v1) {
            break;
        }
    }

    // Assert loop did not time out
    TEST_ASSERT(step < maxSteps, "Simulation timed out before reaching terminal state");

    // For mass ratio 100:1, the exact number of collisions must be 31 (floor(pi * 10))
    TEST_ASSERT(collisionCount == 31, "Collision count for 100:1 mass ratio should equal 31");

    // Final velocity check: Box 2 must be moving faster than Box 1 away from wall
    double final_v1 = box1.getLinearVel().x;
    double final_v2 = box2.getLinearVel().x;

    TEST_ASSERT(final_v2 >= final_v1, "Box 2 must be moving faster than Box 1 at termination");
    TEST_ASSERT(final_v1 > 0.0, "Box 1 must be moving away from wall at termination");

    return true;
}

void registerRigidBodyTests(TestSuite& suite) {
    suite.startSection("RigidBody Module");
    suite.runTest("Finalize Initialization", testFinalizeInitialization);
    suite.runTest("Accumulate Forces and Torque", testAccumulateForcesAndTorque);
    suite.runTest("Piston Kinematic Constraints", testPistonConstraints);
    suite.runTest("Particle Synchronization", testParticleSynchronization);
    suite.runTest("Apply Impulse", testApplyImpulse);
    suite.runTest("AABB Overlap", testAABBOverlap);
    suite.runTest("1D Two-Box Wall Collisions (Pi Test)", test1DTwoBoxWallCollisions);
}