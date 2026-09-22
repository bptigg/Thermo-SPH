#include "test_RigidBody.h"
#include "RigidObject.h"
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
    auto p2 = createParticle(2, Vector2D(2.0, 0.0), 2.0);

    obj.addParticle(p1);
    obj.addParticle(p2);
    obj.finalizeInitialization();

    TEST_ASSERT_NEAR(obj.getTotalMass(), 4.0, 1e-6, "Total mass calculation failed");
    TEST_ASSERT_NEAR(obj.getCenterOfMass().x, 1.0, 1e-6, "Center of Mass X failed");
    TEST_ASSERT_NEAR(obj.getCenterOfMass().y, 0.0, 1e-6, "Center of Mass Y failed");
    TEST_ASSERT_NEAR(obj.getInertia(), 4.0, 1e-6, "Moment of Inertia calculation failed");

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

void registerRigidBodyTests(TestSuite& suite) {
    suite.startSection("RigidBody Module");
    suite.runTest("Finalize Initialization", testFinalizeInitialization);
    suite.runTest("Accumulate Forces and Torque", testAccumulateForcesAndTorque);
    suite.runTest("Piston Kinematic Constraints", testPistonConstraints);
    suite.runTest("Particle Synchronization", testParticleSynchronization);
    suite.runTest("Apply Impulse", testApplyImpulse);
    suite.runTest("AABB Overlap", testAABBOverlap);
}