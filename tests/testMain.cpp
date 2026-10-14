#include "testSuite.h"
#include "test_Vector2D.h"
#include "test_RigidBody.h"
#include "test_Piston.h"
#include "test_Explosion.h"

int main() {
    TestSuite suite;

    registerVector2DTests(suite);
    registerRigidBodyTests(suite);
    registerPistonTests(suite);
    registerExplosionTests(suite);

    return suite.printSummary();
}