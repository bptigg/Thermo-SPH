#include "TestSuite.h"
#include "test_Vector2D.h"
#include "test_RigidBody.h"

int main() {
    TestSuite suite;

    registerVector2DTests(suite);
    registerRigidBodyTests(suite);

    return suite.printSummary();
}