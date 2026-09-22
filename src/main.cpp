#include <iostream>
#include "Vector2D.h"

int main() {
    std::cout << "Thermo-SPH Simulation Starting..." << std::endl;
    
    Vector2D v(1.0, 2.0);
    std::cout << "Vector norm: " << v.norm() << std::endl;

    return 0;
}