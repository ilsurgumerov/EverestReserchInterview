#include <iostream>
#include "internal_point.hpp"

int main() {
    InternalPoint testSolver;
    if (!testSolver.readMpsFile("C:/Users/User/Downloads/test")) {
        std::cout << "Failed to download MPS file\n";
        return 1;
    }

    testSolver.IPM(1000);
    testSolver.printX();
    std::cout << "Optimal value: " << testSolver.c.transpose() * testSolver.x << "\n";

    return 0;
}
