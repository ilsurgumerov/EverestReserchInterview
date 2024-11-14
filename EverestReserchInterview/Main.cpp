#include <iostream>
#include "internal_point.hpp"

int main() {
    setlocale(LC_ALL, "Russian");

    InternalPoint testSolver;
    if (!testSolver.readMpsFile("C:/Users/User/Downloads/test")) {
        std::cout << "Не удалось прочитать MPS файл\n";
        return 1;
    }

    testSolver.IPM(1000);
    testSolver.printX();
    std::cout << "Оптимальное значение: " << testSolver.c.transpose() * testSolver.x << "\n";

    return 0;
}
