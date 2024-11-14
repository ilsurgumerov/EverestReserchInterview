#include <iostream>
#include "internal_point.hpp"

int main() {
    setlocale(LC_ALL, "Russian");

    InternalPoint testSolver;
    if (!testSolver.readMpsFile("C:/Users/User/Downloads/test")) {
        std::cout << "�� ������� ��������� MPS ����\n";
        return 1;
    }

    testSolver.IPM(1000);
    testSolver.printX();
    std::cout << "����������� ��������: " << testSolver.c.transpose() * testSolver.x << "\n";

    return 0;
}
