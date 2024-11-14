#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <coin/CoinMpsIO.hpp>
#include <iostream>
#include <vector>

class InternalPoint {
public:
    // Константы
    static constexpr double EPSILON = 1e-12;
    static constexpr double EPSILONEND = 1e-8;
    static constexpr double SIGMA = 1e-8;
    static constexpr double GAMMA = 0.1;
    static constexpr double INITIAL_X_VALUE = 0.1;

    // Числовые данные
    int numCols;
    int numRows;
    int stepIPM;
    double lambda;

    // Основные матрицы
    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd c;
    Eigen::VectorXd b;
    Eigen::VectorXd x;

    // Вспомогательные матрицы
    Eigen::VectorXd r;
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> D;
    Eigen::SparseMatrix<double> T;
    Eigen::VectorXd q;
    Eigen::VectorXd u;
    Eigen::VectorXd g;
    Eigen::VectorXd s;

    // Конструктор
    InternalPoint();

    // Основные методы
    void IPM(int maxStep);
    bool oneStepIPM();
    bool findLambda();
    bool optimalityCondition();
    bool readMpsFile(const std::string&);

    // Вспомогательные методы
    double norm(const Eigen::VectorXd&);
    void printX();
};
