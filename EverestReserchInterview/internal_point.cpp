#include "internal_point.hpp"

InternalPoint::InternalPoint() : numCols(0), numRows(0), stepIPM(0), lambda(0) {}

bool InternalPoint::readMpsFile(const std::string& filepath) {
    CoinMpsIO m;
    if (m.readMps(filepath.c_str(), "mps") != 0) {
        std::cout << "Ошибка при чтении MPS файла\n";
        return false;
    }

    numCols = m.getNumCols();
    numRows = m.getNumRows();

    c = Eigen::VectorXd(numCols);
    b = Eigen::VectorXd(numRows);
    A.resize(numRows, numCols);

    const double* cCoin = m.getObjCoefficients();
    for (int i = 0; i < numCols; ++i) {
        c(i) = cCoin[i];
    }

    const double* bCoin = m.getRightHandSide();
    for (int i = 0; i < numRows; ++i) {
        b(i) = bCoin[i];
    }

    std::vector<Eigen::Triplet<double>> triplets;
    const CoinPackedMatrix* ACoin = m.getMatrixByCol();
    int k = 0;

    for (int col = 0; col < numCols; ++col) {
        int countInRow = ACoin->getVectorSize(col);
        for (int j = 0; j < countInRow; ++j, ++k) {
            triplets.emplace_back(ACoin->getIndices()[k], col, ACoin->getElements()[k]);
        }
    }

    A.setFromTriplets(triplets.begin(), triplets.end());
    return true;
}

void InternalPoint::IPM(int maxStep) {
    x = Eigen::VectorXd::Constant(numCols, INITIAL_X_VALUE);
    r = b - A * x;

    for (stepIPM = 0; stepIPM < maxStep; ++stepIPM) {
        if (oneStepIPM() || (norm(r) < EPSILONEND && optimalityCondition())) {
            std::cout << stepIPM << " шагов\n";
            break;
        }
    }
}

bool InternalPoint::oneStepIPM() {
    D = x.array().square().matrix().asDiagonal();
    T = A * D * A.transpose();
    q = r + A * D * c;

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(T);
    solver.factorize(T);
    u = solver.solve(q);

    g = c - A.transpose() * u;
    s = (-1) * D * g;

    if (!findLambda()) {
        return true;
    }

    x += lambda * s;
    r = b - A * x;

    return false;
}

bool InternalPoint::findLambda() {
    double lambda_ = 0;
    bool findNegative = false;
    int countZeros = 0;

    for (int i = 0; i < numCols; i++) {
        if (s(i) < 0) {
            lambda_ = (lambda_ == 0) ? -x(i) / s(i) : std::min(lambda_, -x(i) / s(i));
            findNegative = true;
        }
        else if (s(i) < EPSILON) {
            ++countZeros;
        }
    }

    if (!findNegative) // s >= 0
    {
        if (norm(r) < EPSILON)   // r == 0
            if (countZeros == numCols) {    // s == 0
                std::cout << "Решение найдено при s == 0" << "\n";
                return false;
            }
            else {                          // s > 0
                std::cout << "Нет решений" << "\n";
                return false;
            }
        else                     // r != 0
            lambda = 1.;
    }
    else               // s < 0
    {
        if (norm(r) < EPSILON)   // r == 0
            lambda = GAMMA * lambda_;
        else                     // r != 0
            lambda = std::min(1., GAMMA * lambda_);
    }

    return true;
}

bool InternalPoint::optimalityCondition() {
    return (D * g).transpose() * g < SIGMA;
}

void InternalPoint::printX() {
    std::cout << "Вектор X:\n" << x << "\n\n";
}

double InternalPoint::norm(const Eigen::VectorXd& y) {
    return y.cwiseAbs().maxCoeff();
}
