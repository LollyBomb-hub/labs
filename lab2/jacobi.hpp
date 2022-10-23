#ifndef __LABS__JACOBI_METHOD__HPP__

#define __LABS__JACOBI_METHOD__HPP__

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#ifndef TYPE
#error "No type defined earlier"
#endif

// A * x = b
std::vector<TYPE> jacobi(const Eigen::SparseMatrix<TYPE>& A, const std::vector<TYPE>& b, const std::vector<TYPE>& x_0, const TYPE& eps) {
    if (eps < 0) {
        throw std::runtime_error("Could not process negative eps!");
    }
    const auto matrixA = A.toDense();
    const auto size = matrixA.cols();
    if (size != matrixA.rows()) {
        throw std::runtime_error("Matrix A must be square");
    }
    if (size != b.size()) {
        throw std::runtime_error("Vector b must be of cols[or rows] count size!");
    }
    if (size != x_0.size()) {
        throw std::runtime_error("Vector x_0 must be of cols[or rows] count size!");
    }
    const auto q = matrixA.norm();
    TYPE condition;
    if (((1 - q) / q) < 0) {
        condition = eps;
    } else { condition = ((1 - q) / q) * eps; }
    std::cout << "Condition " << condition << ", q = " << q << '\n';
    double norm;
    std::vector<TYPE> result(size);
    std::vector<TYPE> current = x_0;
    do {
        for (int i = 0; i < size; i++) {
            result[i] = b[i];
            for (int g = 0; g < size; g++) {
                if (i != g)
                    result[i] -= matrixA.coeff(i, g) * current[g];
            }
            result[i] /= matrixA.coeff(i, i);
        }
        norm = fabs(current[0] - result[0]);
        for (int h = 0; h < size; h++) {
            if (fabs(current[h] - result[h]) > norm) { norm = fabs(current[h] - result[h]); }
            current[h] = result[h];
        }
    } while (norm > condition);
    return result;
}

#endif