#include <iostream>

#ifdef IS_FLOAT
#define TYPE float
#else
#ifdef IS_LONG_DOUBLE
#define TYPE long double
#else
#define TYPE double
#endif
#endif

#include <Eigen/Sparse>

#include "lab1/Lpq_norm.hpp"
#include "lab1/vector_p_norm.hpp"
#include "lab1/Shatten_norm.hpp"

#include "lab2/jacobi.hpp"

int main() {
    Eigen::SparseMatrix<TYPE> A(3, 3);
    A.coeffRef(0, 0) = 1;
    A.coeffRef(1, 1) = 1;
    A.coeffRef(2, 2) = 1;
    const auto x = jacobi(A, {3, 2, 1}, {0, 0, 0}, 1e-3);
    for (const auto idx : x) {
        std::cout << idx << '\n';
    }
    return 0;
}
