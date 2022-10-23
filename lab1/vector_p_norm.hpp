#ifndef __LABS__VECTOR_P_NORM__HPP__

#define __LABS__VECTOR_P_NORM__HPP__

#include <Eigen/Sparse>

#ifndef TYPE
#error "No type defined earlier"
#endif

TYPE p_norm(const Eigen::SparseMatrix<TYPE> &matrix, const TYPE &p) {
    if (p < 1) {
        throw std::runtime_error("In vector p norm p was less than 1!");
    }
    const auto inverse_p = 1 / p;
    const auto rows = matrix.rows();
    const auto cols = matrix.cols();
    TYPE result = 0;
    for (long long id_col = 0; id_col < cols; id_col++) {
        for (long long id_row = 0; id_row < rows; id_row++) {
            const auto matrix_i_j = fabs(matrix.coeff(id_row, id_col));
            result += pow(matrix_i_j, p);
        }
    }
    return pow(result, inverse_p);
}

#endif