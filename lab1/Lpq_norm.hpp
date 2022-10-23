#ifndef __LABS__LPQ_NORM__HPP__

#define __LABS__LPQ_NORM__HPP__

#include <Eigen/Sparse>

#ifndef TYPE
#error "No type defined earlier"
#endif

TYPE lpq_norm(const Eigen::SparseMatrix<TYPE>& matrix, const TYPE& p, const TYPE& q) {
    if (p < 1) {
        throw std::runtime_error("In Lpq norm p was less than 1!");
    }
    if (q < 1) {
        throw std::runtime_error("In Lpq norm q was less than 1!");
    }
    const auto q_div_p = q / p;
    const auto inverse_q = 1 / q;
    const auto rows = matrix.rows();
    const auto cols = matrix.cols();
    TYPE result = 0;
    for (long long id_col = 0; id_col < cols; id_col++) {
        TYPE column_sum = 0;
        for (long long id_row = 0; id_row < rows; id_row++) {
            const auto matrix_i_j = fabs(matrix.coeff(id_row, id_col));
            column_sum += pow(matrix_i_j, p);
        }
        result += pow(column_sum, q_div_p);
    }
    return pow(result, inverse_q);
}

#endif