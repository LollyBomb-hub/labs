#ifndef __LABS__SHATTEN_NORM__HPP__

#define __LABS__SHATTEN_NORM__HPP__

#include <Eigen/SVD>
#include <Eigen/Sparse>

#ifndef TYPE
#error "No type defined earlier"
#endif

TYPE shatten_norm(const Eigen::SparseMatrix<TYPE>& matrix, const TYPE& p) {
    if (p < 1) {
        throw std::runtime_error("In Shatten norm p was less than 1!");
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> svd;
    svd.compute(matrix.toDense(), Eigen::ComputeThinV | Eigen::ComputeThinU);
    const auto inverse_p = 1 / p;
    const auto& singular_values_vector = svd.singularValues();
    const auto count_values = singular_values_vector.size();
    TYPE result = 0;
    for (long long id_value = 0; id_value < count_values; id_value++) {
        result += pow(singular_values_vector[id_value], p);
    }
    return pow(result, inverse_p);
}

#endif