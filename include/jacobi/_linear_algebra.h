#pragma once
#include "_data_structure.h"
#include "_math.h"

namespace jacobi::linalg {
    template <typename T>
    double norm(const Vector<T>&);
    template <typename T>
    double norm(const Matrix<T>&);
}

template <typename T>
double jacobi::linalg::norm(const Vector<T>& i_V) {
    T sum{};
    for (std::size_t i = 0; i < i_V.size(); ++i) {
        sum += i_V[i] * i_V[i];
    }
    return jacobi::sqrt(sum);
}

template <typename T>
double jacobi::linalg::norm(const Matrix<T>& i_M) {
    T sum{};
    for (std::size_t i = 0; i < i_M.size(); ++i) {
        sum += i_M[i];
    }
    return jacobi::sqrt(sum);
}