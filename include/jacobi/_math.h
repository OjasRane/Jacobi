#pragma once
#include "_data_structure.h"

namespace jacobi {
    template <typename T>
    T abs(const T&);

    double sqrt(double);

    double deg2rad(double);

    double rad2deg(double);

    long long floor(double);
    long long ceil(double);

    double exp(double);

    double ln(double);

    double pow(double, double);

    double sin(double);
    double cos(double);

    template <typename T, typename U>
    Vector<U> map(U (*)(T), const Vector<T>&);
    template <typename T, typename U>
    Vector<U> map(U (*)(T), const Row<T>&);
    template <typename T, typename U>
    Vector<U> map(U (*)(T), const Row<const T>&);
    template <typename T, typename U>
    Matrix<U> map(U (*)(T), const Matrix<T>&);
}

template <typename T>
T jacobi::abs(const T& num) {
    return num > 0 ? num : -num;
}

template <typename T, typename U>
jacobi::Vector<U> jacobi::map(U (*func)(T), const Vector<T>& i_V) {
    Vector<U> r_V(i_V.size());
    for (std::uint64_t i = 0; i < i_V.size(); ++i) {
        r_V[i] = func(i_V[i]);
    }
    return r_V;
}

template <typename T, typename U>
jacobi::Vector<U> jacobi::map(U (*func)(T), const Row<T>& i_R) {
    Vector<U> r_V(i_R.size());
    for (std::uint64_t i = 0; i < i_R.size(); ++i) {
        r_V[i] = func(i_R[i]);
    }
    return r_V;
}

template <typename T, typename U>
jacobi::Vector<U> jacobi::map(U (*func)(T), const Row<const T>& i_R) {
    Vector<U> r_V(i_R.size());
    for (std::uint64_t i = 0; i < i_R.size(); ++i) {
        r_V[i] = func(i_R[i]);
    }
    return r_V;
}

template <typename T, typename U>
jacobi::Matrix<U> jacobi::map(U (*func)(T), const Matrix<T>& i_M) {
    Matrix<U> r_M(i_M.rows(), i_M.cols());
    for (std::uint64_t i = 0; i < i_M.size(); ++i) {
        r_M[i] = func(i_M[i]);
    }
    return r_M;
}
