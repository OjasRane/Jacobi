#pragma once

namespace jacobi {
    template <typename T>
    T abs(const T&);

    double sqrt(double);

    double deg2rad(double);

    double rad2deg(double);

    long long floor(double);

    double exp(double);

    double ln(double);

    double pow(double, double);

    double sin(double);
    double cos(double);
}

template <typename T>
T jacobi::abs(const T& num) {
    return num > 0 ? num : -num;
}
