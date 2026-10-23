#include <cstddef>
#include <cstdint>
#include <bit>
#include <stdexcept>
#include "jacobi/_math.h"
#include "jacobi/_constants.h"

double jacobi::sqrt(const double num) {
    if (num < 0) {
        throw std::domain_error("Cannot compute square root of a negative number");
    }
    double prev = num;
    double sqrt = num / 2.0;
    while (abs(prev - sqrt) > eps) {
        prev = sqrt;
        sqrt = (sqrt * sqrt + num) / (2 * sqrt);
    }
    return sqrt;
}

double jacobi::deg2rad(const double theta) {
    return theta * (pi / 180.0);
}

double jacobi::rad2deg(const double theta) {
    return (180.0 / pi) * theta;
}

long long jacobi::floor(const double num){
    long long floor_num = static_cast<long long>(num);
    if (num < 0 && floor_num !=  num) {
        --floor_num;
    }
    return floor_num;
}

long long jacobi::ceil(const double num) {
    const long long floor_num = floor(num);
    return floor_num != num ? floor_num + 1 : floor_num;
}

double jacobi::exp(const double x) {
    const long long integer_part = floor(x);
    const double frac_part = x - floor(x);
    double exp_int{1.0};
    const bool isSigned{x < 0};
    for (std::uint64_t i = 0; i < integer_part; ++i) {
        exp_int *= e;
    }
    if (isSigned) {
        exp_int = 1 / exp_int;
    }
    double prev = e;
    double term{1.0};
    double exp_frac{term};
    double i{0.0};
    while (abs(prev - exp_frac) > eps) {
        prev = exp_frac;
        ++i;
        term *= frac_part / i;
        exp_frac += term;
    }
    return exp_int * exp_frac;
}

double jacobi::ln(const double num) {
    /*
     * Here for finding $\ln(x)$, first I converted $x=a\times 2^b$ by manipulating bits for double. Refer IEEE 754 double precision floating point format on Wikipedia.
    */
    if (num <= 0) {
        throw std::domain_error("Natural Logarithm does not output real values for 0 and negative real numbers");
    }
    const std::uint64_t bits = std::bit_cast<std::uint64_t>(num);
    const std::uint64_t b = ((bits >> 52) & 0x7FFULL) - 1023 + 1;
    const double a = std::bit_cast<double>((bits & 0xFFFFFFFFFFFFFULL) | 0x3FEULL << 52);
    double lna{0.0};
    double prev{2.0};
    double x = a - 1.0;
    double sign = 1.0;
    double i{1.0};
    while (abs(prev - lna) > eps) {
        prev = lna;
        lna += sign * x / i;
        sign *= -1.0;
        x *= a - 1.0;
        ++i;
    }
    return b * ln2 + lna;
}

double jacobi::pow(const double x, const double y) {
    return exp(y*ln(x));
}

double jacobi::sin(double theta) {
    bool isSigned{false};
    if (theta < 0) {
        isSigned = true;
        theta = -theta;
    }
    while (theta >= 2*pi) {
        theta -= 2*pi;
    }
    if (theta >= pi) {
        isSigned = true;
        theta -= pi;
    }
    double prev{2.0};
    double term{theta};
    double sin{term};
    double i{1.0};
    while (abs(prev - sin) > eps) {
        prev = sin;
        term *= -1 * theta * theta / (2*i*(2*i+1));
        sin += term;
        ++i;
    }
    return isSigned ? -sin : sin;
}

double jacobi::cos(double theta) {
    theta = abs(theta);
    while (theta >= 2*pi) {
        theta -= 2*pi;
    }
    bool sign{false};
    if (theta > pi) {
        theta -= pi;
        sign = true;
    }
    double prev{2.0};
    double term{1.0};
    double i{1.0};
    double cos{term};
    while (abs(prev - cos) > eps) {
        prev = cos;
        term *= -theta * theta / ((2*i-1) * 2*i);
        cos += term;
        ++i;
    }
    return sign ? -cos : cos;
}