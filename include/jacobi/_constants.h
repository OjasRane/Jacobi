#pragma once
#include <bit>

namespace jacobi {
    inline constexpr auto pi = 3.14159265358979323846;
    inline constexpr auto e = 2.718281828459045;
    inline constexpr auto ln2 = 0.6931471805599453094172;
    inline constexpr auto eps = 1e-12;
    inline constexpr double Infinity = std::bit_cast<double>(2047ULL << 52);
    inline constexpr double nan = std::bit_cast<double>((2047ULL << 52) | (1ULL << 51));
    inline constexpr char magic_header[] = "CARL";
    inline constexpr std::streamsize header_length = sizeof(magic_header) - 1;
    inline constexpr char magic_footer[] = "QED";
    inline constexpr std::streamsize footer_length = sizeof(magic_footer) - 1;
}