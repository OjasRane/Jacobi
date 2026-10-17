#pragma once
#include <iostream>
#include <initializer_list>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <fstream>
#include <cstdint>
#include <cstring>
#include "_constants.h"

namespace jacobi {
    template <typename T>
    class Matrix;

    template <typename T>
    class Vector {
    private:
        std::uint64_t m_Length;
        T* m_Data;
    public:
        Vector(std::uint64_t);
        Vector(const Vector&);
        Vector(std::initializer_list<T>);
        Vector& operator=(const Vector&);

        T& operator[](std::uint64_t);
        const T& operator[](std::uint64_t) const;

        Vector operator+(const Vector&) const;
        Vector operator+(const T&) const;

        Vector operator-(const Vector&) const;
        Vector operator-() const;
        Vector operator-(const T&) const;

        T inner_product(const Vector&) const;
        T operator*(const Vector& i_V) const { return inner_product(i_V); }
        Vector operator*(const Matrix<T>& i_M) const;
        Vector operator*(const T&) const;

        Vector operator/(const T&) const;
        Vector operator/(const Vector&) const;

        Vector hadamard_product(const Vector&) const;

        Matrix<T> outer_product(const Vector&) const;

        [[nodiscard]] std::uint64_t size() const { return m_Length; }

        void save(const std::string&, const std::uint8_t = 1) const;

        ~Vector();
    };

    template <typename T>
    Vector<T> operator+(const T& scaler, const Vector<T>& i_V) { return i_V + scaler; }

    template <typename T>
    Vector<T> operator-(const T& scaler, const Vector<T>& i_V) { return (-i_V) + scaler; }

    template <typename T>
    Vector<T> operator*(const T& scaler, const Vector<T>& i_V) { return i_V * scaler; }

    template <typename T>
    Vector<T> operator/(const T&, const Vector<T>&);

    template <typename T>
    Vector<T> zeros(std::uint64_t);
    template <typename T>
    Vector<T> ones(std::uint64_t);

    template <typename T>
    class Matrix {
    private:
        std::uint64_t m_Rows, m_Cols;
        T* m_Data;
    public:
        Matrix(std::uint64_t, std::uint64_t);
        Matrix(const Matrix&);
        Matrix(std::initializer_list<std::initializer_list<T>>);
        Matrix(const Vector<T>&);
        Matrix& operator=(const Matrix&);

        T& operator()(const std::uint64_t row, const std::uint64_t col) { return m_Data[col + row * m_Cols]; }
        const T& operator()(const std::uint64_t row, const std::uint64_t col) const { return m_Data[col + row * m_Cols]; }
        T& operator[](const std::uint64_t n_element) { return m_Data[n_element]; }
        const T& operator[](const std::uint64_t n_element) const { return m_Data[n_element]; }

        Matrix operator+(const Matrix&) const;
        Matrix operator+(const T&) const;

        Matrix operator-(const Matrix&) const;
        Matrix operator-() const;
        Matrix operator-(const T&) const;

        Matrix operator*(const Matrix&) const;
        Matrix operator*(const T&) const;
        Vector<T> operator*(const Vector<T>&) const;

        Matrix operator/(const T&) const;
        Matrix operator/(const Matrix&) const;

        Matrix hadamard_product(const Matrix&) const;

        Matrix transpose() const;

        Vector<T> flatten() const;

        std::uint64_t size() const { return m_Rows * m_Cols; }
        std::uint64_t rows() const { return m_Rows; }
        std::uint64_t cols() const { return m_Cols; }

        void save(const std::string&, const std::uint8_t = 1) const;

        ~Matrix();
    };

    template <typename T>
    Matrix<T> operator+(const T& scaler, const Matrix<T>& i_M) { return i_M + scaler; }

    template <typename T>
    Matrix<T> operator-(const T& scaler, const Matrix<T>& i_M) { return (-i_M) + scaler; }

    template <typename T>
    Matrix<T> operator*(const T& scaler, const Matrix<T>& i_M) { return i_M * scaler; }

    template <typename T>
    Matrix<T> operator/(const T&, const Matrix<T>&);

    template <typename T>
    Matrix<T> zeros(std::uint64_t, std::uint64_t);
    template <typename T>
    Matrix<T> ones(std::uint64_t, std::uint64_t);
    template <typename T>
    Matrix<T> eyes(std::uint64_t);

    template <typename T>
    Vector<T> load_vector(const std::string&);
    template <typename T>
    Matrix<T> load_matrix(const std::string&);
}

template <typename T>
jacobi::Vector<T>::Vector(const std::uint64_t length)
    : m_Length(length) {
    m_Data = new T[m_Length];
}

template <typename T>
jacobi::Vector<T>::Vector(const Vector& i_V)
    : m_Length(i_V.m_Length) {
    m_Data = new T[m_Length];
    for (std::uint64_t i = 0; i < m_Length; ++i) {
        m_Data[i] = i_V.m_Data[i];
    }
}

template <typename T>
jacobi::Vector<T>::Vector(const std::initializer_list<T> values)
    : m_Length(values.size()) {
    m_Data = new T[m_Length];
    std::uint64_t i = 0;
    for (const T& value : values) {
        m_Data[i++] = value;
    }
}

template <typename T>
jacobi::Vector<T>& jacobi::Vector<T>::operator=(const Vector& i_V) {
    if (this != &i_V) {
        T* new_data = new T[m_Length];
        for (std::uint64_t i = 0; i < m_Length; ++i) {
            new_data[i] = i_V.m_Data[i];
        }
        delete[] m_Data;
        m_Data = new_data;
        new_data = nullptr;
        m_Length = i_V.m_Length;
    }
    return *this;
}

template <typename T>
T& jacobi::Vector<T>::operator[](const std::uint64_t idx) {
    return m_Data[idx];
}

template <typename T>
const T& jacobi::Vector<T>::operator[](const std::uint64_t idx) const {
    return m_Data[idx];
}

template <typename T>
jacobi::Vector<T> jacobi::Vector<T>::operator+(const Vector& i_V) const {
    if (m_Length != i_V.m_Length) {
        throw std::length_error("Length of both vectors must be equal");
    }
    Vector r_V(m_Length);
    for (std::uint64_t i = 0; i < m_Length; ++i) {
        r_V.m_Data[i] = m_Data[i] + i_V.m_Data[i];
    }
    return r_V;
}

template <typename T>
jacobi::Vector<T> jacobi::Vector<T>::operator+(const T& scaler) const {
    Vector r_V(m_Length);
    for (std::uint64_t i = 0; i < m_Length; ++i) {
        r_V.m_Data[i] = m_Data[i] + scaler;
    }
    return r_V;
}

template <typename T>
jacobi::Vector<T> jacobi::Vector<T>::operator-(const Vector& i_V) const {
    if (m_Length != i_V.m_Length) {
        throw std::length_error("Length of both vectors must be equal");
    }
    Vector r_V(m_Length);
    for (std::uint64_t i = 0; i < m_Length; ++i) {
        r_V.m_Data[i] = m_Data[i] - i_V.m_Data[i];
    }
    return r_V;
}

template <typename T>
jacobi::Vector<T> jacobi::Vector<T>::operator-() const {
    Vector r_V(m_Length);
    for (std::uint64_t i = 0; i < m_Length; ++i) {
        r_V.m_Data[i] = -m_Data[i];
    }
    return r_V;
}

template <typename T>
jacobi::Vector<T> jacobi::Vector<T>::operator-(const T& scaler) const {
    Vector r_V(m_Length);
    for (std::uint64_t i = 0; i < m_Length; ++i) {
        r_V.m_Data[i] = m_Data[i] - scaler;
    }
    return r_V;
}

template <typename T>
T jacobi::Vector<T>::inner_product(const Vector& i_V) const {
    if (m_Length != i_V.m_Length) {
        throw std::length_error("Length of both vectors must be equal");
    }
    T product{};
    for (std::uint64_t i = 0; i < m_Length; ++i) {
        product += m_Data[i] * i_V.m_Data[i];
    }
    return product;
}

template <typename T>
jacobi::Vector<T> jacobi::Vector<T>::operator*(const Matrix<T>& i_M) const {
    if (m_Length != i_M.rows()) {
        throw std::domain_error("Attempt to multipy non-compatible vector-matrix");
    }
    Vector r_V(m_Length);
    for (std::uint64_t i = 0; i < i_M.cols(); ++i) {
        r_V[i] = T{};
        for (std::uint64_t j = 0; j < m_Length; ++j) {
            r_V[i] += m_Data[j] * i_M(j, i);
        }
    }
    return r_V;
}

template <typename T>
jacobi::Vector<T> jacobi::Vector<T>::operator*(const T& scaler) const {
    Vector r_V(m_Length);
    for (std::uint64_t i = 0; i < m_Length; ++i) {
        r_V.m_Data[i] = m_Data[i] * scaler;
    }
    return r_V;
}

template <typename T>
jacobi::Vector<T> jacobi::Vector<T>::operator/(const T& scaler) const {
    Vector r_V(m_Length);
    for (std::uint64_t i = 0; i < m_Length; ++i) {
        r_V.m_Data[i] = m_Data[i] / scaler;
    }
    return r_V;
}

template <typename T>
jacobi::Vector<T> jacobi::Vector<T>::operator/(const Vector& i_V) const {
    if (m_Length != i_V.m_Length) {
        throw std::length_error("Length of both vectors must be equal");
    }
    Vector r_V(m_Length);
    for (std::uint64_t i = 0; i < m_Length; ++i) {
        r_V.m_Data[i] = m_Data[i] / i_V.m_Data[i];
    }
    return r_V;
}

template<typename T>
jacobi::Vector<T> jacobi::Vector<T>::hadamard_product(const Vector& i_V) const {
    if (m_Length != i_V.m_Length) {
        throw std::length_error("Length of both vectors must be equal");
    }
    Vector r_V(m_Length);
    for (std::uint64_t i = 0; i < m_Length; ++i) {
        r_V.m_Data[i] = m_Data[i] * i_V.m_Data[i];
    }
    return r_V;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const jacobi::Vector<T>& i_V) {
    constexpr int n{10};
    os << '[' << i_V[0];
    for (std::uint64_t i = 1; i < i_V.size(); ++i) {
        os << std::setw(n) << i_V[i];
    }
    os << "]\n";
    return os;
}

template <typename T>
jacobi::Matrix<T> jacobi::Vector<T>::outer_product(const Vector& i_V) const {
    Matrix<T> r_M(m_Length, i_V.size());
    for (std::uint64_t i = 0; i < m_Length; ++i) {
        for (std::uint64_t j = 0; j < i_V.size(); ++j) {
            r_M(i, j) = m_Data[i] * i_V[j];
        }
    }
    return r_M;
}

template <typename T>
void jacobi::Vector<T>::save(const std::string& path, const std::uint8_t version) const {
    std::ofstream file(path, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Could not create file " + path);
    }
    file.write(reinterpret_cast<const char*>(&magic_header), header_length);
    file.write(reinterpret_cast<const char*>(&version), sizeof(version));

    file.put('V');
    file.write(reinterpret_cast<const char*>(&m_Length), sizeof(m_Length));
    file.write(reinterpret_cast<const char*>(m_Data), sizeof(T) * m_Length);

    file.write(reinterpret_cast<const char*>(&magic_footer), footer_length);
    file.close();
}

template <typename T>
jacobi::Vector<T>::~Vector() {
    delete[] m_Data;
    m_Data = nullptr;
}

template<typename T>
jacobi::Vector<T> jacobi::operator/(const T& scaler, const Vector<T>& i_V) {
    Vector<T> r_V(i_V.size());
    for (std::uint64_t i = 0; i < i_V.size(); ++i) {
        r_V[i] = scaler / i_V[i];
    }
    return r_V;
}

template <typename T>
jacobi::Vector<T> jacobi::zeros(std::uint64_t length) {
    Vector<T> r_V(length);
    for (std::uint64_t i = 0; i < length; ++i) {
        r_V[i] = T{};
    }
    return r_V;
}

template <typename T>
jacobi::Vector<T> jacobi::ones(std::uint64_t length) {
    Vector<T> r_V(length);
    for (std::uint64_t i = 0; i < length; ++i) {
        r_V[i] = 1;
    }
    return r_V;
}

template <typename T>
jacobi::Matrix<T>::Matrix(const std::uint64_t rows, const std::uint64_t cols)
    : m_Rows(rows), m_Cols(cols) {
    m_Data = new T[m_Rows * m_Cols];
}

template <typename T>
jacobi::Matrix<T>::Matrix(const Matrix& i_M)
    : m_Rows(i_M.m_Rows), m_Cols(i_M.m_Cols) {
    m_Data = new T[m_Rows * m_Cols];
    for (std::uint64_t i = 0; i < m_Rows * m_Cols; ++i) {
        m_Data[i] = i_M.m_Data[i];
    }
}

template <typename T>
jacobi::Matrix<T>::Matrix(std::initializer_list<std::initializer_list<T>> values)
    : m_Rows(values.size()), m_Cols(values.begin()->size()) {
    m_Data = new T[m_Rows * m_Cols];
    std::uint64_t i = 0;
    for (const auto& row : values) {
        std::uint64_t j = 0;
        for (const auto& element : row) {
            m_Data[i * m_Cols + j] = element;
            ++j;
        }
        ++i;
    }
}

template <typename T>
jacobi::Matrix<T>::Matrix(const Vector<T>& i_V)
    : m_Rows(i_V.size()), m_Cols(1) {
    m_Data = new T[m_Rows * m_Cols];
    for (std::uint64_t i = 0; i < m_Rows * m_Cols; ++i) {
        m_Data[i] = i_V[i];
    }
}

template <typename T>
jacobi::Matrix<T>& jacobi::Matrix<T>::operator=(const Matrix& i_M) {
    if (this != &i_M) {
        if (m_Rows != i_M.m_Rows || m_Cols != i_M.m_Cols) {
            throw std::domain_error("Cannot assign matrices with different shapes");
        }
        T* new_data = new T[i_M.m_Rows * i_M.m_Cols];
        for (std::uint64_t i = 0; i < i_M.m_Rows * i_M.m_Cols; ++i) {
            new_data[i] = i_M.m_Data[i];
        }
        delete[] m_Data;
        m_Data = new_data;
        new_data = nullptr;
        m_Rows = i_M.m_Rows;
        m_Cols = i_M.m_Cols;
    }
    return *this;
}

template <typename T>
jacobi::Matrix<T> jacobi::Matrix<T>::operator+(const Matrix& i_M) const {
    if (m_Rows != i_M.m_Rows || m_Cols != i_M.m_Cols) {
        throw std::length_error("Shape of both matrix must be equal");
    }
    Matrix r_M(m_Rows, m_Cols);
    for (std::uint64_t i = 0; i < m_Rows * m_Cols; ++i) {
        r_M.m_Data[i] = m_Data[i] + i_M.m_Data[i];
    }
    return r_M;
}

template<typename T>
jacobi::Matrix<T> jacobi::Matrix<T>::operator+(const T& scaler) const {
    Matrix r_M(m_Rows, m_Cols);
    for (std::uint64_t i = 0; i < m_Rows * m_Cols; ++i) {
        r_M.m_Data[i] = m_Data[i] + scaler;
    }
    return r_M;
}

template <typename T>
jacobi::Matrix<T> jacobi::Matrix<T>::operator-(const Matrix& i_M) const {
    if (m_Rows != i_M.m_Rows || m_Cols != i_M.m_Cols) {
        throw std::length_error("Shape of both matrix must be equal");
    }
    Matrix r_M(m_Rows, m_Cols);
    for (std::uint64_t i = 0; i < m_Rows * m_Cols; ++i) {
        r_M.m_Data[i] = m_Data[i] - i_M.m_Data[i];
    }
    return r_M;
}

template <typename T>
jacobi::Matrix<T> jacobi::Matrix<T>::operator-() const {
    Matrix r_M(m_Rows, m_Cols);
    for (std::uint64_t i = 0; i < m_Rows * m_Cols; ++i) {
        r_M.m_Data[i] = -m_Data[i];
    }
    return r_M;
}

template<typename T>
jacobi::Matrix<T> jacobi::Matrix<T>::operator-(const T& scaler) const {
    Matrix r_M(m_Rows, m_Cols);
    for (std::uint64_t i = 0; i < m_Rows * m_Cols; ++i) {
        r_M.m_Data[i] = m_Data[i] - scaler;
    }
    return r_M;
}

template <typename T>
jacobi::Matrix<T> jacobi::Matrix<T>::operator*(const Matrix& i_M) const {
    if (m_Cols != i_M.m_Rows) {
        throw std::domain_error("Attempt to multiply non-compatible matrices");
    }
    Matrix r_M(m_Rows, i_M.m_Cols);
    for (std::uint64_t i = 0; i < m_Rows; ++i) {
        for (std::uint64_t j = 0; j < i_M.m_Cols; ++j) {
            r_M(i, j) = T{};
            for (std::uint64_t k = 0; k < m_Cols; ++k) {
                r_M(i, j) += m_Data[i * m_Cols + k] * i_M.m_Data[j + k * m_Rows];
            }
        }
    }
    return r_M;
}

template <typename T>
jacobi::Matrix<T> jacobi::Matrix<T>::operator*(const T& scaler) const {
    Matrix r_M(m_Rows, m_Cols);
    for (std::uint64_t i = 0; i < m_Rows * m_Cols; ++i) {
        r_M.m_Data[i] = scaler * m_Data[i];
    }
    return r_M;
}

template <typename T>
jacobi::Vector<T> jacobi::Matrix<T>::operator*(const Vector<T>& i_V) const {
    if (m_Cols != i_V.size()) {
        throw std::domain_error("Attempt to multiply non-compatible matrix-vector.");
    }
    Vector<T> r_V(m_Rows);
    for (std::uint64_t i = 0; i < m_Rows; ++i) {
        r_V[i] = T{};
        for (std::uint64_t j = 0; j < m_Cols; ++j) {
            r_V[i] += m_Data[j + i * m_Rows] * i_V[j];
        }
    }
    return r_V;
}

template <typename T>
jacobi::Matrix<T> jacobi::Matrix<T>::operator/(const T& scaler) const {
    Matrix r_M(m_Rows, m_Cols);
    for (std::uint64_t i = 0; i < m_Rows * m_Cols; ++i) {
        r_M.m_Data[i] = m_Data[i] / scaler;
    }
    return r_M;
}

template <typename T>
jacobi::Matrix<T> jacobi::Matrix<T>::operator/(const Matrix& i_M) const {
    if (m_Cols != i_M.m_Cols || m_Rows != i_M.m_Rows) {
        throw std::length_error("Shape of both matrix must be equal");
    }
    Matrix r_M(m_Rows, i_M.m_Cols);
    for (std::uint64_t i = 0; i < m_Rows * m_Cols; ++i) {
        r_M.m_Data[i] = m_Data[i] / i_M.m_Data[i];
    }
    return r_M;
}

template <typename T>
jacobi::Matrix<T> jacobi::Matrix<T>::hadamard_product(const Matrix& i_M) const {
    if (m_Rows != i_M.rows() || m_Cols != i_M.cols()) {
        throw std::domain_error("Attempt to multiply non-compatible matrices");
    }
    Matrix r_M(m_Rows, m_Cols);
    for (std::uint64_t i = 0; i < m_Rows * m_Cols; ++i) {
        r_M.m_Data[i] = m_Data[i] * i_M.m_Data[i];
    }
    return r_M;
}

template <typename T>
jacobi::Matrix<T> jacobi::Matrix<T>::transpose() const {
    Matrix r_M(m_Cols, m_Rows);
    for (std::uint64_t i = 0; i < m_Cols; ++i) {
        for (std::uint64_t j = 0; j < m_Rows; ++j) {
            r_M(i, j) = m_Data[i + j * m_Cols];
        }
    }
    return r_M;
}

template <typename T>
jacobi::Vector<T> jacobi::Matrix<T>::flatten() const {
    Vector<T> r_V(m_Rows * m_Cols);
    for (std::uint64_t i = 0; i < m_Rows * m_Cols; ++i) {
        r_V[i] = m_Data[i];
    }
    return r_V;
}

template <typename T>
void jacobi::Matrix<T>::save(const std::string& path, const std::uint8_t version) const {
    std::ofstream file(path, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Could not create file " + path);
    }
    file.write(reinterpret_cast<const char*>(&magic_header), header_length);
    file.write(reinterpret_cast<const char*>(&version), sizeof(version));

    file.put('M');
    file.write(reinterpret_cast<const char*>(&m_Rows), sizeof(m_Rows));
    file.write(reinterpret_cast<const char*>(&m_Cols), sizeof(m_Cols));
    file.write(reinterpret_cast<const char*>(m_Data), sizeof(T) * m_Rows * m_Cols);

    file.write(reinterpret_cast<const char*>(&magic_footer), footer_length);
    file.close();
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const jacobi::Matrix<T>& i_M) {
    constexpr int n{10};
    for (std::uint64_t i = 0; i < i_M.rows(); ++i) {
        os << '[' << i_M(i, 0);
        for (std::uint64_t j = 1; j < i_M.cols(); ++j) {
            os << std::setw(n) << i_M(i, j);
        }
        os << "]\n";
    }
    return os;
}

template <typename T>
jacobi::Matrix<T>::~Matrix() {
    delete[] m_Data;
    m_Data = nullptr;
}

template <typename T>
jacobi::Matrix<T> jacobi::operator/(const T& scaler, const Matrix<T>& i_M) {
    Matrix<T> r_M(i_M.rows(), i_M.cols());
    for (std::uint64_t i = 0; i < i_M.size(); ++i) {
        r_M[i] = scaler / i_M[i];
    }
    return r_M;
}

template <typename T>
jacobi::Matrix<T> jacobi::zeros(std::uint64_t rows, std::uint64_t cols) {
    Matrix<T> r_M(rows, cols);
    for (std::uint64_t i = 0; i < rows; ++i) {
        for (std::uint64_t j = 0; j < cols; ++j) {
            r_M(i, j) = T{};
        }
    }
    return r_M;
}

template <typename T>
jacobi::Matrix<T> jacobi::ones(std::uint64_t rows, std::uint64_t cols) {
    Matrix<T> r_M(rows, cols);
    for (std::uint64_t i = 0; i < rows; ++i) {
        for (std::uint64_t j = 0; j < cols; ++j) {
            r_M(i, j) = 1;
        }
    }
    return r_M;
}

template <typename T>
jacobi::Matrix<T> jacobi::eyes(std::uint64_t order) {
    Matrix<T> r_M(order, order);
    for (std::uint64_t i = 0; i < order; ++i) {
        for (std::uint64_t j = 0; j < order; ++j) {
            r_M(i, j) = i == j ? 1 : T{};
        }
    }
    return r_M;
}

template<typename T>
jacobi::Vector<T> jacobi::load_vector(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Could not open file " + path);
    }
    char footer[footer_length];
    file.seekg(-footer_length, std::ios::end);
    file.read(footer, footer_length);
    char header[header_length];
    file.seekg(0, std::ios::beg);
    file.read(header, header_length);
    if (memcmp(header, magic_header, header_length) != 0 || memcmp(footer, magic_footer, footer_length) != 0) {
        throw std::runtime_error("Not a .jcbi file");
    }
    std::uint8_t version;
    file.read(reinterpret_cast<char*>(&version), 1);
    if (version == 1) {
        char type;
        file.read(&type, 1);
        if (type == 'V') {
            std::uint64_t length{};
            file.read(reinterpret_cast<char*>(&length), sizeof(length));
            Vector<T> r_V(length);

            file.read(reinterpret_cast<char*>(&r_V[0]), sizeof(T) * length);
            // for (std::uint64_t i = 0; i < length; ++i) {
            //     file.read(reinterpret_cast<char*>(&r_V[i]), sizeof(T));
            // }
            return r_V;
        }
        if (type == 'M') {
            std::uint64_t rows{};
            file.read(reinterpret_cast<char*>(&rows), sizeof(rows));
            std::uint64_t cols{};
            file.read(reinterpret_cast<char*>(&cols), sizeof(cols));
            Vector<T> r_V(rows * cols);

            file.read(reinterpret_cast<char*>(&r_V[0]), sizeof(T) * rows * cols);
            // for (std::uint64_t i = 0; i < rows * cols; ++i) {
            //     file.read(reinterpret_cast<char*>(&r_V[i]), sizeof(T));
            // }
            return r_V;
        }
        throw std::runtime_error("Invalid file");
    }
    throw std::runtime_error("Version mismatch");
}

template <typename T>
jacobi::Matrix<T> jacobi::load_matrix(const std::string& path) {
    std::ifstream file(path, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Could not open file " + path);
    }
    char footer[footer_length];
    file.seekg(-footer_length, std::ios::end);
    file.read(footer, footer_length);
    char header[header_length];
    file.seekg(0, std::ios::beg);
    file.read(header, header_length);
    if (memcmp(header, magic_header, header_length) != 0 || memcmp(footer, magic_footer, footer_length) != 0) {
        throw std::runtime_error("Not a .jcbi file");
    }
    std::uint8_t version;
    file.read(reinterpret_cast<char*>(&version), 1);
    if (version == 1) {
        char type;
        file.read(&type, 1);
        if (type == 'M') {
            std::uint64_t rows{};
            file.read(reinterpret_cast<char*>(&rows), sizeof(rows));
            std::uint64_t cols{};
            file.read(reinterpret_cast<char*>(&cols), sizeof(cols));
            Matrix<T> r_M(rows, cols);
            for (std::uint64_t i = 0; i < rows * cols; ++i) {
                file.read(reinterpret_cast<char*>(&r_M[i]), sizeof(T));
            }
            return r_M;
        }
        if (type == 'V') {
            std::uint64_t rows{};
            file.read(reinterpret_cast<char*>(&rows), sizeof(rows));
            std::uint64_t cols{1};
            Matrix<T> r_M(rows, cols);
            for (std::uint64_t i = 0; i < rows; ++i) {
                file.read(reinterpret_cast<char*>(&r_M[i]), sizeof(T));
            }
            return r_M;
        }
        throw std::runtime_error("Invalid file");
    }
    throw std::runtime_error("Version mismatch");
}