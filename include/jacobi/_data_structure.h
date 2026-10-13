#pragma once
#include <iostream>
#include <initializer_list>
#include <iomanip>
#include <stdexcept>

namespace jacobi {
    template <typename T>
    class Matrix;

    template <typename T>
    class Vector {
    private:
        std::size_t m_Length;
        T* m_Data;
    public:
        Vector(std::size_t);
        Vector(const Vector&);
        Vector(std::initializer_list<T>);
        Vector& operator=(const Vector&);

        T& operator[](std::size_t);
        const T& operator[](std::size_t) const;

        Vector operator+(const Vector&) const;
        Vector operator+(const T&) const;

        Vector operator-(const Vector&) const;
        Vector operator-() const;
        Vector operator-(const T&) const;

        T inner_product(const Vector&) const;
        T operator*(const Vector& i_V) const { return inner_product(i_V); }
        Vector operator*(const Matrix<T>& i_M) const;

        Vector operator*(const T&) const;

        Vector hadamard_product(const Vector&) const;

        Matrix<T> outer_product(const Vector&) const;

        [[nodiscard]] std::size_t size() const { return m_Length; }

        ~Vector();
    };

    template <typename T>
    Vector<T> operator+(const T& scaler, const Vector<T>& i_V) { return i_V + scaler; }

    template <typename T>
    Vector<T> operator-(const T& scaler, const Vector<T>& i_V) { return (-i_V) + scaler; }

    template <typename T>
    Vector<T> operator*(const T& scaler, const Vector<T>& i_V) { return i_V * scaler; }

    template <typename T>
    Vector<T> zeros(std::size_t);
    template <typename T>
    Vector<T> ones(std::size_t);

    template <typename T>
    class Matrix {
    private:
        std::size_t m_Rows, m_Cols;
        T* m_Data;
    public:
        Matrix(std::size_t, std::size_t);
        Matrix(const Matrix&);
        Matrix(std::initializer_list<std::initializer_list<T>>);
        Matrix(const Vector<T>&);
        Matrix& operator=(const Matrix&);

        T& operator()(const std::size_t row, const std::size_t col) { return m_Data[col + row * m_Cols]; }
        const T& operator()(const std::size_t row, const std::size_t col) const { return m_Data[col + row * m_Cols]; }
        T& operator[](const std::size_t n_element) { return m_Data[n_element]; }
        const T& operator[](const std::size_t n_element) const { return m_Data[n_element]; }

        Matrix operator+(const Matrix&) const;
        Matrix operator+(const T&) const;

        Matrix operator-(const Matrix&) const;
        Matrix operator-() const;
        Matrix operator-(const T&) const;

        Matrix operator*(const Matrix&) const;
        Matrix operator*(const T&) const;
        Vector<T> operator*(const Vector<T>&) const;

        Matrix hadamard_product(const Matrix&) const;

        Matrix transpose() const;

        Vector<T> flatten() const;

        std::size_t size() const { return m_Rows * m_Cols; }
        std::size_t rows() const { return m_Rows; }
        std::size_t cols() const { return m_Cols; }

        ~Matrix();
    };

    template <typename T>
    Matrix<T> operator+(const T& scaler, const Matrix<T>& i_M) { return i_M + scaler; }

    template <typename T>
    Matrix<T> operator-(const T& scaler, const Matrix<T>& i_M) { return (-i_M) + scaler; }

    template <typename T>
    Matrix<T> operator*(const T& scaler, const Matrix<T>& i_M) { return i_M * scaler; }

    template <typename T>
    Matrix<T> zeros(std::size_t, std::size_t);
    template <typename T>
    Matrix<T> ones(std::size_t, std::size_t);
    template <typename T>
    Matrix<T> eyes(std::size_t);
}

template <typename T>
jacobi::Vector<T>::Vector(const std::size_t length)
    : m_Length(length) {
    m_Data = new T[m_Length];
}

template <typename T>
jacobi::Vector<T>::Vector(const Vector& i_V)
    : m_Length(i_V.m_Length) {
    m_Data = new T[m_Length];
    for (std::size_t i = 0; i < m_Length; ++i) {
        m_Data[i] = i_V.m_Data[i];
    }
}

template <typename T>
jacobi::Vector<T>::Vector(const std::initializer_list<T> values)
    : m_Length(values.size()) {
    m_Data = new T[m_Length];
    std::size_t i = 0;
    for (const T& value : values) {
        m_Data[i++] = value;
    }
}

template <typename T>
jacobi::Vector<T>& jacobi::Vector<T>::operator=(const Vector& i_V) {
    if (this != &i_V) {
        T* new_data = new T[m_Length];
        for (std::size_t i = 0; i < m_Length; ++i) {
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
T& jacobi::Vector<T>::operator[](const std::size_t idx) {
    return m_Data[idx];
}

template <typename T>
const T& jacobi::Vector<T>::operator[](const std::size_t idx) const {
    return m_Data[idx];
}

template <typename T>
jacobi::Vector<T> jacobi::Vector<T>::operator+(const Vector& i_V) const {
    if (m_Length != i_V.m_Length) {
        throw std::length_error("Length of both vectors must be equal");
    }
    Vector r_V(m_Length);
    for (std::size_t i = 0; i < m_Length; ++i) {
        r_V.m_Data[i] = m_Data[i] + i_V.m_Data[i];
    }
    return r_V;
}

template <typename T>
jacobi::Vector<T> jacobi::Vector<T>::operator+(const T& scaler) const {
    Vector r_V(m_Length);
    for (std::size_t i = 0; i < m_Length; ++i) {
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
    for (std::size_t i = 0; i < m_Length; ++i) {
        r_V.m_Data[i] = m_Data[i] - i_V.m_Data[i];
    }
    return r_V;
}

template <typename T>
jacobi::Vector<T> jacobi::Vector<T>::operator-() const {
    Vector r_V(m_Length);
    for (std::size_t i = 0; i < m_Length; ++i) {
        r_V.m_Data[i] = -m_Data[i];
    }
    return r_V;
}

template <typename T>
jacobi::Vector<T> jacobi::Vector<T>::operator-(const T& scaler) const {
    Vector r_V(m_Length);
    for (std::size_t i = 0; i < m_Length; ++i) {
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
    for (std::size_t i = 0; i < m_Length; ++i) {
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
    for (std::size_t i = 0; i < i_M.cols(); ++i) {
        r_V[i] = T{};
        for (std::size_t j = 0; j < m_Length; ++j) {
            r_V[i] += m_Data[j] * i_M(j, i);
        }
    }
    return r_V;
}

template <typename T>
jacobi::Vector<T> jacobi::Vector<T>::operator*(const T& scaler) const {
    Vector r_V(m_Length);
    for (std::size_t i = 0; i < m_Length; ++i) {
        r_V.m_Data[i] = m_Data[i] * scaler;
    }
    return r_V;
}

template<typename T>
jacobi::Vector<T> jacobi::Vector<T>::hadamard_product(const Vector& i_V) const {
    if (m_Length != i_V.m_Length) {
        throw std::length_error("Length of both vectors must be equal");
    }
    Vector r_V(m_Length);
    for (std::size_t i = 0; i < m_Length; ++i) {
        r_V.m_Data[i] = m_Data[i] * i_V.m_Data[i];
    }
    return r_V;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const jacobi::Vector<T>& i_V) {
    constexpr int n{10};
    os << '[' << i_V[0];
    for (std::size_t i = 1; i < i_V.size(); ++i) {
        os << std::setw(n) << i_V[i];
    }
    os << "]\n";
    return os;
}

template <typename T>
jacobi::Matrix<T> jacobi::Vector<T>::outer_product(const Vector& i_V) const {
    Matrix<T> r_M(m_Length, i_V.size());
    for (std::size_t i = 0; i < m_Length; ++i) {
        for (std::size_t j = 0; j < i_V.size(); ++j) {
            r_M(i, j) = m_Data[i] * i_V[j];
        }
    }
    return r_M;
}

template <typename T>
jacobi::Vector<T>::~Vector() {
    delete[] m_Data;
    m_Data = nullptr;
}

template <typename T>
jacobi::Vector<T> jacobi::zeros(std::size_t length) {
    Vector<T> r_V(length);
    for (std::size_t i = 0; i < length; ++i) {
        r_V[i] = T{};
    }
    return r_V;
}

template <typename T>
jacobi::Vector<T> jacobi::ones(std::size_t length) {
    Vector<T> r_V(length);
    for (std::size_t i = 0; i < length; ++i) {
        r_V[i] = 1;
    }
    return r_V;
}

template <typename T>
jacobi::Matrix<T>::Matrix(const std::size_t rows, const std::size_t cols)
    : m_Rows(rows), m_Cols(cols) {
    m_Data = new T[m_Rows * m_Cols];
}

template <typename T>
jacobi::Matrix<T>::Matrix(const Matrix& i_M)
    : m_Rows(i_M.m_Rows), m_Cols(i_M.m_Cols) {
    m_Data = new T[m_Rows * m_Cols];
    for (std::size_t i = 0; i < m_Rows * m_Cols; ++i) {
        m_Data[i] = i_M.m_Data[i];
    }
}

template <typename T>
jacobi::Matrix<T>::Matrix(std::initializer_list<std::initializer_list<T>> values)
    : m_Rows(values.size()), m_Cols(values.begin()->size()) {
    m_Data = new T[m_Rows * m_Cols];
    std::size_t i = 0;
    for (const auto& row : values) {
        std::size_t j = 0;
        for (const auto& element : row) {
            m_Data[i * m_Cols + j++] = element;
        }
        ++i;
    }
}

template <typename T>
jacobi::Matrix<T>::Matrix(const Vector<T>& i_V)
    : m_Rows(i_V.size()), m_Cols(1) {
    m_Data = new T[m_Rows * m_Cols];
    for (std::size_t i = 0; i < m_Rows * m_Cols; ++i) {
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
        for (std::size_t i = 0; i < i_M.m_Rows * i_M.m_Cols; ++i) {
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
    for (std::size_t i = 0; i < m_Rows * m_Cols; ++i) {
        r_M.m_Data[i] = m_Data[i] + i_M.m_Data[i];
    }
    return r_M;
}

template<typename T>
jacobi::Matrix<T> jacobi::Matrix<T>::operator+(const T& scaler) const {
    Matrix r_M(m_Rows, m_Cols);
    for (std::size_t i = 0; i < m_Rows * m_Cols; ++i) {
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
    for (std::size_t i = 0; i < m_Rows * m_Cols; ++i) {
        r_M.m_Data[i] = m_Data[i] - i_M.m_Data[i];
    }
    return r_M;
}

template <typename T>
jacobi::Matrix<T> jacobi::Matrix<T>::operator-() const {
    Matrix r_M(m_Rows, m_Cols);
    for (std::size_t i = 0; i < m_Rows * m_Cols; ++i) {
        r_M.m_Data[i] = -m_Data[i];
    }
    return r_M;
}

template<typename T>
jacobi::Matrix<T> jacobi::Matrix<T>::operator-(const T& scaler) const {
    Matrix r_M(m_Rows, m_Cols);
    for (std::size_t i = 0; i < m_Rows * m_Cols; ++i) {
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
    for (std::size_t i = 0; i < m_Rows; ++i) {
        for (std::size_t j = 0; j < i_M.m_Cols; ++j) {
            r_M(i, j) = T{};
            for (std::size_t k = 0; k < m_Cols; ++k) {
                r_M(i, j) += m_Data[i * m_Cols + k] * i_M.m_Data[j + k * m_Rows];
            }
        }
    }
    return r_M;
}

template <typename T>
jacobi::Matrix<T> jacobi::Matrix<T>::operator*(const T& scaler) const {
    Matrix r_M(m_Rows, m_Cols);
    for (std::size_t i = 0; i < m_Rows * m_Cols; ++i) {
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
    for (std::size_t i = 0; i < m_Rows; ++i) {
        r_V[i] = T{};
        for (std::size_t j = 0; j < m_Cols; ++j) {
            r_V[i] += m_Data[j + i * m_Rows] * i_V[j];
        }
    }
    return r_V;
}

template <typename T>
jacobi::Matrix<T> jacobi::Matrix<T>::hadamard_product(const Matrix& i_M) const {
    if (m_Rows != i_M.rows() || m_Cols != i_M.cols()) {
        throw std::domain_error("Attempt to multiply non-compatible matrices");
    }
    Matrix r_M(m_Rows, m_Cols);
    for (std::size_t i = 0; i < m_Rows * m_Cols; ++i) {
        r_M.m_Data[i] = m_Data[i] * i_M.m_Data[i];
    }
    return r_M;
}

template <typename T>
jacobi::Matrix<T> jacobi::Matrix<T>::transpose() const {
    Matrix r_M(m_Cols, m_Rows);
    for (std::size_t i = 0; i < m_Cols; ++i) {
        for (std::size_t j = 0; j < m_Rows; ++j) {
            r_M(i, j) = m_Data[i + j * m_Cols];
        }
    }
    return r_M;
}

template <typename T>
jacobi::Vector<T> jacobi::Matrix<T>::flatten() const {
    Vector<T> r_V(m_Rows * m_Cols);
    for (std::size_t i = 0; i < m_Rows * m_Cols; ++i) {
        r_V[i] = m_Data[i];
    }
    return r_V;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const jacobi::Matrix<T>& i_M) {
    constexpr int n{10};
    for (std::size_t i = 0; i < i_M.rows(); ++i) {
        os << '[' << i_M(i, 0);
        for (std::size_t j = 1; j < i_M.cols(); ++j) {
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
jacobi::Matrix<T> jacobi::zeros(std::size_t rows, std::size_t cols) {
    Matrix<T> r_M(rows, cols);
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            r_M(i, j) = T{};
        }
    }
    return r_M;
}

template <typename T>
jacobi::Matrix<T> jacobi::ones(std::size_t rows, std::size_t cols) {
    Matrix<T> r_M(rows, cols);
    for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            r_M(i, j) = 1;
        }
    }
    return r_M;
}

template <typename T>
jacobi::Matrix<T> jacobi::eyes(std::size_t order) {
    Matrix<T> r_M(order, order);
    for (std::size_t i = 0; i < order; ++i) {
        for (std::size_t j = 0; j < order; ++j) {
            r_M(i, j) = i == j ? 1 : T{};
        }
    }
    return r_M;
}