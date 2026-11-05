# Matrix

Matrix is a templated data structure in Jacobi.

It is the same concept used in Linear Algebra.

`jacobi::Matrix` uses row-major format, memory is managed on heap with rule of 3 implemented.

**Shape and size of matrix cannot be changed after creating one!**

Various operators are overloaded for ease of use.

## Creating a `jacobi::Matrix`

A matrix can be created using:

1. Providing shape of the matrix.
2. Using initializer list to initialize the matrix.
3. Using `zeros`, `zeros_like`
4. Using `ones`, `ones_like`
5. Using `eyes`, `eyes_like`
6. Using a Vector.
7. Using a `Row` object.

### Providing shape of matrix

To create an uninitialized matrix, pass the desired dimensions as constructor argument and data type as template argument.

Dimensions of a matrix are stored as `std::uint64_t`

```cpp
#include <jacobi/jacobi.h>

int main(){
    jacobi::Matrix<double> mat(10, 10);
}
```

### Using initializer list

This will create an initialized matrix.

```cpp
#include <jacobi/jacobi.h>

int main(){
    jacobi::Matrix mat = {
        { 1, 2, 3 },
        { 3, 4, 5 }    
    };
}
```

### Using `zeros`, `zeros_like`

This will create and initialize all elements of the matrix with 0.

Two `std::uint64_t` are to be provided in function argument for dimensions and typename as template argument for `jacobi::zeros`

Moreover, one can pass another matrix to construct a matrix with same data type and dimensions for `jacobi::zeros_like`.

```cpp
#include <jacobi/jacobi.h>

int main(){
    jacobi::Matrix mat_zeros = jacobi::zeros<int>(3, 3);
    jacobi::Matrix mat_zeros_like = jacobi::zeros_like(mat_zeros);
}
```

### Using `ones`, `ones_like`

This will create and initialize all elements of the matrix with 1.

Two `std::uint64_t` are to be provided in function argument for dimensions and typename as template argument for `jacobi::ones`.

Moreover, one can pass another matrix to construct a matrix with same data type and dimensions for `jacobi::ones_like`.

```cpp
#include <jacobi/jacobi.h>

int main(){
    jacobi::Matrix mat_ones = jacobi::ones<int>(3, 3);
    jacobi::Matrix mat_ones_like = jacobi::ones_like(mat_ones);
}
```

### Using `eyes`, `eyes_like`

This will create an identity matrix of given order.

A `std::uint64_t` is to be provided in function argument and typename as template argument for `jacobi::eyes`

Moreover, one can pass another matrix to construct a matrix with same data type and dimensions for `jacobi::eyes_like`.

```cpp
#include <jacobi/jacobi.h>

int main(){
    jacobi::Matrix mat_ones = jacobi::eyes<int>(3);
    jacobi::Matrix mat_ones_like = jacobi::eyes_like(mat_ones);
}
```

### Using a Vector

One can pass a Vector to construct a matrix.

Note that in Jacobi vectors are considered as columns.

Therefore, if $\mathbf{x}\in\mathbb{R}^n$ then the contructed matrix $X\in\mathbb{R}^{n\times 1}$

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Vector x = { 1, 2, 3 };
    jacobi::Matrix X(x);
    std::cout << X;
}
```

Output

```
[1]
[2]
[3]
```

### Using a `Row` object

A `Row` object is a view into rows of a matrix. Since, it is a row, the resulting matrix constructed from row of a matrix in $\mathbb{R}^{m\times n}$ will belong in $\mathbb{R}^{1\times n}$

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Matrix matA = {
        { 1, 2, 3 },
        { 4, 5, 6 }
    };
    jacobi::Matrix matB(matA.row(0));
    std::cout << matB << std::endl;
}
```

Output

```
[1         2         3]
```

## Displaying a matrix

`operator<<` is overloaded to print matrix to the console. Note that it will also print a trailing new line character.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Matrix mat = {
        { 1, 2, 3 },
        { 4, 5, 6 }
    };
    std::cout << mat;
}
```

Output

```
[1         2         3]
[4         5         6]
```

## Getting size and dimensions of matrix

For a matrix in $\mathbb{R}^{m\times n}$ `size()` method returns $m\times n$, `row()` method returns $m$ and `cols()` method returns $n$.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Matrix mat = {
        { 1, 2, 3 },
        { 4, 5, 6 }
    };
    std::cout << "Size of matrix: " << mat.size() << std::endl;
    std::cout << "Rows in the matrix: " << mat.rows() << std::endl;
    std::cout << "Columns in the matrix: " << mat.cols() << std::endl;
}
```

Output

```
Size of matrix: 6
Rows in the matrix: 2
Columns in the matrix: 3
```

## Accessing elements of the matrix

There are 2 operators overloaded for this:

1. `operator[]`
2. `operator()`

### `operator[]`

Matrix is stored as a single array in row-major format. This operator accepts one `std::uint64_t` as input and that returns the value at that position in the contiguous array.

In simpler words, think of flattening the matrix such that rows are kept next to each other and just accessing one value from that long array of size $mn$.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Matrix mat = {
        { 1, 2, 3 },
        { 4, 5, 6 }
    };
    std::cout << mat[3] << std::endl;
    mat[3] = 10;
    std::cout << mat << std::endl;
}
```

Output

```
4
[1         2         3]
[10         5         6]
```

### `operator()`

This returns a particular element at the position $ij$ by reference.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Matrix mat = {
        { 1, 2, 3 },
        { 4, 5, 6 }
    };
    std::cout << mat(1, 0) << std::endl;
    mat(1, 0) = 10;
    std::cout << mat << std::endl;
}
```

Output

```
4
[1         2         3]
[10         5         6]
```

## `operator+`

The function `operator+` can act differently depending on the input.

### Matrix addition

Two matrices of same dimensions can be added by

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Matrix matA = {
        { 1, 2, 3 },
        { 4, 5, 6 }
    };

    jacobi::Matrix matB = {
        { 7, 8, 9 },
        { 0, 1, 2 }
    };

    std::cout << matA + matB << std::endl;
}
```

Output

```
[8        10        12]
[4         6         8]
```

### Adding scalar to a matrix

Adding scalar to a matrix is not really defined in Linear Algebra, but here in jacobi it performs element-wise addition

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Matrix mat = {
        { 1, 2, 3 },
        { 4, 5, 6 }
    };
    std::cout << mat + 10 << std::endl;
}
```

Output

```
[11        12        13]
[14        15        16]
```

## `operator-`

Works same as `operator+`, just for subtraction.

One can subtract matrix and scalars.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Matrix matA = {
        { 1, 2, 3 },
        { 4, 5, 6 }
    };

    jacobi::Matrix matB = {
        { 7, 8, 9 },
        { 0, 1, 2 }
    };

    std::cout << "matA - matB:\n" << matA - matB << std::endl;
    std::cout << "matA - 10:\n" << matA - 10 << std::endl;
}
```

Output

```
matA - matB:
[-6        -6        -6]
[4         4         4]

matA - 10:
[-9        -8        -7]
[-6        -5        -4]
```

## `operator*`

The function `operator*` can act differently depending on the input.

### Scalar multiplication

Result of scalar multiplication with a matrix can be obtained by

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Matrix mat = {
        { 1, 2, 3 },
        { 4, 5, 6 }
    };

    std::cout << 2 * mat << std::endl;
}
```

Output

```
[2         4         6]
[8        10        12]
```

### Matrix Multiplication

When both inputs for the binary operator `*` are `jacobi::Matrix`, it performs a naive $\mathcal{O}(n^3)$ matrix multiplication.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Matrix matA = {
        { 1, 2, 3 },
        { 4, 5, 6 }
    };
    jacobi::Matrix matB = {
        { 7, 8 },
        { 9, 10 },
        { 11, 12 }
    };

    std::cout << matA * matB << std::endl;
}
```

Output

```
[58        64]
[139       154]
```

### Matrix Vector multiplication

Here the product $Ax$ is evaluated, where $x\in\mathbb{R}^n$ and $A\in\mathbb{R}^{m\times n}.$

The output will be a vector in $\mathbb{R}^m$

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Matrix matA = {
        { 1, 2, 3 },
        { 4, 5, 6 }
    };
    jacobi::Vector vecX = { 1, 2, 3 };

    std::cout << matA * vecX << std::endl;
}
```

Output

```
[14        32]
```

## `operator/`

Operations such as $\frac{A}{\text{Scalar}}$, $\frac{\text{scalar}}{A}$ and $\frac{A}{B}$.

In $\frac{A}{B}$ just element-wise division occurs.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Matrix matA = {
        { 1.0, 2.0, 3.0 },
        { 4.0, 5.0, 6.0 }
    };

    jacobi::Matrix matB = {
        { 1.5, 2.5, 3.5 },
        { 4.5, 5.5, 6.5 }
    };
    
    std::cout << "matA/2:\n" << matA / 2.0 << std::endl;
    std::cout << "2/matA:\n" << 2.0 / matA << std::endl;
    std::cout << "matA/matB:\n" << matA / matB << std::endl;
}
```

Output

```
matA/2:
[0.5         1       1.5]
[2       2.5         3]

2/matA:
[2         1  0.666667]
[0.5       0.4  0.333333]

matA/matB:
[0.666667       0.8  0.857143]
[0.888889  0.909091  0.923077]
```

## Hadamard Product

Hadamard product is element-wise multiplication of a matrix or vector.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Matrix matA = {
        { 1.0, 2.0, 3.0 },
        { 4.0, 5.0, 6.0 }
    };

    jacobi::Matrix matB = {
        { 1.5, 2.5, 3.5 },
        { 4.5, 5.5, 6.5 }
    };
    std::cout << matA.hadamard_product(matB) << std::endl;
}
```

Output

```
[1.5         5      10.5]
[18      27.5        39]
```

## Transpose of a matrix

Transpose of a matrix is an important operation in Linear Algebra.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Matrix matA = {
        { 1.0, 2.0, 3.0 },
        { 4.0, 5.0, 6.0 }
    };
    std::cout << "Transpose of A:\n" << matA.transpose() << std::endl;
}
```

Output

```
Transpose of A:
[1         4]
[2         5]
[3         6]
```

## Flattening a matrix

A matrix can be converted to a Vector by calling the `flatten()` method.

It will convert a matrix in $\mathbb{R}^{m\times n}$ to a vector in $\mathbb{R}^{mn}$ by arranging rows side-by-side.

This generates a copy.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Matrix matA = {
        { 1.0, 2.0, 3.0 },
        { 4.0, 5.0, 6.0 }
    };
    auto vecA = matA.flatten();
    std::cout << "Flattened matrix: " << vecA << std::endl;
}
```

## Accessing a row of matrix

The `row()` method provides a view into a specific row of a Matrix. It does not return a `jacobi::Vector` but rather a different data structure that is not possible for user to construct.

It provides it as reference, so one can edit it whole, use operators like `+=, -=, /=, *=` which returns the same object or `+, -, *, /` which returns a vector.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Matrix matA = {
        { 1.0, 2.0, 3.0 },
        { 4.0, 5.0, 6.0 }
    };
    auto row0 = matA.row(0);
    std::cout << "Row 0: " << row0 << std::endl;
    
    row0[1] = 10.0;
    std::cout << "A:\n" << matA << std::endl;
}
```

Output

```
Row 0: [1         2         3]
A:
[1        10         3]
[4         5         6]
```

## Saving a matrix to disk

For any reason, if needed users can save their matrix on disk, and port it into another language or software.

Jacobi provides a `save` method for Vector and Matrix objects and will output a `*.jcbi` file. For more information refer to the dedicated page on it.

To save a matrix provide a string for file path

```cpp
#include <jacobi/jacobi.h>

int main(){
    jacobi::Matrix mat = {
        { 1, 2, 3 },
        { 4, 5, 6 }
    };
    mat.save("mat.jcbi");
}
```

This will create a `mat.jcbi` file, relative to the compiled binary.

To load it back use `load_matrix` function.

Provide the file name as function argument and data type as template argument.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    auto mat = jacobi::load_matrix<int>("mat.jcbi");
    std::cout << mat << std::endl;
}
```

Output

```
[1         2         3]
[4         5         6]
```