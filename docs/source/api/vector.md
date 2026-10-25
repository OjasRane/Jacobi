# Vector

Vector is a templated data structure in Jacobi (not to be confused by `std::vector`). It is essentially a 1D array.

Memory for `jacobi::Vector` is managed on heap with rule of 3 implemented.

**Length of vector cannot be changed after creating one!**

Various operators are overloaded for easy of use.

## Creating a `jacobi::Vector`

A vector can be created using:

1. Providing a length of vector.
2. Using initializer list to create and initialize the vector.
3. Using `jacobi::zeros`
4. Using `jacobi::ones`
5. Using `Row` object.

### Using a length to create vector

To create an uninitialized vector, pass the desired length as constructor argument and data type as template argument.

Length of vector is stored as `std::uint64_t`

```cpp
#include <jacobi/jacobi.h>

int main(){
    jacobi::Vector<double> vec(10);    
}
```

### Using initializer list

This will create an initialized vector.

```cpp
#include <jacobi/jacobi.h>

int main(){
    jacobi::Vector vec = { 1, 2, 3 };
}
```

### Using `jacobi::zeros`

This will create and initialize all elements of the vector with 0.

A `std::uint64_t` is provided in function argument for length and typename as template argument.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    auto vec = jacobi::zeros<int>(10);
    std::cout << vec;
}
```

Output

```
[0         0         0         0         0         0         0         0         0         0]
```

### Using `jacobi::ones`

This will create and initialize all elements of the vector with 1.

A `std::uint64_t` is provided in function argument for length and typename as template argument.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    auto vec = jacobi::ones<int>(10);
    std::cout << vec;
}
```

Output

```
[1         1         1         1         1         1         1         1         1         1]
```

### Using `Row` object

`Row` object is a view into rows of a matrix (`jacobi::Matrix`).

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Matrix mat = {
        { 1.0, 2.0, 3.0 },
        { 4.0, 5.0, 6.0 }
    };
    
    jacobi::Vector vec = mat.row(0);
    std::cout << vec << std::endl;
}
```

Output

```
[1         2         3]
```

## Displaying a vector

`operator<<` is overloaded to print a vector to the console.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Vector vec = { 1, 2, 3, 4 };
    std::cout << vec;
}
```
Output

```
[1         2         3         4]
```

## Getting size of vector

To get the size of a created vector, there is `size()` method.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Vector vec = { 1, 2, 3, 4 };
    std::cout << vec.size();
}
```
This will print `4` to the console.

## Accessing elements of Vector

To access elements of vector, `operator[]` is used. It returns the element as reference, so the data can be changed.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Vector vec = { 1, 2, 3, 4 };
    vec[2] = 5;
    std::cout << vec;
}
```

Output

```
[1         2         5         4]
```

## `operator+`

The function `operator+` can act differently depending on the input.

### Vector addition

Two vectors of same size can be added by

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Vector vecA = { 1.0, 2.0, 3.0 };
    jacobi::Vector vecB = { 1.0, 2.0, 3.0 };
    
    std::cout << vecA + vecB;
}
```

Output

```
[2         4         6]
```

Here the another vector can also be a `Row`. `Row` is a row from a matrix.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Vector vecA = { 1.0, 2.0, 3.0 };
    auto mat = jacobi::eyes<double>(3);
    
    std::cout << vecA + mat.row(1);
}
```

Output

```
[1         3         3]
```

### Adding scaler to a vector

Adding scaler with vector is not really standard in linear algebra, but here in Jacobi, it will perform element wise addition.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Vector vec = { 1, 2, 3 };
    
    std::cout << vec + 1;
}
```

Output

```
[2         3         4]
```

## `operator-`

Works same as `operator+`, just for subtraction.

One can subtract vectors, scalers and `Row` objects.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Vector vec1 = { 1, 2, 3 };
    jacobi::Vector vec2 = { 2, 5, 6 };
    jacobi::Matrix mat = {
            { 1, 2, 3 },
            { 4, 5, 6 }
        };
    
    std::cout << "vec1 - vec2: " << vec1 - vec2 << std::endl;
    std::cout << "vec1 - 1: " << vec1 - 1 << std::endl;
    std::cout << "vec1 - mat.row(0): " << vec1 - mat.row(0) << std::endl;
}
```

Output

```
vec1 - vec2: [-1        -3        -3]
vec1 - 1: [0         1         2]
vec1 - mat.row(0): [0         0         0]
```

## `operator*`

The function `operator*` can act differently depending on the input.

### Scaler multiplication

Result of scaler multiplication with a vector can be obtained by

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Vector vecA = {1.0, 2.1, 3.2};
    auto vecResult = 2.5 * vecA;
}
```

### Inner product

Inner product of 2 vectors is defined as $x^Tx$, where $x$ is a column vector. The output is a scaler.

We can find inner product by using `operator*`

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Vector vecA = {1.0, 2.1, 3.2};
    jacobi::Vector vecB = {3.5, 3.5, 2.5};
    std::cout << "Inner product is: " << vecA * vecB;
}
```

Output

```
Inner product is: 18.85
```

Additionally one can use the method `inner_product` to calculuate inner product.

Suppose $x_1, x_2$ are 2 vectors of same shape, then $x_1^Tx_2$ can be calculated by

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Vector x1 = { 1.0, 2.0, 3.0 };
    jacobi::Vector x2 = { 2.0, 3.0, 4.0 };
    
    std::cout << "Inner product is: " << x1.inner_product(x2);
}
```

Output

```
Inner product is: 20
```

### Vector Matrix multiplication

Here the product $x^TA$ is evaluated, where $x\in\mathbb{R}^m$ and $A\in\mathbb{R}^{m\times n}.$

The output will be a vector in $\mathbb{R}^{n}$

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Vector xT = { 1.0, 2.0, 3.0 };
    jacobi::Matrix A = {
            { 2.5, 3.5 },
            { 5.2, 1.2 },
            { 2.3, 5.6, 7.8 }
        };
    
    std::cout << xT * A;
}
```

Output

```
[19.8      22.7         0]
```

## `operator/`

Operations such as $\frac{\mathbf{v}}{\text{Scaler}}$ and $\frac{\text{scaler}}{\mathbf{v}}$

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Vector vec = { 1.0, 2.0, 3.0 };
    
    std::cout << "vec/2: " << vec / 2.0 << std::endl;
    std::cout << "2/vec: " << 2.0 / vec << std::endl;
}
```

Output

```
vec/2: [0.5         1       1.5]
2/vec: [2         1  0.666667]
```

## Outer product

Outer product of 2 vectors $x_1,x_2$ is defined as $x_1x_2^T$

If $x_1\in\mathbb{R}^m$ and $x_2\in\mathbb{R}^n$ then $x_1x_2^T\in\mathbb{R}^{m\times n}$

This outputs a `jacobi::Matrix` object.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Vector x1 = { 1.0, 2.0, 1.0 };
    jacobi::Vector x2 = { 4.0, 2.0, 4.0 };
    
    std::cout << "Outer Product:\n" << x1.outer_product(x2);
}
```

Output

```
Outer Product:
[4         2         4]
[8         4         8]
[4         2         4]
```

## Hadamard Product

Hadamard product is element-wise multiplication of a matrix or vector.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Vector x1 = { 1.0, 2.0, 1.0 };
    jacobi::Vector x2 = { 4.0, 2.0, 4.0 };
    
    std::cout << "Hadamard Product: " << x1.hadamard_product(x2);
}
```

Output

```
Hadamard Product: [4         4         4]
```

## Saving a vector to disk

For any reason, if needed users can save their vector on disk, and port it into another language or software.

Jacobi provides a `save` method for Vector and Matrix objects and will output a `*.jcbi` file. For more information refer to the dedicated page on it.

To save a vector provide a string for file path

```cpp
#include <jacobi/jacobi.h>

int main(){
    jacobi::Vector vec = { 1.0, 2.0, 3.0, 4.0, 5.0 };
    vec.save("vec.jcbi");
}
```

This will create a `vec.jcbi` file, relative to the compiled binary.

To load it back use `load_vector` function.

Provide the file name as function argument and data type as template argument.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    auto vec = jacobi::load_vector<double>("vec.jcbi");
    std::cout << vec << std::endl;
}
```

Output

```
[1         2         3         4         5]
```