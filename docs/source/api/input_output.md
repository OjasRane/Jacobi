# I/O

For any reason, users can save their vectors, matrices to disk and load them when needed.

`jacobi::Matrix` and `jacobi::Vector` both have `save` method to write a `.jcbi` file to disk.

The binary format is very simple to parse and already load utility exists in Jacobi to load the file in program.

As of v0.1.0, endianness has not been taken into consideration.

## Saving a file

To save a file, call the `save` method and pass in the argument of relative file path along with name and the extension `.jcbi`.

```cpp
#include <jacobi/jacobi.h>

int main(){
    jacobi::Vector vec = { 1, 2, 3 };
    jacobi::Matrix mat = {
        { 1, 2, 3 },
        { 4, 5, 6 }
    };
    
    vec.save("int_vec.jcbi");
    mat.save("int_mat.jcbi");
}
```

While loading, it is important to know the type of data during saving, so it could be preferred to save the datatype as prefix or suffix in file name.

## Structure of a file

To load a file in other programming language, one will need to write their own loader. Refer to the below structure for it.

### For Vector

1. Magic Header: `CARL` - 4 bytes, taken from the great mathematician Jacobi's full name - Carl Gustav Jacob Jacobi
2. Version: As of v0.1.0 `1` - 1 byte
3. Data structure identifier: A character is placed to identify if a vector was saved or a matrix. `V` is used for vectors - 1 byte
4. Length of vector: Length of vector is stored - 8 bytes
5. Payload: The actual data from the Vector is stored
6. Magic Footer: `QED` - 3 bytes, taken from Quod Erat Demonstrandum to signify the file has ended.

### For Matrix

1. Magic Header: `CARL` - 4 bytes, taken from the great mathematician Jacobi's full name - Carl Gustav Jacob Jacobi
2. Version: As of v0.1.0 `1` - 1 byte
3. Data structure identifier: A character is placed to identify if a vector was saved or a matrix. `M` is used for matrix - 1 byte
4. Dimensions of matrix: Number of rows and columns are stored - 8 bytes each, so effectively 16 bytes
5. Payload: The actual data from the Matrix is stored
6. Magic Footer: `QED` - 3 bytes, taken from Quod Erat Demonstrandum to signify the file has ended.

## Loading file

For loading a file there are 2 functions:

1. `jacobi::load_vector`: Returns a vector
2. `jacobi::load_matrix`: Returns a matrix

### `jacobi::load_vector`

As the name suggests, it will load vector from a file.

If the saved file contained a matrix, then it will still output a vector with length `row * cols`.

In other words if a matrix in $\mathbb{R}^{m\times n}$ was saved, then using `load_vector`, it will output a vector in $\mathbb{R}^{mn}$.

Input parameter is a string of file path and template parameter is the datatype in which the data is to be inferred.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    auto vec = jacobi::load_vector<int>("int_vec.jcbi");
    auto vec_mat = jacobi::load_vector<int>("int_mat.jcbi");
    
    std::cout << "vec: " << vec << std::endl;
    std::cout << "vec_mat: " << vec_mat << std::endl;
}
```

Output

```
vec: [1         2         3]
vec_mat: [1         2         3         4         5         6]
```

### `jacobi::load_matrix`

As the name suggests, it will load matrix from a file.

If the saved file contained a vector, then it will still output a matrix with shape $\text{Length of vector}\times1$.

In other words if a vector in $\mathbb{R}^{n}$ was saved, then using `load_matrix`, it will output a matrix in $\mathbb{R}^{n\times 1}$.

Input parameter is a string of file path and template parameter is the datatype in which the data is to be inferred.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    auto mat = jacobi::load_matrix<int>("int_mat.jcbi");
    auto mat_vec = jacobi::load_matrix<int>("int_vec.jcbi");
    
    std::cout << "mat:\n" << mat << std::endl;
    std::cout << "mat_vec:\n" << mat_vec << std::endl;
}
```

Output

```
mat:
[1         2         3]
[4         5         6]

mat_vec:
[1]
[2]
[3]
```