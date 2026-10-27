# Jacobi

A numerical computing library built from scratch in C++20

## Current implementations

Jacobi is under active development. The functionality listed below is implemented from scratch using numerical algorithms.

Jacobi is a work in progress and is not yet intended for production use. For production applications, consider established alternatives such as `<cmath>` for scalar mathematics and Eigen for linear algebra.

- `Vector<T>`: Templated vector.
- `Matrix<T>`: Templated matrix.
- `abs()`: Returns the absolute value.
- `sqrt()`: Returns the square root of a number.
- `sin()`: Returns the sine value of a number.
- `cos()`: Returns the cosine value of a number.
- `exp()`: Returns $e^x$ for a given number $x$.
- `ln()`: Returns the natural logarithm of a number.
- `floor()`: Returns the floor value of a number.
- `ceil()`: Returns the ceiling value of a number.
- `pow()`: Returns $x^y$.
- `deg2rad()`: Converts degrees to radians.
- `rad2deg()`: Converts radians to degrees.
- `norm()`: Returns $\sqrt{\sum{x_i^2}}$ for a vector $\begin{bmatrix}x_1\\x_2\\\vdots\\x_n\end{bmatrix}$ or a flattened matrix.

For more information, refer to the [documentation site](https://ojasrane.github.io/Jacobi).

More functions will be added as the development of the library moves forward.

## Using Jacobi in your projects

```bash
git clone https://github.com/OjasRane/Jacobi.git
cd Jacobi
```

For Linux/MacOS:

```
mkdir build
cd build
cmake -G Ninja -DCMAKE_BUILD_TYPE=Release ..
ninja
mv libJacobi.a /path/to/your/project
cd ..
cp -r include /path/to/your/project
```

For Windows:

```powershell
New-Item -ItemType Directory -Path build
cd build
cmake -G Ninja -DCMAKE_BUILD_TYPE=Release ..
cmake --build .

# Find the compiled binary name from build/ directory and replace it with libJacobi.lib if needed
Move-Item -Path .\libJacobi.lib -Destination C:\path\to\your\project
cd ..
Copy-Item -Path .\include -Destination C:\path\to\your\project -Recurse
```

### Sample code

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main() {
    jacobi::Matrix matA = {
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 9}
    }; // Created a matrix
    
    auto matB = jacobi::ones<int>(3, 3); // Created an all ones matrix of type int
    auto matC = jacobi::eyes<int>(5); // Created an identity matrix of order 5
    
    std::cout << matA << '\n';
    std::cout << matB << '\n';
    std::cout << matC << '\n';
    return 0;
}
```

# Why the name Jacobi?
This library is named "_Jacobi_" in honor of the great mathematician "_Carl Gustav Jacob Jacobi_".

![Carl Gustav Jacob Jacobi](https://upload.wikimedia.org/wikipedia/commons/9/90/Carl_Jacobi.jpg?utm_source=en.wikipedia.org&utm_campaign=imageinfo&utm_content=thumbnail_unscaled)

By Unknown author - http://www.sil.si.edu/digitalcollections/hst/scientific-identity/explore.htm (reworked), Public Domain, https://commons.wikimedia.org/w/index.php?curid=140344

During my days at junior college, I (Ojas Rane) was confused about how the Jacobian was used for 
change of variables in multivariate calculus as I only had knowledge about change of variables in 
single variable calculus. The path to finding the answer led me through the beautiful world of calculus and numerical methods. 
Even after getting the answer, I could not leave numerical methods behind, and my interest in them eventually found its way into this library.