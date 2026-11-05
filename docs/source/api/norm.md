# Norm

Currently only Euclidean norm is supported.

In linear algebra, norm are defined for vectors. If matrix is passed as argument in `jacobi::linalg::norm` then it will take the norm by assuming the matrix as vector.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    jacobi::Vector vec = { 1, 2, 3 };
    jacobi::Matrix mat = {
        { 1, 2, 3 },
        { 4, 5, 6 },
        { 7, 8, 9 }
    };
    
    std::cout << "||vec||=" << jacobi::linalg::norm(vec) << std::endl;
    std::cout << "||mat||=" << jacobi::linalg::norm(mat) << std::endl;
}
```

Output

```
||vec||=3.74166
||mat||=16.8819
```

Here, $||\text{vec}||=\sqrt{1^2+2^2+3^2}$ and $||\text{mat}||=\sqrt{1^2+2^2+3^2+4^2+5^2+6^2+7^2+8^2+9^2}$.