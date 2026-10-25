# `arange`

This works same as in NumPy's `np.arange`, except, here one has to provide the start point and step explicitly.

Returns evenly spaced values in a given interval.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    auto time = jacobi::arange<double>(0, 2, 0.1);
    std::cout << time << std::endl;
}
```

Output

```
[0       0.1       0.2       0.3       0.4       0.5       0.6       0.7       0.8       0.9         1       1.1       1.2       1.3       1.4       1.5       1.6       1.7       1.8       1.9]
```