# `linspace`

This works same as in NumPy's `np.linspace`, except, here one has to provide start point explicitly.

Returns `N` points in a given intervals.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    auto time = jacobi::linspace<double>(0, 2, 15);
    std::cout << time << std::endl;
}
```

Output

```
[0  0.142857  0.285714  0.428571  0.571429  0.714286  0.857143         1   1.14286   1.28571   1.42857   1.57143   1.71429   1.85714         2]
```