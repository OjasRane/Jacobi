# `abs`

Returns the absolute value of a function.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    int a{-10};
    int abs_a = jacobi::abs(a);
    std::cout << "|a|=" << abs_a << std::endl;
}
```

Output

```
|a|=10
```