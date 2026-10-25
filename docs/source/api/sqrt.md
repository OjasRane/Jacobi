# `sqrt`

Returns the square root of a number.

`sqrt` uses Newton-Raphson method, with initial guess as the half of the original number.

Returns a `double` and during input it will cast the input to a `double`.

The iteration stops after the absolute value of difference between the consequent iterations becomes less that `jacobi::eps`, i.e. `1e-12`.

Is this good? **NO**, I will clean this is next release.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    std::cout << "Square root of 2: " << jacobi::sqrt(2) << std::endl;
}
```

Output

```
Square root of 2: 1.41421
```