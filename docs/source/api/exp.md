# `exp`

Returns the exponent value of a number ($e^x$)

`exp` uses repeated multiplication and Maclaurin series expansion to calculate $e^x$

The part $e^{\floor x}$ is calculated by repeated multiplication and $e^{x-\floor x}$ is calculated using Maclaurin series expansion.

This is done because one can get enough accuracy in lesser number of terms when input is near zero.

The result is obtained by doing $e^{\floor x}\times e^{x-\floor x}$.

Returns a `double` and during input it will cast the input to a `double`.

The iteration stops after the absolute value of difference between the consequent iterations becomes less that `jacobi::eps`, i.e. `1e-12`.

Is this good? **NO**, I will clean this is next release.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    std::cout << "exp(2.5)=" << jacobi::exp(2.5);
}
```

Output

```
exp(2.5)=12.1825
```