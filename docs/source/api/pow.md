# `pow`

Return $x^y$ for inputs $x, y$.

`pow` uses the identity $x^y = e^{y\ln x}$

Returns a `double` and during input it will cast the input to a `double`.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    std::cout << "2^0.5=" << jacobi::pow(2, 0.5) << std::endl;
}
```

Output

```
2^0.5=1.41421
```