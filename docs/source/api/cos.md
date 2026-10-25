# `cos`

Returns cosine value of $x$ ($\cos(x)$). $x$ is in radians.

`cos` uses Maclaurin expansion.

First $x$ is brought between 0 and $\pi$, because it is closer to zero and we can get result with lesser number of terms.

Then appropriate sign is placed.

Returns a `double` and during input it will cast the input to a `double`.

The iteration stops after the absolute value of difference between the consequent iterations becomes less that `jacobi::eps`, i.e. `1e-12`.

Is this good? **NO**, I will clean this is next release.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    std::cout << "cos(1)=" << jacobi::cos(1) << std::endl;
}
```

Output

```
cos(1)=0.540302
```