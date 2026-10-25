# `ln`

Returns the natural logarithm of a number ($\ln x$ or $\log_e x$)

`ln` uses Maclaurin expansion of $\ln(1+x)$.

The radius of convergence for this series is $x\in[-1, 1]$ and using same series for $\ln N$ it becomes $N\in [0, 2]$

To calculate natural logarithm for any number, one needs to express the number in the format $a\space B^b$ where $a\in[0, 2]$ and $\ln B$ is known.

One way is to perform repeated divisions, but this will require expensive division and time.

According to IEEE 754 double precision floating-point representation, the number is represented as $(-1)^s 1.M\times 2^{exp}$

One can extract exponent and $1.M$ by performing bit operations directly from memory. Also in this implementation, I have used $(-1)^s \frac{1.M}{2}\times 2^{exp+1}$.

This brings $a\in[0.5, 1)$. Then using property of logarithm, one can calculate natural logarithm.

$\ln(N) = \ln(a\times 2^b) = b\ln 2 + \ln a$ and $\ln a$ is calculuated using Maclaurin expansion.

Returns a `double` and during input it will cast the input to a `double`.

The iteration stops after the absolute value of difference between the consequent iterations becomes less that `jacobi::eps`, i.e. `1e-12`.

Is this good? **NO**, I will clean this is next release.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    std::cout << "ln(3)=" << jacobi::ln(3) << std::endl;
}
```

Output

```
ln(3)=1.09861
```