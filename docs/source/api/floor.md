# `floor`

Return floor value of a number ($\floor x$)

Returns a `long long` and during input it will cast the input to a `double`.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    std::cout << "Floor of 2.5: " << jacobi::floor(2.5) << std::endl;
}
```

Output

```
Floor of 2.5: 2
```