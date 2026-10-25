# `ceil`

Return ceiling value of a number ($\ceil x$)

Returns a `long long` and during input it will cast the input to a `double`.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    std::cout << "Ceil of 2.5: " << jacobi::ceil(2.5) << std::endl;
}
```

Output

```
Ceil of 2.5: 3
```