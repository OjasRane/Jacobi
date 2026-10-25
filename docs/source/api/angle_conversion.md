# Angle Conversion

Two functions exist for converting from degrees to radians and back.

Returns a `double` and during input it will cast the input to a `double`.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    std::cout << "90 degrees in radians: " << jacobi::deg2rad(90) << std::endl;
    std::cout << "Pi radians to degrees: " << jacobi::rad2deg(jacobi::pi) << std::endl;
}
```

Output
```
90 degrees in radians: 1.5708
Pi radians to degrees: 180
```