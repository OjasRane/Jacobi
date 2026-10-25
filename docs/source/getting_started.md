# Getting Started

Test the if the library is linked by executing the following code

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int main(){
    auto mat = jacobi::eyes<int>(3);
    std::cout << mat;     
}
```

The above code should print an identity matrix of order 3 on the console.