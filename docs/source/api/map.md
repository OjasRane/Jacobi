# `map`

This function maps a function onto a `jacobi::Vector`, `jacobi::Matrix` or `Row` object.

The first argument is a function pointer. For some reason in v0.1.0 lambdas are not working, but other functions work.

Returns the same data structure as input, except for `Row` object, it returns a `jacobi::Vector` object.

```cpp
#include <iostream>
#include <jacobi/jacobi.h>

int add(int x) { return x+1; }

int main(){
    jacobi::Vector vec = { 1, 2, 3 };
    auto mapped_vec = jacobi::map(add, vec);
    std::cout << "add(vec)=" << mapped_vec << std::endl;
}
```

Output

```
add(vec)=[2         3         4]
```