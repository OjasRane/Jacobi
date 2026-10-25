# For Linux or MacOS

Clone the repository

```bash
git clone https://github.com/OjasRane/Jacobi.git
cd Jacobi
```

Build the library

```bash
mkdir build
cd build
cmake -G Ninja ..
ninja
```

A `libJacobi.a` will be created.

Move `libJacobi.a` to the project

```bash
cp libJacobi.a /path/to/your/project/
```

Move include headers to the project

```bash
cd .. # Move to project root
cp -r include /path/to/your/project/
```