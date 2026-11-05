# For Windows

Clone the repository

```bash
git clone https://github.com/OjasRane/Jacobi.git
cd Jacobi
```

Build the library

```powershell
New-Item -ItemType Directory -Path build
cd build
cmake -G Ninja ..
ninja
```

A `libJacobi.lib` will be created.

Move `libJacobi.lib` to the project

```powershell
Move-Item -Path .\libJacobi.lib -Destination C:\path\to\your\project
```

Move include headers to the project

```powershell
cd .. # Move to project root
Copy-Item -Path .\include -Destination C:\path\to\your\project -Recurse
```