# Quadric Simplification

A C++/OpenGL project for CS553 Geometric Modeling at Oregon State University containing an implementation of the paper "Surface simplification using quadric error metrics" by Garland and Heckbert:

> Michael Garland and Paul S. Heckbert. 1997. Surface simplification using quadric error metrics. In Proceedings of the 24th annual conference on Computer graphics and interactive techniques (SIGGRAPH '97). ACM Press/Addison-Wesley Publishing Co., USA, 209–216. https://doi.org/10.1145/258734.258849

This project uses CMake as the build system, depends on an installation of Qt5, and has only been tested on Windows with Visual Studio 17. To build and run the project, from the root directory, use

```
mkdir build
cd build
cmake -DQt5_PATH="<path_to_qt_installation>" ..
```

and

```
cmake --build . --config Debug
.\Debug\QuadricSimplification.exe
```
or alternatively
```
cmake --build . --config Release
.\Release\QuadricSimplification.exe
```

The program accepts triangle meshes in the .ply format, which can be found in the `plyfiles` directory.