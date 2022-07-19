# BEZMAN
Light weight library in C++ for functional composition for Bézier splines.

This small prototyping library is created to test out algorithms to create and modify small microstructures and analytically derive their control points. The objective is the analytical calculation of gradients of microstructures formed through functional composition between splines


## Installation
The Project uses cmake, to facilitate integration in other projects (I recommend cmake version 3.16++). *Project is tested on gcc `10.3.xx` and has also been tested on clang `11.1.xx` and `12.0.xx`.*

Go to your destinated installation directory and create a build and install directory (can also be inside the bezman folder). Step into the build directory.
```
mkdir build install
cd build
```

Run CMake specifying an installation directory
```
cmake -DCMAKE_INSTALL_PREFIX=<link-to-your-install-directory> ..
```
Default build type is debug but if you aim for performace, choose Release instead, which sets compiler optimization flags. GTEST is also activated by default, however it is fetched from the internet, no need to have it preinstalled (can be turned off using -DGOOGLETEST=OFF).

Now install the library by running the command
```
make install
```

That's it. There are no external dependencies, however the `c++17` standard must be supported and new compilers are recommended. To run the unit tests, execute `ctest` (e.g. with the verbose option).

| :grey_question: Changing the standard to `c++20` :grey_question: |
|:---------------------------|
| `std::vector` types are used frequently to store all kinds of information (e.g. control points, Spline groups, etc.). The `c++17` standard library prohibits its use at compile time. The new standard allows for these operations|

## Building an example
There are also some simple examples provided, that show the usage of the library in a bit more detail. To build one of them go inside the example directory. It is recommended to not build directly inside the directory itself, but to provide an additional build folder. To do so run:
```
mkdir build
cd build
```
Run CMake and build the executable
```
cmake ..
make
```
The example files also feature a display script, which is based on the [gustav](https://github.com/tataratat/gustaf) library.

## Building the documentation
If you want to build the documentation, you can do so using [doxygen](https://www.doxygen.nl/index.html), by running these commands:
```
# Go to the documentation directory
cd doc
# Build Documentation
doxygen Doxyfile
```
This will create a new folder named `doxydocs` in the current directory, in which you will find an `index.html` file. Open it with a browser of your choice.

## Python integration
Most of the functionality provided by this library is integrated into the software suite [gustaf](https://github.com/tataratat/gustaf). `Gustaf` aims to provide user-friendly of splines and their visualization and animation, for post-processing and presentation purposes.
