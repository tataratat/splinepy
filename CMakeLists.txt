cmake_minimum_required(VERSION 3.12.0)

project(splinelibpy)

option(COMPILE_SPLINELIB "Compile splinelib together. If ON, cmake 3.19 or higher will be required." ON)
option(VERBOSE_MAKE "Verbose `make` output" ON)

set(C_FLAGS -fPIC)
set(GNU_FLAGS_R -ffast-math -mtune=native -O3)
set(CMAKE_CXX_STANDARD 20)

if(VERBOSE_MAKE)
    set(CMAKE_VERBOSE_MAKEFILE ON)
endif()

if(COMPILE_SPLINELIB)
    set(GOOGLETEST OFF)
    add_subdirectory(third_party/SplineLib)
endif()

add_subdirectory(third_party/pybind11)

add_subdirectory(src)
