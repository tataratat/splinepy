cmake_minimum_required(VERSION 3.12.0)

project(splinelibpy)

option(COMPILE_SPLINELIB "Compile splinelib together. If ON, cmake 3.16 or higher will be required." ON)
option(VERBOSE_MAKE "Verbose `make` output" OFF)

set(SPLINELIBPY_FLAGS -fPIC)
set(GNU_FLAGS_R -ffast-math -mtune=native -O3)
set(CMAKE_CXX_STANDARD 20)

if(VERBOSE_MAKE)
    set(CMAKE_VERBOSE_MAKEFILE ON)
endif()

add_subdirectory(third_party)
add_subdirectory(src)
