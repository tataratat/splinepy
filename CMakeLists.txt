cmake_minimum_required(VERSION 3.12.0)

project(splinepy)

option(COMPILE_SPLINELIB "Compile splinelib together." ON)
option(VERBOSE_MAKE "Verbose `make` output." OFF)
option(MINIMAL "Compile a minimal set of splines." OFF)

# try to minimize runtime dependencies
set(CMAKE_CXX_FLAGS "-static-libstdc++ -static-libgcc")

set(SPLINELIBPY_FLAGS -fPIC)
set(GNU_FLAGS_R -ffast-math -mtune=native -O3)
set(CMAKE_CXX_VISIBILITY_PRESET "hidden")
set(CMAKE_CXX_STANDARD 20)

if(VERBOSE_MAKE)
    set(CMAKE_VERBOSE_MAKEFILE ON)
endif()

add_subdirectory(third_party)
add_subdirectory(src)
