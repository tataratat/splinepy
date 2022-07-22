cmake_minimum_required(VERSION 3.12.0)

project(splinepy VERSION 0.0.2 LANGUAGES CXX)

option(COMPILE_SPLINELIB "Compile splinelib together." ON)
option(VERBOSE_MAKE "Verbose `make` output." OFF)
option(MINIMAL "Compile a minimal set of splines." OFF)
option(COMPILE_PYTHON "Compile python module." ON)
option(ENABLE_WARNINGS "Add warning flags" OFF)

set(SPLINEPY_FLAGS -fPIC)
if(ENABLE_WARNINGS)
    set(SPLINEPY_FLAGS ${SPLINEPY_FLAGS} -Wall)
endif()

set(SPLINEPY_OPTIMIZATION_FLAGS -O3)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  list(APPEND SPLINEPY_OPTIMIZATION_FLAGS -ffast-math -mtune=native)
endif()
set(CMAKE_CXX_VISIBILITY_PRESET "hidden")
set(CMAKE_CXX_STANDARD 20)

if(VERBOSE_MAKE)
    set(CMAKE_VERBOSE_MAKEFILE ON)
endif()

add_subdirectory(third_party)
add_subdirectory(cpp/splinepy)
