cmake_minimum_required(VERSION 3.12.0)

project(splinepy VERSION 0.0.2 LANGUAGES CXX)

option(SPLINEPY_COMPILE_SPLINELIB "Compile splinelib together." ON)
option(VERBOSE_MAKE "Verbose `make` output." OFF)
option(MINIMAL "Compile a minimal set of splines." OFF)
option(SPLINEPY_COMPILE_PYTHON "Compile python module." ON)
option(SPLINEPY_ENABLE_WARNINGS "Add warning flags" OFF)

set(SPLINEPY_FLAGS -fPIC)
set(SPLINEPY_OPTIMIZATION_FLAGS -O3)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  list(APPEND SPLINEPY_OPTIMIZATION_FLAGS -ffast-math)
  set(SPLINEPY_WARNING_FLAGS
        -Wall -Wextra -Wpedantic -Wzero-as-null-pointer-constant -Wno-unused)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(SPLINEPY_WARNING_FLAGS
        -Wall -Wextra -Wmost -Wextra -Wpedantic -Wunreachable-code
        -Wshadow -Wfloat-equal -Weffc++ -Wno-unused-parameter
        -Wno-unused-variable)
endif()

if(SPLINEPY_ENABLE_WARNINGS)
  set(SPLINEPY_FLAGS ${SPLINEPY_FLAGS} ${SPLINEPY_WARNING_FLAGS}>)
endif()

set(CMAKE_CXX_VISIBILITY_PRESET "hidden")
set(CMAKE_CXX_STANDARD 20)

if(VERBOSE_MAKE)
    set(CMAKE_VERBOSE_MAKEFILE ON)
endif()

add_subdirectory(third_party)
add_subdirectory(cpp/splinepy)
